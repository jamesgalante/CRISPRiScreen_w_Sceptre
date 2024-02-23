message("Saving Image")
save.image("format_sceptre_output.rda")
# stop()

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# Install packages
message("Import packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(rtracklayer)
})



# Get input
message("LOADING INPUTS ------------------------------------------")
message("Loading in snakemake inputs")
unformatted_power_analysis_results <- snakemake@input$power_analysis_output
effect_size <- snakemake@params$effect_size

# Create the output directory for other outputs
message("Creating output directory for extra files")
full_path <- snakemake@output[[1]]
base_directory <- dirname(full_path)
outdir <- file.path(base_directory, "formatting_script_figures")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Helper function to save plots
message("Creating helper function for saving plots")
save_plot <- function(plot, filename) {
  ggsave(filename = file.path(outdir, filename), plot = plot, width = 10, height = 8)
}







# There were some operations that dealt with different effect sizes, but we're just dealing with one at a time here.
# Read in the results data frame and do some preliminary filtering
message("CREATING PLOTS ----------------------------------------------")
message("Reading in the results from the power analysis")
results <-  read_tsv(unformatted_power_analysis_results)
results_filt <- results[results$response_id != "response_id",]
results_filt$name <- paste0(results_filt$grna_target,"_",results_filt$response_id)

# Filter results
message("Filtering power analysis results")
power <- results_filt %>%
  filter(!is.na(log_2_fold_change)) %>%
  group_by(rep) %>%
  mutate(pval_adj = p.adjust(p_value, method = "BH")) %>%
  group_by(grna_target, response_id) %>%
  summarize(power = mean(pval_adj < 0.1 & log_2_fold_change < 0),
            .groups = "drop") %>%
  arrange(desc(power))

# Assign effect size to results
power$effect_size = effect_size

# Order the power data frame and add some more columns for investigating power analysis results
message("More power analysis filtering")
power <- power %>% 
  arrange(desc(power)) %>% 
  mutate(power_rank = seq_len(n()),
         power_pct = power_rank / n())

# Plot the Power x Tested Pairs in order of decrease power
message("Plotting the Power x Tested Pairs in order of decrease power")
p1 <- ggplot(power, aes(x = power_pct, y = power, color = as.factor(effect_size))) +
  geom_hline(yintercept = 0.8, lty = "dashed") +
  geom_line(lwd = 1) +
  labs(title = "Power across tested pairs", y = "Power", x = "Tested perturbation-gene pairs",
       color = "Effect size") +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  theme(text = element_text(size = 13))

# Filter for pairs that have 80% power and plot
message("Plotting pairs that have 80% power and plot")
power_80 <- filter(power, power >= 0.8)
p2 <- ggplot(power_80, aes(x = factor(effect_size), fill = as.factor(effect_size))) +
  geom_bar() +
  scale_y_continuous(limits = c(0, nrow(power))) +
  labs(x = "Effect size", y = "Tested perturbation-gene pairs", title = paste0("80% power\n",round(nrow(power_80)/nrow(power),3))) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 13))

# Save plots
message("Saving plots")
save_plot(p1, "power_summary.png")
save_plot(p2, "power_summary_barplot.png")




# Let's filter for unique results and combine those with the power table
message("Filter for unique results and merge")
unique_results_filt <-  unique(results_filt[, c("response_id","grna_target","n_nonzero_trt","disp_outlier_deseq2","average_expression_all_cells")])
power <- merge(unique_results_filt, power)

# Create the `avg_expr` and `cells` columns
message("Create variables for plots")
power$avg_expr <- as.numeric(paste0(power$average_expression_all_cells))
power$cells <- as.numeric(paste0(power$n_nonzero_trt))

# Add columns for power quantiles
power <- power %>% 
  mutate(cells_quant = cut(cells, breaks = quantile(cells,probs = seq(0,1,0.2)), include.lowest = TRUE)) %>% 
  mutate(expr_quant = cut(avg_expr, breaks = quantile(avg_expr, probs = seq(0, 1, 0.1)),
                          include.lowest = TRUE))

# compute average transcripts per expression quantile
avg_expr_quants <- power %>% 
  group_by(expr_quant) %>% 
  summarize(mean_expr = mean(avg_expr))

# compute average power and expression per quantile
binned_power <- power %>%
  group_by(effect_size, cells_quant, expr_quant) %>% 
  summarize(mean_power = mean(power), .groups = "drop") %>% 
  left_join(avg_expr_quants, by = "expr_quant")

# plot binned power
message("Plotting binned power")
p3 <- ggplot(binned_power, 
             aes(x = mean_expr, y = mean_power, color = effect_size, linetype = cells_quant)) +
  geom_hline(yintercept = 0.8, lty = "dashed") +
  geom_line() +
  geom_point() +
  labs(title = "Power vs. gene expression", x = "Average UMI/cell per gene\n(binned in decentiles)",
       y = "Power", linetype = "Pert. cells", color = "Effect size") +
  scale_x_log10() +
  theme_bw() +
  theme(text = element_text(size = 13))

# Save plot with helper function defined above
message("Saving plot")
save_plot(p3, "power_v_expression.png")





# Create a dataframe to combine with the differential expression and power analysis results
# Read in the snakemake inputs
message("CREATING DIST FILE -----------------------------------------------------")
message("Reading in gene_gRNA_group_pairs file")
gene_gRNA_group_pairs <- read_tsv(snakemake@input$gene_gRNA_group_pairs)

# Load gene annotation file
message("Import annotation file")
annot <- import(snakemake@params$annot)

# Parse 'grna_group' to extract 'target_chr', 'target_start', 'target_end', and 'target_type
message("filtering gene_gRNA_group_pairs")
gene_gRNA_group_pairs_filt <- gene_gRNA_group_pairs %>%
  mutate(
    target_chr = sub("^(.*):.*$", "\\1", grna_group),
    target_start = as.numeric(sub("^.*:(.*)-.*$", "\\1", grna_group)),
    target_end = as.numeric(sub("^.*-(.*)$", "\\1", grna_group)),
  ) %>%
  dplyr::rename(gene = response_id)

# Remove all annotations that aren't genes, and that are pseudo genes
message("Filtering annotation file")
annot_filt <- annot[annot$type == "gene",]
annot_filt <- annot_filt[!grepl("PAR_Y", annot_filt$gene_id), ]

# Convert 'annot' to a data frame for merging if it's not already
annot_df <- data.frame(
  gene_id = sub("\\..*$", "", annot_filt$gene_id),
  gene_chr = seqnames(annot_filt),
  gene_tss = start(annot_filt)
)

# Merge to align gene information
message("Creating dist_file")
dist_file <- left_join(gene_gRNA_group_pairs_filt, annot_df, by = c("gene" = "gene_id"))

# Calculate 'dist_to_tss'
dist_file <- dist_file %>%
  mutate(dist_to_tss = abs(gene_tss - ((target_start + target_end) / 2)))






# Combine dist_file, discovery_file, and power
message("CREATING FINAL OUTPUT ----------------------------------------------")
message("Loading discovery results")
discovery_results <- read_tsv(snakemake@input$discovery_results)

# Filter the `power` df and change the name of the `power` column, so that it's compatible with ENCODE pipeline
message("Filtering power")
power_filt <-  power[,c("response_id","grna_target", "disp_outlier_deseq2","average_expression_all_cells", "power","cells")]
colnames(power_filt)[colnames(power_filt) == "power"] <- paste0("power_effect_size_", effect_size*100)

# Merge the `discovery_results` with `power_filt`
message("Merging discovery results and power")
discovery_w_power <-  merge(discovery_results, power_filt)

# Format `discovery_w_power`, add missing columns, and rename columns
message("Formatting")
discovery_w_power$gene_strand <- "."
discovery_w_power$pert_level <- "cre_perts"
discovery_w_power$logFC <- log(2^as.numeric(paste0(discovery_w_power$log_2_fold_change)))
colnames(discovery_w_power)[colnames(discovery_w_power) == "grna_target"] <- "grna_group"
colnames(discovery_w_power)[colnames(discovery_w_power) == "response_id"] <- "gene"

# Merge the `discovery_w_power` with the `dist_file`
message("Merging discovery_w_power and dist_file")
discovery_w_power_w_dist <- merge(discovery_w_power, dist_file)
colnames(discovery_w_power_w_dist)[colnames(discovery_w_power_w_dist) == "grna_group"] <- "perturbation"
colnames(discovery_w_power_w_dist)[colnames(discovery_w_power_w_dist) == "p_value"] <- "pvalue"
colnames(discovery_w_power_w_dist)[colnames(discovery_w_power_w_dist) == "target_chr"] <- "pert_chr"
colnames(discovery_w_power_w_dist)[colnames(discovery_w_power_w_dist) == "target_start"] <- "pert_start"
colnames(discovery_w_power_w_dist)[colnames(discovery_w_power_w_dist) == "target_end"] <- "pert_end"

# Calculate the adjusted p_value, and create empty columns that are required for downstream processing
discovery_w_power_w_dist$pval_adj <- p.adjust(discovery_w_power_w_dist$pvalue, method = "BH")
discovery_w_power_w_dist$ci_low <- NA
discovery_w_power_w_dist$ci_high <- NA
discovery_w_power_w_dist$avg_expr <- NA

message("Subsetting needed columns for final output")
formatted_output <- discovery_w_power_w_dist[,c("perturbation", "gene", "logFC", "ci_high", "ci_low", "pvalue", "pval_adj", 
                                                "pert_chr", "pert_start", "pert_end", "gene_chr", "gene_tss", "gene_strand", 
                                                "dist_to_tss", "pert_level", "target_type", "cells", "avg_expr", 
                                                "disp_outlier_deseq2", paste0("power_effect_size_", effect_size*100))]

# Save output
message("Saving output")
write_tsv(formatted_output, snakemake@output[[1]])
message("Output saved")



# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)