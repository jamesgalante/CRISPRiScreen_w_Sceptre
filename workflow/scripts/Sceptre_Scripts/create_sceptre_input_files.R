## Create the input for the sceptre diff expr and power analysis pipelines from the given files

save.image("create_sceptre_input_files.rda")
# message("Saved Image")
# stop()

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
message("Loading packages")
suppressPackageStartupMessages({
  library(DESeq2)
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyverse)
})

# Load the input files
power_analysis_output <- read_tsv(snakemake@input$power_analysis_output)
guide_targets_w_NT <- read_tsv(snakemake@input$guide_targets_w_NT)



# Create the gene_gRNA_group_pairs file
gene_gRNA_group_pairs <- power_analysis_output[, c("gene", "perturbation", "target_type")]
colnames(gene_gRNA_group_pairs) <- c("response_id", "grna_group", "target_type")



# Create the gRNA_groups_table file
gRNA_groups_table <- guide_targets_w_NT[, c("name", "target_name", "target_type")]
colnames(gRNA_groups_table) <- c("grna_id", "grna_group", "gRNA_type")

# For all rows in `gRNA_groups_table` where the `gRNA_type` is not "enh" or snakemake@params$tss_ctrl_label,
# Make the grna_group == "non-targeting" and the gRNA_type == "negative_targeting"
tss_ctrl_label <- snakemake@params$tss_ctrl_label

# Identify rows where gRNA_type is NA, not "enh", or not tss_ctrl_label
rows_to_update <- which(is.na(gRNA_groups_table$gRNA_type) | 
                          (gRNA_groups_table$gRNA_type != "enh" & 
                             gRNA_groups_table$gRNA_type != tss_ctrl_label))

# Update grna_group and gRNA_type for those rows
gRNA_groups_table$grna_group[rows_to_update] <- "non-targeting"
gRNA_groups_table$gRNA_type[rows_to_update] <- "negative_targeting"


# save perturb_sce with dispersion estimates
message("Writing gene_gRNA_group_pairs.txt to file...")
write_tsv(gene_gRNA_group_pairs, file = snakemake@output[[1]])
message("Done!")
message("Writing gRNA_groups_table.txt to file...")
write_tsv(gRNA_groups_table, file = snakemake@output[[2]])
message("Done!")

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)