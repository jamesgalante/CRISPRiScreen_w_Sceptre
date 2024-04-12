
message("Saving Image")
save.image("split_guide_file.rda")
# stop()


# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages({
  library(tidyverse)
})



# Retrieve parameters from Snakemake
input_file <- snakemake@input$gene_gRNA_group_pairs
batches <- as.numeric(snakemake@params$batches)
output_files <- snakemake@output

# Read the input dataframe
gene_gRNA_group_pairs <- read.csv(input_file, sep="\t", header=TRUE)
gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>%
  dplyr::filter(target_type == "enh")

# Your logic to evenly distribute grna_group across batches
# Calculate the counts of each grna_group
group_counts <- gene_gRNA_group_pairs %>%
  dplyr::count(grna_group) %>%
  dplyr::arrange(desc(n))

# Calculate the target number of unique grna_group for each split
total_unique_groups <- dplyr::n_distinct(gene_gRNA_group_pairs$grna_group)
target_per_split <- ceiling(total_unique_groups / batches)

# Initialize splits
splits <- vector("list", batches)
names(splits) <- paste0("split", seq_len(batches))

# Distribute grna_group to splits trying to even out the number of unique values
for (i in seq_len(nrow(group_counts))) {
  split_counts <- sapply(splits, function(x) dplyr::n_distinct(x$grna_group, na.rm = TRUE))
  split_with_least <- which.min(split_counts)
  
  splits[[split_with_least]] <- dplyr::bind_rows(splits[[split_with_least]], 
                                          gene_gRNA_group_pairs %>%
                                            dplyr::filter(grna_group == group_counts$grna_group[i]))
}

# Write each split to the corresponding output file
for (i in seq_along(splits)) {
  write_tsv(splits[[i]], file=output_files[[i]], col_names = FALSE)
}



# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
