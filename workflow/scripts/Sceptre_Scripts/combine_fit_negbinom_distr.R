## Combine per-chromosome negative binomial distr SingleCellExperiments into one SingleCellExperiment

save.image("combine_fit_negbimom_distr.rda")
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
})



### Load in Data ---------------------------------------------------------------

# Read in all of the single cell objects that have dispersion estimates
message("Reading in chromosome RDS objects")
rds_objects <- lapply(snakemake@input$chr_rds_objects, FUN = readRDS)

# Load in the sce object
message("Loading in full SCE object")
sce <- readRDS(snakemake@input$perturb_sce)


### Retrieve Dispersions -------------------------------------------------------

# Iterate over the RDS objects to update the SCE object
message("Adding dispersion estimates to SCE object")
# Extract relevant values from rds_objects
sep_dipsersion_values <- lapply(rds_objects, function(x) {
  # Extracting the desired data
  data <- data.frame(
    gene_id = rownames(x), # Assuming rownames are accessible like this
    mean = rowData(x)$mean,
    dispersion = rowData(x)$dispersion,
    disp_outlier_deseq2 = rowData(x)$disp_outlier_deseq2
  )
  return(data)
})

# Combining all dataframes into one
all_dispersion_values <- do.call(rbind, sep_dipsersion_values)

### Add Dispersions to SCE -----------------------------------------------------

# Initialize columns to be added
rowData(sce)$mean <- NA
rowData(sce)$dispersion <- NA
rowData(sce)$disp_outlier_deseq2 <- NA

# Convert rowData(sce) to a data.frame for manipulation
sce_row_data <- as.data.frame(rowData(sce))

# Find matching indices between SCE rownames and combined data gene_id
match_indices <- match(rownames(sce), all_dispersion_values$gene_id)

# Only update rows with matching gene_id
non_na_indices <- !is.na(match_indices)

# Update rowData for matched genes
sce_row_data$mean[non_na_indices] <- all_dispersion_values$mean[match_indices[non_na_indices]]
sce_row_data$dispersion[non_na_indices] <- all_dispersion_values$dispersion[match_indices[non_na_indices]]
sce_row_data$disp_outlier_deseq2[non_na_indices] <- all_dispersion_values$disp_outlier_deseq2[match_indices[non_na_indices]]

# Convert back to DataFrame to update rowData(sce)
rowData(sce) <- DataFrame(sce_row_data, row.names = rownames(sce))



### Add other info to SCE ------------------------------------------------------

# We also must recompute "size_factors" in colData as this isn't in the main sce object
message("Creating DESeq2 object")
dds <- DESeqDataSetFromMatrix(countData = assay(sce, "counts"),
                              colData = colData(sce),
                              design = ~ 1)

# compute size factors
message("Computing size factors")
if (snakemake@params$size_factors == "libsize") {
  total_umis <- colSums(assay(dds, "counts"))
  manual_size_factors <- total_umis / mean(total_umis)
  sizeFactors(dds) <- manual_size_factors
} else {
  dds <- estimateSizeFactors(dds, type = snakemake@params$size_factors)
}
# Add size factors to sce object colData
colData(sce)[, "size_factors"] <- sizeFactors(dds)

### Save SCE -------------------------------------------------------------------

# save perturb_sce with dispersion estimates
message("Writing SCE object with dispersion estimates to file...")
saveRDS(sce, file = snakemake@output[[1]])
message("Done!")

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)



