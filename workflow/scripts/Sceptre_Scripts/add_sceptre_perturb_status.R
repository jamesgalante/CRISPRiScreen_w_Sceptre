save.image("add_sceptre_perturb_status.rda")
# message("Saved Image")
# stop()

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(Matrix)
})

# Load Sceptre
message("Loading Sceptre")
devtools::install_github("katsevich-lab/sceptre")
library(sceptre)



# Load in snakemake objects
sceptre_object <- readRDS(snakemake@input$sceptre_object)
sce <- readRDS(snakemake@input$sce)



### Let's add the grna_perts and cre_perts assays to the sce object given the sceptre_assignments
# Add the individual grna perts first
message("Adding the individual grna_perts")
get_grna_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@initial_grna_assignment_list)
}

individual_grna_assignments <- get_grna_assignments(sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_grna_assignments)
nCols <- length(colnames(sce))

# Initialize an empty sparse matrix
sparseMatrix <- Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)

# Assuming sparseMatrix is already initialized correctly
message("Looping through grna_assignments")
for (i in seq_along(individual_grna_assignments)) {
  indices <- individual_grna_assignments[[i]]
  if(length(indices) > 0) { # Check if there are any indices to assign
    sparseMatrix[i, indices] <- 1
  } else {
    # Optionally, handle the case where there are no indices
    # For example, by doing nothing, or logging this case
    cat(sprintf("No indices for row %d\n", i))
  }
}

# Set column names
colnames(sparseMatrix) <- colnames(sce)
rownames(sparseMatrix) <- names(individual_grna_assignments)

# Add to sce object
message("Adding grna_perts to sce object")
altExp(sce, e = "grna_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))










# Now add the individual "cre_perts"
message("Adding the individual cre_perts")
get_cre_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@grna_assignments$grna_group_idxs)
}

individual_cre_assignments <- get_cre_assignments(sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_cre_assignments)
nCols <- length(colnames(sce))

# Initialize an empty sparse matrix
sparseMatrix <- Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)

# Assuming sparseMatrix is already initialized correctly
message("Looping through cre_assignments")
for (i in seq_along(individual_cre_assignments)) {
  print(i)
  # Indices from individual_cre_assignments point to cells_in_use
  indices_to_cells_in_use <- individual_cre_assignments[[i]]

  # Now retrieve the actual indices from sceptre_object@cells_in_use
  indices <- sceptre_object@cells_in_use[indices_to_cells_in_use]

  if(length(indices) > 0) { # Check if there are any indices to assign
    sparseMatrix[i, indices] <- 1
  } else {
    # Optionally, handle the case where there are no indices
    # For example, by doing nothing, or logging this case
    cat(sprintf("No indices for row %d\n", i))
  }
}

# Set column names
colnames(sparseMatrix) <- colnames(sce)
rownames(sparseMatrix) <- names(individual_cre_assignments)

# Add to sce object
message("Adding cre_prets to sce object")
altExp(sce, e = "cre_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))





# Save to output file
message("Saving output to file")
saveRDS(sce, snakemake@output[[1]])




# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)




