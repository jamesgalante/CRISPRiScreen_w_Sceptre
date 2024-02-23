## Run Sceptre for given sample

save.image("sceptre_diff_expr.rda")
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
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(stringr)
  library(Matrix)
  library(data.table)
})

# Load Sceptre
devtools::install_github("katsevich-lab/sceptre")
library(sceptre)

# Load the arguments
mtx_name = snakemake@params$dge
perturb_status_name = snakemake@params$raw_perturb_status
gene_gRNA_group_pairs_name = snakemake@input$gene_gRNA_group_pairs
gRNA_groups_table_name = snakemake@input$gRNA_groups_table
sample_name = snakemake@params$sample_name
sample_number = snakemake@params$sample_number
guides_pre_assigned = snakemake@params$guides_pre_assigned
moi = snakemake@params$moi
side = snakemake@params$side
grna_integration_strategy = snakemake@params$grna_integration_strategy
guide_assignment_method = snakemake@params$guide_assignment_method
do_pos = snakemake@params$do_pos



# Create the output directory for other outputs
full_path <- snakemake@output$final_sceptre_object
base_directory <- dirname(full_path)
outdir <- file.path(base_directory, "outputs")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}


# Read in the gRNA tables
gene_gRNA_group_pairs = read.table(gene_gRNA_group_pairs_name,header=T,sep = '\t')
gRNA_groups_table_ej = read.table(gRNA_groups_table_name,header=T,sep = '\t')

# Read in the dge.txt.gz matrix for the gene counts
gene_matrix <- fread(mtx_name)
gene_matrix <- as(as.matrix(gene_matrix, rownames = 1), "sparseMatrix")

# Read in the raw_perturb_status.rds file
gRNA_matrix <- fread(snakemake@params$raw_perturb_status)
gRNA_matrix <- as(as.matrix(gRNA_matrix, rownames = 1), "sparseMatrix")

#need to add to server version
saveRDS(gRNA_matrix, snakemake@output$gRNA_matrix)
saveRDS(gene_matrix, snakemake@output$gene_matrix)


# I guess this is just always "no". I don't really know what it means - ask Evvie
if (guides_pre_assigned == "yes"){
  gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(perturbation_matrix),]
}
if (guides_pre_assigned == "no"){
  gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(gRNA_matrix),]
}

gene_gRNA_group_pairs1 = gene_gRNA_group_pairs[gene_gRNA_group_pairs$response_id %in% row.names(gene_matrix),]
response_names = row.names(gene_matrix)
colnames(gRNA_groups_table_ej_use)[2] = "grna_target"


# Format cell metadata for Sceptre input if cell_metadata exists
cell_metadata = snakemake@params$cell_metadata
if (!is.null(cell_metadata)) {
  extra_covariates <- data.frame(batch = as.factor(cell_metadata$cell_batches))
  rownames(extra_covariates) <- cell_metadata$cell_barcode
}


#takeout extra covariates if it's 1 sample
if (!is.null(cell_metadata)) {
  sceptre_object <- import_data(
    response_matrix = gene_matrix,
    grna_matrix = gRNA_matrix,
    grna_target_data_frame = gRNA_groups_table_ej_use,
    moi = moi,
    extra_covariates = extra_covariates,
    response_names = response_names
  )
} else {
  sceptre_object <- import_data(
    response_matrix = gene_matrix,
    grna_matrix = gRNA_matrix,
    grna_target_data_frame = gRNA_groups_table_ej_use,
    moi = moi,
    response_names = response_names
  )
}








colnames(gene_gRNA_group_pairs1)[2] = "grna_target"
gene_gRNA_group_pairs1$type = "x"

# Label the positive controls for Sceptre
for (i in 1:nrow(gene_gRNA_group_pairs1)){
  if (gene_gRNA_group_pairs1$target_type[i] == snakemake@params$tss_ctrl_label){
    gene_gRNA_group_pairs1$type[i] = "pos"
  }
}


positive_control_pairs =  gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type == "pos",c("grna_target", "response_id")]
discovery_pairs = gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type != "pos",c("grna_target", "response_id")]


sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = unique(discovery_pairs),
  positive_control_pairs = positive_control_pairs,
  side = side,
  grna_integration_strategy = grna_integration_strategy
)

print(sceptre_object)
#log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero + 1) + log(grna_n_umis + 1) + batch is the default? wanna put it back to what I had




# Helper function to save plots
save_plot <- function(plot, filename) {
  ggsave(filename = file.path(outdir, filename), plot = plot, width = 10, height = 8)
}

# Function to plot GRNA count distributions
plot_grna_count_distributions_fn <- function(sceptre_object) {
  p = plot_grna_count_distributions(sceptre_object, n_grnas_to_plot = 9)
  save_plot(p, "plot_grna_count_distributions.png")
}

# Function to assign GRNAs and plot cell assignment method
assign_grnas_and_plot <- function(sceptre_object, guide_assignment_method) {
  sceptre_object_assigned <- assign_grnas(sceptre_object, method = guide_assignment_method, parallel = FALSE)
  p = plot(sceptre_object_assigned, n_grnas_to_plot = 5)
  save_plot(p, paste0("plot_grna_cell_assignment_method_", guide_assignment_method, ".png"))
  return(sceptre_object_assigned)
}

# Function to plot covariates
plot_covariates_fn <- function(sceptre_object) {
  p = plot_covariates(sceptre_object)
  save_plot(p, "plot_covariates.png")
}

# Function to perform quality control and plot QC results
run_qc_and_plot <- function(sceptre_object) {
  sceptre_object_qc <- run_qc(sceptre_object)
  p = plot(sceptre_object_qc)
  save_plot(p, "plot_qc.png")
  return(sceptre_object_qc)
}

# Function to run calibration check and plot results
run_calibration_check_and_plot <- function(sceptre_object) {
  sceptre_object_calibrated <- run_calibration_check(sceptre_object, parallel = FALSE)
  p = plot(sceptre_object_calibrated)
  save_plot(p, "calibration_check.png")
  calibration_result <- get_result(sceptre_object_calibrated, analysis = "run_calibration_check")
  write.table(calibration_result, file.path(outdir, "calibration_result.txt"), quote = FALSE, sep = '\t', row.names = FALSE)
  return(sceptre_object_calibrated)
}

# Function to run power check and plot results
run_power_check_and_plot <- function(sceptre_object) {
  if (do_pos == "yes") {
    sceptre_object_power_checked <- run_power_check(sceptre_object, parallel = FALSE)
    p = plot(sceptre_object_power_checked)
    save_plot(p, "power_check.png")
    return(sceptre_object_power_checked)
  }
  return(sceptre_object) # Return unmodified object if power check is not performed
}

# Function to run discovery analysis and plot results
run_discovery_analysis_and_plot <- function(sceptre_object) {
  sceptre_object_discovery <- run_discovery_analysis(sceptre_object, parallel = FALSE)
  p = plot(sceptre_object_discovery)
  save_plot(p, "discovery_result.png")
  discovery_result <- get_result(sceptre_object_discovery, analysis = "run_discovery_analysis")
  write.table(discovery_result, snakemake@output$discovery_result, quote = FALSE, sep = '\t', row.names = FALSE)
  return(sceptre_object_discovery)
}

# Sequential execution
plot_grna_count_distributions_fn(sceptre_object)
sceptre_object <- assign_grnas_and_plot(sceptre_object, guide_assignment_method)
plot_covariates_fn(sceptre_object)
sceptre_object <- run_qc_and_plot(sceptre_object)
sceptre_object <- run_calibration_check_and_plot(sceptre_object)
sceptre_object <- run_power_check_and_plot(sceptre_object)
sceptre_object <- run_discovery_analysis_and_plot(sceptre_object)

# Write outputs to directory
write_outputs_to_directory(sceptre_object, outdir)

# Save the final Sceptre object
saveRDS(sceptre_object, paste0(snakemake@output$final_sceptre_object))





# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)

