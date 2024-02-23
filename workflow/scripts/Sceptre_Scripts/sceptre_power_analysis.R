

# save.image("sceptre_power_anal.rda")
# message("Saved Image")
# stop()

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(stringr)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
})


# Load Sceptre
message("Loading Sceptre")
devtools::install_github("katsevich-lab/sceptre")
library(sceptre)



# Read the input arguments into variables
message("Reading in snakemake variables")
guide_file_name <- snakemake@input$guide_file_names
sceptre_object_name <- snakemake@input$sceptre_object_name
full_grna_matrix_name	<- snakemake@input$full_grna_matrix_name
full_response_matrix_name <- snakemake@input$full_response_matrix_name
sce_object_name <- snakemake@input$sce_object_name
n_pert <- snakemake@params$n_pert
effect <- snakemake@params$effect_size
reps <- snakemake@params$reps
center <- snakemake@params$center



# Load in RDS files
message("Loading RDS files")
sceptre_object = readRDS(sceptre_object_name)
full_grna_matrix = readRDS(full_grna_matrix_name)
full_response_matrix = readRDS(full_response_matrix_name)



# Create the output directory for temporary files
message("Creating output directory")
full_path <- snakemake@output[[1]]
base_directory <- dirname(full_path)
outdir <- file.path(base_directory, "outputs")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Add the output name for temporary files
outname = paste0("sceptre_power_sim_effect_size_",effect,"_guide_subset_",basename(guide_file_name), "_n_perts_",n_pert,"_rep_",reps,"_center_effect_size_",center,".txt")


# prepare data =====================================================================================

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
sce <- readRDS(sce_object_name)

# infer perturbation level based on strategy
pert_level <- switch("perCRE", "perGRNA" = "grna_perts", "perCRE" = "cre_perts",
                     stop("incorrect strategy argument"))

# filter for minimum number of cells per perturbation

# perform power simulations ========================================================================
message("Not normalizing transcript counts because this is sceptre.")


# convert 'percentage decrease' effect size to 'relative expression level'
effect_size <- 1 - as.numeric(effect)

# simulate Perturb-seq data and perform differential gene expression tests
message("Performing power simulations.")

# Set variables
message("Hard Code variables")
genes_iter = FALSE
guide_sd = as.numeric(0.13)
norm = "real"


col_names <- colnames(colData(sce)) 
if ("pert" %in% col_names) stop("'pert' cannot be a colData name, please rename.", call. = FALSE)

# check that pert_level is valid
if (!pert_level %in% altExpNames(sce)) {
  stop("pert_level must be one of the altExp names: ", paste(altExpNames(sce), collapse = ", "),
       call. = FALSE)
}

# check that guide_sd is valid
if (!is.numeric(guide_sd) | guide_sd < 0) {
  stop("Invalid 'guide_sd' value. Must be numeric >= 0.", call. = FALSE)
}

# get required functions -------------------------------------------------------------------------
#trying without sampling so that we can just use straight


n_ctrl = FALSE
# get function to generate input data for one perturbation
if (is.numeric(n_ctrl)) {
  pert_input_function <- pert_input_sampled
  n_ctrl <- as.integer(n_ctrl)
} else if (n_ctrl == FALSE) {
  pert_input_function <- pert_input
} else {
  stop("Invalid 'n_ctrl' argument.", call. = FALSE)
}


message("Following sceptre, not sampling control cells")

if (norm == "real") {
  sim_function <- simulate_diff_expr_pert_real
} else if (norm == "sim_nonpert") {
  sim_function <- simulate_diff_expr_pert_sim
} else {
  stop("Invalid 'norm' argument.", call. = FALSE)
}



# Function to get the simulate count matrix
sim_counts_submit <- function(sce, effect_size_mat) {
  
  gene_means = rowData(sce)[, "mean"]
  cell_size_factors = colData(sce)[, "size_factors"]
  gene_ids = names(sce)
  cell_ids = names(cell_size_factors)
  # simulate Perturb-seq count data with parameters from SCE object
  sim_counts <- simulate_tapseq_counts(gene_means = rowData(sce)[, "mean"],
                                       gene_dispersions = rowData(sce)[, "dispersion"],
                                       cell_size_factors = colData(sce)[, "size_factors"],
                                       effect_size_mat = effect_size_mat, gene_ids = gene_ids,cell_ids = cell_ids)
  row.names(sim_counts) = gene_ids
  
  message("Returning sim_counts")
  return(sim_counts)
  
}

# perform simulated DE tests ---------------------------------------------------------------------

#TO DO NOTE: for sceptre, want to make these control cells and stuff exactly the same - come back to that
#get the discovery_pairs to test
discovery_relevant_pairs_all = sceptre_object@discovery_pairs_with_info

#guide_file = read.table(guide_file_name,header = F,sep = '\t')
guide_file = read.table(guide_file_name,header = F,sep = '\t')

#discovery_relevant_pairs_in_guide_subset = discovery_relevant_pairs_all[discovery_relevant_pairs_all$grna_group %in% guide_file$V1,]
discovery_relevant_pairs_in_guide_subset = discovery_relevant_pairs_all[discovery_relevant_pairs_all$grna_group %in% guide_file[[2]],]

guide_targets_sce <- create_guide_targets(sce, pert_level = pert_level)

discovery_relevant_pairs = discovery_relevant_pairs_in_guide_subset[discovery_relevant_pairs_in_guide_subset$response_id %in% rownames(counts(sce)) & discovery_relevant_pairs_in_guide_subset$grna_group %in% rownames(altExp(sce,pert_level)),]

perts = unique(discovery_relevant_pairs$grna_group)


guide_targets = sceptre_object@grna_target_data_frame
colnames(guide_targets)[2] = "target_id"

discovery_results = data.frame()

if (n_pert == "all"){
  n_pert = length(perts)
}

# Creates temp file
# cat(paste0(c("response_id","grna_target","n_nonzero_trt","n_nonzero_cntrl","pass_qc","p_value","log_2_fold_change","significant","rep\n"),collapse = "\t"),
#     file = paste0(outdir,"/tempfile_",outname))


### Basically, I have to get rid of the dispersion estimates that are NA. I added a line that I thought would help
### in power_simulations_fun.R in the simulate_diff_expr_pert_real() function. BUt I still got the error I was dealing with
### Also make sure to fix the for loops etc.

# Runs each power analysis rep for one perturbation
for (pert in perts){
    
  discovery_relevant_pairs_pert = discovery_relevant_pairs[discovery_relevant_pairs$grna_group == pert & discovery_relevant_pairs$pass_qc == "TRUE",]
  
  # Make sure there are pairs in the table
  if (!(nrow(discovery_relevant_pairs_pert) > 0)){
    message(paste0("Skipping ", pert, ". No discovery relevant pairs pert."))
    next
  }
  
  # Generate the pert object
  pert_guides <- guide_targets[guide_targets$target_id == pert, "grna_id"]
  pert_object <-  pert_input_function(pert, sce = sce, pert_level = pert_level)
  
  # get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  
  
  ### Getting a weird error here where a `pert_guide` is not in the `grna_perts` matrix.
  ### I have no idea why this would happen to be honest, but I think whatever happens, it's happening before this script
  ### The sce object that gets imported in this pipeline for instance, also doesn't have the specific `pert_guide`
  ### But it exists in the original perturb_status.txt.gz
  ### Perhaps something occurs when I recombine all the dispersion estimates
  ### Perhaps also have to figure out where grna perts are getting filtered and see what's going on there
  ### For now, I'm just going to add a line to remove all pert_guides that aren't in `grna_perts`
  pert_guides <- pert_guides[pert_guides %in% rownames(grna_perts)]
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts,
                                               pert_guides = pert_guides)
  
  pert_genes = discovery_relevant_pairs_pert$response_id
  
  pert_object = pert_object[pert_genes,]
  
  effect_sizes <- structure(rep(1, nrow(pert_object)), names = rownames(pert_object))
  
  effect_sizes[pert_genes] <- effect_size #does this only for the genes you care about will need to adjust to match sceptre

  
  # Make sure there are non NA dispersion estimates in the pert_object
  if (all(is.na(rowData(pert_object)$dispersion))) {
    message(paste0("Skipping ", pert, ". No valid dispersion values."))
    next
  }
  
  # Remove any rows with NA dispersion values
  filtered_pert_object <- pert_object[!is.na(rowData(pert_object)$dispersion), ]
  
  
  
  for (i in 1:as.numeric(reps)) {

    message("creating simulation for rep,",i)
    
    # Create effect size matrix (sampled from negative binomial distribution around es or 1)
    es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                        gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
    
    # center effect sizes on specified gene-level effect sizes
    if (center == TRUE) {
      es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status,
                                          gene_effect_sizes = effect_sizes)
    }
    
    # Reorder columns
    es_mat_use = es_mat[,colnames(counts(pert_object))]
    
    
    #so now we have the simulated object for the pertrubation with every gene it should be tested with (might need to change with sceptre flow))
    message("Simulating Counts")
    sim_counts = sim_counts_submit(pert_object, effect_size_mat = es_mat_use)
    
    
    # Edge case if there's only one gene for perturbation to keep structure same
    # message("Processing Edge Cases")
    # if (length(pert_genes) > 1) {
    #   full_response_matrix_sim = full_response_matrix[pert_genes,]
    #   shared_cells = intersect(colnames(sim_counts),colnames(full_response_matrix))
    #   full_response_matrix_sim[rownames(sim_counts),shared_cells] = sim_counts[,shared_cells]
    # } else {
    #   sim_counts_mx = as.matrix(sim_counts)
    #   #need to add another row to keep the structure
    #   if (pert_genes[1] != rownames(full_response_matrix)[1]){
    #     alt_gene = rownames(full_response_matrix)[1]
    #   } else {
    #     alt_gene = rownames(full_response_matrix)[2]
    #   }
    #   full_response_matrix_sim = full_response_matrix[c(pert_genes,alt_gene),]
    #   shared_cells = intersect(colnames(sim_counts_mx),colnames(full_response_matrix))
    #   full_response_matrix_sim[pert_genes,shared_cells] = sim_counts_mx[pert_genes,shared_cells]
    # }
    
    
    
    message("Setting up sceptre object for Disovery Analysis")
    #full_response_matrix_sim_sparse = as(as.matrix(full_response_matrix_sim),"RsparseMatrix")
    full_response_matrix_sim_sparse = as(as.matrix(sim_counts),"RsparseMatrix")
    sceptre_object_use = sceptre_object
    sceptre_object_use@response_matrix =  full_response_matrix_sim_sparse
    sceptre_object_use@discovery_pairs_with_info = discovery_relevant_pairs_pert
    
    message("Running discovery analysis")
    sceptre_object2 <- run_discovery_analysis(
      sceptre_object = sceptre_object_use,
      parallel = FALSE
    )
    
    message("Returning discovery results")
    discovery_result <- get_result(
      sceptre_object = sceptre_object2,
      analysis = "run_discovery_analysis"
    )
    
    discovery_result$rep = i
    #write.table(discovery_result, col.names=F, append = T, row.names = F,quote = F,sep = '\t',file = paste0(outdir,"/tempfile_",outname))
    discovery_results = data.frame(rbind(discovery_results,discovery_result))
  }
}

message("Processing output.")
disp_outlier <- data.frame(gene = rownames(rowData(sce)),
                           disp_outlier_deseq2 = rowData(sce)[, "disp_outlier_deseq2"],
                           stringsAsFactors = FALSE)

av_expr = data.frame(rowMeans(counts(sce)))
colnames(av_expr)[1] = "average_expression_all_cells"
av_expr$response_id = rownames(av_expr)

colnames(disp_outlier)[1] = "response_id"


# add to output
discovery_results <- left_join(discovery_results, disp_outlier, by = "response_id")
discovery_results <- left_join(discovery_results, av_expr, by = "response_id")

# save simulation output
message("Saving output to file.")
write_tsv(discovery_results, file = snakemake@output[[1]])


# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
