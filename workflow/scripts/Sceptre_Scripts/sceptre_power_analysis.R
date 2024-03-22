
### Main issues with this script that need to be addressed
### At some point, doubling of response_id, grna_target pairs is occurring - I deal with this in the combine_pwr_analysis script
### But it would be useful to know why that's occurring in the first place
### Also, I had to add this line: pert_guides <- pert_guides[pert_guides %in% rownames(grna_perts)]
### Because some pert_guides were not in the rownames of grna_perts - I have a more detailed description of this below

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
  source(file.path(snakemake@scriptdir, "../R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "../R_functions/power_simulations_fun.R"))
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
n_ctrl <- snakemake@params$n_ctrl
cell_batches <- snakemake@params$cell_batches


# convert 'percentage decrease' effect size to 'relative expression level'
effect_size <- 1 - as.numeric(effect)


# Hard Code Variables
message("Hard Coding variables: guide_sd = as.numeric(0.13)")
guide_sd <- as.numeric(0.13)

# infer perturbation level based on strategy
message("Hard Coding variables: perCRE instead of perGRNA")
pert_level <- switch("perCRE", "perGRNA" = "grna_perts", "perCRE" = "cre_perts",
                     stop("incorrect strategy argument"))





# Load in RDS files
message("Loading RDS files")
sceptre_object <- readRDS(sceptre_object_name)
full_grna_matrix <- readRDS(full_grna_matrix_name)
full_response_matrix <- readRDS(full_response_matrix_name)
sce <- readRDS(sce_object_name)

# Check format of sce object
col_names <- colnames(colData(sce)) 
if ("pert" %in% col_names) stop("'pert' cannot be a colData name, please rename.", call. = FALSE)


# Loading guide file
message("Loading guide file")
guide_file <- read.table(guide_file_name, header = F, sep = '\t')
# Set the column names for the split guide_file
colnames(guide_file) <- c("response_id", "grna_group", "target_type")







# Assigning appropriate sampling function based on status of n_ctrl
message("Assigning pert_input_function with n_ctrl value")
if (is.numeric(n_ctrl)) {
  pert_input_function <- pert_input_sampled
  n_ctrl <- as.integer(n_ctrl)
} else if (n_ctrl == FALSE) {
  pert_input_function <- pert_input
} else {
  stop("Invalid 'n_ctrl' argument.", call. = FALSE)
}


# Function to get the simulate count matrix
sim_counts_submit <- function(sce, effect_size_mat) {
  
  gene_means <- rowData(sce)[, "mean"]
  cell_size_factors <- colData(sce)[, "size_factors"]
  gene_ids <- names(sce)
  cell_ids <- names(cell_size_factors)
  # simulate Perturb-seq count data with parameters from SCE object
  sim_counts <- simulate_tapseq_counts(gene_means = rowData(sce)[, "mean"],
                                       gene_dispersions = rowData(sce)[, "dispersion"],
                                       cell_size_factors = colData(sce)[, "size_factors"],
                                       effect_size_mat = effect_size_mat, gene_ids = gene_ids,cell_ids = cell_ids)
  row.names(sim_counts) <- gene_ids
  
  message("Returning sim_counts")
  return(sim_counts)
  
}




# Get all the perturbation targets in the subsetted guide_file that are represented in the sceptre_object and where at least one pair passes qc
discovery_pairs_which_pass_qc <- sceptre_object@discovery_pairs_with_info[sceptre_object@discovery_pairs_with_info$pass_qc == TRUE,]
perts <- intersect(guide_file$grna_group, discovery_pairs_which_pass_qc$grna_group)



# Set up an empty datafame to store the results
discovery_results <- data.frame()
# Runs each power analysis rep for one perturbation
for (pert in perts){
  
  repeat_loop <- TRUE
  
  # There's a slight chance that the simulated counts for a perturbation cause glm.fit to not converge
  # We include a while loop to redo the for loop pert simulation when this occurs because the p_value > 2
  # Let's make sure we don't run into an infinite error
  number_of_repeats <- 0
  
  while (repeat_loop) {
    
    if (number_of_repeats > 10) {
      message(paste0("Ran into an error with ", pert, ". As it consistently simulated counts which the glm couldn't fit"))
      break
    }
    
    # Get a list of the relevant discovery pairs
    discovery_relevant_pairs_pert <- discovery_pairs_which_pass_qc[discovery_pairs_which_pass_qc$grna_group == pert,]
    
    # Get all the guides that target the current `pert`
    pert_guides <- sceptre_object@grna_target_data_frame[sceptre_object@grna_target_data_frame$grna_target == pert, "grna_id"]
    
    # Get the pert input
    pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level, n_ctrl = n_ctrl, cell_batches = cell_batches)
    
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
      
      
      # So now we have the simulated object for the pertrubation with every gene it should be tested with (might need to change with sceptre flow))
      message("Simulating Counts")
      sim_counts = sim_counts_submit(pert_object, effect_size_mat = es_mat_use)
      
      # Now let's run the discovery analysis
      message("Setting up sceptre object for Disovery Analysis")
      full_response_matrix_sim_sparse <- as(as.matrix(sim_counts),"RsparseMatrix")
      sceptre_object_use <- sceptre_object
      sceptre_object_use@response_matrix <- full_response_matrix_sim_sparse
      sceptre_object_use@discovery_pairs_with_info <- discovery_relevant_pairs_pert
      
      # Fix the `cells_in_use` parameter for indexing when the n_ctrl is not FALSE
      if (is.numeric(n_ctrl)) {
        sceptre_object_use@cells_in_use <- seq(dim(full_response_matrix_sim_sparse)[[2]])
        
        # Only keep the covariate rows that are in the colnames of `full_response_matrix_sim_sparse`
        sceptre_object_use@covariate_data_frame <- sceptre_object_use@covariate_data_frame[rownames(sceptre_object_use@covariate_data_frame) %in% colnames(full_response_matrix_sim_sparse),]
        sceptre_object_use@covariate_matrix <- sceptre_object_use@covariate_matrix[rownames(sceptre_object_use@covariate_matrix) %in% colnames(full_response_matrix_sim_sparse),]
        
        sceptre_object_use@grna_assignments$grna_group_idxs[pert] <- list(seq(sum(colData(pert_object)$pert == 1)))
        #rownames(colData(pert_object))[colData(pert_object)$pert == 1]
      }
      
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
      
      # Check if any p_value is greater than 1
      if (any(discovery_result$p_value > 1)) {
        message(paste0("Repeating analysis for ", pert, " due to p_value > 1."))
        number_of_repeats <- number_of_repeats + 1
        # The loop will repeat
      } else {
        # The loop will not repeat
        repeat_loop <- FALSE
        
        # If everything is good, add to the final results
        discovery_result$rep <- i
        discovery_results <- data.frame(rbind(discovery_results, discovery_result))
      }
    }
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