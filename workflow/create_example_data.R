
# required packages
library(data.table)
library(tidyverse)
library(here)
library(GenomicRanges)

# load data ----------------------------------------------------------------------------------------

# perturbation status matrix
pert <- fread(here("../ENCODE_CRISPR_data/resources/TAPseqChr11/perturb_status.txt"))

# guide targets
guide_targets <- fread(here("../ENCODE_CRISPR_data/resources/TAPseqChr11/guide_targets_hg38.tsv"))

# dge data
dge <- fread(here("../ENCODE_CRISPR_data/resources/TAPseqChr11/dge.txt"))

# process perturbation status ----------------------------------------------------------------------

set.seed(42)

# get coordinates of all screen guides (no positive or negative controls)
guide_coords <- pert %>% 
  select(guide_id = VECTOR) %>% 
  filter(grepl(guide_id, pattern = "^chr")) %>% 
  separate(guide_id, into = c("chr", "start", "end"), sep = ":|-|_", extra = "drop",
           remove = FALSE) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# get guides overlapping example regions
example_region <- GRanges(seqnames = "chr11", ranges = IRanges(273200, 358200), strand = "*")
region_guides <- subsetByOverlaps(guide_coords, example_region)

# get perturbation status for these guides only
region_pert <- pert[pert$VECTOR %in% region_guides$guide_id, ]

# get perturbed cells and a random sample of non-perturbed cells
pert_cells <- colnames(region_pert[, -1])[colSums(region_pert[, -1]) > 0]
non_pert_cells <- colnames(region_pert[, -1])[colSums(region_pert[, -1]) == 0]
sampled_non_pert_cells <- sample(non_pert_cells, size = 10000 - length(pert_cells))

# get pert status for these cells and save to file
cols_output <- c("VECTOR", pert_cells, sampled_non_pert_cells)
pert_output <- region_pert[, colnames(region_pert) %in% cols_output, with = FALSE]

# process guide targets file -----------------------------------------------------------------------

# get guide targets for guides in example region
guide_targets_output <- filter(guide_targets, name %in% pert_output$VECTOR)

# process dge data ---------------------------------------------------------------------------------

# remove any crop-seq vector transcripts to keep gene expression data only
dge_genes <- filter(dge, !grepl(GENE, pattern = "^CROPseq_dCas9_DS_.+$"))

# only retain data on cells from example perturbation status
cols_output <- c("GENE", pert_cells, sampled_non_pert_cells)
dge_output <- dge_genes[, colnames(dge_genes) %in% cols_output, with = FALSE]

# save to output files -----------------------------------------------------------------------------

# save perturbation status and dge data to compressed tsv files
write_tsv(pert_output, file = here("resources/TAPseq_example/perturb_status.txt.gz"))
write_tsv(dge_output, file = here("resources/TAPseq_example/dge.txt.gz"))

# save example guide targets file to tsv
write_tsv(guide_targets_output, file = here("resources/TAPseq_example/guide_targets.tsv"))
