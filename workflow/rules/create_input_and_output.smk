## Rules to create worflow input and output

# function to get samples that require liftover from hg19 to GRCh38
def get_cell_metadata(wildcards):
  cell_metadata = config["samples"][wildcards.sample]["cell_metadata"]
  if cell_metadata is None:
    cell_metadata = []
  return(cell_metadata)

# download gencode annotations
rule download_gencode_annotations:
  output: "resources/{annot}.annotation.gtf.gz"
  params:
    url = lambda wildcards: config["download_urls"][wildcards.annot]
  conda: "../envs/analyze_crispr_screen.yml"
  shell:
    "wget -O {output} {params.url}"
    
# create a SingleCellExperiment object from following input data:
# - gene expression matrix
# - perturbation matrix
# - guide targets file
# - genome annotations
rule create_sce:
  input:
    dge = lambda wildcards: config["samples"][wildcards.sample]["dge"],
    pert_status = lambda wildcards: config["samples"][wildcards.sample]["pert"],
    guide_targets = lambda wildcards: config["samples"][wildcards.sample]["guide_targets"],
    annot = lambda wildcards: config["samples"][wildcards.sample]["annot"],
    cell_metadata = get_cell_metadata
  output: "results/{sample}/perturb_sce.rds"
  log: "results/{sample}/logs/create_sce.log"
  params:
    vector_pattern = lambda wildcards: config["samples"][wildcards.sample]["dge_vector_pattern"]
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G"  # Changed from mem_mb to mem
  script:
    "../scripts/create_sce_object.R"
    
# create output data by combining differential expression and power simulation results
rule create_main_output:
  input:
    diff_expr = "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz",
    power_sim = expand("results/{{sample}}/power_sim/power_{effect}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
                  effect = config["power_simulations"]["effect_sizes"])
  output: "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    size_factors = config["power_simulations"]["size_factors"]
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/create_output.R"
    
# plot differential expression results
rule diff_expr_results:
  input:
    "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz"
  output:
    "results/{sample}/diff_expr_{method}_{strategy}.html"
  params:
    fdr_threshold = config["diff_expr"]["padj_threshold"],
    tss_ctrl_label = lambda wildcards: config["samples"][wildcards.sample]["tss_ctrl_label"]
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/diff_expr_results.Rmd"
    
# plot power analysis results
rule power_analysis:
  input:
    "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz"
  output:
    "results/{sample}/power_analysis_{sd}gStd_{method}_{strategy}.html"
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/power_analysis.Rmd"
   
# format sceptre output for compatibility with ENCODE pipelines
rule format_sceptre_output:
  input:
    power_analysis_output = "results/{sample}/sceptre_pwr_analysis/power_analysis_output.tsv",
    discovery_results = "results/{sample}/sceptre_diff_expr/discovery_result.txt",
    gene_gRNA_group_pairs = "resources/{sample}/Sceptre/gene_gRNA_group_pairs.txt"
  output:
    "results/{sample}/final_sceptre_output/output_0.13gStd_MAST_perCRE.tsv"
  params:
    effect_size = config["sceptre_pwr_anal"]["effect_size"],
    annot = lambda wildcards: config["samples"][wildcards.sample]["annot"],
    gene_format = config["sceptre_pwr_anal"]["gene_format"]
  log: "results/{sample}/logs/format_sceptre_output.log"
  conda: "../envs/sceptre_pwr_env.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/Sceptre_Scripts/format_sceptre_output.R"
