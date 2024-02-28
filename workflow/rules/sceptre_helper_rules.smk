
# Set a rule order, so that the recombined negative binomial distribution is created before the power analysis is done
ruleorder: combine_fit_negbinom_distr_chr > sceptre_power_analysis

# function to get metadata
def get_cell_metadata(wildcards):
  cell_metadata = config["samples"][wildcards.sample]["cell_metadata"]
  if cell_metadata is None:
    cell_metadata = []
  return(cell_metadata)


# Combine each chromosomes neg binom fit into one SingleCellExperiment for Sceptre Input
rule combine_fit_negbinom_distr_chr:
  input:
    perturb_sce = "results/{sample}/perturb_sce.rds",
    chr_rds_objects = expand("resources/{{sample}}/perturb_sce_disp.{chr}.rds", chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]])
  output: "resources/{sample}/recombined_perturb_sce_disp.rds"
  params:
    size_factors = config["power_simulations"]["size_factors"]
  log: "results/{sample}/logs/power_sim/combine_fit_negbinom_distr_chr.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "88G",
    time = "5:00:00"
  script:
    "../scripts/Sceptre_Scripts/combine_fit_negbinom_distr.R"


# These rules are used to create inputs to the Sceptre Diffex and Power Analysis
rule create_sceptre_input_files:
  input:
    power_analysis_output = lambda wildcards: expand("results/" + wildcards.sample + "/output_{sd}gStd_{method}_{strategy}.tsv.gz", sd=["0.13"], method=["MAST"], strategy=["perCRE"], wildcards=wildcards),
    guide_targets_w_NT = lambda wildcards: config["samples"][wildcards.sample]["guide_targets_w_NT"]
  output:
    "resources/{sample}/Sceptre/gene_gRNA_group_pairs.txt", "resources/{sample}/Sceptre/gRNA_groups_table.txt"
  params:
    tss_ctrl_label = lambda wildcards: config["samples"][wildcards.sample]["tss_ctrl_label"]
  log: "results/{sample}/logs/create_sceptre_input_files.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/Sceptre_Scripts/create_sceptre_input_files.R"


# Perform the discovery analysis
rule sceptre_differential_expression:
  input:
    gene_gRNA_group_pairs = "resources/{sample}/Sceptre/gene_gRNA_group_pairs.txt",
    gRNA_groups_table = "resources/{sample}/Sceptre/gRNA_groups_table.txt"
  output:
    gene_matrix = "results/{sample}/sceptre_diff_expr/gene_matrix.rds",
    gRNA_matrix = "results/{sample}/sceptre_diff_expr/gRNA_matrix.rds",
    final_sceptre_object = "results/{sample}/sceptre_diff_expr/final_sceptre_object.rds",
    discovery_result = "results/{sample}/sceptre_diff_expr/discovery_result.txt"
  params:
    dge = lambda wildcards: config["samples"][wildcards.sample]["dge"],
    raw_perturb_status = lambda wildcards: config["samples"][wildcards.sample]["raw_pert"],
    cell_metadata = get_cell_metadata,
    sample_number = config["sceptre_diff_expr"]["sample_number"],
    guides_pre_assigned = config["sceptre_diff_expr"]["guides_pre_assigned"],
    moi = config["sceptre_diff_expr"]["moi"],
    side = config["sceptre_diff_expr"]["side"],
    grna_integration_strategy = config["sceptre_diff_expr"]["grna_integration_strategy"],
    guide_assignment_method = config["sceptre_diff_expr"]["guide_assignment_method"],
    do_pos = config["sceptre_diff_expr"]["do_pos"],
    tss_ctrl_label = lambda wildcards: config["samples"][wildcards.sample]["tss_ctrl_label"]
  log: "results/{sample}/logs/sceptre_differential_expression.log"
  conda: "../envs/R_env.yml"
  resources:
    mem = "88G",
    time = "5:00:00"
  script:
    "../scripts/Sceptre_Scripts/new_sceptre_for_server.R"


# Split the guide file into "N" files for easy parallelization
rule split_guide_file:
    input:
        gene_gRNA_group_pairs = "resources/{sample}/Sceptre/gene_gRNA_group_pairs.txt"
    output:
        temp(expand("resources/{{sample}}/Sceptre/gene_gRNA_group_pairs_split_{split}.txt", split=range(1, config['sceptre_pwr_anal']['batches'] + 1)))
    params:
        batches = config['sceptre_pwr_anal']['batches']
    shell:
        """
        tail -n +2 {input.gene_gRNA_group_pairs} > {input.gene_gRNA_group_pairs}.tmp
        
        total_lines=$(wc -l < {input.gene_gRNA_group_pairs}.tmp)
        batches={params.batches}
        lines_per_batch=$((total_lines / batches))
        extra_lines=$((total_lines % batches))
        if [ $extra_lines -gt 0 ]; then
            lines_per_batch=$(($lines_per_batch + 1))
        fi
        split -d -l $lines_per_batch --numeric-suffixes=1 {input.gene_gRNA_group_pairs}.tmp resources/{wildcards.sample}/Sceptre/gene_gRNA_group_pairs_split_
        for file in resources/{wildcards.sample}/Sceptre/gene_gRNA_group_pairs_split_*; do
            # Remove leading zero for files numbered less than 10
            new_name=$(echo $file | sed 's/split_0/split_/')
            mv "$file" "${{new_name}}.txt"
        done
        
        rm {input.gene_gRNA_group_pairs}.tmp
        """   
    
    
    
    # shell:
    #     """
    #     total_lines=$(wc -l < {input.gene_gRNA_group_pairs})
    #     batches={params.batches}
    #     lines_per_batch=$((total_lines / batches))
    #     extra_lines=$((total_lines % batches))
    #     if [ $extra_lines -gt 0 ]; then
    #         lines_per_batch=$(($lines_per_batch + 1))
    #     fi
    #     split -d -l $lines_per_batch --numeric-suffixes=1 {input.gene_gRNA_group_pairs} resources/{wildcards.sample}/Sceptre/gene_gRNA_group_pairs_split_
    #     for file in resources/{wildcards.sample}/Sceptre/gene_gRNA_group_pairs_split_*; do
    #         # Remove leading zero for files numbered less than 10
    #         new_name=$(echo $file | sed 's/split_0/split_/')
    #         mv "$file" "${{new_name}}.txt"
    #     done
    #     """



# This currently only processes CRISPRi runs that have been separated by chromosome for parallel processing
# This is due to the `recombined_perturb_sce_disp.rds` instead of just `perturb_sce_disp.rds` which is what's created normally
# Add a mechanism, so that this can run with both scenarios
rule sceptre_power_analysis:
  input:
    # Maybe change these names to match the output names of the rules that lead up to it
    sceptre_object_name	= "results/{sample}/sceptre_diff_expr/final_sceptre_object.rds",
    full_grna_matrix_name	= "results/{sample}/sceptre_diff_expr/gRNA_matrix.rds",
    full_response_matrix_name =	"results/{sample}/sceptre_diff_expr/gene_matrix.rds",
    sce_object_name	= "resources/{sample}/recombined_perturb_sce_disp.rds",
    guide_file_names = "resources/{sample}/Sceptre/gene_gRNA_group_pairs_split_{split}.txt"
  output:
    temp("results/{sample}/sceptre_pwr_analysis/power_analysis_output_{split}.tsv")
  params:
    n_pert = config["sceptre_pwr_anal"]["n_pert"],
    effect_size = config["sceptre_pwr_anal"]["effect_size"],
    reps = config["sceptre_pwr_anal"]["reps"],
    center = config["sceptre_pwr_anal"]["center"],
    n_ctrl = config["sceptre_pwr_anal"]["n_ctrl"],
    cell_batches = config["sceptre_pwr_anal"]["cell_batches"]
  log: "results/{sample}/logs/sceptre_power_analysis_{split}.log"
  conda:
    "../envs/sceptre_pwr_env.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/Sceptre_Scripts/sceptre_power_analysis.R"


# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis:
  input:
    expand("results/{{sample}}/sceptre_pwr_analysis/power_analysis_output_{split}.tsv", split=range(1, config['sceptre_pwr_anal']['batches'] + 1))
  output:
    "results/{sample}/sceptre_pwr_analysis/power_analysis_output.tsv"
  log: "results/{sample}/logs/combine_sceptre_power_analysis.log"
  conda:
    "../envs/sceptre_pwr_env.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/Sceptre_Scripts/combine_sceptre_power_analysis.R"








