## Rules to perform differential expression on a Perturb-seq experiment to e.g. map CREs to target 
## genes

# this workflow allows splitting of datasets by chromosome for easy parallelization. whether a
# sample should be split by chromosomes is specified in the config file via the 'split_by_chr'
# entry. rule ordering is used to apply chromosome splitting for samples that are set for splitting.
# wildcard constraints are used to define which samples should be split.

# set rule order so that rules splitting a sample by chromosome are preferred over other rules
ruleorder: combine_de_results > perform_de_tests

# get samples that should be split by chromsome and create wildcards constrain string
config["split_by_chr"].update({"DUMMY": True})  # add dummy sample TODO: find better solution ASAP!
split_samples = list(dict(filter(lambda x: x[1] == True, config["split_by_chr"].items())).keys())
split_samples_wildcards = "|".join(split_samples)

# Rules for differential expression tests ----------------------------------------------------------

# perform differential expression tests for all cis perturbation - gene pairs in the dataset
rule perform_de_tests:
  input: "results/{sample}/perturb_sce.rds"
  output: "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz"
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = config["diff_expr"]["cell_batches"],
    p_adj_method = config["diff_expr"]["padj_method"],
    seed = 20210928
  log: "results/{sample}/logs/diff_expr/diff_expr_{method}_{strategy}.log"
  threads: config["diff_expr"]["threads"]
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "24G",
    time = "6:00:00"
  script:
    "../scripts/differential_expression.R"

# Rules to split dataset by chromosome for DE tests ------------------------------------------------

# extract data for one chromosome
rule extract_chromosome:
  input: "results/{sample}/perturb_sce.rds"
  output: temp("resources/{sample}/perturb_sce.{chr}.rds")
  params:
    rm_zero_cells = True
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/extract_chrom_from_sce.R"

# perform differential expression tests for one chromosome
rule perform_de_tests_chr:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output: temp("results/{sample}/diff_expr/{chr}_output_{method}_{strategy}.tsv.gz")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = config["diff_expr"]["cell_batches"],
    p_adj_method = config["diff_expr"]["padj_method"],
    seed = 20210928
  log: "results/{sample}/logs/diff_expr/diff_expr_{chr}_{method}_{strategy}.log"
  threads: config["diff_expr"]["threads"]
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "24G",
    time = "6:00:00"
  script:
    "../scripts/differential_expression.R"

# run DE tests for all chromosomes and combine into one file
rule combine_de_results:
  input:
    expand("results/{{sample}}/diff_expr/{chr}_output_{{method}}_{{strategy}}.tsv.gz",
      chr = ["chr" + str(i) for i in [*range(1, 23), "X"]])
  output: "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz"
  wildcard_constraints:
    sample = split_samples_wildcards
  params:
    padj_method = config["diff_expr"]["padj_method"]
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/combine_diff_expr_tests.R"
