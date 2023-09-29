## Rules to perform empirical power simulations for Perturb-seq experiment

# this workflow allows splitting of datasets by chromosome for easy parallelization. whether a
# sample should be split by chromosomes is specified in the config file via the 'split_by_chr'
# entry. rule ordering is used to apply chromosome splitting for samples that are set for splitting.
# wildcard constraints are used to define which samples should be split.

# set rule order so that rules splitting a sample by chromosome are preferred over other rules
ruleorder: compute_power_chrs > compute_power

# get samples that should be split by chromsome and create wildcards constrain string
config["split_by_chr"].update({"DUMMY": True})  # add dummy sample TODO: find better solution ASAP!
split_samples = list(dict(filter(lambda x: x[1] == True, config["split_by_chr"].items())).keys())
split_samples_wildcards = "|".join(split_samples)

# Rules for power simulations ----------------------------------------------------------------------

# fit negative binomial distribution to estimate dispersions
rule fit_negbinom_distr:
  input: "results/{sample}/perturb_sce.rds"
  output: temp("resources/{sample}/perturb_sce_disp.rds")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"],
    remove_genes = config["power_simulations"]["remove_genes"],
    size_factors = config["power_simulations"]["size_factors"],
    fit_type = config["power_simulations"]["fit_type"]
  log: "results/{sample}/logs/power_sim/fit_negbinom_distr.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/fit_negbinom_distr.R"

# perform power simulations by submitting one iteration at a time
rule perform_power_simulations:
  input: "resources/{sample}/perturb_sce_disp.rds"
  output:
    temp("results/{sample}/power_sim/rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    norm = config["power_simulations"]["norm"],
    cell_batches = config["diff_expr"]["cell_batches"],
    genes_iter = False
  log: "results/{sample}/logs/power_sim/power_sim_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "24G",
    time = "24:00:00"
  script:
   "../scripts/power_simulations.R"

# compute power
rule compute_power:
  input:
    expand("results/{{sample}}/power_sim/rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/{sample}/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    p_adj_method = config["diff_expr"]["padj_method"],
    pval_threshold = config["diff_expr"]["padj_threshold"]
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/compute_power.R"

# Rules for simulations split by chromosome --------------------------------------------------------

# fit negative binomial distribution to estimate dispersions
rule fit_negbinom_distr_chr:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output: temp("resources/{sample}/perturb_sce_disp.{chr}.rds")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"],
    remove_genes = config["power_simulations"]["remove_genes"],
    size_factors = config["power_simulations"]["size_factors"],
    fit_type = config["power_simulations"]["fit_type"]
  log: "results/{sample}/logs/power_sim/fit_negbinom_distr_{chr}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/fit_negbinom_distr.R"

# perform power simulations by submitting one iteration at a time
rule perform_power_simulations_chr:
  input: "resources/{sample}/perturb_sce_disp.{chr}.rds"
  output:
    temp("results/{sample}/power_sim/{chr}/{chr}_rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    norm = config["power_simulations"]["norm"],
    cell_batches = config["diff_expr"]["cell_batches"],
    genes_iter = False
  log: "results/{sample}/logs/power_sim/{chr}/power_sim_{chr}_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "24G",
    time = "24:00:00"
  script:
   "../scripts/power_simulations.R"

# compute power
rule compute_power_chrs:
  input:
    expand("results/{{sample}}/power_sim/{chr}/{chr}_rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]],
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/{sample}/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz"
  wildcard_constraints:
    sample = split_samples_wildcards
  params:
    p_adj_method = config["diff_expr"]["padj_method"],
    pval_threshold = config["diff_expr"]["padj_threshold"]
  conda: "../envs/analyze_crispr_screen.yml"
  script:
    "../scripts/compute_power.R"
