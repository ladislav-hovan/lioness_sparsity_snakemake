configfile: 'config.yaml'

from snakemake_gpu_manager import GpuManager, allocate_gpus

gpu_manager = GpuManager(config['gpu_ids'])

# TODO: Variables for input and output

rule all:
    input:
        expand(
            'plots/{data_type}_{method}_correlation.png',
            data_type=config['data_types'],
            method=config['methods'],
        )

# TODO: Convert scripts to Python functions, call from here as run block

rule filter_expression_and_priors:
    input:
        'input/expression_tcga_brca_raw_nohead.tsv',
        'input/motif_prior.tsv',
        'input/ppi_prior.tsv',
    output:
        'filtered_input/expression.tsv',
        'filtered_input/motif_prior.tsv',
        'filtered_input/ppi_prior.tsv',
    script:
        'scripts/filter_expression_and_priors.py'

rule sparsify_reads:
    input:
        'filtered_input/expression.tsv'
    output:
        expand(
            'sparse_expression/{{sparsity}}/expression_{repeat}.tsv',
            repeat=range(config['n_repeats']),
        )
    script:
        'scripts/sparsify_reads.py'

# TODO: Make CPU-only compatible via logical statements

rule calculate_lioness_networks:
    input:
        'sparse_expression/{sparsity}/expression_{repeat}.tsv',
        'filtered_input/motif_prior.tsv',
        'filtered_input/ppi_prior.tsv',
    output:
        'lioness_networks/{sparsity}/{repeat}/lioness.tsv'
    resources:
        gpus = 1
    params:
        gpu_manager = gpu_manager
    script:
        'scripts/calculate_lioness_networks.py'

rule calculate_baseline_networks:
    input:
        'filtered_input/expression.tsv',
        'filtered_input/motif_prior.tsv',
        'filtered_input/ppi_prior.tsv',
    output:
        'lioness_networks/baseline/lioness.tsv'
    resources:
        gpus = 1
    params:
        gpu_manager = gpu_manager
    script:
        'scripts/calculate_lioness_networks.py'

rule convert_lioness_tsv_to_feather:
    input:
        'lioness_networks/{path_to_dir}/lioness.tsv'
    output:
        'lioness_networks/{path_to_dir}/lioness.feather',
        'lioness_networks/{path_to_dir}/indegrees.feather',
        'lioness_networks/{path_to_dir}/outdegrees.feather',
    script:
        'scripts/process_lioness_tsv.py'

rule calculate_expression_correlations:
    input:
        'filtered_input/expression.tsv',
        expand(
            'sparse_expression/{{sparsity}}/expression_{repeat}.tsv',
            repeat=range(config['n_repeats']),
        ),
    output:
        'expression_correlations/{sparsity}/pearson.tsv',
        'expression_correlations/{sparsity}/spearman.tsv',
    script:
        'scripts/calculate_correlations.py'

rule calculate_coexpression_error:
    input:
        'filtered_input/expression.tsv',
        expand(
            'sparse_expression/{{sparsity}}/expression_{repeat}.tsv',
            repeat=range(config['n_repeats']),
        ),
    output:
        'coexpression_error/{sparsity}/abs_error.feather'
    script:
        'scripts/calculate_coexpression_error.py'

rule calculate_indegree_correlations:
    input:
        'lioness_networks/baseline/indegrees.feather',
        expand(
            'lioness_networks/{{sparsity}}/{repeat}/indegrees.feather',
            repeat=range(config['n_repeats']),
        ),
    output:
        'indegree_correlations/{sparsity}/pearson.tsv',
        'indegree_correlations/{sparsity}/spearman.tsv',
    script:
        'scripts/calculate_correlations.py'

rule plot_correlation_by_sparsity:
    input:
        expand(
            '{{data_type}}_correlations/{sparsity}/{{method}}.tsv',
            sparsity=config['sparsity_levels'],
        )
    output:
        'plots/{data_type}_{method}_correlation.png'
    script:
        'scripts/plot_correlation_by_sparsity.py'