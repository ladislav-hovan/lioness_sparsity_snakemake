### Configuration file ###
configfile: 'config.yaml'

### Imports and setup ###
import os

global_resources = workflow.global_resources

if 'gpus' in global_resources and global_resources['gpus'] > 0:
    assert global_resources['gpus'] == len(config['gpu_ids']), (
        'The number of GPUs in resources must match the length of provided '
        'GPU ids')

    USE_GPU = True

    import cupy as cp

    gpu_devices = [cp.cuda.Device(i) for i in range(4)]

    from snakemake_gpu_manager import GpuManager, allocate_gpus

    gpu_manager = GpuManager(config['gpu_ids'])
else:
    USE_GPU = False

### Configuration variables ###
input_dir = config['input_dir']

### Input and output ###
SUMMARY_PLOTS = expand(
    os.path.join('plots', '{data_type}_{method}_correlation.png'),
    data_type=config['data_types'],
    method=config['methods'],
)

EXPRESSION_FILE = os.path.join(input_dir, config['expression_file'])
MOTIF_FILE = os.path.join(input_dir, config['motif_file'])
PPI_FILE = os.path.join(input_dir, config['ppi_file'])

# F_ prefix for filtered
F_EXPRESSION_FILE = os.path.join('filtered_input', 'expression.tsv')
F_MOTIF_FILE = os.path.join('filtered_input', 'motif_prior.tsv')
F_PPI_FILE = os.path.join('filtered_input', 'ppi_prior.tsv')

# S_ prefix for sparse
S_EXPRESSION_FILES = expand(
    os.path.join('sparse_expression', '{{sparsity}}', 
        'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)

# TODO: Variables for input and output

### Rules ###
rule all:
    input:
        SUMMARY_PLOTS
    default_target: 
        True

# TODO: Convert scripts to Python functions, call from here as run block

rule filter_expression_and_priors:
    input:
        EXPRESSION_FILE,
        MOTIF_FILE,
        PPI_FILE,
    output:
        F_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    script:
        'scripts/filter_expression_and_priors.py'

rule sparsify_reads:
    input:
        F_EXPRESSION_FILE
    output:
        S_EXPRESSION_FILES
    script:
        'scripts/sparsify_reads.py'

rule calculate_lioness_networks:
    input:
        'sparse_expression/{sparsity}/expression_{repeat}.tsv',
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        'lioness_networks/{sparsity}/{repeat}/lioness.tsv'
    threads:
        1 if USE_GPU else 4
    resources:
        gpus = int(USE_GPU)
    run:
        from scripts.calculate_lioness_networks import create_lioness_networks

        if USE_GPU:
            with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
                # Only one device is yielded here
                with gpu_devices[gpu_ids[0]]:
                    create_lioness_networks(
                        input[0], 
                        input[1], 
                        input[2], 
                        '.', 
                        'gpu', 
                        lioness_options={
                            'start': 1, 
                            'end': int(config['n_networks']),
                            'export_filename': output[0],
                        })
        else:
            create_lioness_networks(
                input[0], 
                input[1], 
                input[2], 
                '.', 
                'cpu', 
                lioness_options={
                    'start': 1, 
                    'end': int(config['n_networks']),
                    'export_filename': output[0],
                    'ncores': threads
                })

rule calculate_baseline_networks:
    input:
        'filtered_input/expression.tsv',
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        'lioness_networks/baseline/lioness.tsv'
    threads:
        1 if USE_GPU else 4
    resources:
        gpus = int(USE_GPU)
    run:
        from scripts.calculate_lioness_networks import create_lioness_networks

        if USE_GPU:
            with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
                # Only one device is yielded here
                with gpu_devices[gpu_ids[0]]:
                    create_lioness_networks(
                        input[0], 
                        input[1], 
                        input[2], 
                        '.', 
                        'gpu', 
                        lioness_options={
                            'start': 1, 
                            'end': int(config['n_networks']),
                            'export_filename': output[0],
                        })
        else:
            create_lioness_networks(
                input[0],
                input[1],
                input[2],
                '.', 
                'cpu', 
                lioness_options={
                    'start': 1, 
                    'end': int(config['n_networks']),
                    'export_filename': output[0],
                    'ncores': threads
                })

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