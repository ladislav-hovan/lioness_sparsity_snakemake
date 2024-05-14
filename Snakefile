### Configuration file ###
configfile: 'config.yaml'
# configfile: os.path.join('tests', 'test_config.yaml')

### Imports and setup ###
import os

global_resources = workflow.global_resources

if 'gpus' in global_resources and global_resources['gpus'] > 0:
    assert global_resources['gpus'] == len(config['gpu_ids']), (
        'The number of GPUs in resources must match the length of provided '
        'GPU ids')

    USE_GPU = True

    import cupy as cp

    gpu_devices = [cp.cuda.Device(i) for i in 
        range(cp.cuda.runtime.getDeviceCount())]

    from lib.gpu_manager import GpuManager, allocate_gpus

    gpu_manager = GpuManager(config['gpu_ids'])
else:
    USE_GPU = False

    gpu_manager = None
    allocate_gpus = None

### Configuration variables ###
input_dir = config['input_dir']

# TODO: Split into multiple snakefiles

### Input and output ###
SUMMARY_CORR_PLOTS = expand(
    os.path.join('plots', '{transform}', '{method}', '{data_type}', 
        '{corr}_correlation.png'),
    transform=config['transformations'],
    method=config['sparsifying_methods'],
    data_type=config['data_types'],
    corr=config['corr_types'],
)
COEXPR_ERROR_PLOTS = expand(
    os.path.join('plots', '{transform}', '{method}_coexpr_error.png'),
    transform=config['transformations'],
    method=config['sparsifying_methods'],
)

EXPRESSION_FILE = os.path.join(input_dir, config['expression_file'])
MOTIF_FILE = os.path.join(input_dir, config['motif_file'])
PPI_FILE = os.path.join(input_dir, config['ppi_file'])

# F_ prefix for filtered
F_EXPRESSION_FILE = os.path.join('filtered_input', 'raw_expression', 
    'expression.tsv')
F_MOTIF_FILE = os.path.join('filtered_input', 'raw_expression', 
    'motif_prior.tsv')
F_PPI_FILE = os.path.join('filtered_input', 'raw_expression', 'ppi_prior.tsv')

# S_ prefix for sparse
# RS_ prefix for resample
S_RS_EXPRESSION_FILES = expand(
    os.path.join('sparse_expression', 'raw_expression', 'resample', 
        '{{sparsity}}', 'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)
S_RS_EXPRESSION_FILE = os.path.join('sparse_expression', 'raw_expression', 
    'resample', '{sparsity}', 'expression_{repeat}.tsv')

# DS_ prefix for downsample
S_DS_EXPRESSION_FILES = expand(
    os.path.join('sparse_expression', 'raw_expression', 'downsample', 
        '{{sparsity}}', 'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)
S_DS_EXPRESSION_FILE = os.path.join('sparse_expression', 'raw_expression',
    'downsample', '{sparsity}', 'expression_{repeat}.tsv')

# C_ prefix for control
C_EXPRESSION_FILES = expand(
    os.path.join('control_expression', 'raw_expression', 
        'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)
C_TRANS_EXPRESSION_FILES = expand(
    os.path.join('control_expression', '{{transform}}', 
        'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)
C_TRANS_EXPRESSION_FILE = os.path.join('control_expression', '{transform}', 
    'expression_{repeat}.tsv')

S_ANY_EXPRESSION_FILE = os.path.join('sparse_expression', '{transform}', 
    '{method}', '{sparsity}', 'expression_{repeat}.tsv')
S_LIONESS_TSV = os.path.join('lioness_networks', '{transform}', 
    '{method}', '{sparsity}', '{repeat}', 'lioness.tsv')
BL_LIONESS_TSV = os.path.join('lioness_networks', '{transform}', 'baseline', 
    'lioness.tsv')
C_LIONESS_TSV = os.path.join('lioness_networks', '{transform}', 'control', 
    '{repeat}', 'lioness.tsv')
ANY_LIONESS_TSV = os.path.join('lioness_networks', '{path_to_dir}', 
    'lioness.tsv')

ANY_EXPRESSION_FILE = os.path.join('{path_to_dir_1}', 
    'raw_expression{path_to_dir_2}', 'expression{repeat}.tsv')
ANY_LOG1P_EXPRESSION_FILE = os.path.join('{path_to_dir_1}', 
    'log1p_rescaled{path_to_dir_2}', 'expression{repeat}.tsv')

ANY_LIONESS_FEATHER = os.path.join('lioness_networks', '{path_to_dir}', 
    'lioness.feather')
ANY_INDEGREE_FEATHER = os.path.join('lioness_networks', '{path_to_dir}', 
    'indegrees.feather')
ANY_OUTDEGREE_FEATHER = os.path.join('lioness_networks', '{path_to_dir}', 
    'outdegrees.feather')

S_ANY_EXPRESSION_FILES = expand(
    os.path.join('sparse_expression', '{{transform}}', '{{method}}', 
        '{{sparsity}}', 'expression_{repeat}.tsv'),
    repeat=range(config['n_repeats']),
)

F_TRANS_EXPRESSION_FILE = os.path.join('filtered_input', '{transform}', 
    'expression.tsv')

EXPR_CORR_PEARSON = os.path.join('expression_correlations', '{transform}', 
    '{method}', '{sparsity}', 'pearson.npy')
EXPR_CORR_SPEARMAN = os.path.join('expression_correlations', '{transform}', 
    '{method}', '{sparsity}', 'spearman.npy')

C_EXPR_CORR_PEARSON = os.path.join('expression_correlations', '{transform}', 
    'control', 'pearson.npy')
C_EXPR_CORR_SPEARMAN = os.path.join('expression_correlations', '{transform}', 
    'control', 'spearman.npy')

COEXPR_ERROR = os.path.join('coexpression_error', '{transform}', '{method}', 
    '{sparsity}', 'abs_error.npy')

BL_INDEGREE_FEATHER = os.path.join('lioness_networks', '{transform}', 
    'baseline', 'indegrees.feather')
C_INDEGREES_FEATHER = expand(
    os.path.join('lioness_networks', '{{transform}}', 'control', '{repeat}', 
        'indegrees.feather'),
    repeat=range(config['n_repeats']),
)
S_INDEGREES_FEATHER = expand(
    os.path.join('lioness_networks', '{{transform}}', '{{method}}', 
        '{{sparsity}}', '{repeat}', 'indegrees.feather'),
    repeat=range(config['n_repeats']),
)

IND_CORR_PEARSON = os.path.join('indegree_correlations', '{transform}', 
    '{method}', '{sparsity}', 'pearson.npy')
IND_CORR_SPEARMAN = os.path.join('indegree_correlations', '{transform}', 
    '{method}', '{sparsity}', 'spearman.npy')

C_IND_CORR_PEARSON = os.path.join('indegree_correlations', '{transform}', 
    'control', 'pearson.npy')
C_IND_CORR_SPEARMAN = os.path.join('indegree_correlations', '{transform}', 
    'control', 'spearman.npy')

BL_LIONESS_FEATHER = os.path.join('lioness_networks', '{transform}', 
    'baseline', 'lioness.feather')
C_LIONESS_FEATHER = expand(
    os.path.join('lioness_networks', '{{transform}}', 'control', '{repeat}', 
        'lioness.feather'),
    repeat=range(config['n_repeats']),
)
S_LIONESS_FEATHER = expand(
    os.path.join('lioness_networks', '{{transform}}', '{{method}}', 
        '{{sparsity}}', '{repeat}', 'lioness.feather'),
    repeat=range(config['n_repeats']),
)

EDGE_CORR_PEARSON = os.path.join('edge_correlations', '{transform}', 
    '{method}', '{sparsity}', 'pearson.npy')
EDGE_CORR_SPEARMAN = os.path.join('edge_correlations', '{transform}', 
    '{method}', '{sparsity}', 'spearman.npy')

C_EDGE_CORR_PEARSON = os.path.join('edge_correlations', '{transform}', 
    'control', 'pearson.npy')
C_EDGE_CORR_SPEARMAN = os.path.join('edge_correlations', '{transform}', 
    'control', 'spearman.npy')

CONTROL_FILES = os.path.join('{data_type}_correlations', '{transform}', 
    'control', '{corr}.npy')
CORR_FILES = expand(
    os.path.join('{{data_type}}_correlations', '{{transform}}', '{{method}}', 
        '{sparsity}', '{{corr}}.npy'),
    sparsity=config['sparsity_levels'],
)
CORR_PLOT = os.path.join('plots', '{transform}', '{method}', '{data_type}', 
    '{corr}_correlation.png')

COEXPR_ERROR_FILES = expand(
    os.path.join('coexpression_error', '{{transform}}', '{{method}}', 
        '{sparsity}', 'abs_error.npy'),
    sparsity=config['sparsity_levels'],
)
COEXPR_ERROR_PLOT = os.path.join('plots', '{transform}', 
    '{method}_coexpr_error.png')

### Rules ###
rule all:
    input:
        SUMMARY_CORR_PLOTS,
        COEXPR_ERROR_PLOTS,
        expand(
            os.path.join('coexpression_networks', '{transform}', 
                '{method}', '{sparsity}', '{repeat}', 'lioness.feather'),
            transform=config['transformations'],
            method=config['sparsifying_methods'],
            sparsity=config['sparsity_levels'],
            repeat=range(config['n_repeats']),
        ),
        expand(
            os.path.join('coexpression_correlations', '{transform}', 
                '{method}', '{sparsity}', 'pearson.npy'),
            transform=config['transformations'],
            method=config['sparsifying_methods'],
            sparsity=config['sparsity_levels'],
        ),
    default_target:
        True

# TODO: Convert scripts to Python functions, call from here as run block
# TODO: Make sure the scripts are command line callable as well
# TODO: Effectively implement one main function that's called with arguments

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

rule resample_reads:
    input:
        F_EXPRESSION_FILE
    output:
        S_RS_EXPRESSION_FILES
    script:
        'scripts/resample_reads.py'

rule downsample_reads:
    input:
        F_EXPRESSION_FILE
    output:
        S_DS_EXPRESSION_FILES
    script:
        'scripts/downsample_reads.py'

rule generate_control_expression:
    input:
        F_EXPRESSION_FILE
    output:
        C_EXPRESSION_FILES
    script:
        'scripts/generate_control_expression.py'

rule mark_zero_expression_genes:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        os.path.join('sparse_expression', '{method}', '{sparsity}', 
            'mark_zero_{repeat}.tsv')
    script:
        'scripts/mark_zero_expression_genes.py'

rule mark_top_half_expression_genes:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        os.path.join('sparse_expression', '{method}', '{sparsity}', 
            'mark_top_half_{repeat}.tsv')
    script:
        'scripts/mark_top_half_expression_genes.py'

rule rescale_filtered_expression_log1p:
    input:
        ANY_EXPRESSION_FILE
    output:
        ANY_LOG1P_EXPRESSION_FILE
    wildcard_constraints:
        # Allow them to be empty as well
        repeat='.*',
        path_to_dir_2='.*',
    script:
        'scripts/rescale_expression_log1p.py'

# TODO: Make sure that changing GPU to CPU does not trigger reruns - but how?
rule calculate_lioness_networks:
    input:
        S_ANY_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        S_LIONESS_TSV,
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        gpus = int(USE_GPU)
    params:
        gpu_manager = gpu_manager
    script:
        f'scripts/calculate_lioness_networks.py'
    
    # Future possible implementation once the run bug is fixed
    # https://github.com/snakemake/snakemake/issues/2350
    # run:
    #     from scripts.calculate_lioness_networks import create_lioness_networks

    #     if USE_GPU:
    #         with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
    #             # Only one device is yielded here
    #             with gpu_devices[gpu_ids[0]]:
    #                 create_lioness_networks(
    #                     input[0], 
    #                     input[1], 
    #                     input[2], 
    #                     '.', 
    #                     'gpu', 
    #                     lioness_options={
    #                         'start': 1, 
    #                         'end': int(config['n_networks']),
    #                         'export_filename': output[0],
    #                     })
    #     else:
    #         create_lioness_networks(
    #             input[0], 
    #             input[1], 
    #             input[2], 
    #             '.', 
    #             'cpu', 
    #             lioness_options={
    #                 'start': 1, 
    #                 'end': int(config['n_networks']),
    #                 'export_filename': output[0],
    #                 'ncores': threads
    #             })

rule calculate_baseline_networks:
    input:
        F_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        BL_LIONESS_TSV
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        gpus = int(USE_GPU)
    params:
        gpu_manager = gpu_manager
    script:
        'scripts/calculate_lioness_networks.py'

rule calculate_control_networks:
    input:
        C_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        C_LIONESS_TSV
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        gpus = int(USE_GPU)
    params:
        gpu_manager = gpu_manager
    script:
        'scripts/calculate_lioness_networks.py'

rule calculate_coexpression_networks:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        os.path.join('coexpression_networks', '{transform}', 
            '{method}', '{sparsity}', '{repeat}', 'lioness.feather')
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        gpus = int(USE_GPU)
    params:
        gpu_manager = gpu_manager
    script:
        f'scripts/calculate_coexpression_networks.py'

rule calculate_baseline_coexpression_matrix:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        os.path.join('coexpression_matrices', '{transform}', 'baseline', 
            'coexpression.feather')
    script:
        'scripts/calculate_coexpression_matrix.py'

rule calculate_coexpression_matrix:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        os.path.join('coexpression_matrices', '{transform}', '{method}', 
            '{sparsity}', 'coexpression_{repeat}.feather')
    script:
        'scripts/calculate_coexpression_matrix.py'

rule convert_lioness_tsv_to_feather:
    input:
        ANY_LIONESS_TSV
    output:
        ANY_LIONESS_FEATHER,
        ANY_INDEGREE_FEATHER,
        ANY_OUTDEGREE_FEATHER,
    script:
        'scripts/process_lioness_tsv.py'

rule calculate_expression_correlations:
    input:
        F_TRANS_EXPRESSION_FILE,
        S_ANY_EXPRESSION_FILES,
    output:
        EXPR_CORR_PEARSON,
        EXPR_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule calculate_expression_correlations_control:
    input:
        F_TRANS_EXPRESSION_FILE,
        C_TRANS_EXPRESSION_FILES,
    output:
        C_EXPR_CORR_PEARSON,
        C_EXPR_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule calculate_coexpression_correlations:
    input:
        os.path.join('coexpression_matrices', '{transform}', 'baseline', 
            'coexpression.feather'),
        expand(
            os.path.join('coexpression_matrices', '{{transform}}', 
                '{{method}}', '{{sparsity}}', 'coexpression_{repeat}.feather'),
            repeat=range(config['n_repeats']),
        ),
    output:
        os.path.join('coexpression_correlations', '{transform}', 
            '{method}', '{sparsity}', 'pearson.npy'),
        os.path.join('coexpression_correlations', '{transform}', 
            '{method}', '{sparsity}', 'spearman.npy'),
    script:
        'scripts/calculate_correlations.py'

rule calculate_coexpression_error:
    input:
        F_TRANS_EXPRESSION_FILE,
        S_ANY_EXPRESSION_FILES,
    output:
        COEXPR_ERROR
    script:
        'scripts/calculate_coexpression_error.py'

rule calculate_indegree_correlations:
    input:
        BL_INDEGREE_FEATHER,
        S_INDEGREES_FEATHER,
    output:
        IND_CORR_PEARSON,
        IND_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule calculate_indegree_correlations_control:
    input:
        BL_INDEGREE_FEATHER,
        C_INDEGREES_FEATHER,
    output:
        C_IND_CORR_PEARSON,
        C_IND_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule calculate_edge_correlations:
    input:
        BL_LIONESS_FEATHER,
        S_LIONESS_FEATHER,
    output:
        EDGE_CORR_PEARSON,
        EDGE_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule calculate_edge_correlations_control:
    input:
        BL_LIONESS_FEATHER,
        C_LIONESS_FEATHER,
    output:
        C_EDGE_CORR_PEARSON,
        C_EDGE_CORR_SPEARMAN,
    script:
        'scripts/calculate_correlations.py'

rule plot_correlation_by_sparsity:
    input:
        CONTROL_FILES,
        CORR_FILES,
    output:
        CORR_PLOT
    script:
        'scripts/plot_correlation_by_sparsity.py'

rule plot_coexpression_error_by_sparsity:
    input:
        COEXPR_ERROR_FILES
    output:
        COEXPR_ERROR_PLOT
    script:
        'scripts/plot_coexpression_error_by_sparsity.py'