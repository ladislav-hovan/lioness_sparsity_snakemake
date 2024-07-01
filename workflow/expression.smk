from scripts.downsample_reads import downsample_reads
from scripts.filter_expression_and_priors import filter_expression_and_priors
from scripts.generate_control_expression import generate_control_expression
from scripts.resample_reads import resample_reads
from scripts.rescale_expression_log1p import rescale_expression_log1p

rule filter_expression_and_priors:
    input:
        EXPRESSION_FILE,
        MOTIF_FILE,
        PPI_FILE,
    output:
        F_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    run:
        filter_expression_and_priors(input[0], input[1], input[2],
            output[0], output[1], output[2])

rule resample_reads:
    input:
        F_EXPRESSION_FILE
    output:
        S_RS_EXPRESSION_FILES
    params:
        n_repeats = config['n_repeats'],
        seed = config['random_seed'] if 'random_seed' in config else None,
    run:
        resample_reads(input[0], output, float(wildcards['sparsity']),
            params['n_repeats'], params['seed'])

rule downsample_reads:
    input:
        F_EXPRESSION_FILE
    output:
        S_DS_EXPRESSION_FILES
    params:
        n_repeats = config['n_repeats'],
        seed = config['random_seed'] if 'random_seed' in config else None,
    run:
        downsample_reads(input[0], output, float(wildcards['sparsity']),
            params['n_repeats'], params['seed'])

rule generate_control_expression:
    input:
        F_EXPRESSION_FILE
    output:
        C_EXPRESSION_FILES
    params:
        n_repeats = config['n_repeats'],
        seed = config['random_seed'] if 'random_seed' in config else None,
    run:
        generate_control_expression(input[0], output, params['n_repeats'],
            params['seed'])

rule mark_zero_expression_genes:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        ZERO_GENES_MARK
    script:
        # TODO: Implement
        'scripts/mark_zero_expression_genes.py'

rule mark_top_half_expression_genes:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        TOP_HALF_GENES_MARK
    script:
        # TODO: Implement
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
    run:
        rescale_expression_log1p(input[0], output[0])