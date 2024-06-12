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