rule calculate_expression_correlations:
    input:
        F_TRANS_EXPRESSION_FILE,
        S_ANY_EXPRESSION_FILES,
    output:
        EXPR_CORR_PEARSON,
        EXPR_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_expression_correlations_control:
    input:
        F_TRANS_EXPRESSION_FILE,
        C_TRANS_EXPRESSION_FILES,
    output:
        C_EXPR_CORR_PEARSON,
        C_EXPR_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

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
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_correlations_control:
    input:
        os.path.join('coexpression_matrices', '{transform}', 'baseline', 
            'coexpression.feather'),
        expand(
            os.path.join('coexpression_matrices', '{{transform}}', 
                'control', 'coexpression_{repeat}.feather'),
            repeat=range(config['n_repeats']),
        ),
    output:
        os.path.join('coexpression_correlations', '{transform}', 
            'control', 'pearson.npy'),
        os.path.join('coexpression_correlations', '{transform}', 
            'control', 'spearman.npy'),
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_network_correlations:
    input:
        os.path.join('coexpression_networks', '{transform}', 'baseline', 
            'lioness.feather'),
        expand(
            os.path.join('coexpression_networks', '{{transform}}', 
                '{{method}}', '{{sparsity}}', '{repeat}', 'lioness.feather'),
            repeat=range(config['n_repeats']),
        ),
    output:
        os.path.join('coexpression_network_correlations', '{transform}', 
            '{method}', '{sparsity}', 'pearson.npy'),
        os.path.join('coexpression_network_correlations', '{transform}', 
            '{method}', '{sparsity}', 'spearman.npy'),
    resources:
        mem_gb = 150
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_network_correlations_control:
    input:
        os.path.join('coexpression_networks', '{transform}', 'baseline', 
            'lioness.feather'),
        expand(
            os.path.join('coexpression_networks', '{{transform}}', 'control', 
                '{repeat}', 'lioness.feather'),
            repeat=range(config['n_repeats']),
        ),
    output:
        os.path.join('coexpression_network_correlations', '{transform}', 
            'control', 'pearson.npy'),
        os.path.join('coexpression_network_correlations', '{transform}', 
            'control', 'spearman.npy'),
    resources:
        mem_gb = 150
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_error:
    input:
        F_TRANS_EXPRESSION_FILE,
        S_ANY_EXPRESSION_FILES,
    output:
        COEXPR_ERROR
    script:
        os.path.join('..', 'scripts', 'calculate_coexpression_error.py')

rule calculate_indegree_correlations:
    input:
        BL_INDEGREE_FEATHER,
        S_INDEGREES_FEATHER,
    output:
        IND_CORR_PEARSON,
        IND_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_indegree_correlations_control:
    input:
        BL_INDEGREE_FEATHER,
        C_INDEGREES_FEATHER,
    output:
        C_IND_CORR_PEARSON,
        C_IND_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_edge_correlations:
    input:
        BL_LIONESS_FEATHER,
        S_LIONESS_FEATHER,
    output:
        EDGE_CORR_PEARSON,
        EDGE_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_edge_correlations_control:
    input:
        BL_LIONESS_FEATHER,
        C_LIONESS_FEATHER,
    output:
        C_EDGE_CORR_PEARSON,
        C_EDGE_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')