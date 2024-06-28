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
        BL_COEXPR_MAT,
        COEXPR_MATS,
    output:
        COEXPR_MAT_CORR_PEARSON,
        COEXPR_MAT_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_correlations_control:
    input:
        BL_COEXPR_MAT,
        C_COEXPR_MATS,
    output:
        C_COEXPR_MAT_CORR_PEARSON,
        C_COEXPR_MAT_CORR_SPEARMAN,
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_network_correlations:
    input:
        BL_COEXPR_NET,
        COEXPR_NETS,
    output:
        COEXPR_NET_CORR_PEARSON,
        COEXPR_NET_CORR_SPEARMAN,
    resources:
        mem_gb = 150
    script:
        os.path.join('..', 'scripts', 'calculate_correlations.py')

rule calculate_coexpression_network_correlations_control:
    input:
        BL_COEXPR_NET,
        C_COEXPR_NETS,
    output:
        C_COEXPR_NET_CORR_PEARSON,
        C_COEXPR_NET_CORR_SPEARMAN,
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