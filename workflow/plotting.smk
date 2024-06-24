rule plot_correlation_by_sparsity:
    input:
        CONTROL_FILES,
        CORR_FILES,
    output:
        CORR_PLOT
    script:
        os.path.join('..', 'scripts', 'plot_correlation_by_sparsity.py')

rule plot_coexpression_error_by_sparsity:
    input:
        COEXPR_ERROR_FILES
    output:
        COEXPR_ERROR_PLOT
    script:
        os.path.join('..', 'scripts', 'plot_coexpression_error_by_sparsity.py')