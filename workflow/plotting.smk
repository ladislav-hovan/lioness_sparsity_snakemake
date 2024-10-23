from scripts.plot_coexpression_error_by_sparsity import (
    plot_coexpression_error_by_sparsity)
from scripts.plot_correlation_by_sparsity import plot_correlation_by_sparsity
from scripts.plot_correlation_by_sparsity_in_groups import plot_correlation_by_sparsity_in_groups

rule plot_correlation_by_sparsity:
    input:
        C_CORR_FILE,
        CORR_FILES,
    output:
        CORR_PLOT
    params:
        sparsity_levels = config['sparsity_levels']
    run:
        plot_correlation_by_sparsity(input,
            [0] + [float(i) for i in params['sparsity_levels']],
            wildcards['corr'], output[0])

rule plot_correlation_by_sparsity_in_groups:
    input:
        os.path.join('{data_type}_correlations', '{transform}',
            'control', '{corr}.npy'),
        expand(
            os.path.join('{{data_type}}_correlations', '{{transform}}', '{{method}}',
                '{sparsity}', '{{corr}}_{{mark}}_t.npy'),
            sparsity=config['sparsity_levels'],
        ),
        expand(
            os.path.join('{{data_type}}_correlations', '{{transform}}', '{{method}}',
                '{sparsity}', '{{corr}}_{{mark}}_f.npy'),
            sparsity=config['sparsity_levels'],
        ),
    output:
        os.path.join('plots', '{transform}', '{method}', '{data_type}',
            '{corr}_correlation_{mark}.png')
    params:
        sparsity_levels = config['sparsity_levels']
    run:
        plot_correlation_by_sparsity_in_groups(input,
            [0] + [float(i) for i in params['sparsity_levels']],
            wildcards['corr'], output[0])

rule plot_coexpression_error_by_sparsity:
    input:
        C_COEXPR_ERROR_FILE,
        COEXPR_ERROR_FILES,
    output:
        COEXPR_ERROR_PLOT
    run:
        plot_coexpression_error_by_sparsity(input,
            [0] + [float(i) for i in params['sparsity_levels']],
            output[0])