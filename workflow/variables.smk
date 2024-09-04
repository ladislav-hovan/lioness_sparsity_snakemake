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

ZERO_GENES_MARK = os.path.join('sparse_expression', '{method}', '{sparsity}',
    'mark_zero_{repeat}.tsv')
TOP_HALF_GENES_MARK = os.path.join('sparse_expression', '{method}',
    '{sparsity}', 'mark_top_half_{repeat}.tsv')

S_ANY_EXPRESSION_FILE = os.path.join('sparse_expression', '{transform}',
    '{method}', '{sparsity}', 'expression_{repeat}.tsv')
S_LIONESS_FEATHER_SINGLE = os.path.join('lioness_networks', '{transform}',
    '{method}', '{sparsity}', '{repeat}', 'lioness.feather')
BL_LIONESS_FEATHER = os.path.join('lioness_networks', '{transform}',
    'baseline', 'lioness.feather')
C_LIONESS_FEATHER_SINGLE = os.path.join('lioness_networks', '{transform}',
    'control', '{repeat}', 'lioness.feather')

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

S_COEXPR_NET_SINGLE = os.path.join('coexpression_networks', '{transform}',
    '{method}', '{sparsity}', '{repeat}', 'lioness.feather')
BL_COEXPR_NET = os.path.join('coexpression_networks', '{transform}',
    'baseline', 'lioness.feather')
C_COEXPR_NET_SINGLE = os.path.join('coexpression_networks', '{transform}',
    'control', '{repeat}', 'lioness.feather')

S_COEXPR_MAT_SINGLE = os.path.join('coexpression_matrices', '{transform}',
    '{method}', '{sparsity}', 'coexpression_{repeat}.feather')
BL_COEXPR_MAT = os.path.join('coexpression_matrices', '{transform}',
    'baseline', 'coexpression.feather')
C_COEXPR_MAT_SINGLE = os.path.join('coexpression_matrices', '{transform}',
    'control', 'coexpression_{repeat}.feather')

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

COEXPR_MATS = expand(
    os.path.join('coexpression_matrices', '{{transform}}',
        '{{method}}', '{{sparsity}}', 'coexpression_{repeat}.feather'),
    repeat=range(config['n_repeats']),
)
C_COEXPR_MATS = expand(
    os.path.join('coexpression_matrices', '{{transform}}',
        'control', 'coexpression_{repeat}.feather'),
    repeat=range(config['n_repeats']),
)

COEXPR_MAT_CORR_PEARSON = os.path.join('coexpression_correlations',
    '{transform}', '{method}', '{sparsity}', 'pearson.npy')
COEXPR_MAT_CORR_SPEARMAN = os.path.join('coexpression_correlations',
    '{transform}', '{method}', '{sparsity}', 'spearman.npy')

C_COEXPR_MAT_CORR_PEARSON = os.path.join('coexpression_correlations',
    '{transform}', 'control', 'pearson.npy')
C_COEXPR_MAT_CORR_SPEARMAN = os.path.join('coexpression_correlations',
    '{transform}', 'control', 'spearman.npy')

COEXPR_NETS = expand(
    os.path.join('coexpression_networks', '{{transform}}',
        '{{method}}', '{{sparsity}}', '{repeat}', 'lioness.feather'),
    repeat=range(config['n_repeats']),
)
C_COEXPR_NETS = expand(
    os.path.join('coexpression_networks', '{{transform}}', 'control',
        '{repeat}', 'lioness.feather'),
    repeat=range(config['n_repeats']),
)

COEXPR_NET_CORR_PEARSON = os.path.join('coexpression_network_correlations',
    '{transform}', '{method}', '{sparsity}', 'pearson.npy')
COEXPR_NET_CORR_SPEARMAN = os.path.join('coexpression_network_correlations',
    '{transform}', '{method}', '{sparsity}', 'spearman.npy')

C_COEXPR_NET_CORR_PEARSON = os.path.join('coexpression_network_correlations',
    '{transform}', 'control', 'pearson.npy')
C_COEXPR_NET_CORR_SPEARMAN = os.path.join('coexpression_network_correlations',
    '{transform}', 'control', 'spearman.npy')

C_CORR_FILE = os.path.join('{data_type}_correlations', '{transform}',
    'control', '{corr}.npy')
CORR_FILES = expand(
    os.path.join('{{data_type}}_correlations', '{{transform}}', '{{method}}',
        '{sparsity}', '{{corr}}.npy'),
    sparsity=config['sparsity_levels'],
)
CORR_PLOT = os.path.join('plots', '{transform}', '{method}', '{data_type}',
    '{corr}_correlation.png')

C_COEXPR_ERROR_FILE = os.path.join('coexpression_error', '{transform}',
    'control', 'abs_error.npy')
COEXPR_ERROR_FILES = expand(
    os.path.join('coexpression_error', '{{transform}}', '{{method}}',
        '{sparsity}', 'abs_error.npy'),
    sparsity=config['sparsity_levels'],
)
COEXPR_ERROR_PLOT = os.path.join('plots', '{transform}',
    '{method}_coexpr_error.png')