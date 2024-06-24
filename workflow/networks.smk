from scripts.calculate_lioness_networks import calculate_lioness_networks

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
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_lioness_networks(input[0], input[1], input[2], 
                    output[0], config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_lioness_networks(input[0], input[1], input[2], output[0], 
                config['n_networks'], None, threads)

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
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_lioness_networks(input[0], input[1], input[2], 
                    output[0], config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_lioness_networks(input[0], input[1], input[2], output[0], 
                config['n_networks'], None, threads)

rule calculate_control_networks:
    input:
        C_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        C_LIONESS_FEATHER_SINGLE
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        gpus = int(USE_GPU)
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_lioness_networks(input[0], input[1], input[2], 
                    output[0], config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_lioness_networks(input[0], input[1], input[2], output[0], 
                config['n_networks'], None, threads)

rule calculate_coexpression_networks:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        os.path.join('coexpression_networks', '{transform}', 
            '{method}', '{sparsity}', '{repeat}', 'lioness.feather')
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    script:
        'scripts/calculate_coexpression_networks.py'

rule calculate_baseline_coexpression_networks:
    input:
        F_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        os.path.join('coexpression_networks', '{transform}', 
            'baseline', 'lioness.feather')
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    script:
        'scripts/calculate_coexpression_networks.py'  

rule calculate_control_coexpression_networks:
    input:
        C_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        os.path.join('coexpression_networks', '{transform}', 
            'control', '{repeat}', 'lioness.feather')
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    script:
        'scripts/calculate_coexpression_networks.py'

rule calculate_baseline_coexpression_matrix:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        os.path.join('coexpression_matrices', '{transform}', 'baseline', 
            'coexpression.feather')
    script:
        'scripts/calculate_coexpression_matrix.py'

rule calculate_control_coexpression_matrix:
    input:
        C_TRANS_EXPRESSION_FILE
    output:
        os.path.join('coexpression_matrices', '{transform}', 'control', 
            'coexpression_{repeat}.feather')
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

# rule convert_lioness_tsv_to_feather:
#     input:
#         ANY_LIONESS_TSV
#     output:
#         ANY_LIONESS_FEATHER,
#         ANY_INDEGREE_FEATHER,
#         ANY_OUTDEGREE_FEATHER,
#     script:
#         'scripts/process_lioness_tsv.py'