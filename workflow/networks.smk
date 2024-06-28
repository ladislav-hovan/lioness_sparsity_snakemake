from scripts.calculate_coexpression_matrix import (
    calculate_coexpression_matrix)
from scripts.calculate_coexpression_networks import (
    calculate_coexpression_networks)
from scripts.calculate_degrees import calculate_degrees
from scripts.calculate_lioness_networks import calculate_lioness_networks

rule calculate_lioness_networks:
    input:
        S_ANY_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        S_LIONESS_FEATHER_SINGLE
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
        BL_LIONESS_FEATHER
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
        S_COEXPR_NET_SINGLE
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_coexpression_networks(input[0], output[0],
                    config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], None, threads)

rule calculate_baseline_coexpression_networks:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        BL_COEXPR_NET
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_coexpression_networks(input[0], output[0],
                    config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], None, threads)

rule calculate_control_coexpression_networks:
    input:
        C_TRANS_EXPRESSION_FILE
    output:
        C_COEXPR_NET_SINGLE
    threads:
        1 if USE_GPU else config['lioness_threads']
    resources:
        mem_gb = 150,
        gpus = int(USE_GPU),
    run:
        if resources['gpus'] > 0:
            with allocate_gpus(gpu_manager,
                resources['gpus']) as gpu_ids:
                calculate_coexpression_networks(input[0], output[0],
                    config['n_networks'], gpu_ids[0], threads)
        else:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], None, threads)

rule calculate_coexpression_matrix:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        S_COEXPR_MAT_SINGLE
    run:
        calculate_coexpression_matrix(input[0], output[0])

rule calculate_baseline_coexpression_matrix:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        BL_COEXPR_MAT
    run:
        calculate_coexpression_matrix(input[0], output[0])

rule calculate_control_coexpression_matrix:
    input:
        C_TRANS_EXPRESSION_FILE
    output:
        C_COEXPR_MAT_SINGLE
    run:
        calculate_coexpression_matrix(input[0], output[0])

rule calculate_degrees:
    input:
        ANY_LIONESS_FEATHER
    output:
        ANY_INDEGREE_FEATHER,
        ANY_OUTDEGREE_FEATHER,
    run:
        calculate_degrees(input[0], output[0], output[1])