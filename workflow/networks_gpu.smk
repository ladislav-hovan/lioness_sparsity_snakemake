from scripts.calculate_coexpression_networks import (
    calculate_coexpression_networks)
from scripts.calculate_lioness_networks import calculate_lioness_networks

rule calculate_lioness_networks:
    input:
        S_ANY_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        S_LIONESS_FEATHER_SINGLE
    resources:
        gpus = 1
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_lioness_networks(input[0], input[1], input[2],
                output[0], config['n_networks'], gpu_ids[0], 1)

rule calculate_baseline_networks:
    input:
        F_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        BL_LIONESS_FEATHER
    resources:
        gpus = 1
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_lioness_networks(input[0], input[1], input[2],
                output[0], config['n_networks'], gpu_ids[0], 1)

rule calculate_control_networks:
    input:
        C_TRANS_EXPRESSION_FILE,
        F_MOTIF_FILE,
        F_PPI_FILE,
    output:
        C_LIONESS_FEATHER_SINGLE
    resources:
        gpus = 1
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_lioness_networks(input[0], input[1], input[2],
                output[0], config['n_networks'], gpu_ids[0], 1)

rule calculate_coexpression_networks:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        S_COEXPR_NET_SINGLE
    resources:
        mem_gb = 150,
        gpus = 1,
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], gpu_ids[0], 1)

rule calculate_baseline_coexpression_networks:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        BL_COEXPR_NET
    resources:
        mem_gb = 150,
        gpus = 1,
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], gpu_ids[0], 1)

rule calculate_control_coexpression_networks:
    input:
        C_TRANS_EXPRESSION_FILE
    output:
        C_COEXPR_NET_SINGLE
    resources:
        mem_gb = 150,
        gpus = 1,
    run:
        with allocate_gpus(gpu_manager, resources['gpus']) as gpu_ids:
            calculate_coexpression_networks(input[0], output[0],
                config['n_networks'], gpu_ids[0], 1)