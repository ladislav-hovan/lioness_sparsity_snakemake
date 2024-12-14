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
    threads:
        config['lioness_threads']
    run:
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
        config['lioness_threads']
    run:
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
        config['lioness_threads']
    run:
        calculate_lioness_networks(input[0], input[1], input[2], output[0],
            config['n_networks'], None, threads)

rule calculate_coexpression_networks:
    input:
        S_ANY_EXPRESSION_FILE
    output:
        S_COEXPR_NET_SINGLE
    threads:
        config['lioness_threads']
    resources:
        mem_gb = 150
    run:
        calculate_coexpression_networks(input[0], output[0],
            config['n_networks'], None, threads)

rule calculate_baseline_coexpression_networks:
    input:
        F_TRANS_EXPRESSION_FILE
    output:
        BL_COEXPR_NET
    threads:
        config['lioness_threads']
    resources:
        mem_gb = 150
    run:
        calculate_coexpression_networks(input[0], output[0],
            config['n_networks'], None, threads)

rule calculate_control_coexpression_networks:
    input:
        C_TRANS_EXPRESSION_FILE
    output:
        C_COEXPR_NET_SINGLE
    threads:
        config['lioness_threads']
    resources:
        mem_gb = 150
    run:
        calculate_coexpression_networks(input[0], output[0],
            config['n_networks'], None, threads)