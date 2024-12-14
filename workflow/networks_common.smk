from scripts.calculate_coexpression_matrix import (
    calculate_coexpression_matrix)
from scripts.calculate_degrees import calculate_degrees

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