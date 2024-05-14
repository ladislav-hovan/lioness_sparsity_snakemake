#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd
import scipy.stats as stats

from lib.functions_io import load_file

# TODO: Make a common function

### Main body ###
baseline_df = load_file(snakemake.input[0])
bl_coexpr_df = baseline_df.T.corr()

coexpr_err = []
for file in snakemake.input[1:]:
    df = load_file(file)
    coexpr_df = df.T.corr()
    diff_df = coexpr_df - bl_coexpr_df
    for pos,row in enumerate(diff_df.columns):
        for col in diff_df.columns[pos+1:]:
            coexpr_err.append(abs(diff_df.loc[row,col]))

np.save(snakemake.output[0], coexpr_err)