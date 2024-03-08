#!/usr/bin/env python

### Imports ###
import pandas as pd
import numpy as np

import scipy.stats as stats

# TODO: Make a common function

### Functions ###
def load_file(
    filename: str,
) -> pd.DataFrame:
    
    extension = filename.split('.')[-1]

    if extension == 'tsv':
        # Expression files
        df = pd.read_csv(filename, sep='\t', index_col=0, header=None)
    elif extension == 'feather':
        # Indegree files
        df = pd.read_feather(filename).set_index('gene')
    else:
        raise ValueError(f'Invalid extension for: {filename}, '
            'expected tsv or feather')

    return df

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

np.savetxt(snakemake.output[0], coexpr_err)