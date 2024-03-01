#!/usr/bin/env python

### Imports ###
import pandas as pd
import numpy as np

import scipy.stats as stats

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

pearson_r = []
spearman_r = []
for file in snakemake.input[1:]:
    df = load_file(file)
    for col in df.columns:
        pearson_r.append(stats.pearsonr(baseline_df[col], 
            df[col])[0])
        spearman_r.append(stats.spearmanr(baseline_df[col], 
            df[col])[0])
        
np.savetxt(snakemake.output[0], pearson_r)
np.savetxt(snakemake.output[1], spearman_r)