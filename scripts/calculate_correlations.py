#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd
import scipy.stats as stats

from lib.functions import load_file

### Main body ###
baseline_df = load_file(snakemake.input[0])

pearson_r = []
spearman_r = []
for file in snakemake.input[1:]:
    df = load_file(file).fillna(0)
    for col in df.columns:
        pearson_r.append(stats.pearsonr(baseline_df[col], 
            df[col])[0])
        spearman_r.append(stats.spearmanr(baseline_df[col], 
            df[col])[0])
        
np.save(snakemake.output[0], pearson_r)
np.save(snakemake.output[1], spearman_r)