#!/usr/bin/env python

### Imports ###
import pandas as pd

### Main body ###
df = pd.read_csv(snakemake.input[0], sep='\t', index_col=['tf', 'gene'])
df.reset_index().to_feather(snakemake.output[0])

ind_df = df.groupby('gene').sum()
ind_df.reset_index().to_feather(snakemake.output[1])

outd_df = df.groupby('tf').sum()
outd_df.reset_index().to_feather(snakemake.output[2])