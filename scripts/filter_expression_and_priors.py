#!/usr/bin/env python

### Imports ###
import pandas as pd

### Main body ###
# Remove genes with zero expression
df = pd.read_csv(snakemake.input[0], sep='\t', header=None, index_col=0)
mask = ((df > 0).sum(axis=1) == df.shape[1])
df_filt = df.loc[mask]

# Ensure correspondence between motif prior and expression
motif_prior = pd.read_csv(snakemake.input[1], sep='\t', 
    names=['tf', 'gene', 'edge'])
df_genes = set(df_filt.index)
motif_prior_genes = set(motif_prior['gene'].unique())
common_genes = df_genes.intersection(motif_prior_genes)
df_common = df_filt.loc[list(common_genes)]
mask = motif_prior['gene'].isin(common_genes)
motif_prior_common = motif_prior.loc[mask]

# Ensure all TFs from PPI prior in the motif prior
ppi_prior = pd.read_csv(snakemake.input[2], sep='\t', 
    names=['tf1', 'tf2', 'edge'])
all_tfs = set(ppi_prior['tf1']).union(set(ppi_prior['tf2']))
motif_prior_tfs = set(motif_prior_common['tf'].unique())
extra_tfs = all_tfs - motif_prior_tfs
mask = ~ppi_prior['tf1'].isin(extra_tfs) & ~ppi_prior['tf2'].isin(extra_tfs)
ppi_prior_common = ppi_prior.loc[mask]

# Save everything properly
df_common.sort_index().to_csv(snakemake.output[0], sep='\t', header=False)
motif_prior_common.sort_values(['tf', 'gene']).to_csv(snakemake.output[1], 
    sep='\t', index=False, header=False)
ppi_prior_common.sort_values(['tf1', 'tf2']).to_csv(snakemake.output[2], 
    sep='\t', index=False, header=False)