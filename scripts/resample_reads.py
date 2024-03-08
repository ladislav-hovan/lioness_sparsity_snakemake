#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

from collections import Counter

### Functions ###
def resample_reads(
    series: pd.Series,
    sparsity: float,
    rng: np.random.Generator,
) -> pd.Series:
    
    target_n = int((1.0 - sparsity) * len(series))
    probs = series.values / series.sum()
    counts = Counter()

    while len(counts) < target_n:
        to_gen = target_n - len(counts)
        selection = rng.choice(series.index, to_gen, p=probs)
        counts.update(selection)

    sparse_series = pd.Series(counts.values(), index=counts.keys()
        ).reindex_like(series).fillna(0).astype(int)

    return sparse_series

### Main body ###
rng = np.random.default_rng(snakemake.config['random_seed'])

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)

for attempt in range(snakemake.config['n_repeats']):
    temp = {}
    for pos,col in enumerate(df.columns):
        temp[pos] = resample_reads(df[col],
            0.01 * float(snakemake.wildcards['sparsity']), rng)
    df_temp = pd.DataFrame({i: temp[i] for i in range(df.shape[1])})
    df_temp.to_csv(snakemake.output[attempt], sep='\t')