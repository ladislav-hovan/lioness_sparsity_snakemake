#!/usr/bin/env python

### Imports ###
import bisect

import numpy as np
import pandas as pd

from collections import Counter

### Functions ###
def downsample_reads(
    series: pd.Series,
    sparsity: float,
    rng: np.random.Generator,
    min_chunk: int = 10000,
) -> pd.Series:
    
    target_n = int((1.0 - sparsity) * len(series))
    orig_counts = Counter({k: v for k,v in series.items()})
    probs = series.values / series.sum()
    counts = Counter()

    while len(counts) < target_n:
        to_gen = target_n - len(counts)
        selection = rng.choice(series.index, max(to_gen, min_chunk), p=probs)
        temp = counts.copy()
        temp.update(selection)
        for k,v in temp.items():
            if v > orig_counts[k]:
                temp[k] = orig_counts[k]
        if len(temp) < target_n:
            counts = temp
        else:
            index = 0
            while len(counts) < target_n:
                counts[selection[index]] += 1
                index += 1
        rem_counts = orig_counts - counts
        probs = np.array([v for v in rem_counts.values()]) / rem_counts.total()

    sparse_series = pd.Series(counts.values(), index=counts.keys()
        ).reindex_like(series).fillna(0).astype(int)

    return sparse_series

### Main body ###
rng = np.random.default_rng(snakemake.config['random_seed'])

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0, header=None)

for attempt in range(snakemake.config['n_repeats']):
    temp = {}
    for pos,col in enumerate(df.columns):
        temp[pos] = downsample_reads(df[col],
            0.01 * float(snakemake.wildcards['sparsity']), rng)
    df_temp = pd.DataFrame({i: temp[i] for i in range(df.shape[1])})
    df_temp.to_csv(snakemake.output[attempt], sep='\t', header=False)