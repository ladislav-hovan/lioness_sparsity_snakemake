#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

from collections import Counter

### Functions ###
def resample_reads_control(
    series: pd.Series,
    rng: np.random.Generator,
) -> pd.Series:
    # Currently too slow to be used

    target_n = series.sum()
    probs = series.values / target_n

    selection = rng.choice(series.index, target_n, p=probs)
    counts = Counter(selection)

    res_series = pd.Series(counts.values(), index=counts.keys()
        ).reindex_like(series).fillna(0).astype(int)

    return res_series

def adjust_reads_control(
    series: pd.Series,
    rng: np.random.Generator,
    proportion: float = 0.1,
) -> pd.Series:

    offsets = series.apply(lambda x: round(rng.normal(0, x * proportion)))
    adj_series = series + offsets
    adj_series = adj_series.apply(lambda x: x if x > 0 else 1)

    return adj_series

### Main body ###
rng = np.random.default_rng(snakemake.config['random_seed'])

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0, header=None)

for attempt in range(snakemake.config['n_repeats']):
    temp = {}
    for pos,col in enumerate(df.columns):
        temp[pos] = adjust_reads_control(df[col], rng)
    df_temp = pd.DataFrame({i: temp[i] for i in range(df.shape[1])})
    df_temp.to_csv(snakemake.output[attempt], sep='\t', header=False)