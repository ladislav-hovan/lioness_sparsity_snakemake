#!/usr/bin/env python

### Imports ###
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from typing import Tuple

### Functions ###
def plot_boxplots(
    data: pd.DataFrame
) -> Tuple[plt.Figure, plt.Axes]:
    
    fig,ax = plt.subplots(figsize=(8,4))

    ax.boxplot(data, labels=data.columns)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Sparsity % of expressed genes')
    corr_name = snakemake.wildcards['corr']
    corr_name = corr_name[0].upper() + corr_name[1:]
    ax.set_ylabel(f'{corr_name} R compared to no sparsity')
    ax.grid(axis='y')

    return fig,ax

### Main body ###
df = pd.DataFrame()
for sparsity,file in zip([0] + snakemake.config['sparsity_levels'], 
    snakemake.input):
    df[sparsity] = np.load(file)
    df[sparsity] = df[sparsity].fillna(df[sparsity].mean())

fig,ax = plot_boxplots(df)
fig.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')