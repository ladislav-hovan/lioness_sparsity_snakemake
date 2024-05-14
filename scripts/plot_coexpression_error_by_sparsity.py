#!/usr/bin/env python

### Imports ###
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from typing import Tuple

# TODO: Make a common function

### Functions ###
def plot_boxplots(
    data: pd.DataFrame
) -> Tuple[plt.Figure, plt.Axes]:
    
    fig,ax = plt.subplots(figsize=(8,4))

    ax.boxplot(data, labels=data.columns)
    ax.set_ylim(0, 2)
    ax.set_xlabel('Sparsity % of expressed genes')
    method_name = snakemake.wildcards['method']
    method_name = method_name[0].upper() + method_name[1:]
    ax.set_ylabel(f'Absolute coexpression error compared to no sparsity')
    ax.grid(axis='y')

    return fig,ax

### Main body ###
df = pd.DataFrame()
for sparsity,file in zip(snakemake.config['sparsity_levels'], snakemake.input):
    df[sparsity] = np.load(file)

fig,ax = plot_boxplots(df)
fig.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')