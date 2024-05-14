#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

### Main body ###
df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0, header=None)

df_mod = df.add(1).apply(np.log)

df_mod.to_csv(snakemake.output[0], sep='\t', header=False)