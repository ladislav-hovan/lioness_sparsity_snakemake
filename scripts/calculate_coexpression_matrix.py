#!/usr/bin/env python

### Imports ###
import pandas as pd

from lib.functions_io import load_file

### Main body ###
df = load_file(snakemake.input[0])
coexpr_df = df.T.corr()
coexpr_df.reset_index().to_feather(snakemake.output[0])