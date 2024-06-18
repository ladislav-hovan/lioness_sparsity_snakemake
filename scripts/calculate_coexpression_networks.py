#!/usr/bin/env python

### Imports ###
import os

import numpy as np
import pandas as pd

from netZooPy.lioness import Lioness
from netZooPy.panda import Panda

use_gpu = bool(snakemake.resources['gpus'])

if use_gpu:
    import cupy as cp

    from lib.gpu_manager import allocate_gpus

### Functions ###
def create_lioness_networks(
    expression: pd.DataFrame,
    motif_prior: pd.DataFrame,
    ppi_prior: pd.DataFrame,
    output_dir: str,
    computing: str,
    lioness_options: dict,
) -> None:

    panda_obj = Panda(expression, motif_prior, ppi_prior, computing=computing,
        save_memory=False, keep_expression_matrix=True)
    lioness_obj = Lioness(panda_obj, computing=computing, save_dir=output_dir,
        save_fmt='npy', output=None, **lioness_options)

### Main body ###
dir_path = os.path.split(snakemake.output[0])[0]

# Prepare LIONESS options
input_options = {
    'expression': snakemake.input[0],
    'motif_prior': None,
    'ppi_prior': None,
    'output_dir': dir_path,
    'computing': 'gpu' if use_gpu else 'cpu',
    'lioness_options': {
        'start': 1,
        'end': int(snakemake.config['n_networks']),
        'ncores': snakemake.threads,
    },
}

# Run LIONESS
if use_gpu:
    gpu_devices = [cp.cuda.Device(i) for i in range(4)]

    with allocate_gpus(snakemake.params['gpu_manager'],
        snakemake.resources['gpus']) as gpu_ids:
        # Only one device is yielded here
        with gpu_devices[gpu_ids[0]]:
            create_lioness_networks(**input_options)
else:
    create_lioness_networks(**input_options)

# Convert the generated .npy file to .feather with the right index
genes = pd.read_csv(snakemake.input[0], sep='\t', header=None, 
    index_col=0).index
mi = pd.MultiIndex.from_product([genes] * 2, names=['gene1', 'gene2'])
npy_path = os.path.join(dir_path, 'lioness.npy')
data = np.load(npy_path)
df = pd.DataFrame(data, index=mi)
df.columns = [str(i + 1) for i in df.columns]
df.reset_index().to_feather(snakemake.output[0])

# Remove the superfluous .npy file
os.remove(npy_path)