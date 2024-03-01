#!/usr/bin/env python

### Imports ###
import cupy as cp
import pandas as pd

from netZooPy.panda import Panda
from netZooPy.lioness import Lioness

from snakemake_gpu_manager import allocate_gpus

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
        save_fmt='npy', **lioness_options)

### Main body ###
gpu_devices = [cp.cuda.Device(i) for i in range(4)]

with allocate_gpus(snakemake.params['gpu_manager'], 
    snakemake.resources['gpus']) as gpu_ids:
    # Only one device is yielded here
    with gpu_devices[gpu_ids[0]]:
        create_lioness_networks(snakemake.input[0], snakemake.input[1], 
            snakemake.input[2], '.', 'gpu', 
            lioness_options={
                'start': 1, 
                'end': int(snakemake.config['n_networks']),
                'export_filename': snakemake.output[0],
            })