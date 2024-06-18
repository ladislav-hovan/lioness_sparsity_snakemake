#!/usr/bin/env python

### Imports ###
import os

import numpy as np
import pandas as pd

from netZooPy.lioness import Lioness
from netZooPy.panda import Panda
from pathlib import Path
from typing import Optional

### Functions ###
def create_lioness_networks(
    expression: pd.DataFrame,
    motif_prior: pd.DataFrame,
    ppi_prior: pd.DataFrame,
    output_dir: str,
    computing: str,
    lioness_options: dict,
) -> None:
    """


    Parameters
    ----------
    expression : pd.DataFrame
        _description_
    motif_prior : pd.DataFrame
        _description_
    ppi_prior : pd.DataFrame
        _description_
    output_dir : str
        _description_
    computing : str
        _description_
    lioness_options : dict
        _description_
    """
    
    panda_obj = Panda(expression, motif_prior, ppi_prior, computing=computing, 
        save_memory=False, keep_expression_matrix=True)
    _ = Lioness(panda_obj, computing=computing, save_dir=output_dir, 
        save_fmt='npy', **lioness_options)


def calculate_lioness_networks(
    expression: Path,
    motif_prior: Path,
    ppi_prior: Path,
    output: Path,
    n_networks: int,
    gpu_id: Optional[int] = None,
    threads: int = 1,
) -> None:
    """
    

    Parameters
    ----------
    expression : Path
        _description_
    motif_prior : Path
        _description_
    ppi_prior : Path
        _description_
    output : Path
        _description_
    n_networks : int
        _description_
    gpu_id : Optional[int], optional
        _description_, by default None
    threads : int, optional
        _description_, by default 1
    """

    dir_path = os.path.split(output)[0]

    # Prepare LIONESS options
    input_options = {
        'expression': expression,
        'motif_prior': motif_prior,
        'ppi_prior': ppi_prior,
        'output_dir': dir_path,
        'computing': 'gpu' if gpu_id is not None else 'cpu',
        'lioness_options': {
            'start': 1,
            'end': n_networks,
            'ncores': threads,
        },
    }

    # Run LIONESS
    if gpu_id is not None:
        from cupy.cuda import Device

        with Device(gpu_id):
            create_lioness_networks(**input_options)
    else:
        create_lioness_networks(**input_options)

    # Convert the generated .npy file to .feather with the right index
    prior = pd.read_csv(motif_prior, sep='\t', header=None,
        names=['tf', 'gene', 'edge'])
    tfs = sorted(prior['tf'].unique())
    genes = sorted(prior['gene'].unique())
    mi = pd.MultiIndex.from_product([genes, tfs], names=['gene', 'tf']
        ).swaplevel()
    npy_path = os.path.join(dir_path, 'lioness.npy')
    data = np.load(npy_path)
    df = pd.DataFrame(data, index=mi)
    df.columns = [str(i + 1) for i in df.columns]
    df.sort_index().reset_index().to_feather(output)

    # Remove the superfluous .npy file
    os.remove(npy_path)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-mp', '--motif-prior', dest='motif_p',
        help='file to load the motif prior from', metavar='FILE')
    parser.add_argument('-pp', '--ppi-prior', dest='ppi_p',
        help='file to load the PPI prior from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the result into', metavar='FILE')
    parser.add_argument('-nn', '--n-networks', dest='n_net', type=int,
        help='number of networks to generate', metavar='INT')
    parser.add_argument('-gid', '--gpu-id', dest='gpu_id', type=int,
        help='ID of the GPU to use', metavar='ID')
    parser.add_argument('-t', '--threads', dest='threads', type=int,
        help='number of threads to use', metavar='INT', default=1)

    parser.parse_args()

    calculate_lioness_networks(
        parser.expression,
        parser.motif_p,
        parser.ppi_p,
        parser.output,
        parser.n_net,
        parser.gpu_id,
        parser.threads,
    )