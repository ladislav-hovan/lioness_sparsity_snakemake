#!/usr/bin/env python

### Imports ###
import pandas as pd

from pathlib import Path

from lib.functions import load_file

### Functions ###
def calculate_degrees(
    lioness_feather: Path,
    ind_feather: Path,
    outd_feather: Path,
) -> None:
    """
    Generates indegree and outdegree files from LIONESS networks saved
    in feather format.

    Parameters
    ----------
    lioness_feather : Path
        The input feather file with the LIONESS networks
    ind_feather : Path
        The output feather file with the calculated indegrees
    outd_feather : Path
        The output feather file with the calculated outdegrees
    """

    df = load_file(lioness_feather)

    ind_df = df.groupby('gene').sum()
    ind_df.reset_index().to_feather(ind_feather)

    outd_df = df.groupby('tf').sum()
    outd_df.reset_index().to_feather(outd_feather)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-l', '--lioness', dest='lioness_feather',
        help='file to load the LIONESS networks from', metavar='FEATHER')
    parser.add_argument('-i', '--indegree', dest='ind_feather',
        help='file to save the indegrees into', metavar='FEATHER')
    parser.add_argument('-o', '--outdegree', dest='outd_feather',
        help='file to save the indegrees into', metavar='FEATHER')

    parser.parse_args()

    calculate_degrees(
        parser.lioness_feather,
        parser.ind_feather,
        parser.outd_feather,
    )