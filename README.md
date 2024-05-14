# LIONESS sparsity testing
This Snakemake pipeline measures the sensitivity of LIONESS to varying 
levels of sparsity.


## Table of Contents
- [LIONESS sparsity testing](#lioness-sparsity-testing)
  - [Table of Contents](#table-of-contents)
  - [General Information](#general-information)
  - [Features](#features)
  - [Setup](#setup)
  - [Usage](#usage)
  - [Project Status](#project-status)
  - [Room for Improvement](#room-for-improvement)
  - [Acknowledgements](#acknowledgements)
  - [Contact](#contact)
  - [License](#license)


## General Information
This pipeline tests and summarises how the output of LIONESS changes 
when the gene expression data is modified to varying levels of sparsity. 
It is intended to automate this process fully, starting from input 
expression file, prior motif network, prior PPI (protein-protein 
interaction) network, and the settings provided in the config.yaml file.

The entire pipeline is implemented using Snakemake. It will utilise
GPUs to compute networks when available, but can also use CPUs only.


## Features
The features already available are:
- Cleanup of expression data and prior networks
- Generation of sparsified expression
- Generation of LIONESS networks
- Summary Pearson and Spearman correlations comparing to the unmodified
  baseline:
  - edge weights
  - indegrees
  - expression
  - coexpression
- Average error in coexpression across sparsity levels


## Setup
The requirements are provided in a `requirements.txt` file.


## Usage
Running a Snakemake pipeline is straightforward:

``` bash
snakemake --cores=10 --resources gpus=2
```

It is assumed that all the input is present and that the settings in the
`config.yaml` file are correct. Relevant settings are:
- `input_dir`: the directory containing the input files
- `expression_file`: name of the gene expression file
- `motif_file`: name of the motif prior file
- `ppi_file`: name of the PPI prior file
- `n_networks`: the number of networks to be generated from the list
of samples, must not exceed the total number of samples in the 
expression file (all samples are always used for background)
- `n_repeats`: the number of realisations of any given sparsity level
- `sparsity_levels`: the levels of sparsity to generate expression and
networks for, in percentage points
- `gpu_ids`: the IDs of GPUs to use, length should match the number of
GPUs specified by `resources`


## Project Status
The project is: _in progress_.


## Room for Improvement
Room for improvement:
- Change the passing of arguments to rules
- Finish defining variables for input and output
- Finish refactoring shared functionality

To do:
- Stratification of results by gene expression levels


## Acknowledgements
Many thanks to the members of the 
[Kuijjer group](https://www.kuijjerlab.org/) 
at NCMM for their feedback and support.

This README is based on a template made by 
[@flynerdpl](https://www.flynerd.pl/).


## Contact
Created by Ladislav Hovan (ladislav.hovan@ncmm.uio.no).
Feel free to contact me!


## License
This project is open source and available under the 
[GNU General Public License v3](LICENSE).