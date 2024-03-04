# LIONESS sparsity testing
This Snakemake pipeline measures the sensitivity of LIONESS to varying levels
of sparsity.


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

The entire pipeline is implemented using Snakemake. In its current
iteration, it depends on GPUs to compute LIONESS networks.


## Features
The features already available are:
- Cleanup of expression data and prior networks
- Generation of sparsified expression
- Generation of LIONESS networks
- Summary Peason and Spearman correlations of both expression and 
LIONESS networks with the unmodified baseline


## Setup
The requirements are provided in a `requirements.txt` file.


## Usage
Running a Snakemake pipeline is straightforward:

``` bash
snakemake --cores=10 --resources gpus=2
```

Currently the pipeline requires at least 1 GPU. It is assumed that all 
the input is present (in the `input/` directory) and that the settings 
in the `config.yaml` file are correct. Relevant settings are:
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
- Define variables for input and output

To do:
- Implement undersampling for sparsified expression generation
- Make the pipeline compatible with CPU-only setups
- Add coexpression and plots that show the effect of sparsity on it


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