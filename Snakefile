### Configuration file ###
configfile: 'config.yaml'
# configfile: os.path.join('tests', 'test_config.yaml')

### Imports and setup ###
import os

from pathlib import Path

global_resources = workflow.global_resources

GPU_GUARD = os.path.join('.snakemake', '.gpu_guard')

if 'gpus' in global_resources and global_resources['gpus'] > 0:
    if not os.path.exists(GPU_GUARD) or (os.path.getmtime(GPU_GUARD) <
        os.path.getmtime(os.path.join('.snakemake', 'log'))):
        assert global_resources['gpus'] == len(config['gpu_ids']), (
            'The number of GPUs in resources must match the length of '
            'provided GPU IDs')
        Path(GPU_GUARD).touch()

    USE_GPU = True

    from lib.gpu_manager import GpuManager, allocate_gpus

    gpu_manager = GpuManager(config['gpu_ids'])
else:
    USE_GPU = False

### Configuration variables ###
input_dir = config['input_dir']

### Input and output variables ###
include: os.path.join('workflow', 'variables.smk')

### Rules ###
rule all:
    input:
        SUMMARY_CORR_PLOTS,
        COEXPR_ERROR_PLOTS,
    default_target:
        True

# TODO: Convert scripts to Python functions, call from here as run block
# TODO: Make sure the scripts are command line callable as well
# TODO: Effectively implement one main function that's called with arguments

include: os.path.join('workflow', 'expression.smk')

include: os.path.join('workflow', 'networks.smk')

include: os.path.join('workflow', 'correlations.smk')

include: os.path.join('workflow', 'plotting.smk')