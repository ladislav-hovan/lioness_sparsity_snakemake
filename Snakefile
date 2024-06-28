### Configuration file ###
configfile: 'config.yaml'
# configfile: os.path.join('tests', 'test_config.yaml')

### Imports and setup ###
import os
import pickle

from pathlib import Path

from lib.functions import get_most_recent_log_time

global_resources = workflow.global_resources

if 'gpus' in global_resources and global_resources['gpus'] > 0:
    from lib.gpu_manager import GpuManager, allocate_gpus

    GPU_GUARD = os.path.join('.snakemake', '.gpu_guard')
    GPU_MANAGER = os.path.join('.snakemake', '.gpu_manager.pkl')

    if (not os.path.exists(GPU_GUARD) or
        get_most_recent_log_time() > os.path.getmtime(GPU_GUARD)):
        print ('Running one-time GPU checks and setting up the manager')

        assert global_resources['gpus'] == len(config['gpu_ids']), (
            'The number of GPUs in resources must match the length of '
            'provided GPU IDs')

        gpu_manager = GpuManager(config['gpu_ids'])

        pickle.dump(gpu_manager, open(GPU_MANAGER, 'wb'))

        Path(GPU_GUARD).touch()
    else:
        gpu_manager = pickle.load(open(GPU_MANAGER, 'rb'))

    USE_GPU = True
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