### Configuration file ###
configfile: 'config.yaml'
# configfile: os.path.join('tests', 'test_config.yaml')

### Imports and setup ###
import os

global_resources = workflow.global_resources

if 'gpus' in global_resources and global_resources['gpus'] > 0:
    assert global_resources['gpus'] == len(config['gpu_ids']), (
        'The number of GPUs in resources must match the length of provided '
        'GPU ids')

    USE_GPU = True

    import cupy as cp

    gpu_devices = [cp.cuda.Device(i) for i in 
        range(cp.cuda.runtime.getDeviceCount())]

    from lib.gpu_manager import GpuManager, allocate_gpus

    gpu_manager = GpuManager(config['gpu_ids'])
else:
    USE_GPU = False

    gpu_manager = None
    allocate_gpus = None

### Configuration variables ###
input_dir = config['input_dir']

### Input and output variables ###
include: 'Snakefile_variables'

### Rules ###
rule all:
    input:
        SUMMARY_CORR_PLOTS,
        COEXPR_ERROR_PLOTS,
        expand(
            os.path.join('coexpression_networks', '{transform}', 
                '{method}', '{sparsity}', '{repeat}', 'lioness.feather'),
            transform=config['transformations'],
            method=config['sparsifying_methods'],
            sparsity=config['sparsity_levels'],
            repeat=range(config['n_repeats']),
        ),
        expand(
            os.path.join('coexpression_correlations', '{transform}', 
                '{method}', '{sparsity}', 'pearson.npy'),
            transform=config['transformations'],
            method=config['sparsifying_methods'],
            sparsity=config['sparsity_levels'],
        ),
    default_target:
        True

# TODO: Convert scripts to Python functions, call from here as run block
# TODO: Make sure the scripts are command line callable as well
# TODO: Effectively implement one main function that's called with arguments

include: 'Snakefile_expression'

include: 'Snakefile_networks'

include: 'Snakefile_correlations'

include: 'Snakefile_plotting'