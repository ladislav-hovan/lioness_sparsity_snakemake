# Code by MuhammedHasan on GitHub
# https://github.com/snakemake/snakemake/issues/281

import random

from contextlib import contextmanager
from pytools.persistent_dict import PersistentDict

class GpuManager:
    """
    Manage GPUs keep track of which GPUs are free and which are 
    allocated
    """

    def __init__(self, gpus):
        self._gpus = PersistentDict('.gpu_manager', container_dir='.snakemake')
        self.free_gpus = set(gpus)
        self.allocated_gpus = set()

    @property
    def free_gpus(self):
        return self._gpus.fetch('free_gpus')

    @free_gpus.setter
    def free_gpus(self, value):
        self._gpus.store('free_gpus', value)

    @property
    def allocated_gpus(self):
        return self._gpus.fetch('allocated_gpus')

    @allocated_gpus.setter
    def allocated_gpus(self, value):
        self._gpus.store('allocated_gpus', value)

    def acquire(self, num_gpus):
        if len(self.free_gpus) < num_gpus:
            raise ValueError(
                f"Cannot allocate {num_gpus} GPUs,"
                f" only {len(self.free_gpus)} available"
            )
        gpus = random.sample(tuple(self.free_gpus), num_gpus)

        self.free_gpus = self.free_gpus.difference(gpus)
        self.allocated_gpus = self.allocated_gpus.union(gpus)

        return gpus

    def release(self, gpus):
        self.free_gpus = self.free_gpus.union(gpus)
        self.allocated_gpus = self.allocated_gpus.difference(gpus)

@contextmanager
def allocate_gpus(gpu_manager, num_gpus):
    """
    Allocate several GPUs from the GPU manager and release them when 
    done.

    Args:
        gpu_manager (GpuManager): the GPU manager
        num_gpus (int): the number of GPUs to allocate
    """

    gpus = gpu_manager.acquire(num_gpus)
    try:
        yield gpus
    finally:
        gpu_manager.release(gpus)