import json
import subprocess
import numpy as np
from loadleveller import jobfile

class Site:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class Bond:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class Lattice:
    def __init__(self, latticedef):
        lattice = json.loads(latticedef)

        self.Lx = lattice['Lx']
        self.Ly = lattice['Ly']
        self.energy_offset = lattice['energy_offset']
        self.uc_spin_count = lattice['uc_spin_count']
        self.bonds = [Bond(**b) for b in lattice['bonds']]
        self.sites = [Site(**s) for s in lattice['sites']]

def load(job, taskname):
    mc_binary = job.jobconfig['mc_binary']

    job_input_filename = job.write_job_input_file(force_overwrite=False)
    res = subprocess.run([mc_binary, 'lattice', job_input_filename, taskname], stdout=subprocess.PIPE)

    return Lattice(res.stdout.decode('utf-8'))
