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


class Mode:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class Model:
    def _load_lattice(self, lattice):
        self.Lx = lattice["Lx"]
        self.Ly = lattice["Ly"]
        self.uc_site_count = lattice["uc_site_count"]
        self.bonds = [Bond(**b) for b in lattice["bonds"]]
        self.sites = [Site(**s) for s in lattice["sites"]]

    def _load_cluster_magnet(self, model):
        self._load_lattice(model)

    def _load_cavity_magnet(self, model):
        self._load_lattice(model)
        self.modes = [Mode(**m) for m in model["modes"]]

    def __init__(self, latticedef):
        model = json.loads(latticedef)
        self.model = model["model"]

        if self.model == "cluster_magnet":
            self._load_cluster_magnet(model)
        elif self.model == "cavity_magnet":
            self._load_cavity_magnet(model)
        else:
            raise Exception("unsupported model {}".format(self.model))


def load(job, taskname, force_overwrite=False):
    mc_binary = job.jobconfig["mc_binary"]

    job_input_filename = job.write_job_input_file(force_overwrite=force_overwrite)
    res = subprocess.run(
        [mc_binary, "lattice", job_input_filename, taskname], stdout=subprocess.PIPE
    )

    return Model(res.stdout.decode("utf-8"))
