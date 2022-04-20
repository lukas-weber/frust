#!/usr/bin/env python3

import numpy as np
import scipy.sparse.linalg as spsl
import os
import argparse
import json
from loadleveller import jobfile

import model
from ed import cluster_magnet
from ed import cavity_magnet


parser = argparse.ArgumentParser(
    description="ED code for the frust project. Automatically calculates the results for a chosen loadleveller jobfile. If there is more than one task in the jobfile, those tasks must only vary by temperature."
)

parser.add_argument("jobfile", type=str, help="Jobscript describing the model.")
parser.add_argument("-o", "--outfile", type=str, help=".json file to save results to")

args = parser.parse_args()


def summarize_tasks(job):
    summarized_tasks = []
    Ts = []
    for task, params in job.tasks.items():

        T = params.pop("T")
        found = False
        for previous in summarized_tasks:
            if params == previous[1]:
                previous[2].append(T)
                found = True
        if not found:
            summarized_tasks.append((task, params, [T]))

    return summarized_tasks


def solve_task(H, num_states=0):
    if num_states == 0:
        H = np.array(H.todense())
        E, psi = np.linalg.eigh(H)
        print("Full diagonalization done.")
    else:
        E, psi = spsl.eigsh(H, num_states, which="SA")
        Emax = E.max()
        psi = psi[:, E < E.max() - 1e-3]
        E = E[E < E.max() - 1e-3]
        print("Sparse diagonalization done ({}/{} states).".format(len(E), num_states))

    Egap = E[E > E.min() + 1e-8].min() - E.min()
    print(
        "{} energies: {:.3g}..{:.3g}, gap={:.3g}".format(len(E), E.min(), E.max(), Egap)
    )

    return E, psi


def to_list(a):
    if type(a) == np.ndarray:
        if len(a) > 1:
            return a.tolist()
        return a[0]
    return a


def get_model(mod):
    if mod.model == "cluster_magnet":
        return cluster_magnet.Model(mod)
    if mod.model == "cavity_magnet":
        return cavity_magnet.Model(mod)
    raise Exception('unsupported model "{}"'.format(mod.model))


jobdir = os.path.dirname(args.jobfile)
if jobdir != "":
    os.chdir(jobdir)

job = jobfile.JobFile("./" + os.path.basename(args.jobfile))

output = []

summarized_tasks = summarize_tasks(job)


for i, (taskname, params, Ts) in enumerate(summarized_tasks):
    print("task {}/{}...".format(i + 1, len(summarized_tasks)))
    mod = model.load(job, taskname, force_overwrite=True)

    ed_model = get_model(mod)

    H = ed_model.hamiltonian()
    print(
        "dimension = 2**{:d} = {}".format(
            int(np.log(H.shape[0]) / np.log(2)), H.shape[0]
        )
    )

    E, psi = solve_task(H, params.get("ed_num_states", 0))

    obs = ed_model.observables(params, np.array(Ts), E, psi)
    print("observables calculated")

    for i, T in enumerate(Ts):
        task_observables = {
            obsname: to_list(obsvalue[i]) for obsname, obsvalue in obs.items()
        }
        output.append({"parameters": {"T": T, **params}, "results": task_observables})

if args.outfile:
    with open(args.outfile, "w") as outfile:
        json.dump(output, outfile, indent=1)
else:
    print(json.dumps(output, indent=4))
