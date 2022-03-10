#!/usr/bin/env python3
import numpy as np
import sys
sys.path.append('../scripts')
import model
import ed.cluster_magnet
from loadleveller import jobfile
import argparse
import shutil
import yaml
import os

parser = argparse.ArgumentParser(description='Test that compares the Hamiltonian of different basis formulations of the triangle square lattice')
parser.add_argument('mcbinary', type=str, help='loadleveller binary')
parser.add_argument('jobfile', type=str, help='Jobscript containing tasks that should yield the same Hamiltonian')

args = parser.parse_args()

jobname = os.path.basename(args.jobfile)
jobdir = 'test_equal_hamiltonian/'

os.makedirs(jobdir, exist_ok=True)
shutil.copy(args.jobfile, jobdir+jobname)
os.chdir(jobdir)

shutil.rmtree(jobname+'.data', ignore_errors=True)

jobconfig = {
    'num_cores': 1,
    'mc_binary': '../' + args.mcbinary,
    'mc_runtime': '24:00:00',
    'mc_checkpoint_time': '24:00:00',
}

with open('jobconfig.yml', 'w') as f:
    yaml.dump(jobconfig, f)

job = jobfile.JobFile('./'+jobname)

H = None
first = True

for task in job.tasks.keys():
    model_data = model.load(job, task)
    mod = ed.cluster_magnet.Model(model_data)

    Hnew = mod.hamiltonian()

    if not first:
        if not np.allclose(Hnew.todense(), H.todense()):
            idcs = np.argwhere(np.logical_not(np.isclose(H.todense(),Hnew.todense())))
            for i, j in zip(idcs[:,0], idcs[:,1]):
                print('{} {}: {} != {}'.format(i, j, H[i,j].real, Hnew[i,j].real))
            sys.exit(1)
     
    H = Hnew
    first = False

