import numpy as np
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import lattice
import argparse
import ed_hamiltonian
from loadleveller import jobfile
import os

parser = argparse.ArgumentParser(description='ED code for the frust project. Automatically calculates the results for a chosen loadleveller jobfile. If there is more than one task in the jobfile, those tasks must only vary by temperature.')

parser.add_argument('jobfile', type=str, help='Jobscript describing the model.')
parser.add_argument('-o', '--outfile', type=str, help='.npz file to save results to')

args = parser.parse_args()

def extract_jobdata(job):
    tasknames = []
    Ts = []
    for task, params in job.tasks.items():
        if params['ed_compare']:
            tasknames.append(task)
            Ts.append(params['T'])
    if len(tasknames) == 0:
        raise Exception('jobfile did not contain any tasks with ed_compare=True')

    return tasknames[0], np.array(Ts)

jobdir = os.path.dirname(args.jobfile)
if jobdir != '':
    os.chdir(jobdir)

job = jobfile.JobFile('./'+os.path.basename(args.jobfile))

taskname, Ts = extract_jobdata(job)

num_states = job.tasks[taskname].get('ed_num_states', 0)
lat = lattice.load(job, taskname, force_overwrite=True)

N = sum(s.nspinhalfs for s in lat.sites)
basis = np.array([
	       [0,	 0 ,  isq_,     0 , -isq_,     0,	       0, 0],
	       [0,	 0 ,     0,  -isq_,     0, isq_ ,	       0, 0],
	       [0, -2*isq6_, isq6_,     0 , isq6_,     0,	       0, 0],
	       [0,	 0 ,	0 , isq6_ ,     0, isq6_, -2*isq6_      , 0]])

Nsite = len(lat.sites)
Nstates = basis.shape[0]

H = sps.dok_matrix(Nstates**Nsites, Nstates**Nsites))



