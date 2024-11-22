#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import filecmp
import yaml
import json
import shutil
import os

parser = argparse.ArgumentParser(description='Test script that runs a Monte Carlo simulation and checks the results against ED (statistically) or a previous “seeded” reference result file (exactly). Seeded simulations are run with reduced sweeps')
parser.add_argument('mcbinary', type=str, help='loadleveller binary')
parser.add_argument('testjob', type=str, help='testjob to simulate')
parser.add_argument('--seeded', action='store_true', help='run in seeded mode')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--result-file', type=str, help='ED or seeded result file')
group.add_argument('--generate', type=str, help='generate ED or seeded result file')


def generate_jobconfig(num_cores, mc_binary):
    jobconfig = {
        'num_cores': num_cores,
        'mc_binary': mc_binary,
        'mc_runtime': '24:00:00',
        'mc_checkpoint_time': '24:00:00',
    }

    with open('jobconfig.yml', 'w') as f:
        yaml.dump(jobconfig, f)
    

args = parser.parse_args()

mc_binary = os.path.abspath(args.mcbinary)

jobname = os.path.basename(args.testjob)
jobdir = 'testjobs/'

os.makedirs(jobdir, exist_ok=True)
shutil.copy(args.testjob, jobdir+jobname)
os.chdir(jobdir)

generate_jobconfig(
    num_cores = multiprocessing.cpu_count()//2,
    mc_binary = mc_binary
)

if not args.seeded: 
    if args.generate:
        subprocess.run(['../../scripts/ed.py', jobname, '-o', args.generate], check=True) 
    else:
        subprocess.run(['loadl', 'run', jobname, '-r'], check=True)
        subprocess.run(['../../scripts/ed_compare.py', args.result_file, jobname + '.results.json', '--cli'], check=True)
else:
    seeded_job = args.testjob + '_seeded'
    subprocess.run(['../../test/gen_seeded_job.sh', args.testjob, jobname + '_seeded'], check=True)
    subprocess.run(['loadl', 'run', jobname + '_seeded', '-sr'], check=True)
    result_file = jobname+'_seeded.results.json'
    if args.generate:
        shutil.copy(result_file, args.generate)
    else:

        def load_without_profiling(filename):
            with open(filename, 'r') as f:
                r = json.load(f)
            return json.dumps([{'parameters': t['parameters'],
                'results': {name: result for name, result in t['results'].items() if not name.startswith('_ll_')},
                'task': t['task']} for t in r], sort_keys=True)

        if load_without_profiling(result_file) != load_without_profiling(args.result_file):
            raise Exception('new and old seeded result file do not match! {} != {}'.format(result_file, args.seeded_result))
