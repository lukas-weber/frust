#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import filecmp
import yaml
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
        'mc_binary': '../'+mc_binary,
        'mc_runtime': '24:00:00',
        'mc_checkpoint_time': '24:00:00',
    }

    with open('jobconfig.yml', 'w') as f:
        yaml.dump(jobconfig, f)
    

args = parser.parse_args()

jobname = os.path.basename(args.testjob)
jobdir = 'testjobs/'

os.makedirs(jobdir, exist_ok=True)
shutil.copy(args.testjob, jobdir+jobname)
os.chdir(jobdir)

generate_jobconfig(
    num_cores = multiprocessing.cpu_count()//2,
    mc_binary = args.mcbinary
)

if not args.seeded: 
    subprocess.run(['loadl', 'run', jobname, '-r'], check=True)
    if args.generate:
        # TODO
        pass
    else:
        subprocess.run(['python3', '../../scripts/ed_compare.py', args.result_file, jobname + '.results.json', '--cli'], check=True)
else:
    seeded_job = args.testjob + '_seeded'
    subprocess.run(['../../test/gen_seeded_job.sh', args.testjob, jobname + '_seeded'], check=True)
    subprocess.run(['loadl', 'run', jobname + '_seeded', '-sr'], check=True)
    if args.generate:
        shutil.copy(jobname+'_seeded.results.json', args.generate)
    else:
        if not filecmp.cmp(jobname+'_seeded.results.json', args.result_file, shallow=False):
            raise Exception('new and old seeded result file do not match! {} != {}'.format(mc_result, args.seeded_result))
