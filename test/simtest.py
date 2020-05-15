#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import yaml
import shutil
import os

parser = argparse.ArgumentParser(description='Test script that runs a couple of simulations and checks them against ed automatically')
parser.add_argument('mcbinary', type=str, help='loadleveller binary')
parser.add_argument('testjob', type=str, help='testjob to simulate')

args = parser.parse_args()

num_cores = multiprocessing.cpu_count()//2

jobname = os.path.basename(args.testjob)
jobdir = 'testjobs/'

os.makedirs(jobdir, exist_ok=True)
shutil.copy(args.testjob, jobdir+jobname)
os.chdir(jobdir)

jobconfig = {
    'num_cores': num_cores,
    'mc_binary': '../' + args.mcbinary,
    'mc_runtime': '24:00:00',
    'mc_checkpoint_time': '24:00:00',
}

with open('jobconfig.yml', 'w') as f:
    yaml.dump(jobconfig, f)


res = subprocess.run(['../../scripts/quicktest.sh', jobname, '--cli'], check=True)
if res.returncode != 0:
    raise Exception('simulation had nonzero return code')

