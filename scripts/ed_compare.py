#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import loadleveller.mcextract as mce
import argparse
import scipy.stats
import sys
import re
import json

parser = argparse.ArgumentParser(description='Compare ed results (generated by ed.py) to mc results (generated by isingsse). Both from the same jobfile, obviously.')

parser.add_argument('edfile', type=str, help='.json file generated by ed.py')
parser.add_argument('mcfile', type=str, help='.results.json file from Monte Carlo')
parser.add_argument('-o', '--observables', type=str, help='Regex to filter observables to be compared.')
parser.add_argument('-c', '--cli', action='store_true', help='cli-only mode. No plots are shown.')
parser.add_argument('-p', '--confidence', type=float, default=1e-4, help='confidence interval for the χ² test to fail.')

args = parser.parse_args()

mc = mce.MCArchive(args.mcfile)
with open(args.edfile, 'r') as edfile:
    ed = json.load(edfile)

mc_obs = set(mc.observables.keys())

obs_names = mc_obs
if args.observables:
    obs_names = list(filter(lambda x: re.fullmatch(args.observables, x) != None, obs_names))

print('Common observables:')
print(*obs_names, sep='\n')
print('\nConfidence level: p = {:.3g}\n'.format(args.confidence))


def match_ed_tasks(mc, ed):
    if mc.num_tasks != len(ed):
        raise ValueError('mc and ed data do not have matching number of tasks')
    matched_tasks = [] 
    for i in range(mc.num_tasks):
        params = {name: values[i] for name, values in mc.parameters.items()}
        match = [task for task in ed if task['parameters'] == params]
        for task in ed:
            if task['parameters'] == params:
                matched_tasks.append(task)
                break
        else:
            raise ValueError('ed data is missing task present in mc')
    return matched_tasks    

failure = 0

def compare(obs_name, Ts, mean, error, edobs):
    global failure

    
    missingvalue = np.isnan(mean)
    error[missingvalue] = 100
    smallerror = error < 1e-10 # roundoff error regime

    mask = np.logical_not(np.logical_or(smallerror, missingvalue))

    maskmean = mean[mask]
    maskerror = error[mask]

    if np.any(smallerror):
        print('WARNING: some values have errors smaller than machine precision. Ergodicity?')
    if len(maskmean) > min(len(mean),3):
        chisq = np.sum((maskmean-edobs[mask])**2/maskerror**2)
        p = scipy.stats.chi2.sf(chisq,len(maskmean))

        print('{} χ²/dof = {:.2g}, p = {:.3g}'.format(obs_name, chisq/len(maskmean), p))
        if p < args.confidence:
            print('\033[91mWARNING: χ² outside of confidence interval.\033[0m')
            failure = 1

    if not args.cli:
        plt.subplot(211)
        plt.ylabel(obs_name)
        plt.plot(Ts, edobs, label='ED')
        plt.errorbar(Ts[mask], maskmean, maskerror, fmt='.', label='MC')
        if np.any(smallerror):
            plt.plot(Ts[smallerror], mean[smallerror], '.', label='MC (small error)')
        plt.subplot(212)
        plt.grid(True)
        plt.errorbar(Ts[mask], maskmean-edobs[mask], maskerror, fmt='.', label='ED-MC')
        if np.any(smallerror):
            plt.plot(Ts[smallerror], mean[smallerror]-edobs[smallerror], '.', label='ED-MC (small error)')
        plt.xlabel('T')
        plt.show()

for obs_name in sorted(obs_names):
    obs = mc.get_observable(obs_name)
    Ts = mc.get_parameter('T')
    ed_idxs = [i for i, task in enumerate(ed) if obs_name in task['results'].keys()]
    if len(ed_idxs) == 0:
        continue
    ed_values = np.array([ed[i]['results'][obs_name] for i in ed_idxs])

    if len(obs.mean.shape) > 1:
        for i in range(obs.mean.shape[1]):
            compare('{}[{}]'.format(obs_name, i), Ts[ed_idxs], obs.mean[ed_idxs, i], obs.error[ed_idxs, i], ed_values[:, i])
    else:
        compare(obs_name, Ts[ed_idxs], obs.mean[ed_idxs], obs.error[ed_idxs], ed_values)
        
sys.exit(failure)
