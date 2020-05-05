import numpy as np
import matplotlib.pyplot as plt
import loadleveller.mcextract as mce
import argparse
import scipy.stats
import sys
import re

parser = argparse.ArgumentParser(description='Compare ed results (generated by ed.py) to mc results (generated by isingsse). Both from the same jobfile, obviously.')

parser.add_argument('edfile', type=str, help='.npz file generated by ed.py')
parser.add_argument('mcfile', type=str, help='.results.yaml file from isingsse')
parser.add_argument('-o', '--observables', type=str, help='Regex to filter observables to be compared.')
parser.add_argument('-c', '--cli', action='store_true', help='cli-only mode. No plots are shown.')
parser.add_argument('-p', '--confidence', type=float, default=1e-3, help='confidence interval for the χ² test to fail.')

args = parser.parse_args()

mc = mce.MCArchive(args.mcfile)
ed = np.load(args.edfile)

cond = {'ed_compare' : True}
Ts = np.array(mc.get_parameter('T',filter=cond))

mc_obs = set(mc.observables.keys())
ed_obs = set(ed.keys())

obs_names = mc_obs & ed_obs
if args.observables:
    obs_names = list(filter(lambda x: re.fullmatch(args.observables, x) != None, obs_names))

print('Common observables:')
print(*obs_names, sep='\n')
print('\nConfidence level: p = {:.3g}\n'.format(args.confidence))

if not args.cli:
    sign = mc.get_observable('Sign', filter=cond)
    plt.errorbar(Ts, sign.mean, sign.error, fmt='.')
    plt.show()


failure = 0

def compare(obs_name, mean, error, edobs):
    global failure
    missingvalue = np.isnan(mean)
    error[missingvalue] = 100
    smallerror = error < 1e-10 # roundoff error regime

    mask = np.logical_not(np.logical_or(smallerror, missingvalue))

    maskmean = mean[mask]
    maskerror = error[mask]


    chisq = np.sum((maskmean-edobs[mask])**2/maskerror**2)
    p = scipy.stats.chi2.sf(chisq,len(maskmean))

    print('{} χ²/dof = {:.2g}, p = {:.3g}'.format(obs_name, chisq/len(maskmean), p))


    if np.any(smallerror):
        print('WARNING: some values have errors smaller than machine precision. Ergodicity?')
    elif p < args.confidence:
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
    obs = mc.get_observable(obs_name, filter=cond)
    edobs = ed[obs_name]

    if len(obs.mean.shape) > 1:
        for i in range(obs.mean.shape[1]):
            compare('{}[{}]'.format(obs_name, i), obs.mean[:,i], obs.error[:,i], edobs[:,i])
    else:
        compare(obs_name, obs.mean, obs.error, edobs.flatten())
        
sys.exit(failure)