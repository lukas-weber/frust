import numpy as np
import loadleveller.mcextract as mce
import matplotlib.pyplot as plt
import argparse
import subprocess
import time

parser = argparse.ArgumentParser(
    description="Benchmark a job by evaluating a Monte Carlo efficiency measure: error at constant time"
)
parser.add_argument("jobfile", type=str, help="job to run")

args = parser.parse_args()

subprocess.run(
    ["loadl", "run", "--restart", "--single", args.jobfile],
    check=True,
    capture_output=True,
)

mc = mce.MCArchive(args.jobfile + ".results.json")

stime = mc.get_observable("_ll_sweep_time")
mtime = mc.get_observable("_ll_measurement_time")
for obsname in ["SignEnergy", "Sign", "SignStagXStagYMag2"]:
    obs = mc.get_observable(obsname)

    stime_error = (
        obs.error
        / np.abs(obs.mean)
        * (
            obs.rebinning_bin_length
            * obs.rebinning_bin_count
            * (stime.mean + mtime.mean)
        )
        ** 0.5
    )

    print("{:20s}: {:.3g}".format(obsname, stime_error.mean() * 60))
