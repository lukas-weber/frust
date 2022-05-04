#!/usr/bin/env python3

import model
import numpy as np
import argparse
import svgwrite
from loadleveller import jobfile

parser = argparse.ArgumentParser(
    description="Draw lattice from a cavity_magnet task into an svg file"
)

parser.add_argument("jobfile", type=str, help="Jobfile describing the model.")
parser.add_argument(
    "-t",
    "--taskname",
    type=str,
    default="task0001",
    help="Taskname of the lattice in question (default task0001).",
)
parser.add_argument("-o", "--outfile", type=str, help="file to save results to")

args = parser.parse_args()

job = jobfile.JobFile(args.jobfile)
lat = model.load(job, args.taskname)

size = 100

poss = np.array([s.pos for s in lat.sites])

bbox = np.min(poss[:, 0]), np.max(poss[:, 0]), np.min(poss[:, 1]), np.max(poss[:, 1])

scale = size / max(bbox[1] - bbox[0], bbox[3] - bbox[2])
ratio = (bbox[1] - bbox[0]) / (bbox[3] - bbox[2])
poss -= np.array([bbox[0], bbox[2]])
poss *= scale

max_mc = np.max(np.abs([b.mode_couplings for b in lat.bonds]))

mode_colors = ["red", "blue"]
d = svgwrite.Drawing(
    args.outfile, size=(size * min(1, ratio), size * min(1, 1 / ratio)), profile="full"
)

for p in poss:
    d.add(d.circle(p, r=0.05 * scale, fill="black"))

for b in lat.bonds:
    diff = poss[b.j] - poss[b.i]
    if np.linalg.norm(diff) > size / 2:
        continue
    d.add(
        d.line(
            start=poss[b.i],
            end=poss[b.j],
            stroke="black",
            stroke_width=0.1,
        )
    )
    for i, mc in enumerate(b.mode_couplings):
        pos = poss[b.i] + (0.5 + 0.4 * (i - (len(b.mode_couplings) - 1) / 2)) * diff
        d.add(d.circle(pos, r=0.1 * scale * abs(mc) / max_mc, fill=mode_colors[i]))


if args.outfile:
    d.save()
else:
    print(d.tostring())
