#!/usr/bin/env python3

import lattice
import numpy as np
import argparse
import svgwrite
from loadleveller import jobfile

parser = argparse.ArgumentParser(description='Draw lattice from a isingsse task into an svg file')

parser.add_argument('jobfile', type=str, help='Jobfile describing the model.')
parser.add_argument('-t', '--taskname', type=str, default='task0001', help='Taskname of the lattice in question (default task0001).')
parser.add_argument('-o', '--outfile', type=str, help='file to save results to')

args = parser.parse_args()

job = jobfile.JobFile(args.jobfile)
lat = lattice.load(job, args.taskname)

d = svgwrite.Drawing(args.outfile, size=(100,100), profile='full')

poss = []

maxJ = max([max(b.J) for b in lat.bonds])


for s in lat.sites:
    idx = np.arange(s.nspinhalfs)
    
    angle = 2*np.pi/s.nspinhalfs * idx+0.5
    r = s.nspinhalfs > 1

    pos = 10*np.array(s.pos)[:,None] - 2*r*np.array([np.cos(angle),np.sin(angle)])
    poss.append(pos)
    #d.add(d.circle([10*s.pos[0],10*s.pos[1]], r = 1, stroke='red'))

    for i in range(s.nspinhalfs):
        d.add(d.circle(pos[:,i].tolist(), r = 0.5, fill='black'))

        if len(s.Jin) != 0:
            idx = 0
            for j in range(s.nspinhalfs):
                if i < j:
                    d.add(d.line(start=pos[:,i], end=pos[:,j], stroke='black', stroke_width=0.1, opacity=s.Jin[idx]/maxJ))
                    idx += 1

for b in lat.bonds:
    nsi = lat.sites[b.i].nspinhalfs
    nsj = lat.sites[b.j].nspinhalfs
    for i in range(nsi): 
       for j in range(nsj):
           len = np.sqrt(np.sum((poss[b.i][:,i]-poss[b.j][:,j])**2))
           if len < 30:
               d.add(d.line(start=poss[b.i][:,i], end=poss[b.j][:,j], stroke='black', stroke_width = 0.1, opacity=b.J[i*nsj+j]/maxJ))

if args.outfile:
    d.save()
else:
    print(d.tostring())
