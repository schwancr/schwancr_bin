#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from msmbuilder import arglib
matplotlib.rcParams['font.size']=22

parser = arglib.ArgumentParser()
parser.add_argument('eigenvalues', help='Eigenvalues (excluding equilibrium)')
parser.add_argument('y_lim', nargs=2, type=float, help='Y-limits')
parser.add_argument('lagtime', type=float, default=1, help='Lag time of MSM')
parser.add_argument('units', default='lagtime', help='Units of lag time')
parser.add_argument('x_label', default='', help='Label for x-axis')
parser.add_argument('plot_red', default=[], nargs='*', type=int, help='plot these timescales in red and dashed. (0-indexed)')
parser.add_argument('output', default='timescales.pdf', help='output filename (pdf)')
parser.add_argument('no_y', default=False, action='store_true', help='Pass this flag to turn of y-labels')

args = parser.parse_args()

vals = np.load(args.eigenvalues)

times = - args.lagtime / np.log(vals)

red_inds = np.array(args.plot_red)
print red_inds
blue_inds = np.arange(len(times))
blue_inds = blue_inds[np.where(1 - np.in1d(blue_inds, red_inds))]

fig = plt.figure(figsize=(6.9,6))
ax = fig.add_subplot(132)

if len(red_inds) >= 1:
    ax.hlines(times[red_inds], 0, 1, lw=3, color='red', linestyle='dashed')
    
ax.hlines(times[blue_inds], 0, 1, lw=3, color='blue')
ax.set_xlim([0,1])
ax.set_ylim(args.y_lim)
ax.set_yscale('log')
if not args.no_y:
    ax.set_ylabel('Implied Timescales (%s)' % args.units)
else:
    ax.set_yticks([])

#ax.set_xlim([0,1])
ax.set_xticks([])
ax.set_xlabel(args.x_label)

fig.savefig(args.output)

