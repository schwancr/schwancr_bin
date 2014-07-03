#!/usr/bin/env python

from msmbuilder import io
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from argparse import ArgumentParser

matplotlib.rcParams['font.size'] = 22

parser = ArgumentParser()
parser.add_argument('-f', dest='amps', help='amplitudes from GetObsAmps.py')
parser.add_argument('-o', dest='output')
parser.add_argument('-l', dest='lagtime', type=float,
    help='lagtime of MSM')
parser.add_argument('-u', dest='units', default='frames',
    help='units of the lagtime')

args = parser.parse_args()

f = io.loadh(args.amps)

vals = f['evals']
amps = f['amplitudes']

ind = np.where(vals > 0)
vals = vals[ind]
amps = amps[ind]

times = - args.lagtime / np.log(vals)

y, x=np.histogram(times, weights=amps, bins=10 ** np.linspace(1, 6, 20))

x = np.vstack([x[:-1], x[1:]]).T.flatten()
y = np.vstack([y, y]).T.flatten() / np.abs(y).max()

figure()
axes((0.18, 0.18, 0.72, 0.72))

fill_between(x, y, lw=2)

xlabel('Timescale (%s)' % args.units)
ylabel('Relative Amplitude')
ylim(-1.2, 1.2)
xscale('log')

hlines([-1, 1], xlim()[0], xlim()[1], color='red', lw=2)

savefig(args.output)
