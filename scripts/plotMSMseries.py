#!/usr/bin/env python
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from msmbuilder import arglib
from msmbuilder import msm_analysis
import numpy as np
from scipy.io import mmread

parser = arglib.ArgumentParser()

parser.add_argument( 'data', help='Data filename which contains the averagevalue of some piece of data for each state in the MSM.')
parser.add_argument( 'tProb', help='tProb.mtx from msmbuilder', default='tProb.mtx' )
parser.add_argument( 'output', help='Output filename to plot the time sereis')
parser.add_argument( 'lagtime', type=float, help='Lag time for your MSM.')
parser.add_argument( 'units', help='Units of your lagtime.')
parser.add_argument( 'num_modes', default=20,type=int, help='Number of eigenmodes to consider.')
parser.add_argument( 'num_points', default=1E4,type=int,help='Number of points to plot in the' )
parser.add_argument( 'starting_pops',default='all',help='Starting populations. Pass "all" to use the uniform distribution')
parser.add_argument( 'equilibrium_pops', default='Populations.dat' )
parser.add_argument( 'y_lbl', default='Ensemble Average', help='Label for y-axis')
parser.add_argument( 'normalize',default=False, action='store_true',help='Pass this flag to normalize to the equilibrium value.' )
args = parser.parse_args()

matplotlib.rcParams['font.size']=22

equil_pops = np.loadtxt( args.equilibrium_pops )
if args.starting_pops != 'all':
    starting_pops = np.loadtxt( args.starting_pops ).flatten()

    if starting_pops.shape[0] == equil_pops.shape[0]:
        starting_pops /= starting_pops.sum() # Normalize just in case.
    else:
        starting_inds = starting_pops.copy()
        starting_pops = np.zeros(equil_pops.shape)
        starting_pops[starting_inds.astype(int)] = 1 / float(len(starting_inds))
else:
    starting_pops = None

try: data = np.loadtxt( args.data )[:,1]
except: data = np.loadtxt(args.data).flatten()
t = mmread( args.tProb )

sol = msm_analysis.calc_expectation_timeseries( t, data, lagtime = args.lagtime, n_modes = args.num_modes, timepoints = args.num_points, init_pop=starting_pops )

N = len(sol)

equil_val = np.dot( equil_pops, data )

h=0.18
axes( ( h, h, 1-2*h, 1-2*h ) )
xi = np.arange( N ) * args.lagtime

if args.normalize:
    plot( xi, sol / equil_val, label='$y(t)$', lw=3 )
    hlines( 1, xlim()[0], xlim()[1], label='$y(\infty)$', linewidth=3, linestyle='--' )
else:
    plot( xi, sol, label='$y(t)$', lw=3 )
    hlines( equil_val, xlim()[0], xlim()[1], label='$y(\infty)$', linewidth=3, linestyle='--' )

xlabel( 'Time (%s)'% args.units, fontsize=27 )
ylabel( args.y_lbl, fontsize=27 )

old_lim = ylim()
yticks( yticks()[0][1::2] )

if sol.mean() < equil_val:
    legend(loc=4)
else:
    legend(loc=1)
ylim(old_lim)
savefig( args.output )

out_ary = np.array( zip( xi, sol ) )
np.save( args.output[:-4] +'.npy', out_ary )
