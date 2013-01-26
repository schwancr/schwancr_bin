#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='traj_dir',default='./Trajectories_AA/Autocorrelations',help='Directory to find the autocorrelation data. [ ./Trajectories_AA/Autocorrelations ]')
parser.add_option('-o',dest='out_FN',default='./AutoFits.dat',help='Output filename [ ./AutoFits.dat ]')
parser.add_option('-p',dest='procs',default=1,type=int,help='Number of processes to run. Be careful, since this is also the number of tajectories open at one time, so you could have a memory issue. [ 1 ]')

options, args = parser.parse_args()
 
import numpy as np
from pyschwancr import dataIO, FitData
import multiprocessing as mp
import os, sys, re

def AnalyzeTraj( trajFN ):
	print "Working on %s " % trajFN
	traj = dataIO.readData( trajFN )
	X = np.arange( len( traj ) )
	timescales = []
	for i in xrange( traj.shape[1] ):
		timescales.append(1. / FitData.ExponFit( X, traj[:,i], LogSample = True )[0] )

	return np.array( timescales )

trajList = dataIO.getTrajList( options.traj_dir, RegEx = '^trj\d+\.npy' )

pool = mp.Pool( options.procs )
result = pool.map_async( AnalyzeTraj, trajList )
sol = np.vstack( result.get() )

np.savetxt( options.out_FN, sol )
print sol.shape
