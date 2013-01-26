#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option( '-t',dest='traj_dir',help=' Directory to look for trajectories in the form trj##.lh5 [ ./Trajectories ]' )
parser.add_option( '-o',dest='out_FN',default='RawRgs.npy',help='Output filename [ ./RawRgs.npy ] ')
parser.add_option( '-p',dest='procs',default=1,type=int, help='Number of processes to run [ 1 ] ')
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Trajectory
from pyschwancr import dataIO, msmTools
import os, sys, re
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import *
import multiprocessing as mp
from time import time

def analyzeTraj( trjFN ):

	print "Working on trajectory %s" % trjFN
	trj = Trajectory.Trajectory.LoadFromLHDF( trjFN )

	return msmTools.calcRg( trj )
	
trajFNs = [ ( int( thing[3:-4] ), thing ) for thing in os.listdir( options.traj_dir ) if re.search( '^trj\d+\.lh5', thing ) ]
trajFNs.sort()
trajFNs = [ os.path.join( options.traj_dir, b ) for (a,b) in trajFNs ]
print "Loaded Data."
maxLength = 0
Rgs = []

pool = mp.Pool( options.procs )
result = pool.map_async( analyzeTraj, trajFNs)
result.wait()

Rgs = concatenate( result.get() )

save( options.out_FN, Rgs )
