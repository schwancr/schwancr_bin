#!/home/schwancr/Installed/epd/bin/python -u
 
from optparse import OptionParser
parser = OptionParser()

parser.add_option('-t',dest='traj_dir',default='./', help='Directory to look for trajectory data')
parser.add_option('-d',dest='metric',help='Distance metric to use in msmbuilder.DistanceMetric [ RMSD for lh5\'s or QNorm for npy\'s ]' )
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of processes to use in multiprocessing [ 1 ]' )
parser.add_option('-o',dest='output',help='Output filename [ METRIC_Drift.npy where METRIC is the input argument ]' )
parser.add_option('-u',dest='stride',type=int,default=1,help='Stride to look at the drift for. [ 1 ]' )
parser.add_option('--dt',dest='tau',type=int,help='Tau to use to calculate the drift.')
options, args = parser.parse_args()
 
from numpy import *
from pyschwancr import dataIO
from msmbuilder import Trajectory
from schwancrtools.metrics_Drift import get_epsilon_neighborhoods
import multiprocessing as mp
import os, sys, re
import pickle

dist = pickle.load( open( options.metric ) )


def AnalyzeTraj( trajFN ):
	# this function analyzes a trajectory to get the drift a given distance metric
	traj = Trajectory.LoadTrajectoryFile( trajFN )[::options.stride]

	print "Working on %s" % trajFN
	
	ptraj = dist.prepare_trajectory( traj )
 
        outAry = get_epsilon_neighborhoods( dist, ptraj, options.tau )

	return outAry

trajFNs = [ ( int( fn[3:-4] ), fn ) for fn in os.listdir( options.traj_dir ) if re.search( '^trj', fn ) ]
trajFNs.sort()
trajFNs = [ os.path.join( options.traj_dir, fn ) for (i, fn) in trajFNs ]
print "Loaded Trajectories."
#print trajFNs
print "Calculating the drifts."
if options.metric.lower() == 'rmsd':
	print "RMSD is already parallelized so only running one process..."
	pool = mp.Pool( 1 )
else:
	pool = mp.Pool( options.procs )

result = pool.map_async( AnalyzeTraj, trajFNs )
result.wait()
sol = result.get()

solVec = []

for thing in sol:
	solVec.extend( thing )

solVec = array( solVec )
if options.output:
	outFN = options.output
else:
	outFN = options.metric + '_drift.npy' 
save( outFN, solVec )
print "Saved output to %s." % outFN 
