#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()

parser.add_option('-t',dest='traj_dir',default='./', help='Directory to look for trajectory data')
parser.add_option('-o',dest='output',default='rmsd_drift.npy',help='Output filename [ rmsd_drift.npy ]' )
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Trajectory, DistanceMetric
import multiprocessing as mp
import os, sys, re


def AnalyzeTraj( trajFN ):
	# this function analyzes a trajectory to get the drift a given distance metric
	traj = Trajectory.Trajectory.LoadFromLHDF( trajFN )['XYZList']
	print "Working on %s" % trajFN
	drifts = []
	for i in range(len( traj )-1 ): # Start at one and end one before the last because I can't average over two values.
		drifts.append( DistanceMetric.RMSD.GetDistance( traj[i], traj[i+1] ) )

	drifts = array( drifts )
	midDrifts = ( drifts[1:] + drifts[:-1] ) / 2.
	
	outAry = concatenate( ( [ drifts[0] ], midDrifts, [ drifts[1] ]) )
	return outAry

trajFNs = [ ( int( fn[3:-4] ), fn ) for fn in os.listdir( options.traj_dir ) if re.search( '^trj', fn ) ]
trajFNs.sort()
trajFNs = [ os.path.join( options.traj_dir, fn ) for (i, fn) in trajFNs ]
print "Loaded Trajectories."
#print trajFNs
print "Calculating the drifts."

sol = []
for fn in trajFNs:
	sol.append( AnalyzeTraj( fn ) )

solVec = []

for thing in sol:
	solVec.extend( thing )

solVec = array( solVec )
outFN = options.output
save( outFN, solVec )
print "Saved output to %s." % outFN 
