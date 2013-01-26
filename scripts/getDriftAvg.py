#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python -u
 
from optparse import OptionParser
parser = OptionParser()

parser.add_option('-t',dest='traj_dir',default='./', help='Directory to look for trajectory data')
parser.add_option('-d',dest='metric',help='Distance metric to use in msmbuilder.DistanceMetric [ RMSD for lh5\'s or QNorm for npy\'s ]' )
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of processes to use in multiprocessing [ 1 ]' )
parser.add_option('-o',dest='output',help='Output filename [ METRIC_Drift.npy where METRIC is the input argument ]' )
parser.add_option('-u',dest='stride',type=int,default=1,help='Stride to look at the drift for. [ 1 ]' )
options, args = parser.parse_args()
 
from numpy import *
from pyschwancr import dataIO
from msmbuilder import Trajectory, DistanceMetric
import multiprocessing as mp
import os, sys, re


#MetricDict = { 'qnorm': DistanceMetric.QNorm, 'rmsd': DistanceMetric.RMSD }
MetricDict = DistanceMetric.MetricDict
if options.metric:
	try: 
		Dist = MetricDict[ options.metric.lower() ]
	except:
		print "%s not a viable metric label. Please input one of: %s" % ( options.metric, str( MetricDict.keys() )[1:-1] )
		exit()
else:
	print "Please input a metric to use of: %s" % str( MetricDict.keys() )[1:-1] 

def AnalyzeTraj( trajFN ):
	# this function analyzes a trajectory to get the drift a given distance metric
	traj = dataIO.readData( trajFN )[::options.stride]
	if options.metric.lower() == 'qnorm':
		traj = traj.astype(uint8)

	print "Working on %s" % trajFN
	drifts = []
	for i in range(len( traj )-1 ): # Start at one and end one before the last because I can't average over two values.
		drifts.append( Dist.GetDistance( traj[i], traj[i+1] ) )
		#print drifts[-1]

	drifts = array( drifts )
	midDrifts = ( drifts[1:] + drifts[:-1] ) / 2.
	
	outAry = concatenate( ( [ drifts[0] ], midDrifts, [ drifts[-1] ]) )
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
pool.close()
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
