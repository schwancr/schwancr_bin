#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]' )
parser.add_option('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments.h5 from msmbuilder [ ./Assignments.h5.Fixed ]' )
parser.add_option('-s',dest='states',action='append',type=int,help='States to check (can give multiple \'-s\' inputs)' )
parser.add_option('-n',dest='num_states',default=10,type=int, help='Number of states to sample from the msm')
parser.add_option('--check-all',dest='check_all',action='store_true',default=False, help='Use this flag to sample all states in the msm.')
parser.add_option('-P',dest='procs',default=1,type=int,help='Number of processes to use')

options, args = parser.parse_args()

from pyschwancr import FitData
from msmbuilder import Project, DistanceMetric
from numpy import *
import multiprocessing as mp
import os
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *

Proj = Project.Project.LoadFromHDF( options.proj_FN )
Ass = Project.Serializer.LoadData( options.ass_FN )
RMSD = DistanceMetric.RMSDMetric()

def CheckStatePWrmsd( stateInd ):
	## This function checks how many gaussians to fit the state PW rmsd's
	which = array( where( Ass == stateInd ) ).T
	N = which.shape[0]

	if N < 20:
		print "Too few conformations (%d)... Skipping state %d" % (N, stateInd )
		return -1

	print "Checking state %d with it's %d conformations" % ( stateInd, N )
	confs = Proj.GetConformations( which )
	prepXYZ = RMSD.PrepareData( confs['XYZList'] )

	rmsds = zeros( (N,N) )
	for i in range( N ):
		rmsds[i,:] = RMSD.GetFastMultiDistance( prepXYZ, prepXYZ, i )

	rmsds_1d = rmsds[ where( eye( N ) == 0 ) ].flatten() # Ignore the diagonal zeros

	a,b,c = hist( rmsds_1d, bins=50 )
	xDist = array( [ ( b[i] + b[i+1] ) / 2. for i in range( len( a ) ) ] )
	yDist = array( a )
	
	# Now I've prepared the data, I need to fit it with multiple gaussians...
	Fits = []
	for i in range(1,100):
		Fits.append( FitData.GaussFit( xDist, yDist, N = i ) )
		if BadFit( Fits[-1], xDist.min(), xDist.max() ):
			return i-1

	return 0
	
def BadFit( params, xMin, xMax ):
	""" 
	Returns true if the fit is bad and false if it is good
	"""

	paramsList = params.reshape( ( params.shape[0] / 3, 3 ) )

	means = paramsList[:,1] 
	stdevs = abs( paramsList[:,2] )

	# Check if any means are outside the range of X	
	if (means < xMin).any():
		return True
	
	if (means > xMax).any():
		return True

	# Check if any of the gaussians are overlapping

	mins = means - stdevs
	maxs = means + stdevs
	#print array( zip( mins, maxs ) )
	minLTmax = zeros( (len(mins),len(mins)) )
	for i in range( len(mins) ):
		for j in range( len(mins) ):
			minLTmax[i,j] = mins[i] < maxs[j]

	# This part is tricky, but think about it and it's true.
	# minLTmax is a 2d array containing 1 if mins[i] < max[j] (row i, col j)
	# then minLTmax.T is actually the negation of maxs[i] < mins[j] (row i, col j) (excluding equals)
	# I don't really care about equals since they are not likely to happen, and if they do, then they will be lumped in with not overlapping

	testMtx = minLTmax * minLTmax.T

	if testMtx[ where( eye( len(mins) ) == 0 ) ].any():
		return True

	return False

def main():
	
	# Need to construct the state list to check. There are three ways to construct it, sooooo I need a hierarchy.
	# I will say the hierarchy is ( left is most prefereable ) --check-all, -s, -n
	
	if options.check_all:
		stateList = arange( Ass.max() + 1 )
	elif options.states:
		stateList = array( options.states )
	else:
		if options.num_states >= ( Ass.max() + 1 ):
			print "Too many states requested ( %d ). There are only %d states in this MSM. Proceeding by sampling all states."
			stateList = arange( Ass.max() + 1 )
		else:
			stateList = random.permutation( arange( Ass.max() + 1 ) )[: options.num_states ]

	
	print "Constructed the list of states."
	print "Will look at these states:"
	print stateList
	
	pool = mp.Pool( options.procs )
	result = pool.map_async( CheckStatePWrmsd, stateList )
	result.wait()
	bestFits = result.get() 
	bestFits = array( bestFits )
	outAry = zip( stateList, bestFits )
	outAry.sort()
	outAry = array( outAry )

	for i in range(10000):
		if os.path.exists( 'BestFits%d.dat' % i ):
			continue
		else:	
			savetxt( 'BestFits%d.dat' % i, outAry, "%d" )
			break
	figure()
	bar( range(1, outAry[:,1].max()+1), bincount( outAry[:,1])[1:] )
	xlabel('Number of gaussians giving the "Best Fit"' )
	ylabel('Frequency')
	text( 0.8*xlim()[1], 0.7 * ylim()[1] ," N = %d " % outAry.shape[0] )
	title('Sum of gaussian fit to Pairwise RMSDs in %d states' % outAry.shape[0] )
	for i in range(100000):
		if os.path.exists( 'BestFits%d.dat' % i ):
			continue
		else:	
			savefig( 'BestFits%d.pdf' % i )
			break


if __name__ == '__main__':
	main()



