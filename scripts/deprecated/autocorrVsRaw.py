#!/usr/bin/env python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-a','--assignments',dest='AsFN',default='./Assignments.Fixed.h5',help='Assignments.h5 file from msmbuilder (You should probably use Assignments.Fixed.h5 from BuildMSM.py')
parser.add_option('-t','--trans-prob',dest='tFN',default='./tProb.mtx',help='Transition probability matrix for the MSM')
parser.add_option('--rd','--raw-data',dest='rawFN',help='Raw data for all the trajectories. Be sure this corresponds to ALL of the trajectories...')
parser.add_option('--md','--msm-data',dest='msmFN',help='State averages for all states in the Assignments file')
parser.add_option('--nt','--num-trajs',dest='nTraj',type=int, default=1,help='Number of raw trajectories to plot')
parser.add_option('-l','--lag-time',dest='lag',type=int,help='Lag time of the particular MSM')
parser.add_option('-p','--proj-info',dest='projFN',default='./ProjectInfo.h5',help='ProjectInfo.h5 generated by msmbuilder')
options, args = parser.parse_args()

from msmbuilder.MSMLib import PropagateModel as pm
from msmbuilder import Serializer
from schwancrtools import autocorrelate
import Correlation
from numpy import *
from pyschwancr.dataIO import readData as rd
from pyschwancr.dataIO import writeData as wd
from scipy.io import mmread
import os
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
import scipy.sparse


def msmAvg( T, N, start, avg ):

#	final, avgVsT = pm( T, N, start, ObservableVector=avg)
# Using my own propogator because msmbuilder's is super slow...

	trajList = zeros(( N+1 , len(start)))
	trajList[0] = start
	if scipy.sparse.issparse(T):
		T = T.tocsr()
	for i in range(1,N+1,1):
		print "\r %5d" % i,
		trajList[i] = trajList[i-1] * T
	
	avgVsT = dot( trajList, avg )
	return avgVsT

def rawObs( trajNum, lens, data ):
	
	trjData = data[ lens[:trajNum].sum() : lens[:trajNum].sum()+lens[trajNum] ]

	return trjData

def plotData( msmList, rawList, trajNums, lag ):
	global options

	if os.path.exists('AutoCorrPDFs'):
		print "Writing data to AutoCorrPDFs. If this is bad, stop now! Data has been saved so you can plot it all manually"
	else:
		print "Writing data to ./AutoCorrPDFs"
		os.mkdir('AutoCorrPDFs')

	for i, traj in enumerate(trajNums):
		figure()
		xRaw = arange(len(rawList[i]))
		xMSM = arange(len(msmList[i])) * lag
		if lag != 1:
			xMSM[0] = 1
		plot(xRaw,rawList[i],label='Raw Data')
		plot(xMSM,msmList[i],label='MSM Average (lag = %d)'%lag)
		xlabel('Time (frames)')
		ylabel('Autocorrelation')
		title('Autocorrelation Function %s' % options.rawFN.split('/')[-1].split('.')[-2] )
		legend()
		xlim([1,xMSM.max() ] )
		ylim([ -0.4, 1.2])
		xscale('log')
		savefig( os.path.join(os.path.abspath('.'),'AutoCorrPDFs','%s_traj%d.pdf'% ( '.'.join( options.rawFN.split('/')[-1].split('.')[:-1] ),traj ) ) )
		
		

def main():

	print "Reading data...",
	print " %s."% options.rawFN,
	raw = rd( options.rawFN )
	print " %s." % options.msmFN,
	msm = rd( options.msmFN )
	msm = msm[:,1] # Use the center column (This is how my data is formated.)
	print " %s." % options.AsFN,
	As = Serializer.LoadData( options.AsFN )
	print " %s." % options.projFN,
	P = Serializer.LoadFromHDF( options.projFN )
	print " %s." % options.tFN
	T = mmread( options.tFN )
	

	lens = [ ( P['TrajLengths'][i], i ) for i in range( P['TrajLengths'].shape[0] ) ]
	lens = sorted(lens,reverse=True)
	
	trajList = [ b for (a,b) in lens[:options.nTraj] ]

	rawFmt = ones( As.shape ) * -1
	tempSum = 0

	print raw.shape

	for i in range( len( P['TrajLengths'] ) ):
		rawFmt[i][ : P['TrajLengths'][i] ] = raw[ tempSum : tempSum + P['TrajLengths'][i] ]
		tempSum += P['TrajLengths'][i]

	# Generate the trajectory data:
	trajData = []
	msmData = []
	N = int( P['TrajLengths'].max() / options.lag ) + 1
	
	for traj in trajList:
		print "\nCalculating trajectory %d" % traj,
		print "Raw data.",
		trajData.append( array(rawObs( traj, P['TrajLengths'], raw ) ) )
		print "MSM Average."
		start = zeros( len(msm) )
		start[ As[ traj ][0] ] += 1
		Cor, Traj, ObsTraj = Correlation.RawMSMCorrelation( T, rawFmt, As, Steps = N, StartingState = As[ traj ][0] ) 
		msmData.append( Cor )

	# Generate the autocorrelations:
	print "Calculating the convolutions ... "
	trajCorr = [ autocorrelate.fft_autocorrelate( thing ) for thing in trajData ]
	msmCorr = [ autocorrelate.fft_autocorrelate( thing ) for thing in msmData ]
	msmCorr = [ thing for thing in msmData ]

	for i in range(1000):
		if not os.path.exists('Autocorr%d'%i):
			outDir = 'Autocorr%d'%i
			os.mkdir( outDir )
			break

	if not outDir:
		print "You have a lot of Autocorr directories..."
		print "Writing output to ./ which could overwrite other data!"
		outDir = './'

	print "Writing the autocorrelation data"
	for index, traj in enumerate( trajList ):
		wd( [ 'autocorr','traj%d'%traj,'msm','lag%d'%options.lag ], msmCorr[index], dir=outDir )
		wd( [ 'autocorr','traj%d'%traj,'raw'], trajCorr[index], dir=outDir )

	plotData( msmCorr, trajCorr, trajList, options.lag )
if __name__ == '__main__':
	main()
