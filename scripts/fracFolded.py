#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder')
parser.add_option('-t',dest='trans_FN',default='./tProb.mtx',help='Transition matrix from msmbuilder')
parser.add_option('-r',dest='raw_FN',help='Raw data to calculate Frac Folded Raw')
parser.add_option('-m',dest='msm_FN',help='MSM data to calculate Frac Folded MSM')
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 file from msmbuilder')
parser.add_option('--fc',type=float,dest='fCut',help='Folded Cutoff')
parser.add_option('--low-is-folded',dest="low_is_folded",default=False,action="store_true",help='If the folded state is low in the rxn coord then pass this flag')
parser.add_option('-l',type=int,dest='lag',help='Lagtime of the particular msm')
parser.add_option('--title',dest='title',help='Title for plot')
options, args = parser.parse_args()

from numpy import *
from msmbuilder import Serializer, Project
from pyschwancr import msmTools, dataIO
from scipy.io import mmread
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *

def main():
	# First load in all the data:
	print "Loading Data..."
	ass = Serializer.LoadData(options.ass_FN)
	tProb = mmread(options.trans_FN)
	Proj = Project.Project.LoadFromHDF(options.proj_FN)
	rawAry = dataIO.readData(options.raw_FN)
	msmAry = dataIO.readData(options.msm_FN)

	if options.low_is_folded:
		testFcn = lambda x: x <= options.fCut
	else:
		testFcn = lambda x: x >= options.fCut
	print "Calculating the raw folded over time..."
	# Now break up the raw data by trajectory:
	sum = 0
	rawTrajs = []
	for i in range( len( Proj['TrajLengths'] ) ):
		rawTrajs.append( rawAry[ sum : sum + Proj['TrajLengths'][i] ] )
		sum += Proj['TrajLengths'][i]

	# Now calculate the steps it takes to fold for each trajectory:
	time2fold_raw = []
	for trj in rawTrajs:
		count = 0
		for frame in trj:
			if testFcn(frame):
				break
			count += 1
		time2fold_raw.append( count )	
	time2fold_raw = array( time2fold_raw )
	rawName = dataIO.writeData( [ "RAW_FractionFolded" ], time2fold_raw )
	# Now calculate the msm fraction folded (using msmTools.calcFracFold)
	x0 = zeros( tProb.shape[0] ) # This start vector is based on the first frame of all the trajectories
	for i in range( ass.shape[0] ):
		x0[ ass[i,0] ] += 1.
	x0 /= float(ass.shape[0])
	print "Defining the folded state and calculating the MSM folded over time"
	# Need to define the folded state. Will use a cutoff, but if there are no states below/above the cutoff then pick the min/max value
	Fstates = []
	for index,stateAvg in enumerate(msmAry[:,1]):
		if testFcn( stateAvg ):
			Fstates.append( index )
	Fstates = array( Fstates )

	if not Fstates.any():
		if options.low_is_folded:
			Fstates = array([ where( msmAry[:,1] == msmAry[:,1].min() ) ] )
		else:
			Fstates = array([ where( msmAry[:,1] == msmAry[:,1].max() ) ] )
	N = time2fold_raw.max() / options.lag
	time2fold_msm = msmTools.calcFracFold( Fstates, tProb, x0, N = N )
	datName = dataIO.writeData( [ "MSM_FractionFolded", str(options.lag) ], time2fold_msm )
	print "Saved data to %s" % datName
	# Now plot everything
	print "Making plot ..."
	hist( time2fold_raw, bins=100, histtype='step',label="Raw Data",cumulative=True,normed=True)
	
	plot( arange( N ) * options.lag , time2fold_msm, label="MSM")
	hlines( 1.0, xmin=0, xmax=time2fold_raw.max(),color='red' )
	xlim([0,time2fold_raw.max()])
	ylim([0,1.25])
	legend()
	xlabel( 'Time (frames)' )
	ylabel( 'Fraction Folded' )
	if options.title:
		title( 'Fraction folded over time (%s)' % options.title)
	else:
		title( 'Fraction folded over time' )
	text( 0.75 * time2fold_raw.max(), 0.2, "N = %d" % ass.shape[0] )
	savefig( "FracFolded_%s_rawVsMsm.pdf" % '.'.join( options.raw_FN.split('/')[-1].split('.')[:-1] ) )
	print "Plot saved to %s" % ("FracFolded_%s_rawVsMsm.pdf" % '.'.join( options.raw_FN.split('/')[-1].split('.')[:-1] ) )	

if __name__ == '__main__':
	main()


