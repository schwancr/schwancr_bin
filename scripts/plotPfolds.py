#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-m',dest='msm_FN',help='FCommittors from TPT')
parser.add_option('-r',dest='raw_FN',help='Raw Pfolds in the format: StateXXX\nFolded = x Unfolded = y Neither = z Total = x+y+z')
parser.add_option('-o',dest='out_FN',help='Output to save plot')

options, args = parser.parse_args()

# This script will read in data about Pfolds and plot the forward committors on the y-axis and the Pfolds calculated (along with their standard deviations on the x-axis
# The two plots will be:
#  1) A plot of all data points
#  2) A plot of the avg plus StdDev

import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from numpy import *
from pyschwancr import dataIO
import re

msm = dataIO.readData(options.msm_FN)
rawFN = open( options.raw_FN, 'r' )

rawList = rawFN.read().split('State')
rawList.pop(0) # There is an initial '' since It splits at the first characters.
rawData = {}

for state in rawList:
	stateList = state.split('\n')
	stateDat = []
	for conf in stateList[1:]:
		m = re.search("Folded\s*=\s*(\d+)\s*Unfolded\s*=\s*(\d+)",conf)
		if m:
			N_unfolded = int( m.group( 2 ) )
			N_folded = int( m.group( 1 ) )
			stateDat.append( N_folded / float( N_folded + N_unfolded ) )
	rawData[ int( stateList[0] ) ] = array( stateDat )

rawDatName = '.'.join( options.raw_FN.split('/')[-1].split('.')[:-1] )
# Now make the first plot: All data points:
x = []
y = []
for key in rawData.keys():
	x.extend( rawData[key] )
	y.extend( [ msm[key] ] * len( rawData[key] ) )
figure()
scatter( x, y )
plot( [0,1],[0,1],'r-')
xlabel( rawDatName )
ylabel( 'Forward Committors from MSM' )
title( 'MSM vs Raw Pfolds' )
xlim([-.01,1])
ylim([0,1])
savefig( options.out_FN )


# Now make the first plot: All data points:
yShort = []
avg = []
std = []
for key in rawData.keys():
	yShort.append( msm[ key ] )
	avg.append( rawData[ key ].mean() )
	std.append( rawData[ key ].std() )

figure()
errorbar( avg,yShort, xerr=std, fmt="o")
plot( [0,1],[0,1],'r-')
xlabel( rawDatName + ' Averages' )
ylabel( 'Forward Committors from MSM' )
title( 'MSM vs Raw Pfolds' )
xlim([0,1])
ylim([0,1])
savefig( options.out_FN[:-4] + '_avg.pdf' )







