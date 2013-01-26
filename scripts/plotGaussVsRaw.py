#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-i','--msm-file',dest='msmFN',help='Flat text file ( 3 columns: i,j,k i=State number, j=average quantity, k=standard deviation')
parser.add_option('-d','--raw-data-file',dest='rawFN',help='The quantity for all conformations in the raw data, numpy save format (.npy)')
parser.add_option('-w','--weights',dest='weightFN',default='./clusterPops.dat',help='Cluster populations in flat text format')
parser.add_option('-o','--output-file',dest='outFN',help='Output pdf with the comparison plot')
parser.add_option('-b','--bin-number',dest='bins',type=int, default=100, help='Number of bins to use in histograms. Pick a divisor of the total number of contacts you distributing over')
parser.add_option('--xl','--x-label',dest='x_lbl',default='Data [ A.U ]',help='X-axis label')
options, args = parser.parse_args()

import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from numpy import *
from pyschwancr import msmTools

figure()
suptitle(options.outFN[:-4],fontsize=14,weight='bold')
# First plot the msm data:

pops = loadtxt(options.weightFN)
msmData = loadtxt(options.msmFN)
raw = load(options.rawFN)

# Number of bins:
N = options.bins

if pops.shape != msmData[:,1].shape:
	print "Weights must be the same length as data! pops: %d data: %d"% (pops.shape[0], msmData[:,1].shape[0])

# Plot the raw data:

subplot(212)
histSol = hist(raw,bins=N,range=(0.,max([1,raw.max()])) )
ylabel('Frequency')
#yticks( [100000, 200000, 300000, 400000, 500000, 600000 ] )
title('Raw Data')
xlabel(options.x_lbl)

newXlim = xlim()
newYlim = ylim()
barWidth = histSol[1][1] - histSol[1][0]
subplot(211)
msmXY = msmTools.distributeStateAvgs( msmData, pops * barWidth, bins=N, xLimit=newXlim )# [0, max([1, raw.max()])] )
bar( msmXY[:,0], msmXY[:,1], align='center', width = barWidth)
ylabel('Frequency')
#print "Raw: %f" % sum( barWidth * histSol[0] )
#print "MSM: %f" % sum( barWidth * msmXY[:,1] )
#print "Actual: %f" % sum( barWidth * pops )
#print "Width: %f" % barWidth
#yticks( [100000, 200000, 300000, 400000, 500000, 600000 ] )
title('MSM %s' % options.msmFN[:-4] )

xlim( newXlim )
ylim( newYlim )
x=xlim()[1]*0.6
y=ylim()[1]*0.8
text(x,y, "Num. States = %d" % len(pops) )

savefig(options.outFN)
