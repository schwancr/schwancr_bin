#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-i','--msm-file',dest='msmFN',help='Flat text file ( 3 columns: i,j,k i=State number, j=average quantity, k=standard deviation')
parser.add_option('-d','--raw-data-file',dest='rawFN',help='The quantity for all conformations in the raw data, numpy save format (.npy)')
parser.add_option('-w','--weights',dest='weightFN',help='Cluster populations in flat text format')
parser.add_option('-o','--output-file',dest='outFN',help='Output pdf with the comparison plot')
parser.add_option('-b','--bin-number',dest='bins',type=int, default=100, help='Number of bins to use in histograms. Pick a divisor of the total number of contacts you distributing over')
parser.add_option('--xl',dest='x_lbl',default='Order Parameter', help='X-Axis label for the plot')

options, args = parser.parse_args()

import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from numpy import *

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

subplot(211)
hist(msmData[:,1],bins=N,weights=pops,range=(0., max([1,raw.max()]) ))
ylabel('Frequency')

#yticks( [100000, 200000, 300000, 400000, 500000, 600000 ] )
title('MSM %s' % options.msmFN[:-4] )

x=0.6
y=ylim()[1]*0.8

text(x,y, "Num. States = %d" % len(pops) )

# Next plot the raw data:


subplot(212)
hist(raw,bins=N,range=(0.,max([1,raw.max()])) )
ylabel('Frequency')
#yticks( [100000, 200000, 300000, 400000, 500000, 600000 ] )
title('Raw Data')
xlabel( options.x_lbl )
newXlim = xlim()
subplot(211)
xlim( newXlim )
savefig(options.outFN)
