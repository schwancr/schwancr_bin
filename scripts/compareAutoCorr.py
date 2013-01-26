#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_FNs',action='append',help='Data files to compare to one another. Will use this regex to get the x-axis: .*_([\d.-]+). Note: This script will open x_FoldTimes.dat and x_UnfoldTimes.dat, where x is the argument passed in.')
parser.add_option('-o',dest='out_FN',help='Output file name (This is the root of all the plots made.')
options, args = parser.parse_args()
 
from numpy import *
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from pyschwancr import dataIO
import os, sys, re
from scipy.optimize import curve_fit

X_values = []
Titles = options.data_FNs 
DatList = []

for fn in Titles:
	DatList.append( dataIO.readData( fn ) )
	
	m = re.search( '.*_([-.\d]+)', fn )
	if m:
		X_values.append( float( m.group(1) ) )
	else:
		print "Filename contains no number !!!! (In the form .*_[-.\d]+)"
		exit()

print "Loaded Data."
# First make the individual plots:

figure()
for i in range( len( Titles ) ):
	plot( DatList[i], label = str( X_values[i] ) + '%' )
	
legend()
xscale('log')
xlabel('Time (frames)')
ylabel('Autocorrelation')
title('Autocorrelations for different topologies')
savefig( options.out_FN + 'AutoCorrs.pdf')
print "Saved autocorrelation plots to %sAutoCorrs.pdf" % options.out_FN

expFit = lambda x, p : exp( - p * x )

FitExps = []
for i in range( len( Titles ) ):
	sol = curve_fit( expFit, arange( len( DatList[i] ) ), DatList[i], p0 = [1] )
	print sol
	FitExps.append( sol[0] )

figure()
scatter( X_values, 1./array( FitExps ) )
xlabel(' Percent of native well depth ')
ylabel('Timescale (from exponential fit)')
title('Exponential timescales for various topologies')
savefig( options.out_FN + 'Timescales.pdf')
print "Saved timescales to %sTimescales.pdf" % options.out_FN


