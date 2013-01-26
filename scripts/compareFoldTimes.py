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
 
FoldTimeList = []
UnfoldTimeList = []
X_values = []
Titles = options.data_FNs 
maxVal = 0
for fn in Titles:
	FoldTimeList.append( dataIO.readData( fn + '_FoldTimes.dat' ) )
	UnfoldTimeList.append( dataIO.readData( fn + '_UnfoldTimes.dat' ) )

	if FoldTimeList[-1].max() > maxVal:
		maxVal = FoldTimeList[-1].max()

	if UnfoldTimeList[-1].max() > maxVal:
		maxVal = UnfoldTimeList[-1].max()

	m = re.search( '.*_([-.\d]+)', fn )
	if m:
		X_values.append( float( m.group(1) ) )
	else:
		print "Filename contains no number !!!! (In the form .*_[-.\d]+)"
		exit()

print "Loaded Data."
# First make the individual plots:



for i in range( len( Titles ) ):
	figure()
	subplot(211)
	hist( FoldTimeList[i], bins=100, color = 'blue', label='Fold Times', range = [0, maxVal ] )
	ylabel( 'Frequency' )
	vlines( FoldTimeList[i].mean(), 0, ylim()[1], color = 'black' )
	text( 0.3 * maxVal, 0.7 * ylim()[1], 'Mean = %.2f\nN = %d' % ( FoldTimeList[i].mean(), len( FoldTimeList[i] ) ) )
	legend()

	subplot(212)
	hist( UnfoldTimeList[i], bins=100, color = 'red', label='Unfold Times', range = [0, maxVal ] )
	vlines( UnfoldTimeList[i].mean(), 0, ylim()[1], color = 'black' )
	text( 0.3 * maxVal, 0.7 * ylim()[1], 'Mean = %.2f\nN = %d' % ( UnfoldTimeList[i].mean(), len( UnfoldTimeList[i] ) ) )
	ylabel( 'Frequency' )
	xlabel( 'Time (frames)' )
	legend()

	subplot(211)
	suptitle( Titles[i] )
	savefig( '%s_UF_dist.pdf' % Titles[i] )
	print "Saved distribution for %s to %s" % ( Titles[i], '%s_UF_dist.pdf' % Titles[i] )

FoldAvgs = array( [ thing.mean() for thing in FoldTimeList ] )
UnfoldAvgs = array( [ thing.mean() for thing in UnfoldTimeList ] )

N_folds = array( [ len( thing ) for thing in FoldTimeList ] )
N_unfolds = array( [ len( thing ) for thing in UnfoldTimeList ] )

figure()
scatter( X_values, FoldAvgs, color = 'blue', label = 'Fold Time' )
scatter( X_values, UnfoldAvgs, color = 'red', marker = 's', label = 'Unfold Time' )
xlabel(u"\u03F5( non-native ) / \u03F5( native ) * 100")
ylabel('Time (frames)')
title( 'Folding/Unfolding times vs. Relative Well Depths' )
legend(loc=2)
savefig( options.out_FN + 'TimeComp.pdf' )
print "Saved Time comparison to %s" % ( options.out_FN + 'TimeComp.pdf' )

figure()
scatter( X_values, N_folds, color = 'blue', label = 'Fold Time' )
scatter( X_values, N_unfolds, color = 'red', marker = 's', label = 'Unfold Time' )
xlabel(u"\u03F5( non-native ) / \u03F5( native ) * 100")
ylabel('Number of transitions')
title( 'Total Folding/Unfolding Transitions vs. Relative Well Depths' )
legend(loc=1)
savefig( options.out_FN + 'NumTransComp.pdf' )
print "Saved Transition number comparison to %s" % ( options.out_FN + 'NumTransComp.pdf' )
