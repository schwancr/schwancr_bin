#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t',dest='tProb',default='./tProb.mtx',help='Transition probability matrix [ ./tProb.mtx ]')
parser.add_option('-p',dest='pops',default='./Populations.dat',help='Populations.dat from msmbuilder [ ./Populations.dat ]')
parser.add_option('-c',dest='cut',type=float,help='Cutoff to use to color states')
parser.add_option('-d',dest='data_FN',help='Data to use to define folded and unfolded states')
parser.add_option('--low-is-folded',dest='low_is_folded',default=False,action='store_true',help='Pass this flag if a low value of your order parameter indicates a folded state (e.g. RMSD)')
parser.add_option('-w',dest='write_dir',default='./',help='Directory to save output to [ ./ ]')
options, args = parser.parse_args()
 
from numpy import *
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import *
from pyschwancr import dataIO, msmTools
import os, sys, re
from scipy.io import mmread

T = mmread( options.tProb )
pops = dataIO.readData( options.pops )
data = dataIO.readData( options.data_FN )[:,1]
cut = options.cut 
print "Loaded Data."

if options.low_is_folded:
   isF = data < cut
else:
   isF = data > cut


degs = msmTools.getDegrees( T )
np.savetxt(os.path.join( options.write_dir,'Degrees.dat' ), degs)
Gs = - log( pops )
Gs = Gs - Gs.max()
plot( degs[ where( 1 - isF ) ], Gs[ where( 1 - isF ) ], 'rs', label='Unfolded States' )
plot( degs[ where( isF ) ], Gs[ where( isF ) ], 'b.', label='Folded States' )
text( 5, ylim()[0] * 0.8, 'Number of States = %d' % degs.shape[0] )
xlabel('Degree of state')
ylabel('Free energy (kT)')
xlim([1,xlim()[1]])
xscale('log')
title('Free energy vs degree of state')
legend()
savefig(os.path.join( options.write_dir,'GvsDeg.png') )

figure()

degDist = bincount( degs )
#plot( degDist[ where( degDist > 0 )], '.' )
hist( degs, bins=10**linspace(0,4) )
xlim([ 0, xlim()[1] ])
title('Degree Distribution')
xlabel('Degree')
ylabel('Frequency')
xscale('log')
savefig(os.path.join( options.write_dir, 'DegDist.png') )
 
