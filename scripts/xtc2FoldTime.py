#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-f',dest='traj_FN',help='A trajectory file that GROMACS can understand.' )
parser.add_option('-s',dest='struc_FN',help='A structure file to use to define the native state. GROMACS must be able to read it.' )
parser.add_option('-o',dest='out_FN',help='Output filename to write Folded and unfolded times. This actually writes two data files with this as the root name.')
parser.add_option('--fc',dest='f_cut',type=float,help='Folded cutoff')
parser.add_option('--uc',dest='u_cut',type=float,help='Unfolded cutoff')

options, args = parser.parse_args()
 
from numpy import *
from pyschwancr import dataIO, msmTools
import os, sys, re
 
# first make the xvg.

os.system('echo "0 0" | g_rms -f %s -s %s -o %s_rmsd.xvg' % ( options.traj_FN, options.struc_FN, options.traj_FN ) )

xvgIn = open( '%s_rmsd.xvg' % options.traj_FN, 'r' )

dat = []
for line in xvgIn:
	if line[0] in [ '#', '@' ]:
		continue
	else:
		dat.append( [ float( i ) for i in line.split() ] )

dat = array( dat )

hist( dat[:,1], bins=100 )
savefig('%s_dist.pdf' % options.traj_FN )
 #Find the correct cutoffs to use... I'm going to simply use the two maxima. 

Folds, Unfolds = msmTools.calcRawFoldTime( dat[:,1], options.f_cut, options.u_cut, low_is_folded = True )

savetxt( options.out_FN[:-4] + '_FoldTimes.dat', Folds )
savetxt( options.out_FN[:-4] + '_UnfoldTimes.dat', Unfolds )
