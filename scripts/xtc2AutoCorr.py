#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-f',dest='traj_FN',help='A trajectory file that GROMACS can understand.' )
parser.add_option('-s',dest='struc_FN',help='A structure file to use to define the native state. GROMACS must be able to read it.' )
parser.add_option('--uc',dest='u_cut',type=float,help='Unfolded cutoff')

options, args = parser.parse_args()
 
from numpy import *
from pyschwancr import dataIO, msmTools
from msmbuilder import autocorrelate
import os, sys, re
 
# first make the xvg.

os.system('echo "0 0" | g_rms -f %s -s %s' % ( options.traj_FN, options.struc_FN ) )

xvgIn = open( 'rmsd.xvg', 'r' )

dat = []
for line in xvgIn:
	if line[0] in [ '#', '@' ]:
		continue
	else:
		dat.append( [ float( i ) for i in line.split() ] )

dat = array( dat )

autoCorr = autocorrelate.fft_autocorrelate( dat[:,1] )

savetxt( options.traj_FN[:-4] + '_autocorr.dat', autoCorr )
