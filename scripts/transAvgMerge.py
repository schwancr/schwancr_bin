#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-f','--input-file',dest='inFN',help='Input (flat-txt) file with assignment data in rows')
parser.add_option('-m','--mapping-file',default='./Mapping.dat',dest='mapFN',help='Mapping file that tells which states in old Assignments.h5 correspond to which in Assignments.Fixed.h5')
parser.add_option('-o','--output-file',dest='outFN',help='Output (flat-txt) file')
parser.add_option('-p',dest='popFN',default='./clusterPops.dat',help='Cluster populations (COUNT NOT EQUILIBRIUM POPS)')
options, args = parser.parse_args()

from numpy import *
import numpy as np

# First load map and check if anything was even trimmed...
M = loadtxt(options.mapFN)
maxState=M.max()
print M.shape, maxState
clusterPops = loadtxt( options.popFN )
print clusterPops.shape
if options.outFN:
	outFN = options.outFN
elif options.inFN.split('.')[-1] in [ 'txt','dat' ]:
# Then we have an extension so output name should use the first part
	outFN = '.'.join( options.inFN.split('.')[:-1] ) + '.Fixed.dat'
else:
	outFN = options.inFN + '.Fixed.dat'

data = loadtxt(options.inFN)

outList = np.zeros( (maxState+1,3) ) 
outList[:,0] = np.arange( maxState+1 )
outList[:,2] = np.ones( maxState+1 ) * -1
Ns = np.zeros( maxState+1 )

for index, line in enumerate(data):
	if M[index] >= 0:
		outList[ M[index] ] += clusterPops[index] * line[1]
		Ns[ M[index] ] += clusterPops[index]

outList[:,1] /= Ns
savetxt( outFN, array(outList) )
