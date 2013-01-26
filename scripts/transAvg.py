#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-f','--input-file',dest='inFN',help='Input (flat-txt) file with assignment data in rows')
parser.add_option('-m','--mapping-file',default='./Mapping.dat',dest='mapFN',help='Mapping file that tells which states in old Assignments.h5 correspond to which in Assignments.Fixed.h5')
parser.add_option('-o','--output-file',dest='outFN',help='Output (flat-txt) file')

options, args = parser.parse_args()

from numpy import *
from pyschwancr import dataIO

# First load map and check if anything was even trimmed...
M = loadtxt(options.mapFN)

if where(M==-1)[0].shape == 0:
	print "No -1 entries in %s ... This probably means nothing was trimmed. Exiting..."%options.mapFN

if options.outFN:
	outFN = options.outFN
elif options.inFN.split('.')[-1] in [ 'txt','dat' ]:
# Then we have an extension so output name should use the first part
	outFN = '.'.join( options.inFN.split('.')[:-1] ) + '.Fixed.dat'
else:
	outFN = options.inFN + '.Fixed.dat'

data = dataIO.readData( options.inFN )

outList = []

for index, line in enumerate(data):
	if M[index] >= 0:
		outList.append( line )

outList = array( outList )

if outList.dtype in [ complex, complex64, complex128 ]:
	save( outFN[:-4]+'.npy', outList )
else:
	savetxt( outFN, array(outList) )
