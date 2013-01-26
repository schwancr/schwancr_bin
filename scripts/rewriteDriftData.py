#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--dir',dest='inDir',help='Input file directory containing files formatted as <metric>_drift_u##.npy')
parser.add_option('-d',dest='metric',default='rmsd',help='Metric used. Will use this to decide the format of the drift files.')
parser.add_option('-o',dest='outFN',help='Output file [ <metric>Data.h5 ]' )

options, args = parser.parse_args()
 
import numpy as np
from pyschwancr import dataIO
import os, sys, re
from msmbuilder import Serializer

regEx = r'^%s_drift_u\d+\.npy' % options.metric.lower()

BeginInd = { 'rmsd' : 12, 'qnorm' : 13, 'dihedral' : 16 }[ options.metric.lower() ]

FileList = dataIO.getTrajList( options.inDir, BeginInd = BeginInd, RegEx = regEx )

Taus = [ int( fn[ BeginInd + len( options.inDir ) + 1 : -4 ] ) for fn in FileList ]

Data = [ np.load( fn ).flatten() for fn in FileList ]

DataAry = np.ones( ( len( Data ), max([ d.shape[0] for d in Data ]) ) ) * -1

for i in xrange( len( Data ) ):
	DataAry[i][:len(Data[i])] = Data[i]

S = Serializer.Serializer( {'Data': DataAry, 'Taus' : Taus } )

if options.outFN:
	outFN = options.outFN
else:
	outFN = '%sData.h5' % options.metric.lower()

S.SaveToHDF( outFN )
print "Saved data to %s" % outFN
