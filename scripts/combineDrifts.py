#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data_FNs',action='append',help='Data file, may have more than one of these entries')
parser.add_option('-c',dest='coefs',action='append',type=float,help='Coefficient corresponding to a data file, used in adding the data together.')

options,args = parser.parse_args()

import numpy as np
from pyschwancr import dataIO

if len( options.data_FNs ) != len( options.coefs ):
   print "Need to enter a coefficient for each filename! There are %d coefficients for %d filenames listed!" % ( len( options.coefs ), len( options.data_FNs ) )
   exit()

Total = []
nameList =[ 'CombinedData' ]

for i in range( len( options.data_FNs ) ):
   fn = options.data_FNs[i] 
   C = options.coefs[i]
   if C % 1: # Nonzero return from mod 1 means this is not an integer, so write as a float
      nameList.append('%.1e'%C)
   else: # It's an integer!
      nameList.append('%d'%int(C))
   nameList.append(fn)
   dat = dataIO.readData( fn )
   if len( dat.shape ) > 1:
      print "Using first column of data... re-shape the data if this doesn't work."
      dat = dat[:,0]
   if len( dat.shape ) == 0:
      dat = np.array([ dat ])

   Total.append( C * dat )


Total = np.array( Total ).sum(axis=0)

name = dataIO.writeData( nameList, Total, txt=False )

print "Wrote the combination to %s" % name


