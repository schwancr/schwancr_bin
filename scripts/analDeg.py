#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='dir',help='Directory to look for data (stateAvg_Qtot.dat.Fixed, Degrees.dat, Populations.dat)')

options, args = parser.parse_args()
 
import numpy as np
from pyschwancr import dataIO
import os, sys, re
 
pops = dataIO.readData( os.path.join( options.dir, 'Populations.dat' ) )
deg = dataIO.readData( os.path.join( options.dir, 'Degrees.dat' ) )
try:
   q = dataIO.readData( os.path.join( options.dir, 'stateAvg_Qtot.dat.Fixed' ) )[:,1]
   isF = q > 0.4
   isU = q <= 0.4
except:
   try:
      r = dataIO.readData( os.path.join( options.dir, 'stateAvg_RMSD.Fixed.dat' ) )[:,1]
      isF = r < 0.4
      isU = r > 0.4
   except:
      print "Need either stateAvg_[Qtot,RMSD].Fixed.dat in the directory"
      exit()

sumDegF = int( deg[ isF ].sum() )
sumDegU = int( deg[ isU ].sum() )

normDegF = ( deg[ isF ] * pops[ isF ] ).sum()
normDegU = ( deg[ isU ] * pops[ isU ] ).sum()

avgDegF = deg[ isF ].mean()
avgDegU = deg[ isU ].mean()

normAvgDegF = ( deg[ isF ] * pops[ isF ] ).sum() / pops[ isF ].sum()
normAvgDegU = ( deg[ isU ] * pops[ isU ] ).sum() / pops[ isU ].sum()

sumDegOverPopF = sumDegF / pops[ isF ].sum()
sumDegOverPopU = sumDegU / pops[ isU ].sum()

print "U (%d - %f): Sum of Degrees %10d -- Degree \dot Pops %20f -- Average Degree %20f -- Normalized Degree \dot Pops %20f -- Sum of Degrees over the sum of populations %20f" % ( len( np.where( isU )[0] ), sum(pops[isU]), sumDegU, normDegU, avgDegU, normAvgDegU, sumDegOverPopU )
print "F (%d - %f): Sum of Degrees %10d -- Degree \dot Pops %20f -- Average Degree %20f -- Normalized Degree \dot Pops %20f -- Sum of Degrees over the sum of populations %20f" % ( len( np.where( isF )[0] ), sum(pops[isF]), sumDegF, normDegF, avgDegF, normAvgDegF, sumDegOverPopF )
