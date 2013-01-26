#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t','--traj-dir',dest='trj_dir',help='Directory to look for xtc\'s')
parser.add_option('--uc','--unfolded-cutoff',type=float,dest='uCut',help='Cutoff to define the unfolded state')
parser.add_option('--fc','--folded-cutoff',type=float,dest='fCut',help='Cutoff to define the folded state')
parser.add_option('-s','--structure-pdb',dest='pdbFN',help='PDB to use as the native state')
parser.add_option('-P','--procs',dest='procs',default=1,type=int, help='Number of procs to use in the calculation')

options, args = parser.parse_args()

from numpy import *
from msmbuilder.xtc import XTCReader, readxtc
import re
import sys
from msmbuilder import Trajectory, Conformation, DistanceMetric
import os
import multiprocessing as mp

RMSD = DistanceMetric.RMSDMetric()
NatState = Conformation.Conformation.LoadFromPDB(options.pdbFN)
NatState['XYZ'] = NatState['XYZ'] 

def analyzeTrajectory(trajFN):

   #print "Analyzing trajectory: %s ..." % trajFN
   coords = readxtc( trajFN, atomindices = range(128) )
   count = 0
   for frame in coords:
      rms = RMSD.GetDistance( frame, NatState['XYZ'] )
 #     print rms
      if rms >= options.uCut:

         return 1
      elif rms <= options.fCut:
         return 0
      else:
         continue

   return -1

def main():
   global options
   trajList = [ os.path.join( options.trj_dir, thing) for thing in os.listdir(options.trj_dir) if thing[-4:] == '.xtc' ]
   output = []

   pool = mp.Pool(options.procs)
   result = pool.map_async(analyzeTrajectory,trajList)
   result.wait()

   sol = array(result.get())
   savetxt('sol.dat',sol)
   folded = len( where( sol == 0 )[0] ) 
   unfolded = len( where( sol == 1 )[0] )
   neither = len( where( sol == -1 )[0] )
   total = len( sol )
   print "Folded = %d Unfolded = %d Neither = %d Total = %d " % ( folded, unfolded, neither, total )

   return 0
if __name__=='__main__':
   main()      

