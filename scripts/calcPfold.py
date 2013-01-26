#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t','--traj-dir',dest='trj_dir',help='Directory to look for xtc\'s')
parser.add_option('--uc','--unfolded-cutoff',type=float,dest='uCut',help='Cutoff to define the unfolded state')
parser.add_option('--fc','--folded-cutoff',type=float,dest='fCut',help='Cutoff to define the folded state')
parser.add_option('--nc','--native-contact-file',dest='nat_contactFN',help='Native contact file with a list of native contacts')
parser.add_option('--co','--cutoff-file',dest='cutoffFN',help='Cutoff file for a cutoff of each native contact that the user can define\n The default will be all contacts have a 6A cutoff.')
parser.add_option('-P','--procs',dest='procs',default=1,type=int, help='Number of procs to use in the calculation')

options, args = parser.parse_args()

from numpy import *
from msmbuilder.xtc import XTCReader, readxtc
import re
import sys
from msmbuilder import Trajectory
import os
import multiprocessing as mp

def getContacts():
   global options

# First parse the contact file:
   contactFile = open(options.nat_contactFN,'r')
   lines = contactFile.readlines()
   contactList = []

   for line in lines[1:]:
      # Skip the first line of the file (it doesn't have contacts on it)
      contactList.append( ( int(line.split()[1]) , int(line.split()[3]) ) )
   
   contactFile.close()
   cutoffFile = open(options.cutoffFN,'r')
   cutoffDict = {}
   for line in cutoffFile.readlines():
      [atomA,atomB] = [ int(thing) for thing in line.split()[:2] ]
      cutoff = float( line.split()[2] )
      if (atomA,atomB) in cutoffDict or (atomB,atomA) in cutoffDict:
         print "(%d,%d) already in the cutoff dictionary! Using the oldest cutoff defined in cutoff file"
      else:
         cutoffDict[(atomA,atomB)] = cutoff*cutoff # Return the cutoff squared so I don't have to sqrt
   return contactList,cutoffDict

def analyzeTrajectory(trajFN):
   global options
   global natContacts
   global cutoff2Dict
# this function will analyze a trajectory for native contacts in the toTestLists

   #print "Analyzing trajectory: %s ..." % trajFN
   coords = readxtc( trajFN, atomindices = range(128) )
   count = 0
   for frame in coords:
      frameSum = 0.
      count += 1
      for contact in natContacts:
         atomA = contact[0]
         atomB = contact[1]

      # Need to check the distance between these atoms in all frames:
         tempArray = []
          
            
         coordA = frame[atomA-1]
         coordB = frame[atomB-1]
 
         diff = [ (thingA - thingB) for (thingA,thingB) in zip(coordA,coordB) ]          

         r2 = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]
         try: cutoff2 = cutoff2Dict[ (atomA,atomB) ]
         except: 
            try: cutoff2 = cutoff2Dict[ (atomB,atomA) ]
            except: 
               print "Cutoff not found for %d,%d, using default = 1.0" % (atomA,atomB)
               cutoff2 = 1
         frameSum += ( r2 <= cutoff2 )
      frameSum /= float(len(natContacts))
      if frameSum <= options.uCut:
         return 1
      elif frameSum >= options.fCut:
         return 0
      else:
         continue

   return -1

def main():
   global options
   global natContacts
   global cutoff2Dict
# Now generate a list of contacts to test. This will really be a list of lists, one for each pair
# And generate the cutoff dictionary that takes in a contact and then finds a cutoff
   natContacts, cutoff2Dict = getContacts()
# Now I have read in all the data I need, I need to read in a trajectory and then analyze it.
   output = []
   trajList = [ os.path.join( options.trj_dir, thing) for thing in os.listdir(options.trj_dir) if thing[-4:] == '.xtc' ]


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

