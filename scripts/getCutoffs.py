#!/usr/bin/env python

# This script will calculate the cutoffs for native contacts such that the given
#   contact is satisfied 85% of the time for all xtc's that are in the given directory

from optparse import OptionParser
import os
import sys
from msmbuilder import Trajectory
from numpy import *

parser = OptionParser()

parser.add_option('-d','--directory-input',default='.',dest='input_dir',help='Directory holding trajectories to analyze')
parser.add_option('-p','--percentage', default=85,type='float',dest='percent',help='Percentage to use in calculating the contacts')
parser.add_option('-c','--contact-file',dest='contactFN',help='Filename with a list of native contacts (from SMOG, UCSD Onuchic Group)')
parser.add_option('-o','--output-file',dest='outputFN',default='cutoffs.txt',help='Filename to write the output to')
parser.add_option('-s','--pdb-file',dest='pdbFN',help='PDB to use to read XTCs')
options, args = parser.parse_args()


def analyzeConf(conf, contactList):
   # This function analyzes a given conformation and returns a numpy array of the distance between residues
   # given in the contactList:

   radii = []
   for atomA,atomB in contactList:
      diff = conf[atomA-1] - conf[atomB-1]
      temp = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]
      temp = sqrt(temp) 
      radii.append( temp )

   return radii

def analyzeTraj(traj,contactList):
   outList = []
   # This script will calculate the distance of each pair listed in the contact list for the whole trajectory
   for frame in traj['XYZList']:
      outList.append(array( analyzeConf(frame,contactList) ) )

   return outList

def main():
   contactFile = open( options.contactFN, 'r' )
   lines = contactFile.readlines()
   contactList = []
   print "Reading contacts from %s ..." % options.contactFN
   for line in lines[1:]:
      contactList.append( ( int(line.split()[1]), int(line.split()[3]) ) )

   # Next generate the list of xtc's:
   print "Finding xtc's in %s ..." % os.path.abspath(options.input_dir)
   dirList = os.listdir(options.input_dir)
 
   xtcList = [ thing for thing in dirList if thing[-4:] == '.xtc' ]

   if not xtcList:
      print "No xtc files in %s ... Exiting ..." % os.path.abspath(options.input_dir)
      exit()

   datList = []

   for trajFile in xtcList:
      traj = Trajectory.Trajectory.LoadFromXTC(os.path.join( options.input_dir,trajFile),PDBFilename=options.pdbFN)
      datList.extend(array( analyzeTraj(traj,contactList) ) )

   datList = array(datList)

   datList = datList.T

# Now the rows correspond to each contact, and the columns represent time snapshots
   cutoffs = []
   for contactLine in datList:
      N = len(contactLine)
      temp = contactLine
      temp.sort()
      COindex = int( round( options.percent / 100. * N , 0 ) )
      cutoffs.append(temp[ COindex ] )

      test = 0
      for item in temp:
         if item < cutoffs[-1]:
            test += 1

      test = float(test) / N * 100.
      if not test - options.percent < 0.1:
         print "Cutoff defined to %5.2f gives a percentage contacting of %5.2f" % (cutoffs[-1], test)

   print "Writing output to %s ..." % options.outputFN
   output = open(options.outputFN,'w')
   for contact,cutoff in zip(contactList,cutoffs):
      output.write("%d %d %5.2f\n" % ( contact[0], contact[1], cutoff ) )

   output.close()

if __name__ == '__main__':
   main()
