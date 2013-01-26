#!/usr/bin/env python


from optparse import OptionParser

parser = OptionParser()

parser.add_option('--nc','--native-contact-file',dest='ncFN',help='File with a list of native contacts (From Onuchic Grp)')
parser.add_option('--rr','--residue-ranges',dest='res_ranges',action='append',help='List of residue ranges to check contacts between\n Number of arguments should be even.\n'+
											'SYNTAX: --rr X,Y --rr W,V.\n This will check for contacts between groups [X,Y] and [W,V].\n'+
											' Note that Y and V ARE included in the two groups.')
parser.add_option('-t','--traj-directory',dest='trj_dir',default='.',help='Directory to look for trajectories in')
parser.add_option('-o','--output-file',dest='outFN',default='NC_calc',help='Output filename')
options, args = parser.parse_args()

import os
from msmbuilder import Serializer
from matplotlib.pyplot import *
from numpy import *
import re

def getNClist(resRanges, ncFN):

   pairs = []

   if len(resRanges) % 2 > 0:
      print "You need to input an even number of ranges with --rr"

   for i in range( len(resRanges) / 2 ):
      a = range(resRanges[2*i][0], resRanges[2*i][1] + 1 ) 
      b = range(resRanges[2*i + 1][0], resRanges[2*i+1][1] + 1 )

      pairs.append( (a,b) )


# Now generate the lists of contacts to calculate for each pair:
   f_contact = open( ncFN, 'r')
   lines = f_contact.readlines()
   NCs = []
   for line in lines[1:]:
      NCs.append( ( int( line.split()[1]), int( line.split()[3] ) ) )
   

   testList = []
   groupList = []
   for (grpA, grpB) in pairs:
      for index, contact in enumerate(NCs):
         if contact[0] in grpA:
            if contact[1] in grpB:
               if contact not in groupList:
                  groupList.append(index)
         elif contact[0] in grpB:
            if contact[1] in grpA:
               if contact not in groupList:
                  groupList.append(index)
      if groupList:
         testList.append(groupList)
      else:
         print " No contacts found in group! "
      groupList = []     
	
   if not testList:
      print "No contacts found between the specified ranges!"
      exit()
   return testList

def readDat(traj, ContactsList ):
	# this function takes in a trajectory numpy array along with a list of contacts and then
   #   outputs the total for each Contacts within ContactsList for each frame in the trajectory
	# So the output array is an N x M array where N = # of frames and M = len(ContactsList) (i.e.
	# number of groups to check
	outList = []
	groups = []
	for frame in traj:
		for contacts in ContactsList:
			N = len(contacts)
			groups.append( sum( frame[ contacts ] ) / float(N) )
		outList.append(array(groups) )
		groups = []
   
	return outList

def main():
	global options
	
	resRanges = [ ( int( thing.split(',')[0] ), int( thing.split(',')[1] ) ) for thing in options.res_ranges ]

 	print "Reading native contacts from %s" % options.ncFN  
	testList = getNClist(resRanges, options.ncFN)
	trajList = [ (int(thing[3:-4]),thing) for thing in os.listdir( options.trj_dir ) if re.match('trj\d+\.npy',thing) ]
	trajList.sort()
	trajList = [ b for (a,b) in trajList ]

	count = 0
	print "Will calculate NC's from %d data files in %s ..." % (len(trajList),os.path.abspath(options.trj_dir))
	sol = []
	for trajFN in trajList:
		print "Calculating for trajectory: %s..." % trajFN
		traj = load(os.path.join( options.trj_dir, trajFN) )
		sol.extend( readDat(traj, testList) )
		count += 1
	print "Saving Data to %s " %options.outFN  
	sol = array(sol)
	save(options.outFN,sol)

	return 
	   

           


if __name__=='__main__':
   main()
