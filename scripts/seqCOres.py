#!/usr/bin/python2.6 -tt -u

from numpy import *
import sys
import re



#'''
#This script reads in a pdb file, then calculates the CO of segments starting from the N terminus.
#
#The output is a txt file with entries corresponding to the number of amino acids considered (from the N-terminus)
#vs the CO calculated.
#'''


def CO(resList):
   # resList is a list of residues, where each residue is a numpy array of coordinates.

   L = len(resList)

   contacts = zeros((L,L))
   resNums = range(L)
   for i in resNums:
      colNums = range(i+1,L)
      for j in colNums:
         contacts[i,j] = isContact(resList[i],resList[j])
   return contacts      

def isContact(resA, resB):
   # This function will decide if atoms in resA and resB are within some cutoff of eachother

   cutoff = 36.
   sol = 0
   for atomA in resA:
      for atomB in resB:
         diff = atomA - atomB
         test = dot(diff,diff)
         if test <= cutoff:
            return 1
         elif test >= 1000:
            return 0

   return 0


def main():

   # first read in the pdb:
   print "Reading structure file: "+sys.argv[1]+"..."

   pdb = open(sys.argv[1],'r')
   out = open(sys.argv[1][:-4]+'resCO.txt','w')
   pdbLines = pdb.readlines()
   
   protein = []
   oldRes = -1
   resTemp = []
   count = 0
   for line in pdbLines:
      if line[:3] == 'TER':
         break
      if line[:4] != 'ATOM':
         continue
      newRes = line[22:26]
      if newRes == oldRes:
         resTemp.append([float(line[30:55].split()[i]) for i in range(3) ])
      else:
         if resTemp:
            protein.append(array(resTemp))
         count = count + len(resTemp)
         resTemp = []
         resTemp.append(array([float(line[30:55].split()[i]) for i in range(3) ]))
      oldRes = newRes
   
   protein.append(array(resTemp)) # In order to add the last residue.

   print "Calculating contacts matrix..."

   contactMat = CO(protein)

   print "Analyzing matrix, and calculating the contact orders"
   pdb.close()
   numRes = len(contactMat)

   sumContacts = 0
   seqDiffSum = 0
   
   for n in range(1,numRes):
      column = n
      for i in range(n):
         seqDiffSum = (n - i) * contactMat[i,n] + seqDiffSum

      sumContacts = sum( contactMat[:,n] ) + sumContacts 
      
      contactOrder = seqDiffSum / float(n+1) / float(sumContacts)

      outLine = str(n+1).rjust(5) + "  %8.5f\n" % contactOrder
      out.write(outLine)

    
   out.close()

if __name__=="__main__":
   main()
