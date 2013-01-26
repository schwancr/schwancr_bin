#!/usr/bin/python2.6 -tt -u

from numpy import *
import sys
import re

from matplotlib.pyplot import *


#'''
#This script reads in a pdb file, then calculates the CO of segments starting from the N terminus.
#
#The output is a txt file with entries corresponding to the number of amino acids considered (from the N-terminus)
#vs the CO calculated.
#'''

# Output is in 4 columns:
# 1) Number of residues considered
# 2) Sequential contact order
# 3) Sequential contact order / total contacts of those residues
# 4) CO per residue

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
            sol = sol + 1
         elif test >= 1000:
            return 0
   return sol


def main():

   # first read in the pdb:
   print "Reading structure file: "+sys.argv[1]+"..."

   pdb = open(sys.argv[1],'r')
   out = open(sys.argv[1][:-4]+'atomCO.txt','w')
   pdbLines = pdb.readlines()
   
   protein = []
   oldRes = -1
   resTemp = []
   count = 0
   readAtoms = False
   for line in pdbLines:
      if line[:3] == 'TER' and readAtoms:
         break
      if line[:4] != 'ATOM':
         continue
      NAtest = line[18:21]
      if NAtest.strip()[0] in ['D','R','d','r']:
         continue
      readAtoms = True
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

#   a = zeros((numRes,numRes))
#
#   for i in range(numRes):
#      for j in range(numRes):
#         if contactMat[i,j] > 0:
#            a[i,j] = 1
#            a[j,i] = 1
#   figure()
#   imshow(a)
#   title(sys.argv[1][:-4])
   savetxt(sys.argv[1][:-4]+'atomCM.txt',contactMat)
   proteinContacts = sum( contactMat )   

   for n in range(1,numRes):
      column = n
      perResCO = 0

      for i in range(n):
         seqDiffSum = (n - i) * contactMat[i,n] + seqDiffSum
      
      for i in range(n):
         perResCO = (n - i) * contactMat[i,n] + perResCO

      for i in range(n+1,numRes):
         perResCO = (i - n) * contactMat[n,i] + perResCO

      sumContacts = sum( contactMat[:,n] ) + sumContacts 
      totalContacts = sum( contactMat[:n,:] )
      
      
      
      contactOrder = { '0': seqDiffSum / float(n+1) / float(sumContacts) , '1': seqDiffSum / float(n+1) / float(totalContacts), '2': perResCO / proteinContacts }

      outLine = str(n+1).rjust(5) + "\t%(0)8.5f\t%(1)8.5f\t%(2)8.5f \n" % contactOrder
      out.write(outLine)

    
   out.close()

if __name__=="__main__":
   main()
