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

#   contacts = zeros((L,L))

#   for i in range(L):
#      for j in range(i+1,L):
#         contacts[i,j] = isContact(resList[i],resList[j])
         
#   contactOrder = 0

#   N = sum(contacts)
#   print N
#   for i in range(L):
#      for j in range(i+1,L):
#         contactOrder = contactOrder + (j - i) * contacts[i,j]

   atoms = []
   getRes = {}
   contactOrder = 0
   atmNum = 0
   for index,thing in enumerate(resList):
      for atom in thing:
         atoms.append(atom)
         #if atmNum == 13:
         #   print atmNum, index, atom, thing
         getRes[atmNum] = index
         atmNum += 1
   #print resList
   #print "------------------------\n",atoms
   N = 0
   #print len(atoms)
   for i in range(len(atoms)):
      for j in range(i+1,len(atoms)):
         if getRes[j] - getRes[i] <= 0:
            continue        
         if linalg.norm(atoms[i]-atoms[j]) <= 6.:
            contactOrder = contactOrder + getRes[j] - getRes[i]
            N += 1
 #  print L
 #  print N
 #  print contactOrder        
 #  print contactOrder / N
 #  print contactOrder / N / L
 #  exit()

   return contactOrder / float(N) / float(L)
   
def isContact(resA, resB):
   # This function will decide if resA and resB are within some cutoff of eachother

   cutoff = 6.
   sol = 0
   for atomA in resA:
      for atomB in resB:
         if linalg.norm(atomA - atomB) <= cutoff:
            sol = sol + 1

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
   for line in pdbLines:
      if line[:3] == 'TER':
         break
      if line[:4] != 'ATOM':
         continue
      newRes = line[22:26]
      if newRes == oldRes:
         #print line
         resTemp.append([float(line[30:55].split()[i]) for i in range(3) ])
      else:
         if resTemp:
            protein.append(array(resTemp))
         count = count + len(resTemp)
         resTemp = []
         resTemp.append(array([float(line[30:55].split()[i]) for i in range(3) ]))
         #print "--------------------NEW RESIDUE"
         #print line
      oldRes = newRes
   
   protein.append(array(resTemp)) # In order to add the last residue.
   count = count + len(resTemp)
   #print "Atoms=",count 
   #for thing in protein:
   #   print "---------------------NEW RESIDUE"
   #   print thing

  # CO(protein)

   pdb.close()
   COs = {}   
   numRes = len(protein)
   for n in range(2,numRes):
      COs[n] = CO(protein[:n])
      print str(n).rjust(5),'\r',
      outLine = str(n).rjust(5) + "  %5.3f\n" % COs[n]
      #print outLine
      out.write(outLine)

    
   out.close()

if __name__=="__main__":
   main()
