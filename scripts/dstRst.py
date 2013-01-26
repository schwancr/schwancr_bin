#!/usr/local/bin/python2.7 -tt

# -------------------------------------------------------------------------------
# This script will assign distance restraints for a ribosome tunnel. The distance
# is calculated away from an interior piece (nascent peptide) specified by atom
# indices. Atoms will be assigned a force constant based on the distance away 
# from the interior atoms. 
# The .gro file, .top file, and ALL .itp files (not including posres files) 
# must be in the current directory.
# -------------------------------------------------------------------------------


import sys
import re
import numpy
import time
# INPUT:
# 1) Filename of .gro file (groFile)
# 2) The radius to test points from the interior (radius) [nm]
# 3) Filename of .top file (topFile)
# 4) Filename of the .gro file containing the interior peptide (interFile)
# OUTPUT:
# 1) dstrst_"...".itp which include the force constants for the particular
#    chain "..."
# 2) restrained_"..."_radiusA.top corresponding to the new topology file.
# 3) restrained_atoms_radiusA.gro corresponding to the atoms restrained.

def main():
  print sys.argv[0]

  if sys.argv[1] == 'help' or sys.argv[1] == '-h':
    print "   This script reads in 4 input arguments and outputs a number of files"
    print " INPUT:\n\t1) Filename of .gro file (groFile)" 
    print " \t2) The radius to test points from the interior in nanometers (radius)"
    print " \t3) Filename of .top file (topFile)"
    print " \t4) Filename of the .gro file containing the interior peptide (interFile)"
    print "  The script determines what atoms in the groFile are greater than"
    print "the radius away from any interior atoms. If such atoms exist, they"
    print "are given position restraints in the new dstrst_...itp files"
    print "\n  Additionally, the script adds a line to each itp file,"
    print "which has the #ifdef DSTRST statements necessary to run with the"
    print "new dstrst files."

    exit()


  timeIn = time.clock()
  groFile = open(sys.argv[1],'r')
  radius = float(sys.argv[2])
  topFile = open(sys.argv[3],'r')
  interFile = open(sys.argv[4],'r')

  groOutTest = open('restrained_atoms_'+str(int(radius*10))+'A.gro','w')
  pdbOutTest = open('pdbOutTest.txt','w')
  atomCoords = {}
  intCoords = []
  
  exclAtms = ['MG  ','NA  ','CL  ','HOH ','K   ','SOL ','CD  ']
  for line in groFile:

    if line[5:9] not in exclAtms:
      
      matchObj = re.search(r'([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)',line[20:])
      # Quit if we don't find anything
      if not matchObj:
        continue
      
      x = float(matchObj.group(1))
      y = float(matchObj.group(2))
      z = float(matchObj.group(3))
            
      atmName = line[9:15]
      atmNum = int(line[15:20])
      
      nameObj = re.search(r'\w',atmName)
      
      # Quit if we didn't find a word character
      if not nameObj:
        continue
      if nameObj.group()=='H':
        continue

#      if atmNum >= beginAtm and atmNum <= endAtm:
#        intCoords.append((x,y,z))
        
      atmType = nameObj.group()
      # Quit if it's a hydrogen, because those do not go in the posre file
      
      atomCoords[atmNum] = (x,y,z)
      
      
  # Now we need to go through each of the atmCoords and determine if the
      # atom is far enough to be included in the dstrst file
      # The atoms to be kept will be stored in a list of indices (keepAtms)
  keepAtms = []
  for line in interFile:
    matchObj = re.search(r'([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)',line[20:])
    if matchObj:
      x = float(matchObj.group(1))
      y = float(matchObj.group(2))
      z = float(matchObj.group(3))
      intCoords.append( (x,y,z) )
  count = 0
  radiusSqr = radius*radius
  for atm in atomCoords:
    count+=1
    percent = round(float(count)/len(atomCoords)*100,2)
    print str(percent)+'%\r',
    coord = atomCoords[atm]
    keep = True
    for testAtm in intCoords:
      diff = ( coord[0] - testAtm[0], coord[1] - testAtm[1], coord[2] - testAtm[2] )
      mag = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]
      if mag <= radiusSqr:
        keep = False
        break
    if keep:
      keepAtms.append(atm)
    
  print "You kept " + str(round(float(len(keepAtms))/len(atomCoords)*100,2))+ "%"
  # Now we must write the actual files. We have a list of the atom numbers
      # But we need to determine which topology each belongs to. 
      # First we get the list of .itp files the topology uses:

  # The following is an output .gro file that gives the locations of all the
  # atoms RESTRAINED.
  # NOTE: H's are not included, since they are not restrained
  groOutTest.write("Cut from greater than " + str(radius*10) +" Angstroms away\n")
  groOutTest.write(str(len(keepAtms))+'\n')
  groFile.seek(0)
  for line in groFile:
    if line[0] == ' ':
      atmTestNum = int(line[15:20])
      if atmTestNum in keepAtms:
        groOutTest.write(line)

  groOutTest.write(line) # Write the last line so the box is included
  outTopFile = open('restrained_'+sys.argv[3],'w')
  
  atomCounted = 0
  for topLine in topFile:
    matchObj = re.search(r'"(\w+\.itp)"',topLine)
    if matchObj:
      print matchObj.group(1)
      # We have the name of the first itp file. The atoms are listed in order
      # of these itp files. So if we go through them, iterating the counter
      # each time we can count through all kept atoms and then output they're
      # chain ID number
      
      # We also need to rewrite the top file
      outTopFile.write(topLine)
#      outTopFile.write('#ifdef DSTRST\n#include "dstrst_'+matchObj.group(1)+'"\n#endif\n')
      outItp = open('dstrst_'+matchObj.group(1),'w')
      inItp = open(matchObj.group(1),'r+')
      outItp.write('[ position_restraints ]\n')
      for itpLine in inItp:
        atomsCheck = re.search(r'\[ atoms \]',itpLine)
        if atomsCheck:
          break
        # We quit iterating if we find the atoms section

      for itpLine in inItp:

        # We need to make sure we haven't finished the atoms:
        bondsCheck = re.search(r'\[ bonds \]',itpLine)
        if bondsCheck:
          outItp.close()
          break

        # We also need to skip over comment lines:
        if itpLine[0] == ';':
          continue

        # We are in the atoms section AND it is not a comment line so we must count
        # an atom that we are looking at:
        
        atom = re.search(r'(\d+)\s+(\w)',itpLine)
        if atom:
          atomCounted+=1
          atomNum = int(atom.group(1))
          atomType = atom.group(2)
 #         if atomType == 'H':
 #           continue
          if atomCounted in keepAtms:
            outputStr = '{0:{align}{fill}5}     1 10000 10000 10000\n'.format(atomNum,align='<',fill=' ')
            outItp.write(outputStr)

        # This means we have reached the end of the atoms in this itp file and must move
        # to the next one
      ifDef = True
      for itpLine in inItp:
        check = re.search('#ifdef\sDSTRST',itpLine)
        if check:
          ifDef = False
        continue
      # We only want to write this if it is not already written.
      if ifDef:
        inItp.write('#ifdef DSTRST\n#include "dstrst_'+matchObj.group(1)+'"\n#endif\n')
      
      inItp.close()
    else:
      outTopFile.write(topLine)
    

  groFile.close()
  topFile.close()
  
  groOutTest.close()
  pdbOutTest.close()
  timeOut = time.clock()
  print 
  print timeOut - timeIn

if __name__ == '__main__':
  main()
