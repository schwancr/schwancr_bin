#!/opt/local/bin/python2.7 -tt

import sys
import re

'''
This script will take in two pdb files: tnl.pdb and rest.pdb

The output pdb (tnlComp.pdb) will be all the atoms in rest.pdb
that are not in tnl.pdb.
'''

tnl = open(sys.argv[1],'r')
rest = open(sys.argv[2],'r')

output = open(sys.argv[1][:-4]+'Comp.pdb','w')

tnlCoords = []
toRemove = []

print "Reading ",tnl.name,"..."

for line in tnl.readlines():
   #if re.search('R[AGUC][N5]',line):
   toRemove.append(line[21:26])
   #tnlCoords.append(line[30:55])

print "Analyizing ",rest.name,"..."

for line in rest.readlines():
   #if line[30:55] in tnlCoords:
   #   continue
   #elif line[21:26] in toRemove:
   #   continue
   if line[21:26] in toRemove:
      continue
   else:
      output.write(line)

print "Done!"
tnl.close()
rest.close()

output.close()

