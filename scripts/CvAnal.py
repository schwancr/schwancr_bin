#!/usr/bin/env python
from numpy import *
from subprocess import *
import sys
import os
import re
'''
This script goes in a directory (sys.argv[1]) that contains a number of .trr and .edr files
with the same prefixing name. 

The script runs g_energy to calculate Cv of those trajectories, and then 
outputs the Cv vs T data to sys.argv[2]
'''


# First generate the list of names:

files = os.listdir(sys.argv[1])
os.putenv("GMX_MAXBACKUP","-1")
edrs = []

for file in files:
   if file[-3:] == 'edr':
      edrs.append(os.path.join( sys.argv[1],file) )


Cvs = []
Temps = []

for file in edrs:
 
   print "Working on file %s " % file

   cmd = 'echo "10 11" | g_energy -f '+file+' -nmol 1 -nconstr 0 2> /dev/null'


   nrgyOut = check_output(cmd, shell=True).split('\n')


   for line in nrgyOut:
      matchObj = re.search('Heat capacity at constant volume Cv: ([\s\d.]*)',line)
      if matchObj:
         try: Cvs.append(float(matchObj.group(1)))
         except: 
            print "Problem in g_energy: Cv = %s " % matchObj.group(1)
            continue

      matchObj2 = re.search('Temperature\s\s\s*([\d.]*)',line)
      if matchObj2:
         try: Temps.append(float(matchObj2.group(1)))
         except: 
            print "Problem in g_energy: T = %s " % matchObj.group(1)
            Cvs.pop(-1)
            continue


sol = array( zip( Temps, Cvs ) )

savetxt( sys.argv[2], sol )

