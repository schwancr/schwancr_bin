#!/usr/bin/env python

import datetime
import sys
from numpy import *
from matplotlib.pyplot import *
from random import random
from mpl_toolkits.mplot3d import Axes3D


now = datetime.datetime.now()

dataFile = open(sys.argv[1],'r')
runs = int(sys.argv[2])
clones = int(sys.argv[3])
col = int(sys.argv[4])
labelList = ['frame','rmsAhx','radius','twist','rise','lenAhx','sstruc','Percent Helix']
# IMPORTANT: The above line needs to be updated everytime getSQLdb.py is updated to include more columns...

print "Warning, dataFile must contain a dictionary!"

data = eval(dataFile.read())

figure()
runColor = []


subplotLoc = { 1 : (1,1),
               2 : (1,2),
               3 : (1,3),
               4 : (2,2),
               5 : (2,3),
               6 : (2,3),
               7 : (3,3),
               8 : (3,3),
               9 : (3,3),
              10 : (3,4),
              11 : (3,4),
              12 : (3,4),
              13 : (3,5),
              14 : (3,5),
              15 : (3,5),
              16 : (4,4) }
if not runs in subplotLoc:
   print "ERROR! You should plot this many runs on multiple subplots..."
   exit()

rows = subplotLoc[runs][0]
columns = subplotLoc[runs][1]



# First get the longest time frame:



maxTime = 0

for i in range(runs):
   for j in range(clones):
      if len(data[(i,j)]) > maxTime:
         maxTime = len(data[(i,j)])
plotSS = {}


xs = range(1,40)
ys = {}

for t in range(maxTime):
   for i in range(runs):
      plotColor = ( random(), random(), random() )
      runColor.append(plotColor)
      SSarray = zeros(39)  # SSarray holds the counts of helical residues. It will be calculated across all clones and normalized by the number of clones at a particular timestep.
      helixCount = 0
      for j in range(clones):
         tempList = array(data[(i,j)])
         if len(tempList) > t: 
         # then we have data for this run,clone, frame
            frame = tempList[t]
            helical = zeros(39) + array( [ 1 if letter in 'THIG' else 0 for letter in frame[6] ] )
            helixCount += 1
            SSarray += helical
      if not helixCount:
         ys[(i,t)] = array([0.]*39)
         continue
      ys[(i,t)] = SSarray / float( helixCount )


# Now plot everything! (For each run.)


subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98,wspace=0.1,hspace=0.1)


for i in range(runs):
   ax = subplot(rows,columns,i+1,projection='3d')
   for t in range(maxTime):
      ax.bar(xs,ys[(i,t)],t/2.,zdir='y',color=list(['gray']*39), alpha=0.5)  
   
   ax.set_xlim3d(0,40)
   ax.set_ylim3d(0,(maxTime+1)/2)
   ax.set_zlim3d(0,1)
   ax.azim = 60
   ax.elev = 30
   ax.set_xlabel('Residue ID')
   ax.set_ylabel('Time (ns)')
   ax.set_zlabel('Percent Helix over all Clones')
   ax.set_title('Run '+str(i))


dbName = sys.argv[1].split('.')[0]

savefig(dbName+'SS.pdf')

