#!/usr/bin/env python

import datetime
import sys
from numpy import *
from matplotlib.pyplot import *
from random import random
import re

now = datetime.datetime.now()

dataFile = open(sys.argv[1],'r')
runs = int(sys.argv[2])
clones = int(sys.argv[3])
col = int(sys.argv[4])
labelList = ['frame','rmsAhx','radius','twist','rise','lenAhx','sstruc','Percent Helix']
# IMPORTANT: The above line needs to be updated everytime getSQLdb.py is updated to include more columns...

HelixMatch = 'H'

print "Warning, dataFile must contain a dictionary!"

data = eval(dataFile.read())

figure()
runColor = []
for i in range(runs):
   plotColor = ( random(), random(), random() )
   runColor.append(plotColor)
   for j in range(clones):
      tempList = array(data[(i,j)])
      if len(tempList):
         # then we have data for this run,clone
         plot( [ float(thing) / 2. for thing in tempList[:,0] ] ,[ float(thing) for thing in tempList[:,col] ], color=plotColor, lw=3, ls=':' )

xlabel('Time (ns)')

ylabel(labelList[col])

dbName = sys.argv[1].split('.')[0]
title(dbName + ' ' + now.strftime("%Y-%m-%d"))

savefig(dbName+'.pdf')

close(figure())


# Now. I want to plot the trajectories' helical properties like this:
  # x-axis helix start, y-axis helix end. This way we might be able to look at the conformations!

# First, I need to find the longest helical part, and then find out where it starts and ends.

longHelixPart = {}

for i in range(runs):
   for j in range(clones):
      tempList = array(data[(i,j)])
      helixPart = []
      if len(tempList):
         for thing in tempList:
            string = thing[6] # This is the dssp string for this frame.
            matchStr = re.split('([\HelixMatch]+)',string)
            HlxList = [ matchItem for matchItem in matchStr if re.search('[\HelixMatch]+',matchItem) ]
            if HlxList:
               HlxLens = [ len(hlxItem) for hlxItem in HlxList ]
               longHlx = HlxList[ HlxLens.index(max(HlxLens)) ]
               HelixStart = string.find(longHlx) + 1
               HelixEnd = HelixStart + len(longHlx) - 1
               helixPart.append( [thing[0], HelixStart, HelixEnd ] )
      longHelixPart[(i,j)] = helixPart

# Now plot the data!

figure()

for i in range(runs):
   plotColor = runColor[i]
   for j in range(clones):
      tempList = array(longHelixPart[(i,j)])
      if len(tempList):
         plot( [ thing[1] for thing in tempList], [ float(thing[2]) - float(thing[1]) for thing in tempList ],'.', color=plotColor,lw=1,ls='-' )


xlabel('Start Helix')
ylabel('End Helix')
title('Longest Helix')


savefig(dbName+'LongHlx.pdf')


