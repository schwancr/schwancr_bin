#!/usr/bin/env python

import datetime
import sys
from numpy import *
from matplotlib.pyplot import *
from random import random

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

