#!/usr/bin/env python

import sys
from numpy import *
from matplotlib.pyplot import *

from random import random

dataFile = open(sys.argv[1],'r')
runs = int(sys.argv[2])
clones = int(sys.argv[3])
col = int(sys.argv[4])


print "Warning, dataFile must contain a dictionary!"

data = eval(dataFile.read())

figure()

for i in range(runs):
   plotColor = ( random(), random(), random() )
   for j in range(clones):
      tempList = array(data[(i,j)])
      if len(tempList):
         # then we have data for this run,clone
         plot( tempList[:,0],tempList[:,col], color=plotColor, lw=3, ls=':' )


