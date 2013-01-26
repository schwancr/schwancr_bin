#!/usr/bin/env python

from numpy import *
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import *

import sys
import re

a = open(sys.argv[1],'r')

labels = {}
data = []

print "Loading Data"
for line in a:
   if line[0] in ['#','@']:
      matchObj = re.search('([xy])axis\s*label\s*\"(.*)\"',line)
      if matchObj:
         labels[matchObj.group(1)] = matchObj.group(2)

      continue
   else:
      data.append(array([ float(thing) for thing in line.split() ]))


data = array(data)
print "Plotting Data"
figure()
# Plot each column against the first
for i in range(1,len(data[0])):
   plot(data[:,0],data[:,i],'.')

xlabel(labels['x'])
ylabel(labels['y'])

title(sys.argv[1][:-4])
print "Saving plots to:"
savefig(sys.argv[1][:-4]+'.png')
print "\tTimetrace: %s" % sys.argv[1][:-4]+'.png'
figure()
hist( data[:,1], bins = 100 )
xlabel( labels['y'] )
ylabel( 'Frequency' )
title( sys.argv[1][:-4] + 'Distribution' )

savefig( sys.argv[1][:-4] + 'Dist.png' )
print "\tDistribution: %s" % sys.argv[1][:-4]+'.png'

savetxt( sys.argv[1][:-4] + '.dat', data )
print "\tData: %s" % sys.argv[1][:-4] + '.dat'
