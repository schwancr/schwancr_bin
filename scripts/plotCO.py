
from numpy import *
from matplotlib.pyplot import *

import sys
import re
import os


files = os.listdir('.')

dataFiles = []

for file in files:
   if re.search('CO\.txt',file):
      dataFiles.append(file)

dataS = []

for file in dataFiles:
   dataS.append( genfromtxt(file) )

n = len(dataFiles)

numFigs = n / 9 + 1

figure(1)
count = 0
for index,data in enumerate(dataS):
   if not index % 9:
      figure(index / 9 + 1)
  #    print "NEW PLOT!"

   
   subplot(3,3, index % 9 + 1 )
  # print index % 9 + 1, index / 9 + 1
   plot(data[:,0],data[:,1],label=dataFiles[index][:4]+'REG')
   plot(data[:,0],data[:,2],label=dataFiles[index][:4]+'NORM')
   plot(data[:,0],data[:,3],label=dataFiles[index][:4]+'PER_RES')
   
   title(dataFiles[index][:4])
   
   count += 1

