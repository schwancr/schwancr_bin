#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p',dest='pcaFN',help='PCA dictionary with vecs and vals')
parser.add_argument('-o',dest='outFN',help='Output filename to save the timescales (in frames) for each mode. These are sorted the same as the vals are in the input pca object')
parser.add_argument('--dt',dest='dt',type=int,help='Time given to TrainPCA_cross.py. In frames')
parser.add_argument('-d',dest='divisor',type=float,default=1,help='Divisor to conver the time in frames to a time unit')
parser.add_argument('--units',dest='units',default='frames',help='Units to use in the xlabel and ylabel')
parser.add_argument('--y-lim',dest='y_lim',type=float,nargs=2,default=None,help='Y-limit to use in the plot')

options = parser.parse_args()
 
from msmbuilder import io
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import *
import os, sys, re

#pca = Serializer.LoadFromHDF( options.pcaFN )
pca = io.loadh( options.pcaFN, deferred=False )

vals = pca['vals']
vals[np.where(vals<10**-8)] = 10**-8
vecs = pca['vecs']

timescales = - options.dt / np.log( vals )
np.savetxt(options.outFN, timescales)

timescales /= options.divisor

subplot(132)

hlines(timescales, 0,1, color='black')
xticks([])
yscale('log')

ylabel('tPCA Timescales (%s)'% options.units)

if options.y_lim != None:
   ylim( options.y_lim )

xlabel( 'PCs from lag of %d %s' % ( options.dt / options.divisor, options.units ) )

savefig( options.outFN[:-4]+'.png')
