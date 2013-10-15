#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-t',dest='T_FN',default='./tProb.mtx',help='Transition probability matrix [ ./tProb.mtx ]')
parser.add_argument('-n',dest='num_vals',default=100,type=int,help='Number of eigenvalues to calculate [ 100 ]')
parser.add_argument('-w',dest='writeFN',default='timescales',help='Base filename to save data and plot to. Will add <num_vals>.npy and <num_vals>.pdf respectively [ timescales ]')
parser.add_argument('-l',dest='lag',default=1,type=int,help='Lag time to build the MSM [ 1 ]')
parser.add_argument('-d',dest='divisor',default=1.,type=float,help='Divisor to convert your data to some time coordinate')
parser.add_argument('-u',dest='units',default='frames',help='Units to write on the y-axis [ frames ]')
parser.add_argument('--x-label',dest='x_lbl',help='x-axis label [ OPTIONAL ]')
parser.add_argument('--title',dest='title',help='Title for plot [ OPTIONAL ]')
parser.add_argument('--y-lim',dest='y_lim',nargs=2,type=float,help='y-limits if you want to choose them',default=None)
parser.add_argument('--print-not-calc',dest='print_not_calc',default=False,action='store_true',help='Pass this flag if you want to print "Not Calculated in the empty space at the bottom of the graph')
parser.add_argument('--color', dest='color', default='blue', help='color of lines.')
parser.add_argument('--no-y', dest='no_yaxis', default=False, action='store_true', help='Pass this flag to remove y-ticks and y-label')
options = parser.parse_args()
 
import matplotlib
matplotlib.use('pdf')
import numpy as np
from msmbuilder import msm_analysis
from matplotlib.pyplot import *
from scipy.io import mmread
import os, sys, re

matplotlib.rcParams['font.size']=22

#import warnings
#warnings.filterwarnings("ignore",category=DeprecationWarning)
if os.path.exists( options.writeFN +'%d.npy'%options.num_vals ):
   print "Found %s, and using these values rather than recalculating them." % ( options.writeFN+'%d.npy'%options.num_vals )
 
   Vals = np.load( options.writeFN+'%d.npy'%options.num_vals )
else:
   T = mmread( options.T_FN )

   Vals,Vecs = msm_analysis.get_eigenvectors( T, options.num_vals+1 )

   Vals = Vals.real[1:]

   np.save( options.writeFN + '%d.npy' % options.num_vals, Vals )
   print "Saved values"

Vals = Vals[np.where(Vals > 0)]
figure()
subplot(132)

Taus = - options.lag / options.divisor / np.log( Vals )
hlines( Taus, 0, 1, color=options.color)
if options.y_lim == None:
   ylim([10**int(np.log10(Taus.min())),10**int(np.log10(Taus.max())+1)])
else:
   ylim( options.y_lim )

yscale('log')
xticks([])
ylabel('Timescales (%s)' % options.units )
if options.x_lbl:
   xlabel(options.x_lbl)
if options.title:
   title(options.title)

if options.print_not_calc:
   text(0.11, ylim()[0] + 10**(np.log10(ylim()[1])-np.log10(ylim()[0]))*0.05, "Not Calculated", fontsize=14)

if options.no_yaxis:
    yticks([])
    ylabel('')

print "Saving plot"
savefig( options.writeFN + '%d.pdf' % options.num_vals )
