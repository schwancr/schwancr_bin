#!/usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i','--input',dest='inputFN',default='./IT.dat',help='Input file to plot (default = IT.dat)')
parser.add_argument('-o','--output',dest='outFN',default='./IT.pdf',help='Output file with plot')
parser.add_argument('-d',dest='divisor',default=1.,type=float,help='Divide everything by this number. For example, if each frame is 200ps, and you want your data to be in ns, then pass 5, and each time will be converted to ns, since 1 ns = 5 frames.')
parser.add_argument('--units',dest='units',default='frames',help='Units to place in the x and y-axes labels.')
parser.add_argument('--line',dest='line',type=float,help='Plot a horizontal line corresponding to the actual time.' )
parser.add_argument('--y-lim',dest='y_lim',type=float,nargs=2,help='y-limit to apply to the log-scale plot')
parser.add_argument('--title',dest='title',nargs='+',default=None,help='Title for the plot')
args = parser.parse_args()

import matplotlib
matplotlib.use('agg')

from matplotlib.pyplot import *
from numpy import *
import os
import re
data = loadtxt(args.inputFN)
data /= args.divisor
name = os.getcwd().split('/')[-1]

Mgrp = re.search('Data_(.*)',name)
if Mgrp:
	name = Mgrp.group(1)
else:
   Mgrp = re.search('Data_(.*)/',args.inputFN)
   if Mgrp:
      name = Mgrp.group(1)
figure()
subplot(211)
scatter(data[:,0],data[:,1])
yscale('log')
xlim([ max( 0,xlim()[0] ), data[:,0].max()+1/args.divisor])
if args.line:
	hlines( args.line / args.divisor, xlim()[0], xlim()[1], color='black',linestyle='--')
#ylim([1,ylim()[1]])
ylabel('Timescale (%s)'%(args.units,))
if args.y_lim:
   ylim( args.y_lim )

subplot(212)

Neigs = len( where( data[:,0] == data[:,0][0] )[0] )
for i in range( Neigs ):
	plot( data[:,0][i::Neigs], data[:,1][i::Neigs] )
xlabel('Lag Time (%s)'%(args.units,))
ylabel('Timescale (%s)'%(args.units,))
xlim([ max( 0,xlim()[0] ), data[:,0].max()+1/args.divisor])
if args.line:
	hlines( args.line / args.divisor, xlim()[0], xlim()[1], color='black',linestyle='--')
#ylim([1,ylim()[1]])
if args.y_lim:
    ylim( args.y_lim )

if args.title==None:
	suptitle( name )
else:
	suptitle( ' '.join( args.title ) )

savefig(args.outFN)

