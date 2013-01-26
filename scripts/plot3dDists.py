#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d',dest='data',action='append',help='Data to plot. Pass in more than one entry here')
parser.add_option('-x',dest='xVal',action='append',type=float,help='X data corresponding to each data file which contains y-values')
parser.add_option('-o',dest='outFN',default='DriftDists.pdf',help='Output filename [ DriftDists.pdf ]') 
parser.add_option('--xl',dest='x_lbl',default='Stride (frames)',help='X-axis label [ "Stride (frames)" ]')
parser.add_option('--yl',dest='y_lbl',default='Drift in Metric',help='Y-axis label [ "Drift in Metric" ]')
parser.add_option('--x-max',dest='x_max',type=float,help='X-axis max to plot [ Max( Xvals ) ]')
parser.add_option('--y-max',dest='y_max',type=float,help='Y-axis max to plot' )
options, args = parser.parse_args()
 
from numpy import *
from msmbuilder import Serializer
import matplotlib
matplotlib.use('pdf')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.pyplot import *
from pyschwancr import dataIO
import os, sys, re
 
Datas = [ dataIO.readData( fn ).flatten() for fn in options.data ]
Means = array( [ dat.mean() for dat in Datas ] )
Vars = array( [ dat.var() for dat in Datas ] ) 
Nbins = 100
Xvals = array( options.xVal )
if options.y_max:
	maxVal = options.y_max
else:
	maxVal = max( [ dat.max() for dat in Datas ] )
minVal = min( [ dat.min() for dat in Datas ] )
minVal = 0

if len(Xvals ) != len( Datas ):
	print "Need to input the same number of x-values and data files... Exiting..."
	exit()

Hists = [ hist( dat, bins= Nbins, range=[minVal,maxVal] ) for dat in Datas ]

Lefts = [ h[1][:-1] for h in Hists ]
Heights = [ h[0] / float(h[0].max()) for h in Hists ]
Modes = array( [ h[1][ h[0].argmax() ] for h in Hists ] )
xind = argsort( Xvals )
figure( figsize=(12,10) )
ax = axes( projection='3d' )
verts = []
for i in xind:
#	print Xvals[i]
	#ax.bar( Lefts[i], Heights[i], zs=Xvals[i], width= maxVal / Nbins, zdir='y',alpha=0.7,edgecolor='none',color=None,fill=False )
#	poly = PolyCollection( [ zip( Lefts[i], Heights[i] )],alpha=0.8 )
	verts.append( zip( Lefts[i], Heights[i] ) )
#	ax.add_collection3d( poly, zs = [ Xvals[i] ] * len( Lefts[i] ), zdir = 'y' )

poly = PolyCollection( verts, alpha=0.7 )

if options.x_max:
	XvalMax = options.x_max
else:
	XvalMax = Xvals.max()

xindRestrict = xind[ where( Xvals[ xind ] <= XvalMax ) ]

ax.plot( Means[xindRestrict], Xvals[xindRestrict], zs=0,zdir='z',color='black',lw=2,label='Averages' )
ax.plot( Modes[xindRestrict], Xvals[xindRestrict], zs=0,zdir='z',color='red',lw=2,label='Modes')
ax.add_collection3d( poly, zs = Xvals[xindRestrict], zdir = 'y' )
ax.legend()
ax.set_xlim3d([0,maxVal])
ax.set_ylim3d([0,XvalMax])
ax.set_zlim3d([0,3])
ax.set_xlabel(options.y_lbl)
ax.set_ylabel(options.x_lbl)
ax.set_zlabel('Relative Frequency')
savefig(options.outFN)
