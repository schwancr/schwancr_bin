#!/usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-d',dest='data_FNs',nargs='+',help='Data filenames to find implied timescale data.')
parser.add_argument('-k',dest='num_states',nargs='+',type=int,help='Number of states associated with each data filename')
parser.add_argument('-o',dest='out_FN',default='ITs2D.pdf',help='Output filename to save the plot.')
parser.add_argument('--divisor',dest='divisor',default=1,type=float,help='Divisor to divide the lag times and timescales by. For example, if each frame corresponded to 200ps in your simulation, if you passed --divisor 5, then the conversion would turn frames into ns')
parser.add_argument('--units',dest='units',default='frames',help='Units to include in the x-axis description for the lag time and in the timescale description')
parser.add_argument('-t',dest='timescale',default=1,type=int,help='Which timescale to plot. [ NOTE: 1-indexed ]')
parser.add_argument('--range',dest='range',nargs=2,type=float,help='Range to use for the colormap scale',default=[0,1000])
parser.add_argument('--font-size=',dest='font_size',type=int,default=16,help='Font size to use for the text in the figure.')
parser.add_argument('--label',dest='label',help='Label to put in the upper right corner of the axes.')
parser.add_argument('--c-ticks',dest='cticks',default=6,type=int,help='Number of ticks to place on the colorbar')
parser.add_argument('--levels',dest='levels',nargs='+',type=float,help='Levels to include in contour plot')
parser.add_argument('--contour',dest='contour',default=False,action='store_true',help='Pass this flag to make a contour plot instead of color plot')
args = parser.parse_args()

import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['font.size'] = args.font_size
matplotlib.rcParams['axes.formatter.limits'] = -1,1
from matplotlib.pyplot import *
import numpy as np
from pyschwancr.custom_tickers import SchwancrFormatter
import re, sys, os

def custom_formatter( number ):
   if number < 1000:
      return str(number)
   return str( number / 1000 ) + 'k'

#myTicker=matplotlib.ticker.ScalarFormatter(useMathText=True,useOffset=False)
myTicker2 = SchwancrFormatter(useMathText=True,useOffset=False,precision=2)
myTicker1 = SchwancrFormatter(useMathText=True,useOffset=False,precision=1)

def tick2_wrap( number ):
   return '$%s$' % myTicker2.format_data( number )

def tick1_wrap( number ):
   return '$%s$' % myTicker1.format_data( number )
Datas = np.array([ np.loadtxt( fn ) / args.divisor for fn in args.data_FNs ])

UniqLags = np.unique( Datas[:,:,0] )

NumTimescales = float( Datas.shape[1] ) / len( UniqLags )

if (NumTimescales - int( NumTimescales ) ) > 0.001:
   print "Something went wrong..."
   exit()

timescale = args.timescale - 1
if not 0 <= timescale < NumTimescales:
   print "Not a valid timescale (%d). There are only %d total timescales in the data." % (timescale, NumTimescales)

Datas = Datas[:,timescale::NumTimescales,1] 

figure(figsize=(8,6))
h=0.18
ax=axes( (h,h,1-2*h,1-2*h) )
if not args.contour:
   im = ax.imshow( Datas, origin='lower',aspect='auto',extent=[UniqLags.min(),UniqLags.max(),-1,Datas.shape[0]],interpolation='nearest',vmin=args.range[0]/args.divisor,vmax=args.range[1]/args.divisor )
else:
   if args.levels != None:
      V = np.array( args.levels ) / args.divisor
   else:
      V = np.linspace( args.range[0], args.range[1],10) / args.divisor
      V=V.round( -2 )[5:]
   print V
   im = ax.contour( Datas,V, origin='lower',extent=[UniqLags.min(),UniqLags.max(),-1,Datas.shape[0]],vmin=args.range[0]/args.divisor,vmax=args.range[1]/args.divisor )
   clabel( im, fmt=tick2_wrap,fontsize=14 )
print Datas.max()

np.savetxt('data.dat', Datas)

if len(UniqLags) > 5:
   subx=len(UniqLags)/3
else:
   subx=1

#yticks(np.arange( Datas.shape[0] ), [ '$%s$'%myTicker.format_data(n) for n in args.num_states ])
#xticks(np.arange( Datas.shape[1] )[::subx],[ '$%s$'%myTicker.format_data(n) for n in UniqLags[::subx] / args.divisor ])
ax.xaxis.set_major_formatter( myTicker2 )
#ax.yaxis.set_major_formatter( myTicker )

ax.set_yticks(np.arange( Datas.shape[0] ))

#ax.set_yticklabels( [ str(n) for n in args.num_states ] ) # [ str(n) for n in args.num_states ])
#ax.set_yticklabels( [r'$%s\times10^%d$' % (a,int(b)) for (a,b) in [ ( '%.1E'%(n,) ).split('E')  for n in args.num_states ] ][::2] ) # [ str(n) for n in args.num_states ])
ax.set_xticks(UniqLags[::subx])#np.arange( Datas.shape[1] )[::subx] )

#ax.set_xticklabels( [ '$%s$'%myTicker.format_data(l) for l in UniqLags[::subx] ],rotation=30 )#[ str(n) for n in UniqLags[::subx] / args.divisor ])
ax.set_yticklabels( [ custom_formatter( n ) for n in args.num_states ] )
ax.set_ylabel('Number of States',fontsize=args.font_size+5)#r'$\textnormal{Number of States}$')
ax.set_xlabel('Lag Time (%s)' % args.units,fontsize=args.font_size+5) #r'$\textnormal{Lag Time (%s)}$'%args.units)

x0,y0 = ax.xaxis.offsetText.get_position()
ax.xaxis.offsetText.set_x(x0*1.2)


oldXlim=xlim()

if not args.contour:
   colorbar_ticks = np.arange( args.range[0], args.range[1] + 100, 100 ) / args.divisor
   colorbar_ticks = colorbar_ticks.astype(float)
#if len(colorbar_ticks) > 20:
#   stride = int( len(colorbar_ticks) / 6)
#   colorbar_ticks = colorbar_ticks[::stride]
#elif 5 < len(colorbar_ticks ) <= 20:
#   colorbar_ticks = colorbar_ticks[::2]

   stride = len(colorbar_ticks) / args.cticks
   colorbar_ticks = colorbar_ticks[::stride]
   c = colorbar(im,ticks=colorbar_ticks, format=myTicker1, fraction=0.05)

   c.set_label('Timescale (%s)' % args.units,fontsize=args.font_size+5 )#r'$\textnormal{Timescale (%s)}$'%args.units)
   c.ax.yaxis.offsetText.set_x( 2 )
#c.ax.yaxis.offsetText.set_y( 0 ) # This doesn't work for some reason... Because it is either bottom or top...
# We can control the "pad" with this:
   c.ax.yaxis.OFFSETTEXTPAD += 15
#c.ax.yaxis.set_ticklabels([ '%.1f'%float(lbl.get_text().split('{')[1].split('}')[0]) for lbl in c.ax.yaxis.get_ticklabels() ])
   c.update_ticks()
   print oldXlim, xlim()
   xlim(oldXlim)

if args.label:
   xloc=xlim()[0] + (xlim()[1]-xlim()[0]) *0.05
   yloc=ylim()[0] + (ylim()[1]-ylim()[0]) *0.9
   text(xloc,yloc,args.label,  backgroundcolor='white') # White background in text object
   #text(xloc,yloc,args.label) # No background in text object

savefig( args.out_FN, dpi=200 )
