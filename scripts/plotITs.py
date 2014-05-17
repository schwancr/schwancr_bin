#!/usr/bin/env python
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument('-d',dest='dataFNs',action='append',help='Data files. Pass more than one')
parser.add_argument('-l',dest='labels',nargs='+',action='append',help='Labels for the legend. There should be one -l for every -d')
parser.add_argument('-o',dest='outFN',default='ITs.pdf',help='Output filename')
parser.add_argument('--divisor',dest='divisor',default=1,type=float,help='Float to divide the data by to convert to a unit you want to plot in')
parser.add_argument('--frame-label',dest='frame_label',help='Label to add to the plot if you want to label it. Will appear in the upper left corner')
parser.add_argument('--units',dest='units',default='frames',help='Units to print on the axes')
parser.add_argument('--slowest',dest='slowest',default=False,action='store_true',help='Pass this flag to only plot the slowest timescale')
parser.add_argument('--y-lim',dest='y_lim',default=None,type=float,nargs=2,help='Two values corresponding to the limits for the y axis')
parser.add_argument('--x-lim',dest='x_lim',default=None,type=float,nargs=2,help='Two values corresponding to the limits for the x axis')
parser.add_argument('--font-size',dest='font_size',default=22, type=int,help='Font size to use.')
parser.add_argument('--top-N',dest='top_N',type=int,default=-1,help='Only plot the top N timescales (-1 means plot all)')
parser.add_argument('--plot-line',dest='plot_line',default=False,action='store_true',help='Pass this flag if you want a line plot instead of a scatter plot')
parser.add_argument('--lin-scale',dest='lin_scale',default=False, action='store_true', help='Pass this flag if you want the y-scaled linearly. Default is log-scale.')
options = parser.parse_args()
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from pyschwancr.custom_tickers import SchwancrFormatter
import re,sys,os

if options.font_size != None:
   matplotlib.rcParams['font.size'] = options.font_size
matplotlib.rcParams['axes.formatter.limits'] = (-1,4)

man_colors = True
ManualColors = [ 'blue', 'red', 'purple','orange','yellow' ]

#myTicker2 = SchwancrFormatter(useMathText=True,useOffset=False,precision=2)
myTicker2 = SchwancrFormatter(useMathText=True,useOffset=False,precision=0)

#matplotlib.rc('text',usetex=True)

colors = matplotlib.colors.Normalize()
colors.autoscale( [-0.15] + range( len(options.dataFNs) ) + [ len(options.dataFNs) -0.85]  )


markerSize = 80
Markers = [ (i,j,0) for j in range(3) for i in range(4+j,6) ]
myCM = matplotlib.cm.jet
colorMap = matplotlib.cm.ScalarMappable(norm=colors,cmap=myCM)

figure(figsize=(8,6))
h=0.18
axes( (h,h,1-2*h,1-2*h) )

gca().xaxis.set_major_formatter( myTicker2 )
gca().yaxis.set_major_formatter( myTicker2 )

if len(options.dataFNs) <= len( ManualColors ):
   man_colors = True

if len(options.labels)!=len(options.dataFNs):
   print "Not the same number of labels as data filenames, will use the filenames in the legend."
   labels=[ '%s'%( '.'.join(fn.split('.')[:-1]) ) for fn in options.dataFNs ] # remove the extension from the filename.
else:
   labels=[ '%s'%(' '.join(lbl) ) for lbl in options.labels ]

for i,fn in enumerate(options.dataFNs):
   try: dat=np.loadtxt(fn)
   except: 
      try: dat=np.load(fn)
      except: 
         print "cannot load with np.load or np.loadtxt"
         exit()

   dat[:, 1][np.where(dat[:, 1] <= 0)] = 1E-10

   if options.slowest: # Need to subsample data by the number of eig values
      numEig = float(len(dat)) / len(np.unique(dat[:,0]))
      if numEig-int(numEig)!=0:
         print numEig, len(dat), len(np.unique(dat[:,0]))
      dat=dat[::int(numEig)]
   dat/=options.divisor
   if options.top_N != -1:
      numEig = float(len(dat)) / len(np.unique(dat[:,0]))
      if numEig-int(numEig)!=0:
         print numEig, len(dat), len(np.unique(dat[:,0]))
      if numEig <= options.top_N:
         pass
      else:
         keep_ind = np.array([ ind for ind in xrange(len(dat)) if (ind % numEig ) < options.top_N ]).astype(int)
         dat=dat[keep_ind]
   
  # if options.slowest: 
  #    plot(dat[:,0],dat[:,1],color=colorMap.to_rgba(i))
   if options.plot_line:
      line_width=3
      a=0.75
      num_vals = len( np.where(dat[:,0] == dat[0,0] )[0] )
      if man_colors:
         #plot(dat[::num_vals,0],dat[::num_vals,1],lw=line_width,color=ManualColors[i],label=labels[i], alpha=a )
         plot(dat[::num_vals,0],dat[::num_vals,1],lw=line_width,color=ManualColors[i], alpha=a )
      else:
         #plot(dat[::num_vals,0],dat[::num_vals,1],lw=line_width,color=colorMap.to_rgba(i),label=labels[i], alpha=a )
         plot(dat[::num_vals,0],dat[::num_vals,1],lw=line_width, color=colorMap.to_rgba(i), alpha=a)
      for start_ind in xrange(1, num_vals):
         if man_colors:
            plot(dat[start_ind::num_vals,0],dat[start_ind::num_vals,1],lw=line_width,color=ManualColors[i], alpha=a )
         else:
            plot(dat[start_ind::num_vals,0],dat[start_ind::num_vals,1],lw=line_width,color=colorMap.to_rgba(i), alpha=a )
      if man_colors:   
         scatter(dat[:,0], dat[:,1],s=markerSize,c=ManualColors[i],marker=Markers[i%len(Markers)],norm=colors,label=labels[i],edgecolor='none',cmap=myCM)
      else:
         scatter(dat[:,0], dat[:,1],s=markerSize,c=colorMap.to_rgba(i),marker=Markers[i%len(Markers)],norm=colors,label=labels[i],edgecolor='none',cmap=myCM)
        
   else:
      scatter(dat[:,0], dat[:,1],s=markerSize,c=np.array([i]*len(dat)),marker=Markers[i%len(Markers)],norm=colors,label=labels[i],edgecolor='none',cmap=myCM)
   
lg=legend(loc=4,prop={'size':options.font_size})  # lower right
#lg=legend(loc=2,prop={'size':options.font_size}) # upper left
#lg.get_frame().set_facecolor('#cdc9c9')
ax = gca()
x0,y0 = ax.xaxis.offsetText.get_position()
#ax.xaxis.offsetText.set_x(x0*1.2)

for i, h in enumerate( lg.legendHandles ):
   if man_colors:
      h.set_color( ManualColors[i] )
   else:
      h.set_color( colorMap.to_rgba(i) )
   #h.set_edgecolor( 'black' )
xlabel('Lag Time (%s)' % options.units, fontsize=options.font_size+5 )
ylabel('Timescale (%s)' % options.units, fontsize=options.font_size+5 )

if options.y_lim:
   ylim( options.y_lim )
else:
   ylim([1,ylim()[1]])

if options.x_lim:
   xlim( options.x_lim )
   xticks( np.linspace( options.x_lim[0], options.x_lim[1], 6) )
else:
   xlim([0,dat[:,0].max()])

if options.frame_label:
   xloc = (xlim()[1] - xlim()[0]) * 0.025 + xlim()[0]
   yloc = (ylim()[1] - ylim()[0]) * 0.6 + ylim()[0]

   text( xloc, yloc, options.frame_label )


if not options.lin_scale:
    yscale('log')

#gca().set_axis_bgcolor('#cdc9c9')
savefig(options.outFN, dpi=200)
