#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-x','--x-file',dest='xFN',help='numpy array saved as .npy to plot on the x-axis')
parser.add_option('-y','--y-file',dest='yFN',help='numpy array saved as .npy to plot on the y-axis')
parser.add_option('-o','--output-PDF',dest='outFN',help='output filename (pdf)\n Default: YvsX.pdf')
parser.add_option('--xl',dest='x_lbl',help='Label for x-axis')
parser.add_option('--yl',dest='y_lbl',help='Label for y-axis')
parser.add_option('-t',dest='title',help='Title for the plot')
parser.add_option('--label',dest='label',help='Label to place in plot')
parser.add_option('--norm-x',dest='normX',help='Pass this to control the length used. For example if I want to plot/interpolate from 0.5 to 1, then I would pass --norm-x 0.5,1')
parser.add_option('--norm-y',dest='normY',help='Pass this to control the length used. For example if I want to plot/interpolate from 0.5 to 1, then I would pass --norm-y 0.5,1')
parser.add_option('--orders',dest='orders',type=int,default=5,help='Number of orders of magnitude to include in the color scale [ 5 ]')
parser.add_option('--weights',dest='weights',help='Weights to use in calculating the destination, if not passed in, a vector of ones is used')
parser.add_option('--font-size',type=int, dest='font_size',help='Font size to use for the xticks and labels and everything... [ None ]')
#parser.add_option('--log-x',dest='logX',default=False,action='store_true',help='Pass this flag to use a log scale on the X axis')
#parser.add_option('--log-y',dest='logY',default=False,action='store_true',help='Pass this flag to use a log scale on the Y axis')
#parser.add_option('-l','--use-log-scale',dest='useLogScale',default=False,action='store_true',help='Use a log scale for the color scale')

options, args = parser.parse_args()

from numpy import *
import matplotlib
from matplotlib.pyplot import *
#from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from pyschwancr import dataIO
from numpy import log10 as l10
import re

#matplotlib.rc('text',usetex=True)
	
if options.font_size != None:
	matplotlib.rcParams['font.size'] = options.font_size

print "IF THE INPUT DATA IS COMPLEX, THE PLOTS WILL ONLY BE THE REAL PART!!!"
X = dataIO.readData( options.xFN ).real
Y = dataIO.readData( options.yFN ).real

if options.xFN[-3:] == '.h5':
    X = X[np.where(X!=-1)]
if options.yFN[-3:] == '.h5':
    Y = Y[np.where(Y!=-1)]

if len( X.shape ) > 1:
	if X.shape[1] == 1:
		print "X formatted strangely... Using X[:,0]"
		X = X[:,0]
	else:
		X = X[:,1]
		print "X formatted strangely... Using X[:,1]"

if len( Y.shape ) > 1:
	if Y.shape[1] == 1:
		print "Y formatted strangely... Using Y[:,0]"
		Y = Y[:,0]
	else:
		Y = Y[:,1]
		print "Y formatted strangely... Using Y[:,1]"

#Nonzeros = where( ( X != 0 ) * ( Y != 0 ) )

#Y = Y[ Nonzeros ]
#X = X[ Nonzeros ]

if X.shape[0] != Y.shape[0]:
	print "Shape of data is not the same!"
	XoverY = X.shape[0] / float( Y.shape[0] )
	YoverX = Y.shape[0] / float( X.shape[0] )
	if XoverY == int( XoverY ): # then subsample X by XoverY
		X = X[::int(XoverY)]
		print "Subsampling X by %d" % int(XoverY)
	elif YoverX == int( YoverX ):
		Y = Y[::int(YoverX)]
		print "Subsampling Y by %d" % int(YoverX)
	else:
		print "Shapes are not compatible (one is not a multiple of the other)"

print X.shape, Y.shape

# remove nans

ind_isNotNan = np.where( ( 1-np.isnan(X) ) * ( 1-np.isnan(Y) ) )

X = X[ ind_isNotNan ]
Y = Y[ ind_isNotNan ]
		


#if options.logX:
#	X = l10(X)
#if options.logY:
#	Y = l10(Y)
# ^^ This doesn't work because there may be negative numbers...
if options.weights:
	pops = dataIO.readData( options.weights )
else:
	print "Using all ones for the weights of each piece of data"
	pops = ones( X.shape[0] )

if re.search( 'RMSD', options.xFN ):
	nx = min( 100, len(unique(X)) )
elif re.search( 'Rg', options.xFN ):
	nx = min( 100, len(unique(X)) )
else:
	nx = len(unique(X))

if re.search( 'RMSD', options.yFN ):
	ny = min( 100, len(unique(Y)) )
elif re.search( 'Rg', options.yFN ):
	ny = min( 100, len(unique(Y)) )
else:
	ny = len(unique(Y))

if nx > 100:
	nx = 100
if ny > 100:
	ny = 100

if options.normX:
	a,b = options.normX.split(',')
	minX = float(a)
	maxX = float(b)
else:
	minX = int( 2 * X.min()  ) / 2.
	maxX = int( 2 * X.max() + 1 ) / 2.

if options.normY:
	a,b = options.normY.split(',')
	minY = float(a)
	maxY = float(b)
else:
	minY = int( 2 * Y.min() ) / 2.
	maxY = int( 2 * Y.max() + 1 ) / 2.
#print nx,ny
dst = zeros((nx,ny))
# dst is going to be the 2D array that I plot with contour

print "Data completely loaded..."
print "Calculating the distribution"

zpd = array([ ( X - minX ) / ( maxX - minX ) * (nx-1), ( Y - minY ) / ( maxY - minY ) * (ny-1) ]).T

for i,pair in enumerate(zpd):
   
	if 0 < pair[0] > nx - 1:
		print "Outside the boundary!"
		continue
	if 0 < pair[1] > ny - 1:
		print "Outside the boundary!"
		continue
	try: dst[ tuple(pair) ] += pops[i]
	except: 
		print "Something didn't work adding pair (%s) to the distribution" % str(pair)
		continue
x = linspace(minX,maxX,nx)
y = linspace(minY,maxY,ny)

print "Interpolating data..."
fit = interpolate.RectBivariateSpline(x,y,dst,kx=1,ky=1,s=0)

print minX, maxX
print minY, maxY
#minX=0 # need to do this so that both origins are the same, since only one value is printed
#minY=0
xi = linspace(minX,maxX,nx*10)
yi = linspace(minY,maxY,ny*10)

smoothDat = fit(xi,yi)
data = - l10( smoothDat.max() ) + l10( smoothDat )
data[ where( data==-inf ) ] = data[ where( data != -inf ) ].min()

savetxt( "%s_vs_%s.dat" % ( options.xFN[:-4].split('/')[-1], options.yFN[:-4].split('/')[-1] ), data.T )

# This makes it log scale!
figure(figsize=(8,6))
h=0.18
axes((h,h,1-2*h,1-2*h))
imshow(data.T,origin='lower',vmin= - options.orders,vmax=0,aspect='auto',cmap='hot_r')
# Y-axis is actually the "primary" axis...
print minX, maxX, minY, maxY
#xticks( xticks()[0][1:-1], [ round(xi[i],2) for i in xticks()[0][1:-1] ] )
#yticks( yticks()[0][1:-1], [ round(yi[i],2) for i in yticks()[0][1:-1] ] )
xLbl = linspace( minX , maxX, (maxX - minX) * 2 + 1 )

yLbl = linspace( minY, maxY, ( maxY - minY ) * 2 + 1 )

if len(yLbl) == 1:
	yLbl = array( [ minY, maxY ] )
if len(xLbl) == 1:
	xLbl = array( [ minX, maxX ] )
#xLbl = arange( 2 * int( maxX + 1 ) )[:-1] / 2.
#yLbl = arange( 2 * int( maxY + 1 ) )[:-1] / 2.
print xLbl, yLbl
xLoc = array([ where( abs(xi - i) == abs(xi - i).min() )[0][0] for i in xLbl ])
yLoc = array([ where( abs(yi - i) == abs(yi - i).min() )[0][0] for i in yLbl ])
print xi[ xLoc ], yi[ yLoc ]

numTicks = 8

if len( xLbl ) > numTicks:
	Xstride = int( len( xLoc ) / float( numTicks ) + 0.5 )
	xLoc = xLoc[::Xstride]
	xLbl = xLbl[::Xstride]

if len( yLbl ) > numTicks:
	Ystride = int( len( yLoc ) / float( numTicks ) + 0.5 )
	yLoc = yLoc[::Ystride]
	yLbl = yLbl[::Ystride]
#xticks( xLoc, [r'$%s$'%str(s) for s in xLbl], fontsize=options.font_size )
#yticks( yLoc, [r'$%s$'%str(s) for s in yLbl], fontsize=options.font_size )

xticks( xLoc[1:], xLbl[1:], fontsize=rcParams['font.size'] ) # Skip the zero, since it overlaps with y's zero
yticks( yLoc, yLbl, fontsize=rcParams['font.size'] )

if options.x_lbl:
	xlabel( options.x_lbl,fontsize=rcParams['font.size']+5 )
else:
	xlabel( options.xFN[:-4], fontsize=rcParams['font.size']+5)

if options.y_lbl:
	ylabel( options.y_lbl, fontsize=rcParams['font.size']+5 )
else:
	ylabel( options.yFN[:-4], fontsize=rcParams['font.size']+5)

if options.title:
	title( options.title, fontsize=rcParams['font.size']+5 )

cbar = colorbar( ticks = -arange( 6 ) )

cbar.set_ticklabels( ['$10^{%d}$'%i for i in -arange(6) ] )
cbar.set_label('Relative Population (Log-Scale)')

if options.label:
   xloc=xlim()[0] + (xlim()[1]-xlim()[0]) *0.05
   yloc=ylim()[0] + (ylim()[1]-ylim()[0]) *0.9
   text(xloc,yloc,options.label)

savefig( options.outFN )

#if options.outFN:
#	if options.outFN[-4:] == '.pdf':
#		filename = options.outFN
#	else:
#		filename = options.outFN+'.pdf'
#	pp = PdfPages(filename)
#else:
#	pp = PdfPages('%sVs%s.pdf'%( options.yFN[:-4].split('/')[-1], options.xFN[:-4].split('/')[-1] ) )
#
#pp.savefig()
#
#pp.close()
