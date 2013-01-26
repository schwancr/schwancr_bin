#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-x','--x-file',dest='xFN',help='numpy array saved as .npy to plot on the x-axis')
parser.add_option('-y','--y-file',dest='yFN',help='numpy array saved as .npy to plot on the y-axis')
parser.add_option('--xc','--x-column',dest='xcol',type=int,help='Column in the x-data to use')
parser.add_option('--yc','--y-column',dest='ycol',type=int,help='Column in the y-data to use')
parser.add_option('-o','--output-PDF',dest='outFN',help='output filename (pdf)\n Default: YvsX.pdf')
#parser.add_option('-l','--use-log-scale',dest='useLogScale',default=False,action='store_true',help='Use a log scale for the color scale')

parser.add_option('-w','--weights',dest='weights',help='File with the population size for each of the clusters.')
options, args = parser.parse_args()

from numpy import *
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from numpy import log10 as l10

pops = loadtxt(options.weights)

if options.xFN[-3:] == 'npy':
	X = load(options.xFN)
else:
	try: X = loadtxt(options.xFN)
	except: print "Cannot read file %s."% options.xFN; exit()
if options.yFN[-3:] == 'npy':
	Y = load(options.yFN)
else:
	try: Y = loadtxt(options.yFN)
	except: print "Cannot read file %s."% options.yFN; exit()

if len( X.shape ) > 1:
	print "X formatted strangely... Using X[:,0]"
	if not options.xcol:
		print "Need to input the column to use as --xc!"
		print "Continuing using column 0. Careful this could give wacky results!"
		X = X[:,0]
	else:
		X = X[:,options.xcol]
	
if len( Y.shape ) > 1:
	print "Y formatted strangely... Need to choose a column"
	if not options.ycol:
		print "Need to input the column to use as --yc!"
		print "Continuing using column 0. Careful this could give wacky results!"
		Y = Y[:,0]
	else:
		Y = Y[:,options.ycol]

nx = 100
ny = 100

print nx,ny

dst = zeros((nx,ny))
# dst is going to be the 2D array that I plot with contour

print "Data completely loaded..."
print "Calculating the distribution"

zpd = array([ X*(nx-1), Y*(ny-1), pops ]).T

for (a,b,p) in zpd:
	dst[ (a,b) ] += p

x = linspace(0,1,nx)
y = linspace(0,1,ny)

print "Interpolating data..."
fit = interpolate.RectBivariateSpline(x,y,dst,kx=1,ky=1,s=0)

xi = linspace(0,1,nx*10)
yi = linspace(0,1,ny*10)

smoothDat = fit(xi,yi)
print "Making the plot..."

data = l10( smoothDat )

data[ where( data==-inf ) ] = 0
# This makes it log scale!
figure()
imshow(data.T,origin='lower',vmin=0,vmax=5,extent=[0,1,0,1])
# Y-axis is actually the "primary" axis...

xlabel(options.xFN[:-4])
ylabel(options.yFN[:-4])
colorbar().set_label('Log10(counts)')
if options.outFN:
	if options.outFN[-4:] == '.pdf':
		filename = options.outFN
		title(options.outFN[:-4])
	else:
		filename = options.outFN+'.pdf'
		title(options.outFN)
	pp = PdfPages(filename)
else:
	title('%sVs%s'%(options.yFN[:-4],options.xFN[:-4]) )
	pp = PdfPages('%sVs%s.pdf'%(options.yFN[:-4],options.xFN[:-4]) )

pp.savefig()

pp.close()
