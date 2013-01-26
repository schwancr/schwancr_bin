#!/usr/bin/env python

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-x','--x-file',dest='xFN',help='numpy array saved as .npy to plot on the x-axis')
parser.add_option('-y','--y-file',dest='yFN',help='numpy array saved as .npy to plot on the y-axis')
parser.add_option('-o','--output-PDF',dest='outFN',help='output filename (pdf)\n Default: YvsX.pdf')
#parser.add_option('-l','--use-log-scale',dest='useLogScale',default=False,action='store_true',help='Use a log scale for the color scale')

options, args = parser.parse_args()

from numpy import *
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from numpy import log10 as l10
	
X = load(options.xFN)
Y = load(options.yFN)

if len( X.shape ) > 1:
	print "X formatted strangely... Using X[:,0]"
	X = X[:,0]

if len( Y.shape ) > 1:
	print "Y formatted strangely... Using Y[:,0]"
	Y = Y[:,0]

nx = len(unique1d(X))
ny = len(unique1d(Y))

print nx,ny

dst = zeros((nx,ny))
# dst is going to be the 2D array that I plot with contour

print "Data completely loaded..."
print "Calculating the distribution"

zpd = array([ X*(nx-1), Y*(ny-1) ]).T

for pair in zpd:
	dst[ tuple(pair) ] += 1

x = linspace(0,1,nx)
y = linspace(0,1,ny)

print "Interpolating data..."
fit = interpolate.RectBivariateSpline(x,y,dst,kx=1,ky=1,s=0)

xi = linspace(0,1,nx*10)
yi = linspace(0,1,ny*10)

smoothDat = fit(xi,yi)
data = l10( smoothDat )

data[ where( data==-inf ) ] = data[ where( data != -inf ) ].min()

savetxt( "%s_vs_%s.dat" % (options.xFN[:-4], options.yFN[:-4]), data.T )

# This makes it log scale!
figure()
imshow(data.T,origin='lower',vmin=0,vmax=5,extent=[0,1,0,1])
# Y-axis is actually the "primary" axis...

xlabel(options.xFN[:-4])
ylabel(options.yFN[:-4])
colorbar().set_label('Log10(counts)')
if options.outFN:
	if options.outFN[:-4] == '.pdf':
		filename = options.outFN
	else:
		filename = options.outFN+'.pdf'
	pp = PdfPages(filename)
else:
	pp = PdfPages('%sVs%s.pdf'%(options.yFN[:-4],options.xFN[:-4]) )

pp.savefig()

pp.close()
