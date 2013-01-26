#!/usr/bin/env python

"""
This script generates a 2D plot of a pathway given two order parameters x and y.
"""
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-x',dest='xFN',help='This is the x-coordinates for all the states')
parser.add_option('-y',dest='yFN',help='This is the y-coordinates for all the states')
parser.add_option('--xc',dest='xcol',type=int,help='The column in the x-data to use')
parser.add_option('--yc',dest='ycol',type=int,help='The column in the y-data to use')
parser.add_option('-d',dest='dataFN',help='Data to plot for all the states.')
parser.add_option('--pos-or-neg',dest='pos_or_neg',help='Pass this flag if you want to change the data to 1 and 0 corresponding to whether it is positive or negative. Like if you are plotting an eigenvector of a transition probability matrix.', default=False, action='store_true')
options,args = parser.parse_args()

import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *
from numpy import *
from pyschwancr.dataIO import readData as rd
from pyschwancr.dataIO import writeData as wd

X = rd( options.xFN ).real
Y = rd( options.yFN ).real
data = rd( options.dataFN ).real
print X.shape, Y.shape, data.shape
if len( X.shape ) > 1:
   print "X formatted strangely... Need to choose a column"
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

X=X[np.where(data!=-1)]
Y=Y[np.where(data!=-1)]
data=data[np.where(data!=-1)]

figure()

eps = 0.001 * np.abs( data ).max()
print eps

if options.pos_or_neg:
   pos_ind = np.where( data < -eps )
   neg_ind = np.where( data > eps )
   if len( pos_ind[0] ) > len( neg_ind[0] ):
      scatter( X[pos_ind], Y[pos_ind], color='black')
      scatter( X[neg_ind], Y[neg_ind], color='red')
   else:
      scatter( X[neg_ind], Y[neg_ind], color='red')
      scatter( X[pos_ind], Y[pos_ind], color='black')
else:
   scatter( X, Y, c=data, vmin=data.min(),vmax=data.max())

   colorbar()

# Reformat the file names so that they look nice
name = '.'.join( options.dataFN.split('/')[-1].split('.')[:-1] )
Xname = '.'.join( options.xFN.split('/')[-1].split('.')[:-1] )
Yname = '.'.join( options.yFN.split('/')[-1].split('.')[:-1] )

xlabel( Xname )
ylabel( Yname )
title( name )
#xlim([0,1])
#ylim([0,1])
outFN = '%s_%sVS%s.pdf' % (name,Yname,Xname)

savefig( outFN )


