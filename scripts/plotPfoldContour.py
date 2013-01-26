#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_option('-p',dest='proj_FN',default='../../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../../ProjectInfo.h5 ]')
parser.add_option('-x',dest='x_FN',help='Data for x-axis (all points)')
parser.add_option('-y',dest='y_FN',help='Data for y-axis (all points)')
parser.add_option('--xl',dest='x_lbl',help='X-axis label')
parser.add_option('--yl',dest='y_lbl',help='Y-axis label')
parser.add_option('-f',dest='fc_FN',default='./FCs.dat',help='Forward Committors for the states in the MSM')
parser.add_option('-o',dest='out_FN',default='FCs_Contour.pdf',help='Output filename')


options, args = parser.parse_args()
 
import numpy as np
from msmbuilder import Serializer
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from pyschwancr import dataIO, msmTools
import os, sys, re
 
print "Loading Data"

Xdat = dataIO.readData( options.x_FN )[:,0]
Ydat = dataIO.readData( options.y_FN )[:,0]

Proj = Serializer.Serializer.LoadFromHDF( options.proj_FN )
Ass = Serializer.LoadData( options.ass_FN )

FCperState = dataIO.readData( options.fc_FN )
print Xdat.shape, Ydat.shape
Xdat = msmTools.reshapeRawData( Xdat, Proj )
Ydat = msmTools.reshapeRawData( Ydat, Proj )

keptAssInd = np.where( Ass != -1 )

Xdat = Xdat[ keptAssInd ]
Ydat = Ydat[ keptAssInd ]
FCdat = FCperState[ Ass[ keptAssInd ] ]
# The above is all the data, but now we need to average the FC's for conformations at the same point in the projection.

if Xdat.shape[0] != Ydat.shape[0] or Ydat.shape[0] != FCdat.shape[0]:
   print "Data is not the same shape! %d, %d, %d" % ( Xdat.shape[0], Ydat.shape[0], FCdat.shape[0] )

nx = len( np.unique( Xdat ) )
ny = len( np.unique( Ydat ) )
Z = np.zeros( (nx, ny) )
Ns = np.zeros( (nx, ny) )

zpd = np.array([ Xdat * (nx-1), Ydat * (ny-1) ]).T
Nrows = len(zpd)
print "Calculating the surface"
for i in xrange( Nrows ):
   #print "\rWorking on point %d / %d" % ( i, Nrows ),
   pair = tuple( zpd[i] )
   Z[ pair ] += FCdat[ i ]
   Ns[ pair ] += 1.
print "here"  
Z[ np.where( Ns != 0 ) ] /= Ns[ np.where( Ns != 0 ) ]


#Z = zip( Xdat, Ydat )
#uniqZ = np.unique( Z )
#FCs = np.zeros( len( uniqZ ) )
#
#totRows = len(uniqZ)
#
#for i in xrange( 100 ):# len( FCdat ) ):
#   print "\rWorking on row %d / %d " % (i, totRows ),
#   ind = np.where( ( Xdat == uniqZ[i,0] ) * ( Ydat  == uniqZ[i,1] ) )
#   FCs[ i ] = FCdat[ ind ].mean()
#   
#
figure()
np.save( options.out_FN[:-4] + '.npy', Z )
contour( Z.T,10, extent = [0,1,0,1] )
contourf( Z.T,10, extent = [0,1,0,1] )
xlabel( options.x_lbl )
ylabel( options.y_lbl )
colorbar( boundaries = [0,1] )
savefig( options.out_FN )
print "Done."
