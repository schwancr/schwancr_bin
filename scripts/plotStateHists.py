#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_argument('-p',dest='proj_FN',default='../ProjectInfo.h5',help='Project info from msmbuilder.')
parser.add_argument('-d',dest='data_FN',help='Data to plot for each state.')
parser.add_argument('-o',dest='out_FN',default='./StateHists.pdf',help='Output filename to save the plots to. (pdf)')
parser.add_argument('--xl',dest='x_lbl',help='X-axis label',nargs='+')

args = parser.parse_args()
 
import numpy as np
from msmbuilder import Serializer
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from pyschwancr import dataIO, msmTools
import os, sys, re

Ass = Serializer.LoadData( args.ass_FN )
Proj = Serializer.LoadFromHDF( args.proj_FN )


if (Ass.max()+1) > 100:
   print "You have %d states... This is going to be a large pdf file..." % (Ass.max()+1)

pp = PdfPages( args.out_FN )

Data = dataIO.readData( args.data_FN )

if len(Data.shape) == 1:
   Data = msmTools.reshapeRawData( Data, Proj )

x0 = 0
x1 = int(Data.max()+1)

if args.x_lbl != None:
   x_lbl = ' '.join( args.x_lbl )
else:
   x_lbl = 'Data in State'

for i in range( Ass.max() + 1 ):

   figure()
   state_dat = Data[ np.where( Ass == i ) ]
   
   hist( state_dat, bins=50, log=True,histtype='step',lw=2,range=(x0,x1) )
   #twinx()
   hist( Data.flatten(),bins=50, log=True,alpha=0.4,range=(x0,x1), histtype='stepfilled',color='gray' )
   yscale('symlog')
   title('State %d'%i)
   xlabel( x_lbl)
   ylabel('Frequency')
 
   pp.savefig()

   close()

pp.close()
