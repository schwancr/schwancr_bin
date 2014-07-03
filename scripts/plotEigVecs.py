#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-t',dest='T_fn',default='./tProb.mtx',help='Transition probability matrix from MSMBuilder.')
parser.add_argument('-k',dest='num_vecs',default=1,type=int,help='Number of dynamic eigenvectors to plot')
parser.add_argument('-d',dest='ord_param',help='Order parameter to use in sorting the eigenvector entries')
parser.add_argument('-o',dest='out_fn',default='Eigenvecs.pdf',help='Output filename (SHOULD BE PDF)')
parser.add_argument('-O',dest='pop_fn',default='./Populations.dat',help='Equilibrium populations from MSMBuilder')

options = parser.parse_args()
 
import numpy as np
from scipy.io import mmread
#from scipy.sparse.linalg import eigs
from msmbuilder import msm_analysis
import matplotlib
matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter
import os, sys, re
 
pdf = PdfPages( options.out_fn )

num_vecs = options.num_vecs
ord_param = np.loadtxt( options.ord_param )

if options.ord_param == 'MacroMapping.dat':
   macro = True
else:
   macro = False
if len( ord_param.shape ) == 2:
   print "2D array inputted. I will use the second column (i.e. data[:,1])"
   ord_param = ord_param[:,1]

T = mmread( options.T_fn )
vals,vecs = msm_analysis.get_eigenvectors( T, num_vecs+1 )
#vals, vecs = eigs( T, k = num_vecs + 1 )

vecs=vecs.real

pi = vecs[:,0] / vecs[:,0].sum()

vecs /= np.sum( np.square( vecs ) / pi.reshape( (-1,1) ), axis=0 )

ord_param_ind = ord_param.argsort()
state_lines = np.array([ np.where( ord_param[ ord_param_ind ] == i )[0][0] for i in np.unique( ord_param ) ])[1:]


for i in range( num_vecs ):
   figure()
   eval = vals[i+1]
   evec = vecs[:,i+1]
   flux = evec #* pops
   flux = flux[ord_param_ind]
   if len(evec)<=10:
      lw=5
   elif 10<len(evec)<=100:
      lw=2
   else:
      lw=1

   pos_ind = np.where( flux >= 0 )[0]
   neg_ind = np.where( flux < 0 )[0]
   # Positive values are black
   vlines( pos_ind, np.zeros( len(pos_ind) ), flux[pos_ind], color='black', lw=lw)
   # Negative values are red
   vlines( neg_ind, np.zeros( len(neg_ind) ), flux[neg_ind], color='red', lw=lw)
   ax=gca()
   ax.yaxis.set_major_formatter(FormatStrFormatter('%.1E'))
   title('Dynamic Eigenvector Number %d' % (i+1))
   xlabel('State Index')
   xlim([ -xlim()[1]*0.05, xlim()[1]*1.05] )
   ylabel('Flux')

   if macro:
      ymin,ymax = ylim()
      vlines( state_lines, ymin,ymax,color='blue' )
      ylim([ ymin,ymax ])
      xlim([-xlim()[1]*0.05, xlim()[1]*1.05])
   pdf.savefig()
 
np.save('eigenvalues.npy',vals) 
   
pdf.close()
