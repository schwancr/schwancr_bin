#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-s',dest='pdbFN',help='PDB Filename to use for the native state')
parser.add_argument('-d',dest='dataFN',help='Data with the state HB maps')
parser.add_argument('-o',dest='outFN',default='StateAvgHBs_vs_N.pdf')
args = parser.parse_args()
 
import numpy as np
from msmbuilder import Trajectory, Serializer
from schwancrtools import metric_HB
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
import os, sys, re

def DrawSquares( CM ):
   nonzero = np.array( np.where( CM ) ).T
   
   for (i,j) in nonzero:
      plot( [j,j,j+1,j+1,j], [i,i+1,i+1,i,i], color='red',lw=2 )
 
print "Loading data..."

Pdb = Trajectory.LoadTrajectoryFile( args.pdbFN )

CMs = np.load( args.dataFN )

pp = PdfPages( args.outFN )

hb = metric_HB.HydrogenBond()


pPdb = hb.prepare_trajectory( Pdb )

num_donors = len( hb.last_donor_ainds )
num_acceptors = len( hb.last_acceptor_ainds )

NatCM = pPdb[0].reshape( ( num_donors, num_acceptors ), order='C' )

for state_ind in xrange( CMs.shape[0] ):
 
   AvgCM = CMs[state_ind]

   figure()
      
   imshow( AvgCM, interpolation='nearest',vmin=0,vmax=1,cmap='gray_r',extent=(0,num_acceptors,num_donors,0) )
   xlabel('Acceptor ID')
   ylabel('Donor ID')
   colorbar().set_label('Average H-Bonding in State')
   gca().grid(True)
   title('State %d' % state_ind )

   DrawSquares( NatCM )

   xlim([ 0, num_acceptors ])
   ylim([ num_donors, 0 ])
   pp.savefig()
   close()

pp.close()


