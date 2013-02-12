#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--hb',dest='in_cm',default='./StateAvgHBs.h5',help='Average HB\'s for each state, along with the corresponding atom indices [ ./StateAvgHBs.h5 ]')
parser.add_argument('--op',dest='out_plot',default='./VecHBs.pdf',help='Output for plotting the HB pattern for each dynamic eigenvector')
parser.add_argument('-s',dest='nat_FN',help='Native state to use to compare the contact maps')
parser.add_argument('-t',dest='tProb',default='./tProb.mtx',help='tProb.mtx for your MSM.')
parser.add_argument('-n',dest='num_vecs',type=int,default=10,help='Number of DYNAMIC eigenvectors to plot')
parser.add_argument('-l',dest='lagtime',type=float,default=None,help='Lagtime to use in calculating the timescale for each eigenvector. These should be in units of whatever you pass to --units')
parser.add_argument('--units',dest='units',default=None,help='Units that your lagtime is in')
args = parser.parse_args()
 
import numpy as np
from scipy.io import mmread
import matplotlib
matplotlib.use('pdf')
from msmbuilder import io, metrics, Trajectory, MSMLib
from schwancrtools import metric_HB
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from pyschwancr import dataIO
import os, sys, re
 
matplotlib.rcParams['font.size']=22
pp = PdfPages( args.out_plot )

HB = metric_HB.HydrogenBond()

input_dict = io.loadh( args.in_cm )
triples = input_dict['donor_h_acceptor_ainds']
CMs = input_dict['HB_maps']

num_acceptors = len( np.unique( triples[:,2] ) )
num_donors = len( np.unique( triples[:,0] ) )

nat_pdb = Trajectory.load_from_pdb( args.nat_FN )
T = mmread( args.tProb )
n_res = np.unique( nat_pdb['ResidueID'] ).shape[0]
CM_pdb = HB.prepare_trajectory( nat_pdb )

CM_pdb = CM_pdb.reshape( (1,num_donors, num_acceptors), order='C')[0]

native_locs = np.array(np.where( CM_pdb )).T

vals,vecs = MSMLib.GetEigenvectors( T, args.num_vecs+1 )
vecs = vecs.real
    
if args.lagtime != None:
   timescales = - args.lagtime / np.log( vals[1:].real ) # NOTE: This means the 1st dyn e-vec is indexed at 0
  
for i in xrange(1,args.num_vecs+1):
   figure(figsize=(12,6))
   #AvgCM = np.array([ CMs[j] * vecs[:,i][j] for j in xrange( len( CMs ) ) ]).sum(axis=0)
   posCM = np.array([ CMs[j] * vecs[:,i][j] for j in np.where( vecs[:,i] < 0 )[0] ]).sum(axis=0)
   negCM = np.array([ CMs[j] * vecs[:,i][j] for j in np.where( vecs[:,i] > 0 )[0] ]).sum(axis=0)
   # Woops... I named those wrong..
   max_range = np.max([ posCM.max(), - posCM.min() ]) #np.abs(AvgCM).max()
   max_range = int( max_range * 2 ) / 2.
   
   #imshow( AvgCM, interpolation='nearest',vmin=-max_range,vmax=max_range,cmap='RdYlBu',extent=(0,num_acceptors,num_donors,0) )
   h=0.1
   axes( (h,2*h,0.5-2*h,1-4*h) )
   imshow( posCM, interpolation='nearest',vmin=-max_range,vmax=max_range,cmap='RdYlBu',extent=(0,num_acceptors,num_donors,0) )
   for (x,y) in native_locs:
      plot( [y,y,y+1,y+1,y], [x,x+1,x+1,x,x], color='black')
   n = 0
   xlim([ 0-n, num_acceptors+n ])
   ylim([ num_donors+n, 0-n ])
   xlabel('Acceptor ID')
   ylabel('Donor ID')
   #colorbar().set_label('Average HBs $\cdot$ Eigenvector')
   gca().grid(True)
   #title('Dynamic Eigenvector %d' % i)
   
   axes( (h+0.4,2*h,(0.5-2*h),1-4*h) )
   im=imshow( negCM, interpolation='nearest',vmin=-max_range,vmax=max_range,cmap='RdYlBu',extent=(0,num_acceptors,num_donors,0) )
   for (x,y) in native_locs:
      plot( [y,y,y+1,y+1,y], [x,x+1,x+1,x,x], color='black')
   n = 0
   xlim([ 0-n, num_acceptors+n ])
   ylim([ num_donors+n, 0-n ])
   xlabel('Acceptor ID')
   ylabel('Donor ID')
   gca().grid(True)
   cbar = colorbar( im, cax=axes( (h+0.75,2*h, 0.02, 1-4*h ) ) ) #.set_label('Average HBs $\cdot$ Eigenvector') 
#   suptitle('Dynamic Eigenvector %d' % i)
   if (args.lagtime != None) and (args.units != None):
      suptitle( 'Eigenvector %d ( %.2f %s )' %( i, timescales[i-1], args.units ), y=0.91 )
   else: 
      suptitle( 'Dynamic Eigenvector %d' % i )
   tx=[-max_range,-max_range/2.,0.0,max_range/2., max_range]
   lbls=[ '%.2f' % t for t in tx ]
   cbar.set_ticks(tx)
   cbar.set_ticklabels(lbls)
   #cbar.ax.yaxis.set_major_formatter('%.2f')
      #imshow( CM_pdb, alpha=0, vmin=0,vmax=1, cmap='hot_r',extent=(0,num_acceptors,num_donors,0) )

    #  for i in range( len( Labels ) ):
          #text( DiffHBs[i][1]+1, DiffHBs[i][0], Labels[i] )
    #      gca().add_patch( Circle( (DiffHBs[i][1]+0.65, DiffHBs[i][0]+0.65 ), 1.5, fill=False ) )
   pp.savefig()

pp.close()
