#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_argument('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_argument('--ot',dest='out_cm',default='./StateAvgCMs.npy',help='Average CM\'s for each state [ ./StateAvgCMs.npy ]')
parser.add_argument('--op',dest='out_plot',default='./StateAvgCMs.pdf',help='Output for plotting each state\'s CM. Will only plot if there are fewer than 20 states [ ./StateAvgCMs.pdf ]')
parser.add_argument('--co',dest='cutoffs',default=None,help='Cutoffs file for residue-residue cutoffs')
parser.add_argument('-s',dest='pdb_FN',help='PDB file for the native state.')
args = parser.parse_args()
 
import numpy as np
import matplotlib
matplotlib.use('pdf')
from msmbuilder import io, Project, metrics, Trajectory
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
import os, sys, re
 
nat_pdb = Trajectory.load_from_pdb( args.pdb_FN )
Proj = Project.load_from( args.proj_FN )
try: Ass = io.loadh(args.ass_FN)['Data']
except: Ass = io.loadh(args.ass_FN)['arr_0']

if args.cutoffs != None:
   cuts = np.loadtxt( args.cutoffs )
   if len(cuts.shape)!=1:
      if cuts.shape[1] == 3:
         contacts=(cuts[:,:2]-1).astype(int)
         cuts = cuts[:,2]
      else:
         print "Cutoff file needs to have residue ID's in the first two columns and the corresponding cutoff in the third."
         exit()
      BoolCont = metrics.BooleanContact( contacts=contacts, cutoff=cuts )
else:
   BoolCont = metrics.BooleanContact( scheme='closest-heavy' )



if ( Ass.max() < 100 ): 
   pp = PdfPages( args.out_plot )
else:
   pp = None

StateAssigns = np.array([ np.where( Ass == i )[0].shape[0] for i in np.unique( Ass[ np.where( Ass >= 0 ) ] )] )
StateAssigns = StateAssigns.reshape( (len(StateAssigns),1) )

CMs_1d = None
loaded_CMs = None
if os.path.exists( args.out_cm ):
   CMs = np.load( args.out_cm )
   n_res = CMs.shape[1]
   loaded_CMs = CMs
   AvgCMs_1d = np.array([ r.flatten(order='C') for r in CMs ]) # THIS DOES NOT WORK!!!
else:
   for traj_ind in xrange( Ass.shape[0] ):
      print "Working on %s" % Proj.traj_filename(traj_ind) 
      Traj = Proj.load_traj( traj_ind )
      pTraj = BoolCont.prepare_trajectory( Traj ).astype(float)
      n_res = np.unique( Traj['ResidueID'] ).shape[0]
      del Traj
      if CMs_1d == None:
         CMs_1d = np.zeros( ( Ass.max()+1, pTraj.shape[1] ) )
 
      if ( (n_res-2) * (n_res-3) / 2 != pTraj.shape[1] ) and ( not args.cutoffs ):
         print "Number of residues is %d, but the shape of pTraj (%d) doesn't fit!" % ( n_res, pTraj.shape[1] )

      for i,state in enumerate( np.unique( Ass[ np.where( Ass >= 0 ) ] ) ):
         CMs_1d[i] += pTraj[ np.where( Ass[traj_ind] == state ) ].sum(axis=0) # remember this '== state' used to be '== i'
   #print "Remember that you changed Ass[traj_ind] == i to Ass[traj_ind] == state without testing. So test it now."

      del pTraj
   AvgCMs_1d = CMs_1d / StateAssigns
   CMs = []

CMs = []
nat_CM1d = BoolCont.prepare_trajectory( nat_pdb )[0]
sum=0
nat_CM = np.zeros( (n_res,n_res) )
for i in xrange( n_res - 3 ):
   nat_CM[i][ (3+i): ] = nat_CM1d[ sum : sum + n_res - 3 -i ]
   sum+= n_res-3-i
nat_CM = (nat_CM+nat_CM.T)/2.
nat_CM += np.eye( n_res, k= 0)
nat_CM += np.eye( n_res, k= 1)
nat_CM += np.eye( n_res, k=-1)
nat_CM += np.eye( n_res, k=-2)
nat_CM += np.eye( n_res, k= 2)

nat_locs = np.array( np.where( nat_CM ) ).T

for state_ind in xrange( AvgCMs_1d.shape[0] ):
   print state_ind
   if loaded_CMs != None:
      AvgCM = loaded_CMs[state_ind]
   else:
      AvgCM = np.zeros( (n_res,n_res) )
      sum=0

      if args.cutoffs == None:
         for i in xrange( n_res - 3):
          #  print i, sum + n_res - 3 - i
            AvgCM[i][ (3+i): ] = AvgCMs_1d[state_ind][ sum : sum + n_res - 3 - i ]
            sum += n_res-3-i
         AvgCM += np.eye( n_res, k= 0  )
         AvgCM += np.eye( n_res, k= 1  )
         AvgCM += np.eye( n_res, k=-1 )
         AvgCM += np.eye( n_res, k= 2 )
         AvgCM += np.eye( n_res, k=-2 )
      else:
         for i in range( len( contacts ) ):
            AvgCM[ contacts[i][0], contacts[i][1] ] = AvgCMs_1d[state_ind][i]

      AvgCM += AvgCM.T

   CMs.append( AvgCM )
   if pp:
      figure()
      imshow( AvgCM, interpolation='nearest',vmin=0,vmax=1,cmap='gray_r',extent=(0,n_res,n_res,0) )
      for (x,y) in nat_locs:
          plot( [y,y,y+1,y+1,y], [x,x+1,x+1,x,x], color='red')
      xlabel('Residue ID')
      gca().grid(True)
      ylabel('Residue ID')
      colorbar().set_label('Average Contacting in State')
      title('State %d (N=%d)' % (np.unique( Ass[np.where( Ass >= 0 )] )[state_ind], StateAssigns[state_ind])  )
      #text(1,n_res-1,'N=%d'%StateAssigns[state_ind])
      pp.savefig()
      close()

if loaded_CMs == None:
   np.save( args.out_cm, np.array(CMs) )


if pp:
   pp.close()
