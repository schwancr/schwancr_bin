#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_argument('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_argument('--ot',dest='out_cm',default='./StateAvgHBs.h5',help='Average HB\'s for each state, along with the corresponding atom indices [ ./StateAvgHBs.h5 ]')
parser.add_argument('--op',dest='out_plot',default='./StateAvgHBs.pdf',help='Output for plotting each state\'s HB. Will only plot if there are fewer than 20 states [ ./StateAvgHBs.pdf ]')
parser.add_argument('-s',dest='nat_FN',help='Native state to use to compare the contact maps')
args = parser.parse_args()
 
import numpy as np
import matplotlib
matplotlib.use('pdf')
from msmbuilder import io, Project, metrics, Conformation, Trajectory
from schwancrtools import metric_HB
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
import os, sys, re
 
Proj = Project.load_from( args.proj_FN )
try: Ass = io.loadh( args.ass_FN )['arr_0']
except: Ass = io.loadh( args.ass_FN )['Data']

HB = metric_HB.HydrogenBond()

pdb = Trajectory.load_from_pdb( Proj.conf_filename )

donorH_ainds = np.where( ( (pdb['AtomNames'] == 'H')|(pdb['AtomNames'] =='HN' ) )&(pdb['ResidueNames']!='PRO') )[0] # _ainds correspond to atom indices in the trajectory
donor_ainds = np.where( (pdb['AtomNames'] == 'N')&(pdb['ResidueNames']!='PRO') )[0]
acceptor_ainds = np.where( pdb['AtomNames'] == 'O' )[0] # This excludes the C-Terminus which has naming issues...


n_res = np.unique( pdb['ResidueID'] ).shape[0]
ppdb = HB.prepare_trajectory(pdb) # Do this to set the donor indices and stuff

atom_indices = np.sort( np.concatenate( ( donor_ainds, donorH_ainds, acceptor_ainds ) ) )

CMs_1d = np.zeros( ( Ass.max()+1, ppdb.shape[1] ) )

if os.path.exists( args.out_cm ):
   print "Data file exists, will use it and just re-plot the data"

   output_dict = io.loadh( args.out_cm )
   triples = output_dict['donor_h_acceptor_ainds']
   CMs = output_dict['HB_maps']
else:
   CMs = None

if ( Ass.max() < 250 ): 
   pp = PdfPages( args.out_plot )
else:
   pp = None

chunk_size=10000

if CMs == None:

   for traj_ind in xrange( Ass.shape[0] ):
      print "Working on %s" % Proj.traj_filename(traj_ind) 
      for chunk_ind, trj_chunk in enumerate(Trajectory.enum_chunks_from_lhdf( Proj.traj_filename(traj_ind), 
                                            ChunkSize=chunk_size, AtomIndices=atom_indices ) ):
         print "chunked"
         ptrj_chunk = HB.prepare_trajectory( trj_chunk ).astype(float)
         ass_chunk = Ass[traj_ind][ chunk_ind*chunk_size : (chunk_ind+1)*chunk_size ] # this behaves as you want at the end of the array
         for i,state in enumerate( np.unique( Ass[ np.where( Ass >= 0 ) ] ) ):
            CMs_1d[i] += ptrj_chunk[ np.where( ass_chunk == state ) ].sum(axis=0)

   StateAssigns = np.array([ np.where( Ass == i )[0].shape[0] for i in np.unique( Ass[ np.where( Ass >= 0 ) ] )] )
   StateAssigns = StateAssigns.reshape( (len(StateAssigns),1) )
   AvgCMs_1d = CMs_1d / StateAssigns

   num_donors = len( HB.last_donor_ainds )
   num_acceptors = len( HB.last_acceptor_ainds )

   CMs = AvgCMs_1d.reshape( (-1,num_donors,num_acceptors),order='C')

   triples = HB.get_angle_list()
   io.saveh(args.out_cm, donor_h_acceptor_ainds=triples, HB_maps=np.array(CMs))
#CMs = [ avg_cm.reshape( (num_donors, num_acceptors ), order='C') for avg_cm in AvgCMs_1d ]

num_acceptors = len( np.unique( triples[:,2] ) )
num_donors = len( np.unique( triples[:,0] ) )

donor_ainds = np.unique( triples[:,0] )
donorH_ainds = np.unique( triples[:,1] )
acceptor_ainds = np.unique( triples[:,2] )


donor_atomnames = pdb['AtomNames'][ donor_ainds ]
donorH_atomnames = pdb['AtomNames'][ donorH_ainds ]
acceptor_atomnames = pdb['AtomNames'][ acceptor_ainds ]

donor_resnames = pdb['ResidueNames'][ donor_ainds ]
donorH_resnames = pdb['ResidueNames'][ donorH_ainds ]
acceptor_resnames = pdb['ResidueNames'][ acceptor_ainds ]

donor_resIDs = pdb['ResidueID'][ donor_ainds ]
donorH_resIDs = pdb['ResidueID'][ donorH_ainds ]
acceptor_resIDs = pdb['ResidueID'][ acceptor_ainds ]

nat_pdb = Trajectory.load_from_pdb( args.nat_FN )

CM_pdb = HB.prepare_trajectory( nat_pdb )

CM_pdb = CM_pdb.reshape( (1,num_donors, num_acceptors), order='C')[0]

native_locs = np.array(np.where( CM_pdb )).T


if pp:

   for state_ind in xrange( CMs.shape[0] ):
      AvgCM = CMs[state_ind]

      DiffHBs = np.array(np.where( (AvgCM - CM_pdb)>0.7 )).T#CMs[9]) > 0.7 )).T
      
      Labels = [ '%d%s %s - %d%s %s' % (donorH_resIDs[i], donorH_resnames[i], donorH_atomnames[i], acceptor_resIDs[j], acceptor_resnames[j], acceptor_atomnames[j] ) for (i,j) in DiffHBs ]
      print state_ind
      print Labels
      figure()
      
      imshow( AvgCM, interpolation='nearest',vmin=0,vmax=1,cmap='gray_r',extent=(0,num_acceptors,num_donors,0) )
      n=10
      for (x,y) in native_locs:
          plot( [y,y,y+1,y+1,y], [x,x+1,x+1,x,x], color='red')
      #xlim([ 0-n, num_acceptors+n ])
      #ylim([ num_donors+n, 0-n ])
      xlabel('Acceptor ID')
      ylabel('Donor ID')
      colorbar().set_label('Average H-Bonding in State')
      gca().grid(True)
      state_ID = np.unique( Ass[np.where( Ass >= 0 )] )[state_ind] 
      title('State %d (N=%d)' % (state_ID, np.where(Ass == state_ID)[0].shape[0]) )
      #imshow( CM_pdb, alpha=0, vmin=0,vmax=1, cmap='hot_r',extent=(0,num_acceptors,num_donors,0) )

    #  for i in range( len( Labels ) ):
          #text( DiffHBs[i][1]+1, DiffHBs[i][0], Labels[i] )
    #      gca().add_patch( Circle( (DiffHBs[i][1]+0.65, DiffHBs[i][0]+0.65 ), 1.5, fill=False ) )
      pp.savefig()
      close()

   pp.close()
