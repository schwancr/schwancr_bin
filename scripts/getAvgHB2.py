#!/usr/bin/env python
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_argument('-a',dest='ass_FN',default='./Assignments.Fixed.h5',help='Assignments from msmbuilder [ ./Assignments.Fixed.h5 ]')
parser.add_argument('--ot',dest='out_cm',default='./StateAvgHBs.h5',help='Average HB\'s for each state, along with the corresponding atom indices [ ./StateAvgHBs.h5 ]')
parser.add_argument('--op',dest='out_plot',default='./StateAvgHBs.pdf',help='Output for plotting each state\'s HB. Will only plot if there are fewer than 20 states [ ./StateAvgHBs.pdf ]')
parser.add_argument('-s',dest='nat_FN',help='Native state to use to compare the contact maps')
parser.add_argument('-w', dest='which', help='which contacts to look for. Should have three columns corresponding to (acceptor, acceptor-hydrogen, donor) atom indices')
args = parser.parse_args()
 
import numpy as np
import matplotlib
matplotlib.use('pdf')
from msmbuilder import io, Project, metrics, Conformation, Trajectory
from schwancrtools import angle
from msmbuilder.geometry import contact
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
import os, sys, re
from time import time
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
sh = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(fmt='%(asctime)s - %(message)s', datefmt="%H:%M:%S")
sh.setFormatter(formatter)
logger.addHandler(sh)
logger.propagate = False

logger.info("start")
Proj = Project.load_from( args.proj_FN )
logger.info("loaded project info")
try: Ass = io.loadh( args.ass_FN )['arr_0']
except: Ass = io.loadh( args.ass_FN )['Data']

pdb = Trajectory.load_from_pdb( Proj.conf_filename )

which = np.loadtxt(args.which).astype(int)

distance_cutoff = 0.32
angle_cutoff = 120

def get_hb(traj):

    # get accH - donor distance:
    dists = contact.atom_distances(traj['XYZList'], atom_contacts=which[:, 1:])

    # get angles
    angles = angle.compute_angles(traj, which)

    hb = ((dists < distance_cutoff) & (angles > angle_cutoff)).astype(int)

    return hb
    
n_res = np.unique( pdb['ResidueID'] ).shape[0]
ppdb = get_hb(pdb) # Do this to set the donor indices and stuff

CMs_1d = np.zeros( ( Ass.max()+1, ppdb.shape[1] ) )

if os.path.exists( args.out_cm ):
   logger.warning("Data file exists, will use it and just re-plot the data")

   output_dict = io.loadh( args.out_cm )
   which = output_dict['which']
   AvgCMs_1d = output_dict['HB_maps']
else:
   AvgCMs_1d = None

if ( Ass.max() < 250 ): 
   pp = PdfPages( args.out_plot )
else:
   pp = None

chunk_size=10000

if AvgCMs_1d == None:

   for traj_ind in xrange( Ass.shape[0] ):
      logger.info("Working on %s" % Proj.traj_filename(traj_ind))
      for chunk_ind, trj_chunk in enumerate(Trajectory.enum_chunks_from_lhdf( Proj.traj_filename(traj_ind), 
                                            ChunkSize=chunk_size) ):
         logger.debug("chunked")
         ptrj_chunk = get_hb( trj_chunk ).astype(float)
         ass_chunk = Ass[traj_ind][ chunk_ind*chunk_size : (chunk_ind+1)*chunk_size ] # this behaves as you want at the end of the array

         for i, ass in enumerate(ass_chunk):
            if ass == -1:
                continue
            CMs_1d[ass] += ptrj_chunk[i]
        
   #StateAssigns = np.array([ np.where( Ass == i )[0].shape[0] for i in np.unique( Ass[ np.where( Ass >= 0 ) ] )] )
   StateAssigns = np.bincount(Ass[np.where(Ass != -1)], minlength=Ass.max() + 1)
   StateAssigns = StateAssigns.reshape( (len(StateAssigns),1) )
   AvgCMs_1d = CMs_1d / StateAssigns

   io.saveh(args.out_cm, which=which, HB_maps=AvgCMs_1d)


uniq_res = np.unique(pdb['ResidueID'])
n_res = uniq_res.shape[0]

acc_res_ids = pdb['ResidueID'][which[:, 0]]
donor_res_ids = pdb['ResidueID'][which[:, 2]]

CMs = np.zeros((len(AvgCMs_1d), n_res, n_res))
CM_pdb = np.zeros((n_res, n_res))
for i in xrange(n_res):
    for j in xrange(n_res):
        if i == j:
            continue
        inds = np.where((acc_res_ids == uniq_res[i]) & (donor_res_ids == uniq_res[j]))[0]
        if len(inds) == 0:
            CMs[:, i, j] = -1
            CM_pdb[i, j] = -1

        else:
            CMs[:, i, j] = AvgCMs_1d[:, inds].sum(axis=1)
            CM_pdb[i, j] = ppdb[:, inds].sum()
            # change this line to get different behavior

native_locs = np.array(np.where(CM_pdb > 0.5)).T
if pp:

   for state_ind in xrange( CMs.shape[0] ):
      AvgCM = CMs[state_ind]

      print AvgCM.shape

      DiffHBs = np.array(np.where( (AvgCM - CM_pdb)>0.7 )).T#CMs[9]) > 0.7 )).T
      
#      Labels = [ '%d%s %s - %d%s %s' % (donorH_resIDs[i], donorH_resnames[i], donorH_atomnames[i], acceptor_resIDs[j], acceptor_resnames[j], acceptor_atomnames[j] ) for (i,j) in DiffHBs ]

      figure()
      
      imshow( AvgCM, interpolation='nearest', vmin=-1, vmax=1, cmap='RdBu',extent=(0,n_res,n_res,0) )
      n=10
      for (x,y) in native_locs:
          plot( [y,y,y+1,y+1,y], [x,x+1,x+1,x,x], color='red')
      xlim(0, n_res)
      ylim(n_res, 0)
      
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
