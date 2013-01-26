#!/usr/bin/env python

"""
This library file is used for defining dimensionality reduction techniques for converting an MSMbuilder project to a lower dimension coordinate
"""

try: import mdp
except: print "Need to install mdp-toolkit in order to use this library"
import numpy as np
from Emsmbuilder import Trajectory
from schwancrtools import metrics_PCA
from Emsmbuilder import metric_LPRMSD as lprmsd
from pyschwancr.msmTools import calc_time_avg
import re, os, sys, copy

def GetTrajList( dir, BeginInd=3, EndInd=-4, RegEx=r'^trj\d+\.lh5' ):
   """
   This function creates a list of data files indexed by integers, and named in some regular way. T
   The files should be named in such a way that contains an integer. 
   
   Inputs:
   1) dir - Directory containing the files
   2) BeginInd [ 3 ] - Index which is the first integer digit in the filename
   3) EndInd [ -4 ] - Index which is one past the last integer digit in the filename
   4) RegEx [ r\'^trj\\d+\\.lh5\' ] - Regular expression to use to filter things to place in the list

   Outputs:
   1) trajList - List of filenames INCLUDING THE DIRECTORY (i.e. os.path.join( dir, fn ) )

   EG)
   
   Trajectories:  
   trj0.lh5  trj1.lh5  trj2.lh5  trj3.lh5 ...

   getTrajList( \'Trajectories\', BeginInd=3, EndInd=-4, RegEx=r\'trj\\d+\\.lh5\' )

   [ Trajectories/trj0.lh5, Trajectories/trj1.lh5, ... ]

   Since fn[3] = 0, fn[-4] = \'.\'

   """
   trajList = [ ( int( fn[ BeginInd : EndInd ] ), fn ) for fn in os.listdir(dir) if re.search( RegEx, fn ) ]
   trajList.sort()
   trajList = [ os.path.join( dir, fn ) for ( i, fn ) in trajList ]

   return trajList


class PCA:
    """
    This class uses PCA to convert trajectories to some number of coordinates.
    """

    def __init__(self, pdbFN=None, input_dir='./Trajectories', prep_with=None,atomindices=None,expl_var=0.95, numComp=None, stride=1, time_avg=1):
        self.input_dir = input_dir
        
        if not numComp:
            self.output_dim = expl_var
        else:
            self.output_dim = numComp

        self.pca = mdp.nodes.PCANode( output_dim = self.output_dim )
        self.stride = stride

        if not atomindices:
            self.atomindices = None
        elif isinstance(atomindices,np.ndarray):
            self.atomindices = atomindices
        else:
            self.atomindices = np.loadtxt( atomindices ).astype(int)

        self.rmsd = lprmsd.LPRMSD( atomindices = self.atomindices )

        if pdbFN:
            self.align2pdb = self.rmsd.prepare_trajectory( lprmsd.LPTraj.LoadFromPDB( pdbFN ) )
            self.use_positions = True
        else:
            self.use_positions = False
            if prep_with:
                self.prep_with = prep_with
            else:
                raise Exception('Must provide one of prep_with or pdbFN')

        self.time_avg = int(time_avg / self.stride)
        if self.time_avg > 1:
            self.do_time_avg = True
        else:
            self.do_time_avg = False

    def train(self):
        
        trajList = GetTrajList( self.input_dir )

        for fn in trajList:
            traj = lprmsd.LPTraj.LoadTrajectoryFile( fn )[::self.stride]
            print "Working on %s" % fn
            #if self.atomindices:
            #    traj.RestrictAtomIndices( self.atomindices )

            if self.use_positions:
                traj = self.rmsd.prepare_trajectory( traj )
                rmsds,traj = self.rmsd.one_to_all_aligned( self.align2pdb, traj, 0 )
                del rmsds
                traj = traj[:,self.atomindices]
            else:
                traj = self.prep_with.prepare_trajectory( traj )

            if self.do_time_avg:
                traj = calc_time_avg( traj, self.time_avg )

            if len(traj.shape) == 3:
                n0,n1,n2 = traj.shape
                traj = traj.reshape( n0, n1*n2 )
            elif len(traj.shape) !=2:
                traj = np.array([ frame.flatten() for frame in traj[1] ]).astype('float64')
            self.pca.train( traj )
            del traj

    def savePCA(self,outFN='pca.npy'):
        self.pca.save(outFN)

