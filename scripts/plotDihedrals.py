#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-s','--pdb',dest='pdbFN',help='PDB to use in finding the residue ID\'s for each of the angles')
parser.add_argument('-p','--pca',dest='pcaFN',help='PCA object. Serliazer object with vecs and vals')
parser.add_argument('-n','--angles',dest='angles',nargs='+',help='Angles used in calculating the PCs. One or more of [ phi, psi, chi, omega ]. Any order is fine, this script will sort them as msmbuilder.geometry.dihedral sorts them.')
parser.add_argument('-o','--out',dest='outFN',help='Output filename (should be PDF) [ DihedralPC1Weights.pdf ]', default='DihedralPC1Weights.pdf')
parser.add_argument('-N',dest='N',default=0,type=int,help='Which eigenvector to look at.')
parser.add_argument('--double',dest='double',default=False,action='store_true',help='Pass this flag if you used msmbuilder.metrics.Dihedrals, which means there is a sin and cosine entry for each angle')
options = parser.parse_args()
 
import numpy as np
from msmbuilder import Serializer, Trajectory
from msmbuilder import metrics
from msmbuilder.geometry import dihedral
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
import os, sys, re
 
pdb = Trajectory.LoadTrajectoryFile( options.pdbFN )

pca = Serializer.LoadFromHDF( options.pcaFN )

decInd = np.argsort( pca['vals'] )[::-1]


v0 = np.abs(pca['vecs'][:,decInd][:,options.N])
if options.double:
    if v0.shape[0] % 2:
        print "There are an odd number of entries, so --double should not be passed here, or something else has gone wrong."
        exit()
 
    n0 = v0.shape[0]
    v0 = v0[:n0/2] + v0[n0/2:]



angles = options.angles

print "Loaded Data."

figure()
count = 0

if 'chi' in angles:
    Rs = [ pdb['ResidueID'][row[0]] for row in dihedral._get_indices_chi( pdb ) ]
    N = len( Rs )
    plot( Rs, v0[ count : count + N ],marker='o',color='blue', label = 'Chi Angles' )
    count += N
    
if 'omega' in angles:
    Rs = [ pdb['ResidueID'][row[0]] for row in dihedral._get_indices_omega( pdb ) ]
    N = len( Rs )
    plot( Rs, v0[ count : count + N ],marker='s',color='cyan', label = 'Omega Angles' )
    count += N

   
if 'phi' in angles:
    Rs = [ pdb['ResidueID'][row[0]] for row in dihedral._get_indices_phi( pdb ) ]
    N = len( Rs )
    plot( Rs, v0[ count : count + N ],marker='^',color='purple', label = 'Phi Angles' )
    count += N

if 'psi' in angles:
    Rs = [ pdb['ResidueID'][row[0]] for row in dihedral._get_indices_psi( pdb ) ]
    N = len( Rs )
    plot( Rs, v0[ count : count + N ],marker='*',color='red', label = 'Psi Angles' )
    count += N

if count != len( v0 ):
    print "You did not input the right angles... Number of residues we looped through here is %d, but the size of the first eigenvector is %d" % ( count, len(v0) )

legend()
xlabel('Residue ID')
ylabel('Norm of Weight in PC1')
ylim([0,1.2*ylim()[1]])
print "Saved Plot to %s" % options.outFN 
savefig( options.outFN )

