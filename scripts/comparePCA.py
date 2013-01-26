#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-d',dest='dataFNs',action='append',help='Data files to compare.')
parser.add_argument('-l',dest='labels',action='append',nargs='+',help='Label for the corresponding data file')
parser.add_argument('-N',dest='num_vecs',type=int,default=10,help='Number of eigenvectors to compare [ 10 ]')
parser.add_argument('-o',dest='outFN',default='comparePCA.pdf',help='Output filename [ comparePCA.pdf ]')
parser.add_argument('--major-label',nargs='+',dest='major_label',default=None,help='Major label for each axis.')
options = parser.parse_args()
 
import numpy as np
from Emsmbuilder import Serializer
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from pyschwancr import dataIO, msmTools
import os, sys, re

outPtr = PdfPages( options.outFN )

PCAs = [ Serializer.LoadFromHDF( fn ) for fn in options.dataFNs ]
Labels = [ ' '.join( lbl ) for lbl in options.labels ]

N = int( options.num_vecs )

HighestEigInd = [ np.argsort( p['vals'].real.astype(float) )[-N:] for p in PCAs ] # use .real.astype... to avoid complexwarning (these are real anyway)

TotVar = [ np.sum( PCAs[i]['vals'].real.astype(float) ) for i in range( len( PCAs ) ) ]
ExplVar = [ np.sum( PCAs[i]['vals'][HighestEigInd[i]].real.astype(float) ) for i in range( len( PCAs ) ) ]
ExplVarProp = [ e / t for (t,e) in zip( TotVar, ExplVar ) ]

SimMat = np.ones( ( len( PCAs ), len( PCAs ) ) ).astype(complex)

for i in range( len( PCAs ) ):
    for j in range( i+1, len( PCAs ) ):
        #SimMat[i,j] = np.linalg.norm( np.sum( PCAs[i]['vecs'][:, HighestEigInd[i] ] * np.conj( PCAs[j]['vecs'][:, HighestEigInd[j]] ) ) / float( N )  )
        SimMat[i,j] = ( np.sum( np.abs( PCAs[i]['vecs'][:, HighestEigInd[i] ] * np.conj( PCAs[j]['vecs'][:, HighestEigInd[j]] ) ) ) / float( N )  ).real
        SimMat[j,i] = SimMat[i,j]

SimMat = SimMat.real.astype(float) # use .real.astype(... to avoid comlpexwarning (these are real anyway)

figure( figsize=(8,8) )

matshow( SimMat, vmin=0,vmax=1,cmap='copper',fignum=False )
xticks( np.arange( len(PCAs) ), Labels, rotation='vertical' )
yticks( np.arange( len(PCAs) ), Labels )

ax=gca()
ax.xaxis.set_label_position('top')
if options.major_label:
    xlabel( ' '.join( options.major_label ) )
    ylabel( ' '.join( options.major_label ) )

colorbar().set_label('Average Similarity of the First %d PCs' % N )

outPtr.savefig()

figure(figsize=(8,8))
width=0.8
bar( np.arange( len( PCAs ) ), ExplVarProp, width=width)

xticks( np.arange( len( PCAs ) ) + width*0.5 , Labels, rotation='vertical' )
ylim([0,1.05])
ylabel('Proportion of Total Variance Explained')

if options.major_label:
    xlabel( ' '.join( options.major_label ) )

hlines( 1,0,len( PCAs ), color='red')
title('Total Variance in top %d PCs' % N )
outPtr.savefig()



outPtr.close()


