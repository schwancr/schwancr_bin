#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-M',dest='mfpt_list',action='append',help='Matrix of MFPTs output from GetMFPTs.py. There can be many of these entries')
parser.add_option('-S',dest='states_list',action='append',help='Filename corresponding to the list of states in the corresponding MFPT matrix')
parser.add_option('-d',dest='data_FN',help='Average value of some order parameter to determine if states are folded or unfolded')
parser.add_option('--cf',dest='cutF',type=float,default=0.4,help='Cutoff to discriminate folded and unfolded. This is the folded cutoff [ 0.4 ]' )
parser.add_option('--cu',dest='cutU',type=float,default=0.6,help='Cutoff to discriminate folded and unfolded. This is the unfolded cutoff [ 0.6 ]' )
parser.add_option('--low-is-folded',dest='low_is_folded',default=False,action='store_true',help='This script assumes that a high value in the order parameter means the state is a folded state. If the opposite is true (like when using RMSD), then pass this flag')
parser.add_option('-o',dest='out_FN',default='MFPT_hists.pdf',help='Output filename [ MFPT_hists.pdf ]')
parser.add_option('-t',dest='titles',default='ABC',help='Labels for the plots, should be one character each [ ABC ]')
options, args = parser.parse_args()

import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from pyschwancr import dataIO
import os, sys, re
 
print "Loading Data"

Ms = []
Sts = []

for mat_fn, state_fn in zip( options.mfpt_list, options.states_list ):
   Ms.append( dataIO.readData( mat_fn ) )
   Sts.append( dataIO.readData( state_fn ) )

Ms = np.hstack( Ms )
Sts = np.concatenate( Sts ).astype(int)

Data = dataIO.readData( options.data_FN )[:,1]
cutF = options.cutF
cutU = options.cutU

if options.low_is_folded:
   isFary = Data <= cutF
   isUary = Data >= cutU
else:
   isFary = Data >= cutF
   isUary = Data <= cutU

rowIsF = np.vstack( [ isFary ] * Ms.shape[1] ).T
colIsF = np.vstack( [ isFary[ Sts ] ] * Ms.shape[0] )
rowIsU = np.vstack( [ isUary ] * Ms.shape[1] ).T
colIsU = np.vstack( [ isUary[ Sts ] ] * Ms.shape[0] )

FF = Ms[ np.where( rowIsF * colIsF ) ]
UF = Ms[ np.where( rowIsU * colIsF ) ]
FU = Ms[ np.where( rowIsF * colIsU ) ]
UU = Ms[ np.where( rowIsU * colIsU ) ]

print Ms[ np.where( Ms != 0 ) ].min()
print Ms.max()

figure(figsize=(8.5,11) )
Dats = [ FF, UF, FU, UU ]
Titles = [ 'F to F (N=%.2E)'% FF[np.where( FF!=0 )].shape[0], 'U to F (N=%.2E)'% UF[np.where( UF!=0 )].shape[0], 'F to U (N=%.2E)'% FU[np.where( FU!=0 )].shape[0], 'U to U (N=%.2E)'% UU[np.where( UU!=0 )].shape[0] ]


for i in range(1,5):
   subplot( 220 + i )
   dat = Dats[i-1]
   hist( dat[np.where( dat!=0 )], bins=10**np.linspace(1,8), range=[10,10**8], normed=True )
   xscale('log')
   title( Titles[i-1] )
   xlabel('MFPT (frames)')
   yticks([])
   ylabel('Frequency (Linear Scale)')
#subplot(221)
#hist(FF[np.where( FF != 0 )], bins=10**np.linspace(1,8), range=[10,10**8],normed=True)#,histtype='step')
#xscale('log')
#title('F to F (N=%.2e)'% FF[np.where( FF!=0 )].shape[0])
##text(10**6, ylim()[1]*0.8,'F to F')
#yticks([])
#subplot(223)
#hist(FU[np.where( FU != 0 )], bins=10**np.linspace(1,8), range=[10,10**8], normed=True)#,histtype='step')
#xscale('log')
#title('F to U (N=%.2e)'% FU[np.where( FU!=0 )].shape[0])
##text(10**6, ylim()[1]*0.8,'F to U')
#yticks([])
#subplot(222)
#hist(UF[np.where( UF != 0 )], bins=10**np.linspace(1,8), range=[10,10**8], normed=True)#,histtype='step')
#xscale('log')
#yticks([])
#title('U to F (N=%.2e)'% UF[np.where( UF!=0 )].shape[0])
#text(10**6, ylim()[1]*0.8,'U to F')
#subplot(224)
#hist(UU[np.where( UU != 0 )], bins=10**np.linspace(1,8), range=[10,10**8], normed=True)#,histtype='step')
#xscale('log')
#yticks([])
#title('U to U (N=%.2e)'% UU[np.where( UU!=0 )].shape[0])
##text(10**6, ylim()[1]*0.8,'U to U')
#
savefig(options.out_FN)

