#!/usr/bin/env python
 
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--a1',dest='ass_FN1',help='First assignments file')
parser.add_argument('--a2',dest='ass_FN2',help='Second assignments file')
parser.add_argument('--norm',dest='norm',type=int,default=0,choices=(0,1,2),help='Normalize to the number of assignments in a particular decomposition. For instance if m is the matrix whos i,j\'th entry corresponds to the conformations assigned to ass1 = i and ass2 = j, then passing --norm 1 will normalize to the states in ass1, and so the rows will be normalized to 1.')
parser.add_argument('-o',dest='out_FN',default='CompareAssignments.pdf',help='Output filename to save the plot to (pdf)')
parser.add_argument('--l1',dest='label1',default='First Assignments', help='Label for first assignments file')
parser.add_argument('--l2',dest='label2',default='Second Assignments', help='Label for second assignments file')

args = parser.parse_args()
 
import numpy as np
from msmbuilder import Serializer
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
import os, sys, re

ass1 = Serializer.LoadData( args.ass_FN1 )
ass2 = Serializer.LoadData( args.ass_FN2 )

num_states1 = ass1.max()+1
num_states2 = ass2.max()+1

if (num_states1 > 20) or (num_states2 > 20):
   print "This matrix is going to be big... (%d,%d), so the picture will likely look bad..." % (num_states1,num_states2)

if ass1.shape != ass2.shape:
   print "The two assignment files do not have the same shape... %s vs %s ...exiting." % ( str(ass1.shape), str(ass2.shape) )

non_neg_ind = np.where( (ass1 >= 0) & (ass2 >= 0) )

ass1 = ass1[non_neg_ind]
ass2 = ass2[non_neg_ind]

mat = np.zeros( (num_states1,num_states2) ).astype(int) 

for (i,j) in zip( ass1, ass2 ):
   mat[i,j] += 1

figure( figsize=(8,6))


if args.norm==0:
   pass
elif args.norm==1:
   mat = mat.astype(float)
   mat /= mat.sum(axis=1).reshape( (-1,1) ) # Reshape so that it divides each row
elif args.norm==2:
   mat = mat.astype(float)
   mat /= mat.sum(axis=0) # Automatically casts to each column
else:
   print "I shouldn't have gotten here..."

maxLog = int( np.log10( mat.max() ) + 1 )

if args.norm==0:
   matshow( np.log10(mat), cmap='hot_r', vmin=0,vmax=maxLog, fignum=1 )
else:
   matshow( np.log10(mat), cmap='hot_r', vmin=-4,vmax=0, fignum=1 )

c = colorbar()

if args.norm==0:
   c.set_ticks( range( maxLog + 1 ) )
   c.set_ticklabels( [ r'$10^{%d}$'%i for i in range( maxLog+1 ) ] )
else:
   c.set_ticks( range( -4, 1 ) )
   c.set_ticklabels( [ r'$10^{%d}$'%i for i in range( -4, 1 ) ] )
ax = gca()

ax.xaxis.set_label_position('top')

xlabel( args.label2 )
ylabel( args.label1 )

if num_states1 <= 20:
   yticks( range( num_states1 ) )

if num_states2 <= 20:
   xticks( range( num_states2 ) )

cbar_label = 'Number of Assignments'

if args.norm == 1:
   cbar_label += ' -- Row Normalized'
elif args.norm == 2:
   cbar_label += ' -- Column Normalized'

c.set_label(cbar_label)

savefig( args.out_FN )

