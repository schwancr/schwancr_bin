#!/usr/bin/env python 

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-q', dest='qdist_FN',default='./Assignments.h5.QDists', help='Assignments.h5.QDist from AssignQDist.py [ ./Assignments.h5.QDists ]')
parser.add_option('-r', dest='cutoff',type=float, help= 'Cutoff used in clustering')

options, args = parser.parse_args()

from numpy import *
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from msmbuilder import Serializer

qdist = Serializer.LoadData( options.qdist_FN )

figure()

hist( qdist[ where( qdist >= 0 ) ].flatten(), bins=100, range=[0,1] )

xlabel('QDist to Generator')
ylabel('Frequency')
xlim([0,1])

if options.cutoff:
	vlines( options.cutoff, 0, ylim()[1], color='red' )

title('QDist Distribution')

savefig('QDist_Dist.pdf')

