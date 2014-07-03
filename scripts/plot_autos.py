#!/usr/bin/env python

from msmbuilder import io
from msmbuilder import msm_analysis
from scipy.io import mmread
from schwancrtools import autocorrelate
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
matplotlib.rcParams['font.size'] = 22

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-t', dest='tProb', default='tProb.mtx')
parser.add_argument('-d', dest='raw_data')
parser.add_argument('-a', dest='avg_msm_data')
parser.add_argument('-o', dest='output', default='Autos.pdf', 
                    help='Output plot (pdf)')
parser.add_argument('-n', dest='num_modes', type=int, default=100)
parser.add_argument('-l', dest='lagtime', type=int, default=1,
                    help='Lag time in units of (--units)')
parser.add_argument('--divisor', dest='divisor', type=float, default=1,
                    help='Divisor to convert frames -> time units in the raw data.')
parser.add_argument('--units', dest='units', default='frames')

args = parser.parse_args()

try:  # should add more cases here...
    msm_data = np.loadtxt(args.avg_msm_data)[:,1]
except:
    print "Cannot load msm data."

tProb = mmread(args.tProb)
try:
    raw_data = io.loadh(args.raw_data)['arr_0']
except:
    raw_data = io.loadh(args.raw_data)['Data']

num_frames = raw_data.shape[1]

num_lagtimes = num_frames / args.lagtime

#msm_acf = msm_analysis.msm_acf(tProb, msm_data, np.arange(num_lagtimes), 
#                               num_modes=args.num_modes)

sampled_traj = msm_analysis.sample(tProb, np.random.randint(tProb.shape[0]), num_lagtimes)
data_traj = msm_data[sampled_traj]
msm_acf = autocorrelate.fft_autocorrelate(data_traj)

raw_acfs = []
for i in xrange(np.max([raw_data.shape[0], 10])):
    
    max_non_neg = np.where(raw_data[i] != -1)[0].max()
    row = raw_data[i][:max_non_neg + 1]

    raw_acfs.append(autocorrelate.fft_autocorrelate(row))

figure()
axes((0.18, 0.18, 0.72, 0.72))

raw_label = 'Raw Data'

for i in xrange(len(raw_acfs)):

    xi = np.arange(len(raw_acfs[i]), dtype=np.float) / args.divisor

    plot(xi, raw_acfs[i], color='blue', label=raw_label)

    raw_label=None

xi_msm = np.arange(len(msm_acf), dtype=np.float) * args.lagtime 

plot(xi_msm, msm_acf, color='red', alpha=0.8, lw=3, label='MSM')
legend()

xscale('symlog')
ylim([-.2,1.2])
xticks(xticks()[0][1:])
yticks(yticks()[0][1:])
xlabel('Time (%s)'%args.units)
ylabel('ACF')

savefig(args.output)

np.save(args.output[:-4] + '.npy', np.vstack((xi_msm, msm_acf)).T)

