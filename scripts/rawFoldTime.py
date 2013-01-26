#!/usr/bin/env python
 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p',dest='proj_FN',default='../ProjectInfo.h5',help='ProjectInfo.h5 from msmbuilder [ ../ProjectInfo.h5 ]')
parser.add_option('-d',dest='data_FN',help='Data to use as a metric for folded and unfolded states' )
parser.add_option('--fc', dest='f_cut',type=float,help='Folded cutoff')
parser.add_option('--uc',dest='u_cut',type=float,help='Unfolded cutoff')
parser.add_option('--low-is-folded',dest='low_is_folded',default=False,action='store_true',help='Pass this flag if a small number means the conformation is folded (i.e. RMSD)')
parser.add_option('-o',dest='out_FN',default='Fold_Unfold_Times.pdf',help='Output file to write to')

options, args = parser.parse_args()
from numpy import *
from msmbuilder import Project
from pyschwancr import dataIO, msmTools
import os, sys, re
import matplotlib
matplotlib.use('pdf')
from matplotlib.pyplot import *
from scipy import optimize
Proj = Project.Project.LoadFromHDF( options.proj_FN )
Data = dataIO.readData( options.data_FN )

# first reshape the data into trajectories.

Lens = Proj['TrajLengths']

Trajs = []
sum = 0
for i in range( len( Lens ) ):
	Trajs.append( Data[ sum : sum + Lens[i] ] )
	sum += Lens[i] 
Folds = []
Unfolds = []

for traj in Trajs:
	(a,b) = msmTools.calcRawFoldTime( traj, options.f_cut, options.u_cut, low_is_folded = options.low_is_folded )
	Folds.extend( a )
	Unfolds.extend( b )

#FoldsDist = bincount( Folds )
#UnfoldsDist = bincount( Unfolds )
figure()
subplot(211)
foldHist = hist( Folds, bins=100, color = 'blue', label='Fold' )
vlines( mean( Folds ), 0, ylim()[1], color = 'black', linewidth=3 )
ylabel('Frequency')
legend()
xFolds = xlim()

subplot(212)
unfoldHist = hist( Unfolds, bins=100, color = 'red', label='Unfold' )
vlines( mean( Unfolds), 0, ylim()[1], color = 'black', linewidth=3 )
ylabel('Frequency')
legend()
xUnfolds = xlim()

xlabel('Fold/Unfold Times (frames)')
suptitle('Distribution of Folding/Unfolding times')

subplot(211)
xlim([ 0, max( xFolds[1], xUnfolds[1] ) ])
text( xlim()[1] * 0.3, ylim()[1] * 0.8, 'Mean = %.2f\nN = %d' % ( mean( Folds ), len( Folds ) ) )
yLimF = ylim()
subplot(212)
xlim([ 0, max( xFolds[1], xUnfolds[1] ) ])
text( xlim()[1] * 0.3, ylim()[1] * 0.8, 'Mean = %.2f\nN = %d' % ( mean( Unfolds ), len( Unfolds ) ) )
yLimU = ylim()
savefig( options.out_FN )


yFold = foldHist[0]
xFold = array( [ ( foldHist[1][i+1] + foldHist[1][i] ) / 2. for i in range( len( foldHist[0] ) ) ] )

yUnfold = unfoldHist[0]
xUnfold = array( [ ( unfoldHist[1][i+1] + unfoldHist[1][i] ) / 2. for i in range( len( unfoldHist[0] ) ) ] )

expFit = lambda p, x : p[0] * exp( - p[1] * x )
powFit = lambda p, x : p[0] * x ** ( - p[1] )

errExp = lambda p, x, y : expFit( p, x ) - y
errPow = lambda p, x, y : powFit( p, x ) - y

foldExp = optimize.leastsq( errExp, x0 = [100,0.001], args = ( xFold, yFold ), maxfev = 1000000 )
foldPow = optimize.leastsq( errPow, x0 = [1,1], args = ( xFold, yFold ), maxfev = 1000000 )
unfoldExp = optimize.leastsq( errExp, x0 = [100,0.001], args = ( xUnfold, yUnfold ), maxfev = 1000000 )
unfoldPow = optimize.leastsq( errPow, x0 = [1,1], args = ( xUnfold, yUnfold ), maxfev = 1000000 )

SStot_F = ( ( yFold - yFold.mean() ) **2 ).sum()
SStot_U = ( ( yUnfold - yUnfold.mean() ) ** 2 ).sum()

SSerr_F_exp = ( ( yFold - expFit( foldExp[0], xFold ) ) ** 2 ).sum()
SSerr_F_pow = ( ( yFold - powFit( foldPow[0], xFold ) ) ** 2 ).sum()

SSerr_U_exp = ( ( yUnfold - expFit( unfoldExp[0], xUnfold ) ) ** 2 ).sum()
SSerr_U_pow = ( ( yUnfold - powFit( unfoldPow[0], xUnfold ) ) ** 2 ).sum()

R2_F_exp = 1 - SSerr_F_exp / SStot_F
R2_F_pow = 1 - SSerr_F_pow / SStot_F
R2_U_exp = 1 - SSerr_U_exp / SStot_U
R2_U_pow = 1 - SSerr_U_pow / SStot_U

figure()
xi = linspace( 1, max(xFolds[1], xUnfolds[1]), 1000 )
subplot(211)
scatter( xFold, yFold, color = 'blue', label='Fold Times' )
plot( xi, expFit( foldExp[0], xi ), color='purple', label='Exponential' )
plot( xi, powFit( foldPow[0], xi ), color='orange', label='Power Law' )
ylabel('Frequency')
xlim([ 0, max( xFolds[1], xUnfolds[1] ) ])
ylim( yLimF )
text(0.3*xlim()[1], ylim()[1]*0.7, u"Exp: R\xb2 = %.4f\nPow: R\xb2 = %.4f" % ( R2_F_exp, R2_F_pow ) )
legend()

subplot(212)
scatter( xUnfold, yUnfold, color = 'red', label = 'Unfold Times' )
plot( xi, expFit( unfoldExp[0], xi ), color='purple', label='Exponential' )
plot( xi, powFit( unfoldPow[0], xi ), color='orange', label='Power Law' )
ylabel('Frequency')
xlim([ 0, max( xFolds[1], xUnfolds[1] ) ])
ylim( yLimU )
text(0.3*xlim()[1], ylim()[1]*0.7, u"Exp: R\xb2 = %.4f\nPow: R\xb2 = %.4f" % ( R2_U_exp, R2_U_pow ) )
legend()

suptitle('Fits of Distribution of Folding/Unfolding Times')
xlabel('Fold/Unfold Times (frames)')

savefig( options.out_FN[:-4] + 'FITS' + options.out_FN[-4:] )
