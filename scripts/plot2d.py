
import numpy as np
from matplotlib.pyplot import *
import argparse
matplotlib.rcParams['font.size'] = 22
matplotlib.rcParams['mathtext.default'] = 'regular'

parser = argparse.ArgumentParser()
parser.add_argument('-x', dest='xfn', help='data for x-axis')
parser.add_argument('-y', dest='yfn', help='data for y-axis')
parser.add_argument('-p', dest='populations', default='Populations.dat', help='populations from msm')
parser.add_argument('-o', dest='output', default='hist.pdf', help='output filename')
parser.add_argument('--xl', dest='xlim', type=float, help='xlimits', nargs=2)
parser.add_argument('--yl', dest='ylim', type=float, help='ylimits', nargs=2)
parser.add_argument('--xlabel', dest='xlabel', help='xlabel')
parser.add_argument('--ylabel', dest='ylabel', help='xlabel')

args = parser.parse_args()

class KDE2D:
    def __init__(self, means, stdevs, weights=None):
        self.means = means
        self.inv_stdevs = 1.0 / stdevs
        self.inv_stdevs2 = 1.0 / np.square(stdevs)

        if weights is None:
            self.weights = np.ones(self.means.shape)

        else:
            self.weights = weights

        d = self.means.shape[1]
        d = float(d)

        prefactors = 1 / np.power( 2. * np.pi, d / 2.0) * np.sqrt(np.product(self.inv_stdevs2, axis=1))
        self.weights = self.weights * prefactors


    def __call__(self, points):
        addl_dims = len(points.shape) - 1
        addl_axes = [1] * addl_dims
        new_shape = list(self.means[0].shape) + addl_axes
        dist2 = np.array([np.square(mu.reshape(new_shape) - points) for mu in self.means])

        new_shape = list(self.inv_stdevs.shape) + addl_axes
        term = - 0.5 * np.sum(dist2 * self.inv_stdevs2.reshape(new_shape), axis=1)

        new_shape = list(self.weights.shape) + addl_axes
        result = np.sum(np.exp(term) * self.weights.reshape(new_shape), axis=0)

        return result

x = np.loadtxt(args.xfn)
y = np.loadtxt(args.yfn)
pops = np.loadtxt(args.populations)

means = np.hstack([x[:, 1:2], y[:, 1:2]])
stdevs = np.hstack([x[:, 2:], y[:, 2:]])

kde = KDE2D(means, stdevs, pops)

xtest, ytest = np.meshgrid(np.linspace(args.xlim[0], args.xlim[1], 100), np.linspace(args.ylim[0], args.ylim[1], 100))

points = np.array([xtest, ytest])

P = kde(points)

axes((0.18, 0.18, 0.72, 0.72))
imshow(P / np.abs(P).max(), extent=args.xlim + args.ylim, origin='lower', 
        aspect='auto', interpolation='nearest', cmap=matplotlib.cm.RdBu,
        vmin=-1, vmax=1)

c = colorbar()
xlabel(args.xlabel)
ylabel(args.ylabel)

savefig(args.output)
