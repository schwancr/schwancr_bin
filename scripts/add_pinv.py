#!/home/schwancr/Installed/epd/bin/python -u
import tables
import sys
import numpy as np

COMPRESSION = tables.Filters(complevel=9, complib='blosc', shuffle=True)

f = tables.openFile(sys.argv[1], 'a')

k = f.root.K[:]

print "calculating pseudoinverse"

pk = np.linalg.pinv(k, 1E-10)

n=k.shape[0] / 2
r = np.zeros(k.shape)
r[:n, n:] = np.eye(n)
r[n:, :n] = np.eye(n)

print "constructing matrix equation"
lhs = pk.dot(pk).dot(k).dot(r).dot(k)

print "calculating eigenvalues"
vals = np.linalg.eigvals(lhs)
print vals

print "adding everything to file"
atom = tables.Atom.from_dtype(lhs.dtype)
node = f.createCArray(where='/', name='pinv_lhs',
        atom=atom, shape=lhs.shape, filters=COMPRESSION)
node[:] = lhs

atom = tables.Atom.from_dtype(vals.dtype)
node = f.createCArray(where='/', name='pinv_vals',
        atom=atom, shape=vals.shape, filters=COMPRESSION)
node[:] = vals

atom = tables.Atom.from_dtype(pk.dtype)
node = f.createCArray(where='/', name='pinv',
        atom=atom, shape=pk.shape, filters=COMPRESSION)
node[:] = pk

f.flush()
f.close()
