
from scipy.linalg.lapack import get_lapack_funcs
from scipy.linalg.misc import LinAlgError
from scipy.linalg.decomp import _make_complex_eigvecs, _I
import numpy


def _geneig(a1, b1, left=False, right=True, overwrite_a=False, 
            overwrite_b=False, return_ab=True):
    ggev, = get_lapack_funcs(('ggev',), (a1, b1))
    cvl, cvr = left, right
    res = ggev(a1, b1, lwork=-1)
    lwork = res[-2][0].real.astype(numpy.int)
    if ggev.typecode in 'cz':
        alpha, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr, lwork,
                                                    overwrite_a, overwrite_b)
        w = alpha / beta
    else:
        alphar, alphai, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr, lwork,
                                                        overwrite_a,overwrite_b)
        w = (alphar + _I * alphai) / beta
        alpha = alphar + _I * alphai
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal ggev'
                                                                    % -info)
    if info > 0:
        raise LinAlgError("generalized eig algorithm did not converge (info=%d)"
                                                                    % info)

    only_real = numpy.logical_and.reduce(numpy.equal(w.imag, 0.0))
    if not (ggev.typecode in 'cz' or only_real):
        t = w.dtype.char
        if left:
            vl = _make_complex_eigvecs(w, vl, t)
        if right:
            vr = _make_complex_eigvecs(w, vr, t)
    if not (left or right):
        return w
    if left:
        if right:
            return w, vl, vr
        return w, vl

    if return_ab:
        return alpha, beta, w, vr

    return w, vr

