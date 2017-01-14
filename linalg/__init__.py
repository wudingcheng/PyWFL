# encoding: utf-8
from __future__ import absolute_import
from . import _linalg
from numpy import iscomplexobj, real


__all__ = ["svd"]


def svd(a):
    """Singular Value Decomposition (SVD)

    Args:
        a (array) : (m, n), can be real or comple

    Returns:
        tuple
            * u (array): (m, n)
            * s (array): (n)
            * v (array): (n, n)

    Notes:
        For the real geophysical data, to get a more reliable results, before do the SVD analysis, please weighted the original data by square root area weight of sqrt(cosine(latitude)) first.  See also sqrootweight(nx) parameters in subroutine :func:`ceof`.
        The original code for the decomposing real/double matrix A to U, W and V is based on "Numerical Recipes for Fortran77", while for the complex/double complex matrix A is based on LINPACK subroutines "csvdc" and "zsvdc".

        For real matrix A:

            .. math:: A = U W V^{*}

        For comple matrix A:

            .. math:: A = U W V^{*}

        here :math:`V^{*}` means conjugate(V)




"""
    if iscomplexobj(a):
        u, s, v = _linalg.py_wfl_svd_c(a)
    else:
        u, s, v = _linalg.py_wfl_svd_d(a)
    s = real(s)
    return u, s, v
