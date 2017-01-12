from __future__ import division
import _harmonic
import numpy as np
import matplotlib.pyplot as plt


def rms(residuals, m, cov=None):
    n = len(residuals)
    if cov is not None:
        w = np.diag(cov)
    else:
        w = np.ones(n)
    chisq = np.dot(residuals**2, w)

    nrms = np.sqrt(chisq / (n - m))
    wrms = np.sqrt((n / (n - m)) * chisq / sum(w))
    return nrms, wrms


def aic(m, n, std):
    return 2 * m + np.log10(std)


def bic(m, n, std):
    return np.log10(std) + np.log10(n) * m


def harmoinc_fit(t, x, polys=1, periods=None, pdef=None, cov=None, timebegin=None, max_order=10):
    if timebegin is None:
        timebegin = t[0]
    harmonic = _harmonic.harmonic
    harmonic.t = t
    harmonic.x = x
    harmonic.residuals = np.zeros(len(t))
    harmonic.periods = periods
    if polys < 0:
        aics = []
        # bics = []
        for i in xrange(max_order):
            harmonic.gen_design_matrix(i)
            harmonic.fit()
            aics.append(aic(harmonic.num_par, harmonic.n, harmonic.std))
            # bics.append(bic(harmonic.num_par, harmonic.n, harmonic.std))
            print i
        polys = aics.index(min(aics))
        # print bics.index(min(bics))
    # harmonic.design_matrix = None

    harmonic.gen_design_matrix(polys)
    # print harmonic.design_matrix
    # for i in range(harmonic.num_par):
    #     print harmonic.design_matrix[:, i]
    harmonic.fit()
    # print harmonic.out, harmonic.std
    # print np.sqrt(np.diag(harmonic.covar))
    # print rms(harmonic.residuals, harmonic.num_par, harmonic.cov)

t = np.arange(5000)
x = 0.001 * t + 0.3 + np.sin(2 * np.pi * t / 365) + np.cos(2 * np.pi * t / 365) + \
    1.2 * np.sin(2 * np.pi * t / 170) + np.cos(2 * np.pi * t / 170)
import cProfile
cProfile.run("harmoinc_fit(t, x, polys=-1, periods=[365])", sort=1)
# for i in range(1, 3):
#     print 'true, sin, i = ', i, np.sin(2 * np.pi * t / i)
#     print 'true, cos, i = ', i, np.cos(2 * np.pi * t / i)
# print _harmonic.harmonic.aic.__doc__
