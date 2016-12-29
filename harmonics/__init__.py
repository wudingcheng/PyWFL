# from _harmonics import harmonics as _harmonics
from _harmonics import harmonics as _harmonics


def cs2kai(cs, cha):
    return _harmonics.py_cs2kai(cs, cha)


def pn(n, x):
    return _harmonics.py_wfl_pn(n, x)


def pn_diff(n, x):
    return _harmonics.py_wfl_pn_diff(n, x)


def pnm(n, theta):
    return _harmonics.py_pnm(n, theta)


def pnmi(n, ts, tn):
    return _harmonics.py_pnmi(n, ts, tn)


def func_exp_fft(func, cord_lon, cord_lat, n, smooth=False, iterative=1):
    """Expand a global function (such as ocean function) to spherical harmonic coefficients (with iterative refinements).

    Args:
        func (array) : Input function (values in grid center) North->South direction in latitude
        cord_lon (array) : Longitude coordinate at lattice (unit: degree 0~360 )
        cord_lat (array) : Latitude coordinate at lattice (unit: degree 90~-90),must be sort from north to 
        n (int) : Maximum spherical harmonic coefficients degree when expand function
        cha (bool) : default False, if True, then a smoothing operator is used for all the output spherical harmonic coefficients c and s. Using smoothing operator can get better coefficients approximation for the grid point value. If cha is False, then no smoothing operator is used. Then all the spherical harmonic coefficients c and s are the coefficients for mean block values (not grid point value) of the global function func.
        iterative (int) : default 1; iterative refinements of the results c and s.  If iterative=1, then do once iterative refinements. if iterative=0, then do not do iterative refinements. if iterative equal other number (such as k), then do k times iterative refinements. In general, iterative=1 is enough; but a second call (iterative=2) to verify convergence can be reassuring.

    Returns:
        tuple
            * c (*array*) : Coefficients of cosine term
            * s (*array*) : Coefficients of sine term
            * diff_percent (float) : Different percent(%) between integrals of f(q,l)2 in the sphere and the sum of K2 (complex harmonic coefficients get from c and s by using subroutine CS2K). If diff_percent is very small, it denotes the coefficients expanded from func are best.

    Notes:
        func_exp_fft has a high speed by using FFT method, and with the iterative refinements and using Parseval theorem to check how much precision of the spherical harmonic coefficients are expanded.
        The maximum order n for func_exp_fft is 1800. if n > 1800, the error will increase quickly.
        > There's some differences for spherical harmonic coefficients c and s between using smoothing operator and not using smoothing operator.
        When expand a global function to spherical harmonic coefficients,  the maximum degree is 180/grid intervals (e.g. for 1*1 degree grid, the maximum degree is 180, but for 2.5*2.5 degree grid, the maximum degree is 72). 
        All spherical harmonic coefficients are fully normalized (see algorithm for details). 

        The definition of the un-normalized Associated Legendre Polynomial is then

        .. math::
            P_{lm}(x) = \\frac{(1 - x^2)^{m/2}}{2^n n!} \\frac{d^{l+m}}{dx^{l+m}}[(x^2 - 1)^n]

        If the normalization factor is defined such that

        .. math::
            N_{lm}^2 = \\sqrt{\\frac{(2 - \\delta_{m0})(2l + 1)(l - m)!}{(l + m)!}}

        and the Associated Legendre Polynomials are normalized by

        .. math::
            \\bar{P}_{lm} = N_{lm} P_{lm}

        then, over a unit sphere :math:`S`

        .. math::
            \\iint_{S}\\left[ \\bar{P}_{lm}(\\sin \\varphi ) \\begin{Bmatrix} \\cos m\\lambda \\\\ \\sin m\\lambda \\end{Bmatrix}\\right]^2 ds = 4 \\pi

        In this convention, the relationship of the spherical harmonic coefficients to the mass distribution function f (or any function) becomes

        .. math::
            \\begin{Bmatrix} \\bar{C}_{lm} \\\\ \\bar{S}_{lm} \\end{Bmatrix}\\ = \\frac{1}{4\\pi} \\iint_S f(\\theta, \\lambda) \\bar{P}_{lm}(\\sin \\varphi)\\begin{Bmatrix} \\cos m\\lambda \\\\ \\sin m\\lambda \\end{Bmatrix}ds


    Reference:
        .. [1] Andre Mainvile,The Altimetry-Gravimetry Problem Using Orthonormal Base Functions, Report No. 373, Department of Geodetic Science and Surveying, the Ohio State University, 1986.
        .. [2] Paul M. K., Recurrence Relations For Integrals Of Associated Legendre Functions, Bulletin Geodesique, VOL.52, P177-190,1978.
        .. [3] Wunsch C, and D. Stammer, The global frequency-wavenumber spectrum of oceanic variability estimated from TOPEX/Poseidon altimetric measurements, J. Geophys. Res., Vol.100, NO. C12, P24895-24910, 1995.  (Eq. 2)

    Examples:
        todo


"""
    cha = 'smooth' if smooth else 'no'
    return _harmonics.py_wfl_func_exp_fft(func, cord_lon, cord_lat, n, cha=cha, iterative=1)


def dpnm(n, theta):
    """
    Get first order derivatives of fully normalized associated Legendre function :math:`dP_{nm}(\cos \\theta)/d\\theta`

    Args:
        n (int) : Maximum degree of fully normalized associated Legendre function (n>=1)
        theta (float) : Co-latitude in degree unit (0-180).

    Returns:
        array : First order derivatives of fully normalized (:math:`4\pi`) associated Legendre function :math:`dP_{nm}(\cos \\theta)/d\\theta` or :math:`P_{nm}^{'}`

    Notes:
        See FUNC_EXP for definition of fully normalized Legendre function.
        check value: see example for this subroutine checking.

    References:
        .. [1] Rummel R, van Gelderen M, Koop R, Schrama E, Sanso` F, Brovelli M, Miggliaccio F, Sacerdote F (1993) Spherical harmonic analysis of satellite gravity gradiometry. Publications on Geodesy, New series, no. 39. Netherlands Geodetic Commission, Delft.  Appendix A-2 (Note, some equations in A-2 have misprints).

    Examples:
        todo
"""

    return _harmonics.py_wfl_dpnm(n, theta)


def dpnm2(n, theta):
    """
    Get first and second order derivatives of fully normalized associated Legendre function :math:`dP_{nm}(\cos \\theta)/d\\theta`, :math:`d^{2}P_{nm}(\cos \\theta)/d\\theta^2)`

    Args:
        n (int) : Maximum degree of fully normalized associated Legendre function (n>=1)
        theta (float) : Co-latitude in degree unit (0-180).

    Returns:
        tuple
            * :math:`P_{nm}^{'}` (array) : First order derivatives of fully normalized (:math:`4\pi`) associated Legendre function :math:`dP_{nm}(\cos \\theta)/d\\theta` or :math:`P_{nm}^{'}`
            * :math:`P_{nm}^{''}` (array) : Second order derivatives of fully normalized (:math:`4\pi`) associated Legendre function :math:`dP_{nm}^{2}(\cos \\theta)/d\\theta^{2}` or :math:`P_{nm}^{''}`

    Notes:
        See FUNC_EXP for definition of fully normalized Legendre function.
        check value: see example for this subroutine checking.

    References:
        .. [1] Rummel R, van Gelderen M, Koop R, Schrama E, Sanso` F, Brovelli M, Miggliaccio F, Sacerdote F (1993) Spherical harmonic analysis of satellite gravity gradiometry. Publications on Geodesy, New series, no. 39. Netherlands Geodetic Commission, Delft.  Appendix A-2 (Note, some equations in A-2 have misprints).

    Examples:
        todo
"""
    return _harmonics.py_wfl_dpnm2(n, theta)

import numpy as np
# import matplotlib.pyplot as plt

data = np.loadtxt('data.txt')
data = data[:, 2].reshape(180, 360).T
lons = np.arange(0.5, 360, 1)
lats = np.arange(-90, 91, 1)[::-1]
# # print data[:, 10]
# # print data.shape, lats.shape
c, s, dif = func_exp_fft(data, lons, lats, 180, iterative=5)
# # plt.imshow(data)
# # plt.show()
print c.shape
print dif
# print c[:, 0]

# help(_harmonics.com)
