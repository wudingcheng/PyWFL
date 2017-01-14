# encoding: utf-8
import _harmonics


def cs2k(c, s, quantum=True):
    """Change real spherical harmonic coefficients to complex spherical harmonic coefficients.

    Args:
        c (array) : c[n, n], Coefficients of cosine term
        s (array) : s[n, n], Coefficients of sine term
        quantum (bool) : default True, then the subroutine will return the fully normalized complex coefficients usually used in quantum mechanics domain; if quantum is False, then the subroutine will return the complex coefficients corresponding to c and s. See algorithm for details.

    Returns:
        array : k[n, 0:2n-1], Complex coefficients of c and s

    Note:
        here c and s are get from subroutine FUNC_EXP_FFT

        if quantum is False, then:

            .. math::
                \\left ( 2 - \\delta_m \\right ) K_{n,m} = C_{n, m} - i S_{n, m} \\quad m \\geq 0, \\left\\{\\begin{matrix} \\delta = 1.0, m \\neq 0 \\\\ \\delta = 2.0, m = 0  \\end{matrix}\\right. \\\\
                K_{n, -m} = (-1)^m conjg(K_{n, m}) \\quad m > 0


        if quantum is True, then:

            .. math::
                K_{n,m} & = & -1^m \\sqrt{4 \\pi} \\left (c_{nm} + is_{nm} \\right ) /  \\sqrt{2 - \\delta_m} \\quad m \\geq 0 \\\\
                k_{n,-m} & = & -1^m conjgK_{n, m}, \\quad m > 0

    Example:
        >>> import numpy as np
        >>> c = np.array([[1.2, 0, 0], [2.4, 3.2, 0], [0.9, 0.8, 9.8]])
        >>> s = np.zeros((3, 3))
        >>> cs2k(c, s).shape
"""
    cha = 'quantum' if quantum else ''
    return _harmonics.py_cs2k(c, s, cha=cha)


def cs2kai(cs, cha='x'):
    """Change degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20,C21,S21) to polar motion and length of day (LOD) domain

    Args:
        cs (float) : One of the degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20 or C21 or S21)  (cs is fully normalized spherical harmonic coefficients, like get from :func:`func_exp_fft`)
        cha (string) : cha is one of follow character, 'X',  'Y' , 'Z','x','y','z'

    Returns:
        float: Polar Motion (unit: mas) or Length of Day (unit: ms)

    Notes
        the transfer pairs are (C20, LOD), (C21, POLAR MOTION X),(S21, POLAR MOTION Y),  and only mass term suit to use this function.

        mas=milli-arcsencod, ms=millisecond

    References:
        .. [1] Yan Haoming, The role of ocean in Earth's Rotation and Gravity field, Doctor thesis, Wuhan, 2005.
"""
    conds = ['x', 'y', 'z']
    if cha.lower() not in conds:
        raise ValueError("cha should be in [{}]".format(' '.join(conds)))
    return _harmonics.py_cs2kai(cs, cha)


def pn(n, x):
    """Compute the Legendre Polynomials Pn(x) using recursion relation equation

    Args:
        n (int) : Maximum degree of Legendre polynomials :math:`P_{n}(x)`
        x (float) : In general, :math:`x = \cos \theta`, :math:`|x| \leq 1`

    Returns:
        array : pnx[n], Value of :math:`P_{n}(x)`

    Notes:
        http://www.efunda.com/math/legendre/index.cfm for more details.
        This subroutine is modified from a Fortran77 function provided by Dr. Wang Hansheng (汪汉胜).

    Examples:
        >>> from numpy import cos, pi
        >>> pn(5, cos(pi / 3))

    See Also:
        pn_diff
"""
    return _harmonics.py_wfl_pn(n, x)


def pn_diff(n, x):
    """Compute the one order derivative of Legendre Polynomials Pn(x) using recursion relation equation. The results are :math:`dP_{n}(\\cos \\theta / d\\theta)`, where :math:`x = \\cos \\theta`, :math:`-1 \\ leq x \\leq 1`

    Args:
        n (int) : Maximum degree of Legendre polynomials :math:`P_{n}(x)`
        x (float) : In general, :math:`x = \\cos \\theta`, :math:`|x| \\leq 1`

    Returns:
        array : pndiff[n], Value of :math:`dP_{n}(x) / d\\theta`

    Notes:
        http://www.efunda.com/math/legendre/index.cfm for more details.
        This subroutine is modified from a Fortran77 function provided by Dr. Wang Hansheng (汪汉胜).

    Examples:
        >>> from numpy import cos, pi
        >>> pn_diff(5, cos(pi / 3))

    See Also:
        pn

"""
    return _harmonics.py_wfl_pn_diff(n, x)


def pnm(n, theta):
    """Get fully normalized associated Legendre function :math:`P_{nm}(\\cos \\theta)`

    Args:
        n (int) : Maximum order number of fully normalized associated Legendre function (:math:`n \geq 1`),
        theta (float) : :math:`\\theta`, co-latitude in degree unit (:math:`0^{\circ}-180^{\circ}`)

    Returns:
        ndarray : p[n+1,n+1], Fully normalized (:math:`4\pi`) associated Legendre function :math:`P_{nm}(\cos \\theta)`


    Notes:
        The best maximum order n for PNM is 1080. If n>1080, the error will increase. Check precision formula is:

            .. math:: \\sum_{m=0, n} \\left ( P_{nm} \\right)^2 = 2n + 1

    Examples:
        >>> from numpy import pi, sin, cos, sqrt
        >>> theta = pi / 3
        >>> pnm(2, 60)
        >>> res = pnm(3, 60)
        >>> sqrt(5) * (3 * cos(pi / 3) ** 2 - 1) / 2 # res[2, 0]
        >>> sqrt(15) * cos(pi / 3) * sin(np.pi / 3) # res[2, 1]
"""
    return _harmonics.py_pnm(n, theta)


def pnmi(n, ts, tn):
    """Calculate :math:`\\int P_{nm}\\left (\\cos \\theta \\right ) \\sin \\theta d\\theta|_{ts}^{tn}`

    Args:
        n (int) : Maximum order number of fully normalized associated Legendre function (n>=1)
        ts (float) : Southern co-latitude (:math:`0^{\circ}-180^{\circ}`)
        tn (float) : Northern co-latitude (:math:`0^{\circ}-180^{\circ}`), ts>tn (such as 35°,30° which means latitude 55° and 60°)

    Returns:
        ndarray: pinmsn[n+1, n+1]

    Notes:
        The maximum order n for PNMI is 1080. If n>1080, the error will increase quickly. Check precision formula is:
            PNMI(ts,t) + PNMI(t,tn) = PNMI(ts,tn) where ts < t <tn

    Examples:
        >>> pnmi(3, 40, 30)

"""
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


def func_sum_fft(lon, lat, c, s, point=True, cmethod='p', depart_lon=0):
    """Get function point or mean value from spherical harmonic coefficients

    Args:
        lon (int) : Number of longitude (center grid number)
        lat (array) : Latitude coordinate (any begin value (such as 88.6 in the first grid) and regular grid after the first one (such as 86.1,83.6,... for 2.5 degree grid), when point is True;  must be grid center latitude when point is False) (unit: degree 90~-90, direction North-->South )
        c (array) : c[n], Coefficients of cosine term
        s (array) : s[n], Coefficients of sine term
        point (boolean) : default True, then output function point value at longitude and latitude intersection; when False, then output mean function value at each grid block. See comment and algorithm in below for details.
        cmethod (string): default 'p', means ... optional in ['p', 'dp']

    Returns:
        array : val[nlon, nlat], Function value at global grid points

    Notes:
        There's some differences for point value and mean value, see Andre Mainvile (1986) for details.
        When expand a global function to spherical harmonic coefficients,  the maximum degree is 180/grid intervals (e.g. for 1*1 degree grid, the maximum degree is 180, but for 2.5*2.5 degree grid, the maximum degree is 72). But from spherical harmonic coefficients to get global function values, we can meet 3 cases like below:

            1. the maximum degree = function grid (2*N=LON), e.g. from 72 degree coefficients to get 2.5*2.5 degree (LON=144) global function value
            2. the maximum degree < function grid (2*N<LON), e.g. from 72 degree coefficients to get 1.25*1.25 degree (LON=288) global function value
            3. the maximum degree > function grid (2*N>LON),  e.g. from 72 degree coefficients to get 5*5 degree (LON=72) global function value

        All the above case can be solved by FUNC_SUM_FFT . The first one is exact, but the last two cases interpolate global function values from coefficients domain, which likes to interpolate in frequency domain and then go back to temporal domain. So be careful to use the last two cases.

        All spherical harmonic coefficients are fully normalized (see FUNC_EXP for details).

        FUNC_SUM uses the popular method to calculate a function point value from spherical harmonic coefficients, but FUNC_SUM_FFT uses the most fast fourier transform method when calculate a function's grid value and mean value, so FUNC_SUM_FFT is more speed than FUNC_SUM. In my computer (2.0G CPU, 512Mb RAM, Windows2000, VISUAL FORTRAN 6.5), the speed of FUNC_SUM_FFT is about n (the maximum degree of coefficients) times quickly than FUNC_SUM and same precision for the coefficients. So the best choice is to use FUNC_SUM_FFT, but not FUNC_SUM.


        All coefficients are normalized (see FUNC_EXP for most comments):

            .. math:: \\begin{Bmatrix} \\bar{C}_{lm} \\\\ \\bar{S}_{lm} \\end{Bmatrix} = \\frac{1}{4\\pi} \\iint_{S} f(\\theta, \\lambda) \\bar{p}_{lm}(\\sin \\phi) \\begin{Bmatrix} \\cos m\\lambda \\\\ \\sin m\\lambda \\end{Bmatrix} ds

        For Point value:

            .. math:: F_{ij} = \\sum_{n=0}^{N} \\sum_{m=0}^{n} \\left( C_{nm} \\cos m \\lambda_j + S_{nm} \\sin m \\lambda_j \\right) \\bar{P}_nm(\\cos \\theta_i)

        For mean value:

            .. math::
                F_{ij} & = & \\sum_{n=0}^{N} \\sum_{m=0}^{n} \\left( C_{nm} \\cos m \\lambda_j + S_{nm} \\sin m \\lambda_j \\right) \\bar{P}_nm(\\cos \\theta_i) \\\\
                \\Delta_{ij} & = & \\Delta \\lambda \\left [ \\cos \\theta_i - \\cos \\theta_{i+1} \\right ] \\\\
                \\bar{I}_{nm}^{i}(\\theta) & = & \\int_{\\theta_i}^{\\theta_{i+1}}\\bar{P}_{nm}\\left ( \\cos \\theta \\right ) \\sin \\theta d\\theta \\\\
                J_{m}^{i}(\\lambda) & = & \\int_{\\lambda_j}^{\\lambda_{j+1}} \\cos m\\lambda d\\lambda = \\frac{1}{m}\\sin \\left ( m\\lambda \\right )|_{j\\Delta \\lambda}^{(j+1) \\Delta \\lambda} \\\\
                K_m^{i}(\\lambda) & = & \\int_{\\lambda_j}^{\\lambda_{j+1}}\\sin m\\lambda d\\lambda = \\frac{-1}{m} \\cos( m\\lambda) |_{j \\Delta \\lambda}^{(j+1)\\Delta \\lambda}

        Example:
            todo
"""
    cha = 'point' if point else 'mean'
    return _harmonics.py_func_sum_fft(lon, depart_lon, lat, c, s, cha=cha, cmethod=cmethod)


def func_sum(lon, lat, c, s, cmethod='p'):
    """Get function point value at grid intersection from harmonic coefficients

    Args:
        lon (array) : lon[nlon], Longitude coordinate at grid lattice or center or between both (unit: degree 0~360 )
        lat (array) : lat[nlat], Latitude coordinate at grid lattice or center or between both (unit: degree 90~-90 )
        c (array) : c[n], Coefficients of cosine term
        s (array) : s[n], Coefficients of sine term
        cmethod (string): default 'p', means ... optional in ['p', 'dp']

    Returns:
        array : val[nlon, nlat], Function value at global grid points

    Notes:
        The maximum order n for FUNC_SUM/FUNC_SUM_FFT is 1800. If n>1800, the error will increase quickly.

        func_sum use popular method (see algorithm for details), so the speed of func_sum is slow, you'd better to using func_sum_fft to get the same precision results but most fast in speed.

        All coefficients are normalized (see func_exp for most comments):

        .. math::
            \\begin{Bmatrix} \\bar{C}_{lm} \\\\ \\bar{S}_{lm} \\end{Bmatrix} = \\frac{1}{4\\pi} \\iint_{S} f(\\theta, \\lambda) \\bar{p}_{lm}(\\sin \\phi) \\begin{Bmatrix} \\cos m\\lambda \\\\ \\sin m\\lambda \\end{Bmatrix} ds \\\\
            F_{ij} = \\sum_{n=0}^{N} \\sum_{m=0}^{n} \\left( C_{nm} \\cos m \\lambda_j + S_{nm} \\sin m \\lambda_j \\right) \\bar{P}_nm(\\cos \\theta_i)

    See Also:
        func_sum_fft
"""
    return _harmonics.py_func_sum(lon, lat, c, s, cmethod=cmethod.upper())


def harmonic(t, x, periods, pdef=None, itrend=2, timebegin=None, method='aic'):
    if timebegin is None:
        timebegin = 0
    if pdef is None:
        return _harmonics.py_wfl_harmonic(t, x, periods, itrend=itrend, timebegin=timebegin, method=method)
    else:
        return _harmonics.py_wfl_harmonic_more(t, x, periods, pdef, itrend=itrend, timebegin=timebegin, method=method)
