import _geophysics


def am2eam(am, cha='xp'):
    """ Change Angular Momentum (AM) to Effective Angular Momentum (EAM)

    Args:
        am (float) : Angular Momentum (unit: :math:`kg.m^2/s`)
        cha (string) : cha eaqual one of follow character, 'xp', 'xc', 'yp', 'yc', 'zp', 'zc'

    Returns:
        float : Effective Angular Momentum (unit: 1.0e-7 rad)

    Notes:a
        cha denote x,y or z coordinate for pressure (p) and current (c) term
        All the parameters about Earth is Eubank's value like below:

            >>> cc  = 7.1236 #earth principal moment of inertia
            >>> ccm = 7.1236
            >>> ca  = 0.003664*cc       #c - a

        barnes's earth parameter (not used here) :

            >>> cc=7.04
            >>> ccm=7.1236
            >>> ca=0.00333*cc

        this function only used in earth rotation domain, which change oceanic am to oceanic eam.

    References:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.

    See Also:
        todo

"""
    cond = ['xp', 'xc', 'yp', 'yc', 'zp', 'zc']
    if cha.lower() not in cond:
        raise ValueError("should should be in [{}]".format(' '.join(cond)))
    return _geophysics.py_am2eam(am, cha=cha.upper())


def geo_distance(phi1, lambda1, phi2, lambda2, ch='angle'):
    """Calculate geocentric angle (or corresponding Earth's surface distance) between two points

    Args:
        phi1 (float) : First point's coordinate latitude, unit: degree, -90~90
        lambda1 (float) : First point's coordinate longitude, unit: degree, 0~360
        phi2 (float) : Second point's coordinate latitude, unit: degree, -90~90
        lambda2 (float) : Second point's coordinate longitude unit: degree, 0~360
        ch (string) : default angle, should be in ['angle', 'distance'], unit radian or meter

    Returns:
        float : Geocentric angle (unit: radian) or corresponding Earth's surface distance (unit: meter).

    Notes:
        .. math::
            \\cos(\\phi) = \\cos( \\theta_1) * \\cos(\\theta_2) + \\sin(\\theta_1) * \\sin(\\theta_2) * \\cos(\\lambda_1 - \\lambda_2)

        :math:`\\phi` is geocentric angle between two points, and :math:`\\theta_1`,:math:`\\theta_2` are colatitude of :math:`\\phi_1` and :math:`\\phi_2` respectively.
        The Earth's radius adopted here is 6371012m.

    Examples:
        todo

    See Also:
        todo
"""

    cond = ['angle', 'distance']
    if ch.lower() not in cond:
        raise ValueError("ch should be in [{}]".format(" ".join(cond)))
    return _geophysics.geo_distance(phi1, lambda1, phi2, lambda2, ch)


def area(dlon, lat1, lat2):
    """Calculate the grid area of the Earth's surface

    Args:
        dlon (float) :Longitude interval (or length) of the grid (unit: degree)
        lat1 (float) :Co-latitude of northern grid border (0 represents Northern Pole, and 180  represents Southern Pole, unit: degree)
        lat2 (float) :Co-latitude of southern grid border (lat2>lat1)

    Return:
        float : Area of grid at the Earth's surface (unit: m**2)

    Notes:
        .. math::  ds=r^2 \\sin( \\theta )d( \\theta )d(\\lambda)

    Examples:
        todo

    See Also
        todo
"""

    lat1, lat2 = sorted([lat1, lat2])
    return _geophysics.py_area(dlon, lat1, lat2)


def cs2kai(cs, cha='x'):
    """ Change degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20,C21,S21) to polar motion and length of day (LOD) domain

    Args:
        cs (float) : One of the degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20 or C21 or S21)  (cs is fully normalized spherical harmonic coefficients, like get from FUNC_EXP_FFT)
        cha (string) : cha is one of follow character, 'x', 'y', 'z'

    Return:
        float : Polar Motion (unit: mas) or Length of Day (unit: ms)

    Notes:
        the transfer pairs are (C20, LOD), (C21, POLAR MOTION X),(S21, POLAR MOTION Y), and only mass term suit to use this function.
        mas=milli-arcsencod, ms=millisecond

    References:
        .. [1] Yan Haoming, The role of ocean in Earth's Rotation and Gravity field, Doctor thesis, Wuhan, 2005.


    Examples:
        todo

    See Also:
        kai2cs
"""
    cond = ['x', 'y', 'z']
    if cha.lower() not in cond:
        raise ValueError("should should be in [{}]".format(' '.join(cond)))
    return _geophysics.py_cs2kai(cs, cha)


def kai2cs(kai, cha='x'):
    """Change polar motion and length of day (LOD) domain to degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20,C21,S21)


    Args:
        kai (float) : Polar Motion (unit: mas) or Length of Day (unit: ms)
        cha (string) : cha is one of follow character, 'x', 'y', 'z'

    Return:
        float : One of the degree 2 spherical harmonic (Stokes) coefficients of the gravity field (C20 or C21 or S21)  (c20,c21,s21 are fully normalized spherical coefficients; the fully normalized (4*pi) Legendre polynomials)

    Notes:
        The transfer pairs are (LOD,C20 ), ( POLAR MOTION X,C21),(POLAR MOTION Y, S21), and only mass term are suit to use this function.
        mas=milli-arcsencod, ms=millisecond
        when convert LOD to c20, we don't consider the global mass change (in other words, we set :math:`\Delta M=0`)

    References:
        .. [1] Yan Haoming, The role of ocean in Earth's Rotation and Gravity field, Doctor thesis, Wuhan, 2005.


    Examples:
        todo

    See Also:
        cs2kai
"""
    cond = ['x', 'y', 'z']
    if cha.lower() not in cond:
        raise ValueError("should should be in [{}]".format(' '.join(cond)))
    return _geophysics.py_kai2cs(kai, cha)


def degree_var(c, s):
    """Get degree variance from spherical harmonic coefficients

    Args:
        c (array) : c(0:n,0:n) Spherical harmonic coefficients (coefficients for cosine term)
        s (array) : s(0:n,0:n) Spherical harmonic coefficients (coefficients for sine term)

    Returns:
        array : degvar(0:n) Degree variance

    Notes:
        .. math::
            degvar(n)=\sum(c(n,m)^2+s(n,m)^2) \qquad m=0,1,2,....,n \quad n=0,1,2, \dots, N

        where :math:`N` the maximum degree

    Examples:
        todo

    See Also:
        todo
"""
    return _geophysics.py_degree_var(c, s)


def j2j3(zlat, dlat, dlon, hi, rho=1000):
    """Calculate j2 and j3 at special grid point

    Args:
        zlat (float): Co-latitude at grid center (unit: degree)
        dlat (float): Interval of latitude (or grid length in latitude direction) (unit: degree)
        dlon (float): Interval of longitude (or grid length in longitude direction) (unit: degree)
        hi (float): Height variation of mass---usually use equivalent water height (unit: meter)
        rho (float): Density of mass--usually water density (1000kg/m^3)  (unit: kg/m^3)

    Returns:
        tuple
            * dj2 (float) : J2 value
            * dj3 (float) : J3 value

    Notes:
        here J2=-C(2,0), and C(2,0) is not fully normalized spherical coefficient, C(2,0)=sqrt(5)*C(2,0)(normalized)

    References:
        .. [1] Dickey J O, et al., Recent Earth Oblateness variations: Unraveling climate and postglacial rebound effects,  Science, Vol. 298, p1975-1977, 2002.

        Please see supplement material of the above paper for the  Equation to calculate J2.

    Examples:
        todo

    See Also:
        todo
"""
    return _geophysics.j2j3(zlat, dlat, dlon, hi, rho)


def eaam_xyz_j2j3(ltrop, slp, pl, u, v, z, hm, method):
    """
    Todo:
        todo
"""
    return _geophysics.py_eaam_xyz_j2j3(ltrop, slp, pl, u, v, z, hm, method)


def load_vert_def():
    return


def load_vert_exp(p):
    """Get the radial elastic deformation of the Earth's surface only estimated from atmospheric surface pressure field globally by experiment method.

    Args:
        p (array) : p[nlon,nlat], atmospheric surface pressure field (sorted by latitude (90~-90) first, then sorted by longitude (0~360). e.g. the coordinate of SIGMA should be (0.5,89.5) (1.5, 89.5),....(359.5,89.5) (0.5, 88.5) (1.5, 88.5).......... (359.5, -89.5) for 1x1 P). Unit of  P is Pascal. (Note: in the ocean, P must be set to 0.0d0)

    Returns:
        tuple
            * u[nlon, nlat]: Radial elastic defomation of the Earth's surface calculated by using experimental method(unit: mm) (grid mean value) , see comment for details. (Note: only 80N~60S have the values, other areas are 0)
            * pbar(nlon, nlat): The average pressure for each grid center point (surrounding area of the point is about 2000km wide), unit: Pascal. (Note: only 80N~60S have the values, other areas are 0)

    Notes:
        :math:`u=-0.9*pbar-0.35(p-pbar)`, here :math:`pbar` is the average (surrounding area of the point is 2000km wide) pressure for the point, and  p is the grid point surface pressure. See reference for details. This result is only the experimental value, the more precise value can be get by using LOAD_VERT_DEF.

    References:
        .. [1] Rabbel, W., J. Zschau, Static deformations and gravity changes at the Earth's surface due to atmospheric loading, J Geophys., 56, 81-99,1985.

"""
    return _geophysics.py_load_vert_exp(p)


def love(n, ch='wang'):
    """ Get load love number from 0 to n (n<=MAX) by cubic spline interpolation method

    Args:
        n (int) : Maximum order of load love. (if n>MAX, then all h(MAX:n),l(MAX:n), and k(MAX:n) will be infinity value, i.e. values of :math:`h_{\infty}`, :math:`l_{\infty}` and :math:`k_{\infty}`)
        ch (string) : optional in ['wang', 'farrell'], for 'wang', MAX=13000, for 'farrell', MAX=10000

    Returns:
        tuple
            * h (array): h(n), Load love number h
            * l (array): l(n), Load love number l
            * k (array): k(n), Load love number k

    References:
        .. [1] FARRELL W. E., Reviews of geophysics and space physics, 10,3,761-797, 1972. Table A2.
        .. [2] Wang Hansheng, etc. New evolve of load love number for SNREI Earth model, Chinese J. Geophys.,vol.36,supp.,182-189,1996 (in Chinese)
"""
    h, l, k = _geophysics.py_love(n, ch=ch)
    return h, l, k


def cs2cs(c, s, cha='g2h'):
    """Change geopotential spherical harmonic coefficients to equivalent water height spherical harmonic coefficients or vice versa

    Args:
        c (array) : c(0:n, 0:n) The spherical harmonic coefficients for geopotential (unit:m? corresponding to geoid height?) or equivalent water height (unit: m)
        s (array) : s(0:n, 0:n) The spherical harmonic coefficients for geopotential (unit:m? corresponding to geoid height?) or equivalent water height (unit: m)
        cha (string) : should be in ['g2h', 'h2g'], default g2h, which means change geopotential spherical harmonic coefficients to equivalent water height spherical harmonic coefficients, if cha='h2g', then convert the equivalent water height spherical harmonic coefficients to geopotential spherical harmonic coefficients

    Returns:
        tuple
            * c (array)
            * s (array)

    Notes:
        here we get the spherical harmonic coefficient as expansion of:

            .. math::
                dh = \\sum\\sum\\left\\{ P_{lm}\\left(\\cos \\theta \\right) \\left[ C_{lm}\\cos(m \\phi) + S_{lm}\\sin(m \\phi) \\right ] \\right \\}

        which is different with surface density defined by Wahr, 1998. The difference is the factor earth semi-radius and water density.

    References:
        .. [1] Wahr J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.


"""
    cond = ['g2h', 'h2g']
    if cha.lower() not in cond:
        raise ValueError("should should be in [{}]".format(' '.join(cond)))
    return _geophysics.py_wfl_cs2cs(c, s, cha=cha.upper())


def cs2vertical(c, s, depart_lon=-1, alat=None, water_density=1e3, frame='ce'):
    return _geophysics.py_wfl_cs2vertical(c, s, depart_lon=depart_lon, alat=alat,
                                          water_density=water_density, frame=frame)


def var_reduction(x, x_correction):
    return _geophysics.py_wfl_var_reduction(x, x_correction)
    # pass


def depth2pressure(depth, lat):
    """convert ocean depth to pressure

    Args:
        depth (float) : ocean depth (unit: meter)
        lat (float) : latitude in degree unit (:math:`-90^{\\circ} - 90^{\\circ})

    Return:
        float: ocean Pressure (unit: Pa).

    Notes:
        This pressure is gauge pressure under 1 atmospheric pressure, which means pressure=0 at sea surface

    References:
        .. [1] Saunders, "Practical Conversion of Pressure to Depth", J. Phys. Oceanog., April 1981.

    Examples:


"""
    return _geophysics.py_wfl_depth2pressure(depth, lat)


def pressure2depth(pressure, lat):
    """convert ocean pressure to depth

    Args:
        pressure (float) : ocean pressure (unit: Pa).
        lat (float) : latitude in degree unit (:math:`-90^{\\circ} - 90^{\\circ})

    Return:
        float : ocean depth (unit: meter)

    Notes:
        This pressure is gauge pressure under 1 atmospheric pressure, which means pressure=0 at sea surface

    References:
        .. [1] Saunders, "Practical Conversion of Pressure to Depth", J. Phys. Oceanog., April 1981.

"""
    return _geophysics.py_wfl_pressure2depth(pressure, lat)


def dm2lod(dm):
    """Calculate global mass change effect on the length of day (LOD) change

    Arg:
        dm (float) : gloabl mass change (unit: kg)

    Return:
        float : global mass balance effect on LOD change (unit: ms)

    References:
        .. [1] Yan H., and B. Chao, Global mass balance effect on Earth's rotation rate, Geophys. Res. Lett., 2010
"""
    return _geophysics.py_wfl_dm2lod(dm)


def fan_filter(nmax, mmax, r_lon, r_lat):
    """Get Gaussian Fan filtering weight function in frequency domain which is used to filter spherical harmonic coefficients

    Args:
        nmax (int): Maximum degree number of Gaussian filtering weight function
        mmax (int): Maximum order number of Gaussian filtering weight function (mmax<=nmax. If mmax<nmax, then all w(i,j)==0.0 for j>mmax)
        r_lon (float): r_lon is the filter average radius in longitude direction (unit: meter), which is the distance on the Earth's surface at which w has dropped to 1/2 it value at alpha=0 (the distance on the Earth's surface = Earth_R*alpha).
        r_lat (float): r_lat is the filter average radius in latitude direction (unit: meter).

    Returns:
        ndarray : w[nmax+1, nmax+1], Gaussian filtering weight function in frequency domain which is used to filter spherical harmonic coefficients

    Notes:
        Fan filter is a expansion of subroutine :func:`gauss_jekeli_fre`, which is the product of two :func:`gauss_jekeli_fre` weight in longitude and latitude direction.

        Here we use the equations of C. Jekeli. So when you use the equation defined by Wahr J. (1998), note the difference for factor 2*pi. E.G., in Wahr J. (1998), eq.(30) now can omit the factor 2*pi, you can get the same results.

        In the frequency domain, the Gaussian filtering weight function can be computed with recursion relations (see Algorithm in :func:`gauss_jekeli` for definition of a):

            .. math::
                W_0 & = 1 \\\\
                W_1 & = \\frac{1 + e^{-2a}}{1 - e^{-2a}} - \\frac{1}{a} \\\\
                W_{n+1} & = -\\frac{2n + 1}{a}W_{n} + W_{n-1}

        Now :math:`W_{ij}` can be multiplied to the spherical harmonic coefficients calculated from :func:`fun_exp_fft`, and then use :func:`func_sum_fft` to get the Gaussian filtering weighted spherical function value in space domain.

    References:
        .. [1] Zhang, Z., B.F. Chao, Y. Lu, and H.Hsu(2009), An effective filtering for GRACE time-variable gravity: Fan filter, Geophys. Res. Lett., 36,L17311, doi:10.1029/2009GL039459.
        .. [2] Jekeli, C., Alternative methods to smooth the Earth's gravity field, Rep.327, Dep. of Geod. Sci. and Surv., Ohio State Univ., Columbus,1981.
        .. [3] Wahr J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.


"""
    return _geophysics.py_wfl_fan_filter(nmax, mmax, r_lon, r_lat)


def green(theta, n_max_love, cha="vd", cha_norm=False):
    """Compute the GREEN functions for surface vertical/horizontal displacement, geoid, gravity change, tilt change, free air gravity, potential, and stress.

    Args:
        theta (array) : theta[n], Radial distance away the observed station
        n_max_love (int) : The maximum order of Love number to be used to calculate Green function. If n_max_love>13000, then there is an interpolation be used to get the value of that Love number. Note: the Love number values are in a file named love.prem.13000.dat, which should be located in c:\wfl\cvf\lib directory, to use WFL_GREEN, you must copy this file in your working dirctory.
        cha (string) : Character to index which Green function can be returned to green_value(n). cha='vd' or 'VD' , cha='hd' or 'HD', cha='gc' or 'GC', cha='dg' or 'DG', cha='ig' or 'IG', cha='dig' or 'DIG', cha='tc' or 'TC', cha='fag' or 'FAG', cha='pc' or 'PC', cha='sc' or 'SC' for Green function of surface Vertical Displacement (VD), surface Horizontal Displacement (HD), Geoid Change (GC), Direct effect of Gravity (DG), Indirect effect of Gravity (IG), Direct + Indirect effect of Gravity (DIG), Tilt Change (TC), Free Air Gravity (FAG), Potential Change (PC), and Stress Change (SC). For more details, please see Comments|Reference.
        cha_norm (bool): If True, then the output Green functions will be normalized by different factors. See Algorithm for more details.

    Return:
        array : green_value[n]

    Notes:
        When cha_norm is True, then the green function is normalized by below factor:
            *1011*Earth_radius*theta (theta is in radian unit) for 'vd', 'hd' and 'gc';
            *1017*Earth_radius*theta (theta is in radian unit) for 'dg', 'ig', 'dig', and 'fag';
            *1011*(Earth_radius*theta)2 (theta is in radian unit) for 'tc' and 'sc';
            *108*Earth_radius*theta (theta is in radian unit) for 'pc'.
        All the unit is in kg, m, s, etc.
        For 'vd' and 'hd', the equation used outside the summation is Earth_radius/Earth_mass
        For 'dg','ig',and 'dig', the equation used outside the summation is g/Earth_mass.
        Other green functions are not used for me, so do it yourself.
        This subroutine is modified from a Fortran77 program provided by Dr. Wang Hansheng.


"""
    cha_norm = 'norm' if cha_norm else 'unnorm'
    conds = ['vd', 'hd', 'gc', 'dg', 'ig', 'dig', 'tc', 'fag', 'pc', 'sc']
    if cha.lower() not in conds:
        raise ValueError('cha should be in [{}]'.format(' '.join(conds)))

    return _geophysics.py_wfl_green(theta, n_max_love, cha, cha_norm=cha_norm)


def green_func_theta(n):
    """Generate angle distance (0.1d-4 ~ 179.5) with smaller step to larger step. The output will be used in Green functions calculation.

    Arg:
        n (int) : # of output (n should be 452, if you want to get all the values)

    Return:
        array : theta[n] angle distance from 0.00001 to 179.5 (unit: degree)

    Notes:
        if input n<=452, the theta[n] will be output; if input n>452, theta[453:]==0.0

"""
    return _geophysics.py_wfl_green_func_theta(n)


def grace_destripe(c, s, n, m, segment=3, fit_order=3):
    pass


def potential_temperture_to_situ(s, theta, pressure, reference_pressure=0):
    """convert potential temperature (theta) to in situ temperature (t) in ocean.

    Args:
        s (float) : salinity (unit: psu)
        theta (float) : potential temperature (unit: deg Celsius) at pressure level pr0
        pressure (float) : in situ pressure level (unit: Pascal)
        reference_pressure (float) : reference pressure level for potential temperature (unit: Pascal). If not presented, pr0==0, which means the sea surface.

    Return:
        float : t, in situ temperature (unit: deg Celsius) at pressure level reference_pressure.

    Notes:
        Potential temperature has been defined classically as the temperature an element of seawater would have if raised adiabatically with no change of salinity to atmospheric pressure (or sea surface with reference_pressure==0). More generally, the potential temperature can be defined as the temperature resulting from an adiabatic displacement to a reference pressure pr that may be greater or less than the initial pressure p. The potential temperature :math:`\\theta` can be computed from the adiabatic lapse rate :math:`\\Gamma`,

        .. math::
            \\theta(s, t, p, p_r) = t + \\int_{p}^{p_r} \\Gamma(s, t^{'}, p)dp^{'}

        by integration along an adiabat, i.e. :math:`\\partial t^{'}/ \\partial p^{'} = \\Gamma^{'}`. The potential temperature can be evaluated by using an empirical formula or by numerical integration of the defining equation.

        Here, the potential temperature :math:`\\theta(s, t, p, p_r)` at reference pressure pr can be computed with sufficient precision using a 4th order Runge-Kutta integration algorithm (Fofonoff, 1977).

        The in situ temperature can be calculated by the inversion of the above process of potential temperature calculation.

        For the atmosphere, the relation between potential temperature and in situ temperature is very simple, and can be acquired by,

        .. math::
            \\theta = T(p_r / p)^{0.286}

        where :math:`p_r` is reference pressure level (usually 1000mb), and p is in situ pressure level, and :math:`T` is in situ temperature. Note, here, the unit of :math:`\theta` and :math:`T` is Kelvin, not deg Celsius.

    References:
        .. [1] Bryden, 1973, polynomial for adiabatic lapse rate ! and runge-kutta 4-th order intergration algorithm.
        .. [2] Bryden H,1973, deep-sea res., 20, 401-408
        .. [3] Rofonoff, N., 1977, Deep-sea res., 24, 489-491

"""
    return _geophysics.py_wfl_in_situ_temp(s, theta, pressure, pr0=reference_pressure)


def situ_temperatrue_to_potential(s, t, pressure, reference_pressure=0):
    """
    Args:
        s (float) : salinity (unit: psu)
        t (float) : in situ temperature (unit: deg Celsius) at pressure level reference_pressure.
        pressure (float) : in situ pressure level (unit: Pascal)
        reference_pressure (float) : reference pressure level for potential temperature (unit: Pascal). If not presented, pr0==0, which means the sea surface.

    Return:
        float : potential temperature (unit: deg Celsius) at pressure level pr0

    See Also:
        :func:`potential_temperture_to_situ`
"""
    return _geophysics.py_wfl_potential_temp(s, t, pressure, pr0=reference_pressure)


def j2j4j6(r, gm, omiga, inv_f):
    """Calculate J2, J4 and J6 from Earth referenced ellipsoid parameters

    Args:
        r (float) : Earth's equatorial radius (in meter), if input r=0.d0, then substitutes it with value of 6378137m.
        gm (float) :Geocentric gravitational constant (in m**3/s**2), if input gm=0.0d0, then substitutes it with value of 3.986005e+14 m**3/s**2.
        omiga (float) :Earth's angular velocity (in Rad/s), if input omig=0.d0, then substitutes it with value of 7.292115D-5 rad/s.
        inv_f (float) :inverse of Earth's flattening (==1/f), if input inv_f=0.d0, then substitutes it with value of 298.257222101.

    Returns:
        tuple:
            * j (array) : j(1:3) corresponding with J2, J4 and J6, respectively.
            * c(array) : c(1:3) corresponding with fully normalized zonal spherical harmonic coefficients C20, C40, and C60, respectively.

    Notes:
        The default input Earth's ellipsoid parameters are adopted from IERS 2003 conventions, see "Table 1.2: Parameters of the Geodetic Reference System GRS80" for details.

        The output values are often used to calculate gravity anomalies from spherical harmonic coefficients geoid modle such as GGM03C, GGM03S, EGM96, and so on.

    reference:
       .. [1] Rapp R. H., A Fortran program for the computation of gravimetric quantities from high degree spherical harmonic expansions, OSU Reports No. 334, 1982.
"""
    return _geophysics.py_wfl_j2j4j6(r, gm, omiga, inv_f)


def lod2torque(dt, lod):
    """Change LOD (length of day)  to torque (z axis)

    Args:
        dt (float) : time interval step of input array
        lod (array) : lod[n], LOD change (unit: ms)

    Returns:
        array: torque_lod[n], The corresponding torque in z axis, unit: Hadley, see :func:`torque` for details

    Notes:
        torque_lod(1),torque_lod(n) are interpolated from torque_lod(2) and torque_lod(n-1), respectively.
"""
    return _geophysics.py_wfl_lod2torque(dt, lod)


def mountain_torque(sp, hm):
    """Using torque method to calculate mountain torques from atmosphere/water to solid Earth

    Args:
        sp : (lon,lat), surface pressure (unit: Pa) or water height (unit: cm*100 or 0.1mm) (centered grid data) (north->south direction (N->S), for 1*1 grid, the coordinate is like (0.5,89.5),(1.5,89.5),(2.5,89.5),.....(359.5,89.5),(0.5,88.5),(1.5,88.5)....)
        hm : (lon,lat), topography/orography data (unit: m) (centered grid data) (N->S, for 1*1 grid, the coordinate is like (0.5, 89.5), (1.5, 89.5),...(359.5, 89.5), (0.5, 88.5), (1.5, 88.5)...)

    Returns:
        tuple:
            * te : (lon,lat), complex output (x+iy), ellipsoidal torque for IB hypothesis (centered grid data, N->S) (unit:Hadley, 1HU=10**18kg.m.m/s/s=N.m)
            * ten : (lon,lat), complex output (x+iy), ellipsoidal torque for NIB hypothesis (centered grid data, N->S) (unit:Hadley)
            * tm_dp : (lon,lat),  complex output (x+iy), mountain torque for IB hypothesis using pressure gradient method (centered grid data, N->S) (unit:Hadley)
            * tmn_dp : (lon,lat), complex output (x+iy), mountain torque for NIB hypothesis using pressure gradient method (centered grid data, N->S) (unit:Hadley)
            * tm_dh : (lon,lat), complex output (x+iy), mountain torque for IB hypothesis using topography gradient method (centered grid data, N->S) (unit:Hadley)
            * tmn_dh : (lon,lat), complex output (x+iy), mountain torque for NIB hypothesis using topography gradient method (centered grid data, N->S) (unit:Hadley)
            * tm3_dp : (lon,lat), z component, mountain torque for IB hypothesis using pressure gradient method (centered grid data, N->S) (unit:Hadley)
            * tm3n_dp : (lon,lat), z component, mountain torque for NIB hypothesis using pressure gradient method (centered grid data, N->S) (unit:Hadley)
            * tm3_dh : (lon,lat), z component, mountain torque for NIB hypothesis using pressure gradient method (centered grid data, N->S) (unit:Hadley)
            * tm3n_dh : (lon,lat), z component, mountain torque for NIB hypothesis using topography gradient method (centered grid data, N->S) (unit:Hadley)

    Notes:
        This is the subset of subroutine :func:`torque`, but with centered data input only for surface pressure/water data.
        For water height data such as land-water, usually use non-IB results. For surface pressure data, usually use IB results.
        For the water data input, if you have water data with cm unit, then you should multiply the water data by 100 as the input SP parameter.

        The mountain torque equations in the three Egger B. J. paper are wrong. The mountain torque equation of Schindelegger M. paper is correct

    References:
        .. [1] Egger B. J. and K. P. Hoinka, Equatorial components of global atmospheric angular momentum: covariance functions, Q. J. R. Meteorol. Soc., Vol. 128, p1137-1157, 2002.
        .. [2] Egger B. J. and K. P. Hoinka, Mountain torques and the equatorial components of global angular momentum, J. Atmos. Sci., p2319-2331, 2000.
        .. [3] Egger B. J., K. Weickmann, and K. P. Hoinka, Angular momentum in the global atmospheric circulation, Rev. Geophys., Vol. 45, RG4007, 2007.
        .. [4] Schindelegger M., J. Bohm and D. Salstein, Seasonal and intraseasonal polar motion variability as deduced from atmospheric torques, J. Geod. Geoinfor., Vol. 1, P89-95, 2012
        .. [5] Schindelegger M., D. Salsteim and J. Bohm, Recent estimates of Earth-atmosphere interaction torques and their use in studying polar motion variability, J. Geophys. Res., Vol. 118, p4586-4598, 2013.
        .. [6] Feldstein S. B., The dynamics of atmospherically driven intraseasonal polar motion, J. Atmos. Sci., Vol. 65, P2290-2307,2008.

    """
    return _geophysics.py_wfl_mountain_torque(sp, hm)


def prem(n, parameter='density'):
    """Acquiring PREM parameters (depth, density, mass, gravity acceleration, pressure, and moment of inertia)

    Args:
        n (int) : The layers of PREM, n must equal 76 here. (n==76)
        parameter (string) : default density
            * 'mass', then  prem_out(:,2) output integral of mass of Earth from Earth's center to surface (unit: kg)
            * 'gravity', then prem_out(:,2) output gravity acceleration of Earth in each layer  (unit:m/s**2)
            * 'pressure', then prem_out(:,2) output pressures of Earth in each layer (unit: Pa)
            * 'inertial', then  prem_out(:,2) output integral of moment of inertia of Earth from Earth's center to surface (unit: kg.m**2)
            * 'density', then  prem_out(:,2) output density of Earth in each layer (unit: kg/m**3)

    Returns:
        2darray : prem_out[n, 2], prem_out[:,0] is the depth of Earth (unit: km), prem_out[:,1] is the other output parameters defined by "prem_parameter" optional input parameter.

    Notes:
        PREM: The isotropic model parameters in the table below are taken from:
            Dziewonski A.M., Anderson D.L.(1981): Preliminary reference Earth model. Phys. Earth planet. Inter., 25, 297.
        Some values, specified in the table II of the paper, which can likely be accurately fitted by cubic splines are left out. Between the given points within each complex block, the parameters are interpolated by means of natural cubic splines. If only one value of a material parameter is specified per a layer, the parameter is assumed constant within the layer. In the input data for complete ray tracing, the values of density are specified only at the interfaces because elsewhere they have no influence upon the high-frequency asymptotics.

"""


def rm_sp_from_sg(fre, fre_idx, spec_sp, spec_sg, m_order=0):
    """Removing surface pressure (solar heating tide) form observed super-conductive gravity

    Args:
        fre (array) : frequency index of power spectrum
        fre_idx (array) : fre_idx[2], frequency index for determination of line equation, then using this line equation to do linear interpolation to get all frequency correction coefficients. Using this coefficients to remove solar heating tide from SG power spectrum.
        spec_sp (array) : Power spectrum of surface pressure at the station.
        spec_sg (array) : Power spectrum of the super-conductive gravity at the station.
        m_order (int) : If larger than 0, then use it as the order of polynomial to fit the back ground noise spectrum. If equal to 0, then use BIC criterion to get the best fit order of polynomial fit, and using this best fit order to fit the back ground noise spectrum.

    Returns:
        array : SG power spectrum removed solar heating tide.
"""
    return _geophysics.py_wfl_rm_sp_from_sg(fre, fre_idx, spec_sp, spec_sg, m_order=m_order)


def search_resonance(fre, spec, n_window, q):
    """Search resonance response from Super-conduction Gravity (SG) power spectrum

    Args:
        fre (array) : Frequency index for the input power spectrum
        spec (array) : Power spectrum to be analysised
        n_window (int) : The parameter "n_window" is a very important parameter when computing the resonance parameter "resonance_parameter". You can change it to different value to see the different results. Typically, n_window should small (maybe around 11, both 5 points from the resonance point) for a true resonance.
        q (float) : Quality factor, typical value of 100~400.

    Return:
        array : Resonance parameter (big value means there is resonance)

    References:
        .. [1] Smylie D. E., J. Hinderer, B. Richter, and B. Ducarme, The porduct spectra of gravity and barometric pressure in Europe, Physics of the Earth and Planetary Interiors, 80, 135-157, 1993. eq. (50-53)

"""
    return _geophysics.py_wfl_search_resonance(fre, spec, n_window, q)


def sg_multistation(alon, alt, sg):
    """Get translational modes (Slichter modes) from multi-station Superconductive Gravimeter (SG) residuals

    Args:
        alon (array) : alon[n], Longitude coordinates of the corresponding SG stations
        alat (array) : alat[n], Latitude coordinates of the corresponding SG stations. Note, alon(i), and alat(i) should be the ith station's coordinates.
        sg (2darray): sg[n, ns], SG residuals time series (n length, n >= 3) at ns stations.

    Returns:
        array : Signals of 3 translational modes from input SG residuals time series (length n, stations ns). sg_out[:,0] for prograde equatorial components, sg_out[:,1] for axial components, and sg_out[:,2] for retrograde equatorial components.

    References:
        .. [1] Courtier, N., et al., Global superconducting gravimeter observations and the search for the translational modes of inner core, Physics of the Earth and Planetary Interiors (PEPI), 117, P3-20, 2000. Eq. (5), (7)

"""
    return _geophysics.py_wfl_sg_multistation(alon, alt, sg)


def soiltemp(amp, period, phase, average_t, depth, year_begin, nday, td=1, unit='day'):
    """Get soil temperature at given depth from surface temperature

    Args:
        amp (float): Amplitude of surface temperature for given period (amp unit: degree)
        period (float) : period of surface temperature, period>2day. (unit: day or hour, see optional parameter unit for details)
        phase (float) : initial phase of the surface temperature (unit: degree). The amp and phase can be accquired from surface temperature time series by using subroutine :func:`harmonic`. The phase definition is the same as in :func:`harmonic`.
        average_t (float) : Multiyear average temperature of surface temperature (unit: degree)
        depth (float) : the depth of soil, depth>0 (unit: mm)
        year_begin (int) : the begin year for output, the time index of first output is Jan 1, year_begin.
        nday (int) : unit of input optional parameter unit (days or hours) from the firs time index
        TD (float) : Thermal Diffusivity. If not presented, the default value is TD=1mm**2/s   (unit: mm**2/s)
        unit (string) : Unit of output time index, default is "day". If not presented or presented and unit='day' (case-insensitive), then the unit of output time index interval is day. If presented and unit='hour' (case-insensitive), then output time index interval is hour. For annual, semi-annual period, and period>2days, you can use unit 'day'. For diurnal, semi-diurnal and period>2 hours, you can use unit 'hour'. The output year and MJD index is in unit year and day with suitable decimals and intervals in day or hour.

    Returns:
        array : soil_t(nday,3) soil_t[:,0] time index (year); soil_t[:,1] time index (MJD); soil_t[:,2] soil temperature at given depth (unit: degree)

    References:
        .. [1] Yan H., W. Chen, Y. Zhu, G. Liu, and M. Zhong, thermal expansion effect on GPS height change, 2009.
"""
    return _geophysics.py_wfl_soiltemp(amp, period, phase, average_t, depth, year_begin, nday, td=td, unit=unit)


def solar_heat_tide(fre, spec, index_fre, n_window, q):
    """Remove solar heating tides effects from Super-conduction Gravity (SG) power spectrum

    Args:
        fre (array) : fre[n], Frequency index for the input power spectrum
        spec(array) : spec[n], Power spectrum to be analysised
        index_fre (array) : index_fre[m], Integer number position of the above solar heating tides frequency in fre(n), e.g. index_fre(1)=83, means fre(83)=1/24 cycle/hr. Note, index_fre(1)>3*n_window/2 , index_fre(m) < n-3*n_window/2 . Also note, index_fre(:) must be monotonic.
        n_window (int) : Window width (must be odd) for removing the solar heating tides effect. e.g. n_window=25, means there are 25/2 points on both side of fre( index_fre(1) ) will be removed, respectively. Note, n_window < n
        Q (float) : Quality factor, typical value of 100~400.

    Returns:
        array : spec_out[n], Power spectrum which the solar heat tides had been removed from the input spec.

    Notes:
        The effects of removing solar heat tides is limited. This is NOT the best method to remove solar heat tides.

    Reference:
        ..[1] Smylie D. E., J. Hinderer, B. Richter, and B. Ducarme, The porduct spectra of gravity and barometric pressure in Europe, Physics of the Earth and Planetary Interiors, 80, 135-157, 1993. eq. (47-49)
"""
    return _geophysics.py_wfl_solar_heat_tide(fre, spec, index_fre, n_window, q)


def stress_torque(ets, nts):
    """Using torque method to calculate stress (both friction and gravity wave) torques from atmosphere to solid Earth

    Args:
        ETS (array) : ETS[lon, lat] : eastward turbulent (or gravity wave) surface stress  (unit: N/m2) (lattice grid data) (N->S) . If unit is (N m-2 s), then divide the data with the time step of the data. For example, if daily data were used, then divide 3600s*24h=86400s to get unit of N/m2 . Other input parameters is the same as here.
        NTS (array) : NTS[lon, lat] : northward turbulent (or gravity wave) surface stress (unit: N/m2) (lattice grid data) (N->S)

    Returns:
        tuple :
            * TF[lon, lat]: complex output (x+iy), friction (or gravity wave) torque (centered grid data, N->S) (unit:Hadley)
            * TF3[lon, lat]: z component, friction (or gravity wave) torque (centered grid data, N->S) (unit:Hadley)

    Notes:
        This is the subset of subroutine :func:`torque`, but with centered data input only for stress data.

    References:
        .. [1] Egger B. J. and K. P. Hoinka, Equatorial components of global atmospheric angular momentum: covariance functions, Q. J. R. Meteorol. Soc., Vol. 128, p1137-1157, 2002.
        .. [2] Egger B. J. and K. P. Hoinka, Mountain torques and the equatorial components of global angular momentum, J. Atmos. Sci., p2319-2331, 2000.
        .. [3] Egger B. J., K. Weickmann, and K. P. Hoinka, Angular momentum in the global atmospheric circulation, Rev. Geophys., Vol. 45, RG4007, 2007.
        .. [4] Schindelegger M., J. Bohm and D. Salstein, Seasonal and intraseasonal polar motion variability as deduced from atmospheric torques, J. Geod. Geoinfor., Vol. 1, P89-95, 2012
        .. [5] Schindelegger M., D. Salsteim and J. Bohm, Recent estimates of Earth-atmosphere interaction torques and their use in studying polar motion variability, J. Geophys. Res., Vol. 118, p4586-4598, 2013.
        .. [6] Feldstein S. B., The dynamics of atmospherically driven intraseasonal polar motion, J. Atmos. Sci., Vol. 65, P2290-2307,2008.

    See Also:
        :func:`mountain_torque`, :func:`torque`
"""
    return _geophysics.py_wfl_stress_torque(ets, nts)


def thermal_vertical_deformation(amp, period, phase, year_begin, nday, pr=0.25, ltec=12e-6, td=1, unit='day'):
    """Using a half-space head conduction model to calculate the vertical change of bed rock

    Args:
        amp (float) : Amplitude of surface temperature for given period (amp unit: degree)
        period (float) : period of surface temperature, period>2day. (unit: day or hour)
        phase (float) : initial phase of the surface temperature (unit: degree). The amp and phase can be accquired from surface temperature time series by usingsubroutine HARMONIC. The phase definition is the same as in HARMONIC.
        year_begin (int) :the begin year for output, the time index of first output is Jan 1, year_begin.
        nday (int) : unit of input optional parameter unit (days or hours) from the firs time index
        pr (float) : Poisson's Ratio. If not presented, the default value is PR=0.25
        ltec (float) : Linear Thermal Expansion Coefficient. If not presented, the default value is LTEC=12.D-6  ( or 12*10-6).  (unit: 1/degree)
        td (float) : Thermal Diffusivity. If not presented, the default value is TD=1mm**2/s   (unit: mm**2/s)
        unit (string) : Unit of output time index, default is "day". If not presented or presented and unit='day' (case-insensitive), then the unit of output time index interval is day. If presented and unit='hour' (case-insensitive), then output time index interval is hour. For annual, semi-annual period, and period>2days, you can use unit 'day'. For diurnal, semi-diurnal and period>2 hours, you can use unit 'hour'. The output year and MJD index is in unit year and day with suitable decimals and intervals in day or hour.

    Returns:
        array : dh(nday,3), dh[:,0] time index (year); dh[:,1] time index (MJD); dh[:,2] vertical change of bed rock induced by thermal expansion (unit: mm)

    References:
        .. [1] Yan H., W. Chen, Y. Zhu, G. Liu, and M. Zhong, thermal expansion effect on GPS height change, submitted to .., 2009.
        .. [2] Dong D., P. Fang, Y. Bock, M. Cheng and S. Miyazaki, 2002, Anatomy of apparent seasonal variations from GPS-derived site position time series, J Geophys Res, Vol. 107, No. B4, 2075, doi:10.1029/2001JB000573

"""
    conds = ['day', 'hour']
    if unit not in conds:
        raise ValueError('unit should be in [{}]'.format(" ".join(conds)))
    return _geophysics.py_wfl_thermalv(amp, period, phase, year_begin, nday, pr=pr, ltec=ltec, td=td, unit=td)


def torque(sp, hm, ets, nts, egws, ngws):
    """Using torque method to calculate torques from atmosphere to solid Earth

    Args:
        sp (array) : sp[lat, lon+1], surface pressure (unit: Pa), (lattice grid data) (north->south direction (N->S), for 1*1 grid, the coordinate is like (0.,90.),(1.,90.),(2.,90.),.....(359.,90),(0.,89.),(1.,89.)....)
        hm (array) : hm[lon, lat], topography/orography data (unit: m) (centered grid data) (N->S, for 1*1 grid, the coordinate is like (0.5, 89.5), (1.5, 89.5),...(359.5, 89.5), (0.5, 88.5), (1.5, 88.5)...)
        ets (array) :[lon, lat+1], eastward turbulent surface stress (unit: n/m2) (lattice grid data) (n->s) . if unit is (n m-2 s), then divide the data with the time step of the data. for example, if daily data were used, then divide 3600s*24h=86400s to get unit of n/m2 . other input parameters is the same as here.
        nts (array) :[lon, lat+1], northward turbulent surface stress (unit: n/m2) (lattice grid data) (n->s)
        egws (array) :[lon, lat+1], eastward gravity wave surface stress (unit: n/m2) (lattice grid data) (n->s)
        ngws (array) :[lon, lat+1], northward gravity wave surface stress (unit: n/m2) (lattice grid data) (n->s)

    Returns:
        tuple :
            * te (array) : [lon,lat], complex output (x+iy), ellipsoidal torque for ib hypothesis (centered grid data, n->s) (unit:hadley, 1hu=10**18kg.m.m/s/s=n.m)
            * ten (array) : [lon,lat], complex output (x+iy), ellipsoidal torque for nib hypothesis (centered grid data, n->s) (unit:hadley)
            * tm_dp (array) : [lon,lat], complex output (x+iy), mountain torque for ib hypothesis using pressure gradient method (centered grid data, n->s) (unit:hadley)
            * tmn_dp (array) : [lon,lat], complex output (x+iy), mountain torque for nib hypothesis using pressure gradient method (centered grid data, n->s) (unit:hadley)
            * tm_dh (array) : [lon,lat], complex output (x+iy), mountain torque for ib hypothesis using topography gradient method (centered grid data, n->s) (unit:hadley)
            * tmn_dh (array) : [lon,lat], complex output (x+iy), mountain torque for nib hypothesis using topography gradient method (centered grid data, n->s) (unit:hadley)
            * tf (array) : [lon,lat], complex output (x+iy), friction torque (centered grid data, n->s) (unit:hadley)
            * tg (array) : [lon,lat], complex output (x+iy), gravity wave torque (centered grid data, n->s) (unit:hadley)
            * tm3_dp (array) : [lon,lat], z component, mountain torque for ib hypothesis using pressure gradient method (centered grid data, n->s) (unit:hadley)
            * tm3n_dp (array) : [lon,lat], z component, mountain torque for nib hypothesis using pressure gradient method (centered grid data, n->s) (unit:hadley)
            * tm3_dh (array) : [lon,lat], z component, mountain torque for ib hypothesis using topography gradient method (centered grid data, n->s) (unit:hadley)
            * tm3n_dh (array) : [lon,lat], z component, mountain torque for nib hypothesis using topography gradient method (centered grid data, n->s) (unit:hadley)
            * tf3 (array) : [lon,lat], z component, friction torque (centered grid data, n->s) (unit:hadley)
            * tg3 (array) : [lon,lat], z component, gravity wave torque (centered grid data, n->s) (unit:hadley)

    References:
         .. [1] Egger B. J. and K. P. Hoinka, Equatorial components of global atmospheric angular momentum: covariance functions, Q. J. R. Meteorol. Soc., Vol. 128, p1137-1157, 2002.
        .. [2] Egger B. J. and K. P. Hoinka, Mountain torques and the equatorial components of global angular momentum, J. Atmos. Sci., p2319-2331, 2000.
        .. [3] Egger B. J., K. Weickmann, and K. P. Hoinka, Angular momentum in the global atmospheric circulation, Rev. Geophys., Vol. 45, RG4007, 2007.
        .. [4] Schindelegger M., J. Bohm and D. Salstein, Seasonal and intraseasonal polar motion variability as deduced from atmospheric torques, J. Geod. Geoinfor., Vol. 1, P89-95, 2012
        .. [5] Schindelegger M., D. Salsteim and J. Bohm, Recent estimates of Earth-atmosphere interaction torques and their use in studying polar motion variability, J. Geophys. Res., Vol. 118, p4586-4598, 2013.
        .. [6] Feldstein S. B., The dynamics of atmospherically driven intraseasonal polar motion, J. Atmos. Sci., Vol. 65, P2290-2307,2008.

    See Also:
        :func:`mountain_torque`, :func:`stress_torque`, :func:`eaam_xyz_j2j3`

"""
    return _geophysics.py_wfl_torque(sp, hm, ets, nts, egws, ngws)


def torque2lod(dt, torque):
    """Change torque (z axis) to dLOD (length of day)  ( lod(i+1)-lod(i) )

    Args:
        dt (float) : Time interval (step) of input array, unit: second
        torque (array) : Torque in z axis, unit: Hadley, see WFL_TORQUE for details

    Returns:
        array: dlod_torque, The corresponding LOD change from input torque, unit: millisecond, ms. Note, this output is == lod(i+1)-lod(i), NOT  lod(i). To get lod(i), please use subroutine :func:`torque2lod_adams`.

    Notes:
        dlod_torque(1:2) are integrated directly, while dlod_torque(3:n) are integrated by Simpson integral method.
"""
    return _geophysics.py_wfl_torque2lod(dt, torque)


def torque2lod_adams(dt, torque, lod0, t=None):
    """Using 4th order Adams-Moulton method to get LOD (length of day) change from torque z

    args:
        dt (float) : time interval (step) of input array, unit: second
        torque (array) : torque[n], torque in z axis, unit: hadley, see wfl_torque for details
        lod (array): lod[3], the first 3 initial values for ordinary differential equation d(lod0)/dt=torque, unit: ms
        t (array) :  t[n], time index of input torque

    returns:
        tuple:
            * lod(n) : the corresponding lod change from input torque, unit: millisecond, ms.
            * lod_poly(n) :polynomial term of output lod(n). due to numerical integral will cause a polynomial term in the output lod(n) results, we use polynomial fit to remove this polynomial term if you present t and lod_poly at the same time in subroutine wfl_torque2lod_adams.


"""
    if t is not None:
        return _geophysics.py_wfl_torque2lod_adams(dt, torque, lod0, t=t, flag=True)
    else:
        return _geophysics.py_wfl_torque2lod_adams(dt, torque, lod0, flag=False)


def vd_station_index(nlat, site_lon, site_lat):
    """Get station index for calculation of vertical displacement of special station by using spherical harmonic method. The output of this subroutine will be used by :func:`load_vert_def` or :func:`cs2vertical`.

    Args:
        nlat (int) :  or degree n, is the center grid # in latitude or the maximum degree n of spherical harmonic coefficients
        site_lon (float) : coordinate longitude (0~360deg) of the station. Unit: degree
        site_lat (float) : coordinate latitude(-90~90deg) of the station. Unit: degree

    Returns:
        tuple
            * depart_lon, alat[nlat] : see details in subroutine :func:`cs2vertical`. Unit: degree
            * idx_lon, idx_lat : Index of the input station location in the output of :func:load_vert_def` or :func:`cs2vertical`. For example, the output of :func:load_vert_def` or :func:`cs2vertical` --> u(idx_lon,idx_lat) is the input station vertical displacement

    Notes:
        If you want to get the vertical deformation of a special GPS station located at (site_lon, site_lat) then do follow:

        1.use global atmosphere/hydrology/ocean pressure data (equivalent water height) and :func:`fun_exp_fft` to get spherical harmonic coefficients c and s of the equivalent water heigh;
        2. use :func:`vd_station_index` to get the output parameters of depart_lon, alat, idx_lon, and idx_lat;
        3. use :func:`cs2vertical` with step 2 output parameters to get vertical deformation u, then u(idx_lon,idx_lat) is the vertical deformation of the above GPS station.

    reference:
        .. [1] Mangiarotti, S., A. Cazenave, Annual vertical crustal motions predicted from surface mass redistribution and observed by space geodesy, J Geophys Res, Vol. 106, B3, P4277-4291, 2001.

    See Also:
        :func:`cs2vertical`, :func:`load_vert_def`, :func:`func_sum`, :func:`func_sum_fft`, :func:`func_exp`, :func:`func_exp_fft`
"""
    return _geophysics.py_wfl_vd_station_index(nlat, site_lon, site_lat)


def wind2tao(windu, windv, air_density=1.3):
    """

    Args:
        windu(float): Wind velocity in East direction (unit: m/s)
        windv(float): Wind velocity in North direction (unit: m/s)
        air_density (float): air density. Default air_density=1.3km/m3

    Returns:
        tuple
            * stressu (float): Wind stress in East direction (unit: kg/(m*s**2) , or N*m-2)
            * stressv (float): Wind stress in North direction (unit: kg/(m*s**2) , or N*m-2)

    Notes:
        This function is used for 10m wind in the ocean area.

    Reference:
        .. [1] Trenberth, et al. (1990): J. Phys. Oceanogr., 20, 1742-1760. Eq. (2)

    See Also:
        :func:`windv2tao`
"""
    return _geophysics.py_wfl_wind2tao(windu, windv, air_density=air_density)


def windv2tao(v):
    """Get wind stress from wind velocity

    Arg:
        w (float): Wind velocity (unit: m/s)

    Return:
        float : tao, Wind stress (unit: kg/(m*s**2) , or N*m-2)

    Notes:
        Wind stress tao has the same direction as wind vector W, the input in the above function is absolute value (:math:`|W|`) of wind vector W, without consideration of wind direction.

        east wind (u) and north wind (v) should be first combine to get W at each grid, after that change :math:`|W|` to tao, and then use the direction of W to get tao(u), tao(v).

    Reference:
        .. [1] Gill A. E., Atmosphere-Ocean Dynamics, Academic press, London, 1982, p29, Eq. (2.4.2)

    Examples:
        >>> windv2tao(4.3) # 0.02542375
        >>> windv2tao(10) # 0.155
    See Also:
        :func:`wfl_wind2tao`
"""
    return _geophysics.py_windv2tao(v)


def x2prograd(a1, fai1, a2, fai2):
    """Change excitation function in Earth's rotation domain to prograde and retrograde format

    Args:
        a1 (float) :The amplitude of excitation function x1
        fai1 (float) :The phase of excitation function x1 (unit: degree)
        a2 (float) :The amplitude of excitation function x2
        fai2 (float) :The phase of excitation function x2 (unit: degree)

    Returns:
        tuple
            *a_pro (float): The amplitude for prograde component
            *fai_pro (float): The phase for prograde component
            *a_retro (float): The amplitude for retrograde component
            *fai_retro (float): The phase for retrograde component

    Notes:
        Equation:  change :math:`x_1 + ix_2 = a_1 \\cos(w_t - \\phi_1) + i a_2 \\cos(w_t - \\phi_2)` to :math:`a_{pro} \\exp(i \\phi_{pro}) \\exp(i w_t) + a_{retro} \\exp(i \\phi_{retro}) \\exp(-i w_t)`, here :math:`x_1` and :math:`x_2` are excitation function used in Earth's rotation domain, and :math:`w` is angular frequency omega. If the model of :math:`x_1` and :math:`x_2` are not like above (maybe sine format, see :func:`harmonic` for transfer to the above format in details )

    See Also:
        :func:`prograd2x0`

"""
    return _geophysics.py_x2prograd(a1, fai1, a2, fai2)


def prograd2x0(a1, fai1, a2, fai2):
    """Change excitation function in Earth's rotation domain to prograde and retrograde format

    Args:
        a1 (float) : The amplitude of excitation function x1
        fai1 (float) : The phase of excitation function x1 (unit: degree)
        a2 (float) : The amplitude of excitation function x2
        fai2 (float) : The phase of excitation function x2 (unit: degree)

    Returns:
        tuple:
            * f1 (float) :The real part for prograde component
            * f2 (float) :The image part for prograde component
            * g1 (float) :The real part for retrograde component
            * g2 (float) :The image part for retrograde component

    See Also:
        :func:`x2prograd`

"""
    return _geophysics.py_x2prograd_0(a1, fai1, a2, fai2)


def angle_dis(fai1, lambda1, fai2, lambda2, ch='radian'):
    """Calculate geocentric angle (or corresponding Earth's surface distance) between two points

    Args:
        fai1, lambda1 (float) : First point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360)
        fai2, lambda2 (float) : Second point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360)
        ch (string) : Must be character. If ch='angle', then return the geocentric angle (unit: radian) to angle_dis_result. If ch='distance', then return the corresponding Earth's surface distance (unit: meter) to angle_dis_result.

    Return:
        angle_dis_result (float), Geocentric angle (unit: radian) or corresponding Earth's surface distance (unit: meter).

    Notes:
        The formula used here (see below for detail) is the same as used in :func:`sphere_len`, but with the different formal.

        .. math::
            \\cos(\\phi)=\\cos(\\theta_{1})\\cos(\\theta_{2})+sin(\\theta_{1})sin(\\theta_{2})\\cos(\\lambda_{1}-\\lambda_{2})

        :math:`\\phi` is geocentric angle between two points, and :math:`\\theta_{1}`,:math:`\\theta_{2}` are colatitude of :math:`\\phi_{1}` and :math:`\\phi_{2}` respectively. The Earth's radius adopted here is 6371012m.

    Examples:
        >>> angle_dis(0., 0., 0., 2.) # 0.0349065850399
        >>> angle_dis(0., 0., 0., 2., ch='distance') # 222390.272168
"""
    return _geophysics.py_angle_dis(fai1, lambda1, fai2, lambda2, ch=ch)


def eaam_xyz_slp(slp, hm, is_sp=False):
    """

    Args:
        slp (array) : slp[lon, lat], Sea level pressure/surface pressure. Sea level pressure if is_sp is not presented. Surface pressure if is_sp is presented, and is_sp==.true.   (North->South) (unit: mbar. Note: 1mbar==1hPa)
        hm (array) : hm[lon, lat], Topography data or elevation at each points (North->South) (unit: m)

    Returns:
        tuple:
            * fp (array) : [lon, lat], Mass term for equatorial component with IB-hypothesis.(unit:10-7)
            * fpn (array) : [lon, lat], Mass term for equatorial component with non IB-hypothesis. (unit:10-7)
            * f3p (array) : [lon, lat], Mass term for axial component with IB-hypothesis.(unit:10-7)
            * f3pn (array) : [lon, lat], Mass term for axial component with non IB-hypothesis.(unit:10-7)
            * pmn (array) : [lon, lat], Total mass of atmosphere. (non-IB) (unit:1020 g, not kg)
            * pm (array) : [lon, lat], Total mass of atmosphere. (IB) (unit:1020 g, not kg)
            * xypn (array) : [lon, lat], Mass center of x and y coordinate (non-IB) (unit: cm)
            * xyp (array) : [lon, lat], Mass center of x and y coordinate (IB) (unit: cm)
            * zpn (array) : [lon, lat], Mass center of z coordinate( non-IB) (unit: cm)
            * zp (array) : [lon, lat], Mass center of z coordinate( IB) (unit: cm)
            * j2n (array) : [lon, lat], J2 value (non-IB) (unit: 10-10)
            * j2 (array) : [lon, lat], J2 value (IB) (unit: 10-10)
            * j3n (array) : [lon, lat], J3 value (non-IB) (unit: 10-10)
            * j3 (array) : [lon, lat], J3 value (IB) (unit: 10-10)

    Notes:
        Note: The better choose is to use Surface Pressure(SP) as the input but not Sea Level Pressure(SLP). Due to here we only use approximate equation (see below) to calculate SP from SLP, this will cause minor error at high altitude area. Especially, the total mass of atmosphere is not correct when use SLP as input, while it's correct when use SP as input. Other outputs are also have minor errors when use SLP as input.

        All input data must be sorted from North to South and at the grid center, order of data may like below (longitude, latitude): (1.25,88.75), (3.75,88.75), ...,(358.75,88.75), (1.25,86.25),(3.75,86.25),.....

        If the data is at grid lattice not at the grid center, then you can use INTERPOLATE2D and INTERPOLATE3D to interpolate the original field first, then  use the interpolation field in EAAM_XYZ_SLP.

        When calculate surface pressure from sea level pressure, we use topography height and an approximate equation like below:
        sp=slp*(1-0.0065*h/288.15)**(9.81/(287.04*0.0065))
        where sp,slp,h means surface pressure, sea level pressure and topography height. If you want to get more accuracy output, please choose EAAM_XYZ_J2J3, but in this case you must provide geopotential height data as input.

        Note: all the parameters about Earth is Eubank's value like below:
        CC   =7.1236 D44       Earth principal moment of inertia
        CCM=7.1236 D44
        CA   =0.003664*CC     C - A

        Barnes's Earth parameter (not used here) :
        CC=7.04D44
        CCM=7.1236 D44
        CA=0.00333*CC

    Reference:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.
"""
    return _geophysics.py_eaam_xyz_slp(slp, hm, is_sp=False)


def eaam_xyz_wind(ltrop, pl, u, v, slp=None, hm=None, z=None, en=False):
    """Get Effective Atmospheric Angular Momentum (EAAM), only for wind-term

    Args:
        ltrop (int) : Lays great and equal (>=) 10mbar(1000->10mbar, include 10mbar. Note: 1mbar==1hPa)
        pl (array) : [lay], Pressure values in each vertical lay (such as: 1000., 850., 700.,...,30.10.,...) (unit: mbar)
        u (array) : [lon,lat,lay], East-ward wind velocities at the each levels (North->South)  (unit: m/s)
        v (array) : [lon,lat,lay], North-ward wind velocities at the each levels (North->South) (unit: m/s)
        slp (array) : [lon,lat], Sea level pressure (North->South) (unit: mbar). The slp,hm,and z must be presented at the same time. If they are not presented, then calculate wind term of effective angular momentum using "BP" method; otherwise using "SP" method to calculate effective angular momentum.
        hm(array) : [lon,lat], Topography data or elevation at each points (North->South) (unit: m)
        z[lon,lat,lay]:Geo-potential height at each level from 1000mb to top (North->South)  (unit: m)
        en (bool) : default False, call :func:`eaam_xyz_wind`, True for :func:`eaam_xyz_wind_en`

    Returns:
        tuple
            * fw (array) : [lon,lat], Wind term for equatorial component (total 1000mbar->top).(unit:10-7)
            * fwl (array) : [lon,lat], Wind term for equatorial component but in upper troposphere alone(10mbar->top).(unit:10-7)
            * f3w (array) : [lon,lat], Wind term for axial component (total 1000mbar->top).(unit:10-7)
            * f3wl (array) : [lon,lat], Wind term for axial component but in upper troposphere alone(10mbar->top).(unit:10-7)

    Notes:
        All input data must be sorted from North to South and at the grid center, order of data may like below (longitude, latitude): (1.25,88.75), (3.75,88.75), ...,(358.75,88.75), (1.25,86.25),(3.75,86.25),.....

        If the data is at grid lattice not at the grid center, then you can use INTERPOLATE2D and INTERPOLATE3D to interpolate the original field first, then  use the interpolation field in EAAM_XYZ_WIND.

        Note: all the parameters about Earth is Eubank's value like below:
        CC   =7.1236 D44       Earth principal moment of inertia
        CCM=7.1236 D44
        CA   =0.003664*CC     C - A

        Barnes's Earth parameter (not used here) :
        CC=7.04D44
        CCM=7.1236 D44
        CA=0.00333*CC

        When calculate wind term or relative angular momentum, use "SP" method which means put the first level (1000mbar) wind as the surface wind in the mountain; use BP method, which do not consider the land surface and calculate all the wind AAM integral from 1000mbar to top.

    Reference:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.


    See Also:
        :func:`eaam_xyz_slp`, :func:`eaam_xyz_j2j3`, :func:`ewam_xyz_j2j3`, :func:`eaam_xyz_wind_en`


"""
    func = _geophysics.py_eaam_xyz_wind_en if en else _geophysics.py_eaam_xyz_wind
    if None in [slp, hm, z]:
        return func(ltrop, pl, u, v, flag=False)
    else:
        return func(ltrop, pl, u, v, flag=True, slp=slp, hm=hm, z=z)


def eam_barnes2eubanks(eam, cha):
    """Change Barnes's EAAM, EWAM or EOAM to Eubanks's EAAM, EWAM or EOAM

    Args:
        eam (float) : Effective angular Momentum (unit: free)
        cha (string): cha eaqual one of follow character, 'xp', 'xc', 'yp', 'yc' , 'zp', 'zc'

    Returns:
        float: eam_results, Effective Angular Momentum (unit: same as input)

    Notes:
        cha denote X,Y or Z coordinate for pressure (P) and current (C) term
        Change Barnes's EAAM or EWAM or EOAM to Eubanks's EAAM or EWAM or EOAM

    See Also:
        :func:`eam2masms`, :func:`masms2eam`, :func:`am2eam`
"""
    conds = ['xp', 'xc', 'yp', 'yc', 'zp', 'zc']
    if cha not in conds:
        raise ValueError("cha should be in [{}]".format(conds))
    return _geophysics.py_eam_barnes2eubanks(eam, cha)


def eam2masms(xin, cha):
    """ Change effective angular momentum(unit: 1.0d-7 rad) to unit mas (for x,y) and ms (for z)

    Args:
        xin (float) : Effective Angular Momentum (unit: 1.0d-7 rad )
        cha (string) : cha is one of follow character, 'x',  'y' , 'z'

    Return:
        float : Polar Motion value (unit: mas) or Length of Day value (unit: ms)

    Notes:
        cha denote X,Y or Z coordinate
        This function only used in the Earth's Rotation domain, which change Effective Angular Momentum to polar motion (unit: mas) or Length of Day (unit: ms).
        mas=milli-arcsencod, ms=millisecond

    Reference:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.


"""
    conds = ['x', 'y', 'z']
    if cha not in conds:
        raise ValueError("cha should be in [{}]".format(conds))
    return _geophysics.py_eam2masms(xin, cha)


def ewam_xyz_j2j3(pre, evp, runoff, hm, grid_status, method='ocean_zero'):
    """Get Effective Water Angular Momentum (EWAM), geocenter, J2 and J3

    Args:
        pre (array) : [lon,lat+1], Global precipitation (North->South) (unit: mm).  Note, pre, evp and runoff must be with size of (lon,lat+1), if grid_status=='center', you can only use (lon,lat) to store the data; if grid_status=='lattice', then you will use (lon,lat+1) to store the data, and in the inner of subroutine, pre, evp and runoff will be interpolated into grid center values.  When output, the value is centered at (lon,lat) and with unit of cm.
        evp (array) : [lon,lat+1], Global evaporation (North->South) (unit: mm). see note above.  When output, the value is centered at (lon,lat) and with unit of cm.
        runoff (array) [lon,lat+1], Land runoff (North->South) (unit: mm). see note above. When output, the value is centered at (lon,lat) and with unit of cm.
        hm (array) : [lon,lat] :Topography data or elevation at each points (unit: m)
        grid_status (string) : status of the grid format, grid_status is case-insensitive. If grid_status=='center', then the input pre, evp and runoff will be taken as grid center data with size of (lon,lat); if grid_status=='lattice', then the input pre, evp and runoff will be taken as grid lattice data with size of (lon,lat+1) and will be interpolated to the grid center values in the inner of this subroutine .
        method (string) : If method=='ocean_zero', then the total water change in the ocean are set to 0.0. If others, then the total water change in ocean are defined as  pre(i,j)-runoff_ocean-evp(i,j), where runoff_ocean is  sum of -runoff(i,j) in land divide ocean area, i.e., the mean runoff in the ocean.

    Returns:
        tuple:
            * fp (array) : [lon,lat], Mass term for equatorial component with IB-hypothesis.(unit:10e-7). Here IB means total water in the ocean will be averaged as a single layer in the ocean.
            * fpn (array) : [lon,lat], Mass term for equatorial component with non IB-hypothesis. (unit:10e-7)
            * f3p (array) : [lon,lat], Mass term for axial component with IB-hypothesis.(unit:10e-7)
            * f3pn (array) : [lon,lat], Mass term for axial component with non IB-hypothesis.(unit:10e-7)
            * pmn (array) : [lon,lat], Total mass of atmosphere. (non-IB) (unit:10e20 g, not kg)
            * pm (array) : [lon,lat], Total mass of atmosphere. (IB) (unit:10e20 g, not kg)
            * xypn (array) : [lon,lat], Mass center of x and y coordinate (non-IB) (unit: cm)
            * xyp (array) : [lon,lat], Mass center of x and y coordinate (IB) (unit: cm)
            * zpn (array) : [lon,lat], Mass center of z coordinate( non-IB) (unit: cm)
            * zp (array) : [lon,lat], Mass center of z coordinate( IB) (unit: cm)
            * j2n (array) : [lon,lat], J2 value (non-IB) (unit: 10e-10)
            * j2 (array) : [lon,lat], J2 value (IB) (unit: 10e-10)
            * j3n (array) : [lon,lat], J3 value (non-IB) (unit: 10e-10)
            * j3 (array) : [lon,lat], J3 value (IB) (unit: 10e-10)

    Notes:
        All input data must be sorted from North to South, order of data may like below (longitude, latitude): (0,90), (1,90), (2,90),...,(359,90), (0,89),(1,89),.....

        Note: all the parameters about Earth is Eubank's value like below:
        CC   =7.1236 D44       Earth principal moment of inertia
        CCM=7.1236 D44
        CA   =0.003664*CC     C - A

        Barnes's Earth parameter (not used here) :
        CC=7.04D44
        CCM=7.1236 D44
        CA=0.00333*CC

    Reference:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.
"""
    return _geophysics.py_ewam_xyz_j2j3(pre, evp, runoff, hm, grid_status, method=method)


def gauss_grid():
    """ Get Gauss grid coordinate 192*94 in longitude*latitude

    Returns:
        tuple
            * gauss_lon (192): Gauss grid longitude coordinate
            * gauss_lat (94): Gauss grid latitude coordinate
            * gauss_lat_edges (95): Gauss grid latitude edges coordinate
"""
    return _geophysics.py_gauss_grid()


def geocenter(zlat_s, zlat_n, zlon_w, zlon_e, rho, r1, r2):
    """Calculate geocenter variation caused by mass change at a grid

    Args:
        zlat_s (float) : Latitude of grid's southern lattice (unit: degree -90~90)
        zlat_n (float) : Latitude of grid's northern lattice (unit: degree -90~90)
        zlon_w (float) : Longitude of grid's western lattice (unit: degree 0~360)
        zlon_e (float) : Longitude of grid's eastern lattice (unit: degree 0~360)
        rho (float) : Density (unit: kg/m**3)
        r1 (float) : Radius of calculated volume (lower one or  benchmark) (unit: m),  see comment for details
        r2 (float) : Radius of calculated volume (higher one or benchmark plus changes) (unit: m),  see comment for details

    Returns:
        tuple
            * x (float) : Geocenter displacement in x direction (unit: m)
            * y (float) : Geocenter displacement in y direction (unit: m)
            * z (float) : Geocenter displacement in z direction (unit: m)

    Notes:
        If r2>r1, it is easy to understand, if r2<r1, it means a geocenter change in opposite direction. Usually, we calculate geocenter from global mass changes, such as from Sea Surface Height (SSH) changes, at this case, we only know the SSH variations or :math:`\\Delta_{ssh}`. To calculate geocenter variations caused by :math:`\\Delta_{ssh}` which may be positive or negative, we can define r1 as the Earth's average radius Re (r1=Re), and :math:`r2=Re+\\Delta_{ssh}`. Here, r2 may greater than r1, or r2 may smaller then r1.
"""
    return _geophysics.py_geocenter(zlat_s, zlat_n, zlon_w, zlon_e, rho, r1, r2)


def masms2eam(xin, cha):
    """Change Polar Motion (unit: mas) or Length of Day (unit: ms) to effective angular momentum (unit: 1.0d-7 rad)
    Args:
        xin (float) : Polar Motion (unit: mas) or Length of Day (unit: ms)
        cha (string) : cha is one of follow character, 'x',  'y' , 'z'

    Return:
        float : Effective Angular Momentum (unit: 1.0d-7 rad )

    Notes:
        cha denote X,Y or Z coordinate
        This function only used in Earth Rotation domain, which change polar motion or Length of Day to EAM.
        mas=milli-arcsencod, ms=millisecond

    Reference:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.


"""
    conds = ['x', 'y', 'z']
    if cha not in conds:
        raise ValueError("cha should be in [{}]".format(conds))
    return _geophysics.py_masms2eam(xin, cha)

def oam_int(rho,r1,r2,zlat_s,zlat_n,zlon_w,zlon_e,u,v):
    """Calculate grid ocean angular momentum use integral calculation

    Args:
        rho (float) : Ocean water's density at grid (unit: kg/m**3)
        r1 (float) : The distance from Earth center to the grid lower edge in vertical (or z) direction (unit: m)
        r2 (float) : The distance from Earth center to the grid higher edge in vertical (or z) direction (unit: m) (r1<=r2)
        zlat_s (float) : Latitude of grid's southern lattice (unit: degree -90~90) (zlat_s<zlat_n, south->north)
        zlat_n (float) : Latitude of grid's northern lattice (unit: degree -90~90) (zlat_s<zlat_n, south->north)
        zlon_w (float) : Longitude of grid's western lattice (unit: degree 0~360) (zlon_w<zlon_e, west->east)
        zlon_e (float) : Longitude of grid's eastern lattice (unit: degree 0~360) (zlon_w<zlon_e, west->east)
        u (float) : Ocean current velocity in latitudinal direction (East is positive) (unit: m/s)
        v (float) : Ocean current velocity in meridional direction (North is positive) (unit: m/s)

    Returns:
        tuple
            * xoamp, yoamp, zoamp (float) Planetary angular momentum (or mass term) in x, y and z direction (unit: kg.m**2/s)
            * xoamc, yoamc, zoamc (float) Relative angular momentum (or current term) in x,y and z direction (unit: kg.m**2/s)

    Notes:
        :func:`oam_int` is more accuracy than simple method :func:`oam_simple`. If you only provide rho but no u and v, then you can only get mass term; if you only provide u and v but no rho, then you can only get current term; if all provide, you can get all output. If you do not provide rho, u, or v, then set them to 0.0d0.


"""
    return _geophysics.py_oam_int(rho,r1,r2,zlat_s,zlat_n,zlon_w,zlon_e,u,v)

def oam_simple(rho,r,zlat,zlon,dlat,dlon,dr,u,v):
    """ Calculate grid ocean angular momentum use simple calculation

    Args:
        rho (float) : Ocean water's density at grid (unit: kg/m**3)
        r1 (float) : The distance from Earth center to the grid lower edge in vertical (or z) direction (unit: m)
        r2 (float) : The distance from Earth center to the grid higher edge in vertical (or z) direction (unit: m) (r1<=r2)
        zlat (float) : Latitude coordinate of grid center (unit: degree -90~90)
        zlon (float) : Longitude coordinate of grid center (unit: degree 0~360)
        dlon (float) : Longitude length (or interval) of grid (unit: degree 0~360)
        dr (float) : The z (or vertical) distance of grid (unit: m)
        u (float) : Ocean current velocity in latitudinal direction (East is positive) (unit: m/s)
        v (float) : Ocean current velocity in meridional direction (North is positive) (unit: m/s)

    Returns:
        tuple
            * xoamp, yoamp, zoamp (float) Planetary angular momentum (or mass term) in x, y and z direction (unit: kg.m**2/s)
            * xoamc, yoamc, zoamc (float) Relative angular momentum (or current term) in x,y and z direction (unit: kg.m**2/s)

    Notes:
        :func:`oam_int` is more accuracy than simple method :func:`oam_simple`. If you only provide rho but no u and v, then you can only get mass term; if you only provide u and v but no rho, then you can only get current term; if all provide, you can get all output. If you do not provide rho, u, or v, then set them to 0.0d0.

"""
    return _geophysics.py_oam_simple(rho,r,zlat,zlon,dlat,dlon,dr,u,v)

def ocean_densiy(s, c, depth, lat):
    """Calculate ocean water density from ocean depth, salinity and temperature field

    Args:
        s (float) :Salinity (unit: psu)
        c (float) :Temperature (unit: degC)
        depth (float) :Ocean depth (unit: m)
        lat (float) : Latitude (unit: deg, -90~90). Gravity variation with latitude is considered here.

    Return:
        float : Ocean water's density (unit: kg/m**3)
"""
    return _geophysics.py_ocnden(s,c,depth,lat=lat)


def pm2xy(dt, pm, cha='mas'):
    """Change Polar Motion data PMx and PMy to excitation fields chi_x,chi_y

    Args:
        dt (float) :  Time interval of vector (or array) pm (unit: days)
        pm (array) :  pm[n], Polar motion data (complex, real part is x component, image part is y component) (unit: mas==arcsecond/1000. ) (the direction of x is Greenwich and y is 90 west)
        cha (string) : If cha=='mas' or cha=='MAS', then chi_xy has unit of mas; if cha=='e-7' or cha=='E-7', then the unit of chi_xy is e-7.


    Return:
        array : chi_xy (n), Excitation fields (complex), real and image part is x and y component respectively (unit: mas or e-7) (the direction of chi_x is Greenwich and chi_y is 90 east)

"""
    conds = ['mas', 'e-7']
    if cha not in conds:
        raise ValueError("cha should be in [{}]".format(conds))
    return _geophysics.py_pm2xy(dt, pm, cha)

def prograd2x_0(f1,f2,g1,g2):
    """Change prograde and retrograde format to excitation function in Earth's rotation domain

    Args:
        f1 (float) :The real part for prograde component
        f2 (float) :The image part for prograde component
        g1 (float) :The real part for retrograde component
        g2 (float) :The image part for retrograde component


    Returns:
        tuple:
            * a1 (float) : The amplitude of excitation function x1
            * fai1 (float) : The phase of excitation function x1 (unit: degree)
            * a2 (float) : The amplitude of excitation function x2
            * fai2 (float) : The phase of excitation function x2 (unit: degree)

    See Also:
        :func:`x2prograd`, :func:`prograd2x0`

"""
    return _geophysics.py_prograd2x_0(f1,f2,g1,g2)


def soda_xyz_st():
    """Get x, y and z coordinate of SODA ocean model for salinity and temperature field

    Returns:
        tuple:
            * x360 (360): x coordinate (or longitude) (unit: degree 0~360)
            * y128 (128): y coordinate (or latitude) (unit: degree -60~60)
            * yedges (129): y coordinate edges (unit: degree -60~60)
            * z20 (20): z coordinate (or depth) (unit: m)
            * zedges (21): z coordinate edges (unit: m)
"""
    return _geophysics.py_soda_xyz_st()


def soda_xyz_uvh():
    """Get x, y and z coordinate of SODA ocean model for oceanic current u,v and SSH field

    Returns:
        tuple:
            * uvhx (360) :x coordinate (or longitude) (unit: degree 0~360)
            * uvhy (128) :y coordinate (or latitude) (unit: degree -60~60)
            * uvhyedges (129) :y coordinate edges (unit: degree -60~60)
"""
    return _geophysics.py_soda_xyz_uvh()


def sphere_len(b1, l1, b2, l2):
    """Calculate spherical length (or distance) between point a1(b1,l1) and a2(b2,l2) on Earth's surface

    Args:
        b1, l1 (float) : First point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360)
        b2, l2 (float) : Second point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360)

    Return:
        float : slen, Spherical length (or distance) between two points on Earth's surface (unit: km )

    Notes:
        The formula used here (see below for detail) is the same as used in ANGLE_DIS, but the different formal.

    .. math::
        [\\sin(\\phi/2)]^{2}=\\sin((\\phi_{1}-\\phi_{2})/2)^{2}+\\cos(\\phi_{1})\\cos(\\phi_{2}) [\\sin((\\lambda_{1}-\\lambda_{2})/2)]^{2}

    :math:`\\phi` is geocentric angle between two points, and (:math:`\\phi_{1}, \\lambda_{1}`) and (:math:`\\phi_{2}, \\lambda_{2}`) are (latitude,longitude) of two points respectively.  The Earth's radius adopted here is 6371.012km.

"""
    return _geophysics.py_sphere_len(b1,l1,b2,l2)

def sphere_len_inv(b, l, slen_bn):
    """Give a point and spherical length depart form the point, calculate increase (or decrease) in latitude and longitude

    Args:
        b, l (float): First point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360)
        slen (float): Spherical length (or distance) on Earth's surface (unit: km )

    Returns:
        tuple : delta_b, delta_l, Increase (or decrease) form first point's coordinate (latitude, longitude) (unit: degree, -90~90, 0~360).

    Notes:
        Because inverse problem is not unique, so we here only give two special case: 1) delta_b=0, and delta_l come from :func:`sphere_len_inv`; 2) delta_l=0 and delta_b come from :func:`sphere_len_inv`. So in subroutine :func:`sphere_len_inv`, we give delta_b and delta_l at the same time, but you must use the two results separately, see example for details.

    Reference:
        YAN Haoming, IB response and it's effect on Earth's rotation

    Examples:
        >>> slen = sphere_len(63, 10, 63, 14)
        >>> print slen
        >>> 201.893584896
        >>> res = sphere_len_inv(63, 10, slen)
        >>> print res
        >>> (1.815669210053464, 3.9999999999999996)
        >>> print sphere_len(63, 10, 63, 10 + res[1]) # correct use of inverse result
        >>> 201.893584896
        >>> print sphere_len(63, 10, 63-res[0], 10) # correct use of inverse result
        >>> 201.893584896
        >>> print sphere_len(63, 10, 63+res[0], 10+res[1]) # misuse !
        >>> 281.02153828
"""
    return _geophysics.py_sphere_len_inv(b, l, slen_bn)


def volume(dlon, lat1, lat2, z1, z2):
    """Calculate grid volume of Earth at different oceanic depth or above Earth's surface.

    Args:
        dlon (float): Longitude interval (or length) of the grid (unit: degree)
        lat1 (float): Co-latitude of northern grid border (0 is north pole, 180 is south pole, unit: degree)
        lat2 (float): Co-latitude of southern grid border (same as lat1, lat2>lat1)
        z1 (float): Ocean depth in grid above (positive, unit: m)
        z2 (float): Ocean depth in grid bottom (positive, z2>z1, unit: m)

    Return:
        float : volume_result, Volume of grid (unit: m**3)

    Notes:
        If calculate mountain volume, then z1 and z2 are negative.
"""
    return _geophysics.py_volume(dlon, lat1, lat2, z1, z2)


def atg(s, t, p):
    """Get adiabatic temperature gradient or adiabatic lapse rate.

    Args:
        s (float) :salinity (unit: psu)
        t (float) :temperature (unit: deg Celsius)
        p0 (float) :pressure level (unit: Pascal)

    Return:
        float : atg,adiabatic temperature gradient (unit: deg Celsius/decibars).

    Notes:
        The adiabatic lapse rate :math:`\\Gamma(s,t,p)` (deg Celsius/decibar) is defined as the change of temperature per unit pressure for an adiabatic change of pressure of an element of seawater. It's assumed that no heat or salt is exchanged with the surroundings so that the pressure change is both adiabatic and isentropic. From thermodynamic considerations, the adiabatic lapse rate :math:`\\Gamma`, a function of pressure, temperature and salinity can be expressed as

        .. math::
            \\Gamma(s,t,p) = T \\left( \\partial V / \\partial t \\right) / C_{p}

        where :math:`T=t+273.15` is absolute temperature (Kelvin), ( :math:`\\partial V / \\partial t` ) (m3/(kg deg Celsius) ) is thermal expansion and :math:`C_{p}` (J/(kg deg Celsius)) specific heat of seawater at constant pressure.

        The lapse rate :math:`\\Gamma` is positive except at low salinities, temperatures and pressures where :math:`\\partial V / \\partial t` is negative. Typical values in the oceanic range are 1~2*104 deg Celsius/decibar.

        Adiabatic lapse rate can be calculated from the equation of state and specific heat (Bryden, 1973) or from direct meansurements.

    References:
        .. [1] BRYDEN,H.,1973,DEEP-SEA RES.,20,401-408

    Examples:
        >>> atg(40, 40, 1e8) # 0.00032559758

"""
    return _geophysics.py_wfl_atg(s, t, p)
