from _geophysics import geophysics as _geophysics


def am2eam(am, cha='xp'):
    """ Change Angular Momentum (AM) to Effective Angular Momentum (EAM)

    Args:
        am (float) : Angular Momentum (unit: :math:`kg.m^2/s`)
        cha (string) : cha eaqual one of follow character, 'xp', 'xc', 'yp', 'yc', 'zp', 'zc'

    Returns:
        float : Effective Angular Momentum (unit: 1.0e-7 rad)

    Notes:
        cha denote x,y or z coordinate for pressure (p) and current (c) term
        All the parameters about Earth is Eubank's value like below:
        >> cc  = 7.1236 #earth principal moment of inertia
        >> ccm = 7.1236
        >> ca  = 0.003664*cc       #c - a

        barnes's earth parameter (not used here) :
        >> cc=7.04
        >> ccm=7.1236
        >> ca=0.00333*cc

        this function only used in earth rotation domain, which change oceanic am to oceanic eam.

    References:
        .. [1] Barnes R. T. H., Atmospheric angular momentum fluctuations, length of day changes and polar motion, Proc. R. Soc. Lond. A, 387,31-73, 1983.

    See Also:
        todo

"""
    cond = ['xp', 'xc', 'yp', 'yc', 'zp', 'zc']
    if cha.lower() not in cond:
        raise ValueError("should should be in [%s]".format(' '.join(cond)))
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
        .. math:: \cos(\phi) = \cos( \theta1_) * \cos(\theta_2) + \sin(\theta_1) * \sin(\theta_2) * \cos(\lambda_1 - \lambda_2) 

        :math:`\phi` is geocentric angle between two points, and :math:`theta_1`,:math:`\theta_2` are colatitude of :math:`phi_1` and :math:`phi_2` respectively.
        The Earth's radius adopted here is 6371012m.

    Examples:
        todo

    See Also:
        todo
"""

    cond = ['angle', 'distance']
    if ch.lower() not in cond:
        raise ValueError("ch should be in [%s]".format(" ".join(cond)))
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
        .. math::  ds=r*r*\sin(\theta)d(\theta)d(\lambda)

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
        raise ValueError("should should be in [%s]".format(' '.join(cond)))
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
        raise ValueError("should should be in [%s]".format(' '.join(cond)))
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
            degvar(n)=\sum(c(n,m)^2+s(n,m)^2) \qquad m=0,1,2,....,n \quad n=0,1,2....N (N the maximum degree)

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
        p (array) : p(nlon,nlat), atmospheric surface pressure field (sorted by latitude (90~-90) first, then sorted by longitude (0~360). e.g. the coordinate of SIGMA should be (0.5,89.5) (1.5, 89.5),....(359.5,89.5) (0.5, 88.5) (1.5, 88.5).......... (359.5, -89.5) for 1x1 P). Unit of  P is Pascal. (Note: in the ocean, P must be set to 0.0d0)

    Returns:
        tuple
            * u(nlon, nlat): Radial elastic defomation of the Earth's surface calculated by using experimental method(unit: mm) (grid mean value) , see comment for details. (Note: only 80N~60S have the values, other areas are 0)
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
            dh = \sum\sum\left\{ P_{lm}\left(\cos \theta \right) \left[ C_{lm}\cos(m \phi) + S_{lm}\sin(m \phi) \right ] \right \}

        which is different with surface density defined by Wahr, 1998. The difference is the factor earth semi-radius and water density.

    References:
        .. [1] Wahr J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.


"""
    cond = ['g2h', 'h2g']
    if cha.lower() not in cond:
        raise ValueError("should should be in [%s]".format(' '.join(cond)))
    return _geophysics.py_wfl_cs2cs(c, s, cha=cha.upper())


def cs2vertical(c, s, water_density=1e3, frame='ce'):
    return _geophysics.py_wfl_cs2vertical(c, s, water_density=1e3, frame=frame)


def var_reduction(x, x_correction):
    return _geophysics.py_wfl_var_reduction(x, x_correction)
    # pass

import numpy as np
import matplotlib.pyplot as plt
cs = np.loadtxt('cs.txt')
# c = cs[:, 2]
c = cs[:, 2].reshape(181, 181).T
s = cs[:, 3].reshape(181, 181).T
print _geophysics.py_wfl_cs2vertical.__doc__
v = cs2vertical(c, s)
print v[50, :]