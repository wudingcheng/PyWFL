# encoding: utf-8
from __future__ import absolute_import
from . import _coordinate

__all__ = ["gd2lambert", "blh2xyz", "xyz2blh", "lambert2gd", "xyz2enu", "gc2gd", "geo2ups", "datum"]


def gd2lambert(lon, lat, cm, ps, pn, ol, ea=None, ef=None):
    """Change geodetic coordinate (longitude and latitude) to two parallel latitudes Lambert projection coordinate (X and Y)

    Args:
        alon,alat (float) :Coordinate of longitude (-180,180) and latitude (-90,90) in degree unit
        cm (float) :Central meridian of Lambert projection (-180,180) in degree unit
        ps (float) : Southern standard parallel (-90,90) in degree unit
        pn (float) : Northern standard parallel (-90,90) in degree unit
        ol (float) :latitude of origin for Lambert projection (-90,90) in degree unit
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            x, y (float) : x and y coordinate of Lambert projection in meter unit

    References:
        .. [1] R. T. WILLIAMS, LAMBERT AND MERCATOR MAP PROJECTIONS IN GEOLOGY AND GEOPHYSICS, COMPUTER AND GEOSCIENCES,  VOL.21,NO.3,P353-364,1995.
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_2lam(lon, lat, cm, ps, pn, ol, ea=ea, ef=ef)


def blh2xyz(b, l, h, ea=None, ef=None):
    """Change WGS84 BLH (latitude, longitude, and height) to X,Y,Z coordinate

    Args:
        b (float) : wgs84 coordinate, latitude, unit: degree
        l (float) : wgs84 coordinate, longitude, unit: degree
        h (float) : wgs84 coordinate, height, unit: meter
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            x, y, z (float) : corresponding x, y, z coordinate in unit meter
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    if not -90 <= b <= 90:
        raise ValueError('Latitude should be between -90 and 90!')
    if not 0 <= l <= 360:
        raise ValueError('Longitude should be between 0 and 360!')
    return _coordinate.py_wfl_blh2xyz(b, l, h, ea=ea, ef=ef)


def xyz2blh(x, y, z, ea=None, ef=None):
    """Transform Cartesian(Geocentric) x,y,z to geodetic coordinates (latitude, longitude, height) based on the exact solution (Borkowski,1989)

    Args:
        x (float): unit meter
        y (float): unit meter
        z (float): unit meter
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            b, l, h (float) : wgs84 coordinate b, l, h are latitude, longitude and height. unit: degree, degree, and meter

    References:
        .. [1] Borkowski, K. M. (1989).  "Accurate algorithms to transform geocentric to geodetic coordinates", *Bulletin Geodesique*, v. 63, pp. 50-56.
        .. [2] Borkowski, K. M. (1987).  "Transformation of geocentric to geodetic coordinates without approximations", Astrophysics and Space Science, v. 139, n. 1, pp. 1-4.
        .. [3] Borkowski, K. M.   Correction in (1988), v. 146, n. 1, p. 201.
        .. [4] Vermeille H., An analytical method to transform geocentric into geodetic coordinates, J Geod (2011) 85:105–117,DOI 10.1007/s00190-010-0419-x
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_map_xyz2blh(x, y, z, ea=ea, ef=ef)


def lambert2gd(x, y, cm, ps, pn, ol, ea=None, ef=None):
    """Change  two parallel latitudes Lambert projection coordinate (X and Y) to geodetic coordinate (longitude and latitude)

    Args:
        x (float) : x coordinate of Lambert projection in meter unit
        y (float) : y coordinate of Lambert projection in meter unit
        cm (float) :Central meridian of Lambert projection (-180,180) in degree unit
        ps (float) : Southern standard parallel (-90,90) in degree unit
        pn (float) : Northern standard parallel (-90,90) in degree unit
        ol (float) :latitude of origin for Lambert projection (-90,90) in degree unit
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            alon,alat (float) :Coordinate of longitude (-180,180) and latitude (-90,90) in degree unit

    References:
        .. [1] R. T. WILLIAMS, LAMBERT AND MERCATOR MAP PROJECTIONS IN GEOLOGY AND GEOPHYSICS, COMPUTER AND GEOSCIENCES,  VOL.21,NO.3,P353-364,1995.
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_invlam(x, y, cm, ps, pn, ol, ea=ea, ef=ef)


def xyz2enu(x, y, z, reference_lon, reference_lat, ea=None, ef=None):
    """Change ECEF(Earth-Center,Earth-Fixed) clestian X,Y,Z coordinate to local E(east),N(north),U(up) coordinate

    Args:
        x (array): unit meter
        y (array): unit meter
        z (array): unit meter
        reference_lon, reference_lat (float) : reference longitude (-180,180) and latitude (-90,90) in degree unit
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        e, n, u (array) : local E(east), N(north), U(up) coordinate (meter)

"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_xyz2enu(x, y, z, ea=ea, ef=ef, s_lat=reference_lat, s_lon=reference_lon)


def gc2gd(x, y, z, ea=None, ef=None):
    """Transform geocentric coordinates to geodetic for a reference ellipsoid of specified form.

    Args:
        x (float): unit meter
        y (float): unit meter
        z (float): unit meter
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            b, l, h (float) : wgs84 coordinate b, l, h are latitude, longitude and height. unit: degree, degree, and meter

    References:
        .. [1] Fukushima, T., "Transformation from Cartesian to geodetic coordinates accelerated by Halley's method", J.Geodesy (2006),79: 689-693
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_map_gc2gd(x, y, z, ea=ea, ef=ef)


def geo2ups(lat, lon, ea=None, ef=None):
    """Change geodetic coordinate (lon, lat ) to UPS (X, Y)
    UPS: Universal Polar Stereographic projection --Polar stereographic (variant A) projection (EPSG 9810)

    Args:
        lat (float) : input geodetic latitude in degree (-90~90)
        lon (float) : input longitude in degree (-180~180)
        ea (float) : Ellipsoid semimajor axis (in meter unit)
        ef (float) : Ellipsoid flattening. ea and ef must be presented at the same time. If no ea and ef parameter are presented, then use WGS 84 value.

    Returns:
        tuple:
            x,y (float): output UPS x(east) and y(north) in meter

    Notes:
        For UPS (x,y) plot, we have

        For North pole case: the original longitude (0) is at S direction, longitude increase anti-colockwise.

        For South pole case: ................................ N ............................. colockwise.

        figure show ::

             north pole case                south pole case
                -/+180                              0
            -90        90                      -90      90
                   0                              -/+180

    References:
        .. [1] Surveying and Positioning Guidance Note Number 7, part 2
            Coordinate Conversions and Transformations including Formulas, p52-56
            Revised - November 2009, available at http://www.epsg.org
"""
    if None in [ea, ef]:
        ea, ef = datum(1)
    return _coordinate.py_wfl_map_geo2ups(lat, lon, ea=ea, ef=ef)


def datum(reference_name):
    """Ellipsoids parameters copied from GEOMAP3.6 MANUAL.

    Args:
        reference_name (int) : the reference_name index list as follows:

    Returns:
        tuple:
            * a (float) : equatorial radius (meters)
            * f (float) : flattening (a number around 0.00335, i.e. around 1/298.)

    Notes:
        1. WGS84
        2. GRS80
        3. WGS72  (below sorted)
        4. Airy1930
        5. Airy1930(modified for Ireland1965)
        6. Australian
        7. Bessel1841
        8. Bessel1841(modified for NGO1948)
        9. Bessel1841(modified for Schwarzeck)
        10.    Clarke1858
        11.    Clarke1866
        12.    Clarke1866(modified for Michigan)
        13.    Clarke1880
        14.    Clarke1880(modified for Arc1950)
        15.    Clarke1880(modified for IGN)
        16.    Clarke1880(modified for Jamaica)
        17.    Clarke1880(modified for Merchich)
        18.    Clarke1880(modified for Palestine)
        19.    Everest1830
        20.    Everest1830(modified for Kalianpur)
        21.    Everest1830(modified for Kertau)
        22.    Everest1830(modified for Timbalai)
        23.    Fischer1960
        24.    Fischer1960(modified for SouthAsia)
        25.    Fischer1968
        26.    GRS67
        27.    Hayford 1909
        28.    Helmert1906
        29.    Hough
        30.    IAG 75 (1980 xi'an coordiante)
        31.    Indonesian
        32.    International1924
        33.    Krassovsky 1948 (1954 beijing coordinate)
        34.    MERIT83
        35.    NewInternational1967
        36.    NWL10D
        37.    NWL9D
        38.    OSU86F
        39.    OSU91A
        40.    Plessis1817
        41.    SouthAmerican
        42.    Sphere (Note: f==0.0 for sphere)
        43.    Struve1860
        44.    Walbeck
        45.    WarOffice
        46.    WGS60
        47.    WGS66
"""
    return _coordinate.py_wfl_map_datum(reference_name)
