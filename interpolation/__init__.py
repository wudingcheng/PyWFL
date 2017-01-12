import _interpolation
import numpy as np


def bilinear(x, y, val, x1, y1):
    """Do bi-linear interpolation on a single grid

    Args:
        x (array) : x(2), coordinate index in x direction, from left to right.
        y (array) : y(2), coordinate index in y direction, from bottom to top.
        val (2d-array) : val(0:4), 4 points values from 1 to 4, counter-clockwise staring from the lower left
        x1 (float) : coordinate index in x direction, :math:`min(x1) \leq x1 \leq  \max(x1)`
        y1 (float) : coordinate index in y direction, :math:`min(y1) \leq y1 \leq  \may(x1)`

    Returns:
        interpolate point (x1, y1) value.

    Note:
        x1, x2 must in the domain of x and y respectively.

    References:
        .. [1] Bilinear interpolation, https://en.wikipedia.org/wiki/Bilinear_interpolation

    Examples:
        >>> bilinear([1, 2], [1, 2], [1, 2, 3, 4], 1.2, 1.2)
"""
    return _interpolation.py_bilinear(x, y, val, x1, y1)


def bilinear_square():
    return _interpolation.py_bilinear_square()


def bicubic(x, y, val, x1, y1, method='ordinary'):
    """Do bi-cubic interpolation on rank 2 dimensions field

    Args:
        x (array): m length, Coordinate index in longitude (or x) direction for original field val
        y (array): n length, Coordinate index in latitude (or y) direction for original field val
        val (array) : (m, n) shape, Input original field
        x1 (array): m1 length, Coordinate index in longitude (or x) direction for interpolate field
        y1 (array): n1 length, Coordinate index in latitude (or x) direction for interpolate field
        Method (string): default 'ordinary', optional in ['ordinary', 'circle', 'spline', 'spline_circle']

    Returns:
        array: (m1, n1) shape, Interpolate field

    Notes:
        x1, y1 must in the domain of x and y respectively.
        * Method, ordinary, ordinary bicubic interploation.
        * Method, circle, x1, y1 must in the domain of x and y respectively, but x1(m2) can less than x1(1)+360 and greater than x(m1). Difference with ordinary: The coordinate index must in global,such as x=(0:359.0) y=(-90:90),so we can take x(m1+1)=x(1) as a continue condition.
        * Method, spline, spline bicubic interploation
        * Method, spline_circle, x1, y1 must in the domain of x and y respectively, but x1(m2) can less than x1(1)+360 and greater than x(m1). Difference with ordinary: The coordinate index must in global,such as x=(0:359.0) y=(-90:90),so we can take x(m1+1)=x(1) as a continue condition.

    Examples:
        todo

"""
    if method == 'ordinary':
        return _interpolation.py_bicubic(x, y, val, x1, y1)
    elif method == 'circle':
        return _interpolation.py_bicubic_circle(x, y, val, x1, y1)
    elif method == 'spline':
        return _interpolation.py_bicubic_spline(x, y, val, x1, y1)
    elif method == 'spline_circle':
        return _interpolation.py_bicubic_spline_circle(x, y, val, x1, y1)


def linear_interpolation(t, x, tint, mask=None):
    """Do linear interpolation on 1 dimension field

    Args:
        t (array) : Time index for original field x
        x (array) : Original field  x
        tint (array) : Time index (must be in the domain of t) for interpolate field xout
        mask (array) : bool array, default None, if mask is not None, all the elements of x corresponding to true elements of mask will be deleted, and use other elements as the original field; if mask is None, then all x elements used as original field.

    Returns:
        array : xout, Linear interpolate field
    """
    return _interpolation.py_linear_int(t, x, tint)


def linear_int_miss(x, mask, method='mean'):
    """Do linear interpolation on 1 dimensions field (interpolate default values), more useful

    Args:
        x (array) : Original field x ( time interval of x must be equal)
        mask (array) : Logical array, if mask==.ture., all the elements of x corresponding to true elements of mask will be taken as default values, and use other elements as the original field to interpolate all default values
        method (string) : ['mean', 'near', 'zero']

    Return:
        array: Linear interpolate default values from original field

    Notes:
        case 1: method='mean' or 'MEAN'. If x1(1) or x1(n1) are true elements of mask, then x1(1) or x1(n1) will be replace by the mean value of x1 (mask=.false.) which exclude true mask elements.
        case 2: method='near' or 'NEAR'. If x1(1) or x1(n1) are true elements of mask, then x1(1) or x1(n1) will be replace by the nearest value of x1 (mask=.false.).
        case 3: method='zero' or 'ZERO'. If x1(1) or x1(n1) are true elements of mask, then x1(1) or x1(n1) will be replace by zero or 0.0d0.
"""

    conds = ['mean', 'near', 'zero']
    if method not in conds:
        raise ValueError('method should be in [%s]'.format(' '.join(conds)))
    return _interpolation.py_linear_int_miss(x, mask, method=method)


def points_average_cycle(x, weight=None, blankvalue=None):
    """Average 9 near grid points as the value of the center grid value.

    Args:
        x (array) : x[nlon, nlat], Input original grid center data, must be global data, (0-360 in longitude and -90-90 in latitude)
        weight (array) : default None, means equal weights
        blackvalue (float) : default None. If not None, then in average process, blankvalue point value will be omitted; if the average grid is blankvalue points, then set the output grid as blankvalue (no average do). If None, then all the points do 9 points average (see Algorithm for details).

    Returns:
        array: [nlon, nlat], average field of the grid center data and other 8 near grid points data

    Notes:
        grids ::
                  |            |            |            |
            --------------------------------------------------
                  |            |            |            |
                  |      1     |      2     |      3     |
                  |            |            |            |
            --------------------------------------------------
                  |            |            |            |
                  |      4     |      A     |      5     |
                  |            |            |            |
            -------------------------------------------------
                  |            |            |            |
                  |      6     |      7     |      8     |
                  |            |            |            |
            -------------------------------------------------
                  |            |            |            |


        Suppose we want to get the 9 points average value of grid center A, then we do below:

        summation 1 to 8 points and A together, then get the average value of the total 9 points as the new value of grid A. We also add area weight when do this average.
"""
    if weight is None:
        weight = np.ones(x.shape)
    flag = False if blankvalue is None else True
    return _interpolation.py_points9_ave_circle(x, weight=weight, blankvalue=blankvalue, flag=flag)
