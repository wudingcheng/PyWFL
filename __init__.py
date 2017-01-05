import _pywfl

def factorial(n):
    """ Get :math:`n!=n \dot (n-1) \cdots 3 \dot 2 \dot 1

    Args:
        n (int): number to be calculated factorial.

    Returns:
        float : get :math:`n!=n \dot (n-1) \cdots 3 \dot 2 \dot 1

    Notes:
        Use recursive method to do factorial.
        FACTORIAL is slower than WFL_FACTORIAL when call it many times.


"""
    return _pywfl.py_factorial(n)

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
        >> bilinear([1, 2], [1, 2], [1, 2, 3, 4], 1.2, 1.2)
""" 
    return _pywfl.py_bilinear(x, y, val, x1, y1)

def bilinear_square():
    return _pywfl.py_bilinear_square()


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
        * Method, circle, 
            x1, y1 must in the domain of x and y respectively, but x1(m2) can less than x1(1)+360 and greater than x(m1) .
            Difference with ordinary: The coordinate index must in global,such as x=(0:359.0) y=(-90:90),so we can take x(m1+1)=x(1) as a continue condition.
        * Method, spline, spline bicubic interploation
        * Method, spline_circle,
            x1, y1 must in the domain of x and y respectively, but x1(m2) can less than x1(1)+360 and greater than x(m1) .
            Difference with ordinary: The coordinate index must in global,such as x=(0:359.0) y=(-90:90),so we can take x(m1+1)=x(1) as a continue condition.

    Examples:
        todo

"""
    if method == 'ordinary':
        return _pywfl.py_bicubic(x, y, val, x1, y1)
    elif method == 'circle':
        return _pywfl.py_bicubic_circle(x, y, val, x1, y1)
    elif method == 'spline':
        return _pywfl.py_bicubic_spline(x, y, val, x1, y1)
    elif method == 'spline_circle':
        return _pywfl.py_bicubic_spline_circle(x, y, val, x1, y1)

def fft():
    pass


if __name__ == '__main__':
    print _pywfl.py_cfftpack.__doc__
    print _pywfl.py_cfftpackb.__doc__