# encoding: utf-8
from __future__ import absolute_import
from . import _fit

__all__ = ["fit", "fitmulti", "fitpoly", "fitpoly_best_degree", "detrend", "lsfit"]


def fit(x, y):
    """ Give a linear fit of a time series by minimizing chi-square method.

    Args:
        x (array): The time index of array y
        y (array): Time series to be fitted (:math:`y=a+bx`)

    Returns:
        tuple:
            * a (float): constant term for :math:`y=a+bx`
            * b (float): linear (or velocity) term for :math:`y=a+bx`
            * sig_a (float): Respective probable uncertainties of constant term a.  (Note: Formal error, unrealistic)
            * sig_b (float): Respective probable uncertainties of linear term b. (Note: Formal error, unrealistic)

    Notes:
        Given a set of data points x(1:n),y(1:n), fit them to a straight line y = a + bx by minimizing chi-square. Returned are a,b and their respective probable uncertainties siga and sigb.

    Examples:
        >>> import numpy as np
        >>> a = np.arange(10)
        >>> b = 1.2 * a + 3.0 + np.random.rand(10)
        >>> fa, fb, sig_a, sig_b = fit(a, b)

"""
    return _fit.py_fitline(x, y)


def fitmulti(x, y):
    """Do multiple linear regression by minimizing chi-square method (or least squares method).

    Args:
        x (array) : x[n, m-1] The multi-variables value (n point, and m-1 variables)
        y (array) : The function value :math:`y(i)=a0+a1*x1(i)+a2*x2(i)+...+a(m-1)*x(m-1)(i)`, where :math:`i=1,n`

    Returns:
        tuple
            * a (array) : a(m) Multiple linear regression coefficients of :math:`y(i)=a0+a1*x1(i)+a2*x2(i)+...+a(m-1)*x(m-1)(i)`
            * covar (matrix): covar(m, m) Covariance matrix (the diagonal is variance) of a
            * chisq (float): Mean standard deviation, see Notes for details
            * s (float): Mean standard deviation, see Notes for details
            * r (float): Coefficient of multiple correlation, see Notes for details
            * v (array): v(m) Partial correlation coefficient, see Notes for details

    Notes:
        Please careful to deal with the relation between regression coefficients and variable x. a(1) corresponding x0 (the constant term), a(2:m) corresponding to x(n, 1:m-1) variables, so covar and v. Also, v(1)=0.0 means no partial correlation coefficient (for constant term) was calculated from FITMULTI.
        The multiple linear regression model is:

        .. math:: y_i = \sum_{j=0}^{m} a_{ij}x_{ij} \qquad i = 1, 2, \dots n

        Residuals sum of squares is:

        .. math:: \chi^2 = \sum_{i=1}^{n} [ y_i - \sum_{j=0}^{m} a_{ij} x_{ij} ]

        Mean standard deviation is:

        .. math::
            s & = & \sqrt{1 - \chi^2 / dyy} \\\\
            dyy & = & \sum_{i=1}^{n} y_i - \bar{y} \\\\
            \bar{y} & = & \sum_{i=1}^{n} y_i / n

        Partial correlation coefficient is:

        .. math::
            v_j & = & \sqrt{1 - \chi^2 / Q_j} \\\\
            Q_j & = & \sum_{i=1}^{n}[ y_i - \sum_{k=1, k \neq j}^{m}a_{ik}x_{ik} ]

    Examples:
        >>> import numpy as np
        >>> a = np.arange(10)
        >>> b = 1.2 * a + 3.0 + np.random.rand(10)
        >>> fitmulti(a, b)
"""
    return _fit.py_fitmulti(x, y)


def fitpoly(x, y, order=2):
    """Give a polynomial fit of a time series by minimizing chi-square method.

    Args:
        x (array): The self-variable of fitted polynomial function
        y (array): Time series to be fitted ( y=f(x) )
        order (int) : default 2, The maximum degree of polynomial (you can use FITPOLY_BEST_DEGREE to get the best degree of polynomial)

    Returns:
        tuple
            * a (m): Polynomial term for :math:`y = a1+a2*x+a3*x**2+...am*x**(m-1)`. Note: After Dec. 19, 2014, the x in the formula should be substituted by x-x(1), see comments for details.
            * covar (m, m): Covariance matrix (the diagonal is variance) of a
            * chisq: :math:`\chi^2`,the postfit residual sums of squares. the postfit residual sums of squares.  chisq=sum [ (y(i)-y_fit(i)) **2 ]  ,i=1,...,n

    Notes:
        Given a set of data points x(1:n), y(1:n) , use chi-square minimization to fit for some or all of the coefficients a(1:m) of a function that depends linearly on a, y = a1+a2*x+a3*x**2+...,The program returns values for a(1:m), :math:`\chi^2` = chisq, and the covariance matrix covar(1:m,1:m).

        the input x(:) was changed to x(:)-x(1) in the subroutine FITPOLY on Dec. 19, 2014. To do this will prevent the problem in x(:) which may cause wrong output values of the subroutine. For example, if x(:) is index of year, i.e., 2000.001-2010.999 with daily interval, using the NON-modified FITPOLY will give wrong results, but the modified FITPOLY will give the correct results (the reason is that we will use x**n in the subroutine, if x value change is small, x**n will similar.). After this modification, to restore the polynomial data using output parameter a manually, you must change x(:) to x(:)-x(1). That is y=a0+a1*(x-x(1))+a2*(x-x(1))**2+... , or simply use subroutine POLYNOMIAL.

    Examples:
        >>> import numpy as np
        >>> a = np.arange(10)
        >>> b = 1.2 * a + 3.0 + np.random.rand(10)
        >>> fitpoly(a, b)
"""
    if order < 2:
        order = fitpoly_best_degree(x, y)
    return _fit.py_fitpoly(x, y, order)


def fitpoly_best_degree(x, y, max_order=10, method='bic'):
    """Give the best fit degree (or order) of a polynomial fit for a time series by using minimizing chi-square method.

    Args:
        x (array) : The time index of array y
        y (array) : Time series to be fitted ( y=f(x) )
        n (int): The number of array x and y
        max_order (int) : The maximum degrees of polynomial fit, default 10
        method(string) : Method to get best fit order, default 'bic' criterion; optional in ['aic', 'bic']

    Returns:
        int : nbest_order, The best order for the polynomial fit

    Notes:
        Given a set of data points x(1:n), y(1:n) , use chi-square minimization to fit for some or all of the coefficients a(1:m) of a function that depends linearly on a, y = a1+a2*x+a3*x**2+...,then use AIC or BIC criterion to get the best fit order.

    Example:
        >>> import numpy as np
        >>> a = np.arange(10)
        >>> b = 1.2 * a + 3.0 + np.random.rand(10)
        >>> fitpoly_best_degree(a, b, method='aic')
"""
    cond = ['aic', 'bic']
    if method not in cond:
        raise ValueError('cond should be in "[%s]"'.format(' '.join(cond)))
    return _fit.py_fitpoly_best_degree(x, y, max_order, ch=method)


def detrend(x, y):
    """Detrend the linear term of a time series

    Args:
        x (array): The time index of array y
        y (array): Array y, time series to be detrend

    Returns:
        array: yout, detrend time series

    Notes:
        Least squares fit
        Gives data x(n),y(n), get detrend data yout(n), the trend is fitted as y=a+bx

    Examples:
        >>> import numpy as np
        >>> a = np.arange(10)
        >>> b = 1.2 * a + 3.0 + np.random.rand(10)
        >>> r = detrend(a, b)
"""
    assert x.shape == y.shape
    return _fit.py_detrend(x, y)


# def harmonic(t, x, pe, itrend=2, timebegin=None, method='bic'):
#     if timebegin is None:
#         timebegin = x[0]
#     if itrend < 0:
#         return _fit.py_wfl_harmonic_polys(t, x, pe, timebegin=timebegin, method=method)
#     else:
#         return _fit.py_wfl_harmonic_fix(t, x, pe, itrend, timebegin=timebegin)

def lsfit(x, y, m, func, sig=None):
    """General linear Least Squares(LS) fits

    Args:
        x (array) : x[n], The time series of self-variable of LS fit function, e.g. function :math:`f(x)`
        y (array) : y[n], Time series to be fitted ( :math:`y=f(x)` )
        m (int) : The number fitted coefficients
        func (function) : External function ( :math:`y=f(x)` ) to be fitted, func must be a subroutine to provide the fit function, it could be considered as the design matrix, this function returns ndarray result[m, n]
        sig (array) : sig[n], default None, means all weight are the same.

    Returns:
        tuple
            * a (array) : a(m), LS fitted coefficients
            * cov (ndarray) : cov[m, m], Covariance matrix (the diagonal is variance) of a
            * chisq (float) : :math:`\chi^2`, the postfit residual sums of squares. :math:`\chi^2 = \sum ( y_i - f(x_i)  )^2, i=1,\dots,n`

    Notes:
        General Linear LS Fit by Gauss-Jordan Elimination.

        Given a set of data points x(1:n), y(1:n) , use chi-square minimization to fit all of the coefficients a(1:m) of a function that depends linearly on a, y = f(x),The program returns values for a(1:m), :math:`\chi^2`, and the covariance matrix covar(1:m,1:m).

        How to write user-defined function? the user-defined function can be regarded as design matrix, see the example


    Examples:
        >>> import numpy as np
        >>>
        >>> # user define function f(x) = 1 + x + x^2
        >>> # the result will return the left design matrix
        >>> def fpoly(x, n):
        >>>     res = np.ones((n, len(x)))
        >>>     for i in range(1, n):
        >>>         res[i, :] = res[i - 1, :] * x
        >>>     return res
        >>>
        >>>
        >>> def func(x):
        >>>     return 1 + x + x**2
        >>>
        >>>
        >>> x = np.arange(10) * 1.0
        >>> y = func(x)
        >>> print lsfit(x, y, 3, fpoly)
        >>>
        >>> with weight example
        >>> yerr = np.random.rand(10)
        >>> y += yerr
        >>> print fpoly(x, 3).shape
        >>> print lsfit(x, y, 3, fpoly, sig=yerr)
"""
    if sig is not None:
        return _fit.lsfit(x, y, m, func, sig=sig)
    else:
        return _fit.lsfit(x, y, m, func)


# def fpoly(x, n):
#     res = np.ones((n, len(x)))
#     for i in range(1, n):
#         res[i, :] = res[i - 1, :] * x
#     return res


# def func(x):
#     return 1 + x + x**2


# if __name__ == '__main__':
#     import matplotlib.pyplot as plt
#     x = np.arange(10) * 1.0
#     y = func(x)
#     print lsfit(x, y, 3, fpoly)
#     yerr = np.random.rand(10)
#     y += yerr
#     print fpoly(x, 3).shape
#     print lsfit(x, y, 3, fpoly, sig=yerr)
    # print _fit.py_wfl_harmonic_polys.__doc__
