# encoding: utf-8
from __future__ import absolute_import
from . import _statistics


__all__ = ["aic", "bic", "chinv", "cumsum", "diff", "geop", "leap_rm", "montonic",
           "montecarlo", "month_mean", "month_mean_idx", "var_percent", "cdff", "cdfnor",
           "cdft", "taylor_diagram"]


def aic(m, residuals):
    """Calculate the A (or minimum) Information Criterion (AIC) in order to select the most suitable parameter of polynomial's order (the minimum AIC order is best)

    Args:
        m (int) : The numbers of normal equation (in HARMONIC subroutine m=n)
        residuals (array) : residuals[n], The standard variation of residual series (different between fit polynomial and original data)
        max_order (int) : The maximum polynomial order to be fitted


    Returns:
        tuple
            * rs (array) : rs[n], the AIC value AIC(q), see algorithm for details
            * order (int) : The best polynomial fit order under AIC criterion

    Notes:
        The fit model is :math:`y(t) = a_i t^{i-1}`, where :math:`i=1,2,\dots` here :math:`i` means :math:`i-1` order polynomial.

        n is between N/3 to 2N/3 for ARMA time series (N is length of time series)

        .. math:: AIC(q) = \\log \\left( \\sigma \\right) + 2q/n

        or

        .. math:: AIC(q) = n \\log \\left( \\sigma \\right) + 2q

        where :math:`\\sigma` is standard deviation of residual series, :func:`q` is ARMA model's order.

    See Also:
        :func:`bic`
"""
    return _statistics.py_aic(m, residuals)


def bic(m, residuals):
    """Calculate the Bayesian extension of the AIC (BIC) in order to select the most suitable parameter of polynomial's order (the minimum BIC order is best)

    Args:
        m (int) : The numbers of normal equation (in HARMONIC subroutine m=n)
        residuals (array) : residuals[n], The standard variation of residual series (different between fit polynomial and original data)
        max_order (int) : The maximum polynomial order to be fitted


    Returns:
        tuple
            * rs (array) : rs[n], the AIC value AIC(q), see algorithm for details
            * order (int) : The best polynomial fit order under AIC criterion

    Notes:
        The fit model is :math:`y(t) = a_i t^{i-1}`, where :math:`i=1,2,\dots` here :math:`i` means :math:`i-1` order polynomial.

        n is between N/3 to 2N/3 for ARMA time series (N is length of time series)

        .. math:: AIC(q) = \\log \\left( \\sigma \\right) + q \\log \\left( n \\right) /n

        or

        .. math:: AIC(q) = n \\log \\left( \\sigma \\right) + q \\log \\left( n \\right)

        where :math:`\\sigma` is standard deviation of residual series, :func:`q` is ARMA model's order.

        The difference between :func:`aic` and func:`bic` is the 2nd term of the equation, in func:`bic` use :math:`log(n)` and in :func:`aic` use 2.0. Generally, :func:`log(n)>>2.0`, so the best order got by func:`bic` is less than by :func:`aic`. When n is infinite, the order determined by :func:`aic` is greater than the true order, but the func:`bic` can get the true order.

"""
    return _statistics.py_bic(m, residuals)


def chinv(alpha, dof):
    """Calculate Chi-square cumulative distribution functions from given significance level or confidence level and Degrees of Freedom (DOF)

    Args:
        alpha (float) : Significance level. alpha typically =0.01,0.05,0.10. Confidence level :math:`p=1-\\alpha \\quad (0.0<p<1.0)`, typically p=0.90, 0.95, 0.99.
        dof (float) : Degrees of Freedom for Chi-square

    Return:
        float : Chi-square value calculated from P and DOF

    Notes:
        :func:`chinv` is modified from DCDFLIB: Available at http://www.netlib.org/random/

"""
    p = 1 - 0.5 * alpha
    return _statistics.py_chinv(p, dof)


def cumsum(x, seed=0.0):
    """cumulative sum on an array, with optional additive seed

    Args:
        x (array) : Array will be  summed cumulatively
        seed (float) : The add constant of x(1)

    Return:
        array : The cumulative sum of array arr, has the same length with x

    Notes:
        Given the rank 1 array arr , returns an array of identical type and size containing the cumulative sums of arr. If the optional argument seed is present, it is added to the first component (and therefore, by cumulation, all components) of the result.

    Examples:
        >>> x = [1, 2, 3, 4, 5]
        >>> cumsum(x) # [  1.   3.   6.  10.  15.]
        >>> cumsum(x, seed=3) # [  4.   6.   9.  13.  18.]

"""
    return _statistics.py_cumsum(x, seed=seed)


def diff(x):
    """Forward differential for a time series

    Arg:
        x (array) : Array x, time series to be forward differential

    Return:
        array : [n-1], Forward differential time series

    Notes:
        forward differential: :math:`xout_k=x_{k+1}-x_k`
        backward differential: :math:`xout_k=x_k-x_{k-1}`
        central differential: :math:`xout_k=x_{k+1}-x_{k-1}`
"""
    return _statistics.py_diff(x)


def geop(first, factor, n):
    """Arguments

    Args:
        first (float) : The first value of geometic progression
        factor (float) : The factor of geometic progression
        n (int) : The number of geometic progression

    Return:
        array : The array of geometic progression

    Notes:
        Returns an array of length n containing a geometric progression whose first value is first and whose multiplier is factor. If first and factor have rank greater than zero, returns an array of one larger rank, with the last subscript having size n and indexing the progression.

        .. math::
            geop_1 & = first \\\\
            geop_k & = geop_{k-1} \dot factor

"""
    return _statistics.py_geop(first, factor, n)


def leap_rm(x, min_leap):
    """Remove leap values from a time series

    Args:
        x (Array): Time series to be remove leap values
        zleap_min : The minimum leap values (positive) taken as leap point
    Return:
        array : Time series remove leap points and remove averages

    Notes:
        Some times, a time series my have some leap (or jump) points, which cause the time series is not continue, so these leap points should be removed from the original field. Then what method should be used to do this? The simple method is to get forward differential first, then check the forward differential values, if the absolute values of forward differential greater than a threshold (here denoted by zleap_min), then we take this point as a leap points. After all leap points are find,  we remove the average between two leap points and get the continue time series.
        If the two continue points are both leap points, the output data may have leap points too. So in this instance, you must combine the two points as one then use :func:`leap_rm`.
"""
    return _statistics.py_leap_rm(x, min_leap)


def montonic(x):
    """Determine if an array is monotonic

    Arg:
        x (array) : Must be a real/double array, one dimensions

    Return:
        bool : If True, then the array x is monotonic.
"""
    return bool(_statistics.py_monotonic(x))


def montecarlo(n, ivars=2):
    """Give critical coherence coefficient at different confidential level (90% 95% and 99%) for multi-variables (:math:`\geq 2`)

    Args:
        n (int) : Number of data (or Degree of Freedom, DOF=n-ivars ) (n>ivars)
        ivars (int) : Total variables (ivars>=2)

    Returns:
        tuple : co90, co90s, co95, co95s, co99, co99s
            * co90, co95, co99 : Critical coherence coefficient of corresponding confidential level (90%,95%,99%)
            * co90s, co95s, co99s : Standard deviation of co90,co95,co99
"""
    return _statistics.py_montecarlo(n, ivars=ivars)


def month_mean(x, year, month, total_month):
    """Get monthly average field from daily field

    Args:
        x (array) : Daily field to be monthly averaged
        year (int) : Year for the first data
        month (int) : Month for the first data
        total_month (int) : Total monthes for daily field x

    Return:
        array: Monthly average field from daily field x

    Notes:
        The first data x[0] must be the field of begin day for any month, the last data x[n-1] must be the end day for any month. That means no omit data in any month for daily field. If the begin day or end day do not correspond with the above law, then return errors.
"""
    return _statistics.py_month_mean(x, year, month, total_month)


def month_mean_idx(t, total_month, year_begin, month_begin, year_end, month_end):
    """Get monthly field index from daily or sub-month field

    Args:
        t (array) : t[n], The time index of original field (at least 2 values in each month), must be monotonic and in unit of Modified Julian Day (MJD).
        total_month : The number of months between begin and end month (n>2*total_month)
        year_begin (int) : Begin year (4 digitals)
        month_beign (int) : Begin month
        year_end (int) : End year (4 digitals)
        month_end (int) : End month

    Returns:
        tuple
            * idx_begin (array) :  idx_begin[total_month], Begin index number of ith month in t
            * idx_end (array) :  idx_end[total_month], End index number of ith month in t

    Notes:
        If we have a field with 10 days interval, we want to get monthly averaged field, what should we do? We know, there are may be 3 or 4 data field in one month, but we do not know exactly in which month we have 3 or 4 data. Then we can find it by MONTH_MEAN_IDX. First, we must have the time index (MJD unit) of 10 days interval field, then we can get the data index (begin and end index) for each month in original field. Now, we can do average easily.  See Example for details.
"""
    return _statistics.py_month_mean_idx(t, total_month, year_begin, month_begin, year_end, month_end)


def var_percent(a, b):
    """Calculate "b can explained variance percent of a "

    Args:
        a (array) : array a
        b (array) : array b,  len(a) == len(b)

    Return:
        float : Percent of variance explained by b for a

    Notes:
        This subroutine only give the one signal variance percent explained by another signal. Please note, the results maybe mis-interpreted. For example, if a and b are the same cosine function with period 365 days but only different with phase initial.

        +--------------------+--------------------------------------------+
        |phase difference    |  variance explained (:func:`var_percent`)  |
        +--------------------+--------------------------------------------+
        |:math:`10^\circ`    |         97%                                |
        +--------------------+--------------------------------------------+
        |:math:`20^\circ`    |         87%                                |
        +--------------------+--------------------------------------------+
        |:math:`30^\circ`    |         71%                                |
        +--------------------+--------------------------------------------+
        |:math:`45^\circ`    |         37%                                |
        +--------------------+--------------------------------------------+
        |:math:`60^\circ`    |          0%                                |
        +--------------------+--------------------------------------------+
        |:math:`90^\circ`    |         -115%                              |
        +--------------------+--------------------------------------------+

        If the amplitude are factor: (e.g. b=0.9b, 0.8b, 0.5b), which means different amplitude and 10deg phase difference.

        +-----------------+-----------------+--------------------------------------------+
        |phase difference | amplitude factor|  variance explained (:func:`var_percent`)  |
        +-----------------+-----------------+--------------------------------------------+
        | :math:`10^\circ`|        0.9      |                  96%                       |
        +-----------------+-----------------+--------------------------------------------+
        | :math:`10^\circ`|        0.8      |                  93%                       |
        +-----------------+-----------------+--------------------------------------------+
        | :math:`10^\circ`|        0.5      |                  73%                       |
        +-----------------+-----------------+--------------------------------------------+

        Thus, the results must be interpreted carefully.

        The algorithm is:

            .. math:: varpercent=1-(var(a-b))/var(a))

    Examples:
        >>> import numpy as np
        >>> t = np.arange(100)
        >>> a = np.cos(t) + 0.00002 * t
        >>> b = 0.25 * a
        >>> var_percent(a, b) # 43.75%
        >>> b = a
        >>> var_percent(a, b) # 100%
"""
    return _statistics.py_varpercent(a, b)


def cdff(p, dof1, dof2):
    """Calculate F cumulative distribution functions from given confidence level P and Degrees of Freedom

    Args:
        p (float) : Confidential level P (0.0<P<1.0), typically p=0.90, 0.95, 0.99
        dof1 (float) : Degrees of freedom of the numerator sum of squares or the first variable.
        dof2 (float) : Degrees of freedom of the denominator sum of squares or the second variable.

    Return:
        float : F cumulative distribution value calculated from p and dof1, dof2.

"""
    return _statistics.py_wfl_cdff(p, dof1, dof2)


def cdfnor(p):
    """Calculate Normal (Gauss) cumulative distribution functions from given significance level or confidence level

    Arg:
        p (float) : Significance level, typically p=0.90, 0.95, 0.99.

    Return:
         float: Normal (Gauss) cumulative distribution value at given significance level

    Notes:
        :func:`cdfnor` is modified from DCDFLIB: Available at http://www.netlib.org/random/
        Normal(1-alpha/2)   where alpha=1-p, and p typically =0.90, 0.95, 0.99
"""
    alpha = 1 - p
    return _statistics.py_wfl_cdfnor(1 - 0.5 * alpha)


def cdft(p, dof):
    """Calculate Student's t cumulative distribution functions from given significance level or confidence level and Degrees of Freedom

    Args:
        p (float) : Significance level, typically p=0.90, 0.95, 0.99.
        dof (float) : Degrees of freedom

    Return:
        float : Student's t cumulative distribution value calculated from p and dof.

    Notes:
        :func:`cdft` is modified from DCDFLIB: Available at http://www.netlib.org/random/

"""
    alpha = 1 - p
    return _statistics.py_wfl_cdft(1 - 0.5 * alpha, dof)


def taylor_diagram(r, f, isnormal=True):
    """Get Taylor diagram parameters

    Args:
        r(n) : the reference time series
        f(n) : the model time series

    Returns:
        array : out[6], output Taylor diagram parameter.

            * out(1:3)==>standard deviation (std) of r(n), std of f(n), correlation coefficient between r(n) and f(n) == cos(fai)
            * out(4:6)==>E' (centered RMS difference), sigma_f*cos(fai)== x value in plot, sigma_f*sin(fai)== y value in plot

    Notes:
        out(1)=std_n(r) !standard deviation by use std_n subroutine
        out(2)=std_n(f)
        out(3)=cor2(r,f) !correlation coefficient by use cor2 subrotine
        out(4)**2=out(1)**2+out(2)**2-2*out(1)*out(2)*out(3)
        out(5)=out(2)*cos(fai)  !cos(fai)==out(3), x coordiante used by Taylor.diagram.bas for plotting in Grapher10
        out(6)=out(2)*sin(fai)   ! y coordinate used by  Taylor.diagram.bas for plotting in Grapher10

    References:
        .. [1] Taylor K. E., Summarizing multiple aspects of model performance in a single diagram, J. Geophys. Res., Vol.106, P7183-7192,2001.
"""
    isnormal = 'y' if isnormal else 'n'
    return _statistics.py_wfl_taylor_diagram(r, f, isnormal=isnormal)


# print geop(1, 1.2, 5)
