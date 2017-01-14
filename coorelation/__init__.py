# encoding: utf-8
from __future__ import absolute_import
from . import _correlation

__all__ = ["cor2", "corfre", "corlag", "arcov", "edof", "fre_response", "mtm_coh", "multicor"]


def cor2(a, b, dof=0):
    """Give correlation coefficient and 95% confidential level between two vectors

    Args:
        a (array) : Array a
        b (array) : Array b
        dof (int) : Degrees of Freedom (DOF), if not equal 0, this value will be used to calculate colev; if 0, DOF will be set to n-2. If you want to do get a more reliable DOF, please use function :func:`edof` to get effect DOF and use it here.

    Returns:
        tuple:
            * co (float) : The correlation coefficient between array a and b
            * colev (float) : The 95% confidential levelarcov of correlation coefficient (Degrees of Freedom is n-2 )

    Notes:
        co is between -1 and 1

    Examples:
        >>> from numpy import cos, arange
        >>> t = arange(100)
        >>> a = cos(t)
        >>> b = cos(t - 0.3)
        >>> cor2(a, b)
"""
    return _correlation.py_cor2(a, b, dof=dof)


def corfre(a, b, dt=1.0):
    """Get  squared coherence amplitude spectrum and phase spectrum in the frequency domain for two vectors

    Args:
        a (array) : Array a
        b (array) : Array b
        dt (float) : The time interval for time index of array a and b

    Returns:
        ndarray: out(i,n-1)--i=0, frequency; i=1, periods; i=2, squared coherency amplitude spectrum; i=3, phase spectrum; i=4, 95% confidential level for amplitude spectrum

    References:
        .. [1] B. F. Chao, correlation of interannual length-of-day variation with El-nino/southern oscillation,1972-1986, J. Geophys. R., Vol.93, pp7709-7715, 1988.

    Examples:
        >>> from numpy import cos, arange
        >>> t = arange(100)
        >>> a = cos(t)
        >>> b = cos(t - 0.3)
        >>> corfre(a, b, dt=1.0)
"""
    return _correlation.py_corfre(dt, a, b)


def corlag(dt, a, b, dof):
    """Give time lagged correlation coefficients of two vectors

    Args:
        dt (float): The time interval for time index of array a and b
        a (array) : Array a
        b (array) : Array b
        dof (float) : Factor of Degrees Of Freedom (DOF) for filtered input vector a and b, if a and b is not filtered, then set  dof=1.0. If array a and b are both filtered by some filter, then the freedom of array a and b will decrease (the value of 95% correlation confidential level will increase), if you know the exact DOF of filtered array a and b, e.g. n1 (n1<n, n is the length of array a, so is array b), then set dof=n1/n will get the exact 95% confidential level; otherwise (set dof=1.0 for filtered a and b) you will get a lower value of 95% confidential level which is not correct.

    Returns:
        ndarray: out(i,-n/2+2:n/2-2)--i=1,time index; i=2, correlation coefficients; i=3, positive 95% confidential level; i=4, negative 95% confidential level

    Notes:
        When you get the maximum correlation coefficient (the most near m=0 point at both positive m and negative m) at out(2, m), if m>0, then array a lead array b; if m<0, then array a lag array b.  See example for details.

        To determine the correct DOF for filtered array a or b by using :func:`vondrak` method, please reference:

        Yan Haoming, Zhong Min, Zhu Yaozhong, Determination of the degree of freedom of digital filtered time series with an application to the correlation analysis between the length of day and the southern oscillation index, Chinese Astronomy and Astrophysics, Vol. 28, No.1, P120-126, 2004.(in English)

        闫昊明，钟敏，朱耀仲，时间序列数字滤波后自由度的确定――应用于日长变化与南方涛动指数的相关分析，天文学报，Vol. 44, No.3, P324-329, 2003. (in Chinese)

        If you want use other filter method to filter the original data, and want to get the correct DOF, please also reference the above papers to do your own Monte Carlo analysis to get the correct DOF of filtered data. Or reference function WFL_EDOF for more details for determine effect DOF.

        Note: the above definition may different with MATLAB etc. software, be aware that the above equation is more accurate. And due to the factor n/(n-lag), the r(lag) may be great than 1.0 when lag is large, please take more attention about it and only use small lag to determine the lag correlation coefficient.

        Sign convention of this routine: if a leads b, i.e., a is shifted to the left of b, then CORLAG  will show a peak at positive lags.

        For lag==0, array a(1),a(2)...a(n) and array b(1),b(2)...b(n) do corelation analysis.
        For lag==1, array a(1),a(2)...a(n-1) and array b(2),b(3)...b(n) do corelation analysis. b is shifted to left by one step.
        For lag==m(m<n/2-2), array a(1),a(2)...a(n-m) and array b(1+m),b(2+m)...b(n) do corelation analysis. b is shifted to left by m steps.
        For lag==-1, array a(2),a(3)...a(n) and array b(1),b(2),...b(n-1) do corelation analysis. a is shifted to right by one step.
        For lag==-m, array a(1+m),a(2+m),..a(n) and array b(1),b(2),...b(n-m) do corelation analysis. a is shifted to right by m steps.




"""
    return _correlation.py_corlag(dt, a, b, dof)


def arcov(x, normalized=True):
    """Calculate the auto-covariance of a time series

    Args:
        x (array) : x[n], Time series to calculate auto-covariance
        normalized (bool) : if True, then x will be normalized first (so mean(x)=0, std(x)=1), and then calculate the auto-covariance; if False, then only get the auto-covariance without normalized.

    Returns:
        array : Auto-covariance

    Notes:
        The autocorrelation function can be acquired by using r(k)=x_arcov(k)/x_arcov(0), where k=0,...,n-1. Also see ARCO or WFL_MULTICOR.

        There two algorithm to calculate autocovariance,  named Algorithm A and Algorithm B, see below for details. Algorithm B is an unbiased estimator of theoretical autocovariance of the ensemble, where Algorithm B is only asymptotically unbiased as the recorder length N tend to infinity. However, it is shown that the biased estimator (Algorithm A) has a smaller mean square error.  It is seen that the two variances of Algorithm A and B coincide when k=0 but as k tends to N the variance of the Algorithm A (biased estimator) tends to zero, whereas the variance of the Algorithm B (unbiased estimator) tends to infinity. It is this behavior which makes the unbiased estimator (Algorithm B) so unsatisfactory.  As a result, here choose Algorithm A to calculate the autocovariance.  See Jenkins and Watts (1968), p174-180, for details.

        Algorithm A:
            autocovar(K)=1/N * ( SUM( (X(T)-X_AVERAGE) * ( X(T+K)-X_AVERAGE)  ) Where T=1,N-K ) and K=0,1,...N-1

        Algorithm B (NOT used):
            autocovar(K)=1/(N-K)* ( SUM( (X(T)-X_AVERAGE) * ( X(T+K)-X_AVERAGE)  ) Where T=1,N-K ) and K=0,1,...N-1

        Reference:
            .. [1] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day,  1968.
"""
    ch = 'y' if normalized else 'n'
    return _correlation.py_wfl_arcov(x, ch=ch)


def edof(x, y):
    """Determine the Effect Degree of Freedom (EDOF) of two time series for correlation analysis

    Args:
        x (array) : Time series to be analysised
        y (array) : Time series to be analysised

    Return:
        int : The effect degree of freedom for time series x and y

    Notes:
        When do correlation analysis, we often use t test to give a 95% confidence level with DOF as n-2. But this method is only suitable for the time series which is not auto-correlated, for the time series which is auto-correlated, we should use EDOF which is smaller than DOF. So this subroutine is used  to calculate the EDOF of two time series and then use it as the DOF to do correlation analysis in subroutine :func:`cor2`.

    References:
        .. [1] Davis R. E., Predictability of sea-surface temperature and sea-level pressure anomalies over the North pacific ocean, J Phys. Oceanogr., Vol.6, 249-266,1976.

"""
    return _correlation.py_wfl_edof(x, y)


def fre_response(dt, x, y, nw, p=0.95):
    """Get Frequency Response for two time series by using Multitaper smoothing method

    Args:
        dt (float) : The time interval for time index of array a and b
        x (array) : Array x
        y (array) : Array y. :math:`y(f)=h(f)*x(f)`  in frequency domain, here :math:`h(f)` is the frequency response function which will be calculated.
        nw (float) : DPSS orders 0 to 2*nw-1 required (2<=nw<n/2). Note, here, nw is real or double type. The typical values of nw is 2, 5/2, 3, 7/2, 4. Thus you can get (nw)Pi-multitapers, e.g. for nw=4, you get 4Pi-multitapers. Smaller nw higher frequency precision, greater nw smaller variance in frequency domain.  See :func:`dpss` for more details.
        p (float) : Confidence level for the F-test (0.<p<1.0)

    Returns:
        tuple:
            * fre_res( n/2, 4) : fre_res(:,1), FFT periods (same unit as dt); fre_res(;,2), frequency response function amplitude spectrum; fre_res(:,3), frequency response function phase spectrum; fre_res(:,4), noise spectrum
            * amp_ci(n/2,2) : 100*P% confidence level confidence intervals for frequency response function amplitude spectrum fre_res(:,2). amp_ci(:,1) is lower confidence intervals, and amp_ci(:,2) is upper confidence intervals.
            * ph_ci(n/2,2) : 100*P% confidence level confidence intervals for frequency response function phase spectrum fre_res(:,3). ph_ci(:,1) is lower confidence intervals, and ph_ci(:,2) is upper confidence intervals.  Note, if the output ph_ci(i,1)==ph_ci(i,2), it means that there's no confidence interval, because the eq. (10.4.4) of Jenkins and Watts(1968) to calculate ph_ci(:,:) is out of the calculation range due to very small squared coherence between x and y. Thus, only at frequency of larger squared coherence value, the phase confidence interval is meaningful.

    Notes:
        Phase confidence intervals may be failure at some frequencies (usually at larger frequencies which close to Nyquist frequency). In this case, ph_ci(i,1)==ph_ci(i,2)==fre_res(i,3). Be careful of this failure.

    References:
        .. [1] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 10 for details of frequency response function)
        .. [2] Bell, B., Percival, D.B. and Walden, A.T. "Calculating Thomson's Spectral Multitapers by Inverse Iteration", J. Comput. and Graph. Stat., 1993.
        .. [3] Percival D. B., and A. T. Walden, Spectral analysis for physical applications: Multitaper and conventional univariate techniques, Cambridge University Press, Cambridge, 583p, 1993.
        .. [4] Thomson, D. J., 1982, Spectral estimation and harmonic analysis: IEEE Proc., v. 70, no. 9, p. 1055-1096.
"""
    return _correlation.py_wfl_fre_response(dt, x, y, nw, p=p)


def mtm_coh(dt, x, y, nw, p=0.95):
    """Get  squared coherence amplitude spectrum and phase spectrum / cross amplitude and phase spectrum in the frequency domain for two vectors by using multitaper(MTM) method

    Args:
        dt (float) : The time interval for time index of array a and b
        x (array) : Array x
        y (array) : Array y. :math:`y(f)=h(f)*x(f)`  in frequency domain, here :math:`h(f)` is the frequency response function which will be calculated.
        nw (float) : DPSS orders 0 to 2*nw-1 required (2<=nw<n/2). Note, here, nw is real or double type. The typical values of nw is 2, 5/2, 3, 7/2, 4. Thus you can get (nw)Pi-multitapers, e.g. for nw=4, you get 4Pi-multitapers. Smaller nw higher frequency precision, greater nw smaller variance in frequency domain.  See :func:`dpss` for more details.
        p (float) : Confidence level for the F-test (0.<p<1.0)

    Returns:
        tuple:
            * coh_spec( n/2, 3):coh_spec(:,1), FFT periods (same unit as dt); coh_spec(;,2), squared coherency amplitude spectrum; coh_spec(:,3), coherency phase spectrum
            * cross_spec( n/2, 3):cross_spec(:,1), FFT periods (same unit as dt); cross_spec(;,2), cross amplitude spectrum of a and b; cross_spec(:,3), cross phase spectrum
            * f_test(n/2):100*P% confidence level F-test value for squared coherency amplitude spectrum

    Notes:
        Similar as subroutine :func:`corfre` bus use multitaper method to calculate the spectrum.

    References:
        .. [1] Bell, B., Percival, D.B. and Walden, A.T. "Calculating Thomson's Spectral Multitapers by Inverse Iteration", J. Comput. and Graph. Stat., 1993.
        .. [2] Percival D. B., and A. T. Walden, Spectral analysis for physical applications: Multitaper and conventional univariate techniques, Cambridge University Press, Cambridge, 583p, 1993.
        .. [3] Thomson, D. J., 1982, Spectral estimation and harmonic analysis: IEEE Proc., v. 70, no. 9, p. 1055-1096.
        .. [4] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 6-9)

    See Also:
        :func:`corlag`, :func:`cor2`, :func:`corfre`
"""
    return _correlation.py_wfl_mtm_coh(dt, x, y, nw, n=len(x), p=p)


def multicor(x, out_type='cor', fft=False, remove_mean=True):
    """Calculate the auto- and cross covariance/correlation of input data (time series)

    Args:
        x (ndarray) : x[n,m], total m time series, each has n points. Time series for calculating auto- and cross covariance/correlation. x was first removed average in the subroutine.
        out_type : out_type='cov' or 'cor' for covariance and correlation output
        fft (bool) : if True, use fft method
        remove_mean (bool) : if True, then first remove mean of input x(n,m) for each column m, if False, no remove mean for the input data x(n,m). If not present  remove_mean, no mean removed from input data x(n,m).

    Returns:
        ndarray : xout
            Auto- and cross covariance for X if out_type='cov';
            auto- and cross correlation for X, if out_type='cor'
            The first index of x(:,:,:) is lags. Thus x_out(0,m,m) is the covariance/correlation of all m series at lag==0. And x_out(:,i,i) is the autocovariance/autocorrelation of the ith series.

            when using FFT method, The first index of output x_out(:,m,m)  is different from not using fft.

    Notes:
        Properties of covariance functions:
            The properties of the auto covariance functions of a real bivariate process are the same as those for the auto covariance functions of a univariate  process, that is,

            .. math::
                \\gamma_{ii}(0) & = Var\\left[ X_i(t) \\right ] = \\sigma^2_{X_i} \\\\
                \\gamma_{ii}(u) & = \\gamma_{ii}(-u)

        Hence the covariance between the two stochastic processed can be described by means of the single cross covariance functions :func:`\\gamma_{12}(u)`, where :func:`-\\infty \\leq u \\leq +\\infty`. Note that although the auto covariance function is an even function of lag, the cross covariance function will **NOT** be an even function in general.

        The cross correlation function:

            In general it may be necessary to study the interactions between two processes with possibly different scales of measurement or different variances. In this situation it is necessary to define the cross correlation function

                .. math::
                    \\rho_{12}(u) = \\gamma_{12}(u) / \\sqrt(\\gamma_{11}(0)\\gamma_{22}(0)) = \\gamma_{12}(u)/(\\sigma_1 \\sigma_2)

            The first property of the cross correlation function is that

                .. math:: | \\rho_{12}(u) | \leq 1

            which follows from the fact that the random variable

                .. math:: Y(t) = \\lambda_1 X_1(t) + \\lambda_2 X_2(t+u)

            has positive variance.
            The second property is that

                ..math:: \\rho_{12}(u) = \\rho_{21}(-u)

        Sign convention of this routine: if x(:,i) lags x(:,j), i.e., x(:,i) is shifted to the right of x(:,j), then WFL_MULTICOR  will show a peak at positive lags.

        Difference between :func:`multicor` with :func:`corlag`:

        The results of :func:`multicor` (for auto and cross correlation only) is different with the results of :func:`corlag` in three parts: first, :func:`corlag` has a factor of n/(n-lag) while :func:`multicor` has the factor 1.0;  second, :func:`corlag` use different variance (each n-lag data) while :func:`multicor` use the same variance cov(0,i,i); and third, the sign convention of func:`corlag` is reverse with that of :func:`multicor`.

    Examples:
        >>> import numpy as np
        >>> x = np.array([[47, 64, 23, 71, 38, 65, 55, 41, 59, 48, 71, 35, 56, 40, 58, 44, 80, 55,
                       37, 74, 51, 58, 50, 60, 44, 57, 50, 45, 25, 59, 50, 71, 56, 74, 50, 58,
                       45, 54, 36, 54, 48, 55, 45, 57, 50, 62, 44, 64, 43, 52, 38, 60, 55, 41,
                       53, 49, 34, 35, 54, 45, 68, 38, 50, 60],
                      [47, 64, 23, 71, 38, 65, 55, 41, 59, 48, 71, 35, 56, 40, 58, 44, 80, 55,
                       37, 74, 51, 58, 50, 60, 44, 57, 50, 45, 25, 59, 50, 71, 56, 74, 50, 58,
                       45, 54, 36, 54, 48, 55, 45, 57, 50, 62, 44, 64, 43, 52, 38, 60, 55, 41,
                       53, 49, 34, 35, 54, 45, 68, 38, 50, 60]]).T
        >>> multicor(x, fft=True)[:, 0, 1]
        >>> multicor(x)[:, 0, 1]

    References:
        .. [1] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications,  P207-208, P322-325, Holden-Day, 1968.

"""
    remove_mean = 'y' if remove_mean else 'n'
    return _correlation.py_wfl_multicor(x, out_type, fft=fft, remove_mean=remove_mean)
