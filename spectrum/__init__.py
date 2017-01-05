# encoding: utf-8
import _spectrum


def arco(x, m=0):
    """Get auto-regressive (AR) process's parameter and autocorrelation coefficients

    Args:
        x (array) : Array or time series to be analysis
        m (int) : Maximum model order given by user (or maximum lag order, e.g. AR(m) )
            if m > 0, use m as the order of AR model; if m <= 0,  then using FPE criterion to get the best order of AR model

    Returns:
        tuple
            * co(0:m) : Array of coefficients of AR model estimated by the subroutine
            * cor(0:m) : Array of autocorrelation function
            * best_order :(Optional) Must be presented with ch at the same time. If presented, return the best order of AR model.

    Notes:
        cor(0)==1.0; cor(1:m)<1.0
        To restore the original time series xx(:) with the AR model coefficients co(0:m), use the following code:

        do j=1,n
            do i=1,m  !here m can be substituted by best_order
                if ( (j-i)>=1) then
                    xx(j)=xx(j)-co(i)*x(j-i)  !note here is -co(i), co(0)==1, and should not be used.
                endif
            enddo
        enddo

    References:
        .. [1] Ding Y., W. Zheng, Astronomical data processing, Nanjing Univ.  Press, Nanjing, 1990. (in Chinese)
        .. [2] 丁月蓉, 郑大伟 , 天文测量数据的处理方法, 南京大学出版社,1990.

    Examples:
        todo

"""
    if m > 0:
        return _spectrum.py_arco(m, x)[:2]
    else:
        return _spectrum.py_arco(1, x, ch='fpe')


def emd(t, x, m=10, std=0.3):
    """Computes Empirical Mode Decomposition (EMD)

    Args:
        t (array) : time index of signal :math:`X`
        x (array) : input signal
        m (int) : number of Intrinsic Mode Functions(IMFs) (last one is residual)
        std (float) : Standard deviation to stop the iteration (0.2<=std<=0.3),default is 0.3

    Returns:
        array : IMFs (first to m-1 column is IMFs, last one (m) is residual)

    Notes:
        Edge conditions: extended symmetrical 2 points of both end maximum values and minimum values respectively.
        EMD_TOL is more speed than EMD, but EMD may more accuracy (personal opinion)

    References:
        .. [1] N. E. Huang et al., "The empirical mode decomposition and the Hilbert spectrum for non-linear and non stationary time series analysis," Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998

    Examples:
        todo
"""
    return _spectrum.py_emd(t, x, m, std=std)


def emd_tol(t, x, m=12):
    """Computes Empirical Mode Decomposition (EMD)

    Args:
        t (array) : time index of signal :math:`X`
        x (array) : input signal
        m (int) : number of Intrinsic Mode Functions(IMFs) (last one is residual)

    Returns:
        array : IMFs (first to m-1 column is IMFs, last one (m) is residual)

    Notes:
        Edge conditions: extended symmetrical 2 points of both end maximum values and minimum values respectively.
        Same as EMD but with different stopping criterion for sifting : at each point : mean amplitude < threshold2*envelope amplitude & mean of boolean array ((mean amplitude)/(envelope amplitude) > threshold) < tolerance & the N.E. HUANG's stopping criterion sd0>0.3
        EMD_TOL is more speed than EMD, but EMD may more accuracy (personal opinion)

    References:
        .. [1] N. E. Huang et al., "The empirical mode decomposition and the Hilbert spectrum for non-linear and non stationary time series analysis," Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998


"""
    return _spectrum.py_emd(t, x, m)


def fft_spec(x, dt=1, p=0.95, lag1_co=0):
    """Get normalized Fourier spectrum, red-noise background spectrum, and 100p% (e.g. p=95%,99%) confidence spectrum of time series.

    Args:
        x (array) : Time series to be analyzed
        dt (float) : Time interval of x
        p (float) : Confidential level (0<p<1.0, typically p=0.95, p=0.99)
        lag1_co (float) : The lag-1 autocorrelations of x, used to output background red-noise spectrum, see FNS for details. If present, lag1_co as the lag-1 autocorrelations; if not present, use ARCO to get it, and then calculate FNS.

    Returns:
        array:　spec(:,4) :　Normalized Fourier spectrum. spec(:,1), period index; spec(:,2), Fourier spectrum of x; spec(:,3) background Red-noise spectrum; spec(:,4), 100p% confidence spectrum

    Reference:
        ..[1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16), (17)

"""
    return _spectrum.py_fft_spec(dt, x, p, lag1_co=lag1_co)


def fft_fre(x, dt=1, amp=True):
    """Get the period, amplitude and instantaneous phase (or real and image part) of a vector by FFT method

    Args:
        x (array) : array (or time series) to be analyzed
        dt (float): time interval of array x
        amp (bool): if true, then get the amplitude and instantaneous phase of x; if false, then get real and image part (FFT result) of x.

    Returns:
        tuple
            * pe (array): period of time series x
            * am (array): Amplitude or real part (FFT result) of time series x corresponding with period pe.
            * ph (array): Instantaneous phase or image part (FFT result) of time series x corresponding with period pe.

    Notes:
        FFT_FRE is almost the same as FFT_PERIOD, but FFT_FRE can get instantaneous phase or it can get the real and image part of time series corresponding with the period analysis. Some times, we need this corresponding relation. E.G. , we have two time series x and y, and we want to know the spectrum of :math:`z=y/x`, then we can get  the spectrum (real and image part) of x and y first, then get :math:`z(\omega)=y(\omega)/x(\omega)`, now we get the spectrum of z.
        When you get the real and image part of a time series, you can get the energy spectrum of x by :math:`\sqrt(real^2+imag^2)`, but this is not the amplitude spectrum.  The energy spectrum and amplitude spectrum is not the same. (see Example for details)

    Examples:
        todo

    See Also:
        fft_period, wave_linear, period_graph

"""
    ch = 'amph' if amp else "reim"
    return _spectrum.py_fft_fre(dt, x, ch=ch)


def fft_preiod(x, dt=1, modifcation=False):
    """ Get the period and amplitude of a vector by FFT method

    Args:
        x (array): x(n), Array (or time series) to be analyzed
        dt (float): Time interval of array x
        modification (boolean): if true, modified amplitude to more accuracy for big dt and small n*dt. For example: if dt=1.0 day, then when get annual and semiannual amplitude, no correction is needed, but when dt=1.0 month, you'd better add the modification (Generally 1% increase in amplitude). More details: Markus Bath, Spectral analysis in geophysics, Elsevier Scientific Publishing Company, New York, 1974.

    Returns:
        tuple
            pe (array): pe(n-1), period of time series x
            am (array): am(n-1), amplitude of time series x corresponding with period pe

    See Also:
        fft_fre, wave_linear, period_graph

"""
    ch = 'y' if modifcation else 'd'
    return _spectrum.py_fft_period(dt, x, ch=ch)


def fns(n, lag, dt=1):
    """ Amplitude of time series x corresponding with period pe

    Args:
        n(int): The number of original time series to get lag
        lag (float): The lag-1 autocorrelations of the original time series, see comments for details
        dt (float): Time interval of the original time series

    Returns:
        tuple
            * per (array): per(1:n/2), The periods of the original time series
            * pk (array): pk(1:n/2), Normalized Fourier Red / White noise spectrum

    Notes:
        If lag is not 0.0, then pk is the normalized Fourier red noise spectrum; if lag is 0.0, then pk==1.0 is the normalized Fourier white noise spectrum. To get lag from a time series, you can use subroutine ARCO (e.g call arco(n,2,x,co,cor), then lag=(cor(1)+sqrt(cor(2)))/2. Note that, to get lag, it's better to pre-whiten (e.g remove known period signals, such as seasonal terms) the original time series first and then to get lag.

    References:
        ..[1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16)

    See Also:
        fft_spec

"""
    return _spectrum.py_fns(n, dt, lag)


def fre_wavanum(c, s, dt=1.0):
    """Get frequency-wavenumber spectrum from time depend spherical harmonic coefficients

    Args:
        c (array) : c(n, n, nt) Coefficients of cosine term
        s (array) : s(n, n, nt) Coefficients of sine term
        dt (float): Time interval

    Returns:
        tuple
            * pe (nt/2):The periods in frequency domain (1/pe is the frequency)
            * fw (nt/2,0:n) : Frequency-wavenumber spectrum
            * f (nt/2): Frequency spectrum (summed on n of fw)
            * w (0:n): Wavenumber spectrum (summed on nt of fw)

    Notes:
        Here, the frequency-wavenumber spectrum means the frequency-dependent degree variance power spectrum. The results of fw(nt/2,0:n), f(nt/2) and w(0:n) is calculated from Eq. 8, 9 and 11 of Wunsch and Stammer, 1995 paper, respectively. If the unit of spatial data obtained from c,s coefficients is mm, then the unit of the above power spectrum is :math:`mm^2`.

    References:
        ..[1] Wunsch C., D. Stammer, The global frequency-wavenumber spectrum of oceanic variability estimated from TOPEX/POSEDION altimetric measurements, J. Geophys. Res., Vol. 100, No. C12, P24895-24910, 1995. (EQ. 8)

    Examples:
        todo

    See Also:
        todo
"""
    return _spectrum.py_fre_wavenum(dt, c, s)


def hilbert(x):
    """Get a vector's Hilbert transform

    Args:
        x (array) : Array (or vector) to do Hilbert transform

    Returns:
        array: Hilbert transform vector of array x

    Notes:
        Hilbert transforms are essential in understanding many modern modulation methods. These transforms effectively phase shift a function by 90 degrees independent of frequency. Of course practical implementations have limitations. For example, the phase shifting of a low frequency implies a long delay, which in turn implies a computational process that maintains a long history of the signal. Hilbert transforms are useful in creating signals with one sided Fourier transforms. Also the concepts of analytic functions and analytic signals will be shown to be related through Hilbert transforms.  (comments from the below reference of Turner)

        Time domain convolution:

            ..math::
                xhil(t) = x(t)/(\pi t) = \frac{1}{\pi}\int_{-\infty}^{\infty} x(u)/(t-u) du

        Frequency domain multiplication:

            ..math::
                xhil(\omega)=(-i\times \textrm{sgn}(\omega))) \times X(\omega)

            where :math:`i*i=-1`, :math:`\textrm{sgn}` is a signum function, and :math:`\textrm{sgn}(\omega)=1, 0, or -1` for :math:`\omega > 0,  \omega=0, \omega < 0`, respectively.

        Hilbert Transform Properties:
            see WFL

    Example:
        todo
"""
    return _spectrum.py_wfl_hilbert(x)


def hilbert_amp(x, dt=1.0):
    """Get a vector's instantaneous envelope, frequency and phase by Hilbert transform  method

    Args:
        x (array): Array (or vector) to do Hilbert transformed
        dt (float): Time interval of array (or vector) x

    Returns:
        tuple
            * amp (n): Instantaneous amplitude of array x by Hilbert transform
            * ph (n): Instantaneous phase of array x by Hilbert transform
            * fre (n): Instantaneous frequency of array x by Hilbert transform

    Notes:

        ..math::
            Amp & = & \sqrt(x_t^2 + xhil_t^2) \\\\
            fre & = & \arctan(xhil/x) \\\\
            ph  & = & \frac{1}{2\pi} \frac{fre_t}{dt}

        here :math:`xhil_t` is the hilbert transform of :math:`x_t`

"""
    return _spectrum.py_hilbert_amp(dt, x)


def period_graph(x, dt, mw, per_low, per_up):
    """Calculate the period graph (or amplitude spectrum) of a time series

    Args:
        x (array): Time series to be used to get period graph (or amplitude spectrum)
        dt (float): Time interval
        mw (int): The number of frequency (or period) to be analyzed
        per_low (float): The lowest analyzed periods
        per_up (float): The highest analyzed periods

    Returns:
        tuple
            * per (array): per(mw), the analyzed periods time series
            * sp (array): sp(mw), The  amplitude spectrum corresponding to per

    References:
        ..[1]  Ding Yuerong, Dawei Zheng, Data processing method for astronomic observed data, Nanjing Unvi. Press, Nanjing, 1990. (in Chinese)

    Examples:
        >> import numpy as np
        >> import matplotlib.pyplot as plt
        >> dt = 18.2621
        >> t = np.arange(400) * dt
        >> x = 0.14 * np.sin(2 * np.pi * t / 435) + 0.10 * np.sin(2 * np.pi * t / 365)
        >> per, sp = period_graph(x, dt, 50, 300, 500)
        >> plt.plot(per, sp)
"""
    return _spectrum.py_pg(dt, x, mw, per_low, per_up)


def wave_linear(y, dt, s0, step, jtot, param, sigma):
    """Calculate wavelet transform of a time series and get amplitude spectrum

    Args:
        y (array) : y(n), time series to be analyzed
        dt (float) : time interval of array
        s0 (float) : begin periods to analysis :math:`(s_0 \leq 2 dt)`
        step (float) : step length of period
        jtot (float) : total number to be analyzed, if the end periods is :math:`ss`, then :math:`jtot=(ss - s0)/ step`

    Returns:
        tuple
            * scale (array) : scale(jtot), fourier period analyzed, it is correspond with amp(n, jtot)
            * amp (array) : amp(n, jtot), amplitude at desired period of input vector y
            * coi (array) : coi(n), cone of influence at time :math:`1 \rightarrow n`

    Notes:
        The results of subroutine WAVE_LINEAR keep the same amplitude of sine/cosine wave at different frequencies but the white noise spectrum is ~1/scale.  See reference of Maraun and Kurths (2004) for details.
        Opinions from Maraun and Kurths (2004):
        In Fourier analysis, any normalization automatically provides the following two features:
            * The Gaussian white noise spectrum is (by definition) flat.
            * Sines of the same amplitude have the same integrated power in the frequency domain.
        In wavelet analysis, we meet difficulties to obtain both. For the factor c(s) in Eq. (1), Torrence and Compo (1998) and Kaiser (1994) suggest different normalizations, which only preserve one of the mentioned features. Torrence suggests c(s) = (Dt/ s)1/2 which preserves a flat white noise spectrum, but sines of equal amplitude exhibit different integrated power proportional to their oscillation scale. This choice equals Kaiser’s normalization with p = 1/2 (Eq. 3.5 in Kaiser, 1994, p. 62). Using the normalization of Kaiser with p = 0, c(s) = (Dt)1/2, sines of equal amplitude also exhibit equal integrated power (when scale is plotted logarithmically). On the other hand, the white noise spectrum is no longer flat but decreases proportional to 1/scale.  Fortunately, the normalization is only relevant for a first inspection by eye. A normalization following  c(s) = (Dt/ s)1/2 emphasizes power on high scales and could lead to misinterpretations. However, when performing a significance test, the significance level already includes the chosen normalization.

        Here, the Morlet wavelet is real part of general complex Morlet wavelet.

        1. Complex Morlet mother wavelet:

            ..math::
                \Psi(t)=\pi^{-1/4}e^{-t^2/2}e^{i \omega_0 t}

            where :math:`\omega_0` is the nondimensional frequency, here taken to be 6 to satisfy the admissibility condition.

        2. Real Morlet mother wavelet (used by WAVE_LINEAR):

            ..math::
                \Psi(t)=e^{-t^2/(2\delta^2)}\cos(\omega_0 t)

            where :math:`\omega_0` and :math:`\delta > 0`

    References:
        ..[1] Torrence, C. and Compo, G.: A practical guide to wavelet analysis, Bull. Amer. Meteor. Soc., 79, 61–78, 1998.
        ..[2] Liu L. T., Basic wavelet theory and its applications in geosciences, Ph. D. Dissertation, WHIGG, Wuhan, 1999 (in Chinese).
        ..[3] Maraun D., and J. Kurths, Cross wavelet analysis: significance testing and pitfalls, Nonlinear Processes in Geophysics, Vol. 11, P505-514, 2004.
"""
    return _spectrum.py_wave_linear(dt, s0, step, jtot, param, sigma, y)

import numpy as np
import matplotlib.pyplot as plt
dt = 18.2621
t = np.arange(400) * dt
x = 0.14 * np.sin(2 * np.pi * t / 435) + 0.10 * np.sin(2 * np.pi * t / 365)

per, sp = period_graph(x, dt, 50, 300, 500)
plt.plot(per, sp)
plt.show()
