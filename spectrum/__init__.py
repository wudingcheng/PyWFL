# encoding: utf-8
from __future__ import absolute_import
from . import _spectrum
import numpy as np


__all__ = ["arco", "emd", "emd_tol", "fft_spec", "fft_fre", "fft_preiod", "fns", "fre_wavanum", "hilbert", "hilbert_amp",
           "period_graph", "wave_linear", "conv", "fft_window", "hilbert_emd_fre", "hilbert_emd_per", "hilbert_spec_fre",
           "hilbert_spec_per", "ceof_ith_eigen", "fft_am_window", "fft_omiga", "fft_spec_window", "mtm", "mtm_overlap",
           "multicor_spec_mtm", "multicor_spectra", "svd_coupled", "wavelet", "wavelet_coherency", "wavelet_coherency_lag",
           "wavelet_signif"]


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


    References:
        .. [1] Ding Y., W. Zheng, Astronomical data processing, Nanjing Univ.  Press, Nanjing, 1990. (in Chinese)
        .. [2] 丁月蓉, 郑大伟 , 天文测量数据的处理方法, 南京大学出版社, 1990.

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
        std (float) : Standard deviation to stop the iteration (0.2<=std<=0.3), default is 0.3

    Returns:
        array : IMFs (first to m-1 column is IMFs, last one (m) is residual)

    Notes:
        Edge conditions: extended symmetrical 2 points of both end maximum values and minimum values respectively.
        EMD_TOL is more speed than EMD, but EMD may more accuracy (personal opinion)

    References:
        .. [1] N. E. Huang et al., "The empirical mode decomposition and the Hilbert spectrum for non-linear and non stationary time series analysis, " Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998

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
        .. [1] N. E. Huang et al., "The empirical mode decomposition and the Hilbert spectrum for non-linear and non stationary time series analysis, " Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998

    Examples:
        todo

"""
    return _spectrum.py_emd(t, x, m)


def fft_spec(x, dt=1, p=0.95, lag1_co=0):
    """Get normalized Fourier spectrum, red-noise background spectrum, and 100p% (e.g. p=95%, 99%) confidence spectrum of time series.

    Args:
        x (array) : Time series to be analyzed
        dt (float) : Time interval of x
        p (float) : Confidential level (0<p<1.0, typically p=0.95, p=0.99)
        lag1_co (float) : The lag-1 autocorrelations of x, used to output background red-noise spectrum, see FNS for details. If present, lag1_co as the lag-1 autocorrelations; if not present, use ARCO to get it, and then calculate FNS.

    Returns:
        array:　spec(:, 4) :　Normalized Fourier spectrum. spec(:, 1), period index; spec(:, 2), Fourier spectrum of x; spec(:, 3) background Red-noise spectrum; spec(:, 4), 100p% confidence spectrum

    Reference:
        .. [1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16), (17)

    Examples:
        todo
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
        :func:`fft_period`, :func:`wave_linear`, :func:`period_graph`

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
        :func:`fft_fre`, :func:`wave_linear`, :func:`period_graph`

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
        If lag is not 0.0, then pk is the normalized Fourier red noise spectrum; if lag is 0.0, then pk==1.0 is the normalized Fourier white noise spectrum. To get lag from a time series, you can use subroutine ARCO (e.g call arco(n, 2, x, co, cor), then lag=(cor(1)+sqrt(cor(2)))/2. Note that, to get lag, it's better to pre-whiten (e.g remove known period signals, such as seasonal terms) the original time series first and then to get lag.

    References:
        ..[1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16)

    See Also:
        :func:`fft_spec`

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
            * fw (nt/2, 0:n) : Frequency-wavenumber spectrum
            * f (nt/2): Frequency spectrum (summed on n of fw)
            * w (0:n): Wavenumber spectrum (summed on nt of fw)

    Notes:
        Here, the frequency-wavenumber spectrum means the frequency-dependent degree variance power spectrum. The results of fw(nt/2, 0:n), f(nt/2) and w(0:n) is calculated from Eq. 8, 9 and 11 of Wunsch and Stammer, 1995 paper, respectively. If the unit of spatial data obtained from c, s coefficients is mm, then the unit of the above power spectrum is :math:`mm^2`.

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

        .. math::
            Amp & = & \\sqrt(x_t^2 + xhil_t^2) \\\\
            fre & = & \\arctan(xhil/x) \\\\
            ph  & = & \\frac{1}{2\\pi} \frac{fre_t}{dt}

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
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> dt = 18.2621
        >>> t = np.arange(400) * dt
        >>> x = 0.14 * np.sin(2 * np.pi * t / 435) + 0.10 * np.sin(2 * np.pi * t / 365)
        >>> per, sp = period_graph(x, dt, 50, 300, 500)
        >>> plt.plot(per, sp)
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
            * coi (array) : coi(n), cone of influence at time :math:`1 \\rightarrow n`

    Notes:
        The results of subroutine WAVE_LINEAR keep the same amplitude of sine/cosine wave at different frequencies but the white noise spectrum is ~1/scale.  See reference of Maraun and Kurths (2004) for details.
        Opinions from Maraun and Kurths (2004):
        In Fourier analysis, any normalization automatically provides the following two features:
            * The Gaussian white noise spectrum is (by definition) flat.
            * Sines of the same amplitude have the same integrated power in the frequency domain.
        In wavelet analysis, we meet difficulties to obtain both. For the factor c(s) in Eq. (1), Torrence and Compo (1998) and Kaiser (1994) suggest different normalizations, which only preserve one of the mentioned features. Torrence suggests c(s) = (Dt/ s)1/2 which preserves a flat white noise spectrum, but sines of equal amplitude exhibit different integrated power proportional to their oscillation scale. This choice equals Kaiser’s normalization with p = 1/2 (Eq. 3.5 in Kaiser, 1994, p. 62). Using the normalization of Kaiser with p = 0, c(s) = (Dt)1/2, sines of equal amplitude also exhibit equal integrated power (when scale is plotted logarithmically). On the other hand, the white noise spectrum is no longer flat but decreases proportional to 1/scale.  Fortunately, the normalization is only relevant for a first inspection by eye. A normalization following  c(s) = (Dt/ s)1/2 emphasizes power on high scales and could lead to misinterpretations. However, when performing a significance test, the significance level already includes the chosen normalization.

        Here, the Morlet wavelet is real part of general complex Morlet wavelet.


        1. Complex Morlet mother wavelet:

            .. math::
                \Psi(t)=\pi^{-1/4}e^{-t^2/2}e^{i \omega_0 t}

            where :math:`\omega_0` is the nondimensional frequency, here taken to be 6 to satisfy the admissibility condition.

        2. Real Morlet mother wavelet (used by WAVE_LINEAR):

            .. math::
                \Psi(t)=e^{-t^2/(2\delta^2)}\cos(\omega_0 t)

            where :math:`\omega_0` and :math:`\delta > 0`

    References:
        .. [1] Torrence, C. and Compo, G.: A practical guide to wavelet analysis, Bull. Amer. Meteor. Soc., 79, 61–78, 1998.
        .. [2] Liu L. T., Basic wavelet theory and its applications in geosciences, Ph. D. Dissertation, WHIGG, Wuhan, 1999 (in Chinese).
        .. [3] Maraun D., and J. Kurths, Cross wavelet analysis: significance testing and pitfalls, Nonlinear Processes in Geophysics, Vol. 11, P505-514, 2004.
"""
    return _spectrum.py_wave_linear(dt, s0, step, jtot, param, sigma, y)


def conv(x, response, method='matlab', isconv=True):
    x = np.asarray(x)
    return _spectrum.py_wfl_conv(x, response, method, isconv)


def fft_window(width, window_type='square', window_parameter=True):
    """Data window used by Fourier spectrum analysis.

    Args:
        width : Length (or data number) of window
        window_type : Type of the data window. window_type is one of below: 'square', 'rectangular', 'bartlett', 'hamming', 'hanning', 'welch'.  Note: 'square'=='rectangular'.
        window_parameter : True for a whole window from index 1 to width, False for a half window from index 1 to width, respectively. See comments for details.

    Returns:
        tuple:
            * window_result (array): output window results
            * bandwidth (float): bandwidth of the window. See comments for details.

    Notes:
        The most common windows and their features are given below. This table can be used to choose the best windowing function for each application.

        +-------+-----------------------------+-----------------------+------------------+---------------------+
        |Window | Best for these Signal Types |  Frequency Resolution | Spectral Leakage |  Amplitude Accuracy |
        +=======+=============================+=======================+==================+=====================+
        |Barlett|  Random                     |Good  　　　　　　　　 |  Fair            |                Fair |
        +-------+-----------------------------+-----------------------+------------------+---------------------+
        |Hanning| Random                      |Good                   |   Good           |                Fair |
        +-------+-----------------------------+-----------------------+------------------+---------------------+
        |Hamming| Random                      |Good                   |  Fair            |                Fair |
        +-------+-----------------------------+-----------------------+------------------+---------------------+
        |Welch  | Random                      |Good                   |  Good            |                Fair |
        +-------+-----------------------------+-----------------------+------------------+---------------------+

        Whole data window equations (here m=window_length, window_parameter=='one'):

            * 'square' or 'rectangular' :math:`w(1:m) = 1.0`
            * 'bartlett': :math:`w(j)= 1 - \\left | \\frac{2j - 2 - m}{m} \\right |  \\quad (j=1, \\dots, m)`
            * 'hamming' : :math:`w(j) = 0.54 - 0.46 \\cos(2\\pi j /m) \\quad (j=1, \\dots, m)`
            * 'hanning' : :math:`w(j) = 0.5 - 0.5 \\cos(2\\pi j /m) \\quad (j=1, \\dots, m)`
            * 'welch'   : :math:`w(j) = 1 - \\left ( \\frac{m^2(j-1)}{4} \\right)^2

        Half data window equations (here m=window_length, window_parameter=='half'):

            * 'square' or 'rectangular':  :math:`w(1:m)=1.0`
            * 'bartlett': :math:`w(j)= 1 - \\left | j/m \\right |   \\quad (j=1, \\dots, m)`
            * 'hamming': :math:`w(j)=0.54+0.46 * \\cos ( \\pi*j/m)  \\quad (j=1, \\dots, m)`
            * 'hanning': :math:`w(j)=0.5+0.5 * \\cos( \\pi*j/m)     \\quad (j=1, \\dots, m)`
            * 'welch': :math:`w(j)=1- (j/m)^2                     \\quad (j=1, \\dots, m)`

        :math:`w(0)==1.0` for half data windows.

        .. math:: bandwidth = 1 / \int_{-\\infty}^{+\\infty} w(u)^2 du

        here, :math:`w(u)` is th window_result, and integral interval is :math:`\\left[ -m, m \\right]`. Note, here bandwidth is defined in :math:`\\left[ -m, m \\right]` but NOT in :math:`\\left[ 0, m \\right]`.

    Referemces:
        ..[1] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 6 and 7)
        ..[2] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, Numerical Recipes in Fortran 77, Second Edition, Syndicate of the University of Cambridge Press, p550, 1992.

    Examples:
        >>> fft_window(10)
        >>> fft_window(10, window='bartlett')
        >>> fft_window(10, window_type='welch')
        >>> fft_window(10, window_type='welch', window_parameter=False)

"""
    window_parameter = 'one' if window_parameter else 'half'
    conds = ['square', 'rectanglular', 'bartlett', 'hamming', 'hanning', 'welch']
    if window_type not in conds:
        raise ValueError("window_type should be in [{}]".format(' '.join(conds)))
    return _spectrum.py_wfl_window(width, window_type, window_parameter)


def hilbert_emd_fre(iunit, t, fre_low, fre_up, emd_imf, nfreline, filename, tlow=None, tup=None):
    """Produce a file which contains instantaneous amplitude and frequency of all EMD's IMFs by Hilbert method. The output file is a contour file which can be plotted in Golden Software Surfer.

    Args:
        iunit (int): File unit of output data
        t (array) : Time index of array (or vector) emd_imf
        fre_low (float) : Lower frequency to be analyzed, If fre_low<=2/(n*dt) then fre_low=2/(n*dt), here dt=t(2)-t(1)
        fre_up (float) : Uper frequency to be analyzed, If fre_up>=1/(2*dt) then fre_up=1/(2*dt), here dt=t(2)-t(1)
        emd_imf (array) : emd_imf[n, m], Input EMD's IMFs (n is number in time domain, m is m_th IMFs)
        nfreline (int) :  Line numbers of frequency between fre_low and fre_up to be analyzed
        filename (string) :  Output contour data filename (add ".dat.grd" to filename after output)


    Notes:
        :func:`hilbert_emd_fre` get an output file which can plot contour in surfer, the results is come form  :func:`hilbert_amp. We get the instantaneous envelope of EMD's IMFs by HILBERT method, then put all the data to a single Golden Software SURFER contour grid file. :func:`hilbert_emd_fre` output time, frequency, envelope for all IMFs.

        :func:`hilbert_emd_fre`, :func:`hilbert_emd_per`, :func:`hilbert_spec_fre` and :func:`hilbert_spec_per`  are similar. :func:`hilbert_emd_fre` and :func:`hilbert_emd_per` output more than one HILBERT amplitude in one single output file, but :func:`hilbert_spec_fre` and :func:`hilbert_spec_per` only output one instantaneous envelope in one single output file.

        :func:`hilbert_emd_fre` not only can be used to output instantaneous envelope of EMD's IMFs, but also can be used by any signal like EMD's IMFs. It means you can analyze two dimension signal x(n,m) where n is time number and m is m_th signal. Specially, :func:`hilbert_spec_fre` and :func:`hilbert_spec_per` is the special case of :func:`hilbert_emd_fre` and :func:`hilbert_emd_per`, in this special case, m==1.

        If the periods you interesting is not long, :func:`hilbert_emd_per` is more suitable than :func:`hilbert_emd_fre`, otherwise you'd better choose :func:`hilbert_emd_fre`.

"""
    if None in [tlow, tup]:
        tlow, tup = t[0], t[-1]
    _spectrum.py_hilbert_emd_fre(iunit, t, fre_low, fre_up, emd_imf, nfreline, filename, tlow=tlow, tup=tup)


def hilbert_emd_per(iunit, t, per_low, per_up, emd_imf, nfreline, filename, tlow=None, tup=None):
    """Produce a file which contains instantaneous amplitude and period of all EMD's IMFs by Hilbert method. The output file is a contour file which can be plotted in Golden Software Surfer.

    Args:
        iunit (int): File unit of output data
        t (array) : Time index of array (or vector) emd_imf
        fre_low (float) : Lower frequency to be analyzed, If fre_low<=2/(n*dt) then fre_low=2/(n*dt), here dt=t(2)-t(1)
        fre_up (float) : Uper frequency to be analyzed, If fre_up>=1/(2*dt) then fre_up=1/(2*dt), here dt=t(2)-t(1)
        emd_imf (array) : emd_imf[n, m], Input EMD's IMFs (n is number in time domain, m is m_th IMFs)
        nfreline (int) :  Line numbers of frequency between fre_low and fre_up to be analyzed
        filename (string) :  Output contour data filename (add ".dat.grd" to filename after output)

    See Also:
        :func:`hilbert_emd_fre`, :func:`hilbert`, :func:`hilbert_amp`, :func:`hilbert_emd_fre`, :func:`hilbert_spec_fre`, :func:`hilbert_spec_per`
"""
    _spectrum.py_hilbert_emd_per(iunit, t, per_low, per_up, emd_imf, nfreline, filename, tlow=tlow, tup=tup)


def hilbert_spec_fre(iunit, t, fre_low, fre_up, x, nfreline, filename):
    """Produce a file which contains instantaneous amplitude and frequency of a vector by Hilbert method. The output file is a contour file which can be plotted in Golden Software Surfer.

    Args:
        iunit (int): File unit of output data
        t (array) : Time index of array (or vector) emd_imf
        fre_low (float) : Lower frequency to be analyzed, If fre_low<=2/(n*dt) then fre_low=2/(n*dt), here dt=t(2)-t(1)
        fre_up (float) : Uper frequency to be analyzed, If fre_up>=1/(2*dt) then fre_up=1/(2*dt), here dt=t(2)-t(1)
        x (array) :  Input vector x
        nfreline (int) :  Line numbers of frequency between fre_low and fre_up to be analyzed
        filename (string) :  Output contour data filename (add ".dat.grd" to filename after output)

    See Also:
        :func:`hilbert_emd_fre`, :func:`hilbert`, :func:`hilbert_amp`, :func:`hilbert_emd_fre`, :func:`hilbert_emd_per`, :func:`hilbert_spec_per`
"""
    _spectrum.py_hilbert_spec_fre(iunit, t, fre_low, fre_up, x, nfreline, filename)


def hilbert_spec_per(iunit, t, per_low, per_up, x, nperline, filename):
    """Produce a file which contains instantaneous amplitude and frequency of a vector by Hilbert method. The output file is a contour file which can be plotted in Golden Software Surfer.

    Args:
        iunit (int): File unit of output data
        t (array) : Time index of array (or vector) emd_imf
        fre_low (float) : Lower frequency to be analyzed, If fre_low<=2/(n*dt) then fre_low=2/(n*dt), here dt=t(2)-t(1)
        fre_up (float) : Uper frequency to be analyzed, If fre_up>=1/(2*dt) then fre_up=1/(2*dt), here dt=t(2)-t(1)
        x (array) :  Input vector x
        nfreline (int) :  Line numbers of frequency between fre_low and fre_up to be analyzed
        filename (string) :  Output contour data filename (add ".dat.grd" to filename after output)

    See Also:
        :func:`hilbert_emd_fre`, :func:`hilbert`, :func:`hilbert_amp`, :func:`hilbert_emd_fre`, :func:`hilbert_emd_per`, :func:`hilbert_emd_per`
"""
    _spectrum.py_hilbert_spec_per(iunit, t, per_low, per_up, x, nperline, filename)


def ceof_ith_eigen(nt, eigenvalues):
    """Determine ith eigenvalus are significant to the error level (only used after all subroutine :func:`ceof`)

    Args:
        nt (int) : Number of original time series (reference nt in :func:`ceof`)
        eigenvalues (array) : Eigenvalues calculate in :func:`ceof`

    Return:
        int : The ith eigenvalues are significant to the error level. Note, if ith==nmodes-1, it has two possible explanations: 1) it's true that the nmode-1 is the last eigenvalues which are significant to the error level; 2) the nmodes is small here, you should use more larger nmodes to see if ith is correct. Because here only input nmodes eigenvalues, so we can get the largest ith is nmodes-1.

    References:
        .. [1] North G., T. Bell, R. Cahalan and F. J. Moeng, Sampling errors in the estimation of empirical orthogonal function, Vol.110, Mon. Wea. Rev., 699-706,1982.
"""
    return _spectrum.py_wfl_ceof_ith_eigen(nt, eigenvalues)


def fft_am_window(dt, x, m, window, overlap, p, lag1_co=-1):
    """Get windowed Fourier amplitude spectrum with suitable window, red-noise background amplitude spectrum, and 100p% (e.g. p=95%,99%) confidence amplitude spectrum of time series.

    Args:
        dt (float) : Time interval of x
        x (array) : x[n], Time series to be analyzed
        m (int) : Number data of window width (m<=n), then the x(n) will be seprate to k=n/m or 2k-1 segments to be analysis with FFT method. See comments for details.
        window_type (string) : Window type used. Can be one of the belows: 'square', 'rectangular', 'bartlett', 'hamming', 'hanning', 'welch'. See :func:`fft_window` for details.
        overlap (bool) :.If overlap, then do overlap analysis, otherwise do non-overlap analysis. See comments for details.
        p (float) :Confidential level (0<p<1.0, typically p=0.95, p=0.99)
        lag1_co (float) : The lag-1 autocorrelations of x, used to output background red-noise amplitude spectrum, see :func:`fns` for details. If presented, lag1_co as the lag-1 autocorrelations; if not present, use :func:`arco` to get it, and then calculate :func:`fns`. Note, if set lag1_co=0.0, the background is white-noise amplitude spectrum.

    Returns:
        tuple:
            * am_window(m/2,2): Windowed Fourier amplitude spectrum of input time series. am_window(:,1), period index; am_window(:,2), windowed Fourier amplitude spectrum. Note: am_window(:,1) is different with am(:,1) output when m/=n.
            * am(n/2,3): Windowed Fourier amplitude spectrum of background (noise). am(:,1), period index; am(:,2) background Red-noise amplitude spectrum; am(:,3), 100p% confidence amplitude spectrum. Note that the period index of am(:,1) and am_window(:,1) is different, so you must plot am_window(:,2) with am_window(:,1) and am(:,2:3) with am(:,1) to get correct plot of windowed amplitude spectrum and corresponding noise amplitude spectrum and amplitude spectrum at confidential level p.

    Notes:
        Here, k=n/m windows are used to filter the original x(n), which can decrease the power leakage in analysis frequency. m is the data number in each segments. If overlap==True, then the original data will be analyzed at 2k-1 segments. Each segments have half window width overlap. E.g., (0,m), (m/2+1,3m/2), (m+1,2m), ... will be analyzed. If overlap==.false., the the original data will be analyzed at k segments, no overlap in the analyzed data. E.g., (0,m), (m+1,2m), (2m+1,3m), ...

        comments from "Numerical Recipes(NR) in Fortran77", P549, modified by WFL Author.

        We might want to obtain the smallest variance from a fixed amount of computation, without regard to the number of data points used. This will generally be the goal when the data are being gathered in real time, with the data-reduction being computer-limited. Alternatively, we might want to obtain the smallest variance from a fixed number of available sampled data points. This will generally be the goal in cases where the data are already recorded and we are analyzing it after the fact.

        In the first situation (smallest spectral variance per computer operation), it is best to segment the data without any overlapping. The first M data points constitute segment number 1; the next M data points constitute segment number 2; and so on, up to segment number K, for a total of KM sampled points. The variance in this case, relative to a single segment, is reduced by a factor K.

        In the second situation (smallest spectral variance per data point), it turns out to be optimal, or very nearly optimal, to overlap the segments by one half of their length. The first M points are segment number 1; the m/2+1 and 3m/2 sets of M points are segment number 2; and so on, up to segment number 2K-1. The reduction in the variance is NOT a factor of K, since the segments are not statistically independent. It can be shown that the variance is instead reduced by a factor of about 9K/11. This is, however, significantly better than the reduction of about K=2 that would have resulted if the same number of data points were segmented without overlapping.

    References:
        .. [1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16), (17)
        .. [2] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 6 and 7)
        .. [3] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, Numerical Recipes in Fortran 77, Second Edition, Syndicate of the University of Cambridge Press, p550, 1992.
        .. [4] Harris, Fredric J. "On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform" in Proceedings of the IEEE Vol. 66, No. 1, January 1978.

"""
    return _spectrum.py_wfl_fft_am_window(dt, x, m, window, overlap, p, lag1_co=lag1_co)


def fft_omiga(n, dt):
    """Get FFT frequency value of a time series

    Args:
        n (int) : number of the time series to be do FFT
        dt (float) : time interval of the time series to be do FFT

    Returns:
        array: omiga[n], corresponding frequency of the time series used for FFT results, unit: 1/(unit of dt)

    Notes:
        omiga(1)==0.0, omiga(i)=(i-1)/((N-1)*dt)  for i=1,N/2+1, while omiga(N/2+1) is the Nyquist positive frequency. The negative frequency is corresponding with the positive frequency but with opposite sign. If N is even, then only one positive Nyquist frequency, if N is odd, then i=N/2+2 is the negative Nyquist frequency.
"""
    return _spectrum.py_wfl_fft_omiga(n, dt)


def fft_spec_window(dt, x, m, window, overlap, p, lag1_co=-1):
    """Get normalized Fourier spectrum with suitable window, red-noise background spectrum, and 100p% (e.g. p=95%,99%) confidence spectrum of time series.

    Args:
        dt (float) : Time interval of x
        x (array) : x[n], Time series to be analyzed
        m (int) : Number data of window width (m<=n), then the x(n) will be seprate to k=n/m or 2k-1 segments to be analysis with FFT method. See comments for details.
        window_type (string) : Window type used. Can be one of the belows: 'square', 'rectangular', 'bartlett', 'hamming', 'hanning', 'welch'. See :func:`fft_window` for details.
        overlap (bool) :.If overlap, then do overlap analysis, otherwise do non-overlap analysis. See comments for details.
        p (float) :Confidential level (0<p<1.0, typically p=0.95, p=0.99)
        lag1_co (float) : The lag-1 autocorrelations of x, used to output background red-noise amplitude spectrum, see :func:`fns` for details. If presented, lag1_co as the lag-1 autocorrelations; if not present, use :func:`arco` to get it, and then calculate :func:`fns`. Note, if set lag1_co=0.0, the background is white-noise amplitude spectrum.

    Returns:
        tuple:
            * spec_window(m/2,2):Normalized windowed Fourier spectrum. spec_window(:,1), period index; spec_window(:,2), Fourier spectrum. Note: spec_window(:,1) is different with spec(:,1) output when m/=n.
            * spec(n/2,3): Normalized Fourier spectrum. spec(:,1), period index; spec(:,2) background Red-noise spectrum; spec(:,3), 100p% confidence spectrum. Note that the period index of spec(:,1) and spec_window(:,1) is different, so you must plot spec_window(:,2) with spec_window(:,1) and spec(:,2:3) with spec(:,1) to get correct plot of windowed spectrum and corresponding noise spectrum and spectrum at confidential level p.

    Notes:
        Here, k=n/m windows are used to filter the original x(n), which can decrease the power leakage in analysis frequency. m is the data number in each segments. If overlap==True, then the original data will be analyzed at 2k-1 segments. Each segments have half window width overlap. E.g., (0,m), (m/2+1,3m/2), (m+1,2m), ... will be analyzed. If overlap==.false., the the original data will be analyzed at k segments, no overlap in the analyzed data. E.g., (0,m), (m+1,2m), (2m+1,3m), ...

        comments from "Numerical Recipes(NR) in Fortran77", P549, modified by WFL Author.

        We might want to obtain the smallest variance from a fixed amount of computation, without regard to the number of data points used. This will generally be the goal when the data are being gathered in real time, with the data-reduction being computer-limited. Alternatively, we might want to obtain the smallest variance from a fixed number of available sampled data points. This will generally be the goal in cases where the data are already recorded and we are analyzing it after the fact.

        In the first situation (smallest spectral variance per computer operation), it is best to segment the data without any overlapping. The first M data points constitute segment number 1; the next M data points constitute segment number 2; and so on, up to segment number K, for a total of KM sampled points. The variance in this case, relative to a single segment, is reduced by a factor K.

        In the second situation (smallest spectral variance per data point), it turns out to be optimal, or very nearly optimal, to overlap the segments by one half of their length. The first M points are segment number 1; the m/2+1 and 3m/2 sets of M points are segment number 2; and so on, up to segment number 2K-1. The reduction in the variance is NOT a factor of K, since the segments are not statistically independent. It can be shown that the variance is instead reduced by a factor of about 9K/11. This is, however, significantly better than the reduction of about K=2 that would have resulted if the same number of data points were segmented without overlapping.

        here the result of :func:`fft_spec_window` is NOT equivalent with the result of "spctrm" in NR.

    References:
        .. [1] Torrence C., G. P. Combo, A Practical Guide to Wavelet Analysis, Bulletin of the American Meteorological Society, Vol. 79, p61-78, 1998.  E.Q. (16), (17)
        .. [2] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 6 and 7)
        .. [3] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, Numerical Recipes in Fortran 77, Second Edition, Syndicate of the University of Cambridge Press, p550, 1992.
"""
    return _spectrum.py_wfl_fft_spec_window(dt, x, m, window, overlap, p, lag1_co=lag1_co)


def mtm(dt, x, nw, spec_type='unity', p=0.95, lag1_co=-1):
    """Get normalized Fourier spectrum by MultiTaper Method (MTM), and red-noise background spectrum, and 100p% (e.g. p=95%,99%) confidence spectrum of time series.

    Args:
        dt (float) : Time interval or sample interval of x
        x (array) : Time series to be analyzed  (real, double, complex, double complex)
        nw (float) : :func:`dpss` orders 0 to 2*nw-1 required (2<=nw<n/2). Note, here, nw is real or double type. The typical values of nw is 2, 5/2, 3, 7/2, 4. Thus you can get (nw)Pi-multitapers, e.g. for nw=4, you get 4Pi-multitapers. Smaller nw higher frequency precision, greater nw smaller variance in frequency domain.
        spec_type (string): weight method for multitaper window. spec_type is one of below (case-insensitive):
            'unity' weighted by: 1/(2*nw-1)* sum s(f)*s(f) k=0,2*nw-1 for sum
            'eigen' weighted by: 1./sum(eigen(:)) * sum s(f)*s(f)*eigen(k) k=0,2*nw-1
            'hires' weighted by: sum s(f)*s(f)/eigen(k) k=0,2*nw-1
            'adapt_pw' weighted by: 1./ sum ( bk/bk/eigen(k)) *sum bk*bk*eigen(k)*s(f)*s(f) k=0,2*nw-1 !see Percival and Walden (1993) eq. 368a,370a
            'adapt_th' (default) weighted by: 1./dk/dk * sum dk*dk*s(f)*s(f) !see Thomson (1982, eq. 5.2, 5.3)
            comment on spec_type: the first 3 weighted methods are appropriate for pure white noise,  and the last two are appropriate for colored noise.
        p (float) : Confidential level (0<p<1.0, typically p=0.95, p=0.99). Note, p and red_noise_spec must be present at same time if you want to get the red_noise_spec, also you can present lag1_co or not when to calculate red_noise_spec.
        lag1_co (float) : The lag-1 autocorrelations of x, used to output background red-noise amplitude spectrum, see :func:`fns` for details. If presented, lag1_co as the lag-1 autocorrelations; if not present, use :func:`arco` to get it, and then calculate :func:`fns`. Note, if set lag1_co=0.0, the background is white-noise amplitude spectrum.

    Returns:
        tuple:
            * f_test(n/2) (array):F-test which compares the variance of an assumed locally white noise background to the variance explained by a periodic signal. See Lees and Park(1994), eq. (6) and (7).
            * x_edof(n/2) (array):Effective degree of freedom in frequency domain.
            * red_noise_spec(n/2,3) (array):Normalized red noise Fourier spectrum. red_noise_spec(:,1), period index; red_noise_spec(:,2), AR(1) red noise Fouried spectrum; red_noise_spec(:,3)100p% confidence spectrum.

"""
    conds = ['unity', 'eigen', 'hires', 'adapt_pw', 'adapt_th']
    if spec_type not in conds:
        raise ValueError('spec_type should be in [{}]'.format(' '.join(conds)))
    return _spectrum.py_wfl_mtm(dt, x, nw, spec_type=spec_type, p=p, lag1_co=lag1_co)


def mtm_overlap(dt, x, m, overlap, nw, spec_type='unity', p=0.95, lag1_co=-1):
    return _spectrum.py_wfl_mtm_overlap(dt, x, m, overlap, nw, spec_type='unity', p=0.95, lag1_co=-1)


def multicor_spec_mtm(dt, x, nw):
    return _spectrum.py_wfl_multicor_spec_mtm(dt, x, nw)


def multicor_spectra(dt, x, max_lag, window_type, out_type, p=0.95):
    """Get windowed FFT spectrum of auto and cross covariance/correlation of input time series.

    Args:
        dt (float) : Time interval of time series of x
        x(n,m) (array) :Time series to be analyzed
        max_lag (int) :Number data of window width (max_lag<=n), w(1)=1.0, w(max_lag)=0.0.
        window_type (string) :Window type used. Can be one of the belows: 'square', 'rectangular', 'bartlett', 'hamming', 'hanning', 'welch'. See :func:`fft_window` for details.
        out_type (string) : Output control. out_type='COV' or 'COR'.
        p(float) :Confidential level (0<p<1.0, typically p=0.95, p=0.99). p, cll, clu must be presented at the same time.

    Returns:
        tuple:
            * pe(n/2) (array): eriods of analyzed time series x(n,m)
            * x_out(n/2,m,m) (array): uto and cross covariance/correlation FFT spectrum of input time series.
                       Auto- and cross covariances FFT spectrum for X if out_type='COV'
                       Auto- and corss correlation FFT spectrum for X, if out_type='COR'
                       The first index is periods( or 1/frequency). Thus x_out(:,i,j) is the covariance/correlation FFT spectrum of
                       the ith and jth series. And x_out(:,i,i) is the autocovariance/autocorrelation FFT spectrum of the ith series.
            * x_out2(n/2,m,m) (array): Phase spectrum. x_out2(:,i,i)==0.0, and x_out2(:,i,j) from -180 to 180 degree. Note, the phase may not continued with the frequency increase, so you should remove the jump by yourself.
            * x_out3(n/2,m,m) (array): Squared coherency spectrum. x_out3(:,i,i)=1.0, x_out3(:,i,j) from 0 to 1.0.  Please also reference :func:`corfre`.
            * cll(n/2,m,m) (array): Lower confidence limits for the auto-spectrum (x_out(:,i,i)) and for the squared coherency spectrum (xout3(:,i,j))
            * clu(n/2, m,m) (array): Upper confidence limits for the auto-spectrum (x_out(:,i,i)) and for the squared coherency spectrum (x_out3(:,i,j))

    Notes:
        If window_type='square' or 'rectangular', max_lag=n-1, out_type='cor', then the result of WFL_MULTICOR_SPECTRA is the same with :func:`fft_spec`. See example 1 for details. If out_type='cov', then var(x) must be divide to get the same result with out_type='cor'.

        The lower and upper confidence limits produce a confidence interval for the auto-spectrum or squared coherency spectrum. The interpretation of this e.g. 95% confidence interval is that, having constructed a random interval which has probability 0.95 of enclosing the true spectrum, we have 95% confidence on this particular occasion that the interval obtained happens to include auto-spectrum or squared coherency spectrum. Thus, the 99% confidence interval > 95% confidence interval. This concept is NOT the same as significance level (or confidence level) such as used in :func:`corfre`.  　

        Difference between :func:`corfre` and :func:`multicor_spectra`:
        use different smoothing method to get squared coherency spectrum. The definition of the squared coherency spectrum is equivalent.

    References:
        .. [1] Jenkins, G. M., and D. G., Watts, Spectral analysis and its applications, Holden-Day, 1968. (chapter 6-9)


"""
    return _spectrum.py_wfl_multicor_spectra(dt, x, max_lag, window_type, out_type, p=p)


def svd_coupled(s, p, kmode):
    """Do SVD decompsion for two coupled given time-space data fields.

    Args:
        S(n,mp) (array): One data set used to do SVD, must be anomalies (you can use RMAVER_NXNT or RMAVER_NTNX to get the anomalies)
        P(n,mq) (array): The other data set used to do SVD, must be anomalies (you can use RMAVER_NXNT or RMAVER_NTNX to get the anomalies)
        kmode (int) : The number of modes to be output. Note: kmode<=MIN(mp,mq)

    Returns:
        tuple
            * lk (kmode)(array) : The computed singular values of covariance matrix (transport(S)* P, or in IMSL operators: S .tx. P), the largest returned first
            * uk(mp,kmode)(array) : The computed singular vectors of array S, the most important returned first  (or spatial pattern for S)
            * vk(mp,kmode)(array) : The computed singular vectors of array P, the most important returned first  (or spatial pattern for P)
            * ak(n,kmode)(array) : The expansion coefficients for array S  (or spatial pattern evolves in time)
            * bk(n,kmode)(array) : The expansion coefficients for array P  (or spatial pattern evolves in time)
            * scfk(kmode)(array) : The percent of Square Covariance Fraction (SCF) explained by each singular value mode
            * cumk(kmode)(array) : The cumulative percent of SCF explained by each singular mode. This is a bit redundant -- it's just the sum of "scfk"
            * sk(n, mp)(array) : The reconstruction of field S use first kmode. Note: must be presented with pk at the same time.
            * pk(n, mp)(array) : The reconstruction of field P use first kmode. Note: must be presented with sk at the same time.

    Notes:
        For the real geophysical data, to get a more reliable results, the two fields such as S and P above should be weighted by square root area weight of sqrt(cosine(latitude)) first.
        See also sqrootweight(nx) parameters in subroutine WFL_CEOF.

        If there is an error like: "WFL ERROR: no convergence in wfl_svd (svdcmp)",please use double precision subroutine of WFL_SVD_COUPLED8.

        If there is an error with "stack overflow", please see here for how to resolve it.

    References:
        .. [1] Mewman M. and P. D. Sardeshmukh, A caveat concerning singular value decomposition, J Climate, Vol. 8, p352-360,1995. (very important)
        .. [2] Bretherton C. S., C. Smith, and J. M. Wallace, An intercomparison of methods for finding coupled patterns in climate data, J Climate, Vol. 5, P541-560,1992.
        .. [3] Björnsson H. and S. A. Venegas, A manual for EOF and SVD analyses of climate data, C2GCR Report No, 97-1, 1997.
        .. [4] Hannachi A., I. T. Jolliffe and D. B. Stephenson, Empirical orthogonal function s and related techniques in atmospheric science: A review, Int. J. Climatol., Vol. 27, p1119-1152, 2007.


"""
    return _spectrum.py_wfl_svd_coupled(s, p, kmode)


def wavelet(dt, s0, dj, jtot, y, param=6.0):
    """Calculate wavelet power spectrum of a time series using complex Molet wavelet (TC98 method)

    Args:
        dt() : Time interval of array y
        s0() : The smallest scale of the wavelet. Typically s0 = 2*dt. Note: for accurate reconstruction computation set s0=dt for Morlet
        dj() : The spacing between discrete scales. Typically dj = 0.25.  A smaller # will give better scale resolution, but be slower.
        jtot() : the # of scales.  Scales range from s0 up to s0*2**[(jtot-1)*dj]. Typically jtot=1+(LOG2(n*dt/s0))/dj
        y (n)() : Time series to be analyzed
        param (float) : Morlet wavelet parameters (default param=6.0). You can change the value of param if you familiar with parameters chosen, otherwise use default.

    Returns
        tuple
            * wave (n, jtot): (array) :complex (real and image part) of Morlet wavelet transform for input vector y versus time and scale.  CABS(wave) gives the wavelet amplitude, ATAN2(AIMAG(wave),REAL(wave)) gives wavelet phase. The wavelet power spectrum is CABS(wave)**2.
            * period(jtot): (array) :the "Fourier" periods (in time units) corresponding to "scale". period corresponding with n index of wave(n,jtot).
            * scale (jtot): (array) :the wavelet scales that were used, it has small difference with period.
            * coi (n): (array) : Cone of influence at time 1->n. the e-folding factor used for the cone of influence.

    Notes:
        The results of subroutine :func:`wavelet` do NOT keep the same amplitude of sine/cosine wave at different frequencies but keep the white noise spectrum flat.  See reference of Maraun and Kurths (2004) for details.

        Opinions from Maraun and Kurths (2004):

        In Fourier analysis, any normalization automatically provides the following two features:

        – The Gaussian white noise spectrum is (by definition) flat.
        – Sines of the same amplitude have the same integrated power in the frequency domain.

        In wavelet analysis, we meet difficulties to obtain both. For the factor c(s) in Eq. (1), Torrence and Compo (1998) and Kaiser (1994) suggest different normalizations, which only preserve one of the mentioned features. Torrence suggests c(s) = (Dt/ s)1/2 which preserves a flat white noise spectrum, but sines of equal amplitude exhibit different integrated power proportional to their oscillation scale. This choice equals Kaiser’s normalization with p = 1/2 (Eq. 3.5 in Kaiser, 1994, p. 62). Using the normalization of Kaiser with p = 0, c(s) = (Dt)1/2, sines of equal amplitude also exhibit equal integrated power (when scale is plotted logarithmically). On the other hand, the white noise spectrum is no longer flat but decreases proportional to 1/scale.  Fortunately, the normalization is only relevant for a first inspection by eye. A normalization following  c(s) = (Dt/ s)1/2 emphasizes power on high scales and could lead to misinterpretations. However, when performing a significance test, the significance level already includes the chosen normalization.

        Here, the Morlet wavelet is complex mother wavelet.

        Complex Morlet mother wavelet:

            .. math:: \\psi(t) = \\pi^{-1/4}\\exp(-t^2/2)\\exp(i \\omega_0 t)

         where :math:`\\omega_0` is the nondimensional frequency, here taken to be 6 to satisfy the admissibility condition.

    References:
        .. [1] Torrence, C. and Compo, G.: A practical guide to wavelet analysis, Bull. Amer. Meteor. Soc., 79, 61–78, 1998.
        .. [2] Maraun D., and J. Kurths, Cross wavelet analysis: significance testing and pitfalls, Nonlinear Processes in Geophysics, Vol. 11, P505-514, 2004.
"""
    return _spectrum.py_wfl_wavelet(dt, s0, dj, jtot, y, param=param)


def wavelet_coherency(time, scale, wave1, wave2, smooth=True, m=8, w=5, tolerance=1e-9):
    """Calculate wavelet squared coherency of two time series using complex Molet wavelet

    Args:
        time(nt) (array):  time index of wavelet transform
        scale(nj) (array):  scale index of wavelet spectrum (accquired from WFL_WAVELET subroutine)
        wave1(nt,nj) (array):  complex (real and image part ) of wavelet transform for time series #1 (accquired from WFL_WAVELET subroutine)
        wave2(nt,nj)  (array):  complex (real and image part ) of wavelet transform for time series #2 (accquired from WFL_WAVELET subroutine)
        smooth (bool):  Logical. Default is True If smooth=True, then smooth the cross wavelet spectrum of wave1 and wave2 in both time and scale direction by normalized Gaussian window and box car window, respectively.
        m (float):   Time direction smooth parameter: the ratio between window length and scale. Default is m=8.w ():  Scale direction smooth parameter, box car smoothing window length (must be odd number, if w is even number, then w=w+1). Default is w=5. if w<3, then w==3.
        tolerance (float):  tolerance for power1*power2, where power1=cabs(wave1)**2 and power2=cabs(wave2)**2. Default is 1.e-9. If power1*power2<tolerance, then wave_coher=0.0.

    Returns:
        tuple:
            * wave_coher (nt, nj) (array): squared wavelet coherency (0.0-1.0)
            * wave_phase(nt,nj) (array): the phase of corresponding squared  wavelet coherency (unit: degree, 0-360)
            * global_coher(nj) (array): the global (mean) squared coherency averaged over all times (0.0-1.0).
            * global_phase(nj) (array): the global (mean) phase averaged over all times (unit: degree, 0-360 ).
            * power1(nt,nj) (array): the wavelet power spectrum of time series #1, corresponding wave1. If smooth=True, then power1 is * smoothed power spectrum.
            * power2(nt,nj) (array): same as power1, but for time series #2.

    Notes:
        Before call WFL_WAVELET_COHERENCY, you must first call WFL_WAVELET twice. see example for details.

        The cross-wavelet spectrum is defined as:

        Wave_Cross(nt,s)=wave1(nt,s)*conj(wave2(nt,s))

        Then, the squared wavelet coherency is defined as::

                                      | < Wave_Cross(nt,ns) > | 2
            R2(nt,s)= -------------------------------------------------
                             < |wave1(nt,s)|2 > < |wave2(nt,s)|2 >

        where, < · > indecates smoothing in both time and scale direction, and |  ·  | indicates the complex module.

        Normalized Gaussian window (to have a total weight of unity) in time direction is : exp( -t2/(2s2)
        Normalized box cars window are used in scale direction.



    References:
        .. [1] Torrence, C. and Webster, P.: Interdecadal Changes in the ENSO Monsoon System, Journal of Climate, 12, 2679–2690, 1999
        .. [2] Torrence, C. and Compo, G.: A practical guide to wavelet analysis, Bull. Amer. Meteor. Soc., 79, 61–78, 1998.
        .. [3] Maraun D., and J. Kurths, Cross wavelet analysis: significance testing and pitfalls, Nonlinear Processes in Geophysics, Vol. 11, P505-514, 2004.
"""
    return _spectrum.py_wfl_wavelet_coherency(time, scale, wave1, wave2, smooth=smooth, m=m, w=w, tolerance=tolerance)


def wavelet_coherency_lag(time, scale, wave1, wave2, smooth=True, m=8, w=5, tolerance=1e-9):
    """Calculate wavelet squared coherency of two time series using complex Molet wavelet

    Args:
        time(nt) (array):  time index of wavelet transform
        scale(nj) (array):  scale index of wavelet spectrum (accquired from WFL_WAVELET subroutine)
        wave1(nt,nj) (array):  complex (real and image part ) of wavelet transform for time series #1 (accquired from WFL_WAVELET subroutine)
        wave2(nt,nj)  (array):  complex (real and image part ) of wavelet transform for time series #2 (accquired from WFL_WAVELET subroutine)
        smooth (bool):  Logical. Default is True If smooth=True, then smooth the cross wavelet spectrum of wave1 and wave2 in both time and scale direction by normalized Gaussian window and box car window, respectively.
        m (float):   Time direction smooth parameter: the ratio between window length and scale. Default is m=8.w ():  Scale direction smooth parameter, box car smoothing window length (must be odd number, if w is even number, then w=w+1). Default is w=5. if w<3, then w==3.
        tolerance (float):  tolerance for power1*power2, where power1=cabs(wave1)**2 and power2=cabs(wave2)**2. Default is 1.e-9. If power1*power2<tolerance, then wave_coher=0.0.

    Returns:
        tuple:
            * wave_coher_lag (-nt/2:nt/2, nj) (array) :lag squared wavelet coherency (0.0-1.0)
            * wave_phase_lag (-nt/2:nt/2,nj) (array) :the phase of corresponding lag squared  wavelet coherency (unit: degree, 0-360)
            * time_idx_lag (-nt/2:nt/2) (array) :lag time index

    Notes:
        Before call :func:`wavelet_coherency_lag`, you must first call WFL_WAVELET twice. see example for details.

        :func:`wavelet_coherency_lag` call :func:`wavelet_coherency` to calculate the global squared coherency with different lags.

        The default for smooth is .false. here, which is different with the default smooth=True in subroutine :func:`wavelet_coherency`

        The cross-wavelet spectrum is defined as:

        Wave_Cross(nt,s)=wave1(nt,s)*conj(wave2(nt,s))

        Then, the squared wavelet coherency is defined as::

                                      | < Wave_Cross(nt,ns) > | 2
            R2(nt,s)= -------------------------------------------------
                             < |wave1(nt,s)|2 > < |wave2(nt,s)|2 >

        where, < · > indecates smoothing in both time and scale direction, and |  ·  | indicates the complex module.



        Normalized Gaussian window (to have a total weight of unity) in time direction is : exp( -t2/(2s2)

        Normalized box cars window are used in scale direction.

        The lag squared coherency is defined as:

        wave_coher_lag(lag, nj)= do :func:`wavelet_coherency` (result of global_coher) using wave1(1:nt-lag,:) and wave2(lag+1:nt,:) for lag=[0,nt/2]

        wave_coher_lag(lag, nj)= do :func:`wavelet_coherency` (result of global_coher) using wave1(lag+1:nt,:) and wave2(1:nt-lag,:) for lag=[-nt/2,-1]


"""
    return _spectrum.py_wfl_wavelet_coherency_lag(time, scale, wave1, wave2, smooth=smooth, m=m, w=w, tolerance=tolerance)


def wavelet_signif(dt, y, s0, dj, scale, period, dof, param=6.0, siglvl=0.05, lag1=-1, isigtest=0):
    return _spectrum.py_wfl_wavelet_signif(dt, y, s0, dj, scale, period, dof, param=param, siglvl=siglvl, lag1=lag1, isigtest=isigtest)
