# encoding: utf-8
"""Filter functions of PyWFL"""
from _filters import filters as _filters

# __all__ = ['fft_filter', 'vondrak', 'gauss_jekeli']


def fft_filter(x, dt, per1, per2, method, window=0.0):
    """
    Args:
        x (ndarray) :Array(or time series) to be analyzed
        dt (float) :Time interval of array x
        per1 (float) :Period for filter
        per2 (float) :Period for filter
        method (string) :Filter method, should be **low**, **high**, **band**, means low pass filter, high pass filter and band pass filter. 
            If method is **low**, all signal period greater than per1 is retained, per2 not used;
            If method is **high**, all signal period smaller than per1 is retained, per2 not used;
            If method is **band**, all signal period from per1 to per2  is retained, per1 < per2;

        window (float) :Window width for filter, :math:`window*dt=T/5` window width (a cosine weighted window), here :math:`T` is window width, such as :math:`T=100days`, if :math:`dt=10days`, then you can set window as 2.0, then 
            filter window:

            .. math::
                w_1(t) & = & (1+\cos(5*\pi*t/T))/2 & \quad &    -T < t < -4T/5 \\\\
                w_2(t) & = & 1 & \quad &                     -4T/5 < t < 4T/5  \\\\
                w_3(t) & = & (1+\cos(5*\pi*t/T))/2 & \quad & 4T/5  < t < T 



    Returns:
        array_alike: Filtered series


    Notes:
        Detials for choose suitable window

    Examples:
        ...

    See Also:
        vondrak
    """
    return _filters.py_fft_filter(x, dt, per1, per2, method, window)


def vondrak(e, t, x):
    """
    Args:
        e (float) : Filter factor
        t (array_alike) : Time index of array x
        x (array_alike) : Time series to be filtered

    Returns:
        tuple : filtered time series and standard error
            * u (list) : Filtered time series
            * em (float) : Standard error between x and u.


    Notes:
        comments here

        For 3 order Vondrak filter, define frequency response function as,

        .. math:: F(e,T)=(1+1/e (2\pi /T)^6)^{-1}=A \quad (0<A<1)
            :label: eq1

        and we can conclude that,

        .. math:: e = (2 \pi /T)^6(A/(1-A))
            :label: eq2

        .. math:: T =2 \pi ((1-A)e/A)^{-1/6}
            :label: eq3

        Then for low pass filter, set :math:`A` close to 1.0( such as :math:`A=0.99` or :math:`A=0.95`) and :math:`T` as the truncate period (which means the the signal of period of :math:`T` can be retained as :math:`(A*100)%`), and we can calculate the corresponding filter factor e by using :eq:`eq2`. By using this :math:`e` in subroutine VONDRAK, we can get the signals for all period greater than T.
        For high pass filter, first do the low pass filter as above, and then use the original data minus the low pass filter data, you can get the high pass filter data.
        For band pass filter, you can do two low pass filter :math:`e1` and :math:`e2` with the different truncate period :math:`T1` and :math:`T2`, then the difference for the two filtered data is the results of band pass filter (see Example below).
        To get more information of parameter e in subroutine VONDRAK, please reference:

    References:
        .. [1] Yan Haoming, Zhong Min, Zhu Yaozhong, Determination of the degree of freedom of digital filtered time series with an application to the correlation analysis between the length of day and the southern oscillation index, Chinese Astronomy and Astrophysics, Vol. 28, No.1, P120-126, 2004.(in English)
        .. [2] 闫昊明，钟敏，朱耀仲，时间序列数字滤波后自由度的确定――应用于日长变化与南方涛动指数的相关分析，天文学报，Vol. 44, No.3, P324-329, 2003. (in Chinese)

    Examples:
        todo

    See Also:
        fft_filter : fft filter function
    """
    return _filters.py_vondrak(e, t, x)


def gauss_jekeli(r, distance, ch='distance'):
    """
    Args:
        r (float): r is the distance on the Earth's surface at which weight function w has dropped to 1/2 it value at :math:`\\alpha=0` (the distance on the Earth's surface = Earth_Radius*alpha). Or we can refer to r as the averaging radius (unit: meter).

        distance (float): Distance from the center point for using Gaussian filtering weight function (unit: meter)

        ch (string): default, 'distance'; optional in ['distance', 'degree']


    Returns:
        float
            If ch=='distance', the unit of distance_from_center should be in meter;
            If ch=='degree', the unit of distance_from_center should be in degree (0°-90° );

    Notes:
        The Gaussian filtering weight function (normalized here so that the global integral of w is 1) in space domain is defined:

        .. math::
            W(\\Psi) & = & \\frac{2 a e^{-a(1-\\cos \\Psi)}}{1 - e^{-2a}} \\\\
            a & = & \\frac{\\ln 2}{1 - \\cos (r/R)}

        here :math:`\\alpha` is the angle distance on the Earth's surface, :math:`a` is the Earth's averaging radius. :math:`r` is the distance on the Earth's surface at which :math:`W(\\alpha)` has dropped to 1/2 it value at :math:`\\alpha=0` (the distance on the Earth's surface = :math:`a*\\alpha`). Or we can refer to :math:`r` as the averaging radius (unit: meter).


    Examples:
        todo


    References:
        .. [1] Jekeli, C., Alternative methods to smooth the Earth's gravity field, Rep.327, Dep. of Geod. Sci. and Surv., Ohio State Univ., Columbus,1981.
        .. [2] Whar, J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.


    See Also:
        gauss_jekeli_fre

"""

    return _filters.py_gauss_jekeli(r, distance, ch)


def gauss_jekeli_fre(n, r):
    """
    Get Gaussian filtering weight function (defined by Jekeli) in frequency domain which is used to filter spherical harmonic coefficients 

    Args:
    n (int): Maximum number of Gaussian filtering weight function

    r (float): r is the distance on the Earth's surface at which w has dropped to 1/2 it value at alpha=0 (the distance on the Earth's surface = Earth_R*alpha). Or we can refer to r as the averaging radius (unit: meter).


    Returns:
        array_alike: w(0:n), Gaussian filtering weight function in frequency domain which is used to filter spherical harmonic coefficients


    Notes:
        Here we use the equations of C. Jekeli. So when you use the equation defined by Wahr J. (1998), note the difference for factor :math:`2*pi`. E.G., in Wahr J. (1998), eq.(30) now can omit the factor :math:`2*pi`, you can get the same results.

        In the frequency domain, the Gaussian filtering weight function can be computed with recursion relations (see Algorithm in GAUSS_JEKELI for definition of a):

        .. math::
            \mathbf{W}_0 & = & 1 \\\\
            \mathbf{W}_1 & = & \\frac{1 + e^{-2a}}{1 - e^{-2a}} - \\frac{1}{1} \\\\
            \mathbf{W}_{n+1} & = & -\\frac{2n+1}{a} \mathbf{W}_n + \mathbf{W}_{n-1}

        Now :math:`w(l)` can multiply the spherical harmonic coefficients calculated from FUN_EXP_FFT, and then use FUNC_SUM_FFT to get the Gaussian filtering weighted spherical function value in space domain. This is also equivalent to do convolution between spherical function value with Gaussian filtering weighted function :math:`w(\alpha)` calculated from GAUSS_JEKELI in space domain.

    References:
        .. [1] Jekeli, C., Alternative methods to smooth the Earth's gravity field, Rep.327, Dep. of Geod. Sci. and Surv., Ohio State Univ., Columbus,1981.
        .. [2] Wahr J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.

    Examples:
        todo

"""
    return _filters.py_gauss_jekeli_fre(n, r)


def gauss_jekeli_fre_nm(n, m_max, r_lon, r_lat):
    """Get non-isotropic Gaussian filtering weight function in frequency domain which is used to filter spherical harmonic coefficients 

    Args:
        n (int): Maximum number (Degree) of Gaussian filtering weight function
        m_order (int): Maximum order used in filter (domain: 1~N), for order grater than m_order, the filter result is 0.
        r_lon (float): The filter distance (in longitude) on the Earth's surface at which w has dropped to 1/2 it value at alpha=0 (the distance on the Earth's surface = Earth_R*alpha). Or we can refer r_lon as the averaging radius (unit: meter).
        r_lat (float): The filter distance (in latitude) on the Earth's surface at which w has dropped to 1/2 it value at alpha=0 (the distance on the Earth's surface = Earth_R*alpha). Or we can refer r_lat as the averaging radius (unit: meter).

    Returns:
        array_alike : w (0:n, 0:n), Gaussian filtering weight function in frequency domain which is used to filter spherical harmonic coefficients

    Notes:
        In the frequency domain, the non-isotropic Gaussian filtering weight function can be computed with recursion relations (see Algorithm in GAUSS_JEKELI for definition of a):

        .. math::
            \mathbf{W}_0 & = & 1 \\\\
            \mathbf{W}_1 & = & \\frac{1 + e^{-2a}}{1 - e^{-2a}} - \\frac{1}{1} \\\\
            \mathbf{W}_{n+1} & = & -\\frac{2n+1}{a} \mathbf{W}_n + \mathbf{W}_{n-1}

        Now :math:`w(l)` can multiply the spherical harmonic coefficients calculated from FUN_EXP_FFT, and then use FUNC_SUM_FFT to get the Gaussian filtering weighted spherical function value in space domain. This is also equivalent to do convolution between spherical function value with Gaussian filtering weighted function :math:`w(\alpha)` calculated from GAUSS_JEKELI in space domain.

    References:
        .. [1] Han S. C., C. k. Shum, C. Jekeli, C. Y. Kuo, C. Wilson, Non-isotropic filtering of GRACE temporal gravity for geohysical signal enhancement, Geophys. J. Int. Vol. 163, P18-25, 2005.
        .. [2] Jekeli, C., Alternative methods to smooth the Earth's gravity field, Rep.327, Dep. of Geod. Sci. and Surv., Ohio State Univ., Columbus,1981.
        .. [3] Wahr J.,M. Molenaar, F. Bryan, Time variablity of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE, J. Geophys. Res. Vol.103, No.b12, P30205-30229,1998.

    See Also:
        gauss_jekeli, gauss_jekeli_fre, func_sum_fft

"""
    return _filters.py_gauss_jekeli_fre_nm(n, m_max, r_lon, r_lat)


def weight_aver(x, npoints, weight_no=True):
    """Do n-points weighted moving (or running) average

    Args:
    x (array) : Time series to be analyzed

    npoints (int) : Number of points to do moving average (npoints must be odd number and >=3)

    weight_no (boolean) : optional, if present, then do running average (which means weight is equal); if not  then do weighted running average

    Returns:
        array_alike : Output data for moving average

    Notes:
        There is a new subroutine RUNNING_MEAN which do the same thing, but can define weight function.
        For 3 points the equation is: 
            .. math:: xout_i = x_{i-1}+3*x_{i} + x_{i+1} / 5.0
        For 5 points the equation is:
            .. math:: xout_i = (x_{i-2}+3*x_{i-1}+5*x_{i}+3*x_{i+1} + x_{i+2} / 13.0
        For npoints the weight is: `1/sum(wt), 3/sum(wt), 5/sum(wt), 7/sum(wt), 9/sum(wt), ... 3/sum(wt), 1/sum(wt)`, where :math:`wt_n=1, 3, 5, 7, 9, 11, ..., 11, 9, 7, 5, 3, 1`

    See Also:
        filters

"""
    return _filters.py_weight_aver(x, npoints, weight_no)


def butterworth(l, fl, fh, iband, n, dt, x, method='no_phase_shift', edge_effect=True):
    """Using Butterworth filter (an Infinite Impulse Response (IIR) filter ) to do signal filter for low pass, high pass, band pass and band stop filter

    Args:
        L (int) : order of normalized Butterworth analog filter. L can be given by first call WFL_BUTTERWORTH_ORDER. See example for details.
        fl (float) : lower cut-off frequency for all filter. (unit: Hz or others, see comments)
        fh (float) : higher cut-off frequency for only band-pass and band-stop filter. (unit: Hz or others, see comments)
        iband (int): Desired digital filter type: 1--low pass   2--high pass    3--band pass    4--band stop
        dt (float) : Sample intervals (unit: second or others, see comments)
        x (array_alike) : Array of input digital signal
        method (string) : if method is presented and method='no_phase_shift' or method='NO_PHASE_SHIFT', then WFL_BUTTERWORTH do twice Butterworth filter. Once forward and once backward. The purpose of doing forward-backward Butterworth filter is that by using this method, the phase shift of forward Butterworth filter will be erase (we can get zero phase shift output). Note, by using forward-backward Butterworth will double the order L of normalized Butterworth analog filter.
        edge_effect (boolean) : must be presented with "method" parameter at the same time. If presented and edge_effect=='NO', then the large edge effects during the filter process will be mostly eliminated (though there are also some small edge effects near both ends of the filtered data ).

    Returns:
        tuple
            * xout (array_alike) : signal array after Butterworth filter
            * fre(n/2) (array_alike) : Frequency of Gain function (unit: Hz or others, see comments). Note, fre, amp, and ph should be presented at the same time.
            * amp(n/2) (array_alike) : Amplitude of Gain function
            * ph(n/2) (array_alike) : Phase of Gain function (unit: degree)

    Notes:
        After the Butterworth-filter, the filtered data has phase shift to compare with the original data if method is not presented. If method is presented and method\='no_phase_shift' or method='NO_PHASE_SHIFT', then the filtered data will have the same (or almost the same) phase with the original data, but the order L will be doubled.
        There are also edge effects for the filtered data if edge_effect parameters are not presented with "NO" value. Be careful when use the data at the beginning and at the end of input signal if you do NOT use this parameter.
        Table 3.1 of Reference: Frequencies of Discrete Fourier Transform (DFT) components in different units

        +---------------------------------------------+-------------+----------------+----------------+
        |                    rad/s                    |     Hz-s    |       rad      |       Hz       |
        +---------------------------------------------+-------------+----------------+----------------+
        |            Symbol :math:`\Omega`            |  :math:`Y`  | :math:`\Omega` |    :math:`f`   |
        +---------------------------------------------+-------------+----------------+----------------+
        |  Frequency of :math:`X_1` :math:`2\pi/N/dt` | :math:`1/N` | :math:`2\pi/N` | :math:`1/N/dt` |
        +---------------------------------------------+-------------+----------------+----------------+
        | Frequency of :math:`X_{n/2}` :math:`\pi/dt` |     0.5     |   :math:`\pi`  | :math:`0.5/dt` |
        +---------------------------------------------+-------------+----------------+----------------+
        |      Sampling Frequency :math:`2\pi/dt`     |     1.0     |  :math:`2\pi`  |  :math:`1/dt`  |
        +---------------------------------------------+-------------+----------------+----------------+

        The unit of fl and fh is in Hz unit. if dt is in second unit, then fl or fh is in Hz unit; if dt is in day unit, then fl or fh is in 1/day unit (or cycle per day). Thus dt can be year or other unit, while fl or fh has the corresponding unit.

        Gain function is defined as:
            .. math::
                    Amplitude\_gain & = & | H(e^{j\omega}) | \\\\
                    Power\_gain\_in\_decibels & = & 10 \log10  | H(e^{j\omega}) |^2 \\\\
                    Phase\_shift & = & \\tan^{-1}{ Im ( | H^(ej\omega) | ) / Re ( | H(e^{j\omega}) | ) }

    References:
        .. [1] Stearns S. D. and R. A. David, Signal processing algorithms using Fortran and c,  PTR Prentice Hall, Englewood Cliffs, New Jersey, 1993. Chapter1-7

    Examples:
        todo

    See Also:
        butterworth_order

"""
    edge_effect = 'yes' if edge_effect else 'no'
    return _filters.py_wfl_butterworth(l, fl, fh, iband, n, dt, x, method, edge_effect)


def butterworth_order(fl, fh, dt):
    """Determine the best order of Butterworth filter (an Infinite Impulse Response (IIR) filter ) 

    Args:
        fl (float) : Minimum order of normalized Butterworth analog filter
        fh (float) : lower cut-off frequency for all filter. (unit: Hz or others, see comments)
        dt (float) : higher cut-off frequency for all filter. (unit: Hz or others, see comments)

    Returns:
        int : Desired digital filter type: 
            1, means low pass, 
            2, means high pass, 
            3, means band pass,
            4, means band stop.


    Notes:
        After the Butterworth-filter, the filtered data has phase shift to compare with the original data if method is not presented. If method is presented and method='no_phase_shift' or method='NO_PHASE_SHIFT', then the filtered data will have the same (or almost the same) phase with the original data. , but the order L will be doubled.

        There are also edge effects for the filtered data if edge_effect parameters are not presented with "NO" value. Be careful when use the data at the beginning and at the end of input signal if you do NOT use this parameter.

        Table 3.1 of Reference: Frequencies of Discrete Fourier Transform (DFT) components in different units

        +---------------------------------------------+-------------+----------------+----------------+
        |                    rad/s                    |     Hz-s    |       rad      |       Hz       |
        +---------------------------------------------+-------------+----------------+----------------+
        |            Symbol :math:`\Omega`            |  :math:`Y`  | :math:`\Omega` |    :math:`f`   |
        +---------------------------------------------+-------------+----------------+----------------+
        |  Frequency of :math:`X_1` :math:`2\pi/N/dt` | :math:`1/N` | :math:`2\pi/N` | :math:`1/N/dt` |
        +---------------------------------------------+-------------+----------------+----------------+
        | Frequency of :math:`X_{n/2}` :math:`\pi/dt` |     0.5     |   :math:`\pi`  | :math:`0.5/dt` |
        +---------------------------------------------+-------------+----------------+----------------+
        |      Sampling Frequency :math:`2\pi/dt`     |     1.0     |  :math:`2\pi`  |  :math:`1/dt`  |
        +---------------------------------------------+-------------+----------------+----------------+

        The unit of fl and fh is in Hz unit. if dt is in second unit, then fl or fh is in Hz unit; if dt is in day unit, then fl or fh is in 1/day unit (or cycle per day). Thus dt can be year or other unit, while fl or fh has the corresponding unit.

        Gain function is defined as:
            .. math::
                    Amplitude\_gain & = & | H(e^{j\omega}) | \\\\
                    Power\_gain\_in\_decibels & = & 10 \log10  | H(e^{j\omega}) |^2 \\\\
                    Phase\_shift & = & \tan^{-1}{ Im [ | H^(ej\omega) | ] / Re [ | H(e^{j\omega}) | ] }

    References:
        .. [1] Stearns S. D. and R. A. David, Signal processing algorithms using Fortran and c,  PTR Prentice Hall, Englewood Cliffs, New Jersey, 1993. Chapter1-7

    Examples:
        todo

    See Also:
        butterworth_order

"""
    return _filters.py_wfl_butterworth_order(fl, fh, dt, iband)


def running_mean(x, m, weight=None):
    """**running\_mean** ( real / double / complex / double complex )
    Do *n*-points weighted moving (or running) average using convolution
    method

    Args:
        x (array_alike) : Time series to be analyzed
        m (int) :Number of points to do moving average (**m** must be odd number and >=3)
        weight (array_alike) : (Optional; Input) If NOT present *weight*, then do running average (which means weight is equal); if present, then do weighted running average. For real/double precision call, the final weight is weight/sum(weight); for complex/double complex call, the final weight is weight/sum(abs(weight)). Thus the weight is normalized in the inner of subroutine WFL\_RUNNING\_MEAN to have a total weight of unity.

    Returns:
        array
            Output data for moving average


    Note
    ----
        weight=1.0 means running average; weight=[1, 3, 5, 7, 9..or anything
        like 8,4, 5,2, ...] means weighted running average, the true weight
        used in the **WFL\_RUNNING\_MEAN** is weight/sum(weight) for
        real/double calling, and weight/sum(abs(weight)) for complex/double
        complex calling.  **WFL\_RUNNING\_MEAN** is the same as
        `WEIGHT\_AVER <weight_aver.htm>`__ but using different method. Note,
        here the weight can be set as your own, while in
        `WEIGHT\_AVER <weight_aver.htm>`__, you only have two type weight.


        xout=CONV(x,weight/sum(weight)) #for real/double

        xout=CONV(x, weight/sum(abs(weight)) ) #for complex/double complex

    Examples:

        　todo

    See Also:
        filters
"""
    if weight is not None:
        if len(weight) != m:
            raise ValueError('weigt length no equal to window m')
        return _filters.py_wfl_running_mean(x, weight)
    else:
        from numpy import ones
        return _filters.py_wfl_running_mean(x, ones(m))
