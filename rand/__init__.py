import _rand


def color_noise(n, beta, iseed=0):
    """Get time series of (color) noise with pre-determined power spectra :math:`P(f)` which is proportional to :math:`1/f^{\\beta}`.

    Args:
        n (int) : number of data to be generated
        beta (float): spectral densities parameter that vary as powers of inverse frequency. More precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\\beta}`. Note, :math:`0<=\\beta<=3`.
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : (color) noise time series (generated from random number) with special spectra :math:`P(f) \propto 1 / f^{\\beta}`

    Notes:
        The **same random seed** will produce the same noise time series.
        Noise (or color noise) will be generated that has spectral densities that vary as powers of inverse frequency, more precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\\beta}` for :math:`\\beta >= 0`.  When :math:`\\beta` is 0 the noise is referred to white noise, when it is 2 it is referred to as Brownian noise, and when it is 1 it normally referred to simply as :math:`1/f` noise which occurs very often in processes found in nature.

    Reference:
        .. [1] J. Timmer and M. Konig, On generating power law noise, Astron. Astrophys., 1995, Vol. 300, P707-710.

        If one plots the log(power) vs log(frequency) for noise with a particular beta value, the slope gives the value of beta. This leads to the most obvious way to generate a noise signal with a particular beta. One generates the power spectrum with the appropriate distribution and then inverse fourier transforms that into the time domain, this technique is commonly called the fBm method and is used to create many natural looking fractal forms. The basic method involves creating frequency components which have a magnitude that is generated from a Gaussian white process and scaled by the appropriate power of f. The phase is uniformly distributed on 0, 2pi. See the code for more precise details. There is a relationship between the value of beta and the fractal dimension D, namely, D = (5 - beta) / 2

    Examples:
        >>> x = color_noise(1024, 2, iseed=1)
        >>> y = color_noise(1024, 2, iseed=1) # x and y are same, because they have same iseed.
        >>> z = color_noise(1024, 1)

"""

    return _rand.py_wfl_noise_color(n, beta, iseed=iseed)


def multi_color_noise(n, beta, pieces, iseed=0):
    """Get time series of multi-color noise with pre-determined power spectra :math:`P(f)` which is proportional to :math:`1/f^{\\beta}`.

    Args:
        n (int) : number of data to be generated
        beta (array): beta(m), spectral densities parameter that vary as powers of inverse frequency. More precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\\beta}`. Note, :math:`0<=\\beta<=3`.
        pieces (array) : The frequency index number of FFT. For example, if we want generate multi-noise color time series, the power spectrum of the time series have the below characteristics: for low frequency beta is 0.5, for middle frequency beta is 1.0, for high frequency beta is 2.0. Then beta(1:m)=beta(1:3)=0.5, 1.0, 2.0; pieces(1:m-1)=pieces(1:2)=100, 300 for n =2048. This means for FFT frequency (1~100)/n we use beta=0.5; for FFT frequency (101~300)/n we use beta =1.0,; for FFT frequency (301~1024+1)/n we use beta=2.0. Note: pieces(1:m-1) must be monotonic and pieces(m-1)<=n/2+1
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : (color) noise time series (generated from random number) with special spectra :math:`P(f) \propto 1 / f^{\\beta}`

    Notes:
        The **same random seed** will produce the same noise time series.
        Noise (or color noise) will be generated that has spectral densities that vary as powers of inverse frequency, more precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\\beta}` for :math:`\\beta >= 0`.  When :math:`\\beta` is 0 the noise is referred to white noise, when it is 2 it is referred to as Brownian noise, and when it is 1 it normally referred to simply as :math:`1/f` noise which occurs very often in processes found in nature.

    Reference:
        .. [1] J. Timmer and M. Konig, On generating power law noise, Astron. Astrophys., 1995, Vol. 300, P707-710.

        If one plots the log(power) vs log(frequency) for noise with a particular beta value, the slope gives the value of beta. This leads to the most obvious way to generate a noise signal with a particular beta. One generates the power spectrum with the appropriate distribution and then inverse fourier transforms that into the time domain, this technique is commonly called the fBm method and is used to create many natural looking fractal forms. The basic method involves creating frequency components which have a magnitude that is generated from a Gaussian white process and scaled by the appropriate power of f. The phase is uniformly distributed on 0, 2pi. See the code for more precise details. There is a relationship between the value of beta and the fractal dimension D, namely, D = (5 - beta) / 2

    Examples:
        >>> import matplotlib.pyplot as plt
        >>> x = multi_color_noise(2048, [0.5, 1.0, 2.0], [100, 300])
        >>> plt.plot(x)
        >>> plt.show()

    See Also:
        :func:`color_noise`
"""
    return _rand.py_wfl_noise_multi_color(n, beta, pieces, iseed=iseed)


def random_mt(n, interval=[0, 1], ran_type=3, iseed=0):
    """Produces uniform distribution random numbers by MT (Mersenne Twister) method.

    Args:
        n (int) : Number of normal distribution
        interval (list) : interval[2], The interval domain for uniform distribution uni
        ran_type (int) : ran_type is one of 1, 2, 3, 4. Default is ran_type==3. This parameter choice use which uniform random to produce normal random data. The uniform random period is: It is proved that the period is 219937-1, and 623-dimensional equidistribution property is assured. [0,1],[0,1) and (0,1) for ran_type= 1,2 and 3 in 32 bit precision, and (0,1) for ran_type=4 in 53 bit precision.
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : Random array

    Notes:
        Give input n then produces standard normal distribution random numbers x(n). If iseed is 0, then get the iseed by the system clock. mean and sig is the mean and standard deviation (or root of variance) of normal distribution which means Normal (mean, sig*sig)

        A random generate based on a C-program for MT19937 with initialization improved 2002/1/26.
        C-program coded by Takuji Nishimura and Makoto Matsumoto,
        FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
        Fortran90 coded base on the above Fortran77 version, and re-interfaced to easy use by Dr. Yan Haoming(2015.4.29)

    References:
        ..[1] http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

    See Also:
        :func:`random_mt_norm`
"""
    return _rand.py_wfl_random_mt(n, iseed, min(interval), max(interval), ran_type)


def random_mt_norm(n, mean=0.0, sig=1.0, ran_type=3, iseed=0):
    """Produces standard normal distribution random numbers  by MT (Mersenne Twister) method

    Args:
        n (int) : Number of normal distribution
        mean (float) : Mean value (or mathematic expectation) of normal distribution
        sig (float) : Standard deviation (or root of variance) of normal distribution.
        ran_type (int) : ran_type is one of 1, 2, 3, 4. Default is ran_type==3. This parameter choice use which uniform random to produce normal random data. The uniform random period is: It is proved that the period is 219937-1, and 623-dimensional equidistribution property is assured. [0,1],[0,1) and (0,1) for ran_type= 1,2 and 3 in 32 bit precision, and (0,1) for ran_type=4 in 53 bit precision.
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : Random array

    Notes:
        Give input n then produces standard normal distribution random numbers x(n). If iseed is 0, then get the iseed by the system clock. mean and sig is the mean and standard deviation (or root of variance) of normal distribution which means Normal (mean, sig*sig)

    See Also:
        :func:`random_mt`
"""
    return _rand.py_wfl_random_mt_norm(n, iseed, mean, sig, ran_type)


def random_norm(n, mean=0.0, sig=1.0, ran_type=1, iseed=0):
    """Produces standard normal distribution random numbers

    Args:
        n (int) : Number of normal distribution
        mean (float) : Mean value (or mathematic expectation) of normal distribution
        sig (float) : Standard deviation (or root of variance) of normal distribution.
        ran_type (int) : ran_type is one of -1, 0, 1, 2, 3. Default is ran_type==1. This parameter choice use which uniform random to produce normal random data. The uniform random period and speed are list as:

            * uniform random period: 3.1*1018 for ran_type=-1 !fastest for serial computer
            * 2.0*1028 for ran_type= 0 !paralle random, fastest
            * 8.5*1037 for ran_type= 1 !20% slower than ran_type=0 (default choose)
            * 8.5*1037 for ran_type= 2 !20% slower than ran_type=1
            * 1.8*1019 for ran_type= 3 !200% slower than ran_type 1

        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : Random array

    Notes:
        Give input n then produces standard normal distribution random numbers x(n). If iseed is 0, then get the iseed by the system clock. mean and sig is the mean and standard deviation (or root of variance) of normal distribution which means Normal (mean, sig*sig)

    See Also:
        :func:`random_mt`, :func:`random_mt_norm`

"""
    return _rand.py_wfl_random_norm(n, iseed, mean, sig, ran_type)

