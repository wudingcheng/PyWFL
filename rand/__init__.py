import _rand


def color_noise(n, beta, iseed=0):
    """Get time series of (color) noise with pre-determined power spectra :math:`P(f)` which is proportional to :math:`1/f^{\beta}`.

    Args:
        n (int) : number of data to be generated
        beta (float): spectral densities parameter that vary as powers of inverse frequency. More precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\beta}`. Note, :math:`0<=\beta<=3`.
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : (color) noise time series (generated from random number) with special spectra :math:`P(f) \propto 1 / f^{\beta}`

    Notes:
        The **same random seed** will produce the same noise time series.
        Noise (or color noise) will be generated that has spectral densities that vary as powers of inverse frequency, more precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\beta}` for :math:`\beta >= 0`.  When :math:`\beta` is 0 the noise is referred to white noise, when it is 2 it is referred to as Brownian noise, and when it is 1 it normally referred to simply as :math:`1/f` noise which occurs very often in processes found in nature.

    Reference:
        .. [1] J. Timmer and M. Konig, On generating power law noise, Astron. Astrophys., 1995, Vol. 300, P707-710.

        If one plots the log(power) vs log(frequency) for noise with a particular beta value, the slope gives the value of beta. This leads to the most obvious way to generate a noise signal with a particular beta. One generates the power spectrum with the appropriate distribution and then inverse fourier transforms that into the time domain, this technique is commonly called the fBm method and is used to create many natural looking fractal forms. The basic method involves creating frequency components which have a magnitude that is generated from a Gaussian white process and scaled by the appropriate power of f. The phase is uniformly distributed on 0, 2pi. See the code for more precise details. There is a relationship between the value of beta and the fractal dimension D, namely, D = (5 - beta) / 2

    Examples:
        >> x = color_noise(1024, 2, iseed=1)
        >> y = color_noise(1024, 2, iseed=1) # x and y are same, because they have same iseed.
        >> z = color_noise(1024, 1)
    
"""

    return _rand.py_wfl_noise_color(n, beta, iseed=iseed)

def multi_color_noise(n, beta, pieces, iseed=0):
    """Get time series of multi-color noise with pre-determined power spectra :math:`P(f)` which is proportional to :math:`1/f^{\beta}`.

    Args:
        n (int) : number of data to be generated
        beta (array): beta(m), spectral densities parameter that vary as powers of inverse frequency. More precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\beta}`. Note, :math:`0<=\beta<=3`.
        pieces (array) : The frequency index number of FFT. For example, if we want generate multi-noise color time series, the power spectrum of the time series have the below characteristics: for low frequency beta is 0.5, for middle frequency beta is 1.0, for high frequency beta is 2.0. Then beta(1:m)=beta(1:3)=0.5, 1.0, 2.0; pieces(1:m-1)=pieces(1:2)=100, 300 for n =2048. This means for FFT frequency (1~100)/n we use beta=0.5; for FFT frequency (101~300)/n we use beta =1.0,; for FFT frequency (301~1024+1)/n we use beta=2.0. Note: pieces(1:m-1) must be monotonic and pieces(m-1)<=n/2+1
        iseed (int) : default 0, use system clock as the random seed, different time with different random array, if not equal to 0, then the same iseed get the same random array.

    Returns:
        array : (color) noise time series (generated from random number) with special spectra :math:`P(f) \propto 1 / f^{\beta}`

    Notes:
        The **same random seed** will produce the same noise time series.
        Noise (or color noise) will be generated that has spectral densities that vary as powers of inverse frequency, more precisely, the power spectra :math:`P(f)` is proportional to :math:`1/f^{\beta}` for :math:`\beta >= 0`.  When :math:`\beta` is 0 the noise is referred to white noise, when it is 2 it is referred to as Brownian noise, and when it is 1 it normally referred to simply as :math:`1/f` noise which occurs very often in processes found in nature.

    Reference:
        .. [1] J. Timmer and M. Konig, On generating power law noise, Astron. Astrophys., 1995, Vol. 300, P707-710.

        If one plots the log(power) vs log(frequency) for noise with a particular beta value, the slope gives the value of beta. This leads to the most obvious way to generate a noise signal with a particular beta. One generates the power spectrum with the appropriate distribution and then inverse fourier transforms that into the time domain, this technique is commonly called the fBm method and is used to create many natural looking fractal forms. The basic method involves creating frequency components which have a magnitude that is generated from a Gaussian white process and scaled by the appropriate power of f. The phase is uniformly distributed on 0, 2pi. See the code for more precise details. There is a relationship between the value of beta and the fractal dimension D, namely, D = (5 - beta) / 2

    Examples:

"""
    return _rand.py_wfl_noise_multi_color(n, beta, pieces, iseed=iseed)

print _rand.py_wfl_random_mt.__doc__
import matplotlib.pyplot as plt
# x = multi_color_noise(2048, [0.5, 1.0, 2.0], [100, 300])
# plt.plot(x)
# plt.show()
