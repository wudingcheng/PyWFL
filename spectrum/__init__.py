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


import numpy as np
data = np.loadtxt(r'C:\Users\wdc\Desktop\gst\sample data\nino3sst.DAT')
# r = arco(data, m=0)
# cor = r[1]
# print r
# # print (cor[1] + np.sqrt(cor[2])) / 2.0
print _spectrum.py_emd.__doc__
