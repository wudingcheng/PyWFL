import _test
print _test.lsfit.__doc__
import numpy as np


def func(x, a):
    return x**a

x = np.linspace(0, 10, 11)
y = func(x, 2)
print x
print y
print _test.pass_func(x, func)
# print _test.lsfit(x, y, 1, func)
