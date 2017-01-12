from ctypes import *
import numpy as np

_add = cdll.LoadLibrary("./add.o")


def add(a, b):
    a = np.require(a, float, ['CONTIGUOUS', 'ALIGNED'])
    b = np.require(b, float, ['CONTIGUOUS', 'ALIGNED'])
    c = np.empty_like(a)
    _add.add(a, b, c)
    return c

print add(1, 2)
