import numpy as np
import cProfile


def cholesky(A, y):
    C = np.dot(A.T, A)
    d = np.dot(A.T, y)
    L = np.linalg.cholesky(C)
    z = np.linalg.solve(L, d)
    x = np.linalg.solve(L.T, z)
    v = np.dot(A, x) - y
    return x


def qr_solve(A, y):
    q, r = np.linalg.qr(A)
    qy = np.dot(q.T, y)
    x = np.linalg.solve(r, qy)
    v = np.dot(A, x) - y
    return x


def adjustment(A, y):
    AtA = np.dot(A.T, A)
    AtY = np.dot(A.T, y)
    Q = np.linalg.inv(AtA)
    x = np.dot(Q, AtY)
    v = np.dot(A, x) - y


def lstsq(A, y):
    return np.linalg.lstsq(A, y)


def curvefit(t, y):


n = 100000
t = np.linspace(0, 20, n)
y = t * 0.0001 + 0.0001 * t**2 + t**3 + 20 * np.sin(2 * np.pi * t / 365) + 30 * np.sin(2 * np.pi * t / 182.5) + 1.2

A = np.zeros((n, 8))
A[:, 0] = 1
A[:, 1] = t
A[:, 2] = np.sin(2 * np.pi * t / 365)
A[:, 3] = np.cos(2 * np.pi * t / 365)
A[:, 4] = np.sin(2 * np.pi * t / 182.5)
A[:, 5] = np.cos(2 * np.pi * t / 182.5)
A[:, 6] = t**2
A[:, 7] = t ** 3


# print np.linalg.lstsq(A, y)
# print np.dot(np.linalg.inv(np.dot(A.T, A)), np.dot(A.T, y))
# print cholesky(A, y)

# cProfile.run("qr_solve(A, y)")
# print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# cProfile.run("cholesky(A, y)")
# print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# cProfile.run("np.linalg.lstsq(A, y)")
# print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# cProfile.run("adjustment(A, y)")
# print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

import timeit

print timeit.timeit(stmt="qr_solve(A, y)", setup="from __main__ import qr_solve, A, y", number=100)
# print timeit.timeit(stmt="cholesky(A, y)", setup="from __main__ import cholesky, A, y", number=100)
print timeit.timeit(stmt="adjustment(A, y)", setup="from __main__ import adjustment, A, y", number=100)
# print timeit.timeit(stmt="lstsq(A, y)", setup="from __main__ import lstsq, A, y", number=100)
