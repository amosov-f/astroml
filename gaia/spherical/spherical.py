import numpy
from math import sin
from math import cos


def fact2(n1, n2):
    res = 1.0
    for i in range(n1, n2 + 1):
        res *= i
    return res


def norm(n, k):
    up = 1.0
    for i in range(n + 1, 2 * n + 1):
        up *= i
    dn = 1.0
    for i in range(1, n - k + 1):
        dn *= i
    up /= dn
    for i in range(1, n + 1):
        up *= 0.5
    return up


def fa(m, n, k):
    up = 1.0
    for i in range(0, 2 * m):
        up *= (n - k - i)
    dn = 1
    for i in range(1, m + 1):
        dn *= 2 * i * (2 * (n - i) + 1)
    if m % 2 == 0:
        return up / dn
    return -up / dn


def pl(n, k, d):
    if n < k:
        return

    x = sin(d)

    f = (n - k) / 2

    if f == 0:
        z = x ** (n - k) * numpy.sqrt(1 - x ** 2) ** k
    else:
        z = x ** (n - k)
        for m in range(1, f + 1):
            z += fa(m, n, k) * x ** (n - k - 2 * m)
        z = z * numpy.sqrt(1 - x ** 2) ** k

    return norm(n, k) * z


def fr(n, k):
    res = numpy.sqrt((2 * n + 1) / (4.0 * numpy.pi))
    if k > 0:
        res *= numpy.sqrt(2.0 / fact2(n - k + 1, n + k))
    return res


def fk(n, k, l, a, d):
    if k == 0:
        return pl(n, 0, d)
    if l == 0:
        return pl(n, k, d) * sin(k * a)
    return pl(n, k, d) * cos(k * a)


def fv(n, k, l, a, d):
    return fr(n, k) * fk(n, k, l, a, d)


def fvj(j, a, d):
    (n, k, l) = indexes(j)
    return fv(n, k, l, a, d)


def indexes(j):
    n = int(numpy.sqrt(j))
    k = j - n ** 2
    if k % 2 == 0:
        l = 1
    else:
        l = 0
    k = (k - l + 1) / 2
    return n, k, l
