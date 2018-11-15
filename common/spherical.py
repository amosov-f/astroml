from numpy import sin, cos, tan, sqrt, pi


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
        return 0.0

    x = sin(d)

    f = (n - k) // 2

    if f == 0:
        z = x ** (n - k) * sqrt(1 - x ** 2) ** k
    else:
        z = x ** (n - k)
        for m in range(1, f + 1):
            z += fa(m, n, k) * x ** (n - k - 2 * m)
        z = z * sqrt(1 - x ** 2) ** k

    return norm(n, k) * z


def fr(n, k):
    res = sqrt((2 * n + 1) / (4.0 * pi))
    if k > 0:
        res *= sqrt(2.0 / fact2(n - k + 1, n + k))
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
    n = int(sqrt(j))
    k = j - n ** 2
    if k % 2 == 0:
        l = 1
    else:
        l = 0
    k = (k - l + 1) // 2
    return n, k, l


def indexj(n, k, l):
# Вычисление индекса j = n ** 2 + 2 * k + k - 1
    return n ** 2 + 2 * k + l - 1


def ta(n, k, l, a, d):
    # Компонент торроидальной функции по A
    if k == 0:
        ta = pl(n, 1, d)
    else:
        if l == 0:
            ta = (-k * tan(d) * pl(n, k, d) + pl(n, k + 1, d)) * sin(k * a)
        else:
            ta = (-k * tan(d) * pl(n, k, d) + pl(n, k + 1, d)) * cos(k * a)
    return fr(n, k) / sqrt(float(n * (n + 1))) * ta


def taj(j, a, d):
  # Компонент торроидальной функции по A от одного индекса
    n, k, l = indexes(j)
    return ta(n, k, l, a, d)


def td(n, k, l, a, d):
  # Компонент торроидальной функции по D
    if k == 0:
        td = 0.0 * d
    else:
        if l == 0:
            td = -k / cos(d) * pl(n, k, d) * cos(k * a)
        else:
            td = +k / cos(d) * pl(n, k, d) * sin(k * a)

    return fr(n, k) / sqrt(float(n * (n + 1))) * td


def tdj(j, a, d):
    # Компонент торроидальной функции по D от одного индекса
    n, k, l = indexes(j)
    return td(n, k, l, a, d)


def sa(n, k, l, a, d):
 # Компонент сфероидальной функции по A
    if k == 0:
        sa = 0.0 * a
    else:
        if l == 0:
            sa = +k / cos(d) * pl(n, k, d) * cos(k * a)
        else:
            sa = -k / cos(d) * pl(n, k, d) * sin(k * a)
    return fr(n, k) / sqrt(float(n * (n + 1))) * sa


def saj(j, a, d):
    # Компонент сфероидальной функции по A от одного индекса
    n, k, l = indexes(j)
    return sa(n, k, l, a, d)


def sd(n, k, l, a, d):
    # Компонент сфероидальной функции по D
    if k == 0:
        sd = pl(n, 1, d)
    else:
        if l == 0:
            sd = (-k * tan(d) * pl(n, k, d) + pl(n, k + 1, d)) * sin(k * a)
        else:
            sd = (-k * tan(d) * pl(n, k, d) + pl(n, k + 1, d)) * cos(k * a)
    return fr(n, k) / sqrt(float(n * (n + 1))) * sd


def sdj(j, a, d):
    # Компонент сфероидальной функции по D от одного индекса
    n, k, l = indexes(j)
    return sd(n, k, l, a, d)