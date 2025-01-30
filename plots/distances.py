import numpy as np
from matplotlib import pyplot as plt

def uniform(x, a, b):
    if x < a or x > b:
        return 0
    return (1 / (b - a))


def inverse(f, x):
    return 1 / x ** 2 * f(1 / x)

def integral(f, a, b, parts=100000):
    sum = 0
    dx = (b - a) / parts
    for i in range(parts):
        x = (i + 0.5) * dx + a
        y = f(x)
        sum += y
    return sum * dx


def exponent(x, mu, sigma):
    return (1 / (sigma * (2 * np.pi) ** 0.5)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def lognormal(x, mu, sigma):
    # if x.any() <= 0:
    #     raise Exception('x should be > 0, but got: ' + str(x))
    return (1 / (x * sigma * (2 * np.pi) ** 0.5)) * np.exp(- (np.log(x) - mu) ** 2 / (2 * sigma ** 2))

def sum_predics(f1, f2, z, f1_min, f1_max):
    res = integral(lambda x: f1(x) * f2(z - x), f1_min, f1_max)
    return res


def solve(sigma):
    plt.clf()

    r1 = 0.1
    r2 = 8.0
    r_real = np.linspace(r1, r2, 1000)

    # u = lambda x: uniform(x, r1, r2)
    u = lambda x: lognormal(x, 0, 0.5)


    inv_u = lambda x: inverse(u, x)
    # inv_u_min = 0
    # inv_u_max = 1 /

    e = lambda x: exponent(x, 0, sigma)

    pi = lambda z: sum_predics(inv_u, e, z, 0, 100)

    # r_observed_f = lambda x: inverse(pi, x)
    r_observed = inverse(pi, r_real)

    # print(integral(r_observed_f, -2, 20))

    plt.title(f'sigma = {sigma} mas')
    plt.plot(r_real, r_observed, label='observed')
    plt.plot(r_real, u(r_real), label='real')

    plt.xlabel('кпк')
    plt.ylabel('доля звезд')

    plt.legend()

    plt.savefig(f'sigma{sigma}.png')
    #
    plt.show()


def main():
    for sigma in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1]:
        solve(sigma)

if __name__ == '__main__':
    main()