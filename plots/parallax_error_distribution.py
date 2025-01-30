import numpy as np
from matplotlib import pyplot as plt


def main():
    mu, sigma = 1, 0.20 # mean and standard deviation
    s = np.random.normal(mu, sigma, 1000)
    bins = np.linspace(0, 1, num=1000)

    # count, bins, ignored = plt.hist(s, 30, density=True)
    # count, bins, ignored = plt.hist(s, 30, density=True)
    sigma = np.linspace(0, 1, num=1000)

    plt.plot(sigma, pizza(sigma))
    plt.xlabel('z')
    plt.ylabel('diff, %')

    # plt.plot(bins, 1 / (sigma * np.sqrt(2 * np.pi)) *
    #          np.exp(- (bins - mu) ** 2 / (2 * sigma ** 2)),
    #          linewidth=2, color='r')

    # plt.plot(bins, reciprocal_normal(bins, mu, sigma))
    plt.show()


def pizza(z):
    return 100 * (1 - (np.sqrt(8 * z**2 + 1) - 1) / (4 * z**2))

def reciprocal_normal(y, mu, sigma):
    return 1 / (np.sqrt(2 * np.pi) * sigma * y ** 2) * np.exp(-0.5 * ((1 / y - mu) / sigma) ** 2)

if __name__ == '__main__':
    main()
