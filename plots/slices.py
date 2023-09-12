from matplotlib import pyplot as plt
import numpy as np


def main():
    # Create a histogram by providing the bin edges (unequally spaced).
    bins = [1, 316, 462, 590, 716, 847, 990, 1157, 1364, 1645, 2058, 2616, 3285, 4100, 5200, 7013]
    x = [2000000] * len(bins)
    mu_x = 200
    sigma_x = 25
    x = []
    for b in bins[:-1]:
        x.extend([b] * 2000000)

    plt.hist(x, bins, histtype='bar', rwidth=0.9)
    plt.title('Гистограмма распределения звезд по выборкам')
    plt.savefig('slices.png')
    # plt.set_title('bar, unequal bins')


if __name__ == '__main__':
    main()
