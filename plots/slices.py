from matplotlib import pyplot as plt
import numpy as np

def main():
    # Create a histogram by providing the bin edges (unequally spaced).
    bins = [3,	208,	300,	386,	474,	571,	687,	835,	1040,	1303,	1594,	1897,	2220,	2582,	3031, 3677]
    x = [400000] * len(bins)
    mu_x = 200
    sigma_x = 25
    x = []
    for b in bins[:-1]:
        x.extend([b] * 400000)

    plt.hist(x, bins, histtype='bar', rwidth=0.9)
    plt.title('Гистограмма распределения звезд по выборкам')
    plt.savefig('slices.png')
    # plt.set_title('bar, unequal bins')

if __name__ == '__main__':
    main()