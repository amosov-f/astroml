from matplotlib import pyplot as plt
import numpy as np

from cartesian.common_cartesian import read_radec, read_and_filter


def main():
    df = read_and_filter()

    bins = []
    slice_size = 2000000

    l1 = []
    l2 = []

    for i in range(0, len(df), slice_size):
        l = 1000 / df.iloc[i].parallax
        ri = min(i + slice_size, len(df) - 1)
        r = 1000 / df.iloc[ri].parallax
        # print(int(l), int(r), int(1000 / df.iloc[(i + ri) // 2].parallax))
        l1.append(int(l))
        l2.append(int(r))
        print(f"{int(l)}-{int(r)}")
        bins.append(l)
    bins.append(1000 / df.iloc[len(df) - 1].parallax)

    for l in l1:
        print(l, end=' ')
    print()
    for l in l2:
        print(l, end=' ')
    # print(1000 / df.iloc[len(df) - 1].parallax)

    # Create a histogram by providing the bin edges (unequally spaced).
    # bins = [1, 316, 462, 590, 716, 847, 990, 1157, 1364, 1645, 2058, 2616, 3285, 4100, 5200, 7013]
    # x = [slice_size] * len(bins)
    # mu_x = 200
    # sigma_x = 25
    x = []
    for b in bins[:-1]:
        x.extend([b] * 2000000)

    plt.hist(x, bins, histtype='bar', rwidth=0.9)
    plt.title('Гистограмма распределения звезд по выборкам')
    plt.savefig('slices.png')
    # plt.set_title('bar, unequal bins')


if __name__ == '__main__':
    main()
