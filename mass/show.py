from common.astro import galaxy_mu
from common.gaia.with_rv import to_galaxy, distance
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import stats

K = 4.74

def prepare_dataset(df):

    l, b, mul, mub, px = to_galaxy(df)
    return pd.DataFrame(data={
        'l': l, 'b': b, 'mul': mul, 'mub': mub, 'px': px, 'r': 1000.0 / px
    })

def main():
    to_center = pd.read_csv('to_center.tsv', sep='\t')
    from_center = pd.read_csv('from_center.tsv', sep='\t')
    from_center['parallax'] = -from_center.parallax
    df = pd.concat([from_center, to_center])

    # solar_pm = 230 / (K * 8500)
    # print(solar_pm)

    df = prepare_dataset(df)

    df['v'] = K * np.abs(df.r) * df['mul']
    # df = df[(df.r > 0) & (df.r < 1000)]#.plot.scatter(x='r', y='mul')

    # x = np.linspace(0, 1000, 10)
    # y = np.linspace(0, 10, 100)
    # H, xedges, yedges = np.histogram2d(df['r'], df['mul'], bins=(x, y))
    dist = [-3000, 8500]
    speed = [-100000, 50000]
    bx = 500
    by = 100
    plt.figure(1)
    plt.hist2d(df.r, df.v, range=[dist, speed], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df.v, bins=bx, range=dist, statistic='median')
    plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label = 'binned statistic of data')
    plt.colorbar()

    mul_range = [-10, 5]
    bx = 500
    plt.figure(2)
    plt.hist2d(df.r, df['mul'], range=[dist, mul_range], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df['mul'], bins=bx, range=dist, statistic='median')
    plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label='binned statistic of data')
    plt.colorbar()

    plt.show()


if __name__ == '__main__':
    main()