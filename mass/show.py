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

    df = prepare_dataset(df)

    # df['mul'] = df['mul'] + (22.6 / (K * (df.r / 1000)))

    df['v'] = K * df['mul'] * np.abs(df.r) / 1000
    # df = df[(df.r > 0) & (df.r < 1000)]#.plot.scatter(x='r', y='mul')

    fig, axes = plt.subplots(nrows=2, ncols=1)


    dist = [-5000, 8000]
    speed = [-100, 50]
    bx = 500
    by = 100

    mul_range = [-10, 5]
    axes[0].set_title('Угловая скорость по галактической долготе')
    h0 = axes[0].hist2d(df.r, df['mul'], range=[dist, mul_range], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df['mul'], bins=bx, range=dist, statistic='median')
    axes[0].hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label='binned statistic of data')
    axes[0].set_ylabel('Угловая скорость [mas/год]')
    plt.colorbar(h0[3], ax=axes[0])

    axes[1].set_title('Абсолютная скорость относительно Солнца')
    h1 = axes[1].hist2d(df.r, df.v, range=[dist, speed], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df.v, bins=bx, range=dist, statistic='median')
    axes[1].hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label = 'binned statistic of data')
    axes[1].set_xlabel('Расстояние [пк]')
    axes[1].set_ylabel('Абсолютная скорость [км/с]')
    plt.colorbar(h1[3], ax=axes[1])


    # fig, axes = plt.subplots(nrows=2, ncols=2)
    # axes[0, 0].set_title('Расстояние до звезд [пк]')
    # axes[1].colorbar()

    plt.show()


if __name__ == '__main__':
    main()