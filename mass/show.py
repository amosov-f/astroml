from common.astro import galaxy_mu
from common.gaia.with_rv import to_galaxy, distance
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats

K = 4.74

def prepare_dataset(df, with_pm_errors=True):
    l, b, mul, mub, px = to_galaxy(df, with_pm_errors)
    return pd.DataFrame(data={
        'l': l, 'b': b, 'mul': mul, 'mub': mub, 'px': px, 'r': 1000.0 / px
    })

def filter(df):
    return df
    # return df[df.parallax > 0]

def main():
    to_center = filter(pd.read_csv('to_center_full.tsv', sep='\t'))
    from_center = filter(pd.read_csv('from_center_full.tsv', sep='\t'))
    from_center['parallax'] = -from_center.parallax
    df = pd.concat([from_center, to_center])

    # solar_pm = 230 / (K * 8500)

    df['r'] = 1000 / df.parallax
    # df = prepare_dataset_full(df)

    print(len(df))

    # df['mul'] = df['mul'] + (22.6 / (K * (df.r / 1000)))

    df['v'] = K * df['mul'] * df.r / 1000
    # df = df[(df.r > 0) & (df.r < 1000)]#.plot.scatter(x='r', y='mul')

    # fig = plt.gcf()
    plt.rcParams["figure.figsize"] = [12, 10]
    plt.rcParams.update({'font.size': 14})

    fig, axes = plt.subplots(nrows=3, ncols=1)


    dist = [-10000, 20000]
    speed = [-500, 100]
    bx = 250
    by = 100

    print(plt.style.available)
    plt.style.use('seaborn-white')


    mul_range = [-7, 4]
    axes[0].set_title('Собственное движение звезд по галактической долготе')
    h0 = axes[0].hist2d(df.r, df['mul'], range=[dist, mul_range], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df['mul'], bins=bx, range=dist, statistic='median')
    axes[0].hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label='binned statistic of data')
    axes[0].set_ylabel('Угловая скорость [мсд/год]')
    plt.colorbar(h0[3], ax=axes[0])

    axes[1].set_title('Поперечная скорость звезд по галактической долготе')
    h1 = axes[1].hist2d(df.r, df.v, range=[dist, speed], bins=(bx, by))
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.r, df.v, bins=bx, range=dist, statistic='median')
    axes[1].hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label = 'binned statistic of data')
    axes[1].set_xlabel('Расстояние [пк]')
    axes[1].set_ylabel('Абсолютная скорость [км/с]')
    plt.colorbar(h1[3], ax=axes[1])

    axes[2].set_title('Количество звезд по галактической долготе')
    h2 = axes[2].hist(df.r, range=dist, bins=bx)
    axes[2].set_xlabel('Расстояние [пк]')
    axes[2].set_ylabel('Количество звезд')



    # fig, axes = plt.subplots(nrows=2, ncols=2)
    # axes[0, 0].set_title('Расстояние до звезд [пк]')
    # axes[1].colorbar()

    # plt.figure(figsize=(20, 10), dpi=100)
    # plt.savefig('speed_curve.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    main()