import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


def main():
    # mag_histogram()
    # parallax_histogram()
    parallax_error()


def parallax_error():
    dataset = pd.read_csv('parallax_error.tsv', sep='\t')
    plt.figure()
    dataset['distance'] = 1000.0 / dataset['parallax']
    dataset['error'] = 100 * dataset['parallax_error'] / dataset['parallax']
    dataset = dataset[((dataset.distance < 2000) & (dataset.distance > 0) & (dataset.error > 0) & (dataset.error < 6))]
    fig = dataset.plot.hexbin(x='distance', y='error', gridsize=50)
    fig.set_xlabel("Расстояние [пк]")
    fig.set_ylabel("Относительная ошибка параллакса [%]")
    plt.savefig('parallax_error.png')


def mag_histogram():
    dataset = pd.read_csv('full.tsv', sep='\t')
    plt.figure()
    min_mag = 7
    max_mag = 16
    step = 1
    bins = np.linspace(min_mag, max_mag, ((max_mag - min_mag) // step) + 1)
    fig = dataset['phot_g_mean_mag'].plot.hist(bins=bins, figsize=(10, 5))
    fig.set_xlabel("Звездная величина [m]")
    fig.set_ylabel("Число звезд")
    plt.savefig('mag_histogram.png')
    print_hist(np.histogram(dataset['phot_g_mean_mag'], bins=bins))


def parallax_histogram():
    dataset = pd.read_csv('full.tsv', sep='\t')
    plt.figure()
    max_dist = 4000
    bins = np.linspace(0, max_dist, (max_dist // 200) + 1)
    dataset['distance'] = 1000.0 / dataset['parallax']
    fig = dataset['distance'].plot.hist(bins=bins, figsize=(10, 5))
    fig.set_xlabel("Расстояние [пк]")
    fig.set_ylabel("Число звезд")
    plt.savefig('parallax_histogram.png')
    print_hist(np.histogram(dataset['distance'], bins=bins))


def print_hist(hist):
    for left, right, v in zip(hist[1][0:-1], hist[1][1:], hist[0]):
        print(f'{left}\t{right}\t{v}')


if __name__ == '__main__':
    main()
