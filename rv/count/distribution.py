import os
from pathlib import Path

import numpy as np
import healpy as hp
import pandas as pd
import scipy
import sklearn
import statsmodels.api as sm
from matplotlib import pyplot as plt
from numpy import sin, cos, sqrt, radians
from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split

from GC2019.gc2019 import read_gc2019
from common.astro import galaxy, galaxy_mu
from common.gaia.with_rv import distance, slices, prepare_galaxy_dataset, read_gaia_with_rv_1500
from common.pandas2 import split



def compute_pixel(nside, l, b):
    theta = np.radians(-b + 90)
    phi = np.radians(l)
    # print(np.min(theta))
    # print(np.max(theta))
    # nan_theta = np.isnan(theta)
    # theta = theta[np.logical_not(nan_theta)]
    # phi = phi[np.logical_not(nan_theta)]
    # nan_phi = np.isnan(phi)
    # theta = theta[np.logical_not(nan_phi)]
    # phi = phi[np.logical_not(nan_phi)]
    # print(nan_theta)
    return hp.ang2pix(nside, theta, phi)


def draw_distribution(df, folder, title):
    nside = 32
    # hist = np.zeros()
    df['pix'] = compute_pixel(nside, np.degrees(df.l.values), np.degrees(df.b.values))
    npix = hp.nside2npix(nside)
    hist, _ = np.histogram(df.pix.values, range=(0, npix), bins=(npix))

    plt.clf()

    hp.mollview(hist, title=title, unit='число скоплений в пикселе',
                coord='G', min=0)
    hp.graticule()
    # plt.show()
    try:
        os.mkdir(folder)
    except:
        pass
    plt.savefig(f'{folder}/galactocentric2_{nside}_{title}.png')


def draw_function(l, b, value, max, folder=None, title=None):

    ## s310, t211, v
    nside = 8

    if folder is not None:
        folder = f"{folder}/{nside}"

    pix = compute_pixel(nside, np.degrees(l.values), np.degrees(b.values))
    npix = hp.nside2npix(nside)
    bin_means, bin_edges, binnumber = scipy.stats.binned_statistic(pix, value, bins=npix, statistic='median')
    # hist, _ = np.histogram(bin_means, range=(0, npix), bins=(npix))

    plt.clf()

    hp.mollview(bin_means, title=title, unit='медианное значение в пикселе',
                coord='G', min=-max, max=+max)
    hp.graticule()

    if folder is not None and title is not None:
        Path(folder).mkdir(parents=True, exist_ok=True)
        plt.savefig(f'{folder}/{title}.png')
    else:
        plt.show()



STEP = 75
SLICE = STEP
MAX = 6000000

def main():
    print(plt.style.available)
    plt.style.use('seaborn-white')

    df = read_gc2019()
    slices(draw_distribution, df, MAX, STEP, SLICE, prepare=False)


if __name__ == '__main__':
    main()