import healpy as hp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def compute_pixel(nside, l, b):
    return hp.ang2pix(nside, np.radians(-b + 90), np.radians(l))


def main():
    nside = 32
    count = np.zeros(hp.nside2npix(nside))

    df2 = pd.read_table('outliers.tsv', dtype=None, index_col=0)
    for indx, row in df2.iterrows():
        pixel = compute_pixel(nside, row['ecl_lon'], row['ecl_lat'])
        count[pixel] += 1

    hp.mollview(map=count, title='Outliers', flip='astro', fig=1, coord='E', min=0)
    plt.show()


if __name__ == '__main__':
    main()
