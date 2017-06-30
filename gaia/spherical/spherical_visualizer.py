# coding=utf-8
from StringIO import StringIO

import healpy as hp
import numpy as np
import pandas
from matplotlib import pyplot as plt
import spherical

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


def flat_for(a, f):
    a = a.reshape(-1)
    for i, v in enumerate(a):
         a[i] = f(v)

if __name__ == '__main__':

    # df = pandas.read_table('../doc/density32.tsv', dtype=None, index_col=0)['count']

    nside = 32

    for j in range(0, 100):
        (n, k, l) = spherical.indexes(j)

        df = []
        for i in range(0, 12 * (nside ** 2)):
            (theta, phi) = hp.pix2ang(nside, i)
            df.append(spherical.fvj(j, phi, - theta + np.pi / 2))

        plt.clf()

        fig = plt.figure(1, figsize=(20, 15))
        hp.mollview(np.array(df), title=u"Spherical function for j=%d (n=%d, k=%d, l=%d)" % (j, n, k, l), flip='astro', fig=1)
        hp.graticule()

        plt.savefig('../../doc/fvj/%d.png' % j, dpi=60)

