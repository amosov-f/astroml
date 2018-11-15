import math
import pandas as pd
import numpy as np
import healpy as hp
import common.spherical
from sklearn import datasets, linear_model
from matplotlib import pyplot as plt
import statsmodels.api as sm

if __name__ == '__main__':
    square_arcmin = 148510660.498

    df = pd.read_table('../../doc/density32.tsv', dtype=None, index_col=0)
    dens = df['count']
    nside = int(math.sqrt(len(dens) / 12))

    ncoeff = 10

    X = np.zeros((len(dens), ncoeff))
    y = np.zeros((len(dens)))
    for ipix in range(0, len(dens)):
        for j in range(0, ncoeff):
            (theta, phi) = hp.pix2ang(nside, ipix)
            (l, b) = (phi, - theta + math.pi / 2)
            X[ipix][j] = spherical.fvj(j, l, b)
            y[ipix] = dens[ipix] / (square_arcmin / len(dens))

    # np.save('../tmp/dens.txt', X)

    print len(X)

    # regr = linear_model.LinearRegression(alpha=0.13, fit_intercept=False)
    res = sm.OLS(y, X).fit()

    # regr.fit(X, y)
    # print regr.coef_


    # res = ols.fit()  # / regr.coef_[1]

    print res.summary()

    # for i in range(0, len(res)):
    #     print "%d %f" % (i, res[i])

    if True:
        last = 0
        for i in range(len(res)):
            if abs(res[i]) > 0.00001:
                last = i
        res = res[:last + 1]

        plt.bar(range(len(res)), res)
        plt.xlabel('j')
        plt.ylabel('coeff')
        plt.xticks(range(len(res)), rotation='90')
        plt.show()

    if False:
        dens_predicted = np.zeros((len(dens)))
        for ipix in range(0, len(dens)):
            dens_predicted[ipix] = regr.predict(X[ipix])

        fig = plt.figure(1, figsize=(10, 7.5))
        hp.mollview(dens_predicted, title="Predicted density", flip='astro', fig=1)

        plt.show()
