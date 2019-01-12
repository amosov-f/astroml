import healpy as hp
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas import DataFrame

from common import spherical
from common.healpix import pixel_centers


def decompose_spherical(df: DataFrame, ncoeff: int):
    X = prepare_features(df, ncoeff)
    y = df.y

    res = sm.OLS(y, X).fit()

    print(res.summary())

    return res


def prepare_features(df: DataFrame, ncoeff: int):
    X = DataFrame()
    for j in range(ncoeff):
        X[str(j)] = spherical.fvj(j, df.l, df.b)
    return X


def show_spherical_decomposition(model):
    coeffs = model.params
    errors = model.bse

    for i in range(len(coeffs)):
        # if abs(coeffs[i]) > 3 * errors[i]:
        print('{}\t{:.2f}\t{:.2f}'.format(i, coeffs[i], errors[i]).replace('.', ','))
    print()
    res = coeffs

    plt.bar(range(len(res)), res)
    plt.xlabel('j')
    plt.ylabel('coeff')
    plt.xticks(range(len(res)), rotation='90')
    plt.show()

def show_spherical_decomposition_on_sphere(model, title):
    nside = 8
    df = pixel_centers(nside)
    X = prepare_features(df, len(model.params))
    y_predicted = model.predict(X)

    plt.figure(2, figsize=(10, 7.5))
    hp.mollview(y_predicted, title=title, flip='astro', fig=2)

    plt.show()