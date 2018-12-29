import time

import healpy as hp
import numpy
import pandas as pd
import sklearn
import statsmodels.api as sm
from matplotlib import pyplot as plt
from numpy import sin, cos, sqrt, radians
from scipy.interpolate import spline
from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split

from common.astro import galaxy, galaxy_mu
from common.gaia.with_rv import distance, slices
from common.pandas2 import split

Pi = 3.1415926535897932384626433832795_8
K = 4.74


def galaxy_apply(row):
    l, b = galaxy(row['ra'], row['dec'])
    return pd.Series([l, b], index=['l', 'b'])


def galaxy_mu_apply(row):
    mul, mub = galaxy_mu(row['pmra'], row['pmdec'], row['l'], row['b'], row['dec'])
    # mul, mub = row['pmra'], row['pmdec']
    return (mul, mub)


def kmul_apply(df):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * sin(l),
                             V=-px * cos(l),
                             W=0,
                             Wx=-sin(b) * cos(l),
                             Wy=-sin(b) * sin(l),
                             Wz=+cos(b),
                             M12=+cos(b) * cos(2 * l),
                             M13=-sin(b) * sin(l),
                             M23=+sin(b) * cos(l),
                             M11=-0.5 * cos(b) * sin(2 * l),
                             M22=+0.5 * cos(b) * sin(2 * l),
                             M33=0.0,
                             y=K * mul))


def kmul_prim_apply(df):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * sin(l),
                             V=-px * cos(l),
                             W=0.0,
                             Wx=-sin(b) * cos(l),
                             Wy=-sin(b) * sin(l),
                             Wz=+cos(b),
                             M12=+cos(b) * cos(2 * l),
                             M13=-sin(b) * sin(l),
                             M23=+sin(b) * cos(l),
                             M11z=-0.5 * cos(b) * sin(2 * l),
                             X=0.0,
                             y=K * mul))


def kmul_cb(df: pd.DataFrame):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * sin(l),
                             V=-px * cos(l),
                             W=0.0,
                             Wx=-sin(b) * cos(l),
                             Wy=-sin(b) * sin(l),
                             Wz=+cos(b),
                             M12=+cos(b) * cos(2 * l),
                             M13=-sin(b) * sin(l),
                             M23=+sin(b) * cos(l),
                             C=-cos(b) * sin(2 * l),
                             K=0.0,
                             y=K * mul))


def kmub_apply(df):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * cos(l) * sin(b),
                             V=+px * sin(l) * sin(b),
                             W=-px * cos(b),
                             Wx=sin(l),
                             Wy=-cos(l),
                             Wz=0.0,
                             M12=-0.5 * sin(2 * b) * sin(2 * l),
                             M13=+cos(2 * b) * cos(l),
                             M23=+cos(2 * b) * sin(l),
                             M11=-0.5 * sin(2 * b) * cos(l) ** 2,
                             M22=-0.5 * sin(2 * b) * sin(l) ** 2,
                             M33=+0.5 * sin(2 * b),
                             y=K * mub))


def kmub_prim_apply(df):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * cos(l) * sin(b),
                             V=+px * sin(l) * sin(b),
                             W=-px * cos(b),
                             Wx=sin(l),
                             Wy=-cos(l),
                             Wz=0.0,
                             M12=-0.5 * sin(2 * b) * sin(2 * l),
                             M13=+cos(2 * b) * cos(l),
                             M23=+cos(2 * b) * sin(l),
                             M11z=-0.25 * sin(2 * b) * cos(2 * l),
                             X=+0.5 * sin(2 * b),
                             y=K * mub))


def kmub_cb(df: pd.DataFrame):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    return pd.DataFrame(dict(U=+px * cos(l) * sin(b),
                             V=+px * sin(l) * sin(b),
                             W=-px * cos(b),
                             Wx=sin(l),
                             Wy=-cos(l),
                             Wz=0.0,
                             M12=-0.5 * sin(2 * b) * sin(2 * l),
                             M13=+cos(2 * b) * cos(l),
                             M23=+cos(2 * b) * sin(l),
                             C=-0.5 * sin(2 * b) * cos(2 * l),
                             K=-0.5 * sin(2 * b),
                             y=K * mub))


def vr_apply(df):
    l, b, mul, mub, px = df.l, df.b, df['mul'], df.mub, df.px
    x = {
        'U': -cos(l) * cos(b),
        'V': -sin(l) * cos(b),
        "W": -sin(b),
        "M12": cos(b) ** 2 * sin(2 * l) / px,
        'M13': sin(2 * b) * cos(l) / px,
        "M23": sin(2 * b) * sin(l) / px,
        "M11": cos(b) ** 2 * cos(l) ** 2 / px,
        "M22": cos(b) ** 2 * sin(l) ** 2 / px,
        "M33": sin(b) ** 2 / px,
        'y': df.radial_velocity
    }
    return pd.DataFrame(x)



def bottlinger_vr(df: pd.DataFrame, l0, R0):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    R = compute_R(l, b, r, l0, R0)
    return pd.DataFrame(dict(UG=cos(b) * cos(l),
                             VG=+cos(b) * sin(l),
                             WG=+sin(b),
                             w01=-R0 * (R - R0) * sin(l - l0) * cos(b),
                             w011=-0.5 * R0 * (R - R0) ** 2 * sin(l - l0) * cos(b),
                             k0=cos(b) ** 2 * r,
                             k01=+(R - R0) * (r * cos(b) - R0 * cos(l - l0)) * cos(b),
                             k011=+0.5 * (R - R0) ** 2 * (r * cos(b) - R0 * cos(l - l0)) * cos(b),
                             y=df.radial_velocity))

def bottlinger_mul(df: pd.DataFrame, l0, R0):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    R = compute_R(l, b, r, l0, R0)
    return pd.DataFrame(dict(UG=-sin(l),
                             VG=+cos(l),
                             w01=-(R - R0) * (R0 * cos(l - l0) - r * cos(b)),
                             w011=-0.5 * (R - R0) ** 2 * (R0 * cos(l - l0) - r * cos(b)),
                             w0=+r * cos(b),
                             k01=+R0 * (R - R0) * sin(l - l0),
                             k011=+0.5 * R0 * (R - R0) ** 2 * sin(l - l0),
                             y=K * r * mul))

def bottlinger_mub(df: pd.DataFrame, l0, R0):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    R = compute_R(l, b, r, l0, R0)
    return pd.DataFrame(dict(UG=-cos(l) * sin(b),
                             VG=-sin(l) * sin(b),
                             WG=+cos(b),
                             w01=+R0 * (R - R0) * sin(l - l0) * sin(b),
                             w011=+0.5 * R0 * (R - R0) ** 2 * sin(l - l0) * sin(b),
                             k0=-cos(b) * sin(b) * r,
                             k01=-(R - R0) * (r * cos(b) - R0 * cos(l - l0)) * sin(b),
                             k011=-0.5 * (R - R0) ** 2 * (r * cos(b) - R0 * cos(l - l0)) * sin(b),
                             y=K * r * mub))


def columns(df: pd.DataFrame):
    return df.l, df.b, df['mul'], df.mub, df.px


def compute_R(l, b, r, l0, R0):
    R2 = (r * cos(l)) ** 2 * R0 * r * cos(b) * cos(l - l0) + R0 ** 2
    return sqrt(R2)


def prepare_mu_vr(dataset: pd.DataFrame):
    kmul = kmul_cb(dataset)
    kmub = kmub_cb(dataset)
    return split(pd.concat([kmul, kmub]).fillna(0))

def prepare_bottlinger(df: pd.DataFrame, l0, R0):
    vr = bottlinger_vr(df, l0, R0)
    mul = bottlinger_mul(df, l0, R0)
    mub = bottlinger_mub(df, l0, R0)
    return split(pd.concat([vr, mul, mub]).fillna(0))


def compute(l, r):
    filter = f'{l} <= 1000 / parallax < {r}'
    print(filter)
    df = pd.read_csv('full.tsv', sep='\t') \
        .query(filter)
    print(df)
    df = df.sample(min(50000, len(df)), random_state=0)

    kmul = df.apply(kmul_apply, axis=1)
    kmub = df.apply(kmub_apply, axis=1)
    # vr = df.apply(vr_apply, axis=1)
    full = pd.concat([kmul, kmub])
    # full = vr
    print(full)

    X = full.iloc[:, :-1]
    y = full.iloc[:, -1]

    print(X)
    print(y)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=0)
    lm = sklearn.linear_model.LinearRegression(fit_intercept=False)
    lm.fit(X_train, y_train)
    # print(lm.score(X_test, y_test))
    y_predicted = lm.predict(X_test)
    print(f'Person for {filter}: {pearsonr(y_test, y_predicted)[0]}')

    smm = sm.OLS(y_train, X_train)
    res = smm.fit()
    print(res.summary())
    print('#################################################')


def solve(l, r):
    pass


STEP = 400000
SLICE = STEP
MAX = 6000000


def compute_pixel(nside, l, b):
    return hp.ang2pix(nside, numpy.radians(-b + 90), numpy.radians(l))


def main():
    start = time.time()
    dataset = pd.read_csv('full.tsv', sep='\t')
    read_ts = time.time()
    print(f'Read finished in {int(read_ts - start)} s')
    dataset = dataset.sort_values(by='parallax', ascending=False)
    sort_ts = time.time()
    print(f'Sort finished in {int(sort_ts - read_ts)} s')

    dists = []
    coefs = {}
    errors = {}
    pearsons = []

    def pizza(df):
        m_dist = distance(df.iloc[len(df) // 2])
        dists.append(m_dist)

        l0 = radians(120)
        R0 = 0.1500

        X, y = prepare_bottlinger(df, l0, R0)

        smm = sm.OLS(y, X)
        st = time.time()
        res = smm.fit()
        # print(f'Model fitting finished in {time.time() - st} s')
        # print(f'Distance from {int(l_dist)}, median {int(m_dist)}, to {int(r_dist)}')
        # print(res.summary())
        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            key = res.params.keys()[i]
            if key not in coefs:
                coefs[key] = []
            if key not in errors:
                errors[key] = []
            coefs[key].append(val)
            errors[key].append(err)
            print(f'{key}={val}±{err}')

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
        lm = sklearn.linear_model.LinearRegression(fit_intercept=False)
        lm.fit(X_train, y_train)
        y_predicted = lm.predict(X_test)
        pearson = pearsonr(y_test, y_predicted)[0]
        pearsons.append(pearson)
        print(f'Pearson: {pearson}')
        print()

    slices(pizza, dataset, MAX, STEP, SLICE)

    spline_dists = numpy.linspace(min(dists), max(dists), 300)

    # figure_lines = [
    #     ['U', 'V', 'W'],
    #     ['Wx', 'Wy', 'Wz'],
    #     ['M12', 'M13', 'M23'],
    #     ['M11', 'M22', 'M33']
    # ]
    # figure_lines = [
    #     ['U', 'V', 'W'],
    #     ['M12', 'M13', 'M23'],
    #     ['M11', 'M22', 'M33']
    # ]
    # figure_lines = [
    #     ['U', 'V', 'W'],
    #     ['Wx', 'Wy', 'Wz'],
    #     ['M12', 'M13', 'M23'],
    #     ['C', 'K']
    # ]
    figure_lines = [
        ['UG', 'VG', 'WG'],
        ['w0', 'w01', 'w011'],
        ['k0', 'k01', 'k011']
    ]

    for figure in figure_lines:
        for label in figure:
            print(label, end='\t')
            for val, err in zip(coefs[label], errors[label]):
                print('{0:.1f}\t{1:.1f}'.format(val, err).replace('.', ','), end='\t')
            print()

    print('Pearson', end='\t')

    for p in pearsons:
        print('{0:.2f}'.format(p).replace('.', ','), end='\t')
        print(end='\t')
    print()

    figures = []
    for i in range(len(figure_lines)):
        figures.append(plt.figure(i).add_subplot(111))

    colors = plt.cm.get_cmap('hsv', len(coefs))
    for i in range(len(figure_lines)):
        for label in figure_lines[i]:
            y = coefs[label]
            figures[i].plot(dists, y, 'ko')
            figures[i].plot(spline_dists, spline(dists, y, spline_dists), color=colors(hash(label) % len(coefs)),
                            label=label)

    for i in range(len(figures)):
        figures[i].set_xlabel('Distance, pc')
        figures[i].legend()

    last_indx = len(figure_lines)
    last_fig = plt.figure(last_indx).add_subplot(111)
    last_fig.plot(dists, pearsons, 'ko', dists, pearsons, 'k')
    last_fig.set_xlabel('Distance, pc')
    last_fig.set_xlabel('Model quantity')

    # plt.xlabel('Distance, pc')
    # plt.ylabel('U, V, W')
    # plt.xticks(range(len(coeff)), rotation='90')
    plt.show()

    # print(f'U={U}, V={V}, W={W}, Wx={Wx}, Wy={Wy}, Wz={Wz}, M12={M12}, M13={M13}, M11={M11}, M22={M22}, M33={M33}')


if __name__ == '__main__':
    main()
    # for r in range(0, 100, 100):
    #     compute(r, r + 100)
