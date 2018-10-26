import time

import numpy
import sklearn

import pandas as pd
from scipy.stats import pearsonr
import statsmodels.api as sm
import math
import sklearn.linear_model as lm
from sklearn.model_selection import cross_validate
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from scipy.interpolate import spline
import healpy as hp

Leo = 282.85948083 # °
L0 = 32.931918056 # °
si = 0.88998807641_8 # sin62.871748611°
ci = 0.45598379779_8 # cos62.871748611°

Pi = 3.1415926535897932384626433832795_8
K = 4.74

def sind(x):
    return math.sin(math.radians(x))

def cosd(x):
    return math.cos(math.radians(x))

def asind(x):
    return math.degrees(math.asin(x))

def acosd(x):
    return math.degrees(math.acos(x))

def atan2d(y, x):
    return math.degrees(math.atan2(y, x))

def galaxy(a, d):
    al=a-Leo
    sa=sind(al)
    ca=cosd(al)
    sd=sind(d)
    cd=cosd(d)
    b = asind(sd * ci - cd * si * sa)
    l = atan2d(sd * si + cd * ci * sa, cd * ca) + L0
    if l < 0:
        l += 360.0
    return l, b

def galaxy_mu(mua, mud, l, b, d):
    cd = cosd(d)
    sfi = si * cosd(l - L0) / cd
    cfi = (cosd(b) * ci - sind(b) * si * sind(l - L0)) / cd
    mul = cfi * mua + sfi * mud
    mub = - sfi * mua + cfi * mud
    return mul, mub

def galaxy_apply(row):
    l, b = galaxy(row['ra'], row['dec'])
    return pd.Series([l, b], index=['l', 'b'])

def galaxy_mu_apply(row):
    mul, mub = galaxy_mu(row['pmra'], row['pmdec'], row['l'], row['b'], row['dec'])
    # mul, mub = row['pmra'], row['pmdec']
    return (mul, mub)


def kmul_apply(row):
    l, b, mul, mub, px = params(row)
    return dict(U=+px * sind(l),
                V=-px * cosd(l),
                W=0,
                Wx=-sind(b) * cosd(l),
                Wy=-sind(b) * sind(l),
                Wz=+cosd(b),
                M12=+cosd(b) * cosd(2 * l),
                M13=-sind(b) * sind(l),
                M23=+sind(b) * cosd(l),
                M11=-0.5 * cosd(b) * sind(2 * l),
                M22=+0.5 * cosd(b) * sind(2 * l),
                M33=0.0,
                y=K * mul)


def kmul_prim_apply(row):
    l, b, mul, mub, px = params(row)
    return series(U=+px * sind(l),
                  V=-px * cosd(l),
                  W=0.0,
                  Wx=-sind(b) * cosd(l),
                  Wy=-sind(b) * sind(l),
                  Wz=+cosd(b),
                  M12=+cosd(b) * cosd(2 * l),
                  M13=-sind(b) * sind(l),
                  M23=+sind(b) * cosd(l),
                  M11z=-0.5 * cosd(b) * sind(2 * l),
                  X=0.0,
                  y=K * mul)


def kmub_apply(row):
    l, b, mul, mub, px = params(row)
    return series(U=+px * cosd(l) * sind(b),
                  V=+px * sind(l) * sind(b),
                  W=-px * cosd(b),
                  Wx=sind(l),
                  Wy=-cosd(l),
                  Wz=0.0,
                  M12=-0.5 * sind(2 * b) * sind(2 * l),
                  M13=+cosd(2 * b) * cosd(l),
                  M23=+cosd(2 * b) * sind(l),
                  M11=-0.5 * sind(2 * b) * cosd(l) ** 2,
                  M22=-0.5 * sind(2 * b) * sind(l) ** 2,
                  M33=+0.5 * sind(2 * b),
                  y=K * mub)


def kmub_prim_apply(row):
    l, b, mul, mub, px = params(row)
    return series(U=+px * cosd(l) * sind(b),
                  V=+px * sind(l) * sind(b),
                  W=-px * cosd(b),
                  Wx=sind(l),
                  Wy=-cosd(l),
                  Wz=0.0,
                  M12=-0.5 * sind(2 * b) * sind(2 * l),
                  M13=+cosd(2 * b) * cosd(l),
                  M23=+cosd(2 * b) * sind(l),
                  M11z=-0.25 * sind(2 * b) * cosd(2 * l),
                  X=+0.5 * sind(2 * b),
                  y=K * mub)

def vr_apply(row):
    l, b, mul, mub, px = params(row)
    x = {
        'U': -cosd(l) * cosd(b),
        'V': -sind(l) * cosd(b),
        "W": -sind(b),
        "M12": cosd(b) ** 2 * sind(2 * l) / px,
        'M13': sind(2 * b) * cosd(l) / px,
        "M23": sind(2 * b) * sind(l) / px,
        "M11": cosd(b) ** 2 * cosd(l) ** 2 / px,
        "M22": cosd(b) ** 2 * sind(l) ** 2 / px,
        "M33": sind(b) ** 2 / px,
        'y': row['radial_velocity']
    }
    return series(**x)


def params(row):
    a, d, mua, mud, px = row['ra'], row['dec'], row['pmra'], row['pmdec'], row['parallax']
    l, b = galaxy(a, d)
    mul, mub = galaxy_mu(mua, mud, l, b, d)
    return l, b, mul, mub, px

def prepare_mu(dataset):
    kmul = dataset.apply(kmul_prim_apply, axis=1)
    kmub = dataset.apply(kmub_prim_apply, axis=1)
    return split(pd.concat([kmul, kmub]))

def prepare_vr(dataset):
    return split(dataset.apply(vr_apply, axis=1))

def fast_apply(df: pd.DataFrame, func):
    new_columns = {}
    for _, row in df.iterrows():
        new_values = func(row)
        for k, v in new_values.items():
            if k not in new_columns:
                new_columns[k] = []
            new_columns[k].append(v)
    return pd.DataFrame(new_columns)

def prepare_mu_vr(dataset: pd.DataFrame):
    start = time.time()
    dataset.iterrows()

    kmul = fast_apply(dataset, kmul_apply)
    kmul_ts = time.time()
    # print(f'kmul computation finished in {int(kmul_ts - start)} s')
    kmub = fast_apply(dataset, kmub_apply)
    kmub_ts = time.time()
    # print(f'kmub computation finished in {int(kmub_ts - kmul_ts)} s')
    vr = fast_apply(dataset, vr_apply)
    vr_ts = time.time()
    # print(f'vr computation finished in {int(vr_ts - kmub_ts)} s')
    return split(pd.concat([kmul, kmub, vr]).fillna(0))


def split(dataset):
    X = dataset.iloc[:, :-1]
    y = dataset.iloc[:, -1]
    return X, y

def series(**kwargs):
    return kwargs

    # items = []
    # index = []
    # for k, v in kwargs.items():
    #     items.append(v)
    #     index.append(k)
    # return pd.Series(items, index=index)


def compute(l, r):
    filter = f'{l} <= 1000 / parallax < {r}'
    print(filter)
    df = pd.read_csv('full.tsv', sep='\t')\
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

def dist(row):
    return 1000 / row['parallax']

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

    for l in range(0, MAX, STEP):
        r = l + SLICE
        print(f'Stars with number from {l} to {r}')
        l_dist = dist(dataset.iloc[l])
        r_dist = dist(dataset.iloc[r])
        m_dist = dist(dataset.iloc[(l + r) // 2])
        dists.append(m_dist)

        nside = 32
        hist = numpy.zeros(hp.nside2npix(nside))
        for _, row in dataset.iloc[l:r].iterrows():
            l, b, mul, mub, px = params(row)
            pix = compute_pixel(nside, l, b)
            hist[pix] += 1

        plt.clf()

        hp.mollview(hist, title=f"Распределение звезд от {int(l_dist)} до {int(r_dist)} пк", unit='число звезд в пикселе', coord='G', min=0, max=500)
        hp.graticule()

        # plt.show()

        plt.savefig(f'count/count_{int(l_dist)}_{int(r_dist)}.png', dpi=150)

        # X, y = prepare_mu_vr(dataset.iloc[l:r])
        # smm = sm.OLS(y, X)
        # st = time.time()
        # res = smm.fit()
        # # print(f'Model fitting finished in {time.time() - st} s')
        # print(f'Distance from {int(l_dist)}, median {int(m_dist)}, to {int(r_dist)}')
        # # print(res.summary())
        # for i in range(len(res.params)):
        #     val = res.params[i]
        #     err = res.bse[i]
        #     key = res.params.keys()[i]
        #     if key not in coefs:
        #         coefs[key] = []
        #     if key not in errors:
        #         errors[key] = []
        #     coefs[key].append(val)
        #     errors[key].append(err)
        #     print(f'{key}={val}±{err}')

        # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
        # lm = sklearn.linear_model.LinearRegression(fit_intercept=False)
        # lm.fit(X_train, y_train)
        # y_predicted = lm.predict(X_test)
        # pearson = pearsonr(y_test, y_predicted)[0]
        # pearsons.append(pearson)
        # print(f'Pearson: {pearson}')
        # print()

    if True:
        return

    spline_dists = numpy.linspace(min(dists), max(dists), 300)

    figure_lines = [
        ['U', 'V', 'W'],
        ['Wx', 'Wy', 'Wz'],
        ['M12', 'M13', 'M23'],
        ['M11', 'M22', 'M33']
    ]

    for figure in figure_lines:
        for label in figure:
            print(label, end='\t')
            for val, err in zip(coefs[label], errors[label]):
                print('{0:.1f}\t{1:.1f}'.format(val, err).replace('.', ','), end='\t')
            print()
    print('Pearson', end='\t')

    for p in pearsons:
        print(int(100 * p), end='\t')
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
            figures[i].plot(spline_dists, spline(dists, y, spline_dists), color=colors(hash(label) % len(coefs)), label=label)

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