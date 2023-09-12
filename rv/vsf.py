import time

import pandas as pd
import statsmodels.api as sm

from common.gaia.with_rv import read_gaia_with_rv_full, SAMPLE, GDR2
from common.gaia.with_rv import slices
from common.pandas2 import split
from common.spherical import taj, saj, tdj, sdj, indexes
from common.spherical_decomposer import decompose_spherical

N = 24
KK = 4.738

def mul_apply(df):
    d = {}
    for j in range(1, N + 1):
        d[f't_{j}'] = taj(j, df.l, df.b)
        d[f's_{j}'] = saj(j, df.l, df.b)
    # d = OrderedDict(sorted(d.items()))
    d['y'] = df['mul'] * KK
    return pd.DataFrame(d)

def mub_apply(df):
    d = {}
    for j in range(1, N + 1):
        d[f't_{j}'] = tdj(j, df.l, df.b)
        d[f's_{j}'] = sdj(j, df.l, df.b)
    # d = OrderedDict(sorted(d.items()))
    d['y'] = df.mub * KK
    return pd.DataFrame(d)


def prepare_ols(dataset):
    mul = mul_apply(dataset)
    mub = mub_apply(dataset)
    return split(pd.concat([mul, mub]))


def to_map(model):
    result = {}
    for i in range(len(model.params)):
        val = model.params[i]
        err = model.bse[i]
        key = model.params.keys()[i]
        result[key] = (val, err)
    return result

STEP = 2000000
SLICE = STEP
MAX = STEP * 15

def main():
    start = time.time()
    # dataset = read_gc2019()#read_gaia_with_rv_in_galactocentric()#read_gc2019()#read_gaia_with_rv()
    dataset = read_gaia_with_rv_full()
    read_ts = time.time()
    print(f'Read finished in {int(read_ts - start)} s')
    # dataset = dataset.sample(min(500000, len(dataset)), random_state=0)
    sort_ts = time.time()
    print(f'Sort finished in {int(sort_ts - read_ts)} s')

    coefs = {}
    errors = {}
    distances = []
    l_dists = []
    r_dists = []

    table = {}

    def put_to_table(k, v1, v2):
        if k not in table:
            table[k] = []
        table[k].append(v1)
        table[k].append(v2)

    t211 = []
    s310 = []
    v310 = []

    def f(slice, dist_l, dist_r):
        X, y = prepare_ols(slice)
        # X2, y2 = prepare_mu_vr(slice)
        r = 1 / slice['px']

        slice.y = slice['radial_velocity'] / r
        rv_model = decompose_spherical(slice, 25)
        # coeffs, errors = show_spherical_decomposition(model, draw=False)


        smm = sm.OLS(y, X)
        # smm2 = sm.OLS(y2, X2)

        res = smm.fit()
        # res2 = smm2.fit()

        print(res.summary())
        vsf = to_map(res)
        # ogorod = to_map(res2)
        rv = to_map(rv_model)

        distances.append(dist_r)
        l_dists.append(dist_l)
        r_dists.append(dist_r)

        t211.append(f'{vsf["t_6"][0]:.1f}')
        s310.append(f'{vsf["s_10"][0]:.1f}')
        v310.append(f'{rv["v_10"][0]:.1f}')
        # t1 = vsf["t_1"][0]
        # t2 = vsf["t_2"][0]
        # t3 = vsf["t_3"][0]
        # s1 = vsf["s_1"][0]
        # s2 = vsf["s_2"][0]
        # s3 = vsf["s_3"][0]
        # s4 = vsf["s_4"][0]
        # s5 = vsf["s_5"][0]
        # s6 = vsf["s_6"][0]
        # s7 = vsf["s_7"][0]
        # s8 = vsf["s_8"][0]
        # w1 = ogorod["Wx"][0]
        # w2 = ogorod["Wy"][0]
        # w3 = ogorod["Wz"][0]
        # U = ogorod['U'][0]
        # V = ogorod['V'][0]
        # W = ogorod['W'][0]
        # C = ogorod['C'][0]
        # K = ogorod['K'][0]
        # M12 = ogorod['M12'][0]
        # M23 = ogorod['M23'][0]
        # M13 = ogorod['M13'][0]
        # v1 = rv['v_1'][0]
        # v2 = rv['v_2'][0]
        # v3 = rv['v_3'][0]
        # v4 = rv['v_4'][0]
        # v5 = rv['v_5'][0]
        # v6 = rv['v_6'][0]
        # v7 = rv['v_7'][0]
        # v8 = rv['v_8'][0]
        #
        # put_to_table('t1', t1, 2.894 * w3)
        # put_to_table('t2', t2, 2.894 * w2)
        # put_to_table('t3', t3, 2.894 * w1)
        # put_to_table('s1', s1, -2.894 * W / r)
        # put_to_table('s2', s2, -2.894 * V / r)
        # put_to_table('s3', s3, -2.894 * U / r)
        # put_to_table('s4', s4, -1.294 * K)
        # put_to_table('s5', s5, 2.242 * M23)
        # put_to_table('s6', s6, 2.242 * M13)
        # put_to_table('s7', s7, 2.242 * M12)
        # put_to_table('s8', s8, 2.242 * C)
        #
        # put_to_table('v1/r', v1 / r, -2.05 * W / r)
        # put_to_table('v2/r', v2 / r, -2.05 * V / r)
        # put_to_table('v3/r', v3 / r, -2.05 * U / r)
        # put_to_table('v4/r', v4 / r, -1.057 * K)
        # put_to_table('v5/r', v5 / r, 1.831 * M23)
        # put_to_table('v6/r', v6 / r, 1.831 * M13)
        # put_to_table('v7/r', v7 / r, 1.831 * M12)
        # put_to_table('v8/r', v8 / r, 1.83 * C)

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

        for i in range(len(rv_model.params)):
            val = rv_model.params[i]
            key = rv_model.params.keys()[i]
            if key not in coefs:
                coefs[key] = []
            coefs[key].append(val)

    slices(f, dataset, MAX, STEP, SLICE, prepare=False)

    keys = list(table.keys())
    # keys = keys[0:len(keys):2] + keys[1:len(keys):2]
    for k in keys:
        print(f'{k}', end='\t')
        for val in table[k]:
            print(f'{val:.1f}'.replace('.', ','), end='\t')
        print()

    print('min', end='\t')
    for d in l_dists:
        print(int(d), end='\t')
    print()
    print('max', end='\t')
    for d in r_dists:
        print(int(d), end='\t')
    print()
    keys = list(coefs.keys())
    keys = keys[0:2 * len(keys) // 3:2] + keys[1:2 * len(keys) // 3:2] + keys[2 * len(keys) // 3:len(keys)]
    for key in keys:
        j = int(key[2:])
        n, k, l = indexes(j)
        print(f'{key[0:2]}{n}{k}{l}', end='\t')
        # for val, err in zip(coefs[k], errors[k]):
        #     print(f'{val:.1f}\t{err:.1f}'.replace('.', ','), end='\t')
        for val in coefs[key]:
            print(f'{val:.1f}'.replace('.', ','), end='\t')
        print()

    # print(f't211 {t211}')
    # print(f's310 {s310}')
    # print(f'v310 {v310}')


if __name__ == '__main__':
    main()