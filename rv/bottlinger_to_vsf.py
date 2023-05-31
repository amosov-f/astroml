import pandas as pd
from common.spherical import sdj, tdj
from pandas import DataFrame
from common.healpix import pixel_centers
from common.pandas2 import split
import statsmodels.api as sm
from scripts.ogorod import bottlinger_mul, bottlinger_mub, K
from rv.vsf import prepare_ols

N = 25

def main():
    l0 = 0
    R0 = 1
    r = 0.5
    for key in ['UG', 'VG', 'WG', 'w01', 'w011', 'k0', 'k01', 'k011']:
        print(key)
        catalog = pixel_centers(16)

        catalog.px = 1 / r
        catalog['mul'] = 0
        catalog['mub'] = 0
        test = bottlinger_mul(catalog, l0, R0)
        if key in test:
            catalog['mul'] = test[key] / K / r
        catalog['mub'] = bottlinger_mub(catalog, l0, R0)[key] / K / r

        X, y = prepare_ols(catalog)

        res = sm.OLS(y, X).fit()
        for k, v in zip(res.params.keys(), res.params.values):
            if abs(v) > 0.01:
                print(f'{k}: {v}')
        print('-----------------------------------')




def main2():
    l0 = 0
    R0 = 1
    r = 1
    for key in ['UG', 'VG', 'WG', 'w01', 'w011', 'k0', 'k01', 'k011']:
        print(key)
        catalog = pixel_centers(16)

        catalog.px = 1 / r
        catalog['mul'] = 0
        catalog['mub'] = 0
        mulx = bottlinger_mul(catalog, l0, R0)
        mulx['l'] = catalog.l
        mulx['b'] = catalog.b
        mubx = bottlinger_mub(catalog, l0, R0)
        mubx['l'] = catalog.l
        mubx['b'] = catalog.b
        pizza = pd.concat([mulx, mubx]).fillna(0)
        df = DataFrame()
        for i in range(1, N + 1):
            df[f's_{i}'] = sdj(i, pizza['l'], pizza['b'])
            df[f't_{i}'] = tdj(i, pizza['l'], pizza['b'])
        df.y = pizza[key]
        # print(df)
        # print(df.y)
        X, y = split(df)
        # print(X)
        # print(y)
        res = sm.OLS(y, X).fit()
        for k, v in zip(res.params.keys(), res.params.values):
            if abs(v) > 0.01:
                print(f'{k}: {v}')
        print('-----------------------------------')
        # for k, v in res.params:
        #     print(k)
        #     print(v)res.params:
        #     print(k)
        #     print(v)
        # print(res.params)




if __name__ == '__main__':
    main()