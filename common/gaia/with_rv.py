import time

from common.astro import galaxy, galaxy_mu
from common.pandas2 import fast_apply, series
import pandas as pd
from numpy import radians


def to_galaxy(df):
    a, d, mua, mud, px = radians(df.ra), radians(df.dec), df.pmra, df.pmdec, df.parallax
    l, b = galaxy(a, d)
    mul, mub = galaxy_mu(mua, mud, l, b, d)
    return l, b, mul, mub, px


def prepare_galaxy_dataset(raw_dataset):
    start = time.time()
    l, b, mul, mub, px = to_galaxy(raw_dataset)
    # raw_dataset['l', 'b', 'mul', 'mub', 'px'] = raw_dataset.apply(lambda row: to_galaxy(row), axis=1)
    print(time.time() - start)
    return pd.DataFrame(data={
        'l': l, 'b': b, 'mul': mul, 'mub': mub, 'px': px, 'radial_velocity': raw_dataset.radial_velocity
    })


def distance(row):
    return 1000 / row['px']


def slices(func, dataset, imax, step, slice):
    dataset = prepare_galaxy_dataset(dataset)
    dataset = dataset.sort_values(by='px', ascending=False)
    for l in range(0, imax, step):
        r = l + slice
        print(f'Stars with number from {l} to {r}')
        l_dist = distance(dataset.iloc[l])
        r_dist = distance(dataset.iloc[r])
        m_dist = distance(dataset.iloc[(l + r) // 2])
        print(f'Distance from {int(l_dist)}, median {int(m_dist)}, to {int(r_dist)}')
        func(dataset.iloc[l:r])
