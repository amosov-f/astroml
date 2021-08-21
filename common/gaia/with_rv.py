from pathlib import Path

import pandas as pd
from numpy import radians, degrees, average, sin, cos, arcsin, arctan2

from common.astro import galaxy, galaxy_mu


def to_galaxy(df, with_pm_errors=True):
    a, d, mua, mud, px = radians(df.ra), radians(df.dec), df.pmra, df.pmdec, df.parallax
    l, b = galaxy(a, d)
    mul, mub = galaxy_mu(mua, mud, l, b, d)
    if with_pm_errors:
        mul_error, mub_error = galaxy_mu(df.pmra_error, df.pmdec_error, l, b, d)
        return l, b, mul, mub, px, mul_error, mub_error
    return l, b, mul, mub, px


def prepare_galaxy_dataset(raw_dataset, with_pm_errors=True):
    if with_pm_errors:
        l, b, mul, mub, px, mul_error, mub_error = to_galaxy(raw_dataset, True)
        return pd.DataFrame(data={
            'l': l, 'b': b, 'mul': mul, 'mub': mub, 'px': px, 'radial_velocity': raw_dataset.radial_velocity,
            'mul_error': mul_error, 'mub_error': mub_error,
            'ra': raw_dataset.ra, 'dec': raw_dataset.dec, 'pmra': raw_dataset.pmra, 'pmdec': raw_dataset.pmdec,
            'parallax': raw_dataset.parallax
        })
    else:
        l, b, mul, mub, px = to_galaxy(raw_dataset, False)
        return pd.DataFrame(data={
            'l': l, 'b': b, 'mul': mul, 'mub': mub, 'px': px, 'radial_velocity': raw_dataset.radial_velocity
        })

def distance(row):
    return 1000 / row['px']


def slices(func, dataset, imax, step, slice, prepare=False):
    if prepare:
        dataset = prepare_galaxy_dataset(dataset, with_pm_errors=False)
    dataset = dataset.sort_values(by='px', ascending=False)
    average_dists = []
    for l in range(0, min(imax, len(dataset)), step):
        r = min(l + slice, len(dataset) - 1)
        print(f'Звезды под номерами от {l} до {r}')
        l_dist = distance(dataset.iloc[l])
        r_dist = distance(dataset.iloc[r])
        m_dist = distance(dataset.iloc[(l + r) // 2])
        average_dist = average(distance(dataset.iloc[l:r]))
        average_dists.append(str(int(average_dist)))
        print(f'На расстоянии от {int(1000 * l_dist)}пк, медиана {int(1000 * m_dist)}пк, среднее {int(1000 * average_dist)}пк до {int(1000 * r_dist)}пк')
        func(dataset.iloc[l:r], str(int(average_dist)))
    print('\t\t'.join(average_dists))


def read_raw_gaia_with_rv():
    file_dir = Path(__file__).parent
    return pd.read_csv(file_dir.joinpath('full_with_errors.tsv'), sep='\t')


def read_gaia_with_rv():
    df = read_raw_gaia_with_rv()
    df = df[df.parallax > 0.66]
    # df = df[df.parallax_error / df.parallax < 0.5]
    return prepare_galaxy_dataset(df)
