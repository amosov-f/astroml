from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance, Galactocentric
from numpy import radians, average

from common.astro import galaxy, galaxy_mu


def to_galaxy(df, with_pm_errors=True):
    a, d, mua, mud, px = radians(df.ra), radians(df.dec), df.pmra, df.pmdec, df.parallax
    l, b = galaxy(a, d)
    mul, mub = galaxy_mu(mua, mud, l, b, d)
    if with_pm_errors:
        mul_error, mub_error = galaxy_mu(df.pmra_error, df.pmdec_error, l, b, d)
        return l, b, mul, mub, px, mul_error, mub_error
    return l, b, mul, mub, px


def prepare_no_error_galaxy_dataset(raw_dataset):
    l, b, mul, mub, px = to_galaxy(raw_dataset, with_pm_errors=False)
    return pd.DataFrame(data={
        'l': l,
        'b': b,
        'mul': mul,
        'mub': mub,
        'px': px,
        'radial_velocity': raw_dataset.radial_velocity
    })


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
        return prepare_no_error_galaxy_dataset(raw_dataset)

def distance(row):
    return 1000 / row['px']


def slices(func, dataset, imax=6000000, step=400000, slice=400000, prepare=False, filter=None):
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
        df = dataset.iloc[l:r]
        if filter:
            df = df[(filter(df))]
        print(f'На расстоянии от {int(l_dist)}пк, медиана {int(m_dist)}пк, среднее {int(average_dist)}пк до {int(r_dist)}пк после фильтра звезд {len(df)}')
        func(df, l_dist, r_dist) #, str(int(average_dist)))
    return average_dists


def read_raw_gaia_with_rv_no_errors():
    file_dir = Path(__file__).parent
    return pd.read_csv(file_dir.joinpath('gdr3_short.csv'))


def read_raw_gaia_with_rv():
    file_dir = Path(__file__).parent
    return pd.read_csv(file_dir.joinpath('full_with_errors.tsv'), sep='\t')


def read_gaia_with_rv_1500():
    df = read_raw_gaia_with_rv()
    df = df[df.parallax > 0.66]
    # df = df[df.parallax_error / df.parallax < 0.5]
    return prepare_galaxy_dataset(df)


def read_gaia_with_rv_full():
    df = read_raw_gaia_with_rv_no_errors()
    return prepare_galaxy_dataset(df, with_pm_errors=False)


def read_gaia_with_rv_xyz():
    df = read_raw_gaia_with_rv()
    # df = df.iloc[0:10]
    df = df[(df.parallax > 0)]
    c = SkyCoord(ra=df.ra.values * u.deg,
                 dec=df.dec.values * u.deg,
                 pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
                 pm_dec=df.pmdec.values * u.mas / u.yr,
                 radial_velocity=df.radial_velocity.values * u.km / u.s,
                 distance=Distance(parallax=df.parallax.values * u.mas))
    g = c.transform_to(Galactocentric)
    return pd.DataFrame(dict(x=np.array(g.cartesian.x),
                             y=np.array(g.cartesian.y),
                             z=np.array(g.cartesian.z),
                             vx=np.array(g.velocity.d_x),
                             vy=np.array(g.velocity.d_y),
                             vz=np.array(g.velocity.d_z)))