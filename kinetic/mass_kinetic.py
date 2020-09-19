import time

import numpy
import pandas
from astropy.coordinates import SkyCoord, Distance
from astroquery.esasky import ESASky
from astroquery.gaia import Gaia
from astropy import units as u
from matplotlib.pyplot import hist
import matplotlib.pyplot as plt

from common.gaia.with_rv import read_gaia_with_rv, read_raw_gaia_with_rv, slices
from common.linear_decomposer import decompose
from scripts.ogorod_fg import compute_ogorod_fg


def match(c1, c2):
    start = time.time()
    print('Matching started')
    idx, d2d, d3d = c1.match_to_catalog_sky(c2)
    print(f'Matching finished in {time.time() - start}s')

    # for i, index in idx.enumerate():
    #     catalog['']
    # numpy.save('idx', idx)
    return idx


def match2(c1, c2):
    return numpy.load('idx.npy')


def abs_vel2(v, sol_v):
    return (numpy.array(v.d_x) + sol_v[0]) ** 2 + (numpy.array(v.d_y) + sol_v[1]) ** 2 + (numpy.array(v.d_z) + sol_v[2]) ** 2


def main():
    # print(ESASky.list_catalogs())
    asu = pandas.read_csv('asu.tsv', sep='\t')
    print(asu)

    gaia_with_rv = read_gaia_with_rv()
    print(gaia_with_rv)

    # gaia_with_rv['px'] = gaia_with_rv['parallax']
    sort = [
        ['U', 'V', 'W'],
        ['A', 'B'],
        ['F', 'G', 'K']
    ]
    coeffs = decompose(gaia_with_rv, compute_ogorod_fg, 100500, 100500, 100500, sort)

    print(coeffs)

    U = coeffs['U'][0]
    V = coeffs['V'][0]
    W = coeffs['W'][0]

    idx = match2(asu, gaia_with_rv)



    asu['gaia_idx'] = idx

    print(asu)

    asu2 = asu.join(gaia_with_rv, on='gaia_idx')

    print(asu2)
    asu2['M'] = asu2.M1 + asu2.M2

    print(asu2)
    for _, row in asu2.iterrows():
        print(row)

    c = SkyCoord(frame="galactic", l=asu2.l.values * u.deg,
                       b=asu2.b.values * u.deg,
                       pm_l_cosb=asu2['mul'].values * u.mas / u.yr,
                       pm_b=asu2.mub.values * u.mas / u.yr,
                       radial_velocity=asu2.radial_velocity.values * u.km / u.s,
                       distance=Distance(parallax=asu2.parallax.values * u.mas, allow_negative=True))

    v2 = abs_vel2(c.velocity, (U, V, W))
    Mv2 = asu2.M.values * numpy.array(v2) / 2

    print(Mv2)

    # asu3 = asu2[(asu2['M2'] is None)]
    # print(asu3)
    plt.hist(Mv2, bins=100, range=(0, 20000))

    plt.show()

    # for row in asu2.iterrows():
    #     # gaia_idx = idx[i]
    #     print(row)
    #     # print(gaia_idx)
    #     print('---')


if __name__ == '__main__':
    main()
