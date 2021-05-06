import time

import numpy
import numpy as np
import pandas
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from common.gaia.with_rv import read_gaia_with_rv
from common.linear_decomposer import compute_coeffs
from rv.count.distribution import draw_distribution
from scripts.ogorod_fg import ogorod_fg_vr, K, compute_mu


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


def abs_vel2(v):
    return numpy.array(v.d_x) ** 2 + numpy.array(v.d_y) ** 2 + numpy.array(v.d_z) ** 2


def compute_absolute(gaia_with_rv):
    c = SkyCoord(frame="galactic",
                 l=gaia_with_rv.l.values * u.deg,
                 b=gaia_with_rv.b.values * u.deg,
                 pm_l_cosb=gaia_with_rv['dmul'].values * u.mas / u.yr,
                 pm_b=gaia_with_rv.dmub.values * u.mas / u.yr,
                 radial_velocity=gaia_with_rv.dvr.values * u.km / u.s,
                 distance=Distance(parallax=gaia_with_rv.px.values * u.mas))
    vx = numpy.array(c.velocity.d_x)
    vy = numpy.array(c.velocity.d_y)
    vz = numpy.array(c.velocity.d_z)
    v = (vx ** 2 + vy ** 2 + vz ** 2) ** 0.5

    vb = numpy.arcsin(- vz / v)
    vl = numpy.arctan2(vy, vx)
    return pandas.DataFrame(dict(l=gaia_with_rv.l, b=gaia_with_rv.b, vb=vb, vl=vl, v=v, dist=numpy.array(c.distance)))

def main():
    gaia_with_rv = read_gaia_with_rv()
    print(gaia_with_rv)

    print(np.min(gaia_with_rv.parallax))
    print(np.max(gaia_with_rv.parallax))

    # gaia_with_rv['px'] = gaia_with_rv['parallax']
    sort = [
        ['U', 'V', 'W'],
        ['A', 'B'],
        ['F', 'G', 'K']
    ]

    _, _, diffs = compute_coeffs(gaia_with_rv, compute_mu)
    # _, _, dmub = compute_coeffs(gaia_with_rv, ogorod_fg_mub)
    _, _, dvr = compute_coeffs(gaia_with_rv, ogorod_fg_vr)

    total = len(diffs)
    print(total)

    gaia_with_rv['dmul'] = diffs.iloc[0:total // 2] / K
    gaia_with_rv['dmub'] = diffs.iloc[total // 2:total] / K
    gaia_with_rv['dvr'] = dvr / gaia_with_rv.parallax#diffs.iloc[2 * total // 3:total]

    print('pizza')

    v = compute_absolute(gaia_with_rv)

    points = numpy.array([v.vl.values, v.vb.values, v.v.values, v.dist.values])
    points = np.transpose(points)

    delta_l = np.deg2rad(3)
    delta_b = np.deg2rad(3)
    delta_v = 10
    delta_d = 10
    max_v = 200
    max_d = 1500

    bins_l = int(2 * np.pi / delta_l)
    bins_b = int(np.pi / delta_b)
    bins_v = int(max_v / delta_v)
    bins_d = int(max_d / delta_d)
    sh = (bins_l, bins_b, bins_v, bins_d)

    H, (le, be, ve, de) = np.histogramdd(points,
                                         bins=sh,
                                         range=[[-np.pi, np.pi], [-np.pi / 2, np.pi / 2], [0, max_v], [0, max_d]])

    totals = np.maximum(H.sum(axis=0, keepdims=1).sum(axis=1, keepdims=1), 1000)

    percents = H / totals

    percents = np.nan_to_num(percents)

    check = percents.sum(axis=0).sum(axis=0)


    indexes = np.dstack(np.unravel_index(np.argsort(percents.ravel()), sh))[0]

    for il, ib, iv, id in list(reversed(indexes))[0:30]:
        vl, vr = ve[iv], ve[iv + 1]
        dl, dr = de[id], de[id + 1]
        print(f'l={le[il]}-{le[il + 1]}')
        print(f'b={be[ib]}-{be[ib + 1]}')
        print(f've={vl}-{vr}')
        print(f'de={dl}-{dr}')

        # stepd = 50
        # stepv = 25

        # mindl = 150
        # maxdl = 200
        # minv = 25
        # maxv = 50

    # for dl in range(mindl, maxdl, stepd):
    #     dr = dl + stepd
        v_fd = v[(dl < v.dist) & (v.dist < dr)]
        # for vl in range(minv, maxv, stepv):
        # vr = vl + stepv
        v_fv = v_fd[(vl < v_fd.v) & (v_fd.v < vr)]

        cl = 154.95
        cb = -35.34
        rad = 5
        center = SkyCoord(frame="galactic", l=cl * u.deg, b=cb * u.deg)
        def in_circle(x):
            c = SkyCoord(frame="galactic", l=x.vl * u.rad, b=x.vb * u.rad)
            sep = center.separation(c)
            return sep.degree < rad

        # v_fv = v_fv[(in_circle(v_fv))]

        draw_distribution(pandas.DataFrame(dict(b=numpy.array(v_fv.vb), l=numpy.array(v_fv.vl))),
                          f'Направления остаточных скоростей для {dl}<d<{dr} {vl}<v<{vr}')

        # draw_distribution(pandas.DataFrame(dict(b=numpy.array(v_fv.b), l=numpy.array(v_fv.l))),
        #                   f'Расположения звезд для сгустка l={cl} b={cb} radius={rad} {dl}<d<{dr} {vl}<v<{vr}')

    # vb = gaia_with_rv.mub
    # vl = gaia_with_rv['mul']

    # vb = vb[~numpy.isnan(vb)]
    # vl = vl[~numpy.isnan(vl)]

    # draw_distribution(pandas.DataFrame(dict(b=gaia_with_rv.b.values, l=gaia_with_rv.l.values)))
    # draw_distribution(pandas.DataFrame(dict(b=vb.values, l=vl.values)))

    # print(c)


if __name__ == '__main__':
    main()
