from astropy.coordinates import SkyCoord, Distance, cartesian_to_spherical
from numpy import cos, sin
from scipy.stats import pearsonr

from common.gaia.with_rv import read_gaia_with_rv_full, slices
from common.linear_decomposer import compute_coeffs
from common.spherical import taj, saj, fvj, tdj, sdj
from rv.count.distribution import draw_distribution, draw_function
from scripts.ogorod_fg import columns, ct211, cs310, cv310, ogorod_fg_mub, ogorod_fg_mul, ogorod_fg_vr
import numpy as np
import pandas as pd
from astropy import units as u

def apply_f(l, b, px):

   # df = df[(df.parallax > 0)]
    c = SkyCoord(frame='galactic',
                 l=l * u.rad,
                 b=b * u.rad,
                 distance=Distance(parallax=px * u.mas, allow_negative=True))

    x = np.array(c.cartesian.x)
    y = np.array(c.cartesian.y)

    x /= 1000
    y /= 1000

    result_x = -x * y
    result_y = y ** 2

    _, lat, lon = cartesian_to_spherical(result_x, result_y, 0)


    # g = c.transform_to(Galactocentric)
    return lat, lon


def main():
    df = read_gaia_with_rv_full()

    df = df.sample(100000)

    l = df.l
    b = df.b
    mul = df['mul']
    mub = df.mub
    r = 1 / df.px
    vr = df.radial_velocity

    def f(df, l_dist, r_dist):
        _, _, dmul = compute_coeffs(df, ogorod_fg_mul)
        _, _, dmub = compute_coeffs(df, ogorod_fg_mub)
        _, _, dvr = compute_coeffs(df, ogorod_fg_vr)

        t211l = taj(6, df.l, df.b)
        s310l = saj(10, df.l, df.b)
        t211b = tdj(6, df.l, df.b)
        s310b = sdj(10, df.l, df.b)
        v310 = fvj(10, df.l, df.b)

        test_x, test_y = apply_f(df.l.values, df.b.values, df.px.values)



        lr = f'{int(l_dist)}-{int(r_dist)}'

        draw_function(df.l, df.b, dmul, 25, "dmul", f"{lr}_dmul")
        draw_function(df.l, df.b, s310l, 1, "dmul", f"{lr}_s310")
        draw_function(df.l, df.b, t211l, 1, "dmul", f"{lr}_t211")

        draw_function(df.l, df.b, dmub, 25, "dmub", f"{lr}_dmub")
        draw_function(df.l, df.b, s310b, 1, "dmub", f"{lr}_s310")
        draw_function(df.l, df.b, t211b, 1, "dmub", f"{lr}_t211")

        draw_function(df.l, df.b, dvr, 25, "dvr", f"{lr}_dvr")
        draw_function(df.l, df.b, v310, 1, "dvr", f"{lr}_v310")

        draw_function(df.l, df.b, test_x, 1000, "test_x", f"{lr}_test_x")
        draw_function(df.l, df.b, test_y, 1000, "test_y", f"{lr}_test_y")

    slices(f, df)


    return

    # pizza = (dmul, dmub)




    t211_l=*(ct211 * cos(2 * b) * cos(l)),
    s310_l=*(cs310 * (-(5 * sin(b) ** 2 - 1) * cos(l))),

    t211_b=*(ct211 * sin(b) * sin(l)),
    s310_b=*(cs310 * ((15 * sin(b) ** 2 - 11) * sin(b) * sin(l))),

    v310=*(cv310 * (5 * sin(b) ** 2 - 1) * cos(b) * sin(l)),

    dmu = np.array([dmul, dmub]).T
    t211 = np.array([t211_l, t211_b]).T
    s310 = np.array([s310_l, s310_b]).T

    # pizza = np.corrcoef(dmu, t211)
    # print(pizza)

    df['t211_corr'] = corr(dmu, t211)
    df['s310_corr'] = corr(dmu, s310)

    # t211_best = df[(df.t211_corr > 0.9)]
    # print(len(t211_best))
    kk = [-1, -0.99, -0.9, -0.8, -0.5, 0.5, 0.8, 0.9, 0.99, 1]

    for i, k1 in enumerate(kk):
        if i == 0 or i == len(kk) - 1:
            continue
        k2 = kk[i + 1]

        def f1(slice, l_dist, r_dist):
            draw_distribution(slice, f"t211_{k1}_{k2}", f"Звезды с корелляцией с t211 {k1} < k < {k2} от {l_dist} до {r_dist}")

        slices(f1, df, filter=lambda df: (k1 < df.t211_corr) & (df.t211_corr < k2))

        def f2(slice, l_dist, r_dist):
            draw_distribution(slice, f"s310_{k1}_{k2}", f"Звезды с корелляцией с s310 {k1} < k < {k2} от {l_dist} до {r_dist}")

        slices(f2, df, filter=lambda df: (k1 < df.s310_corr) & (df.s310_corr < k2))

    # pizza = corr2_coeff(np.array([dmul, dmub]), np.array([t211_l, t211_b]))
    # print(pizza)


def corr(flow1, flow2):
    # flow1 = np.random.uniform(size=(10, 10, 2))  # the 3rd dimension is for the components
    # flow2 = flow1 + np.random.uniform(size=(10, 10, 2))
    flow1_centered = flow1 - np.mean(flow1, axis=(0))
    flow2_centered = flow2 - np.mean(flow2, axis=(0))
    inner_product = np.sum(flow1_centered * flow2_centered, axis=(1))
    r = inner_product / np.sqrt(np.sum(flow1_centered ** 2, axis=(1)) * np.sum(flow2_centered ** 2, axis=(1)))
    return r



if __name__ == '__main__':
    main()