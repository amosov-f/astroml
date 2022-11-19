from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
import astropy.coordinates as coord
from pandas import DataFrame
from scipy import stats

from common.gaia.with_rv import read_raw_gaia_with_rv, read_gaia_with_rv_full
from common.linear_decomposer import compute_coeffs
from scripts.ogorod_fg import ogorod_fg_mul, ogorod_fg_mub, ogorod_fg_vr


def read_gaia_with_rv_residuals():
    df = read_gaia_with_rv_full()

    _, resl, dmul = compute_coeffs(df, ogorod_fg_mul)
    _, resb, dmub = compute_coeffs(df, ogorod_fg_mub)
    _, resvr, dvr = compute_coeffs(df, ogorod_fg_vr)

    c = SkyCoord(frame="galactic",
                 l=df.l.values * u.rad,
                 b=df.b.values * u.rad,
                 pm_l_cosb=dmul.values * u.mas / u.yr,
                 pm_b=dmub.values * u.mas / u.yr,
                 radial_velocity=dvr.values * u.km / u.s,
                 distance=Distance(parallax=df.px.values * u.mas, allow_negative=True))

    return c

def read_gaia_with_rv_redisuals_xyz():
    c = read_gaia_with_rv_residuals()

    c = c.transform_to(coord.Galactocentric)

   # df = df[(df.parallax > 0)]

    # g = c.transform_to(Galactocentric)
    return pd.DataFrame(dict(x=np.array(c.cartesian.x),
                             y=np.array(c.cartesian.y),
                             z=np.array(c.cartesian.z),
                             vx=np.array(c.velocity.d_x),
                             vy=np.array(c.velocity.d_y),
                             vz=np.array(c.velocity.d_z)))

def show_velocity(df, title, subtitle, order, min_x, max_x, min_y, max_y, min_z, max_z, bin_count):
    bins_x = np.linspace(min_x, max_x, bin_count + 1)
    bins_y = np.linspace(min_y, max_y, bin_count + 1)

    df = df[(df.x >= min_x) & (df.x <= max_x) & (df.y >= min_y) & (df.y <= max_y) & (df.z >= min_z) & (df.z <= max_z)]

    print(len(df))

    x_col = df.x
    y_col = df.y

    x_means, x_edge, y_edge, binnumber = stats.binned_statistic_2d(x_col, y_col, df.vx, statistic='mean',
                                                                   bins=(bins_x, bins_y))

    print("x done")
    # x_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(x_col, y_col, df.vx, statistic='std',
    #                                                              bins=(bins_x, bins_y))
    y_means, _, _, _ = stats.binned_statistic_2d(x_col, y_col, df.vy, statistic='mean', bins=(bins_x, bins_y))
    # y_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(x_col, y_col, df.vy, statistic='std',
    #                                                              bins=(bins_x, bins_y))
    print("y done")
    z_means, _, _, _ = stats.binned_statistic_2d(x_col, y_col, df.vz, statistic='mean', bins=(bins_x, bins_y))
    # z_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(x_col, y_col, df.vz, statistic='std',
    #                                                              bins=(bins_x, bins_y))
    print("z done")
    # count, _, _, _ = stats.binned_statistic_2d(x_col, y_col, df.vz, statistic='count', bins=(bins_x, bins_y))

    C = z_means

    x_averages = (x_edge[:-1] + x_edge[1:]) / 2
    y_averages = (y_edge[:-1] + y_edge[1:]) / 2


    X, Y = np.meshgrid(x_averages, y_averages)

    fig, ax = plt.subplots()

    print('min')
    print(np.min(C))
    print('max')
    print(np.max(C))

    np.savetxt('x_points.tsv', x_averages, fmt='%1.0f', delimiter=',')
    np.savetxt('y_points.tsv', y_averages, fmt='%1.0f', delimiter=',')
    np.savetxt('vx_median.tsv', x_means, fmt='%1.3f', delimiter=',')
    # np.savetxt('vx_std.tsv', x_std, fmt='%1.3f', delimiter=',')
    np.savetxt('vy_median.tsv', y_means, fmt='%1.3f', delimiter=',')
    # np.savetxt('vy_std.tsv', y_std, fmt='%1.3f', delimiter=',')
    np.savetxt('vz_median.tsv', z_means, fmt='%1.3f', delimiter=',')
    # np.savetxt('vz_std.tsv', z_std, fmt='%1.3f', delimiter=',')
    # np.savetxt('count.tsv', count, fmt='%1.0f', delimiter=',')

    ax.quiver(X, Y, x_means, y_means, C, width=0.003)

    width_x = 10

    fig.suptitle(title + '\n' + subtitle, fontsize=30)

    fig.set_figwidth(width_x)  # ширина и
    fig.set_figheight(width_x * (max_y - min_y) // (max_x - min_x))  # высота "Figure"

    plt.xlabel('X, кпк')
    plt.ylabel('Y, кпк')

    dir = f'fig/{title}'
    Path(dir).mkdir(parents=True, exist_ok=True)

    path = f'{dir}/{order}_{subtitle}.png'
    fig.savefig(path)

    plt.show()


def main():
    # print(spherical_to_cartesian(100, lat=0, lon=0))
    # return
    df = read_gaia_with_rv_redisuals_xyz()

    # bin_means, bin_edges, binnumber = stats.binned_statistic(df.x, df.vx, bins=bins_x, statistic='mean')
    # plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label='binned statistic of data')

    # plt.hist2d(df.x, df.y, bins=(bins_x, bins_y), cmap='Blues')
    # plt.show()
    min_x = -14000
    max_x = 14000
    min_y = -14000
    max_y = 14000

    show_velocity(df, "pizza", "pasta2", 0, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y, bin_count=10)

    # fig = df.x.plot.hist(bins=bins_x, figsize=(10, 5))
    # fig.set_xlabel("Расстояние [пк]")
    # fig.set_ylabel("Число звезд")
    # plt.savefig('count_to_glactic_center.png')

    return

    print(np.min(df.y))
    print(np.min(df.x))
    df = df[(np.abs(df.y) < 10)]
    print(df)
    # return



if __name__ == '__main__':
    main()