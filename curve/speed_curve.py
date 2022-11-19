from astropy.coordinates import spherical_to_cartesian

from common.gaia.with_rv import read_gaia_with_rv_full, read_gaia_with_rv_xyz
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def main():
    # print(spherical_to_cartesian(100, lat=0, lon=0))
    # return
    df = read_gaia_with_rv_xyz()

    min_x = -16000
    max_x = 0
    bin_size = 200
    bins_x = np.linspace(min_x, max_x, ((max_x - min_x) // bin_size) + 1)
    min_y = - (max_x - min_x) // 2
    max_y = + (max_x - min_x) // 2
    bins_y = np.linspace(min_y, max_y, ((max_y - min_y) // bin_size) + 1)

    x_means, x_edge, y_edge, binnumber = stats.binned_statistic_2d(df.x, df.y, df.vx, statistic='median', bins=(bins_x, bins_y))
    # x_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(df.x, df.y, df.vx, statistic='std', bins=(bins_x, bins_y))
    y_means, _, _, _ = stats.binned_statistic_2d(df.x, df.y, df.vy, statistic='median', bins=(bins_x, bins_y))
    # y_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(df.x, df.y, df.vy, statistic='std', bins=(bins_x, bins_y))
    z_means, _, _, _ = stats.binned_statistic_2d(df.x, df.y, df.vz, statistic='median', bins=(bins_x, bins_y))
    # z_std, x_edge, y_edge, binnumber = stats.binned_statistic_2d(df.x, df.y, df.vz, statistic='std', bins=(bins_x, bins_y))
    # count, _, _, _ = stats.binned_statistic_2d(df.x, df.y, df.vz, statistic='count', bins=(bins_x, bins_y))

    C = z_means

    x_averages = (x_edge[:-1] + x_edge[1:]) / 2
    y_averages = (y_edge[:-1] + y_edge[1:]) / 2

    # bin_means, bin_edges, binnumber = stats.binned_statistic(df.x, df.vx, bins=bins_x, statistic='mean')
    # plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='r', lw=3, label='binned statistic of data')

    # plt.hist2d(df.x, df.y, bins=(bins_x, bins_y), cmap='Blues')
    # plt.show()

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

    fig.set_figwidth(50)  # ширина и
    fig.set_figheight(50)  # высота "Figure"

    plt.show()

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