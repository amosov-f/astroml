from astropy.coordinates import Galactocentric

from cartesian.common_cartesian import read_galactic, read_radec, read_gaia_with_rv_cartesian
from common.gaia.with_rv import read_gaia_with_rv_xyz
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt


def main():
    df = read_gaia_with_rv_cartesian()

    df = df[((df.y < 0.010) & (df.z < 0.010))]

    bin_means, bin_edges, binnumber = stats.binned_statistic(df.x, df.vy, bins=200, range=[-8.0, 8.0],
                                                             statistic='mean')

    x = (bin_edges[:-1] + bin_edges[1:]) / 2

    y = -30.753 -11.174 * x

    y2 = -21.644 -4.774*x -1.811*x**2

    y3 = -20.51 -4.45 * x	+ -1.90 * x**2

    plt.stairs(bin_means, bin_edges, baseline=None, label='Усредненные значения vy звезд')

    plt.plot(x, y, label='Модель Огордникова-Милна (табл. 1, 2)', linestyle='dashed')

    plt.plot(x, y2, label='Квадратичная модель (табл. 5)', linestyle='dashdot')

    plt.plot(x, y3, label='Упрощенная квадратичная модель (табл. 6)', linestyle='dashdot')

    plt.legend()

    plt.xlabel('x, кпк')
    plt.ylabel('vy, км/с')
    plt.title('Скорость движения звезд вдоль оси y\nдля звезд вблизи оси x')
    # plt.show()

    plt.savefig('galactic_vy.png')


if __name__ == '__main__':
    main()