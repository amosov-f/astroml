import os
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np

from common.gaia.with_rv import read_gaia_dr3_with_rv_with_errors
from gaia.yt.gdr3.filter_with_error import read_gdr3


def main():
    df = read_gaia_dr3_with_rv_with_errors()

    Path('dr3_with_rv').mkdir(exist_ok=True)

    plt.rcParams["figure.figsize"] = (6.4*2,4.8*2)
    plt.rcParams.update({'font.size': 22})

    #
    # print(df[(df.parallax_error < 1)].shape[0] / df.shape[0])
    # print(df[(abs(df.parallax_error / df.parallax) < 0.1)].shape[0] / df.shape[0])

    px = 1 / df['parallax']
    plt.hist(px, bins=100, range=(0, 15))
    plt.xlabel('Расстояние, кпк')
    plt.ylabel('Число звезд')
    plt.title('Распределение расстояний\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/distance_dr3.png')

    plt.clf()

    plt.hist(df['parallax'], bins=100, range=(-0.3, 8))
    plt.xlabel('Параллакс, mas')
    plt.ylabel('Число звезд')
    plt.title('Распределение параллаксов\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/parallax_dr3.png')

    print(df[(df.parallax < 0)].shape[0] / df.shape[0])

    plt.clf()

    dv = df['radial_velocity']
    plt.hist(dv, bins=100, range=(-200, 200))
    plt.xlabel('Лучевая скорость, км/с')
    plt.ylabel('Число звезд')
    plt.title('Распределение лучевых скоростей\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/radial_velocity.png')

    plt.clf()

    dv = df['radial_velocity_error']
    plt.hist(dv, bins=100, range=(0, 20))
    plt.xlabel('Ошибка лучевой скорости, км/с')
    plt.ylabel('Число звезд')
    plt.title('Распределение ошибки лучевых скоростей\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/radial_velocity_error.png')

    plt.clf()

    dv = np.abs(100 * df['radial_velocity_error'] / df['radial_velocity'])
    plt.hist(dv, bins=100, range=(0, 100))
    plt.xlabel('Относительная ошибка лучевой скорости, %')
    plt.ylabel('Число звезд')
    plt.title('Распределение относительной ошибки лучевых скоростей\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/radial_velocity_error_sigma.png')

    plt.clf()

    px_error = np.abs(100 * df['pmra_error'] / df['pmra'])
    plt.hist(px_error, bins=100, range=(0, 100))
    plt.xlabel('Относительная точность собственного движения (прямое восхождение), %')
    plt.ylabel('Число звезд')
    plt.title('Распределение относительной точности собственного движения\nзвезд GAIA DR3 with RV (прямое восхождение)')
    plt.savefig('dr3_with_rv/pmra_error_sigma.png')

    plt.clf()

    px_error = (df['pmra_error'] ** 2 + df['pmdec_error'] ** 2) ** 0.5
    plt.hist(px_error, bins=100, range=(0, 0.12))
    plt.xlabel('Точность собственных движений, mas/yr')
    plt.ylabel('Число звезд')
    plt.title('Распределение точности полного собственного движения\nзвезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/pmra_pmdec_error.png')

    plt.clf()

    px_error = 100 * ((df['pmra_error'] ** 2 + df['pmdec_error'] ** 2) ** 0.5) / ((df['pmra'] ** 2 + df['pmdec'] ** 2) ** 0.5)
    plt.hist(px_error, bins=100, range=(0, 5))
    plt.xlabel('Относительная точность собственных движений, %')
    plt.ylabel('Число звезд')
    plt.title('Распределение относительной точности полного собственного\nдвижения звезд GAIA DR3 with RV')
    plt.savefig('dr3_with_rv/pmra_pmdec_error_sigma.png')

    plt.clf()

    px_error = np.abs(100 * df['pmdec_error'] / df['pmdec'])
    plt.hist(px_error, bins=100, range=(0, 100))
    plt.xlabel('Относительная точность собственного движения (склонение), %')
    plt.ylabel('Число звезд')
    plt.title('Распределение относительной точности собственного движения\nзвезд GAIA DR3 with RV (склонение)')
    plt.savefig('dr3_with_rv/pmdec_error_sigma.png')

    plt.clf()

    df = read_gdr3()

    plt.hist(df['phot_g_mean_mag'], bins=15, range=(5, 20))
    plt.xlabel('Звездная величина G, mag')
    plt.ylabel('Число звезд')
    plt.title('Распределение звездной величины G')
    plt.savefig('dr3_with_rv/phot_g_mean_mag.png')

    plt.clf()

    print(np.log10(0.49))

    print(np.log10(df['parallax']))

    print(df['parallax'])

    plt.hist(df['phot_g_mean_mag'] + 5 + np.log10(df['parallax']), bins=15, range=(10, 25))
    plt.xlabel('Абсолютная звездная величина G, mag')
    plt.ylabel('Число звезд')
    plt.title('Распределение абсолютной звездной величины G')
    plt.savefig('dr3_with_rv/phot_g_abs_mag.png')

if __name__ == '__main__':
    main()