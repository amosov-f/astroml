from common.gaia.with_rv import read_gaia_with_rv_1500, read_raw_gaia_with_rv
from matplotlib import pyplot as plt
import numpy as np

def main():
    df = read_raw_gaia_with_rv()

    dv = df['pmra_error']
    plt.hist(dv, bins=100, range=(0, 0.3))
    # plt.show()
    plt.xlabel('Ошибка, mas')
    plt.ylabel('Число звезд')
    plt.title('Распределение абсолютной ошибки собственных движений звезд')
    plt.savefig('pmra_error_abs.png')

    plt.clf()

    dv = df['radial_velocity_error']
    plt.hist(dv, bins=100, range=(0, 10))
    # plt.show()
    plt.xlabel('Ошибка лучевой скорости, км/с')
    plt.ylabel('Число звезд')
    plt.title('Распределение абсолютной ошибки лучевых скоростей звезд')
    plt.savefig('radial_velocity_error_abs.png')

    plt.clf()

    px = df['parallax']
    plt.hist(px, bins=100, range=(-0.1, 4))
    # plt.show()
    plt.xlabel('Параллаксы, mas')
    plt.ylabel('Число звезд')
    plt.title('Распределение параллаксов звезд GAIA DR2 with RV')
    plt.savefig('parallax.png')

    plt.clf()

    dv = np.abs(100 * df['radial_velocity_error'] / df['radial_velocity'])
    plt.hist(dv, bins=100, range=(0, 100), density=True)
    # plt.show()
    plt.xlabel('Относительная ошибка лучевой скорости, %')
    plt.ylabel('Доля звезд')
    plt.title('Распределение относительной ошибки лучевых скоростей звезд')
    plt.savefig('radial_velocity_error.png')

    plt.clf()

    dv = np.abs(100 * df['pmra_error'] / df['pmra'])
    plt.hist(dv, bins=100, range=(0, 100), density=True)
    # plt.show()
    plt.xlabel('Относительная ошибка собственных движений, %')
    plt.ylabel('Доля звезд')
    plt.title('Распределение относительной ошибки собственных движений звезд')
    plt.savefig('pmra_error.png')

    plt.clf()

    M = df['phot_g_mean_mag'] + 5 - 5 * np.log10(1000 / df['parallax'])
    plt.hist(M, bins=100, range=(-5, 10), density=True)
    # plt.show()
    plt.xlabel('Абсолютная звездная величина')
    plt.ylabel('Доля звезд')
    plt.title('Распределение абсолютной звездной величины')
    plt.savefig('M.png')

    plt.clf()

    px_error = df['parallax_error']
    plt.hist(px_error, bins=100, range=(0, 0.3))
    # plt.show()
    plt.xlabel('Точность параллакса, mas')
    plt.ylabel('Число звезд')
    plt.title('Распределение абсолютной точности параллакса')
    plt.savefig('abs_px_error.png')

    plt.clf()

    px_error = np.abs(100 * df['parallax_error'] / df['parallax'])
    plt.hist(px_error, bins=100, range=(-1, 100), density=True)
    # plt.show()
    plt.xlabel('Относительная точность параллакса, %')
    plt.ylabel('Доля звезд')
    plt.title('Распределение относительной точности параллакса')
    plt.savefig('sigma_px_error.png')

if __name__ == '__main__':
    main()