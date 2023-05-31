from matplotlib.colors import LogNorm

from cartesian.common_cartesian import read_gaia_with_rv_cartesian
from matplotlib import pyplot as plt


def main():
    df = read_gaia_with_rv_cartesian()


    h = plt.hist2d(df.x, df.y, bins=100, range=[[-8, 8], [-8, 8]], norm=LogNorm())
    plt.colorbar(h[3])
    plt.title('Гистограмма распределения звезд GAIA DR3 with RV\nв плоскости XY галактической системы координат')
    plt.show()


if __name__ == '__main__':
    main()
