from astropy.coordinates import Galactocentric

from cartesian.common_cartesian import read_galactic, read_radec
from common.gaia.with_rv import read_gaia_with_rv_xyz
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt


def main():
    c = read_radec()
    gc = c.transform_to(Galactocentric)
    df = pd.DataFrame(dict(x=gc.cartesian.x.value,
                             y=gc.cartesian.y.value,
                             z=gc.cartesian.z.value,
                             vx=gc.velocity.d_x.value,
                             vy=gc.velocity.d_y.value,
                             vz=gc.velocity.d_z.value))
    df = df[((df.y < 10) & (df.z < 10))]
    bin_means, bin_edges, binnumber = stats.binned_statistic(df.x, df.vy, bins=200, range=[-16000, 0], statistic='median')
    plt.stairs(bin_means, bin_edges)
    plt.xlabel('Галактическая координата X, пк')
    plt.ylabel('Галактическая скорость Vy, км/с/кпк')
    plt.title('Кривая вращения Галактики на луче зрения от Солнца к Центру')
    plt.show()


if __name__ == '__main__':
    main()