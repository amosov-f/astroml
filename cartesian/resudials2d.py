from scipy import stats

from cartesian.common_cartesian import read_gaia_with_rv_residuals
from curve.sistematic3d import show_velocity
import numpy as np

def show_velocity2d():
    for r in [1, 5, 8]:
        for degree, title in [[1, 'линейной модели'], [2, 'квадратичной модели']]:
            df = read_gaia_with_rv_residuals(degree)

            print(df)

            min_x = -0.1
            max_x = 0.1
            min_y = -r
            max_y = r
            min_z = -2
            max_z = 2
            show_velocity(df, f"Остаточные скорости Vy по z", title + f" до {r} кпк", r,
                          min_x=min_x,
                          max_x=max_x,
                          min_y=min_y,
                          max_y=max_y,
                          min_z=min_z,
                          max_z=max_z,
                          bin_count=20)


def show_velocity_table():
    for d in [1, 2]:
        show_table(d, 5)

def show_table(degree, r):
    print(f'degree={degree}')

    df = read_gaia_with_rv_residuals(degree)

    min_x = -r
    max_x = r
    min_y = -r
    max_y = r
    min_z = -0.1
    max_z = 0.1

    df_x = df[(df.y >= min_y) & (df.y <= max_y) & (df.z >= min_z) & (df.z <= max_z)]
    df_y = df[(df.x >= min_x) & (df.x <= max_x) & (df.z >= min_z) & (df.z <= max_z)]

    step = 0.5
    bins = np.arange(min_x, max_x + step, step)

    xvx_means, _, _ = stats.binned_statistic(df_x.x, df_x.vx, statistic='mean', bins=bins)
    xvy_means, _, _ = stats.binned_statistic(df_x.x, df_x.vy, statistic='mean', bins=bins)
    xvz_means, _, _ = stats.binned_statistic(df_x.x, df_x.vz, statistic='mean', bins=bins)

    yvx_means, _, _ = stats.binned_statistic(df_y.y, df_y.vx, statistic='mean', bins=bins)
    yvy_means, _, _ = stats.binned_statistic(df_y.y, df_y.vy, statistic='mean', bins=bins)
    yvz_means, _, _ = stats.binned_statistic(df_y.y, df_y.vz, statistic='mean', bins=bins)

    print(end=',')
    for l in bins:
        r = l + step
        print(f'{l}', end=',')
    print()
    print('vx', end=',')
    for i in range(len(bins) - 1):
        print(f'{xvx_means[i]:.1f}', end=',')
    print()
    print('vy', end=',')
    for i in range(len(bins) - 1):
        print(f'{xvy_means[i]:.1f}', end=',')
    print()
    print('vz', end=',')
    for i in range(len(bins) - 1):
        print(f'{xvz_means[i]:.1f}', end=',')
    print()
    print()

    for l in bins:
        r = l + step
        print(f'{l}', end=',')
    print()
    for i in range(len(bins) - 1):
        print(f'{yvx_means[i]:.1f}', end=',')
    print()
    for i in range(len(bins) - 1):
        print(f'{yvy_means[i]:.1f}', end=',')
    print()
    for i in range(len(bins) - 1):
        print(f'{yvz_means[i]:.1f}', end=',')
    print()


if __name__ == '__main__':
    show_velocity2d()
    # show_velocity_table()
