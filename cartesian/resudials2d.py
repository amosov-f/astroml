from cartesian.common_cartesian import read_gaia_with_rv_residuals
from curve.sistematic3d import show_velocity

if __name__ == '__main__':
    for degree, title in [[1, 'линейной модели'], [2, 'квадратичной модели']]:
        df = read_gaia_with_rv_residuals(degree)

        print(df)

        min_x = -8
        max_x = 8
        min_y = -8
        max_y = 8
        min_z = -100500
        max_z = 100500
        show_velocity(df, f"Остаточные скорости", title, 0,
                      min_x=min_x,
                      max_x=max_x,
                      min_y=min_y,
                      max_y=max_y,
                      min_z=min_z,
                      max_z=max_z,
                      bin_count=20)