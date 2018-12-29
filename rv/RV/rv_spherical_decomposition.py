from common.spherical_decomposer import decompose_spherical, show_spherical_decomposition, show_spherical_decomposition_on_sphere
import pandas as pd
from common.gaia.with_rv import read_gaia_with_rv


def main():
    df = read_gaia_with_rv()
    df.y = df.radial_velocity
    model = decompose_spherical(df, 50)
    show_spherical_decomposition(model)
    show_spherical_decomposition_on_sphere(model, 'Radial velocity model')


if __name__ == '__main__':
    main()
