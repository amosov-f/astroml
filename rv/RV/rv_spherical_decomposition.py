from GC2019.gc2019 import read_gc2019
from common.spherical_decomposer import decompose_spherical, show_spherical_decomposition, show_spherical_decomposition_on_sphere
import pandas as pd
from common.gaia.with_rv import read_gaia_with_rv, slices
import numpy as np


MAX = 6_000_000
STEP = 75
SLICE = STEP


def main():
    df = read_gc2019()#read_gaia_with_rv()

    # print(f"Average mu error: {np.average(np.sqrt(df['mul_error'] ** 2 + df['mub_error'] ** 2))}")
    # # print(f"Average mub_error: {np.average()}")
    #
    #
    #
    # return
    df['y'] = df['radial_velocity'] / (1000 / df['px']) # км/с/кпк

    ncoeff = 16
    xls = [f"{i}\t" for i in range(ncoeff)]
    def pizza(ring):
        model = decompose_spherical(ring, ncoeff)
        coeffs, errors = show_spherical_decomposition(model, draw=False)
        for i in range(len(coeffs)):
            xls[i] += f"{coeffs[i]:.1f}\t{errors[i]:.1f}\t".replace('.', ',')

        # show_spherical_decomposition_on_sphere(model, 'Модель лучевых скоростей [км/с]')

    slices(pizza, df, MAX, STEP, SLICE, prepare=False)

    # for row in xls:
    #     print(row)



if __name__ == '__main__':
    main()
