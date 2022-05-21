from common.gaia.with_rv import slices
from common.spherical_decomposer import show_spherical_decomposition_on_sphere, decompose_spherical, \
    show_spherical_decomposition
from curve.sistematic3d import read_gaia_with_rv_redisuals_xyz, read_gaia_with_rv_residuals


def velocity_model(x, y, z):
    return - x * y, y ** 2, x ** 2 + y ** 2


def main():
    df = read_gaia_with_rv_residuals()#read_gaia_with_rv()
    df_xyz = read_gaia_with_rv_redisuals_xyz()

    # print(f"Average mu error: {np.average(np.sqrt(df['mul_error'] ** 2 + df['mub_error'] ** 2))}")
    # # print(f"Average mub_error: {np.average()}")
    #
    #
    #
    # return
    # df[]


    df['y'] = df.vr / (1 / df.px) # км/с/кпк

    ncoeff = 16
    xls = [f"{i}\t" for i in range(ncoeff)]
    def pizza(ring, _, _1):
        model = decompose_spherical(ring, ncoeff)
        coeffs, errors = show_spherical_decomposition(model, draw=False)
        for i in range(len(coeffs)):
            xls[i] += f"{coeffs[i]:.1f}\t{errors[i]:.1f}\t".replace('.', ',')

        # show_spherical_decomposition_on_sphere(model, 'Модель лучевых скоростей [км/с]')

    slices(pizza, df)

    # for row in xls:
    #     print(row)



if __name__ == '__main__':
    main()