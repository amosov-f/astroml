# astropy imports
import astropy.coordinates as coord
from astropy.table import QTable
import astropy.units as u
from astroquery.gaia import Gaia

# Third-party imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# gala imports
# import gala.coordinates as gc
# import gala.dynamics as gd
# import gala.potential as gp
# from gala.units import galactic
import pandas as pd
from matplotlib import cm
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm


def main():
    # query_text = '''SELECT TOP 1000000 ra, dec, parallax, pmra, pmdec, radial_velocity,
    # phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
    # FROM gaiadr3.gaia_source
    # WHERE parallax_over_error > 4 AND parallax > 1.0 / 7.906 AND radial_velocity IS NOT null
    # '''
    # job = Gaia.launch_job(query_text)
    # gaia_data = job.get_results()
    # gaia_data.write('gaia_data.fits')
    #
    # gaia_data = QTable.read('gaia_data.fits')

    df = pd.read_csv('gaia_dr3_with_rv2.csv', sep=',')

    # df = df[(df.parallax > 1.0 / 7.906) & (df.parallax_over_error > 4)]
    df = df[(df.parallax_over_error > 4)]

    # df = df[(df.parallax > 10) & (df.parallax_over_error > 10)]

    df = df.sort_values(by='parallax', ascending=False)
    df = df.head(30000000)

    print(df)

    dist = coord.Distance(parallax=df.parallax.values * u.mas)
    print(dist.min(), dist.max())

    # c = coord.SkyCoord(ra=df['ra'],
    #                    dec=df['dec'],
    #                    distance=dist,
    #                    pm_ra_cosdec=df['pmra'],
    #                    pm_dec=df['pmdec'],
    #                    radial_velocity=df['radial_velocity'])

    c = coord.SkyCoord(ra=df.ra.values * u.deg,
                       dec=df.dec.values * u.deg,
                       pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
                       pm_dec=df.pmdec.values * u.mas / u.yr,
                       radial_velocity=df.radial_velocity.values * u.km / u.s,
                       distance=dist)
    g = c.galactic

    # print(c[:4])
    # print('')
    print(c.galactic[:4])

    print(dist.distmod)

    # g = df.phot_g_mean_mag.values

    M_G = df.phot_g_mean_mag.values * u.mag - dist.distmod
    BP_RP = df.phot_bp_mean_mag.values * u.mag - df.phot_rp_mean_mag.values * u.mag

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    #
    # ax.plot(BP_RP.value, M_G.value,
    #         marker='.', linestyle='none', alpha=0.3)
    #
    # ax.set_xlim(-1, 5)
    # ax.set_ylim(16, -5)
    #
    ax.set_xlabel('$G_{BP}-G_{RP}$')
    ax.set_ylabel('$M_G$')

    h = ax.hist2d(BP_RP.value, M_G.value, range=[[0, 3], [-1, 11]], bins=1000, cmap=cm.binary)

    fig.colorbar(h[3], ax=ax)


    x1 = 0.877
    y1 = 2.35
    x2 = 2
    y2 = 7.25

    k = (y2 - y1) / (x2 - x1)
    b = y1 - k * x1

    # Ax1 + By1 + C = 0
    # Ax2 + By2 + C = 0
    #
    # A = (-C -By1) / x1
    #
    # y = (-C -Ax) / B = - A/B x - X/B


    df = pd.DataFrame(dict(x=g.cartesian.x.value / 1000,
                           y=g.cartesian.y.value / 1000,
                           z=g.cartesian.z.value / 1000,
                           vx=g.velocity.d_x.value,
                           vy=g.velocity.d_y.value,
                           vz=g.velocity.d_z.value,
                           lum=k * BP_RP.value + b - M_G.value))

    print(df)

    pizza(df)

    pizza(df[(df.lum < 0)])

    pizza(df[(df.lum > 0)])


    ax.axline((x1, y1), (x2, y2), color='black', label=f'Mg = {k:.2}*(Gbp-Grp) - {abs(b):.2}m')

    ax.text(0.25, 9, 'Звезды\nглавной\nпоследовательности', fontsize=12)
    ax.text(2.4, 1, 'Гиганты', fontsize=12)

    # plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()

    # plt.set_cmap(cm.binary)
    plt.legend()

    plt.show()

    # print(gaia_data)
    # gaia_data.write('gaia_data.fits')


def pizza(df):

    print(len(df))

    X1 = df[['x', 'y', 'z']]
    p = PolynomialFeatures(degree=2).fit(X1)
    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names_out(X1.columns))

    models = []
    for vname, y in [
        ('vx', df.vx),
        ('vy', df.vy),
        ('vz', df.vz),
    ]:

        # print(X)
        # print(y)

        ols = sm.OLS(list(y), X)
        res = ols.fit()

        print(vname, end='\t')

        models.append(res)

        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            col = res.params.keys()[i]
            print(f'\t{val:.3f}±{err:.3f} * {col}', end=',')

        print()

    vx = models[0]
    vy = models[1]
    vz = models[2]

    U = -vx.params[0]
    dU = vx.bse[0]
    V = -vy.params[0]
    dV = vy.bse[0]
    W = -vz.params[0]
    dW = vz.bse[0]

    M11 = vx.params[1]
    dM11 = vx.bse[1]
    M22 = vy.params[2]
    dM22 = vy.bse[2]
    M33 = vz.params[3]
    dM33 = vz.bse[3]

    w1 = (- vy.params[3] + vz.params[2]) / 2
    dw1 = sigma(vy.bse[3], vz.bse[2]) / 2
    w2 = (vx.params[3] - vz.params[1]) / 2
    dw2 = sigma(vx.bse[3], vz.bse[1]) / 2
    w3 = (- vx.params[2] + vy.params[1]) / 2
    dw3 = sigma(vx.bse[2], vy.bse[1]) / 2

    M12 = (vx.params[2] + vy.params[1]) / 2
    dM12 = sigma(vx.bse[2], vy.bse[1]) / 2
    M13 = (vx.params[3] + vz.params[1]) / 2
    dM13 = sigma(vx.bse[3], vz.bse[1]) / 2
    M23 = (vy.params[3] + vz.params[2]) / 2
    dM23 = sigma(vy.bse[3], vz.bse[2]) / 2

    pp(U, dU)
    pp(V, dV)
    pp(W, dW)
    pp(w1, dw1)
    pp(w2, dw2)
    pp(w3, dw3)
    pp(M11, dM11)
    pp(M22, dM22)
    pp(M33, dM33)
    pp(M12, dM12)
    pp(M13, dM13)
    pp(M23, dM23)

    for v in [vx, vy, vz]:
        pp(v.params[4], v.bse[4])  # axx
        pp(v.params[7], v.bse[7])  # ayy
        pp(v.params[9], v.bse[9])  # azz
        pp(v.params[5], v.bse[5])  # axy
        pp(v.params[6], v.bse[6])  # axz
        pp(v.params[8], v.bse[8])  # ayz



def pp(v, err):
    print(f'{v:.3f}±{err:.3f}')

def sigma(err1, err2):
    return (err1 ** 2 + err2 ** 2) ** 0.5

if __name__ == '__main__':
    main()
