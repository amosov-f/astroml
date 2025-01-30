import pandas as pd
from astropy.coordinates import CartesianRepresentation
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm

from cartesian.common_cartesian import *
from common.spherical import taj
from curve.sistematic3d import show_velocity
from numpy import sin, cos


def compute_full_quadric_decomposition_lasso():
    df = read_gaia_with_rv_cartesian()

    X1 = df[['x', 'y', 'z']]

    f = PolynomialFeatures(degree=2)
    p = f.fit(X1)

    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names_out(X1.columns))

    alpha = 5

    f = linear_model.Lasso
    print(f)
    print(alpha)

    print(X.columns)

    for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:
        model = f(alpha=alpha, fit_intercept=False)

        model.fit(X, y)

        print(vname + ' = ', end='\t')
        for i, col in enumerate(X.columns.values):
            if abs(model.coef_[i]) > 0:
                print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

        print()


def compute_full_quadric_decomposition_galactic_coordinates(simple=False, no_radial_velocity=False):
    g = read_galactic()

    df = pd.DataFrame(dict(
        l=np.deg2rad(np.array(g.l)),
        b=np.deg2rad(np.array(g.b)),
        r=np.array(g.distance) / 1000,
        pm_b=g.pm_b,
        pm_l_cosb=g.pm_l_cosb,
        radial_velocity=g.radial_velocity,
    ))

    df = df[(abs(df.r*cos(df.b)*cos(df.l)) < 0.2)]

    print('galcor')
    print(len(df))

    l = df.l
    b = df.b
    r = df.r

    print(min(r))
    print(max(r))
    cosl = cos(l)
    sinl = sin(l)
    cosb = cos(b)
    sinb = sin(b)
    cosl2 = cosl ** 2
    sinl2 = sinl ** 2
    cosb2 = cosb ** 2
    sinb2 = sinb ** 2
    cosl3 = cosl ** 3
    sinl3 = sinl ** 3
    cosb3 = cosb ** 3
    sinb3 = sinb ** 3
    
    sin2l = sin(2*l)
    cos2l = cos(2*l)
    sin2b = sin(2*b)
    cos2b = cos(2*b)

    kmub = pd.DataFrame(dict(U=cosl * sinb / r,
                             V=sinl * sinb / r,
                             W=-cosb / r,
                             w1=sinl,
                             w2=-cosl,
                             w3=0,
                             M12=-0.5 * sin2b * sin2l,
                             M13=cos2b * cosl,
                             M23=cos2b * sinl,

                             # M11=-0.5 * sin2b * cosl2,
                             # M22=-0.5 * sin2b * sinl2,
                             # M33=0.5 * sin2b,

                             C=-0.5 * sin2b * cos2l,
                             K=-0.5 * sin2b,

                             # axx=-r * cosl3 * cosb2 * sinb,
                             # ayy=-r * cosl * sinl2 * cosb2 * sinb,
                             # bxy=-r * cosl * sinl2 * cosb2 * sinb,
                             # azz=-r * cosl * sinb3,
                             # axy=-r * cosl2 * sinl * cosb2 * sinb,
                             # bxx=-r * cosl2 * sinl * cosb2 * sinb,
                             # axz=-r * cosl2 * cosb * sinb2,
                             # ayz=-r * cosl * sinl * cosb * sinb2,
                             # bxz=-r * cosl * sinl * cosb * sinb2,
                             # byy=-r * sinl3 * cosb2 * sinb,
                             # bzz=-r * sinl * sinb3,
                             # cxx=r * cosl2 * cosb3,
                             # cyy=r * sinl2 * cosb3,
                             # czz=r * cosb * sinb2,
                             # cxy=r * cosl * sinl * cosb3,
                             # cxz=r * cosl * cosb2 * sinb,
                             # cyz=r * sinl * cosb2 * sinb,
                             y=K * df.pm_b,
                             ))
    kmulcosb = pd.DataFrame(dict(U=sinl / r,
                                 V=-cosl / r,
                                 # w1=-sinb * cosl,
                                 # w2=-sinb * sinl,
                                 w3=cosb,
                                 M12=cosb * cos2l,
                                 # M13=-sinb * sinl,
                                 # M23=sinb * cosl,

                                 # M11=-0.5 * cosb * sin2l,
                                 # M22=0.5 * cosb * sin2l,
                                 # M33=0,

                                 C=-cosb * sin2l,
                                 K=0.0,

                                 # axx=-r * cosl2 * sinl * cosb2,
                                 # bxy=+r * cosl2 * sinl * cosb2,
                                 # ayy=-r * sinl3 * cosb2,
                                 # azz=-r * sinl * sinb2,
                                 # axy=-r * cosl * sinl2 * cosb2,
                                 # byy=+r * cosl * sinl2 * cosb2,
                                 # axz=-r * cosl * sinl * cosb * sinb,
                                 # byz=+r * cosl * sinl * cosb * sinb,
                                 # ayz=-r * sinl2 * cosb * sinb,
                                 # bxx=r * cosl3 * cosb2,
                                 bzz=r * cosl * sinb2,
                                 # bxz=r * cosl2 * cosb * sinb,
                                 y=K * df.pm_l_cosb,
                                 ))
    if not no_radial_velocity:
        vrr = pd.DataFrame(dict(U=-cosl * cosb / r,
                                V=-sinl * cosb / r,
                                W=-sinb / r,
                                M12=cosb2 * sin2l,
                                M13=sin2b * cosl,
                                M23=sin2b * sinl,
                                M11=cosb2 * cosl2,
                                M22=cosb2 * sinl2,
                                M33=sinb2,
                                axx=r * cosl3 * cosb3,
                                ayy=r * cosl * sinl2 * cosb3,
                                bxy=r * cosl * sinl2 * cosb3,
                                azz=r * cosl * cosb * sinb2,
                                cxz=r * cosl * cosb * sinb2,
                                axy=r * cosl2 * sinl * cosb3,
                                bxx=r * cosl2 * sinl * cosb3,
                                axz=r * cosl2 * cosb2 * sinb,
                                cxx=r * cosl2 * cosb2 * sinb,
                                ayz=r * cosl * sinl * cosb2 * sinb,
                                bxz=r * cosl * sinl * cosb2 * sinb,
                                cxy=r * cosl * sinl * cosb2 * sinb,
                                byy=r * sinl3 * cosb3,
                                bzz=r * sinl * cosb * sinb2,
                                cyz=r * sinl * cosb * sinb2,
                                byz=r * sinl2 * cosb2 * sinb,
                                cyy=r * sinl2 * cosb2 * sinb,
                                czz=r * sinb3,
                                y=df.radial_velocity / r,
                                ))
    print('params computed')
    X = pd.concat([
        kmub,
        kmulcosb,
        # vrr,
    ]).fillna(0)
    y = X[['y']]
    X = X.drop(['y'], axis=1)
    print('fitting model')
    res = sm.OLS(y, X).fit()
    print(res.summary())


def main():
    df = read_gaia_with_rv_cartesian()

    # print(df)

    X1 = df[['x', 'y', 'z']]

    # X1['-1'] = -1
    # X1['z^2'] = X1.z**2
    # X = X1

    #
    # linear_vx = model.fit(X1, df.vx).predict(X1)
    # linear_vy = model.fit(X1, df.vy).predict(X1)
    # linear_vz = model.fit(X1, df.vz).predict(X1)
    #
    # print(df.x)
    #
    # c = Galactic(u=df.x * u.kpc, v=df.y * u.kpc, w=df.z * u.kpc, representation_type="cartesian", U=linear_vx * u.km/u.s, V=linear_vy * u.km/u.s,  W=linear_vz * u.km/u.s)

    p = PolynomialFeatures(degree=2).fit(X1)
    # print(p.get_feature_names(X.columns))

    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names_out(X1.columns))

    precision = 0.01

    alpha = 0.5
    print(alpha)

    print(X.columns)

    models = []

    for vname, y in [
        ('vx', df.vx),
        ('vy', df.vy),
        ('vz', df.vz),
    ]:
        # model = linear_model.Lasso(alpha=alpha, fit_intercept=False)
        # model = make_pipeline(PolynomialFeatures(degree=1), LinearRegression(fit_intercept=False))
        # model = LinearRegression(fit_intercept=False)
        # model.fit(X, y)

        # y_pred = model.predict(X)
        # diff = y - y_pred
        #

        ols = sm.OLS(y, X)
        res = ols.fit()

        # print(res.summary())

        models.append(res)

        # model.coef_[]

        print(vname + ' = ', end='\t')
        # for i, col in enumerate(X.columns.values):
        #     if abs(model.coef_[i]) > precision:
        #         print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')
        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            col = res.params.keys()[i]
            print(f'+ {val:.2f} * {col}', end='\t')

        print()

        print(vname, end='\t')
        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            col = res.params.keys()[i]
            print(f'\t{val:.3f}±{err:.3f} * {col}', end=',')

        print()

    return

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

    print()

    print(f'w1={w1:.3f}±{dw1:.3f}')
    print(f'w2={w2:.3f}±{dw2:.3f}')
    print(f'w3={w3:.3f}±{dw3:.3f}')

    M12 = (vx.params[2] + vy.params[1]) / 2
    dM12 = sigma(vx.bse[2], vy.bse[1]) / 2
    M13 = (vx.params[3] + vz.params[1]) / 2
    dM13 = sigma(vx.bse[3], vz.bse[1]) / 2
    M23 = (vy.params[3] + vz.params[2]) / 2
    dM23 = sigma(vy.bse[3], vz.bse[2]) / 2

    print()

    print(f'M12={M12:.3f}±{dM12:.3f}')
    print(f'M13={M13:.3f}±{dM13:.3f}')
    print(f'M23={M23:.3f}±{dM23:.3f}')

    # print(model.score(X, y))

    print()

    # compute_decomposition()


def sigma(err1, err2):
    return (err1 ** 2 + err2 ** 2) ** 0.5


def compute_decomposition():
    # for i, z in enumerate([-1.5, -1, -0.5, 0, 0.5, 1, 1.5]):
    #     pasta(z, i)
    pasta(None, 0)


def pasta(z, order):
    precision = 0.01

    r2 = 1.5
    dr = 1
    count = 10000

    _x, _y, _z = np.mgrid[-r2:r2 + dr:dr, -r2:r2 + dr:dr, -r2:r2 + dr:dr]
    P = np.array([_x.flatten(), _y.flatten(), _z.flatten()])
    P = [
        np.random.uniform(-r2, r2, count),
        np.random.uniform(-r2, r2, count),
        np.full(count, z) if z is not None else np.random.uniform(-r2, r2, count),
        # np.random.uniform(-r2, r2, count)
    ]

    # print(P)
    P = pd.DataFrame(dict(x=P[0], y=P[1], z=P[2]))
    p = PolynomialFeatures(degree=2).fit(P)

    PX = pd.DataFrame(p.transform(P), columns=p.get_feature_names(P.columns))

    for name, f in [
        ('U', U),
        ('V', V),
        ('W', W),

        ('w1', w1),
        ('w2', w2),
        ('w3', w3),

        ('M12', M12),
        ('M13', M13),
        ('M23', M23),
        ('M11', M11),
        ('M22', M22),
        ('M33', M33),

        ('dw1dr1', dw1dr1),
        ('dw1dr2', dw1dr2),
        ('dw1dr3', dw1dr3),

        ('dw2dr1', dw2dr1),
        ('dw2dr2', dw2dr2),
        ('dw2dr3', dw2dr3),

        ('dw3dr1', dw3dr1),
        ('dw3dr2', dw3dr2),
        ('dw3dr3', dw3dr3),

        ('dM11dr1', dM11dr1),
        ('dM11dr2', dM11dr2),
        ('dM11dr3', dM11dr3),

        ('dM12dr1', dM12dr1),
        ('dM12dr2', dM12dr2),
        ('dM12dr3', dM12dr3),

        ('dM13dr1', dM13dr1),
        ('dM13dr2', dM13dr2),
        ('dM13dr3', dM13dr3),

        ('dM22dr1', dM22dr1),
        ('dM22dr2', dM22dr2),
        ('dM22dr3', dM22dr3),

        ('dM23dr1', dM23dr1),
        ('dM23dr2', dM23dr2),
        ('dM23dr3', dM23dr3),

        ('dM33dr1', dM33dr1),
        ('dM33dr2', dM33dr2),
        ('dM33dr3', dM33dr3),
    ]:
        print(name)
        df = apply(P, f)

        pd.set_option('display.max_rows', None)

        # print(df)
        for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:

            # print(features)
            # model = linear_model.Lasso(alpha=1, fit_intercept=False)
            # model = make_pipeline(PolynomialFeatures(degree=2), LinearRegression(fit_intercept=False))
            # model = LinearRegression(fit_intercept=False)
            model = linear_model.Lasso(alpha=0.001, fit_intercept=False)
            model.fit(PX, y)

            # print(model.)
            print(vname + f' = ', end='\t')
            for i, col in enumerate(PX.columns.values):
                if abs(model.coef_[i]) > precision:
                    print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

            # print(model.score(X, y))
            print()

        min_x = np.min(df.x)
        max_x = np.max(df.x)

        min_y = np.min(df.y)
        max_y = np.max(df.y)

        if z:
            show_velocity(df, name, f"{z}", order, min_x, max_x, min_y, max_y, 10)

        # ax.quiver(X, Y, x_means, y_means, C, width=0.003)

        print()
        print()


if __name__ == '__main__':
    # X, Y, Z = np.mgrid[-1000:1000+1:10, -1000:1000+1:10, -1000:1000+1:10]
    # P = np.array([X.flatten(), Y.flatten(), Z.flatten()])
    # c = SkyCoord(x=P[0], y=P[1], z=P[2], unit='pc', representation_type='cartesian')
    # print(np.deg2rad(np.array(c.galactic.l)))
    # print(np.deg2rad(np.array(c.galactic.b)))
    # print(np.array(c.galactic.distance) / 1000)
    # print(c.galactic)
    # print(P)
    # print('X')
    # print(X)
    # print('Y')
    # print(Y)
    # print('Z')
    # print(Z)
    main()
    # compute_full_quadric_decomposition_lasso()
    # compute_decomposition()
    # compute_full_quadric_decomposition_galactic_coordinates(simple=True, no_radial_velocity=True)
