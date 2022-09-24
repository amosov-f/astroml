import pandas as pd
from astropy.coordinates import CartesianRepresentation
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

from cartesian.common_cartesian import *
from curve.sistematic3d import show_velocity


def main():
    # df = read_gaia_with_rv_cartesian()

    # X1 = df[['x', 'y', 'z']]

    # model = linear_model.LinearRegression()
    # # linear_vx = (df.vx - model.fit(X1, df.vx).predict(X1)).values
    # # linear_vy = (df.vy - model.fit(X1, df.vy).predict(X1)).values
    # # linear_vz = (df.vz - model.fit(X1, df.vz).predict(X1)).values
    #
    # linear_vx = model.fit(X1, df.vx).predict(X1)
    # linear_vy = model.fit(X1, df.vy).predict(X1)
    # linear_vz = model.fit(X1, df.vz).predict(X1)
    #
    # print(df.x)
    #
    # c = Galactic(u=df.x * u.kpc, v=df.y * u.kpc, w=df.z * u.kpc, representation_type="cartesian", U=linear_vx * u.km/u.s, V=linear_vy * u.km/u.s,  W=linear_vz * u.km/u.s)
    c = read_radec()
    g = c.transform_to(Galactocentric)

    df = to_cartesian(g)

    X1 = df[['x', 'y', 'z']]

    p = PolynomialFeatures(degree=2).fit(X1)
    # print(p.get_feature_names(X.columns))

    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names(X1.columns))

    precision = 0.01

    alpha = 0.5
    print(alpha)

    print(X.columns)

    for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:
        model = linear_model.Lasso(alpha=alpha, fit_intercept=False)
    # model = make_pipeline(PolynomialFeatures(degree=1), LinearRegression(fit_intercept=False))
    #     model = LinearRegression(fit_intercept=False)
        model.fit(X, y)

        y_pred = model.predict(X)
        diff = y - y_pred



        # model.coef_[]

        print(vname + ' = ', end='\t')
        for i, col in enumerate(X.columns.values):
            if abs(model.coef_[i]) > precision:
                print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

    # print(model.score(X, y))
        print()

    print()

    # compute_decomposition()



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
    # compute_decomposition()