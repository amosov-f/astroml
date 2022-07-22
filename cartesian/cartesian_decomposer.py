from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

from cartesian.common_cartesian import *
import pandas as pd


def main():
    g = read_pizza()
    df = to_cartesian(g)
    X1 = df[['x', 'y', 'z']]
    p = PolynomialFeatures(degree=2).fit(X1)
    # print(p.get_feature_names(X.columns))

    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names(X1.columns))

    print(X.columns)

    precision = 0.01

    for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:
        model = linear_model.Lasso(alpha=1, fit_intercept=False)
    # model = make_pipeline(PolynomialFeatures(degree=1), LinearRegression(fit_intercept=False))
    #     model = LinearRegression(fit_intercept=False)
        model.fit(X, y)

        print(vname + ' = ', end='\t')
        for i, col in enumerate(X.columns.values):
            if abs(model.coef_[i]) > 0:
                print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

    # print(model.score(X, y))
        print()

    print()

    r2 = 2.95
    dr = 0.100

    _x, _y, _z = np.mgrid[-r2:r2 + dr:dr, -r2:r2 + dr:dr, -r2:r2 + dr:dr]
    P = np.array([_x.flatten(), _y.flatten(), _z.flatten()])
    # print(P)
    P = pd.DataFrame(dict(x=P[0], y=P[1], z=P[2]))
    p = PolynomialFeatures(degree=2).fit(P)

    PX = pd.DataFrame(p.transform(P), columns=p.get_feature_names(P.columns))

    for name, f in [('U', U),
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
        for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:

            # print(features)
            # model = linear_model.Lasso(alpha=1, fit_intercept=False)
            # model = make_pipeline(PolynomialFeatures(degree=2), LinearRegression(fit_intercept=False))
            model = LinearRegression(fit_intercept=False)
            # model = linear_model.Lasso(alpha=0.001, fit_intercept=False)
            model.fit(PX, y)

            # print(model.)
            print(vname + f' = ', end='\t')
            for i, col in enumerate(PX.columns.values):
                if abs(model.coef_[i]) > precision:
                    print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

            # print(model.score(X, y))
            print()

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