from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

from cartesian.common_cartesian import read_gaia_with_rv_cartesian, read_pizza, to_cartesian, w1, w2, w3, M12, M13, M23, \
    M11, M22, M33, U, V, W
import pandas as pd


def main():
    g = read_pizza()
    df = to_cartesian(g)
    X1 = df[['x', 'y', 'z']]
    p = PolynomialFeatures(degree=2).fit(X1)
    # print(p.get_feature_names(X.columns))

    X = pd.DataFrame(p.transform(X1), columns=p.get_feature_names(X1.columns))

    print(X.columns)


    for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:
        # y = df.vx



    # print(features)
        model = linear_model.Lasso(alpha=1, fit_intercept=False)
    # model = make_pipeline(PolynomialFeatures(degree=1), LinearRegression(fit_intercept=False))
    #     model = LinearRegression(fit_intercept=False)
        model.fit(X, y)

        print(vname + ' = ', end='\t')
    # print(model.)
        for i, col in enumerate(X.columns.values):
            print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

    # print(model.score(X, y))
        print()

    print()

    for name, df in [('U', U(g)), ('V', V(g)), ('W', W(g)), ('w1', w1(g)), ('w2', w2(g)), ('w3', w3(g)), ('M12', M12(g)), ('M13', M13(g)), ('M23', M23(g)), ('M11', M11(g)), ('M22', M22(g)), ('M33', M33(g))]:
        print(name)
        for vname, y in [('vx', df.vx), ('vy', df.vy), ('vz', df.vz)]:
            # y = df.vx

            # print(features)
            # model = linear_model.Lasso(alpha=1, fit_intercept=False)
            # model = make_pipeline(PolynomialFeatures(degree=1), LinearRegression(fit_intercept=False))
            model = LinearRegression()
            model.fit(X1, y)

            # print(model.)
            print(vname + f' = {model.intercept_}', end='\t')
            for i, col in enumerate(X1.columns.values):
                print(f'+ {model.coef_[i]:.2f} * {col}', end='\t')

            # print(model.score(X, y))
            print()

        print()
        print()



if __name__ == '__main__':
    main()