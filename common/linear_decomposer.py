from common.gaia.with_rv import slices, distance
import statsmodels.api as sm

from common.pandas2 import split


def prepare_linear_equations(df, compute_equations):
    equations = compute_equations(df)
    return split(equations.fillna(0))


def compute_coeffs(df, compute_equations):
    X, y = prepare_linear_equations(df, compute_equations)

    smm = sm.OLS(y, X)

    res = smm.fit()
    print(res.summary2())

    y_pred = res.predict(X)

    return smm, res, y_pred - y


def decompose(dataset, compute_equations, imax, step, slice, sort_lines):
    dists = []
    coefs = {}
    errors = {}

    def pizza(df):
        m_dist = distance(df.iloc[len(df) // 2])
        dists.append(m_dist)

        _, res, _ = compute_coeffs(df, compute_equations)

        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            key = res.params.keys()[i]
            if key not in coefs:
                coefs[key] = []
            if key not in errors:
                errors[key] = []
            coefs[key].append(val)
            errors[key].append(err)
            print(f'{key}={val}±{err}')

    slices(pizza, dataset, imax, step, slice)

    for figure in sort_lines:
        for label in figure:
            print(label, end='\t')
            for val, err in zip(coefs[label], errors[label]):
                print('{0:.1f}\t{1:.1f}'.format(val, err).replace('.', ','), end='\t')
            print()

    return coefs
