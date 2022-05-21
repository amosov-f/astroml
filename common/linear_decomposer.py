from common.gaia.with_rv import slices, distance
import statsmodels.api as sm
import numpy as np

from common.pandas2 import split


def prepare_linear_equations(df, compute_equations):
    equations = compute_equations(df).fillna(0)
    zero_columns = equations.columns[(equations == 0).all()]
    print("Totally zero columns, removing: " + str(zero_columns))
    equations = equations.loc[:, (equations != 0).any()]
    return split(equations)


def compute_coeffs(df, compute_equations):
    X, y = prepare_linear_equations(df, compute_equations)

    corr = np.corrcoef(X, rowvar=0)


    print('evgen:')
    names = {}
    for i, name in enumerate(list(X)):
        names[i] = name
        print(f'{i}: {name}')

    w, v = np.linalg.eig(corr)

    # print('w')
    # print(w)
    #
    # print('v')
    # print(v)


    for i, name in enumerate(list(X)):
        if w[i] < 0.0001:
            print(f'{name}: {w[i]:.6f}')
            for j, col in enumerate(v[:,i]):
                if abs(col) > 0.1:
                    print(f'\t{names[j]}={col:.6f}')



    smm = sm.OLS(y, X)

    res = smm.fit()
    print(res.summary2())

    y_pred = res.predict(X)

    return smm, res, y_pred - y


def decompose(dataset, compute_equations, imax, step, slice, sort_lines, filter=None):
    dists = []
    coefs = {}
    # errors = {}

    def pizza(df):
        m_dist = distance(df.iloc[len(df) // 2]) / 1000
        dists.append(m_dist)

        _, res, _ = compute_coeffs(df, compute_equations)

        cur_coefs = {}
        cur_errors = {}
        for i in range(len(res.params)):
            val = res.params[i]
            err = res.bse[i]
            key = res.params.keys()[i]
            cur_coefs[key] = val
            cur_errors[key] = err

        # V = cur_coefs.get('V', 0)
        # w1 = cur_coefs.get('w1', 0)
        # w2 = cur_coefs.get('w2', 0)
        # w3 = cur_coefs.get('w3', 0)
        # dw1dr3 = cur_coefs.get('dw1dr3', 0)
        # dw3dr1 = cur_coefs.get('dw3dr1', 0)
        # dM11dr2 =  cur_coefs.get('dM11dr2', 0)
        # dM12dr1 = cur_coefs.get('dM12dr1', 0)
        # dM12dr2 = cur_coefs.get('dM12dr2', 0)
        # dM22dr3 = cur_coefs.get('dM22dr3', 0)
        # dM23dr2 = cur_coefs.get('dM23dr2', 0)
        # dM33dr2 = cur_coefs.get('dM33dr2', 0)
        # dM22dr2 = cur_coefs.get('dM22dr2', 0)
        # dM23dr3 = cur_coefs.get('dM23dr3', 0)
        #
        # cur_coefs['t101'] = 2.894 * w3
        # cur_coefs['s110'] = -2.894 * V / m_dist + (- 1.447 * dw1dr3 + 1.447 * dw3dr1 - 0.289 * dM11dr2 + 0.724 * dM12dr1 - 0.145 * dM12dr2 + 0.868 * dM22dr3 - 0.289 * dM23dr2 - 0.289 * dM33dr2) * m_dist
        #
        # cur_coefs['t211'] = (-1.121 * dw1dr3 + 1.121 * dw3dr1 - 0.374 * dM11dr2 + 0.560 * dM12dr1 + 0.187 * dM12dr2 - 0.374 * dM23dr3 + 0.374 * dM33dr2) * m_dist
        # cur_coefs['s310'] = (-0.126 * dM11dr2 + 0.253 * dM12dr2 - 0.379 * dM22dr2 + 1.012 * dM23dr3 + 0.505 * dM33dr2) * m_dist
        # cur_coefs['v310'] = (-0.109 * dM11dr2 - 0.328 * dM22dr2 + 0.438 * dM33dr2 - 0.219 * dM12dr1 + 0.875 * dM23dr3) * m_dist
        #
        # cur_coefs['dist'] = m_dist
        #
        # print(f"t211 = t_6 = {cur_coefs['t211']} = (-1.121 * {dw1dr3} (dw1dr3) + 1.121 * {dw3dr1} (dw3dr1) - 0.374 * {dM11dr2} (dM11dr2) + 0.560 * {dM12dr1} (dM12dr1) + 0.187 * {dM12dr2} (dM12dr2) - 0.374 * {dM23dr3} (dM23dr3) + 0.374 * {dM33dr2} (dM33dr2)) * {m_dist}")
        # print(f"s310 = s_10 = {cur_coefs['s310']} = (-0.126 * {dM11dr2} (dM11dr2) + 0.253 * {dM12dr2} (dM12dr2) - 0.379 * {dM22dr2} (dM22dr2) + 1.012 * {dM23dr3} (dM23dr3) + 0.505 * {dM33dr2} (dM33dr2)) * {m_dist}")
        # print(f"v310 = v_10 = {cur_coefs['v310']} = (-0.109 * {dM11dr2} (dM11dr2) - 0.328 * {dM22dr2} (dM22dr2) + 0.438 * {dM33dr2} (dM33dr2) - 0.219 * {dM12dr1} (dM12dr1) + 0.875 * {dM23dr3} (dM23dr3)) * {m_dist}")

        # print(f"t211={cur_coefs['t211']}, s310={cur_coefs['s310']}, v310={cur_coefs['v310']}")

        for key, val in cur_coefs.items():
            if key not in coefs:
                coefs[key] = []
            coefs[key].append(val)

        # for key, err in cur_errors.items():
        #     if key not in errors:
        #         errors[key] = []
        #     errors[key].append(err)

            # print(f'{key}={val}±{err}')

    average_dists = slices(pizza, dataset, imax, step, slice, filter=filter)

    print('\t' + '\t\t'.join(average_dists))
    for figure in sort_lines:
        for label in figure:
            if label in coefs:
                print(label, end='\t')
                for val in coefs[label]:
                    print('{0:.1f}'.format(val).replace('.', ','), end='\t')
                print()

    return coefs, average_dists
