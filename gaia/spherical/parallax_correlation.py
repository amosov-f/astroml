import pandas as pd
import statsmodels.api as sm
from matplotlib import pyplot as plt

if __name__ == '__main__':
    tgas_xhip = pd.read_table('tgas_xhip.tsv')
    tgas_xhip = tgas_xhip.sort_values('tgas.parallax')
    tgas = tgas_xhip['tgas.parallax']
    tgas_model = sm.add_constant(tgas)
    xhip = tgas_xhip['xhip.Plx']

    ols = sm.OLS(xhip, tgas_model)
    res = ols.fit()

    print res.summary()

    # xhip_predicted = res.predict(tgas_model)

    plt.scatter(tgas, xhip, s=1)
    # x0, y0 = tgas[0], xhip_predicted[0]
    # x1, y1 = tgas[len(tgas) - 1], xhip_predicted[len(xhip_predicted) - 1]

    # print x0, y0
    # print x1, y1
    # plt.plot([x0, y0], [x1, y1])
    # plt.bar(range(len(coeff)), coeff)
    plt.xlabel('TGAS parallax')
    plt.ylabel('XHIP parallax')
    # plt.xticks(range(len(coeff)), rotation='90')
    plt.show()