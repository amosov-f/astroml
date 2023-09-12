import math
import sys

import numpy
import numpy as np
import pandas as pd

from common.gaia.with_rv import read_gaia_with_rv_full
from common.linear_decomposer import decompose
from numpy import sin, cos, sqrt, radians

from rv.degree_table import  \
    mul_dw1dr1, mul_dw1dr2, mul_dw1dr3, mul_dw2dr1, mul_dw2dr2, mul_dw2dr3, mul_dw3dr1, mul_dw3dr2, mul_dw3dr3, mul_dM11dr1, mul_dM11dr2, \
    mul_dM11dr3, mul_dM12dr1, mul_dM12dr2, mul_dM12dr3, mul_dM13dr1, mul_dM13dr2, mul_dM13dr3, mul_dM22dr1, mul_dM22dr2, \
    mul_dM22dr3, mul_dM23dr1, mul_dM23dr2, mul_dM23dr3, mub_dw1dr1, mub_dw1dr2, mub_dw1dr3, mub_dw2dr1, mub_dw2dr2, mub_dw2dr3, \
    mub_dM11dr1, mub_dM11dr2, mub_dM11dr3, mub_dM12dr1, mub_dM12dr2, mub_dM12dr3, mub_dM13dr1, mub_dM13dr2, mub_dM13dr3, \
    mub_dM22dr1, mub_dM22dr2, mub_dM22dr3, mub_dM23dr1, mub_dM23dr2, mub_dM23dr3, mub_dM33dr1, mub_dM33dr2, mub_dM33dr3, \
    vr_dM11dr1, vr_dM11dr2, vr_dM11dr3, vr_dM22dr1, vr_dM22dr2, vr_dM22dr3, vr_dM33dr1, vr_dM33dr2, vr_dM33dr3, vr_dM12dr1, \
    vr_dM12dr2, vr_dM12dr3, vr_dM13dr1, vr_dM13dr2, vr_dM13dr3, vr_dM23dr1, vr_dM23dr2, vr_dM23dr3

K = 4.74
numpy.set_printoptions(threshold=sys.maxsize)

ct211 = (5 / (8 * math.pi)) ** 0.5
cs310 = - (7 / (128 * math.pi)) ** 0.5
cv310 = (21 / (32 * math.pi)) ** 0.5

def ogorod_fg_mul(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    # R1 = cos(b) * cos(l)
    # R2 = cos(b) * sin(l)
    # R3 = sin(b)

    # la1 = sin(b) * cos(l)
    # la2 = sin(b) * sin(l)
    # la3 = cos(b)
    #
    # p11 = -0.5 * cos(b) * sin(2 * l)
    # p12 = cos(b) * cos(2 * l)
    # p13 = -sin(b) * sin(l)
    #
    # p22 = 0.5 * cos(b) * sin(2 * l)
    # p23 = sin(b) * cos(l)
    #
    # p33 = 0

    return pd.DataFrame(dict(U=sin(l) / r,
                             V=-cos(l) / r,
                             A=cos(b) * cos(2 * l),
                             B=cos(b),
                             # w1=-sin(b) * cos(l),
                             # w2=-sin(b) * sin(l),
                             # w3=+cos(b),
                             F=-r * cos(b) ** 2 * cos(l) ** 3,
                             G=-r * (3 * cos(b) ** 2 * cos(l) - cos(b) ** 2 * cos(l) ** 3),
                             # D=r * cos(b) ** 2 * cos(l) * cos(2 * l),

                             # M12=+cos(b) * cos(2 * l),
                             # M13=-sin(b) * sin(l),
                             # M23=+sin(b) * cos(l),
                             # M11=-0.5 * cos(b) * sin(2 * l),
                             # M22=+0.5 * cos(b) * sin(2 * l),
                             # M33=0.0,

                             # C=-cos(b) * sin(2 * l),
                             # K=0.0,

                             # t211=ct211 * cos(2 * b) * cos(l),
                             # s310=cs310 * (-(5 * sin(b) ** 2 - 1) * cos(l)),

                             # new_mul=cos(b) ** 2 * cos(l),


                             # dw1dr1=r * -mul_dw1dr1.apply(b, l),  # 0
                             # dw1dr2=r * -mul_dw1dr2.apply(b, l),
                             # dw1dr3=r * -mul_dw1dr3.apply(b, l),  # !

                             # dw2dr1=r * -mul_dw2dr1.apply(b, l),
                             # dw2dr2=r * -mul_dw2dr2.apply(b, l),
                             # dw2dr3=r * -mul_dw2dr3.apply(b, l),

                             # dw3dr1=r * mul_dw3dr1.apply(b, l),  # !
                             # dw3dr2=r * mul_dw3dr2.apply(b, l),
                             # dw3dr3=r * mul_dw3dr3.apply(b, l),

                             # dM11dr1=r * mul_dM11dr1.apply(b, l),
                             # dM11dr2=r * mul_dM11dr2.apply(b, l),  # !
                             # dM11dr3=r * mul_dM11dr3.apply(b, l),

                             # dM12dr1=r * mul_dM12dr1.apply(b, l),  # !
                             # dM12dr2=r * mul_dM12dr2.apply(b, l),  # !
                             # dM12dr3=r * mul_dM12dr3.apply(b, l),

                             # dM13dr1=r * mul_dM13dr1.apply(b, l),
                             # dM13dr2=r * mul_dM13dr2.apply(b, l),
                             # dM13dr3=r * mul_dM13dr3.apply(b, l),

                             # dM22dr1=r * mul_dM22dr1.apply(b, l),
                             # dM22dr2=r * mul_dM22dr2.apply(b, l),  # !
                             # dM22dr3=r * mul_dM22dr3.apply(b, l),  # !

                             # dM23dr1=r * mul_dM23dr1.apply(b, l),
                             # dM23dr2=r * mul_dM23dr2.apply(b, l),  # !
                             # dM23dr3=r * mul_dM23dr3.apply(b, l),  # !

                             # zero!!!
                             # dM33dr1=r * (R1 * p33).apply(b, l),
                             # dM33dr2=r * (R2 * p33).apply(b, l), # !
                             # dM33dr3=r * (R3 * p33).apply(b, l),


                             y=K * mul))


def ogorod_fg_mub(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px

    # R1 = cos(b) * cos(l)
    # R2 = cos(b) * sin(l)
    # R3 = sin(b)

    # ta1 = sin(l)
    # ta2 = cos(l)
    #
    # q11 = -0.5 * sin(2 * b) * (cos(l) ** 2)
    # q12 = -0.5 * sin(2 * b) * cos(2 * l)
    # q13 = cos(2 * b) * cos(l)
    #
    # q22 = -0.5 * sin(2 * b) * (sin(l) ** 2)
    # q23 = cos(2 * b) * sin(l)
    # q33 = 0.5 * sin(2 * b)

    return pd.DataFrame(dict(U=cos(l) * sin(b) / r,
                             V=sin(l) * sin(b) / r,
                             W=-cos(b) / r,
                             A=-cos(b) * sin(b) * sin(2 * l),
                             # w1=sin(l),
                             # w2=-cos(l),
                             # w3=0,#cos(l),  # 0?, # TODO!!!
                             F=r * cos(b) ** 2 * sin(b) * sin(l) * cos(l) ** 2,
                             G=r * cos(b) ** 2 * sin(b) * sin(l) ** 3,
                             K=-cos(b) * sin(b),
                             # D=- 0.5 * r * cos(b) * sin(2 * b) * cos(l) * cos(2 * l),

                             # M12=-0.5 * sin(2 * b) * sin(2 * l),
                             # M13=+cos(2 * b) * cos(l),
                             # M23=+cos(2 * b) * sin(l),
                             # M11=-0.5 * sin(2 * b) * cos(l) ** 2,
                             # M22=-0.5 * sin(2 * b) * sin(l) ** 2,
                             # M33=+0.5 * sin(2 * b),

                             # C=-0.5 * sin(2 * b) * cos(2 * l),
                             # K=-0.5 * sin(2 * b),

                             # t211=ct211 * sin(b) * sin(l),
                             # s310=cs310 * ((15 * sin(b) ** 2 - 11) * sin(b) * sin(l)),

                             # new_mub=sin(b) * cos(b) ** 2 * sin(l),

                             # dw1dr1=r * mub_dw1dr1.apply(b, l),  # 0
                             # dw1dr2=r * mub_dw1dr2.apply(b, l),
                             # dw1dr3=r * mub_dw1dr3.apply(b, l),  # !

                             # dw2dr1=r * -mub_dw2dr1.apply(b, l),
                             # dw2dr2=r * -mub_dw2dr2.apply(b, l),
                             # dw2dr3=r * -mub_dw2dr3.apply(b, l),

                             # zero!!!
                             # dw3dr1=r * 0,  # !
                             # dw3dr2=r * 0,
                             # dw3dr3=r * 0,

                             # dM11dr1=r * mub_dM11dr1.apply(b, l),
                             # dM11dr2=r * mub_dM11dr2.apply(b, l),  # !
                             # dM11dr3=r * mub_dM11dr3.apply(b, l),

                             # dM12dr1=r * mub_dM12dr1.apply(b, l),  # !
                             # dM12dr2=r * mub_dM12dr2.apply(b, l),  # !
                             # dM12dr3=r * mub_dM12dr3.apply(b, l),

                             # dM13dr1=r * mub_dM13dr1.apply(b, l),
                             # dM13dr2=r * mub_dM13dr2.apply(b, l),
                             # dM13dr3=r * mub_dM13dr3.apply(b, l),

                             # dM22dr1=r * mub_dM22dr1.apply(b, l),
                             # dM22dr2=r * mub_dM22dr2.apply(b, l),  # !
                             # dM22dr3=r * mub_dM22dr3.apply(b, l),  # !

                             # dM23dr1=r * mub_dM23dr1.apply(b, l),
                             # dM23dr2=r * mub_dM23dr2.apply(b, l),  # !
                             # dM23dr3=r * mub_dM23dr3.apply(b, l),  # !

                             # dM33dr1=r * mub_dM33dr1.apply(b, l),
                             # dM33dr2=r * mub_dM33dr2.apply(b, l),  # !
                             # dM33dr3=r * mub_dM33dr3.apply(b, l),

                             y=K * mub))


def ogorod_fg_vr(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    # R1 = cos(b) * cos(l)
    # R2 = cos(b) * sin(l)
    # R3 = sin(b)
    return pd.DataFrame(dict(U=-cos(b) * cos(l) / r,
                             V=-cos(b) * sin(l) / r,
                             W=-sin(b) / r,
                             A=cos(b) ** 2 * sin(2 * l),
                             F=-r * cos(b) ** 3 * sin(l) * cos(l) ** 2 * r,
                             G=-r * cos(b) ** 3 * sin(l) ** 3 * r,
                             K=cos(b) ** 2 * r,
                             # D=2 * r * (cos(b) ** 3) * (cos(l) ** 2) * sin(l),

                             # M12=cos(b) ** 2 * sin(2 * l),
                             # M13=sin(2 * b) * cos(l),
                             # M23=sin(2 * b) * sin(l),
                             # M11=cos(b) ** 2 * cos(l) ** 2,
                             # M22=cos(b) ** 2 * sin(l) ** 2,
                             # M33=sin(b) ** 2,

                             # v310=cv310 * (5 * sin(b) ** 2 - 1) * cos(b) * sin(l),

                             # dM11dr1=r * vr_dM11dr1.apply(b, l),
                             # dM11dr2=r * vr_dM11dr2.apply(b, l),  # !
                             # dM11dr3=r * vr_dM11dr3.apply(b, l),

                             # dM22dr1=r * vr_dM22dr1.apply(b, l),
                             # dM22dr2=r * vr_dM22dr2.apply(b, l),  # !
                             # dM22dr3=r * vr_dM22dr3.apply(b, l),  # !

                             # dM33dr1=r * vr_dM33dr1.apply(b, l),
                             # dM33dr2=r * vr_dM33dr2.apply(b, l),  # !
                             # dM33dr3=r * vr_dM33dr3.apply(b, l),

                             # dM12dr1=r * 2 * vr_dM12dr1.apply(b, l),  # !
                             # dM12dr2=r * 2 * vr_dM12dr2.apply(b, l),  # !
                             # dM12dr3=r * 2 * vr_dM12dr3.apply(b, l),

                             # dM13dr1=r * 2 * vr_dM13dr1.apply(b, l),
                             # dM13dr2=r * 2 * vr_dM13dr2.apply(b, l),
                             # dM13dr3=r * 2 * vr_dM13dr3.apply(b, l),

                             # dM23dr1=r * 2 * vr_dM23dr1.apply(b, l),
                             # dM23dr2=r * 2 * vr_dM23dr2.apply(b, l),  # !
                             # dM23dr3=r * 2 * vr_dM23dr3.apply(b, l),  # !

                             y=df.radial_velocity / r))


def columns(df: pd.DataFrame):
    return df.l, df.b, df['mul'], df.mub, df.px


def compute_ogorod_fg(df):
    mul = ogorod_fg_mul(df)
    mub = ogorod_fg_mub(df)
    vr = ogorod_fg_vr(df)
    return pd.concat([mul, mub, vr])
    # return vr

def compute_mu(df):
    mul = ogorod_fg_mul(df)
    mub = ogorod_fg_mub(df)
    return pd.concat([mul, mub])


STEP = 2000000
SLICE = STEP
MAX = STEP * 15


def main():
    dataset = read_gaia_with_rv_full()#read_gaia_with_rv_1500()
    print('read')
    # df = df
    sort = [
        ['U', 'V', 'W'],
        ['A', 'B'],

        # ['w1', 'w2', 'w3'],
        # ['M12', 'M13', 'M23', 'M11', 'M22', 'M33'],
        ['F', 'G', 'K'],

        ['dw1dr3', 'dw3dr1'],
        ['dM11dr2', 'dM12dr1', 'dM12dr2', 'dM22dr3', 'dM23dr2', 'dM33dr2', 'dM22dr2', 'dM23dr3'],
        ['t101', 's110', 't211', 's310', 'v310', 'dist']
        # ['new_mul', 'new_mub']
        # ['F', 'G', 'K']
    ]
    # c = [[], [], [], []]
    # decompose(dataset, ogorod_fg_vr, MAX, STEP, SLICE, sort)
    decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort)
    # decompose(dataset, compute_mu, MAX, STEP, SLICE, sort)
    # decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort, filter=lambda df: np.abs(df.b) < math.radians(1))
    # c[0], average_dists = decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort, filter=lambda df: (((0 < df.l) & (df.l < + math.pi / 4)) | ((7 * math.pi / 4 < df.l) & (df.l < 2 * math.pi))))
    # c[1], _ = decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort, filter=lambda df: (+math.pi / 4 < df.l) & (df.l < 3 * math.pi / 4))
    # c[2], _ = decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort, filter=lambda df: (3 * math.pi / 4 < df.l) & (df.l < 5 * math.pi / 4))
    # c[3], _ = decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort, filter=lambda df: (5 * math.pi / 4 < df.l) & (df.l < 7 * math.pi / 4))
    # c[4], _ = decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort)

    # print('\t' + '\t'.join(average_dists))
    # for figure in sort:
    #     for label in figure:
    #         for i in range(len(c)):
    #             print(label + (str(i) if i < 4 else ''), end='\t')
    #             for val in c[i][label]:
    #                 print('{0:.1f}'.format(val).replace('.', ','), end='\t')
    #             print()

    # As = to_center['A']
    # for i in range(len(As) - 1):
    #     A1 = As[i]
    #     A2 = As[i + 1]
    #     r1 = float(average_dists[i]) / 1000
    #     r2 = float(average_dists[i + 1]) / 1000
    #     dA = (A2 - A1) / (r2 - r1)
    #     r = r2 # (r1 + r2) / 2
    #     s310 = -0.26 * dA * r
    #     v310 = -0.22 * dA * r
    #     t211 = 0.37 * dA * r
    #     print(f"{r}\t{s310}\t{v310}\t{t211}".replace('.', ','))

def pasta():
    a = """-0.8
-1.9
-2.5
-3.2
-3.6
-4.5
-5.3
-8.5
-15.0
-20.4
-26.1
-32.1
-38.7
-46.5
-53.4
"""
    pp = [p.strip() for p in a.split("\n")]
    for p in pp:
        print(p.replace('.', ','))


if __name__ == '__main__':
    main()
