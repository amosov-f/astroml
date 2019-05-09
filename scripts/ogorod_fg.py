import pandas as pd

from common.gaia.with_rv import read_gaia_with_rv
from common.linear_decomposer import decompose
from numpy import sin, cos, sqrt, radians


K = 4.74


def ogorod_fg_mul(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    return pd.DataFrame(dict(U=sin(l) / r,
                             V=-cos(l) / r,
                             A=cos(b) * cos(2 * l),
                             B=cos(b),
                             F=-r * cos(b) ** 2 * cos(l) ** 3,
                             G=-r * (3 * cos(b) ** 2 * cos(l) - cos(b) ** 2 * cos(l) ** 3),
                             y=K * mul * cos(b)))


def ogorod_fg_mub(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    return pd.DataFrame(dict(U=cos(l) * sin(b) / r,
                             V=sin(l) * sin(b) / r,
                             W=-cos(b) / r,
                             A=-cos(b) * sin(b) * sin(2 * l),
                             F=r * cos(b) ** 2 * sin(b) * sin(l) * cos(l) ** 2,
                             G=r * cos(b) ** 2 * sin(b) * sin(l) ** 3,
                             K=-cos(b) * sin(b),
                             y=K * mub))


def ogorod_fg_vr(df):
    l, b, mul, mub, px = columns(df)
    r = 1 / px
    return pd.DataFrame(dict(U=-cos(b) * cos(l) / r,
                             V=-cos(b) * sin(l) / r,
                             W=-sin(b) / r,
                             A=cos(b) ** 2 * sin(2 * l),
                             F=-r * cos(b) ** 3 * sin(l) * cos(l) ** 2,
                             G=-r * cos(b) ** 3 * sin(l) ** 3,
                             K=cos(b) ** 2,
                             y=df.radial_velocity / r))


def columns(df: pd.DataFrame):
    return df.l, df.b, df['mul'], df.mub, df.px


def compute_ogorod_fg(df):
    # mul = ogorod_fg_mul(df)
    # mub = ogorod_fg_mub(df)
    vr = ogorod_fg_vr(df)
    return vr


STEP = 400000
SLICE = STEP
MAX = 6000000


def main():
    dataset = read_gaia_with_rv()
    sort = [
        ['U', 'V', 'W'],
        ['A'],
        ['F', 'G', 'K']
    ]
    decompose(dataset, compute_ogorod_fg, MAX, STEP, SLICE, sort)


if __name__ == '__main__':
    main()
