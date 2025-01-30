import numpy
from astropy.coordinates import SkyCoord, Distance, Galactic, Galactocentric
from astropy import units as u
import pandas as pd
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

from common.gaia.with_rv import read_raw_gaia_with_rv_no_errors, read_gaia_dr3_with_rv_with_errors
from numpy import sin, cos
import numpy as np

from rv.degree_table import *


pizza = None

def read_and_filter():
    global pizza
    if pizza is not None:
        return pizza

    df = read_gaia_dr3_with_rv_with_errors()
    df = df[(df.parallax > 0)]
    # df = df[(df.parallax > 1.0 / 7.013) & (np.abs(df.parallax_error / df.parallax) < 0.2)]
    # df = df[(df.parallax > 0) & (np.abs(df.parallax_error / df.parallax) < numpy.sqrt(5) / 9)]
    df = df.sort_values(by='parallax', ascending=False)
    df = (df
          .head(30000000)
          # .sample(1000000)
          )
    print("got dataframe")
    print(len(df))
    print(np.max(1000 / df.parallax))

    pizza = df

    return df

def read_radec():
    df = read_and_filter()
    return SkyCoord(ra=df.ra.values * u.deg,
                 dec=df.dec.values * u.deg,
                 pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
                 pm_dec=df.pmdec.values * u.mas / u.yr,
                 radial_velocity=df.radial_velocity.values * u.km / u.s,
                 distance=Distance(parallax=df.parallax.values * u.mas))

def read_galactic():
    c = read_radec()
    g = c.transform_to(Galactic)
    return g


def convert(g, mul, mub, vr):
    l = np.deg2rad(np.array(g.l))
    b = np.deg2rad(np.array(g.b))
    r = np.array(g.distance) / 1000
    # print('pizza')
    # print(mul(0, 0, 1))
    # print(mub(0, 0, 1))
    # print(vr(0, 0, 1))
    # print(vr(0, 0, 1))
    return SkyCoord(frame="galactic", l=l * u.rad, b=b * u.rad, pm_l_cosb=mul(l, b, r) * u.mas / u.yr,
                    pm_b=mub(l, b, r) * u.mas / u.yr, radial_velocity=vr(l, b, r) * u.km / u.s, distance=g.distance)


def to_cartesian(g, pc=False):
    return pd.DataFrame(dict(x=g.cartesian.x.value / (1000 if pc else 1),
                             y=g.cartesian.y.value / (1000 if pc else 1),
                             z=g.cartesian.z.value / (1000 if pc else 1),
                             vx=g.velocity.d_x.value,
                             vy=g.velocity.d_y.value,
                             vz=g.velocity.d_z.value))

def read_gaia_with_rv_cartesian():
    g = read_galactic()
    df = to_cartesian(g, pc=True)
    # df = df[(abs(df.x) < 1)]
    # df.reset_index(drop=True, inplace=True)
    # print('pasta')
    # print(df)
    return df


def read_gaia_with_rv_residuals(degree):
    df = read_gaia_with_rv_cartesian()

    X = df[['x', 'y', 'z']]

    p = PolynomialFeatures(degree=degree).fit(X)

    X = pd.DataFrame(p.transform(X), columns=p.get_feature_names(X.columns))

    model = LinearRegression(fit_intercept=False)

    return pd.DataFrame(dict(x=df.x,
                             y=df.y,
                             z=df.z,
                             vx=(df.vx - model.fit(X, df.vx).predict(X)).values,
                             vy=(df.vy - model.fit(X, df.vy).predict(X)).values,
                             vz=(df.vz - model.fit(X, df.vz).predict(X)).values))


def read_gaia_with_rv_cubic_residuals():
    df = read_gaia_with_rv_cartesian()

    X = df[['x', 'y', 'z']]

    model = linear_model.LinearRegression()

    return pd.DataFrame(dict(x=df.x,
                             y=df.y,
                             z=df.z,
                             vx=(df.vx - model.fit(X, df.vx).predict(X)).values,
                             vy=(df.vy - model.fit(X, df.vy).predict(X)).values,
                             vz=(df.vz - model.fit(X, df.vz).predict(X)).values))

# VectorPlot[{-0.67 * x^2 + 2.30 * x * y, -1.79 * x^2 + -1.83 * y^2 }, {x, -3, 3}, {y, -3, 3}]

K = 4.74


def apply(P, f):
    c = SkyCoord(x=P.x, y=P.y, z=P.z, unit='kpc', representation_type='cartesian')
    c.representation_type = 'spherical'
    l = np.deg2rad(np.array(c.ra))
    b = np.deg2rad(np.array(c.dec))
    r = np.array(c.distance)
    res = f(l, b, r)
    mul_cosb = res["kmul"] / K
    mub = res["kmub"] / K
    vr = res["vr_r"] * r
    # print(mul_cosb)
    # print(mub)
    # print(vr)
    sc = SkyCoord(frame="galactic",
                  l=l * u.rad,
                  b=b * u.rad,
                  pm_l_cosb=mul_cosb * u.mas / u.yr,
                  pm_b=mub * u.mas / u.yr,
                  radial_velocity=vr * u.km / u.s,
                  distance=c.galactic.distance)
    return to_cartesian(sc)

def U(l, b, r):
    return dict(kmul=sin(l) / r,
                kmub=cos(l) * sin(b) / r,
                vr_r=-cos(b) * cos(l) / r)


def V(l, b, r):
    return dict(kmul=-cos(l) / r,
                kmub=sin(l) * sin(b) / r,
                vr_r=-cos(b) * sin(l) / r)


def W(l, b, r):
    return dict(kmul=0,
                kmub=-cos(b) / r,
                vr_r=-sin(b) / r)


def w1(l, b, r):
    return dict(kmul=-sin(b) * cos(l),
                kmub=sin(l),
                vr_r=0)


def w2(l, b, r):
    return dict(kmul=-sin(b) * sin(l),
                kmub=-cos(l),
                vr_r=0)


def w3(l, b, r):
    return dict(kmul=+cos(b),
                kmub=0,
                vr_r=0)


def M12(l, b, r):
    return dict(kmul=+cos(b) * cos(2 * l),
                kmub=-0.5 * sin(2 * b) * sin(2 * l),
                vr_r=cos(b) ** 2 * sin(2 * l))


def M13(l, b, r):
    return dict(kmul=-sin(b) * sin(l),
                kmub=+cos(2 * b) * cos(l),
                vr_r=sin(2 * b) * cos(l))


def M23(l, b, r):
    return dict(kmul=+sin(b) * cos(l),
                kmub=+cos(2 * b) * sin(l),
                vr_r=sin(2 * b) * sin(l))


def M11(l, b, r):
    return dict(kmul=-0.5 * cos(b) * sin(2 * l),
                kmub=-0.5 * sin(2 * b) * cos(l) ** 2,
                vr_r=cos(b) ** 2 * cos(l) ** 2)


def M22(l, b, r):
    return dict(kmul=+0.5 * cos(b) * sin(2 * l),
                kmub=-0.5 * sin(2 * b) * sin(l) ** 2,
                vr_r=cos(b) ** 2 * sin(l) ** 2)


def M33(l, b, r):
    return dict(kmul=0,
                kmub=+0.5 * sin(2 * b),
                vr_r=sin(b) ** 2)


def dw1dr1(l, b, r):
    return dict(kmul=r * -mul_dw1dr1.apply(b, l),
                kmub=r * mub_dw1dr1.apply(b, l),
                vr_r=0)

def dw1dr2(l, b, r):
    return dict(kmul=r * -mul_dw1dr2.apply(b, l),
                kmub=r * mub_dw1dr2.apply(b, l),
                vr_r=0)


def dw1dr3(l, b, r):
    return dict(kmul=r * -mul_dw1dr3.apply(b, l),
                kmub=r * mub_dw1dr3.apply(b, l),
                vr_r=0)


def dw2dr1(l, b, r):
    return dict(kmul=r * -mul_dw2dr1.apply(b, l),
                kmub=r * -mub_dw2dr1.apply(b, l),
                vr_r=0)


def dw2dr2(l, b, r):
    return dict(kmul=r * -mul_dw2dr2.apply(b, l),
                kmub=r * -mub_dw2dr2.apply(b, l),
                vr_r=0)


def dw2dr3(l, b, r):
    return dict(kmul=r * -mul_dw2dr3.apply(b, l),
                kmub=r * -mub_dw2dr3.apply(b, l),
                vr_r=0)


def dw3dr1(l, b, r):
    return dict(kmul=r * mul_dw3dr1.apply(b, l),
                kmub=0,
                vr_r=0)


def dw3dr2(l, b, r):
    return dict(kmul=r * mul_dw3dr2.apply(b, l),
                kmub=0,
                vr_r=0)


def dw3dr3(l, b, r):
    return dict(kmul=r * mul_dw3dr3.apply(b, l),
                kmub=0,
                vr_r=0)


def dM11dr1(l, b, r):
    return dict(kmul=r * mul_dM11dr1.apply(b, l),
                kmub=r * mub_dM11dr1.apply(b, l),
                vr_r=r * vr_dM11dr1.apply(b, l))


def dM11dr2(l, b, r):
    return dict(kmul=r * mul_dM11dr2.apply(b, l),
                kmub=r * mub_dM11dr2.apply(b, l),
                vr_r=r * vr_dM11dr2.apply(b, l))


def dM11dr3(l, b, r):
    return dict(kmul=r * mul_dM11dr3.apply(b, l),
                kmub=r * mub_dM11dr3.apply(b, l),
                vr_r=r * vr_dM11dr3.apply(b, l))


# !!!
def dM12dr1(l, b, r):
    return dict(kmul=r * mul_dM12dr1.apply(b, l),
                kmub=r * mub_dM12dr1.apply(b, l),
                vr_r=r * 2 * vr_dM12dr1.apply(b, l))

# !!!
def dM12dr2(l, b, r):
    return dict(kmul=r * mul_dM12dr2.apply(b, l),
                kmub=r * mub_dM12dr2.apply(b, l),
                vr_r=r * 2 * vr_dM12dr2.apply(b, l))

# !!!
def dM12dr3(l, b, r):
    return dict(kmul=r * mul_dM12dr3.apply(b, l),
                kmub=r * mub_dM12dr3.apply(b, l),
                vr_r=r * 2 * vr_dM12dr3.apply(b, l))


def dM13dr1(l, b, r):
    return dict(kmul=r * mul_dM13dr1.apply(b, l),
                kmub=r * mub_dM13dr1.apply(b, l),
                vr_r=r * 2 * vr_dM13dr1.apply(b, l))


def dM13dr2(l, b, r):
    return dict(kmul=r * mul_dM13dr2.apply(b, l),
                kmub=r * mub_dM13dr2.apply(b, l),
                vr_r=r * 2 * vr_dM13dr2.apply(b, l))


def dM13dr3(l, b, r):
    return dict(kmul=r * mul_dM13dr3.apply(b, l),
                kmub=r * mub_dM13dr3.apply(b, l),
                vr_r=r * 2 * vr_dM13dr3.apply(b, l))


def dM22dr1(l, b, r):
    return dict(kmul=r * mul_dM22dr1.apply(b, l),
                kmub=r * mub_dM22dr1.apply(b, l),
                vr_r=r * vr_dM22dr1.apply(b, l))


def dM22dr2(l, b, r):
    return dict(kmul=r * mul_dM22dr2.apply(b, l),
                kmub=r * mub_dM22dr2.apply(b, l),
                vr_r=r * vr_dM22dr2.apply(b, l))


def dM22dr3(l, b, r):
    return dict(kmul=r * mul_dM22dr3.apply(b, l),
                kmub=r * mub_dM22dr3.apply(b, l),
                vr_r=r * vr_dM22dr3.apply(b, l))


def dM23dr1(l, b, r):
    return dict(kmul=r * mul_dM23dr1.apply(b, l),
                kmub=r * mub_dM23dr1.apply(b, l),
                vr_r=r * 2 * vr_dM23dr1.apply(b, l))


def dM23dr2(l, b, r):
    return dict(kmul=r * mul_dM23dr2.apply(b, l),
                kmub=r * mub_dM23dr2.apply(b, l),
                vr_r=r * 2 * vr_dM23dr2.apply(b, l))


def dM23dr3(l, b, r):
    return dict(kmul=r * mul_dM23dr3.apply(b, l),
                kmub=r * mub_dM23dr3.apply(b, l),
                vr_r=r * 2 * vr_dM23dr3.apply(b, l))


def dM33dr1(l, b, r):
    return dict(kmul=0,
                kmub=r * mub_dM33dr1.apply(b, l),
                vr_r=r * vr_dM33dr1.apply(b, l))

def dM33dr2(l, b, r):
    return dict(kmul=0,
                kmub=r * mub_dM33dr2.apply(b, l),
                vr_r=r * vr_dM33dr2.apply(b, l))


def dM33dr3(l, b, r):
    return dict(kmul=0,
                kmub=r * mub_dM33dr3.apply(b, l),
                vr_r=r * vr_dM33dr3.apply(b, l))
