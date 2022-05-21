from astropy.coordinates import SkyCoord, Distance, Galactic
from astropy import units as u
import pandas as pd

from common.gaia.with_rv import read_raw_gaia_with_rv_no_errors
from numpy import sin, cos
import numpy as np


def read_pizza():
    df = read_raw_gaia_with_rv_no_errors()
    # df = df[(df.parallax > 0)]
    df = df[(df.parallax > 0.33)]
    print(len(df))
    c = SkyCoord(ra=df.ra.values * u.deg,
                 dec=df.dec.values * u.deg,
                 pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
                 pm_dec=df.pmdec.values * u.mas / u.yr,
                 radial_velocity=df.radial_velocity.values * u.km / u.s,
                 distance=Distance(parallax=df.parallax.values * u.mas))
    g = c.transform_to(Galactic)
    return g

def convert(g, mul, mub, vr):
    l = np.array(g.l)
    b = np.array(g.b)
    r = np.array(g.distance) / 1000
    print(r)
    return SkyCoord(frame="galactic", l=l * u.deg, b=b * u.deg, pm_l_cosb=mul(l, b, r) * u.mas / u.yr, pm_b=mub(l, b, r) * u.mas / u.yr, radial_velocity=vr(l, b, r) * u.km / u.s, distance=Distance(value=r, unit=u.kpc))

def to_cartesian(g):
    return pd.DataFrame(dict(x=g.cartesian.x.value / 1000,
                             y=g.cartesian.y.value / 1000,
                             z=g.cartesian.z.value / 1000,
                             vx=g.velocity.d_x.value,
                             vy=g.velocity.d_y.value,
                             vz=g.velocity.d_z.value))

def read_gaia_with_rv_cartesian():
    g = read_pizza()
    return to_cartesian(g)

# VectorPlot[{-0.67 * x^2 + 2.30 * x * y, -1.79 * x^2 + -1.83 * y^2 }, {x, -3, 3}, {y, -3, 3}]

def U(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: sin(l) / r / 4.74,
                                       mub=lambda l, b, r: cos(l) * sin(b) / r / 4.74,
                                       vr= lambda l, b, r: -cos(b) * cos(l)))

def V(g):
    r = g.distance
    return to_cartesian(convert(g, mul=lambda l, b, r: -cos(l) / r / 4.74,
                                       mub=lambda l, b, r: sin(l) * sin(b) / r / 4.74,
                                       vr= lambda l, b, r: -cos(b) * sin(l)))

def W(g):
    r = g.distance
    return to_cartesian(convert(g, mul=lambda l, b, r: -sin(b) * cos(l)  / 4.74,
                                       mub=lambda l, b, r: -cos(b) / r  / 4.74,
                                       vr= lambda l, b, r: -sin(b)))

def w1(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: 0  / 4.74,
                                       mub=lambda l, b, r: sin(l)  / 4.74,
                                       vr= lambda l, b, r: 0))

def w2(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: -sin(b) * sin(l)  / 4.74,
                                       mub=lambda l, b, r: -cos(l)  / 4.74,
                                       vr= lambda l, b, r: 0))

def w3(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: +cos(b)  / 4.74,
                                       mub=lambda l, b, r: 0  / 4.74,
                                       vr= lambda l, b, r: 0))

def M12(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: +cos(b) * cos(2 * l)  / 4.74,
                                       mub=lambda l, b, r: -0.5 * sin(2 * b) * sin(2 * l)  / 4.74,
                                       vr= lambda l, b, r: cos(b) ** 2 * sin(2 * l) * r))

def M13(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: -sin(b) * sin(l)  / 4.74,
                                       mub=lambda l, b, r: +cos(2 * b) * cos(l)  / 4.74,
                                       vr= lambda l, b, r: sin(2 * b) * cos(l) * r))

def M23(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: +sin(b) * cos(l)  / 4.74,
                                       mub=lambda l, b, r: +cos(2 * b) * sin(l)  / 4.74,
                                       vr= lambda l, b, r: sin(2 * b) * sin(l) * r))

def M11(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: -0.5 * cos(b) * sin(2 * l)  / 4.74,
                                       mub=lambda l, b, r: -0.5 * sin(2 * b) * cos(l) ** 2  / 4.74,
                                       vr= lambda l, b, r: cos(b) ** 2 * cos(l) ** 2 * r))

def M22(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: +0.5 * cos(b) * sin(2 * l)  / 4.74,
                                       mub=lambda l, b, r: -0.5 * sin(2 * b) * sin(l) ** 2  / 4.74,
                                       vr= lambda l, b, r: cos(b) ** 2 * sin(l) ** 2 * r))

def M33(g):
    return to_cartesian(convert(g, mul=lambda l, b, r: 0  / 4.74,
                                       mub=lambda l, b, r: +0.5 * sin(2 * b)  / 4.74,
                                       vr= lambda l, b, r: sin(b) ** 2 * r))