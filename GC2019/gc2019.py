import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import numpy as np

from common.gaia.with_rv import to_galaxy, prepare_galaxy_dataset, read_gaia_with_rv, read_raw_gaia_with_rv
import astropy.units as u
import astropy.coordinates as coord
from numpy import sin, cos

coord.galactocentric_frame_defaults.set('v4.0')

def read_catalog():
    data = []
    oname = []
    radeg = []
    dedeg = []
    dist = []
    hrv = []
    ehrv = []
    pmra = []
    pmde = []
    epmra = []
    epmde = []
    corr = []
    rscale = []
    nstar = []

    file_dir = Path(__file__).parent
    with open(file_dir.joinpath('catalog.dat'), 'r') as f:
        while True:
            name = str(f.read(11 - 1 + 1)).strip()
            if len(name) == 0:
                return pd.DataFrame({
                  'name': name,
                  'oname': oname,
                  'ra': radeg,
                  'dec': dedeg,
                  'dist': dist,
                    'parallax': np.float64(1000) / dist,
                  'radial_velocity': hrv,
                    'ehrv': ehrv,
                    'pmra': pmra,
                    'pmdec': pmde,
                    'pmra_error': epmra,
                    'pmdec_error': epmde,
                    'corr': corr,
                    'rscale': rscale,
                    'nstar': nstar
                })
                                    # # columns=['name', 'oname', 'ra', 'dec', 'dist', 'hrv', 'ehrv', 'pmra', 'pmdec',
                                    # #          'pmra_error', 'pmdec_error', 'corr', 'rscale', 'nstar'],
                                    # # dtype=[np.object_, np.object_, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64,
                                    # #        np.float64, np.float64, np.float64, np.float64, np.int32])
                                    # columns=['ra'],
                                    # dtype=np.dtype(['float']))
            f.read(1)
            oname.append(str(f.read(22 - 13 + 1)).strip())
            f.read(2)
            radeg.append(float(f.read(31 - 25 + 1)))
            f.read(1)
            dedeg.append(float(f.read(39 - 33 + 1)))
            f.read(3)
            try:
                dist.append(float(f.read(47 - 43 + 1)))
            except:
                dist.append(float('nan'))
            f.read(1)
            hrv.append(float(f.read(55 - 49 + 1)))
            f.read(3)
            ehrv.append(float(f.read(63 - 59 + 1)))
            f.read(1)
            pmra.append(float(f.read(71 - 65 + 1)))
            f.read(1)
            pmde.append(float(f.read(79 - 73 + 1)))
            f.read(3)
            epmra.append(float(f.read(87 - 83 + 1)))
            f.read(3)
            epmde.append(float(f.read(95 - 91 + 1)))
            f.read(2)
            corr.append(float(f.read(103 - 98 + 1)))
            f.read(4)
            rscale.append(float(f.read(111 - 108 + 1)))
            f.read(3)
            nstar.append(int(f.read(119 - 115 + 1)))
            f.read(1)


def absolute_velocity(df):
    R = np.float64(1000) / df.px
    k = 4.74
    return pd.DataFrame(data={
        'vx': -sin(df.l) * R * k * df['mul'] * cos(df.b) - sin(df.b) * cos(df.l) * R * k * df.mub + cos(df.b) * cos(df.l) * df.radial_velocity,
        'vy': cos(df.l) * R * k * df['mul'] * cos(df.b) - sin(df.b) * sin(df.l) * R * k * df.mub + cos(df.b) * sin(df.l) * df.radial_velocity,
        'vz': cos(df.b) - sin(df.b) * R * k * df.mub + sin(df.b) * df.radial_velocity
    })


def to_galactocentric_full(df):
    # df = prepare_galaxy_dataset(df)
    r2 = to_galctocentric_row(df)
    return pd.DataFrame(data=r2)


@dataclass(frozen=True)
class Spheric:
    ra: float
    dec: float
    pmra: float
    pmdec: float
    rv: float
    dist: float

    def __str__(self) -> str:
        return f"""ra={self.ra:.2f} deg
dec={self.dec:.2f} deg
pmra={self.pmra:.2f} mas/yr
pmdec={self.pmdec:.2f} mas/yr
rv={self.rv:.2f} km/s
dist={self.dist:.2f} kpc
"""


@dataclass(frozen=True)
class Galactocentric:
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    lon: float
    lat: float
    px: float
    mul: float
    mub: float
    rv: float

    def __str__(self) -> str:
        return f"""x={self.x / 1000:.2f} kpc
y={self.y / 1000:.2f} kpc
z={self.z / 1000:.2f} kpc
vx={self.vx:.2f} km/s
vy={self.vy:.2f} km/s
vz={self.vz:.2f} km/s
lon={self.lon:.2f} deg
lat={self.lat:.2f} deg
px={self.px:.2f} mas
mul={self.mul:.2f} mas/yr
mub={self.mub:.2f} mas/yr
radial_velocity={self.rv:.2f} km/s
"""


def to_galactocentric(ra, dec, pmra, pmdec, rv, dist):
    c = coord.SkyCoord(ra=np.array(ra) * u.deg,
                       dec=np.array(dec) * u.deg,
                       pm_ra_cosdec=np.array(pmra) * u.mas / u.yr,
                       pm_dec=np.array(pmdec) * u.mas / u.yr,
                       radial_velocity=np.array(rv) * u.km / u.s,
                       distance=coord.Distance(value=dist, unit=u.kiloparsec, allow_negative=True))

    g = c.transform_to(coord.Galactocentric)

    return Galactocentric(x=g.cartesian.x.value,
                          y=g.cartesian.y.value,
                          z=g.cartesian.z.value,
                          vx=g.velocity.d_x.value,
                          vy=g.velocity.d_y.value,
                          vz=g.velocity.d_z.value,
                          lon=g.spherical.lon.value,
                          lat=g.spherical.lat.value,
                          px=1000 / g.spherical.distance.value,
                          mul=g.proper_motion[0].value,
                          mub=g.proper_motion[1].value,
                          rv=g.radial_velocity.value)

def to_galctocentric_row(df):
    ra = df['ra']
    dec = df['dec']
    pmra = df['pmra']
    pmdec = df['pmdec']
    rv = df['radial_velocity']
    dist = df['dist']
    gal = to_galactocentric(ra, dec, pmra, pmdec, rv, dist)
    return {
        # 'ra': ra,
        # 'dec': dec,
        # 'pmra': pmra,
        # 'pmdec': pmdec,
        'radial_velocity': gal.rv,
        'galactocentric_x': gal.x,
        'galactocentric_y': gal.y,
        'galactocentric_z': gal.z,
        'galactocentric_vx': gal.vx,
        'galactocentric_vy': gal.vy,
        'galactocentric_vz': gal.vz,
        'galactocentric_lon': gal.lon,
        'galactocentric_lat': gal.lat,
        'galactocentric_px': gal.px,
        'galactocentric_mul': gal.mul,
        'galactocentric_mub': gal.mub,
        'galactocentric_radial_velocity': gal.rv,

        'px': gal.px,
        'l': np.radians(gal.lon),
        'b': np.radians(gal.lat),
        'mul': gal.mul,
        'mub': gal.mub
    }


def read_gc2019():
    df = read_catalog()
    df = df[(df.dist > 0)]
    # df = prepare_galaxy_dataset(df)
    return to_galactocentric_full(df)


def read_gaia_with_rv_in_galactocentric():
    df = read_raw_gaia_with_rv()
    df['dist'] = 1.0 / df['parallax']
    df = df[(df.dist > 0)]
    return to_galactocentric_full(df)


def main():
    #  1- 11  A11   ---       Name    Common name, NGC identifier, etc.
    #   13- 22  A10   ---       OName   Other name
    #   25- 31  F7.3  deg       RAdeg   Right Ascension J2000
    #   33- 39  F7.3  deg       DEdeg   Declination J2000
    #   43- 47  F5.1  kpc       Dist    ? Distance to the cluster
    #                                    (from Harris, 2010arXiv1012.3224H)
    #   49- 55  F7.2  km/s      HRV     Line-of-sight velocity (from Baumgardt et al.,
    #                                    2019, Cat. J/MNRAS/482/5138)
    #   59- 63  F5.2  km/s    e_HRV     Uncertainty in line-of-sight velocity
    #   65- 71  F7.3  mas/yr    pmRA    Proper motion in right ascension, pmRA*cosDE
    #   73- 79  F7.3  mas/yr    pmDE    Proper motion in declination
    #   83- 87  F5.3  mas/yr  e_pmRA    Uncertainty (including systematic) in pmRA
    #   91- 95  F5.3  mas/yr  e_pmDE    Uncertainty (including systematic) in pmDE
    #   98-103  F6.3  ---       corr    Correlation coefficient between uncertainties
    #  108-111  F4.1  arcmin    Rscale  Scale radius of Gaia-detected cluster members
    #  115-119  I5    ---       Nstar   Number of Gaia-detected cluster member stars
    df = read_gc2019()
    # df = df[(df.dist > 0)]
    # # df = prepare_galaxy_dataset(df)
    # start = time.time()
    # to_galctocentric_row(df)
    # for index, row in df.iterrows():
    #     print(to_galctocentric_row(row))
    # print(6000000 * (time.time() - start) / 100 / 3600)
    # df = read_gc2019()
    # df.to_csv('gc2019_with_galactocentric.tsv', sep='\t', encoding='utf-8')
    # for row in df.iterrows():

    # v = absolute_velocity(df)
    # ave_vy = np.average(v.vy)
    # v -= ave_vy
    # print(f'vx={np.average(v.vx):.2f}±{np.std(v.vx, ddof=1):.2f}')
    # print(f'vy={np.average(v.vy):.2f}±{np.std(v.vy, ddof=1):.2f}')
    # print(f'vy={np.average(v.vz):.2f}±{np.std(v.vz, ddof=1):.2f}')


if __name__ == '__main__':
    main()