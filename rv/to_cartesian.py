import sys

import json

import astropy.coordinates as coord
import astropy.units as u


def main():
    for line in sys.stdin:
        row = json.loads(line)
        if not (row['pmra'] and row['pmdec']):
            continue
        c = coord.SkyCoord(ra=row['ra'] * u.deg,
                           dec=row['dec'] * u.deg,
                           pm_ra_cosdec=row['pmra'] * u.mas / u.yr,
                           pm_dec=row['pmdec'] * u.mas / u.yr,
                           radial_velocity=row['radial_velocity'] * u.km / u.s,
                           distance=coord.Distance(parallax=row['parallax'] * u.mas, allow_negative=True))

        g = c.transform_to(coord.Galactic)
        res = {
            'x': g.cartesian.x.value,
            'y': g.cartesian.y.value,
            'z': g.cartesian.z.value,
            'vx': g.velocity.d_x.value,
            'vy': g.velocity.d_y.value,
            'vz': g.velocity.d_z.value
        }
        print(json.dumps(res))

def to_sky():
    pasta = coord.Galactic(representation_type=coord.CartesianRepresentation,
                           u=-10.34848874410621*u.pc,
                           v=-17.554795954741007*u.pc,
                           w=-7.233114343092995*u.pc)
    print(pasta.transform_to(coord.ICRS))

if __name__ == '__main__':
    to_sky()