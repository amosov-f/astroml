from time import time

# import pandas as pd
import numpy
import astropy.coordinates as coord
import astropy.units as u

import yt.wrapper as yt



def main():
    # yt.config["pickling"]["python_binary"] = "python3.4"
    yt.config['proxy']['url'] = 'freud'

    yt.run_map(
        'python3 to_cartesian.py',
        source_table='//home/logsng/amosov-f/gaia/gdr2/with_rv',
        destination_table='//home/logsng/amosov-f/gaia/gdr2/cartesian',
        local_files='to_cartesian.py',
        # yt_files='//home/cocainum/novichok/scripts/defect_detector.py',
        input_format='json',
        output_format='json',
        # memory_limit=yt.common.GB,
        spec={
            'mapper': {
                'layer_paths': [
                    '//home/logsng/amosov-f/layers/astropy.tar.gz',
                ]
            },
            'scheduling_tag_filter': 'porto',
            'job_count': 50_000,
        }
    )

    return

    table = pd.read_csv('gdr2_with_rv.tsv', sep='\t')

    for i, row in table.iterrows():

        print(time())
        c = coord.SkyCoord(ra=row['ra'] * u.deg,
                                         dec=row['dec'] * u.deg,
                                         pm_ra_cosdec=row['pmra'] * u.mas/u.yr,
                                         pm_dec=row['pmdec'] * u.mas/u.yr,
                                         radial_velocity=row['radial_velocity'] * u.km/u.s,
                                         distance=coord.Distance(parallax=row['parallax'] * u.mas, allow_negative=True))
        print(time())
        # print(c)
        # print()

        g = c.transform_to(coord.Galactic)
        # print(g)
        # print()
        #
        # print(g.cartesian)
        # print()
        #
        # print(g.velocity)
        # print()
        print('{}\t{}\t{}\t{}\t{}\t{}'.format(g.cartesian.x.value, g.cartesian.y.value, g.cartesian.z.value, g.velocity.d_x.value, g.velocity.d_y.value, g.velocity.d_z.value))

    # c = coord.CartesianRepresentation.ICRSCoordinates(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree), distance=coord.Distance(770, u.kpc))

    # cart = coord.ICRS(g, )
    # print(cart)

    # print(astropy.coordinates.SkyCoord(row['pmra'] * u.mas, row['pmdec'] * u.rad, radial_velocity=row['radial_velocity']))



if __name__ == '__main__':
    main()