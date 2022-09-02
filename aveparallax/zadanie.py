from common.gaia.with_rv import read_gaia_with_rv_full, distance
import healpy as hp
import numpy as np

from common.linear_decomposer import compute_coeffs
from scripts.ogorod_fg import compute_ogorod_fg, ogorod_fg_vr

import locale
# locale.setlocale(locale.LC_ALL, 'nl_NL')


def compute_pixel(nside, l, b):
    return hp.ang2pix(nside, (-b + np.pi / 2).values, l.values)


def solve(df, max_dist, width):
    print(f'{width}/{max_dist}')

    bins_count = max_dist // width

    df = df[(0 <= df.distance_pc) & (df.distance_pc < max_dist)]
    df['depth'] = np.int32(df.distance_pc / width)
    averages = df.groupby(['depth', 'pixel'], as_index=False).mean()

    row_titles = ['U', 'V', 'W', 'M13', 'M23', 'M12', 'M11', 'M22', 'M33']

    for px in [df.px, 1]:
        df.px = px

        print(df.px)

        columns = []

        for i in range(bins_count):
            ring = averages[(averages.depth == i)]
            _, res, _ = compute_coeffs(ring, ogorod_fg_vr)

            cur_coefs = {}
            cur_errors = {}
            for i in range(len(res.params)):
                val = res.params[i]
                err = res.bse[i]
                key = res.params.keys()[i]
                cur_coefs[key] = val
                cur_errors[key] = err

            column = {}
            for title in row_titles:
                column[title] = (cur_coefs[title], cur_errors[title])
            columns.append(column)

            # print(cur_coefs)
            # print(cur_coefs)

        # print(min(df['depth']))
        # print(max(df['depth']))

        print(end='\t')
        for i in range(bins_count):
            print(f'{i * width}-{(i + 1) * width}', end='\t')
        print()
        for title in row_titles:
            print(title, end='\t')
            for i in range(bins_count):
                print(f'{columns[i][title][0]:.1f}'.replace('.', ','), end='\t')
            print()

        print(end='\t')
        for i in range(bins_count):
            print(f'{i * width}-{(i + 1) * width}', end='\t\t')
        print()
        for title in row_titles:
            print(title, end='\t')
            for i in range(bins_count):
                c = columns[i][title]
                print(f'{c[0]:.3f}\t{c[1]:.3f}'.replace('.', ','), end='\t')
            print()


def main():
    df = read_gaia_with_rv_full()
    print(df.l)
    print(df.b)
    df['pixel'] = compute_pixel(64, df.l, df.b)
    df['distance_pc'] = distance(df)
    solve(df, 400, 20)

    solve(df, 6000, 200)

    # df[]


if __name__ == '__main__':
    main()