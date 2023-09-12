from common.gaia.with_rv import read_gaia_with_rv_full, read_raw_gaia_with_rv_no_errors, slices


def main():
    df = read_gaia_with_rv_full()

    l_dists = []
    r_dists = []
    def pizza(df, l_dist, r_dist):
        l_dists.append(str(int(l_dist)))
        r_dists.append(str(int(r_dist)))

    avgs = slices(pizza, df, imax=30000000, step=2000000, slice=2000000)

    print('min\t' + '\t'.join(l_dists))
    print('max\t' + '\t'.join(r_dists))
    print('avg\t' + '\t'.join(avgs))


if __name__ == '__main__':
    main()