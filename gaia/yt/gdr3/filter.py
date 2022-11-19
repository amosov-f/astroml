import pandas as pd

def main():

    result = []

    i = 0
    chunksize = 10 ** 5
    for df in pd.read_csv('gdr3.csv', chunksize=chunksize):
        print(i)
        df = df[['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'radial_velocity']]
        df = df[df.pmra.notnull()]
        df = df[df.pmdec.notnull()]
        df = df[df.parallax.notnull()]
        df = df[(df.parallax > 0)]
        result.append(df)
        i += 1

    df = pd.concat(result)

    df.to_csv('filtered_notnull.csv', index=False)


if __name__ == '__main__':
    main()