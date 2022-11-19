import pandas as pd

def main():

    result = []

    i = 0
    chunksize = 10 ** 5
    for df in pd.read_csv('gdr3.csv', chunksize=chunksize):
        print(i)
        df = df[['ra', 'dec', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'parallax', 'parallax_error', 'radial_velocity', 'radial_velocity_error']]
        result.append(df)
        i += 1

    df = pd.concat(result)

    df.to_csv('not_filtered_with_error.csv', index=False)


if __name__ == '__main__':
    main()