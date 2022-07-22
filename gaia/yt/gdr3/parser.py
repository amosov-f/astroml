import pandas as pd

def main():

    i = 0
    chunksize = 10 ** 5
    for chunk in pd.read_csv('sample.csv', chunksize=chunksize):
        print(i)
        df = chunk[['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'radial_velocity']]
        df.to_csv('filtered.csv', mode='a', header=i == 0)
        i += 1


if __name__ == '__main__':
    main()