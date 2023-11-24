from pathlib import Path

import pandas as pd


def read_gdr3():
    result = []

    i = 0
    chunksize = 10 ** 5
    for df in pd.read_csv(Path(__file__).parent.joinpath('gdr3.csv'), chunksize=chunksize):
        print(i)
        result.append(df)
        i += 1

    return pd.concat(result)


def main():
    df = read_gdr3()
    df.to_csv('phot.csv', index=False)


if __name__ == '__main__':
    main()
