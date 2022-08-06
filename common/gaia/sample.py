import pandas as pd

def main():
    df = pd.read_csv('gdr3_short_notnull.csv')
    sample = df.sample(300000)
    sample.to_csv('gdr3_short_sample.csv')


if __name__ == '__main__':
    main()