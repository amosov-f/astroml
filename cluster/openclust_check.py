import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

def main():
    found = [
        (179.54, -22.84, 0, 50),
        (-152.65, 31.94, 100, 150),
        (21.19, -13.31, 300, 350),
        (133.54, -2.49, 350, 400),
        (166.42, -23.96, 100, 150),
        (-7.32, 18.12, 100, 150),
        (-1.36, -15.5, 100, 150),
    ]

    df = pd.read_csv('openclust.txt', sep='\t')
    print(df)
    for l, b, lr, rr in found:

        center = SkyCoord(frame="galactic", l=l * u.deg, b=b * u.deg)
        filtered = df[((lr <= df.distance) & (df.distance <= rr))]

        c = SkyCoord(frame="galactic", l=filtered.lii.values, b=filtered.bii.values, unit=u.deg)
        sep = center.separation(c)
        deg = np.min(sep.degree)
        i = np.argmin(sep.degree)
        print(f'{(l, b, lr, rr)} closest {deg} is №{i} {filtered.iloc[i]}')

if __name__ == '__main__':
    main()