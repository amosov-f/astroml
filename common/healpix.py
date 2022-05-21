import healpy as hp
import numpy as np
from pandas import DataFrame


def pixel_centers(nside: int):
    ls = []
    bs = []
    for ipix in range(12 * nside ** 2):
        theta, phi = hp.pix2ang(nside, ipix)
        l, b = (phi, - theta + np.pi / 2)
        ls.append(l)
        bs.append(b)
    return DataFrame({
        'l': ls,
        'b': bs
    })

# def show_distribution(df, nside):
#