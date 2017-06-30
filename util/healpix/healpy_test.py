import numpy as np
import healpy as hp
import pixelfunc

NSIDE = 32


def index_to_decl_ra(index):
    theta, phi = hp.pixelfunc.pix2ang(NSIDE, index)
    return -np.degrees(theta - np.pi / 2.), np.degrees(np.pi * 2. - phi)


def decl_ra_to_index(decl, RA):
    return hp.pixelfunc.ang2pix(NSIDE, np.radians(-decl + 90.), np.radians(360. - RA))


def get_index(dec, ra):
    return pixelfunc.ang2pix(32, np.radians(-dec + 90.), np.radians(360. - ra))


if __name__ == '__main__':
    print get_index(0.005616347754574377, 44.99615220936693)
