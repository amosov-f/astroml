import random

import healpy
import numpy


def ang2pix(nside, theta, phi):
    ns_max = 1048576
    TWOPI = 2. * numpy.pi
    HALFPI = numpy.pi / 2.0
    twothird = 2. / 3.

    SID = "ang2pix_ring:"

    if nside < 1 or nside > ns_max:
        raise Exception(SID + " Nside should be power of 2 >0 and < " + str(ns_max))

    if theta < 0.0 or theta > numpy.pi:
        raise Exception(SID + " Theta out of range [0,pi]")

    z = numpy.cos(theta)
    za = abs(z)

    if phi >= TWOPI:
        phi -= TWOPI

    if phi < 0.:
        phi += TWOPI  # phi in [0, 2pi]
    tt = phi / HALFPI  # in [0, 4]

    # tt = BitManipulation.MODULO(phi, TWOPI) / HALFPI; # in [0, 4]
    nl2 = 2 * nside
    nl4 = 4 * nside
    ncap = nl2 * (nside - 1)  # number of pixels in the north polar cap
    npix = 12 * nside * nside
    if za < twothird:  # equatorial region
        jp = long(nside * (0.5 + tt - 0.75 * z))  # index of ascending
        # edge line
        jm = long(nside * (0.5 + tt + 0.75 * z))  # index of descending
        # edge line

        ir = nside + 1 + jp - jm  # in[1, 2n+1]
        kshift = 0
        if ir % 2 == 0:
            kshift = 1  # 1 if ir even, 0 otherwise
        ip = (jp + jm - nside + kshift + 1) / 2 + 1  # in[1, 4n]
        if ip > nl4:
            ip -= nl4
        ipix1 = ncap + nl4 * (ir - 1) + ip

    else:  # North and South polar caps
        tp = tt - long(tt)
        tmp = numpy.sqrt(3.0 * (1.0 - za))
        jp = long(nside * tp * tmp)  # increasing edge line index
        jm = long(nside * (1.0 - tp) * tmp)  # decreasing edge index

        ir = jp + jm + 1  # ring number counted from closest pole
        ip = long(tt * ir) + 1  # in[1, 4 * ir]
        if ip > 4 * ir:
            ip -= 4 * ir

        ipix1 = 2 * ir * (ir - 1) + ip
        if z <= 0.0:
            ipix1 = npix - 2 * ir * (ir + 1) + ip

    pix = ipix1 - 1  # in [0, npix - 1]

    return pix


if __name__ == '__main__':
    for i in range(0, 10000000):
        pow = random.randint(1, 19)
        nside = 2 << (pow - 1)
        theta = random.uniform(0, numpy.pi)
        phi = random.uniform(0, 2 * numpy.pi)

        healpy_pix = healpy.ang2pix(nside, theta, phi)
        my_pix = ang2pix(nside, theta, phi)

        if healpy_pix != my_pix:
            raise Exception("%d %f %f -> %d vs %d" % (nside, theta, phi, healpy_pix, my_pix))
