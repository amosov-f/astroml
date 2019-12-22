import unittest

from GC2019.gc2019 import Spheric, to_galactocentric


class Test(unittest.TestCase):
    def setUp(self) -> None:
        print('========')

    def test1(self):
        star = Spheric(ra=0, dec=0, pmra=0, pmdec=0, rv=0, dist=0)
        print(star)
        res = to_galactocentric(star)
        print(res)

    def test_galactic_center(self):
        ra = self.to_degrees(17, 45.6)
        dec = -28.94
        star = Spheric(ra=ra, dec=dec, pmra=0, pmdec=0, rv=0, dist=8.5)
        print(star)
        res = to_galactocentric(star)
        print(res)

    def test_north_pool(self):
        ra = self.to_degrees(12, 51.4)
        dec = 27.13
        star = Spheric(ra=ra, dec=dec, pmra=0, pmdec=0, rv=100000, dist=0)
        print(star)
        res = to_galactocentric(star)
        print(res)

    def to_degrees(self, h, m):
        return h * (360 / 24) + m * (360 / (24 * 60))