import unittest

from math import sin, cos, sqrt, pi

import numpy as np

from common.spherical import indexes, pl, fk, fr


class TestSpherical(unittest.TestCase):
    def test_indexes(self):
        self.assertEqual(indexes(45), (6, 5, 0))

    def test_pl(self):
        for d in np.arange(-np.pi / 2, np.pi / 2, 0.01):
            self.assertAlmostEqual(pl(3, 2, d), 15 * sin(d) * cos(d) ** 2)

    def test_fv(self):
        for d in np.arange(-np.pi / 2, np.pi / 2, 0.01):
            for a in np.arange(0, 2 * np.pi, 0.01):
                self.assertAlmostEqual(fk(3, 2, 0, a, d), 15 * sin(d) * cos(d) ** 2 * sin(2 * a))
                self.assertAlmostEqual(fk(3, 2, 1, a, d), 15 * sin(d) * cos(d) ** 2 * cos(2 * a))

    def test_fr(self):
        self.assertEqual(fr(3, 2), sqrt(7 / (4 * pi)) * sqrt(2.0 / (5 * 4 * 3 * 2)))
