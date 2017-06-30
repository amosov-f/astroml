import unittest

from math import sin, cos, sqrt, pi

import numpy as np

import spherical


class TestSpherical(unittest.TestCase):
    def indexes_test(self):
        self.assertEqual(spherical.indexes(45), (6, 5, 0))

    def pl_test(self):
        for d in np.arange(-np.pi / 2, np.pi / 2, 0.01):
            self.assertAlmostEqual(spherical.pl(3, 2, d), 15 * sin(d) * cos(d) ** 2)

    def fv_test(self):
        for d in np.arange(-np.pi / 2, np.pi / 2, 0.01):
            for a in np.arange(0, 2 * np.pi, 0.01):
                self.assertAlmostEqual(spherical.fk(3, 2, 0, a, d), 15 * sin(d) * cos(d) ** 2 * sin(2 * a))
                self.assertAlmostEqual(spherical.fk(3, 2, 1, a, d), 15 * sin(d) * cos(d) ** 2 * cos(2 * a))

    def fr_test(self):
        self.assertEqual(spherical.fr(3, 2), sqrt(7 / (4 * pi)) * sqrt(2.0 / (5 * 4 * 3 * 2)))
