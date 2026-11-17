"""Self-tests for tools/speed_difference.py.

    cd tools && python3 -m unittest -v test_speed_difference

Hand-made fields with known answers: the counts above 0.25, 0.5 and 1 m/s, the column interpolation to
other heights (exact for a field linear in z), the coarsening mean, the near-ground mask and the column maximum.
"""
import unittest

import numpy as np

import speed_difference as sd


class SpeedDifference(unittest.TestCase):
    def test_counts_are_strictly_above_each_threshold(self):
        d = np.array([0.1, 0.25, 0.3, -0.6, 0.5, 1.0, 1.2, -2.0])
        self.assertEqual(sd.counts(d), [6, 4, 2])

    def test_at_height_is_exact_for_a_linear_profile(self):
        zr = np.cumsum(np.full((2, 3, 5), 2.0), axis=2) + np.arange(3)[None, :, None]
        fr = 3.0 + 0.5 * zr
        zn = zr[:, :, :4] + 0.7
        np.testing.assert_allclose(sd.at_height(zn, zr, fr), 3.0 + 0.5 * zn, rtol=0, atol=1e-12)

    def test_coarsen_is_the_two_by_two_mean(self):
        a = np.arange(4 * 6 * 2, dtype=float).reshape(4, 6, 2)
        c = sd.coarsen(a, 2)
        self.assertEqual(c.shape, (2, 3, 2))
        self.assertAlmostEqual(c[1, 2, 1], np.mean(a[2:4, 4:6, 1]))

    def test_near_ground_and_column_max(self):
        z = np.broadcast_to(np.array([1.0, 50.0, 150.0, 250.0]), (2, 2, 4)).copy()
        m = sd.near_ground(z, 200.0)
        self.assertEqual(m[0, 0].tolist(), [True, True, True, False])
        d = np.zeros((2, 2, 4))
        d[1, 0, 3] = 5.0    # above the layer: not counted
        d[1, 0, 1] = -0.7
        self.assertEqual(sd.column_max(d, m)[1, 0], 0.7)

    def test_speed_is_horizontal(self):
        self.assertEqual(sd.speed({"u": np.array([3.0]), "v": np.array([4.0])})[0], 5.0)


if __name__ == "__main__":
    unittest.main()
