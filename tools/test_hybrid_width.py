"""Self-tests for tools/hybrid_width.py.

    cd tools && python3 -m unittest -v test_hybrid_width

The envelope against a direct loop over every node's square (cut and periodic), its closed form on a ramp, the
calibration ridges (the range each gives) and the input errors.
"""
import contextlib
import io
import os
import tempfile
import unittest

import numpy as np

import hybrid_width as hw

RIDGE = ["--domain", "-4096", "0", "4096", "40", "--dy", "40", "--periodic", "y", "--max-cells", "8"]


def direct(h, px, py, perx, pery):
    """h_s by a loop over each node's square (Source/HybridTerrain.cpp, written out)."""
    nx, ny = h.shape

    def square(i, j):
        ii = [(i + s) % nx if perx else i + s for s in range(-px, px + 1) if perx or 0 <= i + s < nx]
        jj = [(j + s) % ny if pery else j + s for s in range(-py, py + 1) if pery or 0 <= j + s < ny]
        return ii, jj

    m = np.empty_like(h)
    for i in range(nx):
        for j in range(ny):
            ii, jj = square(i, j)
            m[i, j] = h[np.ix_(ii, jj)].min()
    e = np.empty_like(h)
    for i in range(nx):
        for j in range(ny):
            ii, jj = square(i, j)
            e[i, j] = m[np.ix_(ii, jj)].mean()
    return np.minimum(e, h)


def run(args):
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            code = hw.main(args)
        except SystemExit as e:
            code = e.code
    return code, out.getvalue()


class Envelope(unittest.TestCase):
    def test_matches_a_direct_loop_cut_and_periodic(self):
        rng = np.random.default_rng(7)
        h = rng.uniform(0, 50, (13, 9))
        for px, py, perx, pery in ((2, 1, False, False), (3, 2, True, False), (1, 3, False, True), (2, 2, True, True)):
            np.testing.assert_allclose(hw.lower_envelope(h, px, py, perx, pery), direct(h, px, py, perx, pery),
                                       rtol=0, atol=1e-12)

    def test_zero_width_is_the_terrain_and_the_envelope_never_rises_above_it(self):
        h = np.random.default_rng(3).uniform(0, 10, (8, 8))
        np.testing.assert_array_equal(hw.lower_envelope(h, 0, 0), h)
        self.assertTrue((hw.lower_envelope(h, 2, 2) <= h).all())

    def test_ramp_interior_is_lowered_by_p_cells(self):
        # h = s x: the minimum over [i - p, i + p] is h(i - p), its mean over the same window is h(i - p) too.
        s, dx, p = 0.3, 10.0, 3
        x = dx * np.arange(40)
        h = (s * x)[:, None]
        e = hw.lower_envelope(h, p, 0)
        np.testing.assert_allclose(e[2 * p:-2 * p, 0], s * (x[2 * p:-2 * p] - p * dx), rtol=0, atol=1e-12)


class Calibration(unittest.TestCase):
    def test_gentle_ridge_range_is_terrain_following_to_128_m(self):
        code, out = run(["--agnesi", "200", "500", "ridge_x", "--dx", "64"] + RIDGE)
        self.assertEqual(code, 0, out)
        self.assertIn("Suggested: terrain.hybrid_half_width = 0   (range 0 to 128 m)", out)

    def test_steep_ridge_at_dx_32_suggests_two_cells(self):
        code, out = run(["--agnesi", "400", "250", "ridge_x", "--dx", "32"] + RIDGE)
        self.assertEqual(code, 0, out)
        self.assertIn("Suggested: terrain.hybrid_half_width = 64   (range 64 to 64 m)", out)

    def test_steep_ridge_at_dx_64_has_no_width_and_names_the_spacing(self):
        code, out = run(["--agnesi", "400", "250", "ridge_x", "--dx", "64"] + RIDGE)
        self.assertEqual(code, 2, out)
        self.assertIn("No width meets both rules", out)
        self.assertIn("about dx = 32 m", out)

    def test_one_cell_in_y_without_periodicity_steps_by_dx_as_gust_does(self):
        # mesh.n_cell y = 1 and not periodic: two nodes in y, but GUST's default counts only dx.
        code, out = run(["--agnesi", "400", "250", "ridge_x", "--dx", "32", "--domain", "-4096", "0", "4096", "40",
                         "--dy", "40", "--max-cells", "8"])
        self.assertEqual(code, 0, out)
        self.assertIn("Suggested: terrain.hybrid_half_width = 64   (range 64 to 64 m)", out)

    def test_deck_gives_the_same_answer_as_the_arguments(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "inputs")
            with open(p, "w") as f:
                f.write("mesh.n_cell = 256 1 50\nmesh.prob_lo = -4096 0 0\nmesh.prob_hi = 4096 40 1000\n"
                        "mesh.is_periodic = 0 1 0\nterrain.type = agnesi   # a ridge\nterrain.hill_shape = ridge_x\n"
                        "terrain.hill_height = 400.\nterrain.hill_half_width = 250.\n")
            code, out = run(["--deck", p, "--max-cells", "8"])
        self.assertEqual(code, 0, out)
        self.assertIn("Suggested: terrain.hybrid_half_width = 64 ", out)


class Errors(unittest.TestCase):
    def test_each_input_error_exits_3_naming_the_input(self):
        base = ["--agnesi", "200", "500", "ridge_x", "--dx", "64"]
        cases = [
            (base + RIDGE[:-4] + ["--periodic", "z"], "--periodic"),
            (base + ["--domain", "0", "0", "100", "40", "--dy", "40"], "whole number of cells"),
            (["--dx", "64", "--domain", "0", "0", "128", "64"], "--terrain or --agnesi"),
            (["--agnesi", "200", "500", "cone", "--dx", "64", "--domain", "0", "0", "128", "64"], "--agnesi shape"),
            (base + RIDGE + ["--gap-fraction", "0"], "--gap-fraction"),
            (base + RIDGE + ["--min-cells", "-1"], "--min-cells"),
            (["--agnesi", "200", "500", "ridge_x", "--dx", "-1", "--domain", "0", "0", "128", "64"], "--dx"),
            (["--deck", "/nonexistent/inputs"], "--deck"),
        ]
        for args, word in cases:
            code, out = run(args)
            self.assertEqual(code, 3, (args, out))
            self.assertIn(word, out)


if __name__ == "__main__":
    unittest.main()
