"""Self-tests for tools/nest_placement.py.

    cd tools && python3 -m unittest -v test_nest_placement

Analytic terrains with known answers: a channel through the inflow and outflow faces at 0 and 35 degrees (and a
wide one at 35 degrees that the outer cells follow), the
same site scaled down four times (the thresholds follow dx), an escarpment crossing a face, a second feature
that only some wind directions see, a channel that bends beyond the outflow face (a note, not a flag), a site
where no position is clear, the margin from the outer faces, a smooth hill, the exit status of a check, and
the input checks (each aborts with status 3 and names the input).
"""
import io
import math
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from types import SimpleNamespace

import numpy as np

import nest_placement as npl

DOMAIN = (-9600.0, -8000.0, 9600.0, 8000.0)
NEST = (-3000.0, -3000.0, 3000.0, 3000.0)


def args(dx=200.0, scale=1.0, **kw):
    a = SimpleNamespace(dx=dx, ratio=2, domain=tuple(scale * v for v in DOMAIN), direction=[270.0], angle=20.0,
                        turn=20.0, misfit=None, prominence=None, step_rise=None, step_width=None, bend_length=None,
                        downstream=None, margin=None, search_step=None, top=5, plot=None, nest=None, recommend=None)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def site(fn, a, scale=1.0):
    """Terrain fn(x, y) on the tool's fine grid (dx / 4); x, y in the unscaled frame, h scaled."""
    step = a.dx / (2 * a.ratio)
    x0, y0, x1, y1 = a.domain
    xs = np.arange(x0, x1 + 0.5 * step, step)
    ys = np.arange(y0, y1 + 0.5 * step, step)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return npl.Site(xs, ys, scale * fn(X / scale, Y / scale), a)


def channel(theta_deg, depth=300.0, sigma=400.0):
    """A Gaussian channel through the origin whose axis is theta degrees from the x axis."""
    t = math.radians(theta_deg)

    def f(X, Y):
        d = -X * math.sin(t) + Y * math.cos(t)
        return 400.0 - depth * np.exp(-(d / sigma) ** 2)
    return f


def quiet(fn, *a):
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        return fn(*a)


def faults(ev, name):
    return next(f for f in ev["faces"] if f["name"] == name)["faults"]


class NestPlacement(unittest.TestCase):
    def test_channel_along_the_wind_is_clear_and_oblique_is_flagged(self):
        a = args()
        along = quiet(site, channel(0.0), a).evaluate(NEST)
        self.assertEqual(along["flagged"], {})
        self.assertAlmostEqual(along["cover"], 6000.0 / 19200.0, delta=0.02)   # the nest's share of the channel
        oblique = quiet(site, channel(35.0), args()).evaluate(NEST)
        self.assertEqual(sorted(nm for nm, _, _ in oblique["flagged"][270.0]), ["east", "west"])
        self.assertIn("cuts a channel at 3", faults(oblique, "west")[0])

    def test_an_oblique_channel_the_outer_cells_follow_is_clear(self):
        # 35 deg, but 2 km wide: the outer mesh's misfit stays under 1 m, so the face is not flagged.
        ev = quiet(site, channel(35.0, sigma=2000.0), args()).evaluate(NEST)
        west = next(f for f in ev["faces"] if f["name"] == "west")
        self.assertTrue(any(s.startswith("channel 3") for s in west["feats"]), west)
        self.assertEqual(ev["flagged"], {})

    def test_thresholds_follow_dx(self):
        # The oblique channel scaled down 4 times (lengths and heights) with dx 50 m: the same verdict, and
        # every threshold a quarter of its value at dx 200 m.
        big = args()
        quiet(site, channel(35.0), big)
        a = args(dx=50.0, scale=0.25)
        s = quiet(site, channel(35.0), a, 0.25)
        self.assertAlmostEqual(s.relief, 75.0, delta=1.0)
        self.assertEqual((big.misfit, big.prominence, big.step_rise, big.margin), (5.0, 25.0, 25.0, 2000.0))
        self.assertEqual((a.misfit, a.prominence, a.step_rise, a.margin), (1.25, 6.25, 6.25, 500.0))
        for name in ("misfit", "prominence", "step_rise", "step_width", "bend_length", "downstream", "margin"):
            self.assertAlmostEqual(getattr(a, name), 0.25 * getattr(big, name), delta=1e-3 * getattr(big, name), msg=name)
        ev = s.evaluate(tuple(0.25 * v for v in NEST))
        self.assertEqual(sorted(nm for nm, _, _ in ev["flagged"][270.0]), ["east", "west"])

    def test_escarpment_crossing_a_face_is_flagged(self):
        t = math.radians(50.0)
        scarp = lambda X, Y: 300.0 + 200.0 * np.tanh((-(X + 3000.0) * math.sin(t) + Y * math.cos(t)) / 250.0)
        ev = quiet(site, scarp, args()).evaluate(NEST)
        west = faults(ev, "west")
        self.assertTrue(any("cuts a step" in f for f in west), west)
        self.assertIn("west", [nm for nm, _, _ in ev["flagged"][270.0]])

    def test_a_feature_seen_only_by_other_directions(self):
        # Channel along x (clear for westerlies) plus a ridge 30 deg off the y axis that crosses only the
        # south and north faces, 1.1 km from the corners: side faces for 270, inflow and outflow for 0 and 180.
        t = math.radians(60.0)
        ridge = lambda X, Y: channel(0.0)(X, Y) + 250.0 * np.exp(-((X * math.sin(t) - Y * math.cos(t)) / 400.0) ** 2)
        s = quiet(site, ridge, args(direction=list(npl.ALL_DIRECTIONS)))
        ev = s.evaluate(NEST)
        self.assertNotIn(270.0, ev["flagged"])
        self.assertIn(0.0, ev["flagged"])
        self.assertIn(180.0, ev["flagged"])
        self.assertEqual(sorted(nm for nm, _, _ in ev["flagged"][0.0]), ["north", "south"])
        self.assertTrue(npl.verdict(ev, 16).startswith("move the nest"))

    def test_a_bend_beyond_the_outflow_face_is_a_note(self):
        # Channel along y = 0 that turns north 600 m beyond the east face (x = 3600).
        def bend(X, Y):
            d = np.where(X <= 3600.0, np.abs(Y), np.hypot(X - 3600.0, np.minimum(Y, 0.0)))
            d = np.where((X > 3600.0) & (Y > 0.0), np.minimum(np.abs(X - 3600.0), np.hypot(X - 3600.0, Y)), d)
            return 400.0 - 300.0 * np.exp(-(d / 400.0) ** 2)
        ev = quiet(site, bend, args()).evaluate(NEST)
        east = next(f for f in ev["faces"] if f["name"] == "east")
        self.assertEqual(east["faults"], [])
        self.assertTrue(any("bends" in n for n in east["notes"]), east)
        self.assertEqual(ev["flagged"], {})

    def test_no_clear_position_warns_and_returns_2(self):
        t = math.radians(45.0)
        waves = lambda X, Y: 300.0 + 150.0 * np.sin(2 * math.pi * (X * math.cos(t) + Y * math.sin(t)) / 2000.0)
        a = args(recommend=(4000.0, 4000.0))
        s = quiet(site, waves, a)
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code, ok = npl.recommend(a, s)
        self.assertEqual((code, ok), (2, []))
        self.assertIn("WARNING: no position", err.getvalue())
        self.assertIn("Nothing is recommended", out.getvalue())

    def test_recommend_keeps_the_margin_and_check_notes_it(self):
        a = args(recommend=(4000.0, 4000.0))
        s = quiet(site, channel(0.0), a)
        with redirect_stdout(io.StringIO()):
            code, ok = npl.recommend(a, s)
        self.assertEqual(code, 0)
        self.assertTrue(ok)
        X0, Y0, X1, Y1 = DOMAIN
        for ev in ok:
            x0, y0, x1, y1 = ev["nest"]
            self.assertGreaterEqual(min(x0 - X0, y0 - Y0, X1 - x1, Y1 - y1), 2000.0 - 1e-6)
            self.assertEqual(ev["edge"], [])
        near = s.evaluate((-8800.0, -3000.0, -2000.0, 3000.0))
        self.assertEqual(len(near["edge"]), 1)
        self.assertIn("west face 800 m (4 outer cells)", near["edge"][0])

    def test_smooth_hill_needs_no_nest(self):
        hill = lambda X, Y: 100.0 * np.exp(-(X ** 2 + Y ** 2) / 5000.0 ** 2)
        a = args()
        out = io.StringIO()
        with redirect_stdout(out):
            s = site(hill, a)
            npl.summary(s)
        self.assertFalse(s.high.any())
        self.assertIn("a nest adds little here", out.getvalue())
        self.assertEqual(s.evaluate(NEST)["flagged"], {})

    def test_face_kind_uses_60_degrees(self):
        # West face (outward normal -x): the wind is inflow within 60 deg of the normal, a side face beyond.
        west = (-1.0, 0.0)
        self.assertEqual(npl.face_kind(west, 270.0), "inflow")
        self.assertEqual(npl.face_kind(west, 270.0 + 55.0), "inflow")
        self.assertEqual(npl.face_kind(west, 270.0 + 65.0), "side")
        self.assertEqual(npl.face_kind(west, 90.0 - 55.0), "outflow")
        self.assertEqual(npl.face_kind(west, 90.0 + 65.0), "side")

    def test_check_exit_status(self):
        # 1 when a checked nest has a flagged face, 0 when every nest is clear.
        with redirect_stdout(io.StringIO()):
            self.assertEqual(npl.check(args(nest=[NEST]), site(channel(35.0), args())), 1)
            self.assertEqual(npl.check(args(nest=[NEST]), site(channel(0.0), args())), 0)

    def test_terrain_smaller_than_the_domain_aborts(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.xyz")
            xs = np.arange(-4000.0, 4001.0, 200.0)
            X, Y = np.meshgrid(xs, xs, indexing="ij")
            np.savetxt(path, np.c_[X.ravel(), Y.ravel(), 0 * X.ravel()])
            err = io.StringIO()
            with self.assertRaises(SystemExit) as cm, redirect_stderr(err):
                npl.load_terrain(path, DOMAIN, 50.0)
        self.assertEqual(cm.exception.code, npl.INPUT_ERROR)
        self.assertIn("--domain", err.getvalue())
        self.assertIn("reaches beyond", err.getvalue())

    def test_directions_and_ranges_are_checked(self):
        dom = ["--domain", *map(str, DOMAIN)]
        base = ["--terrain", "t.xyz", *dom, "--dx", "200", "--nest", *map(str, NEST)]
        self.assertEqual(npl.parse(base + ["--direction", "all"]).direction, npl.ALL_DIRECTIONS)
        self.assertEqual(npl.parse(base + ["--direction", "250", "270"]).direction, [250.0, 270.0])
        rec = ["--terrain", "t.xyz", *dom, "--dx", "200", "--recommend", "4000", "4000"]
        self.assertEqual(npl.parse(rec + ["--search-step", "400", "--top", "3"]).top, 3)
        bad = [
            (base + ["--direction", "360"], "--direction"), (base + ["--direction", "west"], "--direction"),
            (base + ["--step-rise", "0"], "--step-rise"), (base + ["--margin", "-1"], "--margin"),
            (base + ["--misfit", "-5"], "--misfit"), (base[:-7] + ["--dx", "0", *base[-5:]], "--dx"),
            (base + ["--ratio", "1"], "--ratio"), (base + ["--angle", "90"], "--angle"),
            (base + ["--turn", "0"], "--turn"),
            (["--terrain", "t.xyz", "--domain", "9600", "-8000", "-9600", "8000", "--dx", "200", "--nest",
              *map(str, NEST)], "--domain"),
            (["--terrain", "t.xyz", "--domain", "-9600", "-8000", "9700", "8000", "--dx", "200", "--nest",
              *map(str, NEST)], "--domain"),
            (base[:-5] + ["--nest", "-3000", "-3000", "9800", "3000"], "--nest"),
            (base[:-5] + ["--nest", "-3000", "-3000", "3100", "3000"], "outer nodes"),
            (base + ["--recommend", "4000", "4000"], "--recommend"),
            (base[:-5], "--nest"),
            (base + ["--search-step", "400"], "--search-step"), (base + ["--top", "3"], "--top"),
            (rec + ["--search-step", "0"], "--search-step"), (rec + ["--top", "0"], "--top"),
            (rec[:-2] + ["4100", "4000"], "--recommend"), (rec[:-2] + ["0", "4000"], "--recommend"),
        ]
        for argv, name in bad:
            err = io.StringIO()
            with self.assertRaises(SystemExit) as cm, redirect_stderr(err):
                npl.parse(argv)
            self.assertEqual(cm.exception.code, npl.INPUT_ERROR, argv)
            self.assertIn(name, err.getvalue().splitlines()[-1], argv)

if __name__ == "__main__":
    unittest.main()
