"""Tests of tools/gust_library/ensemble.py, tools/gust_library/grids.py and the ensemble and output options of
tools/make_gust_library.py (CTest gust_ensemble_selftest).

The libraries here are written by hand (index.json, grid.npz, one .npz per sector) with fields chosen so that
the answers are closed forms.

    - two sectors with known fields, a two-point direction distribution: the mean, the spread, the
      percentiles, the exceedance and the reversal equal the closed form (test_two_point_distribution)
    - the ASCII grid header (cell size, corner, no data) and a round trip of the values
      (test_grid_header_and_round_trip)
    - the library tab of the window: tools/test_gust_window.py (LibraryTabTests)

    cd tools && python3 -m unittest -v test_gust_ensemble
"""
import contextlib
import io
import json
import math
import os
import sys
import tempfile
import unittest

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import make_gust_library  # noqa: E402
from gust_case import geo  # noqa: E402
from gust_library import Library, LibraryError, ensemble, grids, index  # noqa: E402

NX, NY, DX = 6, 4, 25.0
X = -62.5 + DX * np.arange(NX)          # column centres, case frame (m)
Y = -37.5 + DX * np.arange(NY)
HEIGHTS = [10.0, 40.0]
REF = (12.5, 12.5)                      # a column centre
REVERSED = (4, 1)                       # the column where the wind from 90 degrees blows back towards the east


def f0(x, y):
    """Speed per unit reference speed of the sector from 0 degrees (u = 0, v = -f0: towards the south)."""
    return 1.0 + 0.004 * x


def f90(x, y):
    """Speed per unit reference speed of the sector from 90 degrees (u = -f90, v = 0: towards the west)."""
    return 0.6 + 0.01 * y


def write_library(case_dir, case_resolved=None, y=Y, unconverged=()):
    """A two-sector library (0 and 90 degrees) in case_dir/wind_library; stored per unit reference speed. The
    fields at the reference point are not 1 (f0 = 1.05, f90 = 0.725 at REF): the tests compare with the fields,
    not with the speed given. y: the column centres in y; unconverged: directions marked as not converged."""
    folder = os.path.join(case_dir, index.FOLDER)
    os.makedirs(folder)
    ny = len(y)
    gx, gy = np.meshgrid(X, y, indexing="ij")
    np.savez(os.path.join(folder, "grid.npz"), x=X, y=y, ground=np.zeros((NX, ny)))
    sectors = []
    for d in (0.0, 90.0):
        u = np.zeros((len(HEIGHTS), NX, ny))
        v = np.zeros_like(u)
        for n, h in enumerate(HEIGHTS):
            scale = 1.0 + 0.01 * (h - HEIGHTS[0])
            if d == 0.0:
                v[n] = -scale * f0(gx, gy)
            else:
                u[n] = -scale * f90(gx, gy)
                u[n][REVERSED] = +0.3 * scale            # a recirculation: blows towards the east
        k = np.full_like(u, 0.01)
        name = "sector_%03d.npz" % d
        np.savez(os.path.join(folder, name), u=u, v=v, w=np.zeros_like(u), k=k)
        sectors.append(dict(direction=d, folder="dir_%03d" % d, file=name, converged=d not in unconverged,
                            status="steady" if d not in unconverged else "stalled",
                            reference_direction=d, reference_speed=1.0))
    with open(os.path.join(folder, index.INDEX), "w") as fh:
        json.dump(dict(format_version=index.FORMAT_VERSION, case="twosector", heights=HEIGHTS, reference_point=REF,
                       reference_height=HEIGHTS[0], grid="grid.npz", sectors=sectors), fh)
    if case_resolved is not None:
        with open(os.path.join(case_dir, "case_resolved.yaml"), "w") as fh:
            yaml.safe_dump(case_resolved, fh)


def domain(**given):
    d = dict(center=None, center_lat_lon=None, utm_zone=None, half_width=[100.0], west=None, east=None, south=None,
             north=None, top=1000.0, origin="center")
    d.update(given)
    return {"case": {"name": "twosector"}, "domain": d}


def quiet_main(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = make_gust_library.main(argv)
    return rc, out.getvalue(), err.getvalue()


class EnsembleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.case = os.path.join(self.tmp.name, "case")
        write_library(self.case)
        self.lib = Library(self.case)
        gx, gy = np.meshgrid(X, Y, indexing="ij")
        self.s0 = 6.0 * f0(gx, gy)                    # member speeds at 10 m
        self.s90 = 10.0 * f90(gx, gy)
        self.s90[REVERSED] = 10.0 * 0.3

    def tearDown(self):
        self.tmp.cleanup()

    def test_two_point_distribution(self):
        """from 0 at 6 m/s (weight 1) and from 90 at 10 m/s (weight 3): weights 0.25 and 0.75."""
        members = [(0.0, 6.0, 1.0), (90.0, 10.0, 3.0)]
        st = ensemble.ensemble(self.lib, members, 10.0, percentiles=(10.0, 25.0, 30.0, 90.0), thresholds=(6.5,))
        w0, w90 = 0.25, 0.75
        s0, s90 = self.s0, self.s90
        np.testing.assert_allclose(st["mean_speed"], w0 * s0 + w90 * s90, rtol=1e-13)
        np.testing.assert_allclose(st["std_speed"], math.sqrt(w0 * w90) * np.abs(s0 - s90), rtol=1e-12, atol=1e-13)
        u90 = -s90.copy()
        u90[REVERSED] = +s90[REVERSED]
        np.testing.assert_allclose(st["mean_u"], w90 * u90, rtol=1e-13)
        np.testing.assert_allclose(st["mean_v"], -w0 * s0, rtol=1e-13)
        frm = np.degrees(np.arctan2(w90 * -u90, w0 * s0))          # from: atan2(-u, -v) clockwise from north
        np.testing.assert_allclose(st["mean_direction"], np.mod(frm, 360.0), rtol=1e-12)
        lo, hi = np.minimum(s0, s90), np.maximum(s0, s90)
        w_lo = np.where(s0 <= s90, w0, w90)
        self.assertTrue((s0 < s90).any() and (s0 > s90).any())     # the order of the two changes across the grid
        for p in (10.0, 25.0, 30.0, 90.0):
            np.testing.assert_array_equal(st[ensemble.percentile_name(p)], np.where(p / 100.0 <= w_lo, lo, hi),
                                          err_msg="p%g" % p)
        np.testing.assert_allclose(st["pexceed_6p5"], w0 * (s0 > 6.5) + w90 * (s90 > 6.5), rtol=0, atol=1e-15)
        expect = np.zeros((NX, NY))
        expect[REVERSED] = w90
        np.testing.assert_allclose(st["preversal"], expect, rtol=0, atol=1e-15)

    def test_blended_member_and_speed_scaling(self):
        """A member between sectors is Library.wind's blend; members of one direction share it, scaled."""
        st = ensemble.ensemble(self.lib, [(45.0, 4.0, 1.0), (45.0, 8.0, 1.0)], 40.0, percentiles=(50.0, 100.0))
        w = self.lib.wind(45.0, 1.0, 40.0)
        s1 = np.hypot(w["u"], w["v"])
        np.testing.assert_allclose(st["mean_speed"], 6.0 * s1, rtol=1e-13)
        np.testing.assert_allclose(st["std_speed"], 2.0 * s1, rtol=1e-12)
        np.testing.assert_allclose(st["p50"], 4.0 * s1, rtol=1e-13)    # cumulative 0.5 reached by the slower one
        np.testing.assert_allclose(st["p100"], 8.0 * s1, rtol=1e-13)

    def test_points_are_statistics_of_the_members_at_the_point(self):
        """At a point between columns the statistics are those of the interpolated members."""
        members = [(0.0, 6.0, 0.5), (90.0, 10.0, 0.5)]
        x, y = 0.0, 0.0                                            # midway between four columns
        pt = ensemble.at_points(self.lib, members, 10.0, [(x, y), REF], percentiles=(50.0,), thresholds=())
        a, b = 6.0 * f0(x, y), 10.0 * f90(x, y)                    # both fields are linear: exact
        self.assertAlmostEqual(pt["mean_speed"][0], 0.5 * (a + b), places=12)
        self.assertAlmostEqual(pt["p50"][0], min(a, b), places=12)
        self.assertAlmostEqual(pt["mean_speed"][1], 0.5 * (6.0 * f0(*REF) + 10.0 * f90(*REF)), places=12)

    def test_member_refusals(self):
        for members, word in (([], "no members"), ([(0.0, 5.0, -1.0)], "weights must be 0 or more"),
                              ([(0.0, 5.0, 0.0)], "sum to 0"), ([(0.0, 0.0, 1.0)], "above 0"),
                              ([(0.0, float("nan"), 1.0)], "not a finite"), ([(0.0, 5.0)], "not (direction")):
            with self.assertRaises(LibraryError) as cm:
                ensemble.normalise(members)
            self.assertIn(word, str(cm.exception))
        with self.assertRaisesRegex(LibraryError, "--direction-weights: 1 weights for 2 directions"):
            ensemble.product([0.0, 90.0], [1.0], [5.0], None)
        with self.assertRaisesRegex(LibraryError, "--speed-weights: the weights sum to 0"):
            ensemble.product([0.0], None, [5.0], [0.0])
        for p, t, word in (([0.0], [], "(0, 100]"), ([101.0], [], "(0, 100]"), ([], [-1.0], "0 m/s or more"),
                           ([50.0, 50.0], [], "given twice")):
            with self.assertRaises(LibraryError) as cm:
                ensemble.ensemble(self.lib, [(0.0, 5.0, 1.0)], 10.0, p, t)
            self.assertIn(word, str(cm.exception))

    def test_product_and_file(self):
        m = ensemble.product([0.0, 90.0], [1.0, 3.0], [5.0, 10.0], [2.0, 2.0])
        self.assertEqual(m, [(0.0, 5.0, 0.125), (0.0, 10.0, 0.125), (90.0, 5.0, 0.375), (90.0, 10.0, 0.375)])
        path = os.path.join(self.tmp.name, "members.txt")
        with open(path, "w") as fh:
            fh.write("# direction speed weight\n0 5 1\n\n90, 10   # weight 1 left out\n450 7 2\n")
        self.assertEqual(ensemble.from_file(path), [(0.0, 5.0, 0.25), (90.0, 10.0, 0.25), (90.0, 7.0, 0.5)])
        with open(path, "w") as fh:
            fh.write("0 5 1 2\n")
        with self.assertRaisesRegex(LibraryError, "members.txt:1: .* a member is 'direction speed"):
            ensemble.from_file(path)

    def test_percentile_of_many_equal_members(self):
        """n members of equal weight with speeds 1..n: the P-th percentile is the speed P n / 100 (the cumulative
        weight of n terms rounds by about n eps: 1.2e-12 short of 0.4 at n = 3e5 and 1.0e-12 short of 0.95 at
        n = 1e6, beyond a fixed tolerance of 1e-12)."""
        ps = (10.0, 40.0, 50.0, 90.0, 95.0)
        for n in (1000, 100000, 300000, 1000000):
            speeds = np.arange(1.0, n + 1.0)
            members = ensemble.normalise([(270.0, s, 1.0) for s in speeds])
            st = ensemble.statistics(speeds[:, None], np.zeros((n, 1)), members, percentiles=ps)
            for p in ps:
                self.assertEqual(st[ensemble.percentile_name(p)][0], p * n / 100.0, "n %d p%g" % (n, p))

    def test_zero_weight_member_is_never_a_percentile(self):
        """A member of weight 0 is slower than the others: no percentile, however small, returns its speed."""
        members = [(270.0, 1.0, 0.0), (270.0, 2.0, 1.0)]
        st = ensemble.statistics(np.array([[1.0], [2.0]]), np.zeros((2, 1)), members, percentiles=(1e-10, 50.0))
        self.assertEqual(st["p1e-10"][0], 2.0)
        self.assertEqual(st["p50"][0], 2.0)

    def test_direction_just_below_zero_is_zero(self):
        """-1e-20 % 360 rounds to 360.0; the member's direction is 0."""
        self.assertEqual(ensemble.normalise([(-1e-20, 5.0, 1.0)]), [(0.0, 5.0, 1.0)])
        self.assertEqual(ensemble.normalise([(-90.0, 5.0, 1.0)]), [(270.0, 5.0, 1.0)])

    def test_directions_that_cancel_have_no_mean_direction(self):
        """Two opposite members: the vector mean is round-off (4e-16 of the mean speed), its direction NaN."""
        members = [(0.0, 5.0, 0.5), (180.0, 5.0, 0.5)]
        st = ensemble.statistics(np.zeros((2, 1)), np.array([[-5.0], [5.0 + 8e-15]]), members)
        self.assertGreater(abs(st["mean_v"][0]), 0.0)
        self.assertTrue(np.isnan(st["mean_direction"][0]))
        st = ensemble.statistics(np.zeros((2, 1)), np.array([[-5.0], [4.0]]), members)
        self.assertAlmostEqual(st["mean_direction"][0], 0.0, places=12)

    def test_member_on_an_unconverged_sector_is_named(self):
        case = os.path.join(self.tmp.name, "stalled")
        write_library(case, unconverged=(90.0,))
        lib = Library(case)
        with self.assertRaises(LibraryError) as cm:
            ensemble.ensemble(lib, [(0.0, 5.0, 1.0), (45.0, 5.0, 1.0)], 10.0, source="--members m.txt")
        text = str(cm.exception)
        self.assertTrue(text.startswith("--members m.txt: member 2 (from 45 degrees): the sector dir_090"), text)
        self.assertNotIn("--direction", text)


class GridTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_grid_header_and_round_trip(self):
        """the header, the north row first, NaN as no data, the values back."""
        gx, gy = np.meshgrid(X, Y, indexing="ij")
        field = 3.0 + 0.01 * gx - 0.002 * gy ** 2
        field[2, NY - 1] = np.nan
        g = grids.Georeference((500000.0, 4100000.0), 13, False)
        path = os.path.join(self.tmp.name, "f.asc")
        files = grids.write_grid(path, X, Y, field, g)
        self.assertEqual(files, [path, os.path.join(self.tmp.name, "f.prj")])
        with open(path) as fh:
            head = [fh.readline().split() for _ in range(7)]
        self.assertEqual([h[0] for h in head[:6]], ["ncols", "nrows", "xllcorner", "yllcorner", "cellsize",
                                                   "NODATA_value"])
        self.assertEqual((int(head[0][1]), int(head[1][1])), (NX, NY))
        self.assertEqual(float(head[2][1]), 500000.0 + X[0] - 0.5 * DX)       # 499925
        self.assertEqual(float(head[3][1]), 4100000.0 + Y[0] - 0.5 * DX)      # 4099950
        self.assertEqual(float(head[4][1]), DX)
        self.assertEqual(float(head[5][1]), -9999.0)
        north = np.where(np.isnan(field[:, -1]), -9999.0, field[:, -1])
        np.testing.assert_allclose([float(a) for a in head[6]], north, rtol=1e-8)          # north row, west first
        self.assertEqual(head[6][2], "-9999")
        header, back, bx, by = grids.read_grid(path)
        np.testing.assert_allclose(back, field, rtol=1e-8)
        self.assertTrue(np.isnan(back[2, NY - 1]))
        np.testing.assert_allclose(bx, X + 500000.0)
        np.testing.assert_allclose(by, Y + 4100000.0)
        with open(files[1]) as fh:
            prj = fh.read()
        self.assertIn('PROJCS["WGS_1984_UTM_Zone_13N"', prj)
        self.assertIn('PARAMETER["Central_Meridian",-105.0]', prj)
        self.assertIn('PARAMETER["False_Northing",0.0]', prj)
        south = grids.Georeference((0.0, 0.0), 19, True).prj()
        self.assertIn('Zone_19S', south)
        self.assertIn('PARAMETER["False_Northing",10000000.0]', south)
        self.assertIn('PARAMETER["Central_Meridian",-69.0]', south)

    def test_no_prj_without_zone_and_copied_prj(self):
        path = os.path.join(self.tmp.name, "f.asc")
        self.assertEqual(grids.write_grid(path, X, Y, np.zeros((NX, NY)), grids.Georeference()), [path])
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "f.prj")))
        grids.write_grid(path, X, Y, np.zeros((NX, NY)), grids.Georeference((0, 0), 13), prj_text="PROJCS[\"a\"]")
        with open(os.path.join(self.tmp.name, "f.prj")) as fh:
            self.assertEqual(fh.read(), "PROJCS[\"a\"]\n")                 # a line end is added

    def test_copied_prj_that_is_one_of_the_outputs(self):
        """--prj naming the .prj of the first grid: every grid gets the text it had before the export."""
        case = os.path.join(self.tmp.name, "case")
        write_library(case, domain(center_lat_lon=[43.40, -113.03]))
        lib = Library(case)
        out = os.path.join(self.tmp.name, "out")
        os.makedirs(out)
        mine = os.path.join(out, "w_vel.prj")
        with open(mine, "w") as fh:
            fh.write("PROJCS[\"my projection\"]\n")
        w = lib.wind(0.0, 5.0, 10.0)
        fields = dict(speed=np.hypot(w["u"], w["v"]), u=w["u"], v=w["v"], k=w["k"])
        files, _, note = grids.export_grids(lib, out, "w", fields, ("speed", "u", "v"), extra=("k",), prj=mine)
        self.assertEqual(note, "")
        prjs = [f for f in files if f.endswith(".prj")]
        self.assertEqual(len(prjs), 3)
        for f in prjs:
            with open(f) as fh:
                self.assertEqual(fh.read(), "PROJCS[\"my projection\"]\n", f)

    def test_columns_that_are_not_square_are_resampled(self):
        """Columns of 25 m by 50 m: cells of 25 m from the domain's south-west corner, bilinear inside, the edge
        value beyond the outermost centres."""
        y = 2.0 * Y                                                  # -75, -25, 25, 75
        xs, ys, note = grids.grid_of(X, y)
        np.testing.assert_array_equal(xs, X)
        np.testing.assert_allclose(ys, -87.5 + 25.0 * np.arange(8), rtol=0, atol=1e-12)
        self.assertEqual(note, "columns of 25 m by 50 m resampled onto 6 x 8 square cells of 25 m (bilinear)")
        gx, gy = np.meshgrid(X, y, indexing="ij")
        field = 2.0 + 0.03 * gx - 0.05 * gy
        cells = grids.resample(X, y, field, xs, ys)
        cx, cy = np.meshgrid(xs, np.clip(ys, y[0], y[-1]), indexing="ij")
        np.testing.assert_allclose(cells, 2.0 + 0.03 * cx - 0.05 * cy, rtol=1e-14, atol=1e-14)
        self.assertEqual(grids.grid_of(X, Y)[2], "")
        square = field[:, :NY]
        self.assertIs(grids.resample(X, Y, square, X, Y), square)

    def test_cells_refused(self):
        path = os.path.join(self.tmp.name, "f.asc")
        with self.assertRaisesRegex(LibraryError, "square cells"):
            grids.write_grid(path, X, 2.0 * Y, np.zeros((NX, NY)), grids.Georeference())
        x = X.copy()
        x[-1] += 1.0
        with self.assertRaisesRegex(LibraryError, "spacing in x is not uniform"):
            grids.write_grid(path, x, Y, np.zeros((NX, NY)), grids.Georeference())
        with self.assertRaisesRegex(LibraryError, "a field of shape"):
            grids.write_grid(path, X, Y, np.zeros((NY, NX)), grids.Georeference())

    def test_direction_from_convention(self):
        """Wind blowing towards the east (u > 0) comes from 270; towards the south (v < 0) from 0; calm is NaN."""
        d = grids.direction_from(np.array([2.0, 0.0, -1.0, 0.0, 0.0]), np.array([0.0, -3.0, 0.0, 1.0, 0.0]))
        np.testing.assert_allclose(d[:4], [270.0, 0.0, 90.0, 180.0], atol=1e-12)
        self.assertTrue(np.isnan(d[4]))

    def test_georeference_from_case_resolved(self):
        """Centre in latitude and longitude: UTM origin and zone; south-west origin; projected centre; none."""
        lat, lon = 43.40, -113.03
        e, n, zone, south = geo.latlon_to_utm(lat, lon)
        for origin, shift in (("center", (0.0, 0.0)), ("south-west corner", (-100.0, -100.0))):
            case = os.path.join(self.tmp.name, origin.replace(" ", "_"))
            write_library(case, domain(center_lat_lon=[lat, lon], origin=origin))
            g = grids.georeference(Library(case))
            self.assertEqual((g.zone, g.south), (zone, south))
            self.assertAlmostEqual(g.origin[0], float(e) + shift[0], places=6)
            self.assertAlmostEqual(g.origin[1], float(n) + shift[1], places=6)
        case = os.path.join(self.tmp.name, "projected")
        write_library(case, domain(center=[345000.0, 4800000.0]))
        g = grids.georeference(Library(case))
        self.assertEqual((g.origin, g.zone, g.prj()), ((345000.0, 4800000.0), None, None))
        case = os.path.join(self.tmp.name, "bare")
        write_library(case)
        g = grids.georeference(Library(os.path.join(case, index.FOLDER)))
        self.assertEqual((g.origin, g.zone), ((0.0, 0.0), None))
        self.assertIn("not found", g.source)


class CommandLineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.case = os.path.join(self.tmp.name, "case")
        write_library(self.case, domain(center_lat_lon=[43.40, -113.03]))
        self.out = os.path.join(self.tmp.name, "out")
        self.points = os.path.join(self.tmp.name, "points.txt")
        with open(self.points, "w") as fh:
            fh.write("12.5 12.5 mast A\n0 0\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_ensemble_writes_grids_and_points(self):
        rc, out, err = quiet_main(["ensemble", self.case, "--height", "10", "--directions", "0", "90",
                                   "--direction-weights", "1", "3", "--speeds", "6", "10", "--speed-weights", "1", "0",
                                   "--thresholds", "6.5", "--grids", self.out, "--points", self.points,
                                   "--out", os.path.join(self.tmp.name, "st")])
        self.assertEqual(rc, 0, err)
        self.assertIn("ensemble: 4 members, 2 directions", out)
        names = sorted(os.listdir(self.out))
        stem = "twosector_ensemble_10m_agl"
        for s in ("vel", "ang", "std_speed", "p10", "p50", "p90", "pexceed_6p5", "preversal"):
            self.assertIn("%s_%s.asc" % (stem, s), names)
            self.assertIn("%s_%s.prj" % (stem, s), names)
        st = dict(np.load(os.path.join(self.tmp.name, "st.npz")))
        _, vel, _, _ = grids.read_grid(os.path.join(self.out, stem + "_vel.asc"))
        np.testing.assert_allclose(vel, st["mean_speed"], rtol=1e-8)
        _, ang, _, _ = grids.read_grid(os.path.join(self.out, stem + "_ang.asc"))
        np.testing.assert_allclose(ang, st["mean_direction"], rtol=1e-8)
        with open(os.path.join(self.out, stem + "_points.csv")) as fh:
            rows = [line.rstrip("\n").split(",") for line in fh]
        self.assertEqual(rows[0][:6], ["name", "x", "y", "easting", "northing", "height"])
        self.assertEqual(rows[1][0], "mast A")
        self.assertEqual(rows[2][0], "point2")
        col = rows[0].index("mean_speed")
        self.assertAlmostEqual(float(rows[1][col]), 0.25 * 6.0 * f0(*REF) + 0.75 * 6.0 * f90(*REF), places=7)   # 10 m/s has weight 0

    def test_query_grids_direction_from_and_units(self):
        """A wind from 90 (towards the west): ang 90 away from the recirculation, 270 in it; vel in mph."""
        rc, out, err = quiet_main(["query", self.case, "--direction", "90", "--speed", "10", "--height", "10",
                                   "--grids", self.out, "--speed-units", "mph", "--points", self.points])
        self.assertEqual(rc, 0, err)
        stem = "twosector_90deg_10mps_10m_agl"
        _, ang, _, _ = grids.read_grid(os.path.join(self.out, stem + "_ang.asc"))
        expect = np.full((NX, NY), 90.0)
        expect[REVERSED] = 270.0
        np.testing.assert_allclose(ang, expect, atol=1e-9)
        _, vel, _, _ = grids.read_grid(os.path.join(self.out, stem + "_vel.asc"))
        gx, gy = np.meshgrid(X, Y, indexing="ij")
        speed = 10.0 * f90(gx, gy)
        speed[REVERSED] = 3.0
        np.testing.assert_allclose(vel, speed * 3600.0 / 1609.344, rtol=1e-8)
        with open(os.path.join(self.out, stem + "_points.csv")) as fh:
            rows = [line.rstrip("\n").split(",") for line in fh]
        self.assertAlmostEqual(float(rows[1][rows[0].index("speed")]), 10.0 * f90(*REF), places=7)  # CSV in m/s
        self.assertAlmostEqual(float(rows[1][rows[0].index("direction")]), 90.0, places=7)
        e, n, _, _ = geo.latlon_to_utm(43.40, -113.03)
        self.assertAlmostEqual(float(rows[1][rows[0].index("easting")]), float(e) + REF[0], places=2)

    def test_speed_units_apply_to_every_speed_grid(self):
        """--speed-units kph: the mean speed, the spread and the percentiles in km/h; probabilities unchanged."""
        st_path = os.path.join(self.tmp.name, "st")
        rc, out, err = quiet_main(["ensemble", self.case, "--height", "10", "--directions", "0", "90", "--speeds",
                                   "6", "10", "--percentiles", "50", "--thresholds", "6", "--grids", self.out,
                                   "--speed-units", "kph", "--out", st_path])
        self.assertEqual(rc, 0, err)
        st = dict(np.load(st_path + ".npz"))
        stem = os.path.join(self.out, "twosector_ensemble_10m_agl_")
        for name, factor in (("vel", 3.6), ("std_speed", 3.6), ("p50", 3.6), ("pexceed_6", 1.0),
                             ("preversal", 1.0)):
            _, g, _, _ = grids.read_grid(stem + name + ".asc")
            key = {"vel": "mean_speed"}.get(name, name)
            np.testing.assert_allclose(g, st[key] * factor, rtol=1e-8, atol=1e-12, err_msg=name)

    def test_columns_that_are_not_square_on_the_command_line(self):
        case = os.path.join(self.tmp.name, "tall")
        write_library(case, domain(center_lat_lon=[43.40, -113.03]), y=2.0 * Y)
        rc, out, err = quiet_main(["query", case, "--direction", "0", "--speed", "10", "--height", "10",
                                   "--grids", self.out])
        self.assertEqual(rc, 0, err)
        self.assertIn("resampled onto 6 x 8 square cells of 25 m", out)
        header, vel, _, _ = grids.read_grid(os.path.join(self.out, "twosector_0deg_10mps_10m_agl_vel.asc"))
        self.assertEqual(vel.shape, (6, 8))
        self.assertEqual(float(header["cellsize"]), 25.0)
        np.testing.assert_allclose(vel, 10.0 * f0(X, 0.0)[:, None] * np.ones((1, 8)), rtol=1e-8)   # f0 is linear in x
        _, ang, _, _ = grids.read_grid(os.path.join(self.out, "twosector_0deg_10mps_10m_agl_ang.asc"))
        np.testing.assert_allclose(ang, 0.0, atol=1e-9)

    def test_output_options_that_would_be_ignored(self):
        plain_file = os.path.join(self.tmp.name, "plain_file")
        with open(plain_file, "w") as fh:
            fh.write("x\n")
        for argv, word in ((["--prj", plain_file], "--prj: only used with --grids"),
                           (["--speed-units", "kph"], "--speed-units: only used with --grids"),
                           (["--prefix", "p"], "--prefix: only used with --grids or --points"),
                           (["--grids", plain_file], "is a file, not a folder")):
            for command in (["query", self.case, "--direction", "0", "--speed", "5"],
                            ["ensemble", self.case, "--directions", "0", "--speeds", "5"]):
                rc, out, err = quiet_main(command + ["--height", "10"] + argv)
                self.assertEqual(rc, 2, command + argv)
                self.assertIn(word, err)
                self.assertNotIn("wrote", out)

    def test_refusals(self):
        members = os.path.join(self.tmp.name, "m.txt")
        with open(members, "w") as fh:
            fh.write("0 5\n")
        for argv, word in ((["--members", members, "--directions", "0"], "not both"),
                           ([], "give --members FILE"),
                           (["--directions", "0", "--speeds", "5", "--percentiles", "0"], "(0, 100]"),
                           (["--directions", "0", "--speeds", "5", "--direction-weights", "-1"], "0 or more"),
                           (["--directions", "0", "--speeds", "5", "--points", members + "x"], "cannot read"),
                           (["--directions", "0", "--speeds", "5", "--grids", self.out, "--prj", "nope.prj"],
                            "--prj: nope.prj not found")):
            rc, out, err = quiet_main(["ensemble", self.case, "--height", "10"] + argv)
            self.assertEqual(rc, 2, argv)
            self.assertIn(word, err)


if __name__ == "__main__":
    unittest.main()
