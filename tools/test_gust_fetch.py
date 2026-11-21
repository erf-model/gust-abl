"""Tests of the upstream roughness per direction (tools/gust_case/fetch.py, wind.upstream_roughness = auto),
CTest gust_fetch_selftest.

The numbers pinned here are closed forms: the integral of exp(-x/L) over a step map, the geometric mean of a
half-plane map, the flux weights of two inflow faces, and the log law's u* of each direction's z0.

    cd tools && python3 -m unittest -v test_gust_fetch
"""
import math
import os
import re
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gust_case import CaseError, build, validate, wind  # noqa: E402
from gust_case import fetch  # noqa: E402
from gust_case.deck import read_deck  # noqa: E402
from gust_case.terrain import Frame  # noqa: E402


def frame(half_x, half_y):
    return Frame((0.0, 0.0), half_x, half_x, half_y, half_y)


def grid_points(x, y, z0_of):
    gx, gy = np.meshgrid(x, y)
    return gx.ravel(), gy.ravel(), z0_of(gx, gy).ravel()


def step_closed_form(z0_near, z0_far, s, length):
    """Fetch-weighted ln z0 when z0_near holds for x < s upwind of the face and z0_far beyond, out to REACH L."""
    reach = fetch.REACH * length
    near = 1.0 - math.exp(-s / length)
    far = math.exp(-s / length) - math.exp(-reach / length)
    return (near * math.log(z0_near) + far * math.log(z0_far)) / (1.0 - math.exp(-reach / length))


def base_case(map_path, **wind_keys):
    w = {"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": "auto"}
    w.update(wind_keys)
    return {"case": {"name": "t"}, "domain": {"half_width": [2000.0], "top": 1000.0},
            "mesh": {"cell_size": 250.0, "layers": 16, "first_layer": 8.0},
            "wind": w, "surface": {"roughness": 0.1, "regional_roughness_map": {"file": map_path}}}


class BinnedMapTests(unittest.TestCase):
    def test_bin_size_is_the_raster_spacing(self):
        x = 312.5 + 40.0 * np.arange(37)
        y = -900.0 + 40.0 * np.arange(23)
        binned = fetch.BinnedMap(*grid_points(x, y, lambda gx, gy: 0.1 + 0.0 * gx), "k")
        self.assertAlmostEqual(binned.h, 40.0, places=9)
        self.assertEqual((binned.nx, binned.ny), (37, 23))

    def test_a_point_takes_its_own_bin_and_outside_is_not_covered(self):
        x = 10.0 * np.arange(5)
        binned = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: np.exp(0.1 * gx + gy)), "k")
        self.assertAlmostEqual(float(binned.sample(np.array([20.0]), np.array([30.0]))[0]), 32.0, places=12)
        self.assertAlmostEqual(float(binned.sample(np.array([24.9]), np.array([25.1]))[0]), 32.0, places=12)
        self.assertTrue(np.isnan(binned.sample(np.array([-5.1]), np.array([0.0]))[0]))
        self.assertTrue(np.isnan(binned.sample(np.array([0.0]), np.array([45.1]))[0]))

    def test_an_empty_bin_is_filled_from_its_neighbours(self):
        x = 10.0 * np.arange(5)
        gx, gy, z0 = grid_points(x, x, lambda gx, gy: np.where(gx < 20.0, 1.0, math.e ** 2))
        keep = ~((gx == 20.0) & (gy == 20.0))
        binned = fetch.BinnedMap(gx[keep], gy[keep], z0[keep], "k")
        # the eight neighbours of (20, 20): three at x = 10 (ln 0), three at x = 30 and two at x = 20 (ln 2)
        self.assertAlmostEqual(float(binned.sample(np.array([20.0]), np.array([20.0]))[0]), 10.0 / 8.0, places=12)

    def test_a_line_of_points_is_refused(self):
        with self.assertRaisesRegex(CaseError, "k.file: the points lie on a line"):
            fetch.BinnedMap(np.arange(5.0), np.zeros(5), np.ones(5), "k")


class UpstreamTests(unittest.TestCase):
    def test_half_plane_map_gives_each_side_its_own_z0(self):
        """UT-13.4.1: z0 0.03 m west of the domain, 1 m east: the wind from 270 gets 0.03 m, from 90 gets 1 m."""
        x = np.arange(-8000.0, 8000.1, 100.0)
        binned = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: np.where(gx < 0.0, 0.03, 1.0)), "k")
        f = frame(2000.0, 2000.0)
        ground = fetch.uniform_ground(0.1)
        west = fetch.upstream_roughness(binned, f, 270.0, 1000.0, ground)
        east = fetch.upstream_roughness(binned, f, 90.0, 1000.0, ground)
        self.assertAlmostEqual(west.z0, 0.03, places=12)
        self.assertAlmostEqual(east.z0, 1.0, places=12)
        self.assertEqual([face.face for face in west.faces], ["xlo"])
        self.assertEqual([face.face for face in east.faces], ["xhi"])
        self.assertAlmostEqual(west.faces[0].mismatch(), math.log(0.3), places=12)
        self.assertAlmostEqual(west.cover, 1.0, places=12)

    def test_uniform_map_returns_its_z0_for_every_sector(self):
        """UT-13.4.2."""
        x = np.arange(-9000.0, 9000.1, 250.0)
        binned = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: 0.2 + 0.0 * gx), "k")
        for direction in np.arange(0.0, 360.0, 22.5):
            r = fetch.upstream_roughness(binned, frame(2000.0, 1000.0), float(direction), 1000.0, fetch.uniform_ground(0.2))
            self.assertAlmostEqual(r.z0, 0.2, places=12, msg="direction %g" % direction)
            for face in r.faces:
                self.assertAlmostEqual(face.mismatch(), 0.0, places=12)

    def test_step_map_weight_against_the_closed_form(self):
        """UT-13.4.3: z0_near within s of the face, z0_far beyond: the weight exp(-x/L) integrated in closed form.
        Map points every h upwind of the face; a step half-way between two points lies on a segment edge, so the
        discrete sum is the integral."""
        x_lo, h = -200.0, 5.0
        near, far = 0.5, 0.02
        for length in (60.0, 100.0, 250.0):
            for steps in (3, int(0.2 * fetch.REACH * length / h), int(0.7 * fetch.REACH * length / h)):
                s = (steps + 0.5) * h
                x = x_lo - h * np.arange(int(fetch.REACH * length / h) + 3)[::-1]
                y = np.arange(-200.0, 200.1, h)
                binned = fetch.BinnedMap(*grid_points(x, y, lambda gx, gy: np.where(x_lo - gx < s, near, far)), "k")
                segments = int(min(fetch.MAX_SEGMENTS, max(fetch.SEGMENTS_PER_FETCH * fetch.REACH,
                                                           math.ceil(2.0 * fetch.REACH * length / binned.h))))
                self.assertEqual(segments, int(round(2.0 * fetch.REACH * length / h)), "segments of h/2")
                r = fetch.upstream_roughness(binned, frame(200.0, 200.0), 270.0, length, fetch.uniform_ground(0.1))
                self.assertAlmostEqual(math.log(r.z0), step_closed_form(near, far, s, length), places=12,
                                       msg="L = %g, s = %g" % (length, s))

    def test_two_inflow_faces_are_weighted_by_the_flow_through_them(self):
        """A 4000 x 2000 m domain, wind from 225: xlo (2000 m) and ylo (4000 m) take the flow at the same normal
        speed, so the direction's ln z0 is (2000 ln z0_xlo + 4000 ln z0_ylo) / 6000."""
        x = np.arange(-10000.0, 10000.1, 100.0)
        # the line y = x + 1000 runs upwind from the corner (-2000, -1000): 0.05 m above it (upwind of xlo), 0.8 m below
        binned = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: np.where(gy - gx > 1000.0, 0.05, 0.8)), "k")
        r = fetch.upstream_roughness(binned, frame(2000.0, 1000.0), 225.0, 300.0, fetch.uniform_ground(0.1))
        by_face = {face.face: face for face in r.faces}
        self.assertEqual(sorted(by_face), ["xlo", "ylo"])
        expected = (2000.0 * math.log(by_face["xlo"].z0_upwind) + 4000.0 * math.log(by_face["ylo"].z0_upwind)) / 6000.0
        self.assertAlmostEqual(math.log(r.z0), expected, places=12)
        # transects next to the corner sample bins that straddle the line
        self.assertGreater(by_face["ylo"].z0_upwind, 0.7)
        self.assertLess(by_face["xlo"].z0_upwind, 0.06)

    def test_turned_faces_are_weighted_by_the_turned_flux(self):
        """turning -15 degrees from 270: xlo (4010 m) and yhi (4000 m) are the inflow faces, weighted by the turned
        direction's normal components cos 15 and sin 15; the transects run west along the wind asked for. The map is
        a south of y = 2005 m and b north of it, so xlo's transects see a and yhi's see b."""
        a, b, h = 0.02, 0.5, 20.0
        f = Frame((0.0, 0.0), 2000.0, 2000.0, 2000.0, 2010.0)
        x = np.arange(-3100.0, 2100.1, h)
        y = np.arange(-2100.0, 2100.1, h)
        binned = fetch.BinnedMap(*grid_points(x, y, lambda gx, gy: np.where(gy > 2005.0, b, a)), "k")
        r = fetch.upstream_roughness(binned, f, 270.0, 200.0, fetch.uniform_ground(0.1), turning=-15.0)
        by_face = {face.face: face for face in r.faces}
        self.assertEqual(sorted(by_face), ["xlo", "yhi"])
        self.assertAlmostEqual(by_face["xlo"].z0_upwind, a, places=12)
        self.assertAlmostEqual(by_face["yhi"].z0_upwind, b, places=12)
        c, s = math.cos(math.radians(15.0)), math.sin(math.radians(15.0))
        expected = (4010.0 * c * math.log(a) + 4000.0 * s * math.log(b)) / (4010.0 * c + 4000.0 * s)
        self.assertAlmostEqual(math.log(r.z0), expected, places=12)
        self.assertEqual(r.turning, -15.0)

    def test_the_ground_at_the_face_comes_from_the_inner_map(self):
        x = np.arange(-6000.0, 6000.1, 100.0)
        regional = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: 0.3 + 0.0 * gx), "k")
        inner_x = np.arange(-2000.0, 2000.1, 50.0)
        inner = fetch.BinnedMap(*grid_points(inner_x, inner_x, lambda gx, gy: np.where(gy < 0.0, 0.01, 0.04)), "i")
        r = fetch.upstream_roughness(regional, frame(2000.0, 2000.0), 270.0, 500.0, inner.sample)
        # 40 points along xlo, half below y = 0: the geometric mean of 0.01 and 0.04
        self.assertAlmostEqual(r.faces[0].z0_ground, math.sqrt(0.01 * 0.04), places=12)
        self.assertAlmostEqual(r.faces[0].mismatch(), math.log(0.3 / 0.02), places=12)

    def test_an_inner_map_starting_inside_the_face_gives_its_edge_bin(self):
        x = np.arange(-6000.0, 6000.1, 100.0)
        regional = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: 0.3 + 0.0 * gx), "k")
        inner_x = np.arange(-1930.0, 1930.1, 50.0)           # first point 70 m inside the faces
        # 0.05 m west of x = 0, 0.2 m east: a wrapped index would read the east column
        inner = fetch.BinnedMap(*grid_points(inner_x, inner_x, lambda gx, gy: np.where(gx < 0.0, 0.05, 0.2)), "i")
        self.assertTrue(np.isnan(inner.sample(np.array([-2000.0]), np.array([0.0]))[0]))
        r = fetch.upstream_roughness(regional, frame(2000.0, 2000.0), 270.0, 500.0, inner.sample_nearest_edge)
        self.assertAlmostEqual(r.faces[0].z0_ground, 0.05, places=12)

    def test_an_unknown_ground_leaves_out_the_mismatch(self):
        x = np.arange(-6000.0, 6000.1, 100.0)
        regional = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: 0.3 + 0.0 * gx), "k")
        r = fetch.upstream_roughness(regional, frame(2000.0, 2000.0), 270.0, 500.0,
                                     lambda px, py: np.full(np.shape(px), np.nan))
        self.assertIsNone(r.faces[0].z0_ground)
        self.assertIsNone(r.faces[0].mismatch())
        self.assertIn("xlo: upwind 0.3 m, ground not known along the face", "\n".join(fetch.report_lines([r])))

    def test_a_map_of_two_spacings_is_refused_as_uneven(self):
        """A 10 m raster over the domain inside a 400 m raster: the mean spacing leaves the coarse part's bins
        empty; the message says so rather than asking for a larger map."""
        fine = np.arange(-2000.0, 2000.1, 10.0)
        coarse = np.arange(-9200.0, 9200.1, 400.0)
        cx, cy, cz = grid_points(coarse, coarse, lambda gx, gy: 0.1 + 0.0 * gx)
        outside = (np.abs(cx) > 2000.0) | (np.abs(cy) > 2000.0)
        fx, fy, fz = grid_points(fine, fine, lambda gx, gy: 0.1 + 0.0 * gx)
        binned = fetch.BinnedMap(np.concatenate([fx, cx[outside]]), np.concatenate([fy, cy[outside]]),
                                 np.concatenate([fz, cz[outside]]), "k")
        with self.assertRaisesRegex(CaseError, "within the map's extent, where its points are further apart"):
            fetch.upstream_roughness(binned, frame(2000.0, 2000.0), 270.0, 1000.0, fetch.uniform_ground(0.1),
                                     key="k")

    def test_a_short_map_is_noted_or_refused(self):
        length = 1000.0
        f = frame(2000.0, 2000.0)
        for reach, refused in ((3500.0, False), (500.0, True)):
            x = np.arange(-2000.0 - reach, 2000.1, 50.0)
            binned = fetch.BinnedMap(*grid_points(x, x, lambda gx, gy: 0.1 + 0.0 * gx), "k")
            if refused:
                with self.assertRaisesRegex(CaseError, "k.file covers 41 % of the fetch weight upwind of the wind "
                                                       "from 270 degrees.*wind.fetch_length = 1000 m"):
                    fetch.upstream_roughness(binned, f, 270.0, length, fetch.uniform_ground(0.1), key="k")
            else:
                r = fetch.upstream_roughness(binned, f, 270.0, length, fetch.uniform_ground(0.1), key="k")
                # covered out to 3500 m (+ half a bin): (1 - exp(-3.525)) / (1 - exp(-5))
                self.assertAlmostEqual(r.cover, (1.0 - math.exp(-3525.0 / length)) / (1.0 - math.exp(-5.0)), places=2)
                self.assertAlmostEqual(r.z0, 0.1, places=12)


class CaseTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)
        x = np.arange(-12000.0, 12000.1, 200.0)
        gx, gy, z0 = grid_points(x, x, lambda gx, gy: np.where(gx < 0.0, 0.03, 1.0))
        self.map = os.path.join(self.dir, "regional.xyz")
        np.savetxt(self.map, np.column_stack([gx, gy, z0]), fmt="%.1f %.1f %.6g")

    def build(self, raw, dry_run=False):
        case, given = validate(raw)
        return build(case, given, self.dir, check=False, dry_run=dry_run)

    def refused(self, raw):
        with self.assertRaises(CaseError) as caught:
            self.build(raw, dry_run=True)
        return caught.exception.problems

    def test_sweep_decks_take_each_directions_z0_and_ustar(self):
        raw = base_case(self.map, direction=None, directions=[90.0, 270.0])
        result = self.build(raw)
        for label, z0 in (("dir_270", 0.03), ("dir_090", 1.0)):
            deck = read_deck(os.path.join(result.case_dir, label, "inputs"))
            self.assertAlmostEqual(float(deck["bc.inflow_z0"][0]), z0, places=9)
            self.assertAlmostEqual(float(deck["bc.inflow_ustar"][0]), wind.friction_velocity(8.0, 10.0, z0), places=9)
            self.assertEqual(deck["bc.surface_z0"], ["0.1"])       # the ground inside the domain is unchanged
        self.assertIn("wind from 270: z0 = 0.03 m, map covers 100.0 % of the weight", result.summary)
        self.assertIn("xhi: upwind 1 m, ground 0.1 m, ln mismatch +2.303", result.summary)
        self.assertIn("L = 1000 m", result.summary)               # the default: 100 x wind.height

    def test_geostrophic_column_deck_takes_the_directions_z0(self):
        for drive in ("geostrophic_inflow", "geostrophic_column"):
            with self.subTest(drive=drive):
                self._geostrophic(drive)

    def _geostrophic(self, drive):
        raw = base_case(self.map, drive=drive, direction=None, directions=[90.0, 270.0])
        raw["atmosphere"] = {"latitude": 45.0}
        raw["case"]["name"] = drive
        result = self.build(raw)
        # the faces turn with the drag law: a y face is an inflow face too and its transects cross the map's step
        # at x = 0, so z0 lies between the two sides; the decks take the summary's value and its faces
        listed = re.findall(r"wind from (\d+): z0 = ([0-9.e+-]+) m", result.summary)
        self.assertEqual([d for d, _ in listed], ["90", "270"])
        for label, (_, text) in zip(("dir_090", "dir_270"), listed):
            z0 = float(text)
            column = read_deck(os.path.join(result.case_dir, label, "column", "inputs"))
            terrain = read_deck(os.path.join(result.case_dir, label, "inputs"))
            inflow = sorted(face for face in ("xlo", "xhi", "ylo", "yhi") if terrain["bc." + face] == ["inflow"])
            self.assertEqual(len(inflow), 2, "a turned face is an inflow face")
            block = result.summary.split("wind from %s:" % label[4:].lstrip("0"))[1].split("wind from")[0]
            self.assertEqual(sorted(re.findall(r"^\s+([xy](?:lo|hi)): upwind", block, re.M)), inflow)
            self.assertAlmostEqual(float(column["bc.surface_z0"][0]) / z0, 1.0, places=5)   # the summary's 6 digits
            self.assertEqual(terrain["bc.inflow_z0"], column["bc.surface_z0"])
            self.assertTrue(0.03 < z0 < 1.0)
            z0 = float(column["bc.surface_z0"][0])
            ustar = wind.friction_velocity(8.0, 10.0, z0)
            start = wind.geostrophic_speed(ustar, z0, 45.0)
            self.assertAlmostEqual(math.hypot(float(column["forcing.ug"][0]), float(column["forcing.vg"][0])), start,
                                   places=6)

    def test_fetch_length_is_given_to_the_transects(self):
        result = self.build(base_case(self.map, direction=0.0, fetch_length=400.0), dry_run=True)
        self.assertIn("L = 400 m", result.summary)

    def test_land_cover_regional_map(self):
        cover = os.path.join(self.dir, "cover.asc")
        codes = np.where(np.arange(121)[None, :] < 60, 130, 70) * np.ones((121, 1), dtype=int)
        with open(cover, "w", encoding="utf-8") as out:
            out.write("ncols 121\nnrows 121\nxllcorner -12100\nyllcorner -12100\ncellsize 200\n")
            np.savetxt(out, codes, fmt="%d")
        raw = base_case(cover, direction=None, directions=[90.0, 270.0])
        result = self.build(raw)
        west = read_deck(os.path.join(result.case_dir, "dir_270", "inputs"))
        east = read_deck(os.path.join(result.case_dir, "dir_090", "inputs"))
        self.assertAlmostEqual(float(west["bc.inflow_z0"][0]), 0.07, places=9)     # shrubland (130)
        self.assertAlmostEqual(float(east["bc.inflow_z0"][0]), 1.0, places=9)      # needle-leaved forest (70)

    def test_refusals_name_their_keys(self):
        table = os.path.join(self.dir, "inflow.dat")
        with open(table, "w", encoding="utf-8") as out:
            out.write("0 0 0 300\n1000 10 0 300\n")
        cases = [
            (dict(surface={"roughness": 0.1}), "wind.upstream_roughness = auto needs surface.regional_roughness_map.file"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0}),
             "surface.regional_roughness_map.file has no effect without wind.upstream_roughness = auto"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": 0.1, "fetch_length": 500.0}),
             "wind.fetch_length has no effect without wind.upstream_roughness = auto"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": "auto", "drive": "couette"}),
             "wind.upstream_roughness = auto: not with wind.drive = couette"),
            (dict(wind={"direction": 270.0, "upstream_roughness": "auto", "drive": "geostrophic_inflow"},
                  atmosphere={"latitude": 45.0, "geostrophic_wind": [10.0, 0.0]}),
             "wind.fetch_length: required with wind.upstream_roughness = auto"),
            (dict(wind={"direction": 270.0, "upstream_roughness": "auto", "profile": "table", "table": table,
                        "fetch_length": 500.0}),
             "wind.upstream_roughness = auto: not with wind.profile = table"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": "auto",
                        "profile": "power_law", "power_law_exponent": 0.14}),
             "wind.upstream_roughness = auto: not with wind.profile = power_law"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": "auto",
                        "fetch_length": 0.0}), "wind.fetch_length = 0.0: must be a number > 0"),
            (dict(wind={"speed": 8.0, "height": 10.0, "direction": 270.0, "upstream_roughness": "automatic"}),
             "wind.upstream_roughness = 'automatic': must be a number > 0, or auto"),
        ]
        for change, text in cases:
            raw = base_case(self.map)
            for section, value in change.items():
                raw[section] = value
            try:
                problems = self.refused(raw)
            except AssertionError:
                self.fail("%s: not refused (expected '%s')" % (change, text))
            self.assertTrue(any(text in p for p in problems), "%s: expected '%s' in %s" % (change, text, problems))

    def test_nodata_in_a_regional_raster_is_left_out(self):
        cover = os.path.join(self.dir, "cover_nodata.asc")
        codes = np.where(np.arange(121)[None, :] < 60, 130, 70) * np.ones((121, 1), dtype=int)
        codes[25:30, :] = -9999                                # rows 6100 .. 7100 m north of the centre hold no data
        with open(cover, "w", encoding="utf-8") as out:
            out.write("ncols 121\nnrows 121\nxllcorner -12100\nyllcorner -12100\ncellsize 200\n"
                      "NODATA_value -9999\n")
            np.savetxt(out, codes, fmt="%d")
        result = self.build(base_case(cover), dry_run=True)
        self.assertIn("points with no value left out", result.summary)
        self.assertIn("wind from 270: z0 = 0.07 m", result.summary)

    def test_a_map_that_does_not_reach_upwind_is_refused(self):
        x = np.arange(-2200.0, 12000.1, 200.0)
        gx, gy, z0 = grid_points(x, np.arange(-12000.0, 12000.1, 200.0), lambda gx, gy: 0.1 + 0.0 * gx)
        np.savetxt(self.map, np.column_stack([gx, gy, z0]), fmt="%.1f %.1f %.6g")
        with self.assertRaisesRegex(CaseError, "surface.regional_roughness_map.file covers"):
            self.build(base_case(self.map), dry_run=True)


class LibrarySummaryTests(unittest.TestCase):
    def test_roughness_per_direction_is_shown_as_a_range(self):
        from gust_case import library
        entry = {"direction": 270.0, "roughness": 0.03, "latitude": 45.0, "geostrophic_wind": {"speed": 15.0,
                 "direction": 290.0}, "wind_at_height": None, "turning_deg": None, "turning_from_stress_deg": -18.0,
                 "ustar": 0.5, "drag_law": None}
        lib = {"reference_wind": {"speed": 8.0, "height": 10.0}, "directions": [entry, dict(entry, direction=90.0,
                                                                                             roughness=1.0)]}
        self.assertIn("z0 per direction 0.03 .. 1 m", library.summary_text(lib))
        lib["directions"][1]["roughness"] = 0.03
        self.assertIn("z0 = 0.03 m", library.summary_text(lib))


if __name__ == "__main__":
    unittest.main()
