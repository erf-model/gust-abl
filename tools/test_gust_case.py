"""Tests of tools/gust_case and tools/make_gust_case.py (CTest gust_case_selftest).

The numbers pinned here come from outside the package: closed forms worked
by hand, pyproj for the UTM projection, amrTerrain's land-cover table (five classes
moved to Wieringa's 1992 Davenport classes), and
the decks of the repository (Tests/Real_Terrain/WindNinja2019/inputs_bsb_srtm,
Tests/Real_Terrain/Askervein_CNBL).

    cd tools && python3 -m unittest -v test_gust_case
"""
import contextlib
import io
import json
import math
import os
import re
import struct
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import make_gust_case  # noqa: E402
import tiff_to_xyz  # noqa: E402
from gust_case import CaseError, apply_override, build, key_reference_yaml, key_table_rst, resolve, validate  # noqa: E402
from gust_case import atmosphere, geo, grid, library, roughness, schema, terrain, wind  # noqa: E402
from gust_case.deck import Deck, format_value, read_deck  # noqa: E402


def base_case():
    """Flat ground, wind from the west: the smallest complete case."""
    return {"case": {"name": "t"}, "domain": {"half_width": [2000.0], "top": 1000.0},
            "mesh": {"cell_size": 250.0, "layers": 16, "first_layer": 8.0},
            "wind": {"speed": 8.0, "height": 10.0, "direction": 270.0},
            "surface": {"roughness": 0.03}}


def changed(case, **changes):
    """A copy of ``case`` with section__key=value changes (None removes the key)."""
    import copy
    out = copy.deepcopy(case)
    for path, value in changes.items():
        node = out
        parts = path.split("__")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        if value is None:
            node.pop(parts[-1], None)
        else:
            node[parts[-1]] = value
    return out


def write_hill(path, n=81, cell=100.0, x0=500000.0, y0=4000000.0, base=1200.0, height=100.0, nodata_at=None):
    """An ESRI ASCII grid of a Gaussian hill at the centre of an (n cell) square."""
    x = x0 + cell * (np.arange(n) + 0.5)
    y = y0 + cell * (np.arange(n) + 0.5)
    gx, gy = np.meshgrid(x, y)
    half = 0.5 * n * cell
    z = base + height * np.exp(-((gx - x0 - half) ** 2 + (gy - y0 - half) ** 2) / (2.0 * 400.0 ** 2))
    if nodata_at is not None:
        z[nodata_at] = -9999.0
    with open(path, "w", encoding="utf-8") as out:
        out.write("ncols %d\nnrows %d\nxllcorner %.1f\nyllcorner %.1f\ncellsize %g\nNODATA_value -9999\n"
                  % (n, n, x0, y0, cell))
        np.savetxt(out, z[::-1], fmt="%.3f")
    return (x0 + half, y0 + half)


def write_tiff(path, values, x_west, y_north, pixel, dtype="<f4", nodata=None, point=False):
    """A little-endian, single-strip GeoTIFF with a tie point, a pixel scale and a GeoKeyDirectory.

    dtype is any of the reader's integer or float types; nodata writes GDAL_NODATA; point marks the raster
    PixelIsPoint (the tie point is the first pixel's centre)."""
    values = np.asarray(values, dtype=dtype)
    ny, nx = values.shape
    data = values.tobytes()
    sample = {"u": 1, "i": 2, "f": 3}[values.dtype.kind]
    doubles = struct.pack("<3d", pixel, pixel, 0.0) + struct.pack("<6d", 0.0, 0.0, 0.0, x_west, y_north, 0.0)
    keys = struct.pack("<12H", 1, 1, 0, 2, 1024, 0, 1, 1, 1025, 0, 1, 2 if point else 1)
    text = b"" if nodata is None else (repr(nodata).encode() + b"\0")
    entries = [(256, 4, 1, nx), (257, 4, 1, ny), (258, 3, 1, 8 * values.dtype.itemsize), (259, 3, 1, 1),
               (262, 3, 1, 1), (273, 4, 1, "data"), (277, 3, 1, 1), (278, 4, 1, ny), (279, 4, 1, len(data)),
               (339, 3, 1, sample), (33550, 12, 3, "scale"), (33922, 12, 6, "tie"), (34735, 3, 12, "keys")]
    if text:
        entries.append((42113, 2, len(text), "text"))
    ifd_size = 2 + 12 * len(entries) + 4
    at = {"scale": 8 + ifd_size}
    at["tie"] = at["scale"] + 24
    at["keys"] = at["scale"] + len(doubles)
    at["text"] = at["keys"] + len(keys)
    at["data"] = at["text"] + len(text)
    ifd = struct.pack("<H", len(entries))
    for tag, kind, count, value in entries:
        value = at.get(value, value)
        if kind == 3 and count == 1:
            ifd += struct.pack("<HHIHH", tag, kind, count, value, 0)
        else:
            ifd += struct.pack("<HHII", tag, kind, count, value)
    ifd += struct.pack("<I", 0)
    with open(path, "wb") as out:
        out.write(b"II" + struct.pack("<HI", 42, 8) + ifd + doubles + keys + text + data)


def reference_defaults():
    """{key: default tokens} of Tests/inputs_reference ('' where the key has none)."""
    defaults = {}
    with open(os.path.join(REPO, "Tests", "inputs_reference"), encoding="utf-8") as handle:
        for line in handle:
            match = re.match(r"# ([a-z_]+\.[A-Za-z0-9_]+) =(.*)$", line)
            if match:
                defaults[match.group(1)] = match.group(2).split("#", 1)[0].split()
    return defaults


def same_tokens(a, b, rel=1.0e-3):
    """Deck values agree: numbers to ``rel``, words exactly."""
    if len(a) != len(b):
        return False
    for p, q in zip(a, b):
        try:
            if not math.isclose(float(p), float(q), rel_tol=rel, abs_tol=1.0e-12):
                return False
        except ValueError:
            if p != q:
                return False
    return True


class CaseDir(unittest.TestCase):
    """A scratch folder per test, and helpers that build cases in it."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

    def build(self, raw, **options):
        case, given = validate(raw)
        options.setdefault("check", False)
        return build(case, given, self.dir, **options)

    def deck(self, raw, rel="inputs"):
        result = self.build(raw)
        return read_deck(os.path.join(result.case_dir, rel)), result

    def refused(self, raw):
        """The problems of a case that must be refused."""
        with self.assertRaises(CaseError) as caught:
            self.build(raw, dry_run=True)
        return caught.exception.problems


class WindTests(unittest.TestCase):
    def test_direction_from_becomes_direction_towards(self):
        # from the west blows along +x; from the south-west (222) is the WindNinja decks' 48 degrees
        for met, towards in ((270.0, 0.0), (222.0, 48.0), (180.0, 90.0), (90.0, 180.0), (0.0, 270.0), (360.0, 270.0)):
            self.assertAlmostEqual(wind.towards_from_meteorological(met), towards, places=12)
            self.assertAlmostEqual(wind.meteorological_from_towards(towards), met % 360.0, places=12)

    def test_faces_of_the_eight_compass_directions(self):
        side = "slip_wall"
        expected = {   # from: (xlo, xhi, ylo, yhi)
            0.0: (side, side, "outflow", "inflow"),
            45.0: ("outflow", "inflow", "outflow", "inflow"),
            90.0: ("outflow", "inflow", side, side),
            135.0: ("outflow", "inflow", "inflow", "outflow"),
            180.0: (side, side, "inflow", "outflow"),
            225.0: ("inflow", "outflow", "inflow", "outflow"),
            270.0: ("inflow", "outflow", side, side),
            315.0: ("inflow", "outflow", "outflow", "inflow"),
        }
        for met, faces in expected.items():
            got = wind.lateral_faces(wind.towards_from_meteorological(met))
            self.assertEqual((got["xlo"], got["xhi"], got["ylo"], got["yhi"]), faces, "wind from %g" % met)
        periodic = wind.lateral_faces(0.0, "periodic")
        self.assertEqual((periodic["ylo"], periodic["yhi"]), ("periodic", "periodic"))

    def test_geostrophic_faces_follow_the_turning_wind(self):
        """UT-13.1.1: with a geostrophic drive no face is a slip wall: the faces parallel to a wind along an axis
        take the side the wind aloft crosses (clockwise turn in the north, counter-clockwise in the south)."""
        turn_n = wind.face_turning(0.564336173406487, 0.03, 45.0)
        turn_s = wind.face_turning(0.564336173406487, 0.03, -45.0)
        # half of atan2(4.5, ln(0.5643 / (1.0284e-4 * 0.03)) - 1.8) = 11.78 degrees (the column turns 18.1)
        self.assertAlmostEqual(turn_n, -0.5 * math.degrees(math.atan2(4.5, math.log(
            0.564336173406487 / (1.0284451246e-4 * 0.03)) - 1.8)), places=8)
        self.assertAlmostEqual(turn_n, -11.7830, places=3)
        self.assertEqual(turn_s, -turn_n)
        north = {   # from: (xlo, xhi, ylo, yhi) at 45 N
            0.0: ("outflow", "inflow", "outflow", "inflow"),      # aloft from about 18 degrees: towards -x
            90.0: ("outflow", "inflow", "inflow", "outflow"),
            180.0: ("inflow", "outflow", "inflow", "outflow"),
            270.0: ("inflow", "outflow", "outflow", "inflow"),    # aloft from about 288: towards -y
            240.0: ("inflow", "outflow", "inflow", "outflow"),    # off the axes: as without turning
        }
        for met, faces in north.items():
            got = wind.lateral_faces(wind.towards_from_meteorological(met), "slip_wall", turn_n)
            self.assertEqual((got["xlo"], got["xhi"], got["ylo"], got["yhi"]), faces, "wind from %g, 45 N" % met)
        got = wind.lateral_faces(wind.towards_from_meteorological(270.0), "slip_wall", turn_s)
        self.assertEqual((got["ylo"], got["yhi"]), ("inflow", "outflow"), "wind from 270, 45 S")
        for met in range(0, 360, 5):
            got = wind.lateral_faces(wind.towards_from_meteorological(met), "slip_wall", turn_n)
            self.assertNotIn("slip_wall", got.values(), "wind from %g" % met)

    def test_drag_law_inverse(self):
        for z0, lat in ((0.03, 45.0), (1.0, -30.0), (0.0002, 60.0)):
            for ustar in (0.05, 0.3, 1.2):
                G = wind.geostrophic_speed(ustar, z0, lat)
                self.assertAlmostEqual(wind.drag_law_ustar(G, z0, lat) / ustar, 1.0, places=12)
        with self.assertRaises(ValueError):
            wind.drag_law_ustar(0.0, 0.03, 45.0)

    def test_log_law(self):
        # 0.41 * 8 / ln(10.03 / 0.03)
        self.assertAlmostEqual(wind.friction_velocity(8.0, 10.0, 0.03), 0.564336173406487, places=12)
        self.assertAlmostEqual(wind.log_law_speed(0.564336173406487, 10.0, 0.03), 8.0, places=12)
        # the lid of inputs_bsb_srtm: u* 0.5963 over z0 0.01 at 2806 m is (12.212, 13.562) m/s to 4 digits
        self.assertAlmostEqual(wind.log_law_speed(0.5963, 2806.0, 0.01) / math.hypot(12.212, 13.562), 1.0, places=3)

    def test_drag_law(self):
        # f = 2 (2 pi / 86400 s) sin(45 degrees), as GUST's forcing.rotation_period default; G = 0.3/0.41 * hypot(ln(0.3 / (f * 0.03)) - 1.8, 4.5)
        self.assertAlmostEqual(wind.coriolis_parameter(45.0), 1.0284451246e-4, places=13)
        self.assertLess(wind.coriolis_parameter(-45.0), 0.0)
        self.assertAlmostEqual(wind.geostrophic_speed(0.3, 0.03, 45.0), 7.814101451, places=8)
        self.assertAlmostEqual(wind.geostrophic_speed(0.3, 0.03, -45.0), 7.814101451, places=8)
        with self.assertRaises(ValueError):
            wind.geostrophic_speed(0.3, 0.03, 0.0)

    def test_sweeps_and_labels(self):
        self.assertEqual(wind.sweep_directions({"start": 0, "stop": 315, "step": 45}),
                         [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
        self.assertEqual(wind.sweep_directions([222, 48]), [222.0, 48.0])
        self.assertEqual([wind.direction_label(d) for d in (0, 48, 222, 22.5, 360)], ["000", "048", "222", "022p5", "000"])
        # to 0.001 degree: no empty fraction from a rounding to the next degree, and no wrap to '360'
        self.assertEqual([wind.direction_label(d) for d in (22.9996, 22.0004, 22.1236, 359.9996, 0.001, -0.5)],
                         ["023", "022", "022p124", "000", "000p001", "359p5"])


class GridTests(unittest.TestCase):
    def test_cell_counts_of_the_butte_deck(self):
        # inputs_bsb_srtm: 18690 x 20160 m in 272 x 296 cells, three levels with boxes of 64
        self.assertEqual(grid.horizontal_cells(18690.0, 68.4), 272)
        self.assertEqual(grid.horizontal_cells(20160.0, 68.4), 296)
        self.assertEqual(grid.hierarchy_levels(272, 296, 64), 3)
        self.assertEqual(grid.horizontal_cells(10.0, 50.0), 8)
        self.assertEqual([grid.horizontal_cells(n * 100.0, 100.0) for n in (52, 60, 68)], [56, 64, 72])   # halves go up

    def test_hierarchy_stops_where_the_counts_do(self):
        self.assertEqual(grid.hierarchy_levels(40, 40, 64), 2)     # 40/4 = 10 cells are too few
        self.assertEqual(grid.hierarchy_levels(72, 72, 36), 3)
        self.assertEqual(grid.hierarchy_levels(72, 72, 18), 2)     # the box size is not divisible by 4
        self.assertEqual(grid.hierarchy_levels(24, 64, 64), 1)     # 24/2 = 12 cells are too few
        self.assertEqual(grid.hierarchy_levels(256, 256, 64, wanted=4), 4)

    def test_stretched_layers(self):
        # inputs_bsb_srtm: 40 layers from 4 m reach 2806 m with a ratio of 1.11678
        ratio = grid.stretch_ratio_for_top(4.0, 40, 2806.0)
        self.assertAlmostEqual(ratio, 1.11678, places=5)
        self.assertAlmostEqual(grid.stretched_top(4.0, ratio, 40), 2806.0, places=8)
        self.assertEqual(grid.stretch_ratio_for_top(25.0, 40, 1000.0), 1.0)
        self.assertEqual(grid.stretched_top(25.0, 1.0, 40), 1000.0)
        # 1 + 2 + 4 + 8 = 15
        self.assertEqual(grid.stretched_top(1.0, 2.0, 4), 15.0)
        self.assertEqual(grid.layer_edges(1.0, 2.0, 4), [0.0, 1.0, 3.0, 7.0, 15.0])
        self.assertEqual(grid.layers_for_ratio(1.0, 2.0, 15.0), 4)
        self.assertEqual(grid.layers_for_ratio(1.0, 2.0, 15.5), 5)
        self.assertEqual(grid.layers_for_ratio(4.0, 1.11678, 2806.0), 40)
        with self.assertRaises(ValueError):
            grid.stretch_ratio_for_top(30.0, 40, 1000.0)   # 40 layers of 30 m pass the top


class GeoTests(unittest.TestCase):
    def test_utm_against_pyproj(self):
        # pyproj 3.6.1, EPSG:4326 -> 32612 and 32740
        east, north, zone, south = geo.latlon_to_utm(43.40, -113.02)
        self.assertEqual((zone, south), (12, False))
        self.assertAlmostEqual(east, 336423.239, places=3)
        self.assertAlmostEqual(north, 4807217.318, places=3)
        east, north, zone, south = geo.latlon_to_utm(-21.1, 55.5)
        self.assertEqual((zone, south), (40, True))
        self.assertAlmostEqual(east, 344197.575, places=3)
        self.assertAlmostEqual(north, 7666050.711, places=3)

    def test_zones(self):
        self.assertEqual([geo.utm_zone(lon) for lon in (-180.0, -113.02, -0.1, 0.0, 55.5, 179.9)], [1, 12, 30, 31, 40, 60])
        # the central meridian of zone 31 is 3 E: easting 500000 there, northing 0 on the equator
        east, north, _, _ = geo.latlon_to_utm(0.0, 3.0)
        self.assertAlmostEqual(east, 500000.0, places=6)
        self.assertAlmostEqual(north, 0.0, places=6)


class AtmosphereTests(unittest.TestCase):
    def test_inversion_sounding_of_the_askervein_cnbl_deck(self):
        # 300 K to 700 m, 9 K over 300 m, 0.03 K/m to 1600 m
        heights, theta = atmosphere.inversion_sounding(300.0, 1600.0, 700.0, 300.0, 9.0, 0.03)
        self.assertEqual(heights, [0.0, 700.0, 1000.0, 1600.0])
        for got, want in zip(theta, [300.0, 300.0, 309.0, 327.0]):
            self.assertAlmostEqual(got, want, places=10)
        self.assertAlmostEqual(atmosphere.gradient_at_top(heights, theta), 0.03, places=12)
        self.assertEqual(atmosphere.inversion_sounding(300.0, 1000.0), ([0.0], [300.0]))
        self.assertEqual(atmosphere.inversion_sounding(300.0, 1000.0, lapse_rate=0.003), ([0.0, 1000.0], [300.0, 303.0]))

    def test_refusals(self):
        with self.assertRaisesRegex(CaseError, "atmosphere.inversion"):
            atmosphere.inversion_sounding(300.0, 1000.0, 800.0, 200.0, 5.0)
        with self.assertRaisesRegex(CaseError, "first height must be 0"):
            atmosphere.check_sounding([10.0, 2000.0], [300.0, 305.0], 1000.0)
        with self.assertRaisesRegex(CaseError, "increase strictly"):
            atmosphere.check_sounding([0.0, 500.0, 500.0, 2000.0], [300.0] * 4, 1000.0)
        with self.assertRaisesRegex(CaseError, "below the top"):
            atmosphere.check_sounding([0.0, 900.0], [300.0, 305.0], 1000.0)


class TerrainTests(CaseDir):
    def test_esri_ascii_rows_run_from_north_to_south(self):
        path = os.path.join(self.dir, "g.asc")
        with open(path, "w", encoding="utf-8") as out:
            out.write("ncols 3\nnrows 2\nxllcorner 100\nyllcorner 200\ncellsize 10\nNODATA_value -9999\n"
                      "1 2 3\n4 5 -9999\n")
        raster = terrain.read_esri_ascii(path, "terrain.file")
        self.assertEqual(list(raster.x), [105.0, 115.0, 125.0])
        self.assertEqual(list(raster.y), [205.0, 215.0])
        self.assertEqual(raster.values[1].tolist(), [1.0, 2.0, 3.0])      # the northern row was first in the file
        self.assertEqual(raster.values[0][:2].tolist(), [4.0, 5.0])
        self.assertTrue(np.isnan(raster.values[0][2]))

    def test_esri_ascii_cell_centre_header(self):
        path = os.path.join(self.dir, "g.asc")
        with open(path, "w", encoding="utf-8") as out:
            out.write("ncols 2\nnrows 1\nxllcenter 100\nyllcenter 200\ncellsize 10\n7 8\n")
        raster = terrain.read_esri_ascii(path, "terrain.file")
        self.assertEqual((list(raster.x), list(raster.y)), ([100.0, 110.0], [200.0]))

    def geotiff_values(self):
        values = np.arange(12, dtype=float).reshape(3, 4) + 0.5      # row 0 is the northern one
        path = os.path.join(self.dir, "g.tif")
        write_tiff(path, values, x_west=1000.0, y_north=5000.0, pixel=20.0)
        return path, values

    def check_geotiff(self, raster, values):
        self.assertEqual(list(raster.x), [1010.0, 1030.0, 1050.0, 1070.0])
        self.assertEqual(list(raster.y), [4950.0, 4970.0, 4990.0])
        self.assertEqual(raster.values.tolist(), values[::-1].tolist())

    def test_geotiff_with_the_built_in_reader(self):
        path, values = self.geotiff_values()
        with mock.patch.dict(sys.modules, {"rasterio": None}):      # as if rasterio were not installed
            self.check_geotiff(terrain.read_geotiff(path, "terrain.file"), values)

    def test_geotiff_integer_with_no_data(self):
        """int16 heights read exactly; the GDAL_NODATA pixel becomes NaN."""
        values = np.array([[100, -9999, 3], [-7, 250, 32767]])
        path = os.path.join(self.dir, "i.tif")
        write_tiff(path, values, x_west=0.0, y_north=40.0, pixel=20.0, dtype="<i2", nodata=-9999)
        got, x0, y0, sx, sy = tiff_to_xyz.read_tiff(path)
        self.assertEqual((x0, y0, sx, sy), (0.0, 40.0, 20.0, 20.0))
        self.assertEqual(np.isnan(got).sum(), 1)
        self.assertEqual(sorted(got[~np.isnan(got)].tolist()), [-7.0, 3.0, 100.0, 250.0, 32767.0])

    def test_geotiff_pixel_is_point(self):
        """PixelIsPoint: the tie point is the north-west pixel's centre, not its corner."""
        path, values = self.geotiff_values()
        write_tiff(path, values, x_west=1010.0, y_north=4990.0, pixel=20.0, point=True)
        with mock.patch.dict(sys.modules, {"rasterio": None}):
            self.check_geotiff(terrain.read_geotiff(path, "terrain.file"), values)

    def test_geotiff_with_rasterio(self):
        try:
            import rasterio  # noqa: F401
        except ImportError:
            self.skipTest("rasterio is not installed")
        path, values = self.geotiff_values()
        self.check_geotiff(terrain.read_geotiff(path, "terrain.file"), values)

    def test_border_weight(self):
        frame = terrain.Frame((0.0, 0.0), 1000.0, 1000.0, 1000.0, 1000.0)
        x = np.array([0.0, -800.0, -800.0, -950.0, -900.0, 700.0, 1000.0])
        y = np.array([0.0, 0.0, -800.0, 0.0, 0.0, 0.0, 0.0])
        # flat for 100 m from each face, then 200 m of slope: (d - 100) / 200 per side, multiplied
        weight = terrain.border_weight(x, y, frame, 100.0, 200.0)
        for got, want in zip(weight, [1.0, 0.5, 0.25, 0.0, 0.0, 1.0, 0.0]):
            self.assertAlmostEqual(got, want, places=12)
        step = terrain.border_weight(x, y, frame, 100.0, 0.0)
        self.assertEqual(step.tolist(), [1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0])

    def test_prepared_terrain(self):
        path = os.path.join(self.dir, "hill.asc")
        center = write_hill(path)
        frame = terrain.Frame(center, 2000.0, 2000.0, 2000.0, 2000.0)
        plain = terrain.prepare(path, "auto", "auto", frame, 400.0, None, "lowest", 0.0, 0.0, "lowest")
        self.assertEqual(plain.lowest, 0.0)
        self.assertAlmostEqual(plain.highest, 100.0, places=2)      # the middle pixel is on the hill's top
        self.assertAlmostEqual(plain.spacing, 100.0, places=6)
        self.assertLessEqual(plain.x.min(), -2000.0)
        self.assertGreaterEqual(plain.x.max(), 2000.0)
        kept = terrain.prepare(path, "auto", "auto", frame, 400.0, None, "none", 0.0, 0.0, "lowest")
        self.assertAlmostEqual(kept.lowest, 1200.0, places=2)
        level = terrain.prepare(path, "auto", "auto", frame, 400.0, None, 1000.0, 0.0, 0.0, "lowest")
        self.assertAlmostEqual(level.lowest, 200.0, places=2)
        # a flat border at 1300 m in the file's datum: 100 m above the plain, so the edge points sit at 100 m
        flat = terrain.prepare(path, "auto", "auto", frame, 400.0, None, 1200.0, 300.0, 500.0, 1300.0)
        edge = (np.abs(flat.x) > 1750.0) | (np.abs(flat.y) > 1750.0)
        self.assertTrue(np.allclose(flat.z[edge], 100.0))
        self.assertAlmostEqual(flat.flat_level, 100.0, places=10)
        thinned = terrain.prepare(path, "auto", "auto", frame, 400.0, 200.0, "lowest", 0.0, 0.0, "lowest")
        self.assertAlmostEqual(thinned.spacing, 200.0, places=6)

    def test_data_that_does_not_cover_the_domain_is_refused(self):
        path = os.path.join(self.dir, "hill.asc")
        center = write_hill(path)
        frame = terrain.Frame((center[0] + 3000.0, center[1]), 2000.0, 2000.0, 2000.0, 2000.0)
        with self.assertRaisesRegex(CaseError, r"terrain.file = .* does not cover the domain on the east"):
            terrain.prepare(path, "auto", "auto", frame, 400.0, None, "lowest", 0.0, 0.0, "lowest")

    def test_missing_data_in_the_domain_is_refused(self):
        path = os.path.join(self.dir, "hill.asc")
        center = write_hill(path, nodata_at=(40, 40))
        frame = terrain.Frame(center, 2000.0, 2000.0, 2000.0, 2000.0)
        with self.assertRaisesRegex(CaseError, "terrain.file = .*1 of the .* points in the domain have no data"):
            terrain.prepare(path, "auto", "auto", frame, 400.0, None, "lowest", 0.0, 0.0, "lowest")

    def test_geographic_raster(self):
        # 0.001 degrees: 111 m in latitude; the hill's top is at 43.415 N, 113.005 W
        lon = -113.05 + 0.001 * np.arange(90)
        lat = 43.37 + 0.001 * np.arange(90)
        glon, glat = np.meshgrid(lon, lat)
        z = 1500.0 + 150.0 * np.exp(-(((glon + 113.005) * 80000.0) ** 2 + ((glat - 43.415) * 111000.0) ** 2) / (2 * 500.0 ** 2))
        path = os.path.join(self.dir, "geo.asc")
        with open(path, "w", encoding="utf-8") as out:
            out.write("ncols 90\nnrows 90\nxllcenter -113.05\nyllcenter 43.37\ncellsize 0.001\n")
            np.savetxt(out, z[::-1], fmt="%.3f")
        raw = changed(base_case(), domain={"center_lat_lon": [43.415, -113.005], "half_width": [1500.0], "top": 1500.0},
                      terrain={"file": path}, mesh__cell_size=125.0)
        result = self.build(raw)
        data = np.loadtxt(os.path.join(result.case_dir, "terrain.xyz"))
        top = data[np.argmax(data[:, 2])]
        self.assertLess(math.hypot(top[0], top[1]), 120.0)        # the hill's top is at the centre of the domain
        self.assertAlmostEqual(data[:, 2].max(), 150.0, delta=3.0)
        self.assertEqual(data[:, 2].min(), 0.0)
        problems = self.refused(changed(raw, domain={"half_width": [1500.0], "top": 1500.0}))
        self.assertTrue(any("domain.center_lat_lon, or terrain.coordinates = projected" in p for p in problems), problems)


class RoughnessTests(CaseDir):
    def test_esa_cci_table(self):
        # amrTerrain, src/backend/backendInterface.py: lccs_class and roughness_length, with the five changes of
        # roughness.py's docstring (class 0 removed; 12, 160, 170 and 200-202 to Wieringa's 1992 Davenport classes)
        codes = [10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 110, 120, 121, 122, 130,
                 140, 150, 151, 152, 153, 160, 170, 180, 190, 200, 201, 202, 210, 220]
        z0 = [0.1, 0.1, 0.25, 0.2, 0.3, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.15, 0.3, 0.07,
              0.07, 0.07, 0.07, 0.05, 0.05, 0.05, 0.05, 0.05, 1.0, 1.0, 0.2, 0.8, 0.005, 0.005, 0.005, 1e-4, 0.0012]
        self.assertEqual(roughness.class_table("esa_cci"), dict(zip(codes, z0)))
        with self.assertRaisesRegex(CaseError, "holds the classes 0, "):      # no data is not a land cover
            roughness.classes_to_roughness([0, 10], roughness.class_table("esa_cci"))

    def test_classes_become_roughness(self):
        table = roughness.class_table({1: 0.03, "2": 0.5})
        self.assertEqual(roughness.classes_to_roughness([1, 2, 2.0, 1], table).tolist(), [0.03, 0.5, 0.5, 0.03])
        self.assertEqual(roughness.classes_to_roughness([1, 9], table, default=0.1).tolist(), [0.03, 0.1])
        with self.assertRaisesRegex(CaseError, r"surface.roughness_map.file holds the classes 9, "):
            roughness.classes_to_roughness([1, 9], table)
        with self.assertRaisesRegex(CaseError, "finite and > 0"):
            roughness.class_table({1: 0.0})

    def test_land_cover_raster_becomes_the_roughness_file(self):
        hill = os.path.join(self.dir, "hill.asc")
        center = write_hill(hill)
        cover = os.path.join(self.dir, "cover.asc")
        codes = np.where(np.arange(81)[None, :] < 40, 130, 70) * np.ones((81, 1), dtype=int)
        with open(cover, "w", encoding="utf-8") as out:
            out.write("ncols 81\nnrows 81\nxllcorner 500000\nyllcorner 4000000\ncellsize 100\n")
            np.savetxt(out, codes, fmt="%d")
        raw = changed(base_case(), domain={"center": list(center), "half_width": [2000.0], "top": 1500.0},
                      terrain={"file": hill}, surface={"roughness_map": {"file": cover}},
                      wind__upstream_roughness=0.07)
        deck, result = self.deck(raw)
        self.assertEqual(deck["bc.surface_z0_file"], ["roughness.xyz"])
        self.assertNotIn("bc.surface_z0", deck)
        data = np.loadtxt(os.path.join(result.case_dir, "roughness.xyz"))
        self.assertTrue(np.all(data[data[:, 0] < -100.0, 2] == 0.07))     # shrubland (130) in the west
        self.assertTrue(np.all(data[data[:, 0] > 100.0, 2] == 1.0))       # needle-leaved forest (70) in the east


class SchemaTests(CaseDir):
    def test_unknown_keys_are_refused_with_a_suggestion(self):
        with self.assertRaises(CaseError) as caught:
            validate(changed(base_case(), wind__sped=8.0, mesh__layer=3, weather={}))
        text = "\n".join(caught.exception.problems)
        self.assertIn("wind.sped", text)
        self.assertIn("wind.speed", text)        # the suggestion
        self.assertIn("mesh.layer", text)
        self.assertIn("weather", text)

    def test_ranges_abort_and_name_the_key(self):
        for change, key in (({"wind__speed": 0.0}, "wind.speed"), ({"wind__direction": 400.0}, "wind.direction"),
                            ({"mesh__layers": 0}, "mesh.layers"), ({"surface__roughness": -1.0}, "surface.roughness"),
                            ({"wind__drive": "jet"}, "wind.drive"), ({"mesh__layers": 2.5}, "mesh.layers"),
                            ({"domain__top": "high"}, "domain.top"), ({"solver__cfl_max": 0.0}, "solver.cfl_max")):
            with self.assertRaises(CaseError, msg=key) as caught:
                validate(changed(base_case(), **change))
            self.assertTrue(any(p.startswith(key + " =") for p in caught.exception.problems),
                            "%s: %s" % (key, caught.exception.problems))

    def test_required_keys(self):
        with self.assertRaises(CaseError) as caught:
            validate({"wind": {"speed": 8.0}})
        text = "\n".join(caught.exception.problems)
        self.assertIn("domain.top", text)
        self.assertIn("mesh.cell_size", text)

    def test_defaults_are_filled_in_and_given_keys_remembered(self):
        case, given = validate(base_case())
        self.assertEqual(case["solver"]["max_iterations"], 20000)
        self.assertEqual(case["atmosphere"]["theta"], "neutral")
        self.assertIn("wind.speed", given)
        self.assertNotIn("solver.max_iterations", given)

    def test_yaml_reader(self):
        path = os.path.join(self.dir, "c.yaml")
        with open(path, "w", encoding="utf-8") as out:
            out.write("wind:\n  speed: 8.\n  speed: 9.\n")
        with self.assertRaisesRegex(CaseError, "speed.*twice|twice.*speed"):
            schema.load_yaml(path)
        with open(path, "w", encoding="utf-8") as out:
            out.write("domain: {half_width: 2000, top: 1e3}\nmesh: {cell_size: 250}\nmodel: {tke_floor: 1e-6}\n")
        case, _ = validate(schema.load_yaml(path))      # PyYAML reads 1e3 and 1e-6 as text
        self.assertEqual((case["domain"]["top"], case["model"]["tke_floor"]), (1000.0, 1.0e-6))
        with self.assertRaisesRegex(CaseError, "no such file|cannot read"):
            schema.load_yaml(os.path.join(self.dir, "absent.yaml"))

    def test_overrides(self):
        raw = base_case()
        apply_override(raw, "wind.speed=12")
        apply_override(raw, "domain.half_width=[1000, 3000]")
        apply_override(raw, "gust_inputs.steady.cfl0=2.")
        self.assertEqual(raw["wind"]["speed"], 12)
        self.assertEqual(raw["domain"]["half_width"], [1000, 3000])
        self.assertEqual(raw["gust_inputs"], {"steady.cfl0": 2.0})
        with self.assertRaisesRegex(CaseError, "section.key=value"):
            apply_override(raw, "wind.speed")

    def test_key_table_lists_every_key(self):
        table = key_table_rst()
        for name, _ in schema.iter_keys():
            self.assertIn("``%s``" % name, table)

    def test_documented_key_table_is_current(self):
        path = os.path.join(REPO, "docs", "case_setup_keys.inc")
        with open(path, encoding="utf-8") as handle:
            self.assertEqual(handle.read(), key_table_rst(),
                             "docs/case_setup_keys.inc is not what tools/make_gust_case.py --write-key-table writes")

    def test_yaml_reference_is_current(self):
        path = os.path.join(REPO, "Tests", "inputs_yaml_reference")
        with open(path, encoding="utf-8") as handle:
            self.assertEqual(handle.read(), key_reference_yaml(), "Tests/inputs_yaml_reference is not what "
                             "tools/make_gust_case.py --write-yaml-reference writes")

    def test_yaml_reference_uncommented_is_yaml_with_every_key_and_default(self):
        lines = key_reference_yaml().splitlines(True)
        self.assertTrue(all(line.startswith("#") or not line.strip() for line in lines))
        data = yaml.safe_load("".join(line[2:] if line.startswith("# ") else line for line in lines))
        for name, key in schema.iter_keys():
            node = data
            for part in name.split("."):
                node = node[part[:-2]][0] if part.endswith("[]") else node[part]
            if key.required or key.default is None:
                self.assertIsNone(node, name)
            elif key.default == "":
                self.assertEqual(node, "", name)
            elif isinstance(key.default, str) and not isinstance(node, str):
                self.assertEqual(" ".join("%g" % v for v in (node if isinstance(node, list) else [node])),
                                 key.default, name)
            else:
                self.assertEqual(node, list(key.default) if isinstance(key.default, tuple) else key.default, name)
        self.assertEqual(sorted(data), sorted(schema.SPEC))


class PlanTests(CaseDir):
    def test_cross_key_refusals_name_their_keys(self):
        inversion = {"base": 500.0, "thickness": 100.0, "strength": 5.0}
        cases = [
            ({"domain__west": 1000.0}, "domain.half_width and domain.west"),
            ({"domain__half_width": None}, "domain: the extent is required"),
            ({"domain__center": [0.0, 0.0], "domain__center_lat_lon": [43.0, -113.0]}, "domain.center and domain.center_lat_lon"),
            ({"wind__directions": [0.0, 90.0]}, "wind.direction"),
            ({"wind__direction": None}, "wind.direction"),
            ({"surface__roughness": None}, "surface.roughness"),
            ({"wind__profile": "table"}, "wind.table"),
            ({"wind__power_law_exponent": 0.2}, "wind.power_law_exponent"),
            ({"wind__drive": "geostrophic_inflow"}, "atmosphere.latitude"),
            ({"atmosphere__latitude": 45.0}, "atmosphere.latitude has no effect"),
            ({"wind__drive": "geostrophic_inflow", "atmosphere__latitude": 45.0, "atmosphere__geostrophic_wind": [10.0, 0.0]},
             "wind.speed and atmosphere.geostrophic_wind"),
            ({"wind__drive": "geostrophic_inflow", "atmosphere__latitude": 45.0, "wind__side_faces": "periodic"},
             "wind.side_faces has no effect with wind.drive"),
            ({"wind__drive": "geostrophic_column", "atmosphere__latitude": 45.0, "wind__side_faces": "slip_wall"},
             "wind.side_faces has no effect with wind.drive"),
            ({"wind__drive": "couette", "atmosphere__theta": "capped", "atmosphere__inversion": inversion,
              "atmosphere__lapse_rate": 0.003}, "wind.drive = couette needs atmosphere.theta = neutral"),
            ({"wind__drive": "couette", "wind__upstream_roughness": 0.1}, "wind.upstream_roughness = surface.roughness"),
            ({"atmosphere__theta": "capped"}, "atmosphere.theta = capped"),
            ({"atmosphere__theta": "surface_flux"}, "atmosphere.surface_heat_flux"),
            ({"atmosphere__surface_heat_flux": 0.02}, "atmosphere.surface_heat_flux"),
            ({"atmosphere__inversion": {"base": 500.0}}, "atmosphere.inversion"),
            ({"atmosphere__inversion": dict(inversion, base=950.0)}, "atmosphere.inversion"),
            ({"mesh__first_layer": 100.0}, "mesh.first_layer"),
            ({"mesh__stretch_ratio": 1.1, "mesh__first_layer": None}, "mesh.stretch_ratio"),
            ({"mesh__stretch_ratio": 1.5}, "domain.top"),
            ({"mesh__hierarchy_levels": 3, "mesh__max_grid_size": 30}, "mesh.hierarchy_levels"),
            ({"terrain__border": {"flat_width": 100.0}}, "terrain.border"),
            ({"terrain__smoothing": "stf"}, "terrain.smoothing"),
            ({"nests": [{"name": "a", "center": [0.0, 0.0], "half_width": [1900.0]}]}, "nests[0]"),
            ({"nests": [{"name": "a b", "center": [0.0, 0.0], "half_width": [500.0]}]}, "nests[0].name"),
            ({"nests": [{"name": "a", "half_width": [500.0]}]}, "nests[0].center"),
            ({"nests": [{"name": "a", "center": [0.0, 0.0], "half_width": [500.0]}], "wind__side_faces": "periodic"},
             "wind.side_faces = periodic"),
            ({"gust_inputs": {"nonsense": 1}}, "gust_inputs"),
            ({"case__name": "../up"}, "case.name = '../up'"),
            ({"case__name": ""}, "case.name = ''"),
            ({"case__name": "a/b"}, "case.name = 'a/b'"),
            ({"case__name": ".."}, "case.name = '..'"),
            ({"wind__direction": 225.0, "wind__side_faces": "periodic"}, "wind.side_faces has no effect"),
            ({"solver__initial_blend": [1.0, 0.1]}, "solver.initial_blend has no effect"),
            ({"nests": [{"name": "a", "center": [0.0, 0.0], "half_width": [500.0], "blend_width": 100.0}]},
             "nests[0].blend_width has no effect"),
        ]
        for change, text in cases:
            problems = self.refused(changed(base_case(), **change))
            self.assertTrue(any(text in p for p in problems), "%s: expected '%s' in %s" % (change, text, problems))

    def test_every_problem_is_reported_at_once(self):
        problems = self.refused(changed(base_case(), surface__roughness=None, wind__direction=None, mesh__first_layer=100.0))
        self.assertGreaterEqual(len(problems), 3)

    def test_grid_of_the_plan(self):
        case, given = validate(changed(base_case(), domain__half_width=[2000.0, 3000.0], mesh__cell_size=110.0))
        plan = resolve(case, given, self.dir)
        # 4000 / 110 = 36.4 and 6000 / 110 = 54.5 go to the nearest multiples of 8
        self.assertEqual((plan.nx, plan.ny, plan.layers), (40, 56, 16))
        self.assertAlmostEqual(plan.dy, 6000.0 / 56, places=10)
        self.assertAlmostEqual(grid.stretched_top(8.0, plan.stretch_ratio, 16), 1000.0, places=8)
        self.assertAlmostEqual(plan.ustar, 0.564336173406487, places=12)

    def test_steep_stretch_writes_its_ratio(self):
        """80 layers from 0.14 m to 3000 m: GUST's own Newton (from 1.03) fails, so the deck carries the ratio."""
        from gust_input_checker.context import gust_newton_converges
        raw = changed(base_case(), domain__top=3000.0, mesh__layers=80, mesh__first_layer=0.14)
        deck, result = self.deck(raw)
        ratio = float(deck["mesh.stretch_ratio"][0])
        self.assertEqual(float(deck["mesh.dz0"][0]), 0.14)
        self.assertAlmostEqual(ratio, 1.100725926, places=8)
        self.assertAlmostEqual(0.14 * (ratio ** 80 - 1.0) / (ratio - 1.0), 3000.0, places=6)
        self.assertFalse(gust_newton_converges(0.14, 3000.0, 80))

    def test_stretch_of_one_is_written_uniform(self):
        """48 layers of 6.4 m reach 307.2 m with no stretch: uniform layers, no mesh.dz0 (GUST refuses
        dz0 above top / n by any rounding)."""
        raw = changed(base_case(), domain__top=307.2, mesh__layers=48, mesh__first_layer=6.4)
        deck, result = self.deck(raw)
        self.assertNotIn("mesh.dz0", deck)
        self.assertNotIn("mesh.stretch_ratio", deck)
        self.assertEqual(deck["mesh.n_cell"][2], "48")

    def test_uniform_layers_past_the_top_are_refused(self):
        """first_layer 10 m, ratio 1: 10 layers pass a 95 m top; a 100 m top is accepted."""
        raw = changed(base_case(), domain__top=95.0, mesh__layers=None, mesh__first_layer=10.0, mesh__stretch_ratio=1.0)
        problems = self.refused(raw)
        self.assertTrue(any("pass it even without stretching" in p for p in problems), problems)
        deck, result = self.deck(changed(raw, domain__top=100.0))
        self.assertEqual(deck["mesh.n_cell"][2], "10")

    def test_patch_lies_on_parent_cell_faces(self):
        raw = changed(base_case(), mesh__cell_size=125.0,
                      nests=[{"name": "mid", "center": [130.0, -70.0], "half_width": [610.0, 480.0]}])
        case, given = validate(raw)
        plan = resolve(case, given, self.dir)
        patch = plan.patches[0]
        # faces at multiples of the parent's 125 m: 130 -/+ 610 -> -500, 750; -70 -/+ 480 -> -500, 375
        self.assertEqual((patch.x_lo, patch.x_hi, patch.y_lo, patch.y_hi), (-500.0, 750.0, -500.0, 375.0))
        self.assertEqual((patch.nx, patch.ny), (20, 14))        # twice the parent's 10 x 7 cells


class BuilderTests(CaseDir):
    def test_flat_deck(self):
        deck, result = self.deck(base_case())
        self.assertEqual(deck["mesh.n_cell"], ["16", "16", "16"])
        self.assertEqual(deck["mesh.prob_lo"], ["-2000.", "-2000.", "0."])
        self.assertEqual(deck["mesh.prob_hi"], ["2000.", "2000.", "1000."])
        self.assertEqual(deck["terrain.type"], ["flat"])
        self.assertEqual((deck["bc.xlo"], deck["bc.xhi"], deck["bc.ylo"], deck["bc.yhi"]),
                         (["inflow"], ["outflow"], ["slip_wall"], ["slip_wall"]))
        self.assertAlmostEqual(float(deck["bc.inflow_ustar"][0]), 0.564336173406487, places=9)
        self.assertEqual(deck["bc.inflow_direction"], ["0."])
        self.assertEqual((deck["abl.drive"], deck["abl.theta"]), (["inflow_only"], ["neutral"]))
        self.assertEqual(deck["forcing.gravity"], ["0."])
        self.assertNotIn("convection.hybrid_height", deck)       # no blend over flat ground
        self.assertEqual(sorted(os.listdir(result.case_dir)),
                         [".gust_case_data", ".gust_case_files", "case_resolved.yaml", "case_summary.txt", "inputs", "run.sh"])
        self.assertTrue(os.access(os.path.join(result.case_dir, "run.sh"), os.X_OK))

    def test_tool_decks_check_clean(self):
        """The column and nest decks the tool writes pass the deck checker: the column table and the
        parent checkpoint are files run.sh makes, given to the checker as made by the run."""
        raw = changed(base_case(), wind__drive="geostrophic_column", atmosphere__latitude=45.0,
                      nests=[{"name": "mid", "center": [0.0, 0.0], "half_width": [750.0]}])
        result = self.build(raw, check=True)
        rels = [rel for rel, _, _, _ in result.decks]
        self.assertEqual(rels, ["inputs", "column/inputs", "nest_mid/inputs"])
        self.assertEqual((result.errors, result.warnings, result.findings), (0, 0, []))

    def butte(self):
        """The mesh, wind and set-up of inputs_bsb_srtm, over a stand-in hill of the same relief scale."""
        x, y = np.meshgrid(np.arange(-400.0, 19200.0, 100.0), np.arange(-400.0, 20700.0, 100.0))
        z = 1500.0 + 786.0 * np.exp(-((x - 9345.0) ** 2 + (y - 10080.0) ** 2) / (2.0 * 1500.0 ** 2))
        np.savetxt(os.path.join(self.dir, "butte.xyz"), np.column_stack([x.ravel(), y.ravel(), z.ravel()]), fmt="%.3f")
        speed = 0.5963 / 0.41 * math.log(10.01 / 0.01)       # the deck's u* as a speed at 10 m
        return {"case": {"name": "bsb"},
                "domain": {"center": [9345.0, 10080.0], "half_width": [9345.0, 10080.0], "top": 2806.0,
                           "origin": "southwest"},
                "terrain": {"file": "butte.xyz"},
                "mesh": {"cell_size": 68.4, "layers": 40, "first_layer": 4.0},
                "wind": {"speed": speed, "height": 10.0, "direction": 222.0, "drive": "couette"},
                "surface": {"roughness": 0.01},
                "atmosphere": {"damping": {"start": 2000.0}},
                "model": {"blend_height": 1965.0}}

    def test_butte_deck_is_reproduced(self):
        ours, _ = self.deck(self.butte())
        theirs = read_deck(os.path.join(REPO, "Tests", "Real_Terrain", "WindNinja2019", "inputs_bsb_srtm"))
        defaults = reference_defaults()
        # the differences that are meant: the names of the set-up (the deck predates them), the terrain file (the
        # deck's comes on the command line), and the ratio (ours leaves GUST to find it from the top)
        meant = {"abl.drive": ["couette"], "abl.theta": ["neutral"], "terrain.file": ["terrain.xyz"]}
        for key in sorted(set(ours) | set(theirs)):
            if key in ours and key in theirs:
                self.assertTrue(same_tokens(ours[key], theirs[key]), "%s: %s, the deck has %s" % (key, ours[key], theirs[key]))
            elif key in ours:
                self.assertEqual(ours[key], meant.get(key), "%s is not in the deck" % key)
            elif key == "mesh.stretch_ratio":
                self.assertAlmostEqual(grid.stretch_ratio_for_top(4.0, 40, 2806.0), float(theirs[key][0]), places=5)
            else:
                self.assertTrue(same_tokens(theirs[key], defaults[key]),
                                "%s = %s of the deck is absent here and is not GUST's default %s"
                                % (key, theirs[key], defaults[key]))

    def test_blend_height_follows_the_relief(self):
        raw = self.butte()
        del raw["model"]
        deck, result = self.deck(raw)
        data = np.loadtxt(os.path.join(result.case_dir, "terrain.xyz"))
        self.assertAlmostEqual(float(deck["convection.hybrid_height"][0]), 2.5 * (data[:, 2].max() - data[:, 2].min()), places=2)
        self.assertEqual(deck["turb.realizability"], ["0.309"])
        self.assertEqual(deck["turb.nut_relax"], ["0.5"])
        plain, _ = self.deck(changed(raw, model={"lee_setting": "none"}, case__name="plain"))
        for key in ("turb.realizability", "turb.nut_relax", "convection.hybrid_height", "convection.hybrid_sigma_max"):
            self.assertNotIn(key, plain)

    def test_sweep_writes_one_folder_per_direction(self):
        raw = changed(base_case(), wind__direction=None, wind__directions=[0.0, 222.0, 22.5])
        result = self.build(raw)
        self.assertEqual(sorted(d for d in os.listdir(result.case_dir) if d.startswith("dir_")),
                         ["dir_000", "dir_022p5", "dir_222"])
        north = read_deck(os.path.join(result.case_dir, "dir_000", "inputs"))
        self.assertEqual((north["bc.yhi"], north["bc.ylo"], north["bc.inflow_direction"]), (["inflow"], ["outflow"], ["270."]))
        south_west = read_deck(os.path.join(result.case_dir, "dir_222", "inputs"))
        self.assertEqual((south_west["bc.xlo"], south_west["bc.ylo"], south_west["bc.inflow_direction"]),
                         (["inflow"], ["inflow"], ["48."]))
        with open(os.path.join(result.case_dir, "run.sh"), encoding="utf-8") as handle:
            script = handle.read()
        for folder in ("dir_000", "dir_022p5", "dir_222"):
            self.assertIn("cd %s\n" % folder, script)

    def test_couette_lid(self):
        deck, _ = self.deck(changed(base_case(), wind__drive="couette", wind__direction=222.0))
        lid = [float(v) for v in deck["bc.zhi_wall_velocity"]]
        speed = 0.564336173406487 / 0.41 * math.log(1000.03 / 0.03)      # the log law at the top
        self.assertEqual(deck["bc.zhi"], ["no_slip_wall"])
        self.assertAlmostEqual(lid[0], speed * math.cos(math.radians(48.0)), places=6)
        self.assertAlmostEqual(lid[1], speed * math.sin(math.radians(48.0)), places=6)
        self.assertEqual(lid[2], 0.0)
        self.assertEqual(deck["turb.max_geom_lscale"], ["1000000."])

    def test_capped_theta(self):
        raw = changed(base_case(), atmosphere={"theta": "capped", "lapse_rate": 0.003,
                                               "inversion": {"base": 500.0, "thickness": 100.0, "strength": 5.0}})
        deck, _ = self.deck(raw)
        self.assertEqual(deck["bc.theta_profile_z"], ["0.", "500.", "600.", "1000."])
        self.assertTrue(same_tokens(deck["bc.theta_profile"], ["300.", "300.", "305.", "306.2"], rel=1.0e-12))
        self.assertAlmostEqual(float(deck["bc.zhi_theta_grad"][0]), 0.003, places=12)
        self.assertEqual(deck["abl.theta"], ["capped"])
        self.assertEqual(deck["bc.rayleigh_zdamp"], ["700."])        # 0.7 of the top
        self.assertNotIn("forcing.gravity", deck)
        self.assertEqual(deck["bc.rayleigh_fields"], ["w", "theta"])

    def test_geostrophic_inflow_takes_the_column_runs_profile(self):
        raw = changed(base_case(), wind__drive="geostrophic_inflow", atmosphere={"latitude": 45.0})
        result = self.build(raw)
        deck = read_deck(os.path.join(result.case_dir, "inputs"))
        column = read_deck(os.path.join(result.case_dir, "column", "inputs"))
        self.assertEqual((deck["bc.inflow_profile"], deck["bc.inflow_file"]), (["file"], ["column/column_guess.dat"]))
        self.assertEqual(deck["abl.drive"], ["geostrophic_inflow"])
        self.assertNotIn("forcing.ug", deck)           # run.sh passes the column run's values
        # from 270 at 45 N the wind aloft crosses y towards -y: no slip walls (UT-13.1.1)
        self.assertEqual([deck["bc." + f][0] for f in ("xlo", "xhi", "ylo", "yhi")],
                         ["inflow", "outflow", "outflow", "inflow"])
        self.assertNotIn("mesh.is_periodic", deck)
        self.assertEqual(column["mesh.n_cell"], ["4", "4", "16"])
        self.assertEqual(column["mesh.is_periodic"], ["1", "1", "0"])
        self.assertEqual((column["init.column_target_u"], column["init.column_target_v"], column["init.column_target_height"]),
                         (["8."], ["0."], ["10."]))
        # the column starts from the drag law's speed for u* = 0.5643 m/s, z0 = 0.03 m at 45 N, along the wind
        start = 0.564336173406487 / 0.41 * math.hypot(math.log(0.564336173406487 / (1.0284451246e-4 * 0.03)) - 1.8, 4.5)
        self.assertAlmostEqual(float(column["forcing.ug"][0]), start, places=5)
        self.assertEqual(column["forcing.vg"], ["0."])
        with open(os.path.join(result.case_dir, "run.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('( cd column && $MPI "$GUST" inputs > log.txt )', script)
        self.assertIn('$MPI "$GUST" inputs forcing.ug="$ug" forcing.vg="$vg" "$@" > log.txt', script)

    def test_given_geostrophic_wind_goes_into_both_decks(self):
        raw = changed(base_case(), wind__drive="geostrophic_inflow", wind__speed=None, wind__height=None,
                      atmosphere={"latitude": 45.0, "geostrophic_wind": [10.0, -2.0]})
        result = self.build(raw)
        for rel in ("inputs", os.path.join("column", "inputs")):
            deck = read_deck(os.path.join(result.case_dir, rel))
            self.assertEqual((deck["forcing.ug"], deck["forcing.vg"]), (["10."], ["-2."]), rel)
        self.assertNotIn("init.column_target_u", read_deck(os.path.join(result.case_dir, "column", "inputs")))
        # the faces follow G's direction (-11.31 degrees) turned back by half the drag law's angle (u* 0.3835 m/s):
        # towards +0.88, so y-low is an inflow face, whatever wind.direction (270) says
        deck = read_deck(os.path.join(result.case_dir, "inputs"))
        self.assertEqual([deck["bc." + f][0] for f in ("xlo", "xhi", "ylo", "yhi")],
                         ["inflow", "outflow", "inflow", "outflow"])

    def test_given_geostrophic_wind_sets_the_faces(self):
        """UT-13.1.1: G = (10, 0) at 45 N: the wind near the ground comes from G turned back 12.2 degrees
        (counter-clockwise), so it enters through y-low; at 45 S through y-high. wind.direction does not matter."""
        for lat, ylo, yhi in ((45.0, "inflow", "outflow"), (-45.0, "outflow", "inflow")):
            for direction in (270.0, 300.0):
                raw = changed(base_case(), wind__drive="geostrophic_inflow", wind__speed=None, wind__height=None,
                              wind__direction=direction, atmosphere={"latitude": lat, "geostrophic_wind": [10.0, 0.0]})
                deck = read_deck(os.path.join(self.build(raw, overwrite=True).case_dir, "inputs"))
                self.assertEqual([deck["bc." + f][0] for f in ("xlo", "xhi", "ylo", "yhi")],
                                 ["inflow", "outflow", ylo, yhi], "latitude %g, direction %g" % (lat, direction))

    def test_southern_hemisphere_target_faces(self):
        """UT-13.1.1: 8 m/s at 10 m from 270 at 45 S: the wind aloft turns counter-clockwise (towards +y)."""
        raw = changed(base_case(), wind__drive="geostrophic_column", atmosphere={"latitude": -45.0})
        deck = read_deck(os.path.join(self.build(raw).case_dir, "inputs"))
        self.assertEqual([deck["bc." + f][0] for f in ("xlo", "xhi", "ylo", "yhi")],
                         ["inflow", "outflow", "inflow", "outflow"])

    def test_patch_deck(self):
        hill = os.path.join(self.dir, "hill.asc")
        center = write_hill(hill)
        raw = changed(base_case(), domain={"center": list(center), "half_width": [2000.0], "top": 1500.0},
                      terrain={"file": hill}, mesh__cell_size=125.0,
                      nests=[{"name": "crest", "center": list(center), "half_width": [600.0]}])
        result = self.build(raw)
        parent = read_deck(os.path.join(result.case_dir, "inputs"))
        patch = read_deck(os.path.join(result.case_dir, "nest_crest", "inputs"))
        self.assertNotIn("io.write_checkpoint", parent)      # GUST's default: the patch needs the checkpoint
        self.assertEqual(patch["mesh.n_cell"], ["20", "20", "16"])       # 1250 m of 62.5 m cells
        self.assertEqual((patch["mesh.prob_lo"], patch["mesh.prob_hi"]), (["-625.", "-625.", "0."], ["625.", "625.", "1500."]))
        self.assertEqual({patch["bc." + f][0] for f in ("xlo", "xhi", "ylo", "yhi")}, {"nest"})
        self.assertEqual(patch["init.first_guess"], ["nest"])
        self.assertEqual(patch["terrain.file"], ["../terrain_nest_crest.xyz"])
        self.assertLess(float(patch["terrain.nest_blend_width"][0]), 0.5 * 1250.0)
        self.assertNotIn("bc.nest_parent", patch)      # run.sh names the parent's last checkpoint
        self.assertNotIn("mesh.hierarchy_levels", patch)
        for key in ("mesh.dz0", "bc.surface_z0", "bc.theta_profile", "turb.realizability"):
            self.assertEqual(patch[key], parent[key], key)
        with open(os.path.join(result.case_dir, "run.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('( cd nest_crest && $MPI "$GUST" inputs bc.nest_parent="../$parent" "$@" > log.txt )', script)
        # the patch's ground is the parent's data: the same elevations inside the patch
        data = np.loadtxt(os.path.join(result.case_dir, "terrain_nest_crest.xyz"))
        self.assertAlmostEqual(data[:, 2].max(), np.loadtxt(os.path.join(result.case_dir, "terrain.xyz"))[:, 2].max(), places=3)

    @unittest.skipIf(os.name == "nt", "run.sh is a bash script")
    def test_run_script_hands_the_newest_checkpoint_to_the_patch(self):
        """A stale chk200000 of an earlier run is ignored, and chk100000 beats chk99999 (it sorts first as text)."""
        import subprocess
        raw = changed(base_case(), nests=[{"name": "mid", "center": [0.0, 0.0], "half_width": [750.0]}])
        result = self.build(raw)
        os.mkdir(os.path.join(result.case_dir, "chk200000"))             # an earlier run's checkpoint
        fake = os.path.join(self.dir, "fake_gust")
        with open(fake, "w") as out:                                     # parent: two checkpoints; patch: its arguments
            out.write('#!/bin/bash\nif [ "$(basename "$PWD")" = nest_mid ]; then echo "$@" > args.txt; '
                      'else mkdir -p chk99999 chk100000; fi\n')
        os.chmod(fake, 0o755)
        env = dict(os.environ, GUST=fake, MPI="")
        run = subprocess.run(["bash", os.path.join(result.case_dir, "run.sh")], env=env, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        with open(os.path.join(result.case_dir, "nest_mid", "args.txt")) as handle:
            self.assertEqual(handle.read().split(), ["inputs", "bc.nest_parent=../chk100000"])

    def test_gust_inputs_add_replace_and_remove_keys(self):
        deck, _ = self.deck(changed(base_case(), gust_inputs={"steady.cfl0": 2.0, "io.plot_file": "test_plt",
                                                              "steady.cfl_max": 8.0, "simplec.w_vnn_cap": None}))
        self.assertEqual((deck["steady.cfl0"], deck["io.plot_file"]), (["2."], ["test_plt"]))
        self.assertEqual(deck["steady.cfl_max"], ["8."])         # replaces the value of solver.cfl_max
        self.assertNotIn("simplec.w_vnn_cap", deck)              # null removes a key

    def test_existing_folder_needs_overwrite(self):
        self.build(base_case())
        with self.assertRaisesRegex(CaseError, "--overwrite"):
            self.build(base_case())
        self.build(base_case(), overwrite=True)

    def test_dry_run_writes_nothing(self):
        result = self.build(base_case(), dry_run=True)
        self.assertFalse(os.path.exists(result.case_dir))
        self.assertIn("16 x 16 x 16", result.summary)
        self.assertEqual([rel for rel, _, _, _ in result.decks], ["inputs"])


class DeckTests(unittest.TestCase):
    def test_values_are_written_as_gust_reads_them(self):
        self.assertEqual([format_value(v) for v in (1, 1.0, 0.5, True, False, "lusgs", [1, 2.0, 0.25], 1.0e-6, 1.0e6)],
                         ["1", "1.", "0.5", "1", "0", "lusgs", "1 2. 0.25", "1e-06", "1000000."])

    def test_a_key_is_set_once(self):
        deck = Deck()
        deck.set("a.b", 1)
        with self.assertRaises(KeyError):
            deck.set("a.b", 2)
        deck.replace("a.b", 2)
        self.assertEqual(deck.get("a.b"), "2")
        deck.remove("a.b")
        self.assertNotIn("a.b", deck)


class CommandLineTests(CaseDir):
    def run_tool(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = make_gust_case.main(list(argv))
        return status, out.getvalue(), err.getvalue()

    def case_file(self):
        path = os.path.join(self.dir, "case.yaml")
        with open(path, "w", encoding="utf-8") as out:
            out.write(make_gust_case.EXAMPLE)
        return path

    def test_example_is_a_valid_case(self):
        status, out, _ = self.run_tool("--example")
        self.assertEqual((status, out), (0, make_gust_case.EXAMPLE))
        status, _, err = self.run_tool(self.case_file(), "--quiet", "--no-check")
        self.assertEqual((status, err), (0, ""))
        deck = read_deck(os.path.join(self.dir, "flat_west", "inputs"))
        self.assertEqual(deck["mesh.n_cell"], ["80", "80", "32"])

    def test_options_replace_the_files_values(self):
        status, _, err = self.run_tool(self.case_file(), "--quiet", "--no-check", "--name", "other", "--direction", "180",
                                       "--speed", "5", "--cell-size", "250", "--layers", "16", "--max-iterations", "7",
                                       "--set", "solver.cfl_max=8", "--set", "gust_inputs.steady.cfl0=2.")
        self.assertEqual((status, err), (0, ""))
        deck = read_deck(os.path.join(self.dir, "other", "inputs"))
        self.assertEqual(deck["mesh.n_cell"], ["16", "16", "16"])
        self.assertEqual((deck["bc.ylo"], deck["bc.inflow_direction"]), (["inflow"], ["90."]))
        self.assertAlmostEqual(float(deck["bc.inflow_ustar"][0]), 0.41 * 5.0 / math.log(10.03 / 0.03), places=9)
        self.assertEqual((deck["steady.max_iterations"], deck["steady.cfl_max"], deck["steady.cfl0"]), (["7"], ["8."], ["2."]))

    def test_directions_option_makes_a_sweep(self):
        status, _, _ = self.run_tool(self.case_file(), "--quiet", "--no-check", "--directions", "0", "90", "--cell-size", "250")
        self.assertEqual(status, 0)
        self.assertEqual(sorted(d for d in os.listdir(os.path.join(self.dir, "flat_west")) if d.startswith("dir_")),
                         ["dir_000", "dir_090"])
        with self.assertRaises(SystemExit) as caught, contextlib.redirect_stderr(io.StringIO()) as err:
            make_gust_case.main([self.case_file(), "--direction", "180", "--directions", "0", "90"])
        self.assertEqual(caught.exception.code, 2)
        self.assertIn("--direction and --directions", err.getvalue())

    def test_overwrite_removes_the_previous_write(self):
        """A sweep, then one direction with --overwrite: the dir_* decks go, a run output stays."""
        case = self.case_file()
        self.assertEqual(self.run_tool(case, "--quiet", "--no-check", "--directions", "0", "90")[0], 0)
        folder = os.path.join(self.dir, "flat_west")
        with open(os.path.join(folder, "dir_090", "log.txt"), "w") as out:        # a run of the first write
            out.write("run\n")
        self.assertEqual(self.run_tool(case, "--quiet", "--no-check", "--direction", "270")[0], 2)     # no --overwrite
        self.assertEqual(self.run_tool(case, "--quiet", "--no-check", "--direction", "270", "--overwrite")[0], 0)
        self.assertFalse(os.path.exists(os.path.join(folder, "dir_000")))
        self.assertEqual(os.listdir(os.path.join(folder, "dir_090")), ["log.txt"])
        self.assertTrue(os.path.isfile(os.path.join(folder, "inputs")))

    def test_refused_case_exits_2_and_writes_nothing(self):
        status, out, err = self.run_tool(self.case_file(), "--speed", "-3", "--set", "mesh.layer=4")
        self.assertEqual(status, 2)
        self.assertIn("wind.speed = -3", err)
        self.assertIn("mesh.layer", err)
        self.assertIn("nothing was written", err)
        self.assertFalse(os.path.exists(os.path.join(self.dir, "flat_west")))
        self.assertEqual(self.run_tool()[0], 2)

    def test_dry_run_prints_the_deck(self):
        status, out, _ = self.run_tool(self.case_file(), "--dry-run")
        self.assertEqual(status, 0)
        self.assertIn("==== inputs (not written)", out)
        self.assertIn("bc.inflow_profile", out)
        self.assertFalse(os.path.exists(os.path.join(self.dir, "flat_west")))

    def test_list_keys(self):
        status, out, _ = self.run_tool("--list-keys")
        self.assertEqual(status, 0)
        for name, _ in schema.iter_keys():
            self.assertIn(name, out)


class LibraryTests(CaseDir):
    """library.json from the column runs (make_gust_case.py --collect). The column outputs are written here
    with a profile linear in ln((z + z0) / z0) in u and v, W(z) = c1 L(z) + c2, so the wind at any height has a
    closed form; c2 != 0 makes the direction turn with height."""

    Z0, HEIGHT = 0.03, 10.0
    L10 = math.log((HEIGHT + Z0) / Z0)
    C1, C2 = complex(0.4 / L10, -8.0 / L10), complex(-0.4, 0.0)  # W(10 m) = (0, -8): 8 m/s from 0 degrees (north)
    G = complex(-4.0, -11.0)                                     # the column's geostrophic wind for the wind from 0
    LEVELS = (2.0, 6.0, 12.0, 30.0, 200.0, 1500.0)

    def W(self, z, turn):
        L = math.log((z + self.Z0) / self.Z0)
        return (self.C1 * L + self.C2) * complex(math.cos(math.radians(turn)), math.sin(math.radians(turn)))

    def write_column(self, folder, turn, ustar=0.4, adjusted=True, inv_L=0.0, heat_flux=0.0, scale=1.0, row=None):
        column = os.path.join(folder, "column")
        os.makedirs(column, exist_ok=True)
        with open(os.path.join(column, "column_guess.dat"), "w", encoding="utf-8") as out:
            out.write("# the converged column (init.first_guess = column): z u v theta k; u* %.6e m/s, 12 steps; "
                      "surface heat flux %.6e K m/s, 1/L %.6e 1/m, surface theta 300.000000 K\n"
                      % (ustar, heat_flux, inv_L))
            out.write("0 0 0 300 0.1\n")
            for z in self.LEVELS:
                w = self.W(z, turn) * scale
                out.write("%.17g %.17g %.17g 300 0.1\n" % (z, w.real, w.imag))
            if row is not None:
                out.write(row)
        g = self.G * complex(math.cos(math.radians(turn)), math.sin(math.radians(turn)))
        with open(os.path.join(column, "log.txt"), "w", encoding="utf-8") as out:
            if adjusted:
                # GUST prints G to 6 digits in this line and to 9 in the restart lines below it (App.cpp)
                out.write("gust: column guess: target wind reached after 3 forcing adjustments; the run takes geostrophic "
                          "wind (%.6g, %.6g) m/s and pressure gradient (0, 0) m/s^2\n" % (g.real, g.imag))
                out.write("gust: column guess: a restart from this run's checkpoints needs these inputs in place of the "
                          "init.* keys:\n  forcing.ug = %.9g\n  forcing.vg = %.9g\n" % (g.real, g.imag))
            else:
                out.write("gust: an unfinished column run\n  forcing.ug = 1\n  forcing.vg = 2\n")
        return g

    def sweep(self, directions=(0.0, 135.0), latitude=45.0, overwrite=False):
        raw = changed(base_case(), wind__direction=None,
                      wind__directions=list(directions) if isinstance(directions, (list, tuple)) else directions,
                      wind__drive="geostrophic_column", atmosphere={"latitude": latitude})
        return self.build(raw, overwrite=overwrite).case_dir

    def one(self, **write):
        """A one-direction sweep (from 0) with its column written; returns (case folder, its library entry)."""
        case_dir = self.sweep(directions=[0.0], overwrite=True)
        self.write_column(os.path.join(case_dir, "dir_000"), 0.0, **write)
        return case_dir, library.collect(case_dir)["directions"][0]

    def test_wind_at_height_g_and_turning_closed_form(self):
        case_dir = self.sweep(directions=[0.0, 135.0])
        # from 135 the flow turns by towards(135) - towards(0) = 135 - 270 = -135 degrees
        for direction, turn in ((0.0, 0.0), (135.0, -135.0)):
            self.write_column(os.path.join(case_dir, "dir_" + wind.direction_label(direction)), turn)
        lib = library.collect(case_dir)
        w10 = self.W(self.HEIGHT, 0.0)
        self.assertAlmostEqual(abs(w10 - complex(0.0, -8.0)), 0.0, places=12)
        # by hand: W(10 m) points to -y, G to (-4, -11): G is atan(4/11) clockwise of it
        turning = -math.degrees(math.atan2(4.0, 11.0))
        # the surface stress: the wind at the first centre (2 m), W(2) = c1 L(2) + c2
        w2 = self.W(2.0, 0.0)
        stress = math.degrees(math.atan2(self.G.imag, self.G.real) - math.atan2(w2.imag, w2.real))
        for e, direction, turn in zip(lib["directions"], (0.0, 135.0), (0.0, -135.0)):
            g = self.G * complex(math.cos(math.radians(turn)), math.sin(math.radians(turn)))
            self.assertEqual(e["direction"], direction)
            self.assertAlmostEqual(e["wind_at_height"]["speed"], 8.0, places=12)
            self.assertAlmostEqual(library._angle(e["wind_at_height"]["direction"] - direction), 0.0, places=9)
            self.assertAlmostEqual(e["geostrophic_wind"]["ug"], g.real, delta=1.0e-8 * abs(g))   # the 9 digits
            self.assertAlmostEqual(e["geostrophic_wind"]["speed"], abs(self.G), delta=1.0e-8 * abs(self.G))
            self.assertAlmostEqual(e["turning_deg"], turning, delta=1.0e-6)
            self.assertAlmostEqual(e["turning_from_stress_deg"], library._angle(stress), delta=1.0e-6)
            self.assertTrue(e["geostrophic_wind"]["adjusted_by_column"])
            self.assertEqual((e["roughness"], e["latitude"]), (0.03, 45.0))     # from the column deck
        # G from 0 points to (-4, -11): it comes from 20.0 degrees clockwise of north (east of north)
        self.assertAlmostEqual(lib["directions"][0]["geostrophic_wind"]["direction"],
                               math.degrees(math.atan2(4.0, 11.0)), places=6)
        self.assertEqual([e["folder"] for e in lib["directions"]], ["dir_000", "dir_135"])

    def test_below_the_first_cell_centre_the_log_law(self):
        rows = [[0.0, 0, 0, 300, 0.1]] + [[z, self.W(z, 0).real, self.W(z, 0).imag, 300, 0.1] for z in self.LEVELS]
        # GUST's wind_at below zc[0] (neutral): W(z1) ln((z + z0)/z0) / ln((z1 + z0)/z0)
        f = math.log((1.0 + self.Z0) / self.Z0) / math.log((2.0 + self.Z0) / self.Z0)
        u, v = library.wind_at_height(rows, 1.0, self.Z0)
        self.assertAlmostEqual(u, f * self.W(2.0, 0).real, places=12)
        self.assertAlmostEqual(v, f * self.W(2.0, 0).imag, places=12)
        with self.assertRaisesRegex(CaseError, "stratified"):           # psi_m is not evaluated here
            library.wind_at_height(rows, 1.0, self.Z0, inv_L=0.01)
        u, v = library.wind_at_height(rows, 12.0, self.Z0, inv_L=0.01)  # a level: no psi_m needed above z1
        self.assertAlmostEqual(complex(u, v), self.W(12.0, 0), places=12)
        for height in (0.0, 1600.0):
            with self.assertRaisesRegex(CaseError, "outside"):
                library.wind_at_height(rows, height, self.Z0)

    def test_drag_law_numbers_for_the_columns_ustar(self):
        d = self.one(ustar=0.4)[1]["drag_law"]
        # by hand: f = 2 (2 pi / 86400) sin 45, G = u*/kappa sqrt((ln(u*/(f z0)) - 1.8)^2 + 4.5^2)
        self.assertAlmostEqual(d["speed"], 10.673988565180437, places=9)
        self.assertAlmostEqual(d["turning_deg"], -24.28671580017051, places=9)   # clockwise in the north

    def test_stratified_column_has_no_drag_law(self):
        for header in ({"inv_L": 0.004}, {"heat_flux": -0.01}):
            e = self.one(**header)[1]
            self.assertIsNone(e["drag_law"])
            self.assertIsNotNone(e["turning_deg"])

    def test_southern_hemisphere_turns_the_other_way(self):
        case_dir = self.sweep(directions=[0.0], latitude=-30.0)
        self.write_column(os.path.join(case_dir, "dir_000"), 0.0)
        self.assertGreater(library.collect(case_dir)["directions"][0]["drag_law"]["turning_deg"], 0.0)

    def test_latitude_from_center_lat_lon(self):
        # atmosphere.latitude left out: the case tool takes domain.center_lat_lon's, and so must the library
        raw = changed(base_case(), wind__direction=None, wind__directions=[0.0], wind__drive="geostrophic_column",
                      domain__center_lat_lon=[40.0, -105.0])
        case_dir = self.build(raw).case_dir
        self.write_column(os.path.join(case_dir, "dir_000"), 0.0)
        e = library.collect(case_dir)["directions"][0]
        self.assertEqual(e["latitude"], 40.0)
        self.assertIsNotNone(e["drag_law"])

    def test_start_stop_step_folders_and_single_direction(self):
        case_dir = self.sweep(directions={"start": 0.0, "stop": 45.0, "step": 22.5})
        for d in (0.0, 22.5, 45.0):
            self.write_column(os.path.join(case_dir, "dir_" + wind.direction_label(d)), -d)
        lib = library.collect(case_dir)
        self.assertEqual([e["folder"] for e in lib["directions"]], ["dir_000", "dir_022p5", "dir_045"])
        self.assertEqual([e["direction"] for e in lib["directions"]], [0.0, 22.5, 45.0])
        raw = changed(base_case(), wind__direction=0.0, wind__drive="geostrophic_column", atmosphere={"latitude": 45.0},
                      case__name="single")
        single = self.build(raw).case_dir
        self.write_column(single, 0.0)
        e = library.collect(single)["directions"][0]
        self.assertEqual((e["folder"], e["direction"]), (".", 0.0))
        self.assertAlmostEqual(e["turning_deg"], -math.degrees(math.atan2(4.0, 11.0)), delta=1.0e-6)

    def test_given_geostrophic_wind_read_from_the_column_deck(self):
        raw = changed(base_case(), wind__drive="geostrophic_column", wind__speed=None, wind__height=None,
                      atmosphere={"latitude": 45.0, "geostrophic_wind": [3.0, -9.0]})
        case_dir = self.build(raw).case_dir
        self.write_column(case_dir, 0.0, adjusted=False)
        e = library.collect(case_dir)["directions"][0]
        self.assertEqual((e["folder"], e["geostrophic_wind"]["ug"], e["geostrophic_wind"]["vg"]), (".", 3.0, -9.0))
        self.assertFalse(e["geostrophic_wind"]["adjusted_by_column"])
        self.assertIsNone(e["wind_at_height"])                           # no wind.height: no reference wind
        self.assertIsNone(e["turning_deg"])

    def test_write_replaces_the_summary_section(self):
        case_dir = self.one()[0]
        summary = os.path.join(case_dir, "case_summary.txt")
        with open(summary, "a", encoding="utf-8") as out:                # a mention inside a line is not the section
            out.write("note: see %s after run.sh\n" % library.SUMMARY_MARK)
        for _ in range(2):
            library.write(case_dir, library.collect(case_dir))
        with open(summary, encoding="utf-8") as fh:
            text = fh.read()
        self.assertEqual(text.count("\n" + library.SUMMARY_MARK + "\n"), 1)
        self.assertIn("note: see %s after run.sh\n" % library.SUMMARY_MARK, text)
        self.assertTrue(text.startswith("Case t"))                     # the build's summary is kept
        with open(os.path.join(case_dir, "library.json"), encoding="utf-8") as fh:
            self.assertEqual(json.load(fh)["format_version"], 1)

    def test_refusals(self):
        case_dir = self.sweep(directions=[0.0])
        folder = os.path.join(case_dir, "dir_000")
        with self.assertRaisesRegex(CaseError, "run run.sh first"):      # not run yet
            library.collect(case_dir)
        for write, why in (({"adjusted": False}, "target wind reached"),    # an unfinished column: not the deck's G
                           ({"scale": 1.01}, "not the 8 m/s"),               # a table from another run
                           ({"row": "2000 nan 0 300 0.1\n"}, "not finite"),
                           ({"row": "1500 1 0 300 0.1\n"}, "do not increase"),
                           ({"row": "2000 1 0\n"}, "expected 5")):
            self.write_column(folder, 0.0, **write)
            with self.assertRaisesRegex(CaseError, why):
                library.collect(case_dir)
        self.write_column(folder, 0.0)
        with self.assertRaises(ValueError):                              # NaN never reaches the json
            library.write(case_dir, {"x": float("nan")})
        flat = self.build(changed(base_case(), case__name="inflow")).case_dir
        with self.assertRaisesRegex(CaseError, "no column run"):
            library.collect(flat)

    def test_run_script_collects_for_geostrophic_drives_only(self):
        with open(os.path.join(self.sweep(), "run.sh"), encoding="utf-8") as fh:
            script = fh.read()
        self.assertIn('--collect . --quiet || echo', script)               # a failed collect leaves the runs' status
        self.assertLess(script.index("dir_135"), script.index("--collect"))   # after the last run
        self.assertLess(script.index("rm -f library.json"), script.index("dir_000"))   # the old one goes first
        with open(os.path.join(self.build(changed(base_case(), case__name="inflow")).case_dir, "run.sh"),
                  encoding="utf-8") as fh:
            self.assertNotIn("library.json", fh.read())

    def test_overwrite_removes_the_previous_library(self):
        case_dir = self.one()[0]
        library.write(case_dir, library.collect(case_dir))
        raw = changed(base_case(), wind__direction=None, wind__directions=[0.0], wind__drive="geostrophic_column",
                      atmosphere={"latitude": 45.0})
        self.build(raw, overwrite=True)
        self.assertFalse(os.path.exists(os.path.join(case_dir, "library.json")))

    def test_command_line_collect(self):
        case_dir = self.one()[0]
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = make_gust_case.main(["--collect", case_dir])
        self.assertEqual((status, err.getvalue()), (0, ""))
        self.assertIn(library.SUMMARY_MARK, out.getvalue())
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            self.assertEqual(make_gust_case.main(["--collect", self.dir]), 2)   # not a case folder
        for extra in (["case.yaml"], ["--overwrite"], ["--speed", "8"]):
            with contextlib.redirect_stderr(io.StringIO()) as e, self.assertRaises(SystemExit) as caught:
                make_gust_case.main(["--collect", case_dir] + extra)
            self.assertEqual(caught.exception.code, 2)
            self.assertIn("not with", e.getvalue())


if __name__ == "__main__":
    unittest.main()
