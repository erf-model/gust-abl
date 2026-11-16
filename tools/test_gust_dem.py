"""Tests of the terrain download, tools/gust_case/dem.py and tools/fetch_terrain.py (CTest gust_dem_selftest).

No test touches the network: the tiles are made here, gzip-compressed, and
served by a stand-in for urllib's urlopen. Each tile holds the closed form

    z = 2 i + j - 1500,   i = 1200 (lat - 57),   j = 1200 (lon + 8)   (i, j rounded to whole points)

so every height in the written grid can be checked against the point's
latitude and longitude, with no ties in the rounding; the form goes below
0 m in the south-west, which tests the sea-floor option.

    cd tools && python3 -m unittest -v test_gust_dem
"""
import contextlib
import gzip
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fetch_terrain  # noqa: E402
from gust_case import CaseError, validate  # noqa: E402
from gust_case import dem, terrain  # noqa: E402
from gust_case.plan import _resolve_domain  # noqa: E402

N = 1201          # a 3 arc-second tile keeps the tests small (2.9 MB per tile)


def form(lat, lon):
    return 2.0 * np.rint(1200.0 * (lat - 57.0)) + np.rint(1200.0 * (lon + 8.0)) - 1500.0


def make_tile(lat0, lon0, n=N, void_at=None):
    """The .hgt bytes of the tile with south-west corner (lat0, lon0): rows north to south."""
    lat = lat0 + 1.0 - np.arange(n) / (n - 1.0)
    lon = lon0 + np.arange(n) / (n - 1.0)
    z = form(lat[:, None], lon[None, :]).astype(">i2")
    if void_at is not None:
        z[void_at] = dem.VOID
    return z.tobytes()


class Server:
    """A stand-in for urlopen: serves gzip tiles by URL and counts the requests."""

    def __init__(self, tiles=None, broken=()):
        self.tiles = tiles or {}
        self.broken = set(broken)
        self.requests = []

    def __call__(self, url, timeout):
        self.requests.append(url)
        if url in self.broken:
            return io.BytesIO(b"\x1f\x8b not a whole gzip stream")
        for (a, o), data in self.tiles.items():
            if url == dem.tile_url(a, o):
                return io.BytesIO(gzip.compress(data, compresslevel=1))
        import urllib.error
        raise urllib.error.URLError("404 for %s" % url)


class TmpDir(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name
        self.cache = os.path.join(self.dir, "cache")

    def tearDown(self):
        self._tmp.cleanup()


class TileTests(unittest.TestCase):
    def test_tile_names(self):
        self.assertEqual(dem.tile_name(57, -8), "N57W008")
        self.assertEqual(dem.tile_name(-1, -1), "S01W001")
        self.assertEqual(dem.tile_name(0, 0), "N00E000")
        self.assertEqual(dem.tile_name(-34, 151), "S34E151")
        self.assertEqual(dem.tile_url(57, -8),
                         "https://s3.amazonaws.com/elevation-tiles-prod/skadi/N57/N57W008.hgt.gz")

    def test_tiles_for_box(self):
        self.assertEqual(dem.tiles_for_box(57.1, 57.3, -7.5, -7.2), [(57, -8)])
        self.assertEqual(dem.tiles_for_box(56.9, 57.3, -8.2, -7.2), [(56, -9), (56, -8), (57, -9), (57, -8)])
        # a box ending on a whole degree does not need the next tile: its edge row is in the tile below
        self.assertEqual(dem.tiles_for_box(57.5, 58.0, -7.5, -7.0), [(57, -8)])
        # a tile whose only point in the box is its shared edge row or column is not fetched
        self.assertEqual(dem.tiles_for_box(57.5, 58.0001, -7.5, -7.0), [(57, -8)])
        self.assertEqual(dem.tiles_for_box(56.9999, 57.5, -8.0001, -7.0), [(57, -8)])
        self.assertEqual(dem.tiles_for_box(57.5, 58.0 + 1.0 / 3600.0, -7.5, -7.0), [(57, -8), (58, -8)])
        self.assertEqual(dem.tiles_for_box(56.99995, 57.00005, -7.5, -7.4), [(57, -8)])
        self.assertEqual(dem.tiles_for_box(-0.5, 0.5, -0.5, 0.5), [(-1, -1), (-1, 0), (0, -1), (0, 0)])

    def test_tiles_for_box_refusals(self):
        with self.assertRaisesRegex(CaseError, "180th meridian"):
            dem.tiles_for_box(10.0, 11.0, 179.5, -179.5)
        with self.assertRaisesRegex(CaseError, "latitudes"):
            dem.tiles_for_box(11.0, 10.0, 0.0, 1.0)
        with self.assertRaisesRegex(CaseError, "at most 16"):
            dem.tiles_for_box(10.0, 15.0, 0.0, 5.0)

    def test_read_hgt(self):
        z = dem.read_hgt(make_tile(57, -8))
        self.assertEqual(z.shape, (N, N))
        self.assertEqual(int(z[0, 0]), form(58.0, -8.0))          # first row north, first column west
        self.assertEqual(int(z[-1, -1]), form(57.0, -7.0))
        with self.assertRaisesRegex(CaseError, "not a 3601 x 3601 or 1201 x 1201"):
            dem.read_hgt(b"\0" * (2 * 1000 * 1000))


class MosaicTests(TmpDir):
    def tiles(self, corners, **kw):
        files = {}
        os.makedirs(self.cache, exist_ok=True)
        for a, o in corners:
            path = os.path.join(self.cache, dem.tile_name(a, o) + ".hgt")
            with open(path, "wb") as out:
                out.write(make_tile(a, o, **kw))
            files[(a, o)] = path
        return files

    def test_four_tiles_match_the_closed_form(self):
        box = (56.8, 57.2, -8.3, -7.6)
        files = self.tiles(dem.tiles_for_box(*box))
        lon, lat, z = dem.mosaic(box, files)
        step = 1.0 / (N - 1)
        self.assertTrue(np.all(np.diff(lat) > 0) and np.all(np.diff(lon) > 0))
        np.testing.assert_allclose(np.diff(lat), step, rtol=1e-9)
        np.testing.assert_allclose(np.diff(lon), step, rtol=1e-9)
        self.assertGreaterEqual(lat[0], box[0])
        self.assertLess(lat[0] - box[0], step)
        self.assertLessEqual(lat[-1], box[1])
        self.assertLess(box[1] - lat[-1], step)
        # every point, across the seams at 57 N and 8 W
        np.testing.assert_array_equal(z, form(lat[:, None], lon[None, :]))

    def test_box_edges_on_nodes_keep_those_nodes(self):
        # 0.7 N is node 3240 of the mosaic from 2 S, and -2 + 3240/1200 rounds above 0.7: the cut is by index
        box = (-1.2, 0.7, -0.5, -0.4)
        lon, lat, z = dem.mosaic(box, self.tiles(dem.tiles_for_box(*box)))
        self.assertAlmostEqual(lat[-1], 0.7, places=12)
        self.assertAlmostEqual(lat[0], -1.2, places=12)
        self.assertAlmostEqual(lon[0], -0.5, places=12)
        self.assertAlmostEqual(lon[-1], -0.4, places=12)
        self.assertEqual(z.shape, (2281, 121))
        # nodes whose (edge - corner) * 1200 rounds just below (north, east) or just above (south, west) the index
        box = (57.005, 57.01, -7.9975, -7.995)
        lon, lat, z = dem.mosaic(box, self.tiles([(57, -8)]))
        np.testing.assert_allclose([lat[0], lat[-1], lon[0], lon[-1]], box, atol=1e-12)
        self.assertEqual(z.shape, (7, 4))

    def test_voids_become_nan(self):
        box = (57.4, 57.6, -7.6, -7.4)            # holds the tile's middle point
        lon, lat, z = dem.mosaic(box, self.tiles([(57, -8)], void_at=(N // 2 + 1, N // 2 + 1)))
        self.assertEqual(int(np.sum(np.isnan(z))), 1)

    def test_mixed_resolutions_refused(self):
        files = self.tiles([(57, -8)])
        files.update({k: v for k, v in self.tiles([(57, -9)], n=3601).items()})
        with self.assertRaisesRegex(CaseError, "different resolutions"):
            dem.mosaic((57.1, 57.3, -8.3, -7.8), files)


class DownloadTests(TmpDir):
    def test_written_grid_reads_back_as_the_closed_form(self):
        server = Server({(57, -8): make_tile(57, -8)})
        out = os.path.join(self.dir, "dem.asc")
        box = (57.4, 57.7, -7.9, -7.5)       # the form changes sign in this box
        result = dem.download(box, out, self.cache, sea_floor="keep", opener=server, date="2026-11-15")
        r = terrain.read_esri_ascii(out, "terrain.file")
        self.assertTrue(r.is_geographic())
        np.testing.assert_allclose(r.values, form(r.y[:, None], r.x[None, :]), atol=0)
        self.assertEqual(r.values.shape, result.shape)
        self.assertEqual(result.lowest, float(form(r.y[0], r.x[0])))
        self.assertLess(result.lowest, 0.0)
        with open(result.sources, encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("registry.opendata.aws/terrain-tiles", text)
        self.assertIn("SRTM terrain data courtesy of the U.S. Geological Survey", text)
        self.assertIn("2026-11-15", text)
        self.assertIn("Tiles: N57W008", text)

    def test_sea_floor_zero(self):
        server = Server({(57, -8): make_tile(57, -8)})
        out = os.path.join(self.dir, "dem.asc")
        box = (57.4, 57.7, -7.9, -7.5)       # the form changes sign in this box
        result = dem.download(box, out, self.cache, sea_floor="zero", opener=server)
        r = terrain.read_esri_ascii(out, "terrain.file")
        expect = form(r.y[:, None], r.x[None, :])
        self.assertEqual(result.sea_points, int(np.sum(expect < 0)))
        self.assertGreater(result.sea_points, 0)
        self.assertLess(result.sea_points, r.values.size)
        np.testing.assert_array_equal(r.values, np.maximum(expect, 0.0))
        self.assertTrue(any("set to 0 m" in n for n in result.notes))

    def test_cache_and_requests(self):
        server = Server({(57, -8): make_tile(57, -8)})
        box = (57.4, 57.7, -7.9, -7.5)       # the form changes sign in this box
        dem.download(box, os.path.join(self.dir, "a.asc"), self.cache, opener=server)
        dem.download(box, os.path.join(self.dir, "b.asc"), self.cache, opener=server)
        self.assertEqual(server.requests, [dem.tile_url(57, -8)])     # the second download reads the cache

    def test_cut_download_leaves_nothing_cached(self):
        server = Server(broken={dem.tile_url(57, -8)})
        with self.assertRaisesRegex(CaseError, "N57W008.hgt.gz"):
            dem.download((57.1, 57.2, -7.9, -7.8), os.path.join(self.dir, "d.asc"), self.cache, opener=server)
        self.assertEqual(os.listdir(self.cache), [])
        self.assertFalse(os.path.exists(os.path.join(self.dir, "d.asc")))

    def test_failed_write_leaves_nothing_cached(self):
        server = Server({(57, -8): make_tile(57, -8)})
        with mock.patch("gust_case.dem.os.replace", side_effect=OSError("no space left on device")):
            with self.assertRaisesRegex(CaseError, "no space left on device"):
                dem.fetch_tile(57, -8, self.cache, opener=server)
        self.assertEqual(os.listdir(self.cache), [])

    def test_failed_grid_write_keeps_the_old_grid(self):
        server = Server({(57, -8): make_tile(57, -8)})
        out = os.path.join(self.dir, "dem.asc")
        with open(out, "w", encoding="utf-8") as old:
            old.write("the grid in use")
        dem.fetch_tile(57, -8, self.cache, opener=server)          # cached, so only the grid's write fails
        with mock.patch("gust_case.dem.os.replace", side_effect=PermissionError(13, "Permission denied")):
            with self.assertRaisesRegex(CaseError, "cannot write .*dem.asc: Permission denied"):
                dem.download((57.4, 57.7, -7.9, -7.5), out, self.cache, opener=server)
        with open(out, encoding="utf-8") as old:
            self.assertEqual(old.read(), "the grid in use")
        self.assertEqual(sorted(os.listdir(self.dir)), sorted(["dem.asc", os.path.basename(self.cache)]))

    def test_damaged_cached_tile_names_the_file(self):
        os.makedirs(self.cache, exist_ok=True)
        with open(os.path.join(self.cache, "N57W008.hgt"), "wb") as bad:
            bad.write(b"cut")
        with self.assertRaisesRegex(CaseError, "N57W008.hgt holds 3 bytes.*delete it to fetch the tile again"):
            dem.download((57.4, 57.7, -7.9, -7.5), os.path.join(self.dir, "d.asc"), self.cache, opener=Server())

    def test_cache_that_is_not_a_folder_is_refused(self):
        not_folder = os.path.join(self.dir, "file")
        open(not_folder, "w").close()
        with self.assertRaisesRegex(CaseError, "N57W008.hgt.gz"):
            dem.download((57.4, 57.7, -7.9, -7.5), os.path.join(self.dir, "d.asc"), os.path.join(not_folder, "sub"),
                         opener=Server({(57, -8): make_tile(57, -8)}))

    def test_missing_tile_names_the_url(self):
        with self.assertRaisesRegex(CaseError, "skadi/N57/N57W008.hgt.gz: 404"):
            dem.download((57.1, 57.2, -7.9, -7.8), os.path.join(self.dir, "d.asc"), self.cache, opener=Server())

    def test_refusals(self):
        with self.assertRaisesRegex(CaseError, "sea_floor = 'clamp'"):
            dem.download((57.1, 57.2, -7.9, -7.8), os.path.join(self.dir, "d.asc"), self.cache, sea_floor="clamp")
        with self.assertRaisesRegex(CaseError, "extension .asc"):
            dem.download((57.1, 57.2, -7.9, -7.8), os.path.join(self.dir, "d.tif"), self.cache)


def askervein_case(**domain):
    d = {"center_lat_lon": [57.18855, -7.37267], "half_width": [3000.0], "top": 1500.0}
    d.update(domain)
    raw = {"case": {"name": "t"}, "domain": d, "terrain": {"file": "dem.asc"},
           "mesh": {"cell_size": 50.0, "layers": 20, "first_layer": 4.0},
           "wind": {"speed": 8.0, "height": 10.0, "direction": 210.0}, "surface": {"roughness": 0.03}}
    return validate(raw)


class CaseBoxTests(TmpDir):
    """The download box holds every point the builder reads (terrain.prepare's cover check passes)."""

    def prepare(self, box, case, given):
        server = Server({(57, -8): make_tile(57, -8)})
        out = os.path.join(self.dir, "dem.asc")
        dem.download(box, out, self.cache, opener=server)
        frame = _resolve_domain(case, given, [])
        return terrain.prepare(out, "auto", "auto", frame, 200.0, 25.0, "lowest", 0.0, 0.0, "lowest")

    def test_box_covers_the_domain(self):
        for extra in ({}, {"half_width": [4000.0, 1500.0]}, {"origin": "southwest"}):
            case, given = askervein_case(**extra)
            box = dem.box_for_case(case, given)
            t = self.prepare(box, case, given)
            self.assertGreater(len(t.x), 1000)

    def test_a_smaller_box_is_refused_by_the_builder(self):
        # the cover check is what the box is for: a box half as wide leaves the domain's faces out
        case, given = askervein_case()
        s, n, w, e = dem.box_for_case(case, given)
        lat0, lon0 = 57.18855, -7.37267
        small = (lat0 - 0.5 * (lat0 - s), lat0 + 0.5 * (n - lat0), lon0 - 0.5 * (lon0 - w), lon0 + 0.5 * (e - lon0))
        with self.assertRaisesRegex(CaseError, "does not cover the domain on the west"):
            self.prepare(small, case, given)

    def test_margin_from_the_cell_counts(self):
        case, given = askervein_case()
        s, n, _, _ = dem.box_for_case(case, given)
        # 6000 m / 50 m = 120 cells; margin 4 cells = 200 m; reach 1.1 sqrt(2) 3200 m, plus one arc-second
        reach = 1.1 * np.sqrt(2.0) * 3200.0 / terrain.METERS_PER_DEGREE + 1.0 / 3600.0
        self.assertAlmostEqual(n - 57.18855, reach, places=12)
        self.assertAlmostEqual(57.18855 - s, reach, places=12)

    def test_pole_and_antimeridian_refused(self):
        case, given = askervein_case(center_lat_lon=[57.0, 179.99])
        with self.assertRaisesRegex(CaseError, "across the 180th meridian"):
            dem.box_for_case(case, given)
        case, given = askervein_case(center_lat_lon=[89.99, 10.0])
        with self.assertRaisesRegex(CaseError, "past a pole"):
            dem.box_for_case(case, given)

    def test_needs_center_lat_lon(self):
        raw = {"case": {"name": "t"}, "domain": {"center": [0.0, 0.0], "half_width": [3000.0], "top": 1500.0},
               "mesh": {"cell_size": 50.0, "layers": 20}, "wind": {"speed": 8.0, "height": 10.0, "direction": 0.0},
               "surface": {"roughness": 0.03}}
        case, given = validate(raw)
        with self.assertRaisesRegex(CaseError, "domain.center_lat_lon"):
            dem.box_for_case(case, given)


class CommandLineTests(TmpDir):
    def run_tool(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = fetch_terrain.main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_box_or_case_file(self):
        with self.assertRaises(SystemExit):
            self.run_tool([])

    def test_refusal_exits_2(self):
        path = os.path.join(self.dir, "case.yaml")
        with open(path, "w", encoding="utf-8") as out:
            out.write("case: {name: t}\ndomain: {center: [0, 0], half_width: 1000., top: 500.}\n"
                      "mesh: {cell_size: 50., layers: 10}\nwind: {speed: 5., height: 10., direction: 0.}\n"
                      "surface: {roughness: 0.03}\n")
        code, _, err = self.run_tool([path, "--cache", self.cache])
        self.assertEqual(code, 2)
        self.assertIn("domain.center_lat_lon", err)


if __name__ == "__main__":
    unittest.main()
