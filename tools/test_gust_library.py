"""Tests of tools/gust_library and tools/make_gust_library.py (CTest gust_library_selftest).

The plotfiles here are written by hand in AMReX's format (one or several boxes,
double or single precision), with fields chosen so that the answers are closed
forms: a field linear in the height above a sloping ground, uniform winds from
known directions.

    cd tools && python3 -m unittest -v test_gust_library
"""
import contextlib
import io
import json
import math
import os
import struct
import sys
import tempfile
import unittest

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import make_gust_library  # noqa: E402
from gust_library import Library, LibraryError, SampleError, build, index, plotfile, sample  # noqa: E402

NX, NY = 7, 5
DX = 50.0
LO = (-175.0, -125.0)                      # prob_lo x, y: the centre is (0, 0)
ETA = (0.0, 4.0, 10.0, 20.0, 40.0, 80.0)   # node heights above the ground, the same in every column
SLOPE_X, SLOPE_Y, BASE = 0.1, -0.05, 30.0  # a plane ground: z_g = BASE + SLOPE_X x + SLOPE_Y y
BOXES = [((0, 0, 0), (2, 4, 4)), ((3, 0, 0), (6, 1, 4)), ((3, 2, 0), (6, 4, 4))]   # uneven split


def ground_at(x, y, sloped=True):
    return BASE + SLOPE_X * x + SLOPE_Y * y if sloped else BASE + 0.0 * x


def _fab(fh, lo, hi, comps, single, nodal=False, ghost=0):
    """One FAB over lo..hi (ghost cells included): header line, then the components, x fastest.

    comps: arrays over the domain, padded by ``ghost`` points on every side."""
    order = "4 3 2 1" if single else "8 7 6 5 4 3 2 1"
    fmt = "(8, (64 11 52 0 1 12 0 1023))" if not single else "(8, (32 8 23 0 1 9 0 127))"
    typ = "(1,1,1)" if nodal else "(0,0,0)"
    head = "FAB (%s,(%d, (%s)))((%d,%d,%d) (%d,%d,%d) %s) %d\n" % (fmt, 4 if single else 8, order, *lo, *hi, typ, len(comps))
    off = fh.tell()
    fh.write(head.encode())
    for a in comps:
        sub = a[lo[0] + ghost:hi[0] + 1 + ghost, lo[1] + ghost:hi[1] + 1 + ghost, lo[2] + ghost:hi[2] + 1 + ghost]
        fh.write(struct.pack("<%d%s" % (sub.size, "f" if single else "d"), *sub.transpose(2, 1, 0).ravel()))
    return off


def _multifab(prefix, boxes, comps, single, nodal=False, ghost=0):
    """A VisMF MultiFab; with ghost > 0 each FAB holds that many ghost points of -999 around its box."""
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    data = os.path.basename(prefix) + "_D_00000"
    offs = []
    padded = [np.pad(a, ghost, constant_values=-999.0) for a in comps]
    with open(os.path.join(os.path.dirname(prefix), data), "wb") as fh:
        for lo, hi in boxes:
            if nodal:
                hi = tuple(h + 1 for h in hi)
            offs.append(_fab(fh, tuple(v - ghost for v in lo), tuple(v + ghost for v in hi), padded, single, nodal,
                             ghost))
    with open(prefix + "_H", "w") as fh:
        fh.write("1\n0\n%d\n0\n(%d 0\n" % (len(comps), len(boxes)))
        for lo, hi in boxes:
            hi_w = tuple(h + 1 for h in hi) if nodal else hi
            fh.write("((%d,%d,%d) (%d,%d,%d) %s)\n" % (*lo, *hi_w, "(1,1,1)" if nodal else "(0,0,0)"))
        fh.write(")\n%d\n" % len(boxes))
        for off in offs:
            fh.write("FabOnDisk: %s %d\n" % (data, off))


def write_plotfile(path, u, v, w, k, z_nodes, boxes=BOXES, single=False):
    """A GUST-like plotfile: u, v, w, p, theta, tke, z and Level_0/z_phys_nd."""
    nx, ny, nz = u.shape
    z_cc = 0.125 * sum(z_nodes[i:i + nx, j:j + ny, m:m + nz] for i in (0, 1) for j in (0, 1) for m in (0, 1))
    names = ["u", "v", "w", "p", "theta", "tke", "z"]
    comps = [u, v, w, np.zeros_like(u), np.full_like(u, 300.0), k, z_cc]
    os.makedirs(os.path.join(path, "Level_0"), exist_ok=True)
    with open(os.path.join(path, "Header"), "w") as fh:
        fh.write("HyperCLaw-V1.1\n%d\n%s\n3\n0\n0\n" % (len(names), "\n".join(names)))
        fh.write("%g %g 0\n%g %g 1000\n\n" % (LO[0], LO[1], LO[0] + nx * DX, LO[1] + ny * DX))
        fh.write("((0,0,0) (%d,%d,%d) (0,0,0))\n0\n%g %g %g\n" % (nx - 1, ny - 1, nz - 1, DX, DX, 1000.0 / nz))
    _multifab(os.path.join(path, "Level_0", "Cell"), boxes, comps, single)
    _multifab(os.path.join(path, "Level_0", "z_phys_nd"), boxes, [z_nodes], single, nodal=True, ghost=1)


def mesh(sloped=True, nx=NX, ny=NY):
    """(z_nodes, x and y of the column centres, cell-centre heights above the ground)."""
    xn = LO[0] + np.arange(nx + 1) * DX
    yn = LO[1] + np.arange(ny + 1) * DX
    X, Y = np.meshgrid(xn, yn, indexing="ij")
    z_nodes = ground_at(X, Y, sloped)[:, :, None] + np.asarray(ETA)[None, None, :]
    xc, yc = 0.5 * (xn[1:] + xn[:-1]), 0.5 * (yn[1:] + yn[:-1])
    eta_c = 0.5 * (np.asarray(ETA[1:]) + np.asarray(ETA[:-1]))
    agl = np.broadcast_to(eta_c, (nx, ny, len(eta_c))).copy()
    return z_nodes, xc, yc, agl


def uniform_run(folder, direction_from, speed, k=0.5, status="steady", sloped=False, single=False, nan=False):
    """A run folder whose wind is from ``direction_from`` everywhere, ``speed`` x h / 10 m at h above the ground
    (``speed`` at 10 m, 5 x ``speed`` at 50 m: linear in h, so sampled exactly); k uniform."""
    z_nodes, xc, yc, agl = mesh(sloped)
    towards = math.radians(270.0 - direction_from)
    shape = agl.shape
    u = speed * math.cos(towards) * agl / 10.0
    v = speed * math.sin(towards) * agl / 10.0
    if nan:
        u[1, 1, 1] = float("nan")
    os.makedirs(folder, exist_ok=True)
    write_plotfile(os.path.join(folder, "plt00010"), u, v, np.zeros(shape), np.full(shape, k), z_nodes, single=single)
    with open(os.path.join(folder, "inputs"), "w") as fh:
        fh.write("bc.surface_z0 = 0.03\nbc.inflow_z0 = 0.05\n")
    write_log(folder, status)


def write_log(folder, status="steady", plt="plt00010"):
    with open(os.path.join(folder, "log.txt"), "w") as fh:
        fh.write("  Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r\n")
        fh.write("    10  1.600e+01  lusgs        1.0e-05  2.0e-06  3.0e-06  4.0e-06  5.0e-10  2.0e-09  9.870e-05\n")
        if status:
            fh.write("gust: %s after 10 iterations (10 over the run)\n" % status)
        if plt:
            fh.write("gust: wrote %s\ngust: wrote turb00010\n" % plt)


def write_case(case_dir, directions):
    os.makedirs(case_dir, exist_ok=True)
    with open(os.path.join(case_dir, "case_resolved.yaml"), "w") as fh:
        yaml.safe_dump({"case": {"name": "t"}, "wind": {"directions": list(directions), "speed": 8.0}}, fh)


def quiet_main(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = make_gust_library.main(argv)
    return rc, out.getvalue(), err.getvalue()


class PlotfileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_uneven_boxes_and_precision(self):
        """Three uneven boxes, double and single precision: every value lands at its (i, j, k)."""
        z_nodes, xc, yc, agl = mesh()
        idx = np.indices(agl.shape)
        u = 100.0 * idx[0] + 10.0 * idx[1] + idx[2]
        for single in (False, True):
            p = os.path.join(self.dir, "plt_%d" % single)
            write_plotfile(p, u, -u, 0 * u, 0 * u + 1, z_nodes, single=single)
            cells, hdr = plotfile.read_cells(p, ["u", "v"])
            np.testing.assert_array_equal(cells["u"], u)
            np.testing.assert_array_equal(cells["v"], -u)
            np.testing.assert_allclose(plotfile.read_nodes(p), z_nodes, rtol=1e-6 if single else 0, atol=0)

    def test_missing_field_and_box_gap(self):
        z_nodes, xc, yc, agl = mesh()
        p = os.path.join(self.dir, "plt")
        write_plotfile(p, agl, agl, agl, agl, z_nodes, boxes=BOXES[:2])
        with self.assertRaisesRegex(ValueError, "no field speed"):
            plotfile.read_cells(p, ["speed"])
        with self.assertRaisesRegex(ValueError, "do not cover the domain"):
            plotfile.read_cells(p, ["u"])


class SampleTests(unittest.TestCase):
    def test_linear_in_height_above_sloping_ground(self):
        """UT: a field a + b h + c x + d y, h the height above a sloping ground, sampled at fixed heights."""
        z_nodes, xc, yc, agl = mesh(sloped=True)
        X, Y = np.meshgrid(xc, yc, indexing="ij")
        a, b, c, d = 2.0, 0.07, 0.003, -0.002
        z_cc = 0.125 * sum(z_nodes[i:i + NX, j:j + NY, m:m + agl.shape[2]] for i in (0, 1) for j in (0, 1) for m in (0, 1))
        h_cc = sample.heights_above_ground(z_cc, z_nodes)
        np.testing.assert_allclose(h_cc, agl, atol=1e-12)
        np.testing.assert_allclose(sample.ground(z_nodes), ground_at(X, Y), atol=1e-12)
        f = a + b * h_cc + (c * X + d * Y)[:, :, None]
        got = sample.at_heights(f, h_cc, [2.0, 10.0, 33.0, 60.0])
        for n, h in enumerate([2.0, 10.0, 33.0, 60.0]):
            np.testing.assert_allclose(got[n], a + b * h + c * X + d * Y, atol=1e-12)

    def test_heights_refused(self):
        z_nodes, xc, yc, agl = mesh()
        for bad, msg in ((1.0, "below the first cell centre"), (61.0, "above the last cell centre"),
                         (0.0, "not a finite height"), (float("nan"), "not a finite height")):
            with self.assertRaisesRegex(SampleError, "--heights: .*" + msg):
                sample.at_heights(agl, agl, [bad])
        with self.assertRaisesRegex(SampleError, "given twice"):
            sample.at_heights(agl, agl, [10.0, 10.0])

    def test_point(self):
        """Bilinear at a point; exact for a + c x + d y; refused outside the outer centres."""
        _, xc, yc, _ = mesh()
        X, Y = np.meshgrid(xc, yc, indexing="ij")
        f = 1.0 + 0.01 * X - 0.02 * Y
        self.assertAlmostEqual(float(sample.at_point(xc, yc, f, 12.0, -31.0)), 1.0 + 0.12 + 0.62, places=12)
        self.assertAlmostEqual(float(sample.at_point(xc, yc, f, xc[-1], yc[-1])), f[-1, -1], places=12)
        with self.assertRaisesRegex(SampleError, "--reference-point: .*outside the column centres"):
            sample.at_point(xc, yc, f, xc[-1] + 1.0, 0.0)


class LibraryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.case = os.path.join(self.tmp.name, "case")

    def tearDown(self):
        self.tmp.cleanup()

    def sweep(self, runs):
        """runs: [(direction, speed, k, status)]."""
        write_case(self.case, [r[0] for r in runs])
        from gust_case.wind import direction_label
        for d, s, k, st in runs:
            uniform_run(os.path.join(self.case, "dir_" + direction_label(d)), d, s, k, st)

    def test_blend_between_two_sectors(self):
        """UT: uniform winds from 0 (6 m/s, k 0.36) and 10 degrees (10 m/s, k 2): the answer at 5 degrees."""
        self.sweep([(0.0, 6.0, 0.36, "steady"), (10.0, 10.0, 2.0, "steady")])
        idx = build(self.case, [10.0])
        self.assertEqual([s["direction"] for s in idx["sectors"]], [0.0, 10.0])
        self.assertAlmostEqual(idx["sectors"][0]["reference_speed"], 6.0, places=12)
        self.assertAlmostEqual(idx["sectors"][1]["reference_direction"], 10.0, places=9)
        w = Library(self.case).wind(5.0, 8.0, 10.0)
        towards = math.radians(265.0)
        np.testing.assert_allclose(w["u"], 8.0 * math.cos(towards), atol=1e-12)
        np.testing.assert_allclose(w["v"], 8.0 * math.sin(towards), atol=1e-12)
        k = 0.5 * (0.36 / 36.0 + 2.0 / 100.0) * 8.0 ** 2      # k per unit reference speed squared, blended
        np.testing.assert_allclose(w["k"], k, rtol=1e-12)
        self.assertEqual(w["weights"], [(0.0, 0.5), (10.0, 0.5)])
        w = Library(self.case).wind(0.5, 4.0, 10.0)          # half a degree from a sector: blended
        self.assertEqual(w["weights"], [(0.0, 0.95), (10.0, 0.05)])
        w = Library(self.case).wind(2.5, 4.0, 10.0)
        self.assertEqual(w["weights"], [(0.0, 0.75), (10.0, 0.25)])
        np.testing.assert_allclose(np.hypot(w["u"], w["v"]), 4.0, rtol=1e-12)

    def test_convention_west_wind(self):
        """UT: a run with u > 0, v = 0 is wind from 270 (the west); querying 270 gives u = speed."""
        self.sweep([(270.0, 7.0, 0.5, "steady"), (280.0, 7.0, 0.5, "steady")])
        idx = build(self.case, [10.0])
        self.assertAlmostEqual(idx["sectors"][0]["reference_direction"], 270.0, places=9)
        w = Library(self.case).wind(270.0, 5.0, 10.0)
        self.assertEqual(w["weights"], [(270.0, 1.0)])
        np.testing.assert_allclose(w["u"], 5.0, atol=1e-12)
        np.testing.assert_allclose(w["v"], 0.0, atol=1e-12)

    def test_wrap_past_north(self):
        """UT: sectors 350, 0 and 10; 355 blends 350 and 0, 5 blends 0 and 10, 357.5 weights 0.25 and 0.75."""
        self.sweep([(350.0, 8.0, 0.5, "steady"), (0.0, 8.0, 0.5, "steady"), (10.0, 8.0, 0.5, "steady")])
        build(self.case, [10.0])
        lib = Library(self.case)
        self.assertEqual(lib.wind(355.0, 8.0, 10.0)["weights"], [(350.0, 0.5), (0.0, 0.5)])
        self.assertEqual(lib.wind(-5.0, 8.0, 10.0)["weights"], [(350.0, 0.5), (0.0, 0.5)])
        self.assertEqual(lib.wind(5.0, 8.0, 10.0)["weights"], [(0.0, 0.5), (10.0, 0.5)])
        w = lib.wind(357.5, 8.0, 10.0)
        self.assertEqual(w["weights"], [(350.0, 0.25), (0.0, 0.75)])
        # The blend of the unit vectors with weights 0.25 and 0.75, scaled to 8 m/s (it points 0.0024 degrees
        # from 357.5: a vector blend is on the bisector only midway).
        e = 0.25 * np.array([math.cos(math.radians(-80.0)), math.sin(math.radians(-80.0))]) + 0.75 * np.array([0.0, -1.0])
        e *= 8.0 / np.hypot(*e)
        np.testing.assert_allclose(w["u"], e[0], atol=1e-12)
        np.testing.assert_allclose(w["v"], e[1], atol=1e-12)
        with self.assertRaisesRegex(LibraryError, "180 degrees apart"):
            lib.wind(180.0, 8.0, 10.0)

    def test_unconverged_sector_refused_by_name(self):
        """UT: a sector without GUST's steady verdict is refused when the blend uses it, and only then."""
        self.sweep([(260.0, 8.0, 0.5, "steady"), (270.0, 8.0, 0.5, "not steady"), (280.0, 8.0, 0.5, "steady")])
        idx = build(self.case, [10.0])
        self.assertEqual([s["converged"] for s in idx["sectors"]], [True, False, True])
        lib = Library(self.case)
        with self.assertRaisesRegex(LibraryError, "dir_270 .*did not reach GUST's steady verdict \\(not steady\\)"):
            lib.wind(265.0, 8.0, 10.0)
        with self.assertRaisesRegex(LibraryError, "dir_270"):
            lib.wind(270.0, 8.0, 10.0)
        lib.wind(260.0, 8.0, 10.0)

    def test_index_fields(self):
        """The index records the log's numbers; the npz holds the fields over the reference speed."""
        self.sweep([(90.0, 4.0, 0.16, "steady"), (100.0, 4.0, 0.16, "steady")])
        idx = build(self.case, [10.0, 50.0], reference_point=(25.0, -25.0))
        with open(os.path.join(self.case, "wind_library", "index.json")) as fh:
            self.assertEqual(json.load(fh), idx)
        s = idx["sectors"][0]
        self.assertEqual((s["iterations"], s["last_r"], s["status"], s["surface_z0"], s["inflow_z0"]),
                         (10, 9.87e-05, "steady", 0.03, 0.05))
        self.assertEqual(s["plotfile"], os.path.join("dir_090", "plt00010"))
        self.assertEqual(idx["reference_point"], [25.0, -25.0])
        with np.load(os.path.join(self.case, "wind_library", s["file"])) as data:
            self.assertEqual(data["u"].shape, (2, NX, NY))
            np.testing.assert_allclose(data["u"][0], -1.0, atol=1e-12)
            np.testing.assert_allclose(data["u"][1], -5.0, atol=1e-12)
            np.testing.assert_allclose(data["k"], 0.01, rtol=1e-12)

    def test_subset_and_refusals(self):
        self.sweep([(0.0, 8.0, 0.5, "steady"), (20.0, 8.0, 0.5, "steady"), (40.0, 8.0, 0.5, "steady")])
        idx = build(self.case, [10.0], directions=[40.0, 0.0])
        self.assertEqual([s["direction"] for s in idx["sectors"]], [0.0, 40.0])
        with self.assertRaisesRegex(LibraryError, "--directions: 30 is not a direction of the sweep"):
            build(self.case, [10.0], directions=[30.0])
        with self.assertRaisesRegex(LibraryError, "--heights: 1 m is below the first cell centre"):
            build(self.case, [1.0])
        with self.assertRaisesRegex(LibraryError, "--reference-point: .*outside"):
            build(self.case, [10.0], reference_point=(1e4, 0.0))
        lib = Library(self.case)
        with self.assertRaisesRegex(LibraryError, "--height: 20.0 m is not a height of the library"):
            lib.wind(10.0, 8.0, 20.0)
        with self.assertRaisesRegex(LibraryError, "--speed: "):
            lib.wind(10.0, 0.0, 10.0)
        with self.assertRaisesRegex(LibraryError, "--direction: "):
            lib.wind(float("nan"), 8.0, 10.0)
        for bad in (True, "8", None):
            with self.assertRaisesRegex(LibraryError, "--speed: "):
                lib.wind(10.0, bad, 10.0)
        with self.assertRaisesRegex(LibraryError, "--direction: "):
            lib.wind(np.True_, 8.0, 10.0)
        with self.assertRaisesRegex(LibraryError, "--height: "):
            lib.wind(10.0, 8.0, "10")
        # numpy scalars, single precision included, are numbers
        w = lib.wind(np.float32(20.0), np.float32(8.0), np.float32(10.0))   # the library holds 0 and 40
        self.assertEqual(w["weights"], [(0.0, 0.5), (40.0, 0.5)])
        np.testing.assert_allclose(np.hypot(w["u"], w["v"]), 8.0, rtol=1e-12)
        lib.wind(np.int64(10), 8, np.float64(10.0))
        s32 = np.float32(8.3)                     # k scales by the speed squared in double precision
        np.testing.assert_allclose(lib.wind(0.0, s32, 10.0)["k"], 0.5 / 64.0 * float(s32) ** 2, rtol=1e-14)

    def test_rebuild_replaces_whole_library(self):
        """A refused build leaves the previous library; a build of fewer sectors leaves no stale sector file."""
        self.sweep([(0.0, 8.0, 0.5, "steady"), (20.0, 8.0, 0.5, "steady")])
        idx = build(self.case, [10.0])
        lib_dir = os.path.join(self.case, "wind_library")
        before = sorted(os.listdir(lib_dir))
        os.remove(os.path.join(self.case, "dir_020", "log.txt"))
        with self.assertRaisesRegex(LibraryError, "dir_020"):
            build(self.case, [10.0, 40.0])
        self.assertEqual(sorted(os.listdir(lib_dir)), before)
        self.assertEqual(Library(self.case).index, idx)
        self.assertFalse(os.path.exists(lib_dir + ".tmp"))   # refused before anything is written
        build(self.case, [10.0], directions=[0.0])
        self.assertEqual(sorted(os.listdir(lib_dir)), ["grid.npz", "index.json", "sector_000.npz"])
        self.assertFalse(os.path.exists(lib_dir + ".old"))

    def test_library_files_checked(self):
        """A library of another format version, or a sector of the wrong shape, is refused."""
        self.sweep([(0.0, 8.0, 0.5, "steady"), (20.0, 8.0, 0.5, "steady")])
        build(self.case, [10.0])
        lib_dir = os.path.join(self.case, "wind_library")
        with np.load(os.path.join(lib_dir, "sector_020.npz")) as d:
            np.savez(os.path.join(lib_dir, "sector_020.npz"), **{f: d[f][:, :-1] for f in ("u", "v", "w", "k")})
        lib = Library(self.case)
        lib.wind(0.0, 8.0, 10.0)
        with self.assertRaisesRegex(LibraryError, "sector_020.npz: arrays of shape"):
            lib.wind(10.0, 8.0, 10.0)
        p = os.path.join(lib_dir, "index.json")
        with open(p) as fh:
            data = json.load(fh)
        data["format_version"] = 2
        with open(p, "w") as fh:
            json.dump(data, fh)
        with self.assertRaisesRegex(LibraryError, "format_version 2"):
            Library(self.case)
        with open(p, "w") as fh:
            fh.write("{")
        with self.assertRaisesRegex(LibraryError, "not a library index"):
            Library(self.case)

    def test_runs_refused(self):
        write_case(self.case, [0.0, 10.0])
        uniform_run(os.path.join(self.case, "dir_000"), 0.0, 8.0)
        with self.assertRaisesRegex(LibraryError, "dir_010.*log.txt: not found"):
            build(self.case, [10.0])
        uniform_run(os.path.join(self.case, "dir_010"), 10.0, 8.0)
        write_log(os.path.join(self.case, "dir_010"), status=None)
        with self.assertRaisesRegex(LibraryError, "no verdict line"):
            build(self.case, [10.0])
        write_log(os.path.join(self.case, "dir_010"), plt=None)
        with self.assertRaisesRegex(LibraryError, "no plotfile plt<step>"):
            build(self.case, [10.0])
        uniform_run(os.path.join(self.case, "dir_010"), 10.0, 8.0, nan=True)
        with self.assertRaisesRegex(LibraryError, "not finite"):
            build(self.case, [10.0])
        uniform_run(os.path.join(self.case, "dir_010"), 10.0, 8.0, sloped=True)
        with self.assertRaisesRegex(LibraryError, "the mesh differs"):
            build(self.case, [10.0])
        with open(os.path.join(self.case, "case_resolved.yaml"), "w") as fh:
            yaml.safe_dump({"case": {"name": "t"}, "wind": {"direction": 10.0}}, fh)
        with self.assertRaisesRegex(LibraryError, "no wind.directions"):
            build(self.case, [10.0])

    def test_log_reading(self):
        """The verdict, r before it, and the plotfile written after it: notes after r, Newton rows, a level-1
        row, a plotfile before the verdict, a hierarchy plotfile and the deck's io.plot_file."""
        folder = os.path.join(self.tmp.name, "run")
        uniform_run(folder, 0.0, 8.0)
        os.rename(os.path.join(folder, "plt00010"), os.path.join(folder, "flow00010"))
        with open(os.path.join(folder, "inputs"), "a") as fh:
            fh.write("io.plot_file = flow   # the deck's prefix\n")
        with open(os.path.join(folder, "log.txt"), "w") as fh:
            fh.write("gust: wrote flow00005\n"
                     "     5  1.600e+01  lusgs  1e-5 1e-5 1e-5 1e-5 1e-9 1e-9  4.000e-03  confidence 0.95\n"
                     "     6        inf  newton 1e-5 1e-5 1e-5 1e-5 1e-9 1e-9  3.000e-04\n"
                     "     7          -  newton 1e-5 1e-5 1e-5 1e-5 1e-9 1e-9  2.000e-04  Newton: 2 inner\n"
                     "  level 1:     7  1.0e+01  lusgs 1e-5 1e-5 1e-5 1e-5 1e-9 1e-9  9.000e-01\n"
                     "gust: steady after 7 iterations (7 over the run)\n"
                     "     8  1.600e+01  lusgs  1e-5 1e-5 1e-5 1e-5 1e-9 1e-9  1.000e-07\n"
                     "gust: wrote flow00010\ngust: wrote flow_level1_00010\ngust: wrote plt00011\n")
        status, it, r, plt = index.read_log(folder)
        self.assertEqual((status, it, r, os.path.basename(plt)), ("steady", 7, 2.0e-4, "flow00010"))
        self.assertIsNone(index.table_r("  level 1:     7  1.0e+01  lusgs 1 1 1 1 1 1  9.0e-01"))
        self.assertIsNone(index.table_r("  Iter  CFL  Line solver  u-res v-res w-res tke-res theta-res p-res r"))
        with open(os.path.join(folder, "log.txt"), "w") as fh:
            fh.write("gust: wrote flow00010\ngust: steady after 7 iterations (7 over the run)\n")
        with self.assertRaisesRegex(LibraryError, "no plotfile flow<step> written after the verdict"):
            index.read_log(folder)
        with open(os.path.join(folder, "log.txt"), "w") as fh:    # a plotfile after an earlier verdict only
            fh.write("gust: not steady after 5 iterations (5 over the run)\ngust: wrote flow00010\n"
                     "gust: steady after 7 iterations (12 over the run)\n")
        with self.assertRaisesRegex(LibraryError, "no plotfile flow<step> written after the verdict"):
            index.read_log(folder)
        for line, msg in (("gust: unsteady run ended after 100 steps\n", "an unsteady run"),
                          ("gust: joint run: outer + 1 nest\n", "two-way joint run")):
            with open(os.path.join(folder, "log.txt"), "w") as fh:
                fh.write(line + "gust: steady after 7 iterations (7 over the run)\ngust: wrote flow00010\n")
            with self.assertRaisesRegex(LibraryError, msg):
                index.read_log(folder)

    def test_plotfile_refusals(self):
        """One column of cells, a VisMF header of another version, and a header cut short are refused."""
        z_nodes, xc, yc, agl = mesh(nx=1)
        folder = os.path.join(self.tmp.name, "run")
        os.makedirs(folder)
        write_plotfile(os.path.join(folder, "plt00010"), agl, agl, agl, agl, z_nodes, boxes=[((0, 0, 0), (0, 4, 4))])
        write_log(folder)
        with self.assertRaisesRegex(LibraryError, "inputs: not found"):
            index.sector(folder, [10.0], None)
        open(os.path.join(folder, "inputs"), "w").close()
        with self.assertRaisesRegex(LibraryError, "1 x 5 x 5 cells; the library needs at least 2"):
            index.sector(folder, [10.0], None)
        uniform_run(folder, 0.0, 8.0)
        cell_h = os.path.join(folder, "plt00010", "Level_0", "Cell_H")
        with open(cell_h) as fh:
            text = fh.read()
        with open(cell_h, "w") as fh:
            fh.write("2" + text[1:])
        with self.assertRaisesRegex(ValueError, "VisMF header version 2"):
            index.sector(folder, [10.0], None)
        with open(cell_h, "w") as fh:
            fh.write("\n".join(text.split("\n")[:6]))
        with self.assertRaisesRegex(ValueError, "Cell_H: not a complete VisMF MultiFab"):
            index.sector(folder, [10.0], None)

    def test_command_line(self):
        self.sweep([(0.0, 8.0, 0.5, "steady"), (10.0, 8.0, 0.5, "steady")])
        rc, out, err = quiet_main(["build", self.case, "--heights", "10", "40"])
        self.assertEqual(rc, 0, err)
        self.assertIn("2 sectors, heights 10, 40 m", out)
        self.assertIn(os.path.join(self.case, "wind_library"), out)
        npz = os.path.join(self.tmp.name, "w")
        rc, out, err = quiet_main(["query", self.case, "--direction", "5", "--speed", "6", "--height", "40",
                                   "--point", "10", "20", "--out", npz])
        self.assertEqual(rc, 0, err)
        self.assertIn("speed 24 m/s from 5.0000 degrees", out)      # 6 m/s at 10 m, 4 x at 40 m
        self.assertIn("wrote %s.npz" % npz, out)                     # numpy adds the suffix; the message says so
        with np.load(npz + ".npz") as data:
            self.assertEqual(data["u"].shape, (NX, NY))
        rc, out, err = quiet_main(["query", self.case, "--direction", "5", "--speed", "6", "--height", "40",
                                   "--out", os.path.join(self.tmp.name, "no_such_dir", "w.npz")])
        self.assertEqual(rc, 2)
        self.assertIn("no_such_dir", err)
        rc, out, err = quiet_main(["query", self.case, "--direction", "5", "--speed", "6", "--height", "11"])
        self.assertEqual(rc, 2)
        self.assertIn("not a height of the library", err)
        rc, out, err = quiet_main(["build", self.case, "--heights", "0.5"])
        self.assertEqual(rc, 2)
        self.assertIn("--heights", err)


if __name__ == "__main__":
    unittest.main()
