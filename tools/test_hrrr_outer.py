"""Tests of tools/hrrr_outer.py and tools/hrrr_outer_case.py on synthetic HRRR crops and checkpoints (numpy only).

    cd tools && python3 -m unittest -v test_hrrr_outer
"""
import math
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

import hrrr_outer as ho
import hrrr_outer_case as hc
from gust_case import geo
from gust_library import plotfile

LAT0, LON0 = 36.362, -97.405
R_HRRR, LAT_T, LON_C = 6371229.0, 38.5, -97.5     # HRRR's sphere and Lambert conformal projection
NL = math.sin(math.radians(LAT_T))
FL = math.cos(math.radians(LAT_T)) * math.tan(math.pi / 4 + math.radians(LAT_T) / 2) ** NL / NL


E_WGS84 = math.sqrt(0.00669437999014)


def conformal(lat):
    """Conformal latitude (radians) of a geodetic latitude (degrees): a sphere map of it is conformal on WGS84."""
    p = np.radians(lat)
    es = E_WGS84 * np.sin(p)
    return 2 * np.arctan(np.tan(np.pi / 4 + p / 2) * ((1 - es) / (1 + es)) ** (E_WGS84 / 2)) - np.pi / 2


def geodetic(chi):
    p = chi
    for _ in range(12):
        es = E_WGS84 * np.sin(p)
        p = 2 * np.arctan(np.tan(np.pi / 4 + chi / 2) * ((1 + es) / (1 - es)) ** (E_WGS84 / 2)) - np.pi / 2
    return np.degrees(p)


def lambert(lat, lon):
    """HRRR's Lambert conformal projection made conformal on the ellipsoid (through the conformal latitude), so that
    directions measured with gust_case.geo's UTM are exact; the cone constant, and so the grid rotation, is HRRR's."""
    rho = R_HRRR * FL / np.tan(np.pi / 4 + conformal(lat) / 2) ** NL
    t = NL * np.radians(np.asarray(lon) - LON_C)
    return rho * np.sin(t), -rho * np.cos(t)


def lambert_inv(x, y):
    rho = np.hypot(x, y)
    lat = geodetic(2 * np.arctan((R_HRRR * FL / rho) ** (1 / NL)) - np.pi / 2)
    return lat, LON_C + np.degrees(np.arctan2(x, -y)) / NL


H_LEV = np.array([10., 20., 60., 150., 300., 600., 1000., 1500., 2200., 3000., 3800.])


def crop(lat0=LAT0, lon0=LON0, n=10, fxx=1, cloud=None, wind=None):
    """A crop on HRRR's grid, 2n+1 columns a side: theta_v = 300 + 0.004 h exactly (p = 1e5, q = 0.01, T scaled by
    1/(1 + 0.61 q)), tke = 0.5 + 0.001 h, grid wind u = 5 + 0.002 h + 1e-4 x, v = 2 + 0.001 h (x east of the centre
    along the grid, m) or ``wind`` everywhere, ground rising 1 m per km east, cloud water ``cloud`` = (level, value)
    in the centre column."""
    x0, y0 = lambert(lat0, lon0)
    i = np.arange(-n, n + 1) * 3000.0
    xx, yy = np.meshgrid(x0 + i, y0 + i)
    lat, lon = lambert_inv(xx, yy)
    orog = 300.0 + 1e-3 * (xx - x0)
    h = H_LEV[:, None, None] * np.ones_like(xx)[None]
    one = np.ones_like(h)
    q = 0.01
    u, v = (5.0 + 0.002 * h + 1e-4 * (xx - x0)[None], 2.0 + 0.001 * h) if wind is None else \
        (wind[0] * one, wind[1] * one)
    d = {"lat": lat, "lon": lon, "sfc_orog": orog, "sfc_fsr": 0.05 * np.ones_like(xx),
         "nat_gh": orog[None] + h, "nat_pres": 1e5 * one, "nat_q": q * one, "nat_t": (300.0 + 0.004 * h) / (1 + 0.61 * q),
         "nat_u": u, "nat_v": v, "nat_w": 0 * one, "nat_tke": 0.5 + 0.001 * h,
         "nat_clwmr": 0 * one, "init": np.array("2023-05-10T19:00Z"), "fxx": np.array(fxx)}
    if cloud is not None:
        d["nat_clwmr"][cloud[0], n, n] = cloud[1]
    return d


# --------------------------------------------------------------------------------------------- a synthetic checkpoint
NX, NY, NZ, HALF_X, HALF_Y, TOP = 6, 4, 5, 6000.0, 4000.0, 3000.0
SPLIT = 3          # two boxes in x: cells 0..2 and 3..5
GHOST = 1
SENTINEL = 777.0


def _ground(x, y):
    return 10.0 + 2e-3 * x + 1e-3 * y


def _write_mf(folder, name, typ, ncomp, fill):
    """A VisMF MultiFab with two boxes (split in x), one ghost, little-endian doubles; fill(i, j, k) of global
    indices gives the values, SENTINEL outside the domain."""
    boxes = []
    for lo, hi in ((0, SPLIT - 1), (SPLIT, NX - 1)):
        hi_x = hi + (typ[0] if hi == NX - 1 else 0)
        boxes.append(((lo, 0, 0), (hi_x, NY - 1 + typ[1], NZ - 1 + typ[2])))
    lines = ["1", "1", str(ncomp), str(GHOST), f"({len(boxes)} 0"]
    lines += [f"(({b[0][0]},{b[0][1]},{b[0][2]}) ({b[1][0]},{b[1][1]},{b[1][2]}) ({typ[0]},{typ[1]},{typ[2]}))"
              for b in boxes]
    lines += [")", str(len(boxes))]
    data = b""
    offsets = []
    nmax = (NX + typ[0] - 1, NY + typ[1] - 1, NZ + typ[2] - 1)
    for blo, bhi in boxes:
        flo = [v - GHOST for v in blo]
        fhi = [v + GHOST for v in bhi]
        offsets.append(len(data))
        hdr = (f"FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))(({flo[0]},{flo[1]},{flo[2]}) "
               f"({fhi[0]},{fhi[1]},{fhi[2]}) ({typ[0]},{typ[1]},{typ[2]})) {ncomp}\n")
        k, j, i = np.meshgrid(np.arange(flo[2], fhi[2] + 1), np.arange(flo[1], fhi[1] + 1),
                              np.arange(flo[0], fhi[0] + 1), indexing="ij")
        inside = (i >= 0) & (j >= 0) & (k >= 0) & (i <= nmax[0]) & (j <= nmax[1]) & (k <= nmax[2])
        a = np.where(inside, fill(i, j, k), SENTINEL)
        data += hdr.encode() + np.ascontiguousarray(np.broadcast_to(a, (ncomp,) + a.shape), "<f8").tobytes()
    lines += [f"FabOnDisk: {name}_D_00000 {o}" for o in offsets]
    lines += ["", f"{len(boxes)},{ncomp}"] + ["0," * ncomp] * len(boxes) + ["", f"{len(boxes)},{ncomp}"] + \
        ["0," * ncomp] * len(boxes)
    with open(os.path.join(folder, name + "_H"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(folder, name + "_D_00000"), "wb") as f:
        f.write(data)


def node_z(i, j, k):
    dx, dy = 2 * HALF_X / NX, 2 * HALF_Y / NY
    zg = _ground(-HALF_X + dx * i, -HALF_Y + dy * j)
    with np.errstate(invalid="ignore"):      # k = -1 ghosts are NaN: a reader that used them would fail the tests
        return zg + (np.asarray(k) / NZ) ** 1.3 * (TOP - zg)


def make_checkpoint(path):
    lvl = os.path.join(path, "Level_0")
    os.makedirs(lvl)
    with open(os.path.join(path, "Header"), "w") as f:
        f.write("gust-checkpoint 2\nreal_bytes 8\nstep 0\n")
        f.write(f"n_cell {NX} {NY} {NZ}\n")
        f.write("prob_lo %s %s %s\n" % ((-HALF_X).hex(), (-HALF_Y).hex(), (0.0).hex()))
        f.write("prob_hi %s %s %s\n" % (HALF_X.hex(), HALF_Y.hex(), TOP.hex()))
    _write_mf(lvl, "z_phys_nd", (1, 1, 1), 1, node_z)
    for name, typ in (("u", (1, 0, 0)), ("v", (0, 1, 0)), ("w", (0, 0, 1)), ("theta", (0, 0, 0)),
                      ("tke", (0, 0, 0)), ("p", (0, 0, 0))):
        _write_mf(lvl, name, typ, 1, lambda i, j, k: -1.0 + 0 * i)


def read_back(chk, name, typ):
    """Independent reader (tools/gust_library/plotfile.py): the valid boxes, assembled to (nx, ny, nz)."""
    _, boxes = plotfile._multifab(os.path.join(chk, "Level_0", name), [0], 1)
    return plotfile._assemble(boxes, (NX + typ[0], NY + typ[1], NZ + typ[2]), [0])[0]


class Rotation(unittest.TestCase):
    """UT-16.2.1: HRRR grid wind -> UTM equals the closed form at two longitudes (east and west of the UTM central
    meridian, and either side of HRRR's 97.5 W). The smallest missing rotation (HRRR's alone at 97.4 W) is 1e-3 rad."""

    def test_closed_form(self):
        for lat, lon in ((36.362, -97.405), (40.0, -104.5), (33.0, -100.9)):
            zone = geo.utm_zone(lon)
            x, y = lambert(lat, lon)
            e0, n0, _, _ = geo.latlon_to_utm(lat, lon, zone=zone, south=False)
            for ug, vg, dxy in ((1.0, 0.0, (1.0, 0.0)), (0.0, 1.0, (0.0, 1.0))):
                la, lo = lambert_inv(x + dxy[0], y + dxy[1])          # one metre along the grid axis
                e1, n1, _, _ = geo.latlon_to_utm(la, lo, zone=zone, south=False)
                ex, ey = e1 - e0, n1 - n0
                s = math.hypot(ex, ey)
                u, v = ho.grid_to_utm(np.array(ug), np.array(vg), lat, lon, zone)
                self.assertAlmostEqual(float(u), ex / s, delta=2e-6, msg=(lat, lon, ug, vg))
                self.assertAlmostEqual(float(v), ey / s, delta=2e-6, msg=(lat, lon, ug, vg))

    def test_field_rotates_the_crop(self):
        """HrrrField turns the crop's grid wind: grid-x wind at 104.5 W points along grid x as seen in UTM."""
        lat, lon = 40.0, -104.5
        zone = geo.utm_zone(lon)
        x, y = lambert(lat, lon)
        e0, n0, _, _ = geo.latlon_to_utm(lat, lon, zone=zone, south=False)
        la, lo = lambert_inv(x + 1.0, y)
        e1, n1, _, _ = geo.latlon_to_utm(la, lo, zone=zone, south=False)
        ref = math.atan2(n1 - n0, e1 - e0)
        f = ho.HrrrField(crop(lat, lon, wind=(3.0, 0.0)), lat, lon)
        o, hh = np.zeros((1, 1)), 100.0 * np.ones((1, 1, 1))
        u, v = f.at("u", o, o, hh).item(), f.at("v", o, o, hh).item()
        self.assertAlmostEqual(math.hypot(u, v), 3.0, delta=1e-9)
        self.assertAlmostEqual(math.atan2(v, u), ref, delta=2e-6)

    def test_rotation_is_large_enough_to_matter(self):
        # at 104.5 W the two rotations add to about 4.7 degrees: a missing rotation is 0.08 of the speed
        u, v = ho.grid_to_utm(np.array(1.0), np.array(0.0), 40.0, -104.5, geo.utm_zone(-104.5))
        self.assertGreater(abs(math.degrees(math.atan2(float(v), float(u)))), 4.0)


class Field(unittest.TestCase):
    def test_closed_form_theta_and_tke(self):
        f = ho.HrrrField(crop(), LAT0, LON0)
        x = np.array([[-7000.0, 0.0, 12345.0]])
        y = np.array([[3000.0, -9000.0, 500.0]])
        h = np.array([12.0, 47.3, 812.0, 2999.0])[:, None, None] * np.ones((1, 1, 3))
        np.testing.assert_allclose(f.at("theta", x, y, h), 300.0 + 0.004 * h, atol=1e-9)
        np.testing.assert_allclose(f.at("tke", x, y, h), 0.5 + 0.001 * h, atol=1e-9)

    def test_log_law_below_first_level(self):
        """Below the first level (10 m) the speed is S ln(h/z0)/ln(10/z0) at fixed direction, w goes to 0 and theta
        is held; checked on the 2 m table rows (4 and 6 m) where no interpolation enters."""
        d = crop()
        d["nat_w"][0] = 0.3
        f = ho.HrrrField(d, LAT0, LON0)
        o = np.zeros((1, 1))
        for h in (4.0, 6.0):
            hh = h * np.ones((1, 1, 1))
            u, v = f.at("u", o, o, hh).item(), f.at("v", o, o, hh).item()
            u10, v10 = f.at("u", o, o, 10 * np.ones((1, 1, 1))).item(), f.at("v", o, o, 10 * np.ones((1, 1, 1))).item()
            self.assertAlmostEqual(math.hypot(u, v), math.hypot(u10, v10) * math.log(h / .05) / math.log(10 / .05),
                                   delta=1e-9)
            self.assertAlmostEqual(math.atan2(v, u), math.atan2(v10, u10), delta=1e-12)
            w10 = f.at("w", o, o, 10 * np.ones((1, 1, 1))).item()      # omega 0.3 Pa/s -> -0.3 / (rho g) m/s
            self.assertAlmostEqual(w10, -0.3 * 287.04 * 300.04 / (1e5 * 9.81), delta=2e-4)
            self.assertAlmostEqual(f.at("w", o, o, hh).item(), w10 * h / 10, delta=1e-12)
            self.assertAlmostEqual(f.at("theta", o, o, hh).item(), 300.04, delta=1e-9)

    def test_outside_crop_refused(self):
        f = ho.HrrrField(crop(n=4), LAT0, LON0)
        with self.assertRaisesRegex(ho.HrrrOuterError, "outside the crop"):
            f.at("u", np.array([[20000.0]]), np.array([[0.0]]), np.ones((1, 1, 1)))


class ForecastHour(unittest.TestCase):
    """Mutant: an f00 crop accepted."""

    def test_f00_refused(self):
        with self.assertRaisesRegex(ho.HrrrOuterError, "fxx = 0: the forecast hour must be 1 or later"):
            ho.HrrrField(crop(fxx=0), LAT0, LON0)

    def test_fetch_refuses_f00_before_download(self):
        with self.assertRaisesRegex(ho.HrrrOuterError, "--fxx = 0"):
            ho.fetch("2023-05-10T19:00Z", 0, LAT0, LON0, 30.0, "unused", "unused.npz")

    def test_f01_accepted(self):
        self.assertIn("f01 (valid 2023-05-10T20:00Z)", ho.HrrrField(crop(fxx=1), LAT0, LON0).label)


class FetchRegex(unittest.TestCase):
    """fetch downloads every field the tool reads: each crop key's name in HRRR's .idx inventory matches the regex
    (cloud water is CLMR in wrfnat; cfgrib names it clwmr). Mutant: CLMR left out of NAT_REGEX."""
    NAT = {"nat_u": "UGRD", "nat_v": "VGRD", "nat_w": "VVEL", "nat_gh": "HGT", "nat_t": "TMP", "nat_pres": "PRES",
           "nat_q": "SPFH", "nat_tke": "TKE", "nat_clwmr": "CLMR"}
    SFC = {"sfc_orog": ":HGT:surface:", "sfc_fsr": ":SFCR:surface:", "sfc_ishf": ":SHTFL:surface:"}

    def test_every_crop_field_matched(self):
        import re
        for key in ho.CROP_KEYS:
            if key.startswith("nat_"):
                for lev in (1, 30):
                    self.assertRegex(f"12:345:d=2023051019:{self.NAT[key]}:{lev} hybrid level:1 hour fcst:",
                                     ho.NAT_REGEX, key)
                self.assertIsNone(re.search(ho.NAT_REGEX, f":{self.NAT[key]}:31 hybrid level:"), key)
            elif key.startswith("sfc_"):
                self.assertRegex(f"7:890:d=2023051019{self.SFC[key]}1 hour fcst:", ho.SFC_REGEX, key)

    def test_skin_temperature_matched(self):
        """The skin temperature (hrrr.surface = surface_temperature) is TMP at the surface, not only at 2 m.
        Mutant: :(HGT|TMP|PRES):surface back to :(HGT|PRES):surface."""
        self.assertRegex("7:890:d=2023051019:TMP:surface:1 hour fcst:", ho.SFC_REGEX)
        self.assertRegex("7:890:d=2023051019:TMP:2 m above ground:1 hour fcst:", ho.SFC_REGEX)


class Cloud(unittest.TestCase):
    """UT-16.2.3: an hour with cloud in the domain is refused by name."""

    def test_cloud_refused_by_name(self):
        f = ho.HrrrField(crop(cloud=(6, 2e-5)), LAT0, LON0)      # level 6: 1000 m
        with self.assertRaisesRegex(ho.HrrrOuterError,
                                    r"HRRR 2023-05-10T19:00Z f01 \(valid 2023-05-10T20:00Z\): cloud water 0.02 g/kg "
                                    r"at 1000 m .*this hour is refused"):
            f.check_cloud(10000.0, 3000.0)

    def test_cloud_above_top_or_below_limit_accepted(self):
        ho.HrrrField(crop(cloud=(10, 2e-5)), LAT0, LON0).check_cloud(10000.0, 3000.0)   # 3800 m
        ho.HrrrField(crop(cloud=(6, 5e-6)), LAT0, LON0).check_cloud(10000.0, 3000.0)

    def test_to_chk_refuses_cloudy_hour(self):
        tmp = tempfile.mkdtemp(prefix="ut_hrrr_outer_")
        try:
            make_checkpoint(os.path.join(tmp, "chk"))
            with self.assertRaisesRegex(ho.HrrrOuterError, "cloud water"):
                ho.write_checkpoint(ho.HrrrField(crop(cloud=(4, 1e-4)), LAT0, LON0), os.path.join(tmp, "chk"))
            # refused before anything is written
            np.testing.assert_array_equal(read_back(os.path.join(tmp, "chk"), "theta", (0, 0, 0)), -1.0)
        finally:
            shutil.rmtree(tmp)


class Checkpoint(unittest.TestCase):
    """UT-16.2.2: a checkpoint written and read back holds the given fields (two boxes, ghosts, sloped ground)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ut_hrrr_outer_")
        self.chk = os.path.join(self.tmp, "chk")
        make_checkpoint(self.chk)
        self.field = ho.HrrrField(crop(), LAT0, LON0)
        self.stats = ho.write_checkpoint(self.field, self.chk)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _agl_cells(self):
        i, j, k = np.meshgrid(np.arange(NX + 1), np.arange(NY + 1), np.arange(NZ + 1), indexing="ij")
        z = node_z(i, j, k)
        zc = sum(z[a:a + NX, b:b + NY, c:c + NZ] for a in (0, 1) for b in (0, 1) for c in (0, 1)) / 8
        zg = sum(z[a:a + NX, b:b + NY, 0] for a in (0, 1) for b in (0, 1)) / 4
        return zc - zg[:, :, None]

    def test_theta_and_tke(self):
        h = self._agl_cells()
        np.testing.assert_allclose(read_back(self.chk, "theta", (0, 0, 0)), 300.0 + 0.004 * h, atol=1e-9)
        np.testing.assert_allclose(read_back(self.chk, "tke", (0, 0, 0)), 0.5 + 0.001 * h, atol=1e-9)

    def test_w_zero_on_ground_and_p_kept(self):
        w = read_back(self.chk, "w", (0, 0, 1))
        np.testing.assert_array_equal(w[:, :, 0], 0.0)
        np.testing.assert_array_equal(read_back(self.chk, "p", (0, 0, 0)), -1.0)

    def test_u_v_w_on_their_faces(self):
        """u on x faces, v on y faces, w on z faces, each at its face's height above the ground (the mean of the
        face's four nodes minus the mean of their ground nodes); the field varies in x and height, so a face
        taken for a cell or a height for another fails."""
        dx, dy = 2 * HALF_X / NX, 2 * HALF_Y / NY
        i, j, k = np.meshgrid(np.arange(NX + 1), np.arange(NY + 1), np.arange(NZ + 1), indexing="ij")
        z = node_z(i, j, k)
        h = z - z[:, :, :1]
        xn, yn = -HALF_X + dx * np.arange(NX + 1), -HALF_Y + dy * np.arange(NY + 1)
        xc, yc = 0.5 * (xn[1:] + xn[:-1]), 0.5 * (yn[1:] + yn[:-1])
        faces = {"u": ((1, 0, 0), xn, yc, 0.25 * (h[:, :-1, :-1] + h[:, 1:, :-1] + h[:, :-1, 1:] + h[:, 1:, 1:])),
                 "v": ((0, 1, 0), xc, yn, 0.25 * (h[:-1, :, :-1] + h[1:, :, :-1] + h[:-1, :, 1:] + h[1:, :, 1:])),
                 "w": ((0, 0, 1), xc, yc, 0.25 * (h[:-1, :-1] + h[1:, :-1] + h[:-1, 1:] + h[1:, 1:]))}
        for name, (typ, xs, ys, hf) in faces.items():
            X, Y = np.meshgrid(xs, ys, indexing="ij")
            want = np.transpose(self.field.at(name, X.T, Y.T, np.transpose(hf, (2, 1, 0))), (2, 1, 0))
            if name == "w":
                want[:, :, 0] = 0.0
            np.testing.assert_allclose(read_back(self.chk, name, typ), want, atol=1e-12, err_msg=name)
        self.assertGreater(np.ptp(read_back(self.chk, "u", (1, 0, 0))[:, 1, 2]), 1.0)     # varies along x

    def test_ghosts_outside_domain_kept(self):
        for flo, fhi, a in ho.chk_fabs(self.chk, "theta"):
            self.assertEqual(float(a[0, 0, 1, 1]), SENTINEL)           # k = -1: below the ground
            self.assertEqual(float(a[0, -1, 1, 1]), SENTINEL)          # k = NZ: above the top

    def test_balance_of_a_uniform_wind(self):
        """u = 1, v = 2, w = 0: each face's flux is the speed times the face's area, 2 HALF (TOP - mean ground) on
        the linear ground (west 8000 x 3002, east 8000 x 2978, south 12000 x 2994, north 12000 x 2986 m2)."""
        for name, val in (("u", 1.0), ("v", 2.0), ("w", 0.0)):
            for _, _, a in ho.chk_fabs(self.chk, name, mode="r+"):
                a[:] = val
                a.flush()
        r = ho.balance(self.chk)
        self.assertAlmostEqual(r["west"], 8000.0 * 3002.0, delta=1e-6)
        self.assertAlmostEqual(r["east"], -8000.0 * 2978.0, delta=1e-6)
        self.assertAlmostEqual(r["south"], 2.0 * 12000.0 * 2994.0, delta=1e-6)
        self.assertAlmostEqual(r["north"], -2.0 * 12000.0 * 2986.0, delta=1e-6)
        self.assertEqual(r["top_outflow"], 0.0)
        with self.assertRaisesRegex(ho.HrrrOuterError, "--cells"):
            ho.balance(self.chk, 2)

    def test_refused_top_leaves_the_checkpoint(self):
        """A crop whose highest level lies below the domain top is refused before any field is written."""
        d = crop()
        d["nat_gh"] = d["sfc_orog"][None] + 0.7 * (d["nat_gh"] - d["sfc_orog"][None])     # top level 2660 m
        before = read_back(self.chk, "u", (1, 0, 0)).copy()
        with self.assertRaisesRegex(ho.HrrrOuterError, r"the domain top \(3000 m\) lies above HRRR's highest"):
            ho.write_checkpoint(ho.HrrrField(d, LAT0, LON0), self.chk)
        np.testing.assert_array_equal(read_back(self.chk, "u", (1, 0, 0)), before)


class Score(unittest.TestCase):
    """The score table on a small plotfile (tools/test_gust_library.py's writer, three uneven boxes, sloped
    ground): a nest equal to HRRR scores 0, a nest 1 m/s faster everywhere scores bias = rms = 1 m/s."""

    def _run(self, extra_speed, verdict="steady", margin=0.0, w_spike=None):
        import test_gust_library as tgl
        f = ho.HrrrField(crop(), LAT0, LON0)
        zn, xc, yc, agl = tgl.mesh(sloped=True)
        X, Y = np.meshgrid(xc, yc, indexing="ij")
        T = lambda a: np.transpose(a, (2, 1, 0))
        g = {n: T(f.at(n, X.T, Y.T, T(agl))) for n in ("u", "v", "w", "tke")}
        s = np.hypot(g["u"], g["v"])
        fac = (s + extra_speed) / s
        tmp = tempfile.mkdtemp(prefix="ut_hrrr_outer_")
        self.addCleanup(shutil.rmtree, tmp)
        plt = os.path.join(tmp, "plt00042")
        w = g["w"].copy()
        if w_spike is not None:
            w[0] = w_spike                 # the first column of cells in x, inside any margin
        tgl.write_plotfile(plt, g["u"] * fac, g["v"] * fac, w, g["tke"], zn)
        log = os.path.join(tmp, "log.txt")
        with open(log, "w") as fh:
            fh.write(f"gust: {verdict} after 42 iterations (42 over the run)\n"
                     "gust: stage account: 42 iterations, 1 s, 0 AU\n")
        return ho.score(f, plt, margin=margin, log=log, low=30.0), g

    def test_equal_fields_score_zero(self):
        sc, g = self._run(0.0)
        for key in ("bias", "rms", "low_speed_bias", "low_vector_rms", "low_dir_diff"):
            self.assertAlmostEqual(sc[key], 0.0, delta=1e-9, msg=key)
        self.assertAlmostEqual(sc["k_ratio"], 1.0, delta=1e-12)
        self.assertAlmostEqual(sc["w_rms"], sc["w_rms_hrrr"], delta=1e-12)
        self.assertEqual((sc["iterations"], sc["steady"], sc["interior_cells"]), (42, True, 7 * 5 * 5))

    def test_faster_nest(self):
        sc, g = self._run(1.0)
        for key in ("bias", "rms", "low_speed_bias", "low_vector_rms"):
            self.assertAlmostEqual(sc[key], 1.0, delta=1e-9, msg=key)
        # the direction of the mean vector moves 2e-4 degrees: +1 m/s weighs the slower, more turned cells more
        self.assertAlmostEqual(sc["low_dir_diff"], 0.0, delta=1e-3)
        self.assertEqual(len(sc["levels"]), 5)
        self.assertAlmostEqual(sc["levels"][0]["agl"], 2.0, delta=1e-12)


    def test_not_steady_and_diverged(self):
        self.assertIs(self._run(0.0, "not steady")[0]["steady"], False)
        self.assertIs(self._run(0.0, "diverged")[0]["steady"], False)

    def test_w_max_over_the_interior(self):
        """w_max, like w_rms, covers the cells inside the margin: a spike in the first column of cells is left
        out with a margin of one cell and counted without one."""
        import test_gust_library as tgl
        _, xc, _, _ = tgl.mesh(sloped=True)
        dx = xc[1] - xc[0]
        self.assertLess(self._run(0.0, margin=dx, w_spike=50.0)[0]["w_max"], 1.0)
        self.assertEqual(self._run(0.0, margin=0.0, w_spike=50.0)[0]["w_max"], 50.0)


class Forcing(unittest.TestCase):
    """The geostrophic wind and the surface flux of the hrrr section (variants C, D, F) against closed forms."""

    def test_geostrophic_closed_form(self):
        """p = 1e5 + a X + b Y at every level (X, Y: UTM metres from the centre): G = (-b, a) / (rho f), rho from
        p / (R_d T_v) with T_v = 300 + 0.004 h, as a mean over 200..1000 m every 100 m."""
        d = crop()
        f0 = ho.HrrrField(d, LAT0, LON0)
        a, b = 1.0e-3, -2.0e-3                                   # Pa/m
        nk = d["nat_gh"].shape[0]
        d["nat_pres"] = np.broadcast_to(1e5 + a * f0.X + b * f0.Y, (nk,) + f0.X.shape).copy()
        fld = ho.HrrrField(d, LAT0, LON0)
        f = 2 * ho.OMEGA * math.sin(math.radians(LAT0))
        want = np.mean([(-b, a) for h in range(200, 1001, 100)], axis=0)
        tv = np.mean([300.0 + 0.004 * h for h in range(200, 1001, 100)])
        rho_inv = np.mean([ho.RD * (300.0 + 0.004 * h) / 1e5 for h in range(200, 1001, 100)])
        ug, vg = fld.geostrophic(10000.0, LAT0)
        np.testing.assert_allclose((ug, vg), want * rho_inv / f, rtol=2e-4)
        self.assertGreater(tv, 300.0)

    def test_geostrophic_refused_in_steep_terrain(self):
        d = crop()
        d["sfc_orog"] = d["sfc_orog"].copy()
        d["sfc_orog"][10, 12] += 400.0                           # a hill 6 km east of the centre
        d["nat_gh"] = d["nat_gh"].copy()
        d["nat_gh"][:, 10, 12] += 400.0
        with self.assertRaisesRegex(ho.HrrrOuterError, "geostrophic wind: 200 m above the box-mean ground"):
            ho.HrrrField(d, LAT0, LON0).geostrophic(10000.0, LAT0)

    def test_surface_flux_closed_form(self):
        d = crop()
        one = np.ones_like(d["sfc_orog"])
        d.update({"sfc_sp": 1e5 * one, "sfc_t2m": 300.0 * one, "sfc_sh2": 0.01 * one, "sfc_ishf": 200.0 * one,
                  "sfc_slhtf": 100.0 * one})
        rho = 1e5 / (287.04 * 300.0 * 1.0061)
        want = 200.0 / (rho * 1004.0) * 1.0061 + 0.61 * 300.0 * 100.0 / (rho * 2.5e6)
        self.assertAlmostEqual(ho.HrrrField(d, LAT0, LON0).surface_flux(10000.0), want, delta=1e-12)

    def test_surface_theta_closed_form(self):
        """theta_v,s = T_skin (1e5 / p_s)^0.2857 (1 + 0.61 q_2m), column by column. Mutant: q_2m left out."""
        d = crop()
        one = np.ones_like(d["sfc_orog"])
        tskin = 300.0 + 0.001 * (np.arange(one.size).reshape(one.shape))
        d.update({"sfc_sp": 9e4 * one, "sfc_sh2": 0.01 * one, "sfc_t": tskin})
        got = ho.HrrrField(d, LAT0, LON0).surface_theta()
        want = (tskin * (1e5 / 9e4) ** 0.2857 * 1.0061).ravel()
        np.testing.assert_allclose(got, want, rtol=1e-14)

    def test_surface_theta_refusals(self):
        d = crop()
        one = np.ones_like(d["sfc_orog"])
        d.update({"sfc_sp": 1e5 * one, "sfc_sh2": 0.01 * one})
        with self.assertRaisesRegex(ho.HrrrOuterError, r"no skin temperature \(sfc_t\)"):
            ho.HrrrField(d, LAT0, LON0).surface_theta()
        d["sfc_t"] = 400.0 * one
        with tempfile.TemporaryDirectory(prefix="ut_hrrr_tsurf_") as tmp:
            with self.assertRaisesRegex(ho.HrrrOuterError, r"must be finite and in \[150, 350\] K"):
                ho.HrrrField(d, LAT0, LON0).write_surface_theta(os.path.join(tmp, "t.xyz"), 10000.0)
            self.assertFalse(os.path.exists(os.path.join(tmp, "t.xyz")))


class Top(unittest.TestCase):
    def test_top_refusals(self):
        fld = ho.HrrrField(crop(), LAT0, LON0)
        with self.assertRaisesRegex(ho.HrrrOuterError, r"--top \(3700 m\) must be > 0 and at most 3600 m"):
            ho.theta_profile(fld, 10000.0, 3700.0)
        with self.assertRaisesRegex(ho.HrrrOuterError, r"--half-width \(0.0\) must be > 0"):
            ho.theta_profile(fld, 0.0, 3000.0)
        self.assertEqual(ho.theta_profile(fld, 10000.0, 3000.0)[0][-1], 3000.0)

    def test_cloud_check_needs_the_crop_to_reach(self):
        """The crop is 60 km a side; the cloud check of a raw half-width of 28 km reaches 31 km, outside it."""
        with self.assertRaisesRegex(ho.HrrrOuterError, "lies outside the crop"):
            ho.HrrrField(crop(), LAT0, LON0).check_cloud(28000.0, 3000.0)

    def test_command_line_refusal_exits_1(self):
        tmp = tempfile.mkdtemp(prefix="ut_hrrr_outer_")
        self.addCleanup(shutil.rmtree, tmp)
        path = os.path.join(tmp, "c.npz")
        np.savez(path, **crop())
        import contextlib
        import io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = ho.main(["profile", path, "--lat", str(LAT0), "--lon", str(LON0), "--half-width", "5000",
                          "--top", "4000"])
        self.assertEqual(rc, 1)
        self.assertIn("hrrr_outer.py: --top (4000 m)", err.getvalue())


class Meshes(unittest.TestCase):
    def test_awaken(self):
        m = hc.plan_meshes(12500.0, 125.0, 4, 2000.0)
        self.assertEqual(m, {"nest": (12500.0, 125.0), "parent": (26000.0, 500.0), "hrrr_raw": (28000.0, 1000.0)})

    def test_faces_align_in_a_sweep(self):
        for cell in (50.0, 100.0, 125.0, 200.0):
            for n8 in range(4, 40):
                half = 4.0 * cell * n8
                for ratio in (2, 4, 8, 16):
                    if (4 * n8) % ratio:
                        continue
                    m = hc.plan_meshes(half, cell, ratio, 1500.0)
                    (hp, cp), (hr, cr) = m["parent"], m["hrrr_raw"]
                    self.assertGreaterEqual(hp, 2 * half)
                    self.assertGreaterEqual(hr, hp + 1500.0)
                    for a, b, c in ((hp, half, cp), (hr, hp, cr)):
                        self.assertAlmostEqual((a - b) / c, round((a - b) / c), places=9)
                    for a, c in ((hp, cp), (hr, cr)):
                        self.assertEqual(round(2 * a / c) % 8, 0)

    def test_refusals(self):
        with self.assertRaisesRegex(ho.HrrrOuterError, "multiple of 8"):
            hc.plan_meshes(12250.0, 125.0, 4, 2000.0)
        with self.assertRaisesRegex(ho.HrrrOuterError, "hrrr.parent_ratio"):
            hc.plan_meshes(500.0 * 5, 125.0, 16, 2000.0)


class Chain(unittest.TestCase):
    """build_chain on a small flat case: three decks and one run.sh, the scheme written where it belongs."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ut_hrrr_outer_")
        self.addCleanup(shutil.rmtree, self.tmp)
        np.savez(os.path.join(self.tmp, "hour.npz"), **crop())
        with open(os.path.join(self.tmp, "case.yaml"), "w") as fh:
            fh.write("case: {name: c}\n"
                     f"domain: {{center_lat_lon: [{LAT0}, {LON0}], half_width: 1000., top: 3000.}}\n"
                     "mesh: {cell_size: 125., layers: 20, first_layer: 10.}\n"
                     "surface: {roughness: 0.05}\n"
                     "gust_inputs: {steady.max_iterations: 77}\n"
                     "hrrr: {file: hour.npz, theta: frozen}\n")

    def _deck(self, level):
        out = {}
        for line in open(os.path.join(self.tmp, "c", level, "inputs")):
            line = line.split("#", 1)[0]
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
        return out

    def test_decks(self):
        info = hc.build_chain(os.path.join(self.tmp, "case.yaml"))
        self.assertEqual(info["meshes"], {"nest": (1000.0, 125.0), "parent": (2000.0, 500.0),
                                          "hrrr_raw": (4000.0, 1000.0)})
        for level in ("hrrr_raw", "parent"):
            d = self._deck(level)
            stage = [k for k in d if k.startswith(hc.STAGE_NAMESPACES) or k in hc.STAGE_KEYS]
            self.assertEqual(stage, [], level)
        p, n = self._deck("parent"), self._deck("nest")
        self.assertEqual(p["bc.nest_parent"], "../hrrr_raw/chk00000")
        self.assertEqual(float(p["terrain.nest_blend_width"]), 2000.0)
        for face in ("xlo", "xhi", "ylo", "yhi"):
            self.assertEqual((p["bc." + face], n["bc." + face]), ("nest", "nest"))
        want = {"steady.stage": "simplec", "bc.nest_parent": "../parent/chk00000", "bc.nest_ratio": "4",
                "init.first_guess": "nest", "forcing.coriolis": "1", "simplec.theta": "frozen",
                "steady.max_iterations": "77"}
        self.assertEqual({k: n.get(k) for k in want}, want)
        self.assertEqual(float(n["turb.max_geom_lscale"]), 100.0)
        self.assertAlmostEqual(float(n["forcing.latitude"]), LAT0)
        self.assertAlmostEqual(float(n["forcing.theta0"]), 300.04, delta=1e-9)    # theta_v held below 10 m
        self.assertAlmostEqual(float(n["bc.zhi_theta_grad"]), 0.004, delta=1e-6)
        self.assertNotIn("forcing.gravity", n)
        run = open(os.path.join(self.tmp, "c", "run.sh")).read()
        order = [run.index(t) for t in ("\nstep hrrr_raw", "\" to-chk ", "\nstep parent", "\nstep nest", "\" score ")]
        self.assertEqual(order, sorted(order))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "c", "nest", "run.sh")))

    def test_rerun_clears_the_parent_checkpoint(self):
        """run.sh's clean-up line removes the parent checkpoint a previous run left (a directory) and succeeds."""
        hc.build_chain(os.path.join(self.tmp, "case.yaml"))
        case = os.path.join(self.tmp, "c")
        clean = [ln for ln in open(os.path.join(case, "run.sh")).read().splitlines() if ln.startswith("rm ")]
        self.assertEqual(len(clean), 1)
        os.makedirs(os.path.join(case, "parent", "chk00000", "Level_0"))
        open(os.path.join(case, "parent", "chk00000", "Header"), "w").close()
        r = subprocess.run(["bash", "-e", "-c", clean[0]], cwd=case, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(os.path.exists(os.path.join(case, "parent", "chk00000")))

    def test_surface_temperature_map(self):
        """hrrr.surface = surface_temperature: the nest deck names the map, and the map holds HRRR's columns within
        the nest's half-width plus 6 km in the case frame, with theta_v,s of each. Mutant: the map written in the
        crop's grid order but with the columns of the whole crop."""
        d = crop()
        one = np.ones_like(d["sfc_orog"])
        d.update({"sfc_sp": 1e5 * one, "sfc_sh2": 0.01 * one, "sfc_t": 305.0 + 1e-4 * d["sfc_orog"]})
        np.savez(os.path.join(self.tmp, "hour.npz"), **d)
        yaml = open(os.path.join(self.tmp, "case.yaml")).read().replace("theta: frozen", "theta: frozen, "
                                                                        "surface: surface_temperature")
        open(os.path.join(self.tmp, "case.yaml"), "w").write(yaml)
        info = hc.build_chain(os.path.join(self.tmp, "case.yaml"))
        n = self._deck("nest")
        self.assertEqual((n["bc.surface_stability"], n["bc.surface_temperature_file"]),
                         ("surface_temperature", hc.SURFACE_THETA_FILE))
        for level in ("hrrr_raw", "parent"):
            self.assertNotIn("bc.surface_temperature_file", self._deck(level))
        pts = np.loadtxt(os.path.join(self.tmp, "c", "nest", hc.SURFACE_THETA_FILE))
        fld = ho.HrrrField(d, LAT0, LON0)
        m = fld.box(7000.0)
        self.assertEqual(len(pts), int(m.sum()))
        self.assertEqual(info["surface_temperature"][0], len(pts))
        np.testing.assert_allclose(pts[:, 0], fld.X.ravel()[m], atol=0.005)
        np.testing.assert_allclose(pts[:, 1], fld.Y.ravel()[m], atol=0.005)
        want = (305.0 + 1e-4 * d["sfc_orog"].ravel()[m]) * 1.0061
        np.testing.assert_allclose(pts[:, 2], want, atol=5e-4)

    def test_strip_stage_keys(self):
        path = os.path.join(self.tmp, "deck")
        with open(path, "w") as fh:
            fh.write("# a deck\nmesh.n_cell = 8 8 8\nturb.tke_floor = 1e-6\nsteady.stage = simplec   # x\n"
                     "steady.cfl_max = 4\nforcing.theta0 = 300\nconvection.hybrid_height = 9\nsimplec.w_vnn_cap = 2\n"
                     "io.turb_plot_file = turb\nio.plot_file = plt\nabl.theta = 300\nabl.drive = speed\n"
                     "bc.xlo = nest\n\n")
        removed = hc.strip_stage_keys(path)
        self.assertEqual(len(removed), 9)
        self.assertEqual(open(path).read(), "# a deck\nmesh.n_cell = 8 8 8\nio.plot_file = plt\nbc.xlo = nest\n\n")


class CaseSection(unittest.TestCase):
    def test_unknown_key_and_ranges(self):
        with tempfile.NamedTemporaryFile(suffix=".npz") as f:
            base = os.path.dirname(f.name)
            for sec, msg in (({"file": f.name, "ratio": 4}, "hrrr.ratio"),
                             ({"file": f.name, "parent_ratio": 3}, "hrrr.parent_ratio"),
                             ({"file": f.name, "theta": "capped"}, "hrrr.theta"),
                             ({"file": f.name, "surface": "skin"}, r"hrrr.surface \(skin\) must be neutral, "
                              "heat_flux or surface_temperature"),
                             ({"file": f.name, "raw_margin": 0}, "hrrr.raw_margin"),
                             ({"file": f.name, "cloud_mixing_ratio": -1}, "hrrr.cloud_mixing_ratio"),
                             ({"file": f.name, "cloud_mixing_ratio": "1e-5"}, r"hrrr.cloud_mixing_ratio \('1e-5'\) "
                              "must be a number"),
                             ({"file": f.name, "parent_ratio": 4.0}, r"hrrr.parent_ratio \(4.0\)"),
                             ({}, "hrrr.file")):
                with self.assertRaisesRegex(ho.HrrrOuterError, msg):
                    hc.read_hrrr_section({"hrrr": sec}, base)

    def test_scheme_refusals(self):
        raw = {"domain": {"center_lat_lon": [LAT0, LON0], "half_width": 12500., "top": 2000.},
               "wind": {"speed": 5}, "mesh": {"hierarchy_levels": 2}, "gust_inputs": {"bc.nest_ratio": 2}}
        with self.assertRaises(ho.HrrrOuterError) as cm:
            hc.check_case(raw)
        for name in ("domain.top (2000.0) must be at least 3000", "wind: set by", "mesh.hierarchy_levels",
                     "gust_inputs.bc.nest_ratio"):
            self.assertIn(name, str(cm.exception))
        raw = {"domain": {"center_lat_lon": [LAT0, LON0], "half_width": 12500., "top": 4000., "utm_zone": 15},
               "output": {"prefix": "a_"}}
        with self.assertRaises(ho.HrrrOuterError) as cm:
            hc.check_case(raw)
        for name in ("domain.top (4000.0) must be at most 3600", "domain.utm_zone (15): HRRR is placed in the "
                     "zone of the centre (14)", "output.prefix"):
            self.assertIn(name, str(cm.exception))
        hc.check_case({"domain": {"center_lat_lon": [LAT0, LON0], "half_width": 12500., "top": 3000.,
                                  "utm_zone": 14}})


if __name__ == "__main__":
    unittest.main()
