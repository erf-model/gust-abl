#!/usr/bin/env python3
"""HRRR as the outer level of a GUST nest: crop, checkpoint, decks and score.

    hrrr_outer.py fetch --cycle INIT_TIME --fxx N --lat LAT --lon LON --half-km H --save DIR --out CROP.npz
    hrrr_outer.py to-chk CROP.npz CHK --lat LAT --lon LON [--cloud-mixing-ratio Q]
    hrrr_outer.py balance CHK [--cells M]
    hrrr_outer.py profile CROP.npz --lat LAT --lon LON --half-width M [--top M]
    hrrr_outer.py score CROP.npz PLOTFILE --lat LAT --lon LON [--log LOG] [--margin M] [--json OUT]
    hrrr_outer.py case CASE.yaml [--output-dir DIR] [--overwrite]

The chain (``case`` writes its decks and a run.sh that runs it):

1. ``hrrr_raw``: a GUST run that only writes its first guess, on a coarse mesh larger than the parent. Its
   checkpoint is the skeleton: ``to-chk`` overwrites u, v, w, theta and tke with HRRR's at every point.
2. ``parent``: ``init.first_guess = nest`` from the raw checkpoint, which applies GUST's mass-consistent adjustment;
   again first guess only.
3. ``nest``: the steady run, with the parent's field on its four lateral faces (``bc.nest_ratio``).

The scheme is fixed: a forecast hour of 1 or later (the f00 analysis carries almost no vertical motion), a domain
top of at least 3000 m (HRRR's horizontal convergence is balanced near 3 km, not below), the adjusted parent, the
virtual potential temperature, Coriolis on, the length-scale cap at 100 m, and an hour with cloud water anywhere in
the raw domain below its top refused by name (GUST is dry).

HRRR's fields. The crop holds HRRR's native (hybrid) levels [k, j, i] and surface fields [j, i] (``fetch`` writes
it). Winds are rotated from HRRR's grid (Lambert conformal, true latitude 38.5 N, central longitude 97.5 W) to east
and north, then to the UTM grid of the domain centre (meridian convergence c = atan(tan(lon - lon_cm) sin lat)).
theta = T (1e5 / p)^0.2857 (1 + 0.61 q) is the virtual potential temperature; w = -omega Rd Tv / (p g). Each column
is put on a table of heights above HRRR's own ground, 0 to 3600 m every 2 m; below HRRR's first level u and v follow
the neutral log law with HRRR's z0 at fixed direction, w goes linearly to 0, theta and k keep the first level's
value. A point of GUST takes the table at its height above GUST's ground, interpolated bilinearly in HRRR's grid
indices (the map from indices to UTM inverted by Newton's method). Needs numpy; ``fetch`` needs herbie, xarray and
cfgrib.
"""
import argparse
import datetime as dt
import glob
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hrrr_regime  # noqa: E402
from gust_case import geo  # noqa: E402

G, RD, KAPPA_P, OMEGA = 9.81, 287.04, 0.2857, 7.292e-5
CP, LV = 1004.0, 2.5e6
DH, HTOP = 2.0, 3600.0
# the height table: finer near the ground, where the log law curves, then every DH m to HTOP
HG = np.concatenate([[0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5], np.arange(DH, HTOP + 0.5 * DH, DH)])
MIN_TOP = 3000.0
CLOUD_MIXING_RATIO = 1e-5          # kg/kg; the classifier's cloud threshold (hrrr_regime.THRESHOLDS)
FIELDS = ("u", "v", "w", "theta", "tke")
CROP_KEYS = ("lat", "lon", "nat_u", "nat_v", "nat_w", "nat_gh", "nat_t", "nat_pres", "nat_q", "nat_tke",
             "nat_clwmr", "sfc_orog", "sfc_fsr", "init", "fxx")


class HrrrOuterError(Exception):
    """A refused input or hour; the message names it."""


# ---------------------------------------------------------------------------------------------------- winds

def earth_to_utm(u_e, v_n, lat, lon, zone):
    """East and north wind components to the UTM grid of ``zone`` (x along grid east, y along grid north)."""
    lon_cm = 6.0 * zone - 183.0
    c = np.arctan(np.tan(np.radians(np.asarray(lon) - lon_cm)) * np.sin(np.radians(lat)))
    return u_e * np.cos(c) - v_n * np.sin(c), u_e * np.sin(c) + v_n * np.cos(c)


def grid_to_utm(ug, vg, lat, lon, zone):
    """HRRR grid-relative wind to the UTM grid of ``zone``."""
    u_e, v_n = hrrr_regime.rotate_to_earth(ug, vg, lon)
    return earth_to_utm(u_e, v_n, lat, lon, zone)


def check_site(lat0, lon0):
    if not (math.isfinite(lat0) and 21.0 <= lat0 <= 53.0):
        raise HrrrOuterError(f"--lat ({lat0}) must lie in HRRR's domain, 21 to 53")
    if not (math.isfinite(lon0) and -135.0 <= lon0 <= -60.0):
        raise HrrrOuterError(f"--lon ({lon0}) must lie in HRRR's domain, -135 to -60")


def check_fxx(fxx, what="--fxx"):
    if int(fxx) != fxx or fxx < 1:
        raise HrrrOuterError(f"{what} = {fxx}: the forecast hour must be 1 or later; the f00 analysis carries almost "
                             "no vertical motion")


# ---------------------------------------------------------------------------------------------------- the field

class HrrrField:
    """HRRR's u, v, w, theta (virtual) and tke of a crop as functions of (x, y, height above ground).

    x, y: UTM (zone of the centre) minus the centre's easting and northing, the frame of a gust_case case with
    domain.origin = center."""

    def __init__(self, crop, lat0, lon0):
        check_site(lat0, lon0)
        d = dict(np.load(crop)) if isinstance(crop, str) else crop
        missing = [k for k in CROP_KEYS if k not in d]
        if missing:
            raise HrrrOuterError(f"crop: no {', '.join(missing)} (write it with 'hrrr_outer.py fetch')")
        self.fxx = int(d["fxx"])
        self.init = str(d["init"])
        check_fxx(self.fxx, f"crop {self.init}: fxx")
        valid = hrrr_regime.parse_time(self.init, "init") + dt.timedelta(hours=self.fxx)
        self.valid = hrrr_regime.fmt_time(valid)
        self.label = f"HRRR {self.init} f{self.fxx:02d} (valid {self.valid})"
        lat = np.asarray(d["lat"], float)
        lon = np.asarray(d["lon"], float)
        lon = np.where(lon > 180.0, lon - 360.0, lon)
        self.shape = lat.shape
        self.zone = geo.utm_zone(lon0)
        e, n, _, _ = geo.latlon_to_utm(lat, lon, zone=self.zone, south=False)
        e0, n0, _, _ = geo.latlon_to_utm(lat0, lon0, zone=self.zone, south=False)
        self.X, self.Y = e - e0, n - n0
        nk = np.asarray(d["nat_u"]).shape[0]
        f = lambda k: np.asarray(d[k], float).reshape(nk, -1)
        p, T, q = f("nat_pres"), f("nat_t"), f("nat_q")
        u, v = grid_to_utm(f("nat_u"), f("nat_v"), lat.ravel()[None], lon.ravel()[None], self.zone)
        tv = T * (1.0 + 0.61 * q)
        w = -f("nat_w") * RD * tv / (p * G)
        thv = T * (1e5 / p) ** KAPPA_P * (1.0 + 0.61 * q)
        tke = np.maximum(f("nat_tke"), 1e-4)
        self.orog = np.asarray(d["sfc_orog"], float).ravel()
        self.z0 = np.asarray(d["sfc_fsr"], float).ravel()
        self.h_levels = f("nat_gh") - self.orog[None]
        if (np.diff(self.h_levels, axis=0) <= 0).any() or (self.h_levels[0] <= self.z0).any():
            raise HrrrOuterError(f"{self.label}: level heights must rise from above z0 in every column")
        self.clw = f("nat_clwmr")
        self.d = d
        h = self.h_levels
        self.hg = HG
        self.h_top = float(self.h_levels[-1].min())
        self.tab = {}
        for name, a3 in (("u", u), ("v", v), ("w", w), ("theta", thv), ("tke", tke)):
            g = np.empty((h.shape[1], self.hg.size))
            for c in range(h.shape[1]):
                g[c] = np.interp(self.hg, h[:, c], a3[:, c])
                below = self.hg < h[0, c]
                if name in ("u", "v"):
                    zz = np.maximum(self.hg[below], self.z0[c] * 1.0001)
                    g[c, below] = a3[0, c] * np.log(zz / self.z0[c]) / math.log(h[0, c] / self.z0[c])
                elif name == "w":
                    g[c, below] = a3[0, c] * self.hg[below] / h[0, c]
            self.tab[name] = g
        nj, ni = self.shape
        jj, ii = np.mgrid[0:nj, 0:ni]
        A = np.column_stack([np.ones(self.X.size), self.X.ravel(), self.Y.ravel()])
        self._ci = np.linalg.lstsq(A, ii.ravel().astype(float), rcond=None)[0]
        self._cj = np.linalg.lstsq(A, jj.ravel().astype(float), rcond=None)[0]

    # the map from grid indices to (x, y) is bilinear in each cell; invert it by Newton's method
    def locate(self, x, y):
        x, y = np.asarray(x, float).ravel(), np.asarray(y, float).ravel()
        nj, ni = self.shape
        fi = self._ci[0] + self._ci[1] * x + self._ci[2] * y
        fj = self._cj[0] + self._cj[1] * x + self._cj[2] * y
        X, Y = self.X, self.Y
        for _ in range(10):
            i0 = np.clip(np.floor(fi).astype(int), 0, ni - 2)
            j0 = np.clip(np.floor(fj).astype(int), 0, nj - 2)
            a, b = fi - i0, fj - j0
            corners = []
            for F in (X, Y):
                f00, f10, f01, f11 = F[j0, i0], F[j0, i0 + 1], F[j0 + 1, i0], F[j0 + 1, i0 + 1]
                val = (1 - a) * (1 - b) * f00 + a * (1 - b) * f10 + (1 - a) * b * f01 + a * b * f11
                da = (1 - b) * (f10 - f00) + b * (f11 - f01)
                db = (1 - a) * (f01 - f00) + a * (f11 - f10)
                corners.append((val, da, db))
            (xv, xa, xb), (yv, ya, yb) = corners
            rx, ry = x - xv, y - yv
            det = xa * yb - xb * ya
            fi = fi + (yb * rx - xb * ry) / det
            fj = fj + (-ya * rx + xa * ry) / det
        tol = 1e-6
        out = ~((fi >= -tol) & (fi <= ni - 1 + tol) & (fj >= -tol) & (fj <= nj - 1 + tol))   # NaN is out
        if out.any():
            k = int(np.argmax(out))
            raise HrrrOuterError(f"{self.label}: the point ({x[k]:.0f}, {y[k]:.0f}) m lies outside the crop; fetch a "
                                 "larger one (--half-km)")
        i0 = np.clip(np.floor(fi).astype(int), 0, ni - 2)
        j0 = np.clip(np.floor(fj).astype(int), 0, nj - 2)
        return i0, j0, fi - i0, fj - j0

    def at(self, name, x2, y2, h3):
        """name at horizontal points x2, y2 (same shape S) and heights above ground h3 (shape (nz,) + S)."""
        x2 = np.asarray(x2, float)
        h3 = np.asarray(h3, float)
        i0, j0, a, b = self.locate(x2, y2)
        ni = self.shape[1]
        hh = h3.reshape(h3.shape[0], -1)
        hg = self.hg
        k = np.clip(np.searchsorted(hg, hh, side="right") - 1, 0, hg.size - 2)
        t = np.clip((hh - hg[k]) / (hg[k + 1] - hg[k]), 0.0, 1.0)
        tab = self.tab[name]
        out = np.zeros_like(hh)
        for dj, di, wgt in ((0, 0, (1 - a) * (1 - b)), (0, 1, a * (1 - b)), (1, 0, (1 - a) * b), (1, 1, a * b)):
            col = ((j0 + dj) * ni + i0 + di)[None, :]
            out += wgt[None, :] * ((1 - t) * tab[col, k] + t * tab[col, k + 1])
        return out.reshape(h3.shape)

    def box(self, half):
        """Columns whose (x, y) lie within the square of half-width half (m)."""
        return (np.abs(self.X.ravel()) <= half) & (np.abs(self.Y.ravel()) <= half)

    def mean_profile(self, name, half):
        m = self.box(half)
        if not m.any():
            raise HrrrOuterError(f"{self.label}: no HRRR column within {half:.0f} m of the centre")
        return self.tab[name][m].mean(axis=0)

    def check_top(self, top, what="the domain top"):
        """Refuse a top (m above ground) above the height table (HTOP) or above HRRR's highest level."""
        if not (math.isfinite(top) and 0 < top <= HTOP):
            raise HrrrOuterError(f"{what} ({top:g} m) must be > 0 and at most {HTOP:.0f} m (the top of the "
                                 "height table)")
        if top > self.h_top:
            raise HrrrOuterError(f"{what} ({top:g} m) lies above HRRR's highest fetched level "
                                 f"({self.h_top:.0f} m above ground)")

    def check_cloud(self, half, top, limit=CLOUD_MIXING_RATIO):
        """Refuse the hour if cloud water passes limit (kg/kg) in a column within half + 3 km, below top (m).
        The crop must reach that far."""
        if not (math.isfinite(limit) and limit > 0):
            raise HrrrOuterError(f"--cloud-mixing-ratio ({limit}) must be > 0")
        r = half + 3000.0
        self.locate(np.array([-r, r, -r, r]), np.array([-r, -r, r, r]))
        m = self.box(half + 3000.0)       # and the HRRR columns just outside, whose cells reach in
        low = self.h_levels[:, m] <= top
        c = np.where(low, self.clw[:, m], 0.0)
        if c.max() > limit:
            k, n = np.unravel_index(np.argmax(c), c.shape)
            raise HrrrOuterError(f"{self.label}: cloud water {1e3 * c.max():.3g} g/kg at "
                                 f"{self.h_levels[k, m][n]:.0f} m above ground inside the domain (limit "
                                 f"{1e3 * limit:.3g} g/kg below {top:.0f} m); this hour is refused: GUST is dry")

    def geostrophic(self, half, lat0, z_lo=200.0, z_hi=1000.0):
        """Geostrophic wind (UTM components) from a plane fit of p at constant height over the box, mean over
        z_lo..z_hi above the mean ground (the 'x, y' of the fit are UTM, so G is in the same frame)."""
        d = self.d
        nk = self.h_levels.shape[0]
        fl = lambda k: np.asarray(d[k], float).reshape(nk, -1)
        p, gh, T, q = fl("nat_pres"), fl("nat_gh"), fl("nat_t"), fl("nat_q")
        m = self.box(half)
        f = 2.0 * OMEGA * math.sin(math.radians(lat0))
        g0 = self.orog[m].mean()
        lowest = float(gh[0, m].max())
        if g0 + z_lo < lowest:
            raise HrrrOuterError(f"{self.label}: geostrophic wind: {z_lo:.0f} m above the box-mean ground "
                                 f"({g0 + z_lo:.0f} m) lies below HRRR's first level in a column of the box "
                                 f"({lowest:.0f} m); the plane fit of p would use extrapolated pressures")
        A = np.column_stack([self.X.ravel()[m], self.Y.ravel()[m], np.ones(m.sum())])
        out = []
        for h in np.arange(z_lo, z_hi + 1.0, 100.0):
            pz = np.array([math.exp(np.interp(g0 + h, gh[:, i], np.log(p[:, i]))) for i in np.nonzero(m)[0]])
            tv = np.array([np.interp(g0 + h, gh[:, i], T[:, i] * (1 + 0.61 * q[:, i])) for i in np.nonzero(m)[0]])
            rho = (pz / (RD * tv)).mean()
            c = np.linalg.lstsq(A, pz, rcond=None)[0]
            out.append((-c[1] / (rho * f), c[0] / (rho * f)))
        return tuple(np.mean(out, axis=0))

    def surface_flux(self, half):
        """Box-mean kinematic virtual heat flux w'theta_v', K m/s (upward positive): the sensible heat flux H and
        the latent heat flux LE, H / (rho c_p) (1 + 0.61 q) + 0.61 theta LE / (rho L_v), with rho = p_s / (R_d T_v)
        at 2 m (GUST's theta here is theta_v)."""
        d, m = self.d, self.box(half)
        f = lambda k: np.asarray(d[k], float).ravel()[m]
        sp, t2, q2 = f("sfc_sp"), f("sfc_t2m"), f("sfc_sh2")
        rho = sp / (RD * t2 * (1.0 + 0.61 * q2))
        theta = t2 * (1e5 / sp) ** KAPPA_P
        return float((f("sfc_ishf") / (rho * CP) * (1.0 + 0.61 * q2) + 0.61 * theta * f("sfc_slhtf") / (rho * LV))
                     .mean())

    def surface_theta(self):
        """Per column, the ground's virtual potential temperature theta_v,s = T_skin (1e5 / p_s)^0.2857
        (1 + 0.61 q_2m), K: HRRR's skin temperature (TMP at the surface), its surface pressure and its 2 m specific
        humidity (the skin's own humidity is not in HRRR's output). GUST's theta is theta_v."""
        if "sfc_t" not in self.d:
            raise HrrrOuterError(f"{self.label}: the crop has no skin temperature (sfc_t); fetch it again with "
                                 "'hrrr_outer.py fetch'")
        f = lambda k: np.asarray(self.d[k], float).ravel()
        return f("sfc_t") * (1e5 / f("sfc_sp")) ** KAPPA_P * (1.0 + 0.61 * f("sfc_sh2"))

    def pblh(self):
        """Per column, HRRR's boundary layer height z_i (HPBL, m above ground): the depth in Beljaars's convective
        velocity w* = (g / theta_0 z_i q)^(1/3) (bc.surface_wstar = beljaars)."""
        if "sfc_blh" not in self.d:
            raise HrrrOuterError(f"{self.label}: the crop has no boundary layer height (sfc_blh); fetch it again "
                                 "with 'hrrr_outer.py fetch'")
        return np.asarray(self.d["sfc_blh"], float).ravel()

    def write_pblh(self, path, half, lo=10.0, hi=10000.0):
        """Write the columns within half (m) of the centre as 'x y z_i' lines (bc.surface_pblh_file); returns
        (count, min, max, mean). A value outside [lo, hi] m (bc.surface_pblh's range) refuses the hour."""
        return self._write_column_map(path, half, self.pblh(), "boundary layer height", "m", lo, hi, "{:.1f}")

    def _write_column_map(self, path, half, values, what, unit, lo, hi, fmt):
        m = self.box(half)
        if not m.any():
            raise HrrrOuterError(f"{self.label}: no HRRR column within {half:.0f} m of the centre")
        v = values[m]
        if not (np.isfinite(v).all() and (v >= lo).all() and (v <= hi).all()):
            raise HrrrOuterError(f"{self.label}: the {what} ({np.nanmin(v):.2f} to {np.nanmax(v):.2f} {unit}) "
                                 f"must be finite and in [{lo:g}, {hi:g}] {unit}")
        with open(path, "w") as out:
            for x, y, t in zip(self.X.ravel()[m], self.Y.ravel()[m], v):
                out.write(f"{x:.2f} {y:.2f} {fmt.format(t)}\n")
        return int(m.sum()), float(v.min()), float(v.max()), float(v.mean())

    def write_surface_theta(self, path, half, lo=150.0, hi=350.0):
        """Write the columns within half (m) of the centre as 'x y theta_v,s' lines (bc.surface_temperature_file);
        returns (count, min, max, mean). A value outside [lo, hi] K refuses the hour."""
        return self._write_column_map(path, half, self.surface_theta(), "skin theta_v", "K", lo, hi, "{:.3f}")


# ---------------------------------------------------------------------------------------------------- checkpoints

def chk_header(chk):
    """n_cell, prob_lo, prob_hi of a GUST checkpoint Header (hex floats)."""
    h = {}
    path = os.path.join(chk, "Header")
    if not os.path.isfile(path):
        raise HrrrOuterError(f"{chk}: no Header (not a GUST checkpoint)")
    with open(path) as f:
        first = f.readline().split()
        if first[:1] != ["gust-checkpoint"]:
            raise HrrrOuterError(f"{path}: not a GUST checkpoint header")
        for line in f:
            s = line.split()
            if s and s[0] in ("n_cell", "prob_lo", "prob_hi"):
                h[s[0]] = [int(x) if s[0] == "n_cell" else float.fromhex(x) for x in s[1:4]]
    if len(h) != 3:
        raise HrrrOuterError(f"{path}: no n_cell, prob_lo or prob_hi")
    return h


def chk_fabs(chk, name, mode="r"):
    """(lo, hi, memmap [ncomp, nz, ny, nx]) of every FAB of <chk>/Level_0/<name>, ghosts included; with mode r+
    the arrays write through to the files. Only native little-endian doubles (GUST's checkpoints) are taken."""
    base = os.path.join(chk, "Level_0", name)
    if not os.path.isfile(base + "_H"):
        raise HrrrOuterError(f"{base}_H: no such file")
    with open(base + "_H") as f:
        entries = re.findall(r"FabOnDisk:\s+(\S+)\s+(\d+)", f.read())
    out = []
    for fname, off in entries:
        fp = os.path.join(chk, "Level_0", fname)
        with open(fp, "rb") as f:
            f.seek(int(off))
            line = f.readline().decode()
            start = f.tell()
        m = re.search(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \((\d),(\d),(\d)\)\)\s+(\d+)", line)
        if not m or "(8, (8 7 6 5 4 3 2 1))" not in line:
            raise HrrrOuterError(f"{fp}: FAB header {line.strip()!r}: little-endian doubles expected")
        g = [int(x) for x in m.groups()]
        lo, hi, nc = np.array(g[0:3]), np.array(g[3:6]), g[9]
        if nc != 1:
            raise HrrrOuterError(f"{fp}: {nc} components; one expected")
        n = hi - lo + 1
        out.append((lo, hi, np.memmap(fp, dtype="<f8", mode=mode, offset=start, shape=(nc, n[2], n[1], n[0]))))
    return out


def chk_global(chk, name, shape):
    """The domain part of a checkpoint MultiFab as one array of the given (nz, ny, nx) shape."""
    g = np.full(shape, np.nan)
    top = np.array(shape[::-1]) - 1
    for flo, fhi, a in chk_fabs(chk, name):
        i0, j0, k0 = np.maximum(flo, 0)
        i1, j1, k1 = np.minimum(fhi, top)
        if (np.array([i1, j1, k1]) < [i0, j0, k0]).any():
            continue
        g[k0:k1 + 1, j0:j1 + 1, i0:i1 + 1] = a[0, k0 - flo[2]:k1 - flo[2] + 1, j0 - flo[1]:j1 - flo[1] + 1,
                                                i0 - flo[0]:i1 - flo[0] + 1]
    if np.isnan(g).any():
        raise HrrrOuterError(f"{chk}: {name} does not cover the domain")
    return g


def _avg(a, ax):
    n = a.shape[ax]
    return 0.5 * (np.take(a, range(1, n), ax) + np.take(a, range(0, n - 1), ax))


def write_checkpoint(field, chk, cloud_limit=CLOUD_MIXING_RATIO):
    """Overwrite u, v, w, theta and tke of a GUST checkpoint in place with HRRR's (p and the heights are kept).

    Every point inside the domain takes HRRR at its own x, y and height above the checkpoint's ground (the mean of
    the surrounding nodes minus the mean of their ground nodes); w = 0 on the ground faces. Ghost points outside the
    domain keep the first guess. Returns {field: (min, max, mean)} and the divergence of the result."""
    hd = chk_header(chk)
    nx, ny, nz = hd["n_cell"]
    lo, hi = np.array(hd["prob_lo"]), np.array(hd["prob_hi"])
    half = max(abs(lo[0]), abs(hi[0]), abs(lo[1]), abs(hi[1]))
    field.check_top(hi[2], f"{chk}: the domain top")
    field.check_cloud(half, hi[2], cloud_limit)
    dx, dy = (hi[0] - lo[0]) / nx, (hi[1] - lo[1]) / ny
    zn = chk_global(chk, "z_phys_nd", (nz + 1, ny + 1, nx + 1))
    xn, yn = lo[0] + dx * np.arange(nx + 1), lo[1] + dy * np.arange(ny + 1)
    xc, yc = _avg(xn, 0), _avg(yn, 0)
    agl = zn - zn[0:1]
    pts = {"u": (xn, yc, _avg(_avg(agl, 0), 1)), "v": (xc, yn, _avg(_avg(agl, 0), 2)),
           "w": (xc, yc, _avg(_avg(agl, 1), 2)), "theta": (xc, yc, _avg(_avg(_avg(agl, 0), 1), 2))}
    pts["tke"] = pts["theta"]
    glob_ = {}
    stats = {}
    for name in FIELDS:            # every field first: a refusal leaves the checkpoint as it was
        xs, ys, h3 = pts[name]
        x2, y2 = np.meshgrid(xs, ys)
        g = field.at(name, x2, y2, h3)
        if name == "w":
            g[0] = 0.0
        glob_[name] = g
        stats[name] = (float(g.min()), float(g.max()), float(g.mean()))
    for name in FIELDS:
        g = glob_[name]
        top = np.array(g.shape[::-1]) - 1
        for flo, fhi, a in chk_fabs(chk, name, mode="r+"):
            i0, j0, k0 = np.maximum(flo, 0)
            i1, j1, k1 = np.minimum(fhi, top)
            if (np.array([i1, j1, k1]) < [i0, j0, k0]).any():
                continue
            a[0, k0 - flo[2]:k1 - flo[2] + 1, j0 - flo[1]:j1 - flo[1] + 1, i0 - flo[0]:i1 - flo[0] + 1] = \
                g[k0:k1 + 1, j0:j1 + 1, i0:i1 + 1]
            a.flush()
    # the discrete divergence in flux form on the terrain-following cells (the face areas follow the terrain;
    # w is the vertical velocity through the horizontal projection of the top and bottom faces)
    dzc = np.diff(zn, axis=0)
    dzx, dzy, dz = _avg(dzc, 1), _avg(dzc, 2), _avg(_avg(dzc, 1), 2)
    div = (np.diff(glob_["u"] * dzx, axis=2) * dy + np.diff(glob_["v"] * dzy, axis=1) * dx
           + np.diff(glob_["w"], axis=0) * dx * dy) / (dx * dy * dz)
    stats["divergence_rms"] = float(np.sqrt((div ** 2).mean()))
    return stats


def balance(chk, m=0):
    """Volume fluxes through the faces of a checkpoint's velocity, m cells in from each lateral face."""
    hd = chk_header(chk)
    nx, ny, nz = hd["n_cell"]
    lo, hi = np.array(hd["prob_lo"]), np.array(hd["prob_hi"])
    if not (isinstance(m, int) and 0 <= m < min(nx, ny) // 2):
        raise HrrrOuterError(f"--cells ({m}) must be an integer from 0 to {min(nx, ny) // 2 - 1}")
    dx, dy = (hi[0] - lo[0]) / nx, (hi[1] - lo[1]) / ny
    zn = chk_global(chk, "z_phys_nd", (nz + 1, ny + 1, nx + 1))
    u = chk_global(chk, "u", (nz, ny, nx + 1))
    v = chk_global(chk, "v", (nz, ny + 1, nx))
    w = chk_global(chk, "w", (nz + 1, ny, nx))
    dzc = np.diff(zn, axis=0)
    dzx, dzy = _avg(dzc, 1), _avg(dzc, 2)
    s, t = slice(m, ny - m), slice(m, nx - m)
    fw = (u[:, s, m] * dzx[:, s, m]).sum() * dy
    fe = (u[:, s, nx - m] * dzx[:, s, nx - m]).sum() * dy
    fs = (v[:, m, t] * dzy[:, m, t]).sum() * dx
    fn = (v[:, ny - m, t] * dzy[:, ny - m, t]).sum() * dx
    net = fw - fe + fs - fn
    entering = sum(max(q, 0.0) for q in (fw, -fe, fs, -fn))
    return {"west": fw, "east": -fe, "south": fs, "north": -fn, "net_lateral_inflow": net,
            "net_fraction": net / entering if entering > 0 else float("nan"),
            "top_outflow": float(w[nz, s, t].sum() * dx * dy)}


# ---------------------------------------------------------------------------------------------------- profile

PROFILE_HEIGHTS = (0., 10., 40., 90., 175., 300., 450., 650., 850., 1000., 1250., 1500., 2000., 2500., 3000.)


def theta_profile(field, half, top):
    """Box-mean theta_v at PROFILE_HEIGHTS up to the top (and the top), and its gradient over the last 100 m."""
    field.check_top(top, "--top")
    if not (math.isfinite(half) and half > 0):
        raise HrrrOuterError(f"--half-width ({half}) must be > 0")
    p = field.mean_profile("theta", half)
    z = [h for h in PROFILE_HEIGHTS if h < top - 1.0] + [float(top)]
    th = [float(np.interp(h, field.hg, p)) for h in z]
    grad = float((np.interp(top, field.hg, p) - np.interp(top - 100.0, field.hg, p)) / 100.0)
    return z, th, grad


# ---------------------------------------------------------------------------------------------------- score

def _read_plotfile(plt):
    from gust_library import plotfile
    try:
        c, hdr = plotfile.read_cells(plt, ["u", "v", "w", "tke", "z"])
        zn = plotfile.read_nodes(plt)
    except ValueError as err:
        raise HrrrOuterError(str(err))
    T = lambda a: np.transpose(a, (2, 1, 0))
    return {k: T(a) for k, a in c.items()}, T(zn), hdr


def _wrap(a):
    return (np.asarray(a) + 180.0) % 360.0 - 180.0


def _from_dir(u, v):
    return float(np.degrees(np.arctan2(-u, -v)) % 360.0)


def score(field, plt, margin=2500.0, log=None, low=500.0):
    """The nest against HRRR on the nest's cells, more than margin (m) from the lateral faces."""
    if not (math.isfinite(margin) and margin >= 0):
        raise HrrrOuterError(f"--margin ({margin}) must be >= 0")
    c, zn, hdr = _read_plotfile(plt)
    nz, ny, nx = c["u"].shape
    lo, hi = hdr["prob_lo"], hdr["prob_hi"]
    dx, dy = (hi[0] - lo[0]) / nx, (hi[1] - lo[1]) / ny
    xc, yc = lo[0] + dx * (np.arange(nx) + 0.5), lo[1] + dy * (np.arange(ny) + 0.5)
    ground = 0.25 * (zn[0, :-1, :-1] + zn[0, 1:, :-1] + zn[0, :-1, 1:] + zn[0, 1:, 1:])
    agl = c["z"] - ground[None]
    ix = (xc - lo[0] > margin) & (hi[0] - xc > margin)
    iy = (yc - lo[1] > margin) & (hi[1] - yc > margin)
    if not ix.any() or not iy.any():
        raise HrrrOuterError(f"--margin ({margin}) leaves no interior cell")
    X, Y = np.meshgrid(xc[ix], yc[iy])
    sub = lambda a: a[:, iy][:, :, ix]
    h = sub(agl)
    u, v, w, k = (sub(c[n]) for n in ("u", "v", "w", "tke"))
    H = {n: field.at(n, X, Y, h) for n in ("u", "v", "w", "tke")}
    S, SH = np.hypot(u, v), np.hypot(H["u"], H["v"])
    dS = S - SH
    dV2 = (u - H["u"]) ** 2 + (v - H["v"]) ** 2
    levels = []
    for kk in range(nz):
        um, vm, uh, vh = u[kk].mean(), v[kk].mean(), H["u"][kk].mean(), H["v"][kk].mean()
        levels.append({"level": kk, "agl": float(h[kk].mean()), "speed": float(S[kk].mean()),
                       "speed_hrrr": float(SH[kk].mean()), "bias": float(dS[kk].mean()),
                       "rms": float(np.sqrt((dS[kk] ** 2).mean())), "vector_rms": float(np.sqrt(dV2[kk].mean())),
                       "dir": _from_dir(um, vm), "dir_hrrr": _from_dir(uh, vh),
                       "k": float(k[kk].mean()), "k_hrrr": float(H["tke"][kk].mean())})
    m = h <= low
    out = {"hour": field.label, "plotfile": os.path.abspath(plt), "margin": margin, "interior_cells": int(dS.size),
           "bias": float(dS.mean()), "rms": float(np.sqrt((dS ** 2).mean())),
           "low_top": low, "low_speed_bias": float(dS[m].mean()), "low_vector_rms": float(np.sqrt(dV2[m].mean())),
           "low_dir_diff": float(_wrap(_from_dir(u[m].mean(), v[m].mean()) - _from_dir(H["u"][m].mean(),
                                                                                    H["v"][m].mean()))),
           "w_rms": float(np.sqrt((w ** 2).mean())), "w_max": float(np.abs(w).max()),
           "w_rms_hrrr": float(np.sqrt((H["w"] ** 2).mean())),
           "k_ratio": float(k[m].sum() / max(H["tke"][m].sum(), 1e-30)), "levels": levels,
           "iterations": None, "steady": None}
    if log:
        with open(log, errors="replace") as f:
            t = f.read()
        it = re.findall(r"stage account: (\d+) iterations", t)
        out["iterations"] = int(it[-1]) if it else None
        v = re.findall(r"^gust: (steady|not steady|diverged) after", t, re.M)
        out["steady"] = (v[-1] == "steady") if v else None
    return out


def score_text(s):
    lines = [f"{s['hour']}; {os.path.basename(s['plotfile'])}; {s['interior_cells']} cells more than "
             f"{s['margin']:.0f} m from the faces",
             " level  AGL(m)  speed nest  HRRR   bias   rms  vec rms   dir nest  HRRR    k nest   HRRR"]
    for r in s["levels"]:
        lines.append(f"{r['level']:5d} {r['agl']:7.0f} {r['speed']:10.2f} {r['speed_hrrr']:6.2f} {r['bias']:6.2f} "
                     f"{r['rms']:5.2f} {r['vector_rms']:8.2f} {r['dir']:10.0f} {r['dir_hrrr']:5.0f} "
                     f"{r['k']:9.3f} {r['k_hrrr']:6.3f}")
    st = {None: "unknown (no log)", True: "yes", False: "no"}[s["steady"]]
    lines += [f"all heights: speed bias {s['bias']:.2f}, rms {s['rms']:.2f} m/s",
              f"below {s['low_top']:.0f} m: speed bias {s['low_speed_bias']:.2f} m/s, vector rms "
              f"{s['low_vector_rms']:.2f} m/s, direction nest - HRRR {s['low_dir_diff']:.1f} degrees, "
              f"k / k_HRRR {s['k_ratio']:.2f}",
              f"w: rms {s['w_rms']:.3f} (HRRR {s['w_rms_hrrr']:.3f}), max |w| {s['w_max']:.2f} m/s",
              f"iterations {s['iterations']}, steady {st}"]
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------------- fetch

NAT_REGEX = r":(UGRD|VGRD|HGT|TMP|SPFH|PRES|VVEL|TKE|CLMR|CLWMR):([1-9]|[1-2][0-9]|30) hybrid level"
SFC_REGEX = r":(SHTFL|LHTFL|FRICV|HPBL|SFCR):|:(HGT|TMP|PRES):surface|:(TMP|SPFH):2 m"


def fetch(cycle, fxx, lat0, lon0, half_km, save, out):
    """Download HRRR's wrfnat and wrfsfc subsets of one cycle and forecast hour (byte ranges, Herbie) and save the
    crop around (lat0, lon0) with a margin for HRRR's rotation against the UTM grid."""
    check_site(lat0, lon0)
    check_fxx(fxx)
    if not (math.isfinite(half_km) and half_km > 0):
        raise HrrrOuterError(f"--half-km ({half_km}) must be > 0")
    init = hrrr_regime.parse_time(cycle, "--cycle")
    import warnings
    warnings.filterwarnings("ignore")
    from herbie import Herbie
    res, sl = {}, None
    when = init.strftime("%Y-%m-%d %H:%M")
    for product, regex, prefix in (("nat", NAT_REGEX, "nat_"), ("sfc", SFC_REGEX, "sfc_")):
        dss = Herbie(when, model="hrrr", product=product, fxx=fxx, save_dir=save, verbose=False).xarray(
            regex, remove_grib=False)
        for ds in dss if isinstance(dss, list) else [dss]:
            if sl is None:
                lat = ds.latitude.values
                lon = np.where(ds.longitude.values > 180, ds.longitude.values - 360, ds.longitude.values)
                d2 = (lat - lat0) ** 2 + ((lon - lon0) * math.cos(math.radians(lat0))) ** 2
                j, i = np.unravel_index(np.argmin(d2), d2.shape)
                n = int(math.ceil(math.sqrt(2.0) * half_km / 3.0)) + 3
                if j - n < 0 or i - n < 0 or j + n + 1 > lat.shape[0] or i + n + 1 > lat.shape[1]:
                    raise HrrrOuterError(f"--half-km ({half_km}): the crop passes the edge of the HRRR domain")
                sl = (slice(j - n, j + n + 1), slice(i - n, i + n + 1))
                res["lat"], res["lon"] = lat[sl], lon[sl]
            for name in ds.data_vars:
                a = ds[name].values
                res[prefix + name] = a[(slice(None),) + sl] if a.ndim == 3 else a[sl]
    res["init"] = np.array(hrrr_regime.fmt_time(init))
    res["fxx"] = np.array(int(fxx))
    res["source"] = np.array(f"HRRR {init:%Y-%m-%d %H} UTC f{fxx:02d} (wrfnat, wrfsfc), NOAA, via Herbie")
    np.savez_compressed(out, **res)
    return res


# ---------------------------------------------------------------------------------------------------- command line

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def site(p):
        p.add_argument("--lat", type=float, required=True, help="latitude of the domain centre, 21 to 53")
        p.add_argument("--lon", type=float, required=True, help="longitude of the domain centre, -135 to -60")
    f = sub.add_parser("fetch", help="download one HRRR hour and save its crop")
    f.add_argument("--cycle", required=True, help="initial time of the HRRR run, YYYY-MM-DDTHH:MMZ")
    f.add_argument("--fxx", type=int, required=True, help="forecast hour, 1 or later")
    site(f)
    f.add_argument("--half-km", type=float, required=True, help="half-width of the square to cover, km")
    f.add_argument("--save", required=True, help="Herbie download folder")
    f.add_argument("--out", required=True, help="crop to write (.npz)")
    c = sub.add_parser("to-chk", help="overwrite a checkpoint's u, v, w, theta, tke with HRRR's")
    c.add_argument("crop")
    c.add_argument("chk")
    site(c)
    c.add_argument("--cloud-mixing-ratio", type=float, default=CLOUD_MIXING_RATIO,
                   help="cloud water above which an hour is refused, kg/kg (%g)" % CLOUD_MIXING_RATIO)
    b = sub.add_parser("balance", help="volume fluxes through a checkpoint's faces")
    b.add_argument("chk")
    b.add_argument("--cells", type=int, default=0, help="cells in from each lateral face (0)")
    p = sub.add_parser("profile", help="box-mean theta_v profile, geostrophic wind and surface flux")
    p.add_argument("crop")
    site(p)
    p.add_argument("--half-width", type=float, required=True, help="half-width of the box, m")
    p.add_argument("--top", type=float, default=MIN_TOP, help="domain top, m (%g)" % MIN_TOP)
    s = sub.add_parser("score", help="a nest plotfile against HRRR")
    s.add_argument("crop")
    s.add_argument("plotfile")
    site(s)
    s.add_argument("--log", help="the nest's log (iterations and the steady verdict)")
    s.add_argument("--margin", type=float, default=2500.0, help="width left out beside each face, m (2500)")
    s.add_argument("--json", help="write the score as JSON")
    k = sub.add_parser("case", help="decks of the chain from a gust_case YAML with an hrrr section")
    k.add_argument("yaml")
    k.add_argument("--output-dir", help="where the case folder is made (default: case.output_dir)")
    k.add_argument("--overwrite", action="store_true", help="replace the files of an earlier write")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "fetch":
            res = fetch(a.cycle, a.fxx, a.lat, a.lon, a.half_km, a.save, a.out)
            print(f"wrote {a.out}: {res['lat'].shape[0]} x {res['lat'].shape[1]} columns, "
                  f"{res['nat_u'].shape[0]} levels")
        elif a.cmd == "to-chk":
            st = write_checkpoint(HrrrField(a.crop, a.lat, a.lon), a.chk, a.cloud_mixing_ratio)
            for n in FIELDS:
                print(f"{n:6s} min {st[n][0]:9.4f}  max {st[n][1]:9.4f}  mean {st[n][2]:9.4f}")
            print(f"divergence rms {st['divergence_rms']:.3e} 1/s")
        elif a.cmd == "balance":
            r = balance(a.chk, a.cells)
            print(f"{a.chk}: flux in through W {r['west']:.3e}, E {r['east']:.3e}, S {r['south']:.3e}, "
                  f"N {r['north']:.3e} m3/s; net lateral inflow {r['net_lateral_inflow']:.3e} = "
                  f"{r['net_fraction']:.3f} of entering; top outflow {r['top_outflow']:.3e}")
        elif a.cmd == "profile":
            fld = HrrrField(a.crop, a.lat, a.lon)
            z, th, grad = theta_profile(fld, a.half_width, a.top)
            ug, vg = fld.geostrophic(a.half_width, a.lat)
            print(json.dumps({"hour": fld.label, "theta_z": z, "theta": th, "top_gradient": grad,
                              "geostrophic": [ug, vg], "kinematic_heat_flux": fld.surface_flux(a.half_width)},
                             indent=1))
        elif a.cmd == "score":
            sc = score(HrrrField(a.crop, a.lat, a.lon), a.plotfile, a.margin, a.log)
            print(score_text(sc))
            if a.json:
                with open(a.json, "w") as fh:
                    json.dump(sc, fh, indent=1)
                    fh.write("\n")
        elif a.cmd == "case":
            import hrrr_outer_case
            hrrr_outer_case.main_case(a.yaml, a.output_dir, a.overwrite)
    except (HrrrOuterError, getattr(sys.modules.get("hrrr_outer"), "HrrrOuterError", HrrrOuterError),
            hrrr_regime.RegimeError) as err:
        print(f"hrrr_outer.py: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
