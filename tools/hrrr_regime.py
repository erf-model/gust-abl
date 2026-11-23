#!/usr/bin/env python3
"""Put an HRRR hour into a stratification class for a GUST run with an HRRR outer level.

    hrrr_regime.py fetch VALID_TIME... --lat LAT --lon LON [--half-km 12.5] [--grad-half-km 40] [--fxx 1]
                         --save DIR --out-dir DIR [--jobs N]
    hrrr_regime.py column CROP.npz --valid VALID_TIME --lat LAT --lon LON [--half-km 12.5] [--grad-half-km 40]
                          [--fxx 1] --out COL.json
    hrrr_regime.py classify COL.json [--before COL.json] [--after COL.json] [--relief H] [--thresholds FILE]
    hrrr_regime.py catalogue COL.json... [--jobs N] [--thresholds FILE] [--out CSV]
    hrrr_regime.py thresholds

A column is the mean over a square box of HRRR native (hybrid) levels around a point: heights above the HRRR ground,
the earth-relative wind, the virtual potential temperature, the largest cloud water in the box at each level, the
horizontal gradient of the virtual potential temperature at constant height, and the surface fields (sensible and
latent heat flux, u*, z0, PBL height, 2 m temperature and humidity, surface pressure, the relief of the HRRR ground
in the box). ``fetch`` downloads an hour with Herbie (byte-range subsets of the wrfnat and wrfsfc files; needs herbie,
xarray and cfgrib) and writes its column; ``column`` makes the column from a crop saved earlier with the same keys.
``classify`` and ``catalogue`` need numpy only.

Classes (the first rule that holds wins):

  C9 cloud water above the cloud threshold anywhere in the box below the domain top
  C7 evening transition: q_v above the flux threshold at t-1 and not at t, or at t and not at t+1
  C8 morning transition: q_v not above the flux threshold at t-1 and above it at t, or the same between t and t+1
  C1 convective: L < 0, -L <= the neutral length and -z_i/L above the convective ratio
  C2 weakly unstable: L < 0, -L <= the neutral length, -z_i/L up to the convective ratio
  C5 stable with a low-level jet: q_v below minus the flux threshold and a jet in the speed profile, whatever L
  C3 neutral: |L| above the neutral length
  C6 very stable (0 < L <= the neutral length): z_i/L above the very-stable ratio, or the bulk Richardson number
     of the lowest layer above its limit
  C4 weakly stable: the other stable hours

Modifiers, added to any class: Fr (U / (N H) < 1 over the relief H), B (geostrophic shear from the thermal wind above
its threshold), T (change of the layer-mean wind over an hour above its threshold). Transitions and T need the columns
of the hours before or after; without them the result says so and they are not tested.

q_v is the virtual surface heat flux w'theta_v' = H / (rho c_p) (1 + 0.61 q) + 0.61 theta E / (rho L_v), from HRRR's
sensible (H) and latent (E) heat flux, upward positive. L = -u*^3 theta_v / (kappa g q_v) with HRRR's u* and
kappa = 0.41 (GUST's turb.kappa). z_i is HRRR's PBL height. The thresholds are starting values, to be set from a
catalogue of hours; ``thresholds`` prints them and a JSON (or YAML) file replaces any of them.
"""
import argparse
import csv
import datetime as dt
import json
import math
import os
import sys

import numpy as np

G, RD, CP, LV, KAPPA_P, KAPPA_VK, OMEGA, R_EARTH = 9.81, 287.04, 1005.0, 2.5e6, 0.2857, 0.41, 7.292e-5, 6.371e6
CLASSES = ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")
MODIFIERS = ("Fr", "B", "T")

# name: (default, unit, lower bound, upper bound or None, lower bound included); each is checked when read.
THRESHOLDS = {
    "flux_eps": (0.005, "K m/s", 0.0, None, True),          # |q_v| up to this: weak flux, Ri_b decides (C3, C6)
    "neutral_length": (500.0, "m", 0.0, None, False),        # |L| above this is neutral
    "convective_ratio": (10.0, "-", 0.0, None, False),       # -z_i/L above this is convective
    "very_stable_ratio": (1.0, "-", 0.0, None, False),       # z_i/L above this is very stable
    "rib_top": (100.0, "m", 0.0, None, False),               # top of the lowest layer for the bulk Richardson number
    "rib_crit": (0.25, "-", 0.0, None, False),               # its limit: above this is very stable
    "n_depth": (500.0, "m", 0.0, None, False),               # layer above z_i for N
    "jet_nose_max": (1000.0, "m", 0.0, None, False),         # highest jet nose
    "jet_depth": (1000.0, "m", 0.0, None, False),            # depth above the nose searched for the fall-off
    "jet_drop": (2.0, "m/s", 0.0, None, False),              # fall-off of speed above the nose
    "jet_drop_frac": (0.1, "-", 0.0, 1.0, False),            # the same as a fraction of the nose speed
    "cloud_top": (3000.0, "m", 0.0, None, False),            # domain top: cloud above it is ignored
    "cloud_mixing_ratio": (1e-5, "kg/kg", 0.0, None, False), # cloud water above this is a cloud
    "thermal_depth": (1000.0, "m", 0.0, None, False),        # layer above z_i for the geostrophic shear
    "geo_shear": (2.0, "m/s per km", 0.0, None, False),      # geostrophic shear above this is baroclinic
    "tendency_top": (1000.0, "m", 0.0, None, False),         # top of the layer for the wind tendency
    "tendency": (1.0, "m/s per h", 0.0, None, False),        # change of the layer-mean wind above this is unsteady
    "froude_min_depth": (100.0, "m", 0.0, None, False),      # least depth of the layer for U and N in Fr
    "froude": (1.0, "-", 0.0, None, False),                  # U / (N H) below this is low Froude number
}


class RegimeError(Exception):
    pass


def default_thresholds():
    return {k: v[0] for k, v in THRESHOLDS.items()}


def check_thresholds(t):
    """Every key known, every value finite and in range; aborts naming the key."""
    out = default_thresholds()
    for k, v in t.items():
        if k not in THRESHOLDS:
            raise RegimeError(f"thresholds: unknown key '{k}'")
        _, unit, lo, hi, closed = THRESHOLDS[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            raise RegimeError(f"thresholds: '{k}' must be a finite number ({unit}), got {v!r}")
        if v < lo or (v == lo and not closed) or (hi is not None and v > hi):
            rng = ("[" if closed else "(") + f"{lo}, " + (f"{hi}]" if hi is not None else "inf)")
            raise RegimeError(f"thresholds: '{k}' = {v} is outside {rng} ({unit})")
        out[k] = float(v)
    return out


def _read_text(path, what):
    try:
        with open(path) as f:
            return f.read()
    except OSError as e:
        raise RegimeError(f"{what} {path}: {e.strerror or e}")


def read_thresholds(path):
    if path is None:
        return default_thresholds()
    text = _read_text(path, "--thresholds")
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise RegimeError(f"--thresholds {path}: reading YAML needs PyYAML; give the thresholds as JSON instead")
        try:
            t = yaml.safe_load(text) or {}
        except yaml.YAMLError as e:
            raise RegimeError(f"--thresholds {path}: not YAML ({e})")
    else:
        try:
            t = json.loads(text)
        except ValueError as e:
            raise RegimeError(f"--thresholds {path}: not JSON ({e})")
    if not isinstance(t, dict):
        raise RegimeError(f"{path}: the thresholds must be a mapping of names to numbers")
    return check_thresholds(t)


# ---------------------------------------------------------------------------- the column

COLUMN_KEYS = ("z", "u", "v", "theta_v", "cloud_max", "dthdx", "dthdy")
SURFACE_KEYS = ("shf", "lhf", "ustar", "z0", "hpbl", "t2m", "q2m", "psfc", "relief")


def rotate_to_earth(u, v, lon):
    """HRRR grid-relative wind (Lambert conformal, true latitude 38.5, central longitude -97.5) to east and north."""
    a = np.radians(math.sin(math.radians(38.5)) * (np.asarray(lon) + 97.5))
    return u * np.cos(a) + v * np.sin(a), -u * np.sin(a) + v * np.cos(a)


def _pick(d, *names):
    for n in names:
        if n in d:
            return np.asarray(d[n], dtype=float)
    raise RegimeError(f"crop: none of {names} found")


def column_from_crop(d, lat0, lon0, half_km, valid, fxx, source="", grad_half_km=40.0):
    """Box mean of a crop (nat_* on hybrid levels [k, j, i], sfc_* [j, i], lat, lon).

    Profiles and surface fields are means over the square of half-width half_km; the horizontal gradient of theta_v is
    a plane fit over the larger square of half-width grad_half_km, which must lie inside the crop."""
    check_box(lat0, lon0, half_km, grad_half_km)
    lon0 = _wrap_lon(lon0)
    lat, lon = np.asarray(d["lat"], float), np.asarray(d["lon"], float)
    lon = np.where(lon > 180, lon - 360, lon)
    x = R_EARTH * math.cos(math.radians(lat0)) * np.radians(lon - lon0)
    y = R_EARTH * np.radians(lat - lat0)
    m = (np.abs(x) <= half_km * 1e3) & (np.abs(y) <= half_km * 1e3)
    if m.sum() < 4:
        raise RegimeError(f"crop: {m.sum()} HRRR cells inside the {half_km} km box; at least 4 are needed")
    mg = (np.abs(x) <= grad_half_km * 1e3) & (np.abs(y) <= grad_half_km * 1e3)
    if mg[0].any() or mg[-1].any() or mg[:, 0].any() or mg[:, -1].any():
        raise RegimeError(f"crop: the {grad_half_km} km gradient box ('--grad-half-km') reaches the edge of the crop")
    p, T, q = _pick(d, "nat_pres"), _pick(d, "nat_t"), _pick(d, "nat_q")
    ug, vg = _pick(d, "nat_u"), _pick(d, "nat_v")
    gh, clw = _pick(d, "nat_gh"), _pick(d, "nat_clwmr", "nat_clmr")
    orog = _pick(d, "sfc_orog")
    u_e, v_n = rotate_to_earth(ug, vg, lon[None])
    thv = T * (1e5 / p) ** KAPPA_P * (1 + 0.61 * q)
    h = gh - orog[None]
    nk = p.shape[0]
    xs, ys = x[mg], y[mg]
    A = np.column_stack([np.ones(xs.size), xs, ys])
    if np.linalg.matrix_rank(A) < 3:
        raise RegimeError(f"crop: {xs.size} cells in the gradient box do not define a plane")
    col = {k: [] for k in COLUMN_KEYS}
    for k in range(nk):
        col["z"].append(h[k][m].mean())
        col["u"].append(u_e[k][m].mean())
        col["v"].append(v_n[k][m].mean())
        col["theta_v"].append(thv[k][m].mean())
        col["cloud_max"].append(clw[k][m].max())
        # gradient at constant height: the level's gradient minus d(theta)/dz times the level's slope
        cth = np.linalg.lstsq(A, thv[k][mg], rcond=None)[0]
        cz = np.linalg.lstsq(A, gh[k][mg], rcond=None)[0]
        col["dthdx"].append(cth[1])
        col["dthdy"].append(cth[2])
        col["_gx"] = col.get("_gx", []) + [cz[1]]
        col["_gy"] = col.get("_gy", []) + [cz[2]]
    z = np.array(col["z"])
    th = np.array(col["theta_v"])
    dthdz = np.gradient(th, z)
    col["dthdx"] = list(np.array(col["dthdx"]) - dthdz * np.array(col.pop("_gx")))
    col["dthdy"] = list(np.array(col["dthdy"]) - dthdz * np.array(col.pop("_gy")))
    surf = {
        "shf": _pick(d, "sfc_ishf", "sfc_sshf", "sfc_shtfl")[m].mean(),
        "lhf": _pick(d, "sfc_slhtf", "sfc_lhtfl")[m].mean(),
        "ustar": _pick(d, "sfc_fricv")[m].mean(),
        "z0": _pick(d, "sfc_fsr", "sfc_sfcr")[m].mean(),
        "hpbl": _pick(d, "sfc_blh", "sfc_hpbl")[m].mean(),
        "t2m": _pick(d, "sfc_t2m")[m].mean(),
        "q2m": _pick(d, "sfc_sh2")[m].mean(),
        "psfc": _pick(d, "sfc_sp")[m].mean(),
        "relief": float(orog[m].max() - orog[m].min()),
    }
    out = {"valid_time": valid, "fxx": int(fxx), "lat": float(lat0), "lon": float(lon0), "half_km": float(half_km),
           "grad_half_km": float(grad_half_km), "cells": int(m.sum()), "source": source}
    out.update({k: [float(a) for a in col[k]] for k in COLUMN_KEYS})
    out["surface"] = {k: float(v) for k, v in surf.items()}
    check_column(out, "column")
    return out


def _wrap_lon(lon):
    return (lon + 180.0) % 360.0 - 180.0


def check_box(lat0, lon0, half_km, grad_half_km):
    """Range checks of the column's site and boxes (fetch and column refuse before any download or read)."""
    if not (math.isfinite(lat0) and 21.0 <= lat0 <= 53.0):
        raise RegimeError(f"--lat ({lat0}) must lie in HRRR's domain, 21 to 53")
    if not (math.isfinite(lon0) and -135.0 <= _wrap_lon(lon0) <= -60.0):
        raise RegimeError(f"--lon ({lon0}) must lie in HRRR's domain, -135 to -60 (or 225 to 300)")
    if not (math.isfinite(half_km) and half_km > 0):
        raise RegimeError(f"--half-km ({half_km}) must be > 0")
    if not (math.isfinite(grad_half_km) and grad_half_km >= half_km):
        raise RegimeError(f"--grad-half-km ({grad_half_km}) must be >= --half-km ({half_km})")


def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def check_column(c, name):
    for k in ("valid_time", "lat", "surface") + COLUMN_KEYS:
        if k not in c:
            raise RegimeError(f"{name}: missing '{k}'")
    for k in COLUMN_KEYS:
        if not isinstance(c[k], list) or not all(_num(x) for x in c[k]):
            raise RegimeError(f"{name}: '{k}' must be a list of numbers")
    if not isinstance(c["surface"], dict):
        raise RegimeError(f"{name}: 'surface' must be a mapping")
    for k in SURFACE_KEYS:
        if not _num(c["surface"].get(k)):
            raise RegimeError(f"{name}: surface '{k}' missing or not a number")
    if not _num(c["lat"]):
        raise RegimeError(f"{name}: 'lat' must be a number")
    n = len(c["z"])
    if n < 3:
        raise RegimeError(f"{name}: at least 3 levels are needed, got {n}")
    for k in COLUMN_KEYS:
        a = np.asarray(c[k], float)
        if a.shape != (n,):
            raise RegimeError(f"{name}: '{k}' has {a.size} values, 'z' has {n}")
        if not np.all(np.isfinite(a)):
            raise RegimeError(f"{name}: '{k}' is not finite everywhere")
    if not np.all(np.diff(c["z"]) > 0) or c["z"][0] <= 0:
        raise RegimeError(f"{name}: 'z' must be positive and increase upward")
    for k in SURFACE_KEYS:
        if not math.isfinite(c["surface"][k]):
            raise RegimeError(f"{name}: surface '{k}' is not finite")
    for k in ("t2m", "psfc"):
        if c["surface"][k] <= 0:
            raise RegimeError(f"{name}: surface '{k}' must be > 0")
    if c["surface"]["ustar"] < 1e-3:
        raise RegimeError(f"{name}: surface 'ustar' must be >= 1e-3 m/s")
    if not 0 <= c["surface"]["q2m"] < 0.1:
        raise RegimeError(f"{name}: surface 'q2m' must be in [0, 0.1) kg/kg")
    if not np.all(np.asarray(c["theta_v"], float) > 0):
        raise RegimeError(f"{name}: 'theta_v' must be > 0")
    if c["surface"]["hpbl"] <= 0:
        raise RegimeError(f"{name}: surface 'hpbl' must be > 0")
    if not -90 < c["lat"] < 90 or c["lat"] == 0:
        raise RegimeError(f"{name}: 'lat' must be in (-90, 90) and not 0 (the thermal wind needs f != 0)")
    parse_time(c["valid_time"], name)


def parse_time(s, name="time"):
    for fmt in ("%Y-%m-%dT%H:%MZ", "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M", "%Y-%m-%d %H"):
        try:
            return dt.datetime.strptime(s, fmt)
        except (ValueError, TypeError):
            pass
    raise RegimeError(f"{name}: time '{s}' is not YYYY-MM-DDTHH:MMZ")


def fmt_time(t):
    return t.strftime("%Y-%m-%dT%H:%MZ")


# ---------------------------------------------------------------------------- diagnostics

def virtual_flux(s):
    """w'theta_v' [K m/s] from sensible and latent heat flux [W/m2, upward positive]."""
    q = s["q2m"]
    tv = s["t2m"] * (1 + 0.61 * q)
    rho = s["psfc"] / (RD * tv)
    theta = s["t2m"] * (1e5 / s["psfc"]) ** KAPPA_P
    return s["shf"] / (rho * CP) * (1 + 0.61 * q) + 0.61 * theta * s["lhf"] / (rho * LV)


def obukhov_length(s, qv):
    if qv == 0:
        return math.inf
    thv = s["t2m"] * (1e5 / s["psfc"]) ** KAPPA_P * (1 + 0.61 * s["q2m"])
    return -s["ustar"] ** 3 * thv / (KAPPA_VK * G * qv)


def _interp(z, a, zz, name):
    if zz < z[0] or zz > z[-1]:
        raise RegimeError(f"column: height {zz:.0f} m for {name} is outside the column ({z[0]:.0f} to {z[-1]:.0f} m)")
    return float(np.interp(zz, z, a))


def _layer_mean(z, a, z0, z1):
    zz = np.linspace(z0, z1, 201)
    f = np.interp(zz, z, a)
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(zz)) / (z1 - z0))


def bulk_richardson(c, ztop):
    z, th = np.asarray(c["z"]), np.asarray(c["theta_v"])
    u, v = np.asarray(c["u"]), np.asarray(c["v"])
    if ztop <= z[0]:
        raise RegimeError(f"thresholds: 'rib_top' = {ztop} m is not above the lowest level ({z[0]:.0f} m)")
    dth = _interp(z, th, ztop, "rib_top") - th[0]
    du = _interp(z, u, ztop, "rib_top") - u[0]
    dv = _interp(z, v, ztop, "rib_top") - v[0]
    thm = 0.5 * (_interp(z, th, ztop, "rib_top") + th[0])
    sh2 = du * du + dv * dv
    num = G / thm * dth * (ztop - z[0])
    if sh2 == 0:
        return math.copysign(math.inf, num) if num != 0 else 0.0
    return num / sh2


def buoyancy_frequency(c, z0, z1, name):
    z, th = np.asarray(c["z"]), np.asarray(c["theta_v"])
    t0, t1 = _interp(z, th, z0, name), _interp(z, th, z1, name)
    return G / (0.5 * (t0 + t1)) * (t1 - t0) / (z1 - z0)   # N^2 [1/s^2]


def jet(c, t):
    """Nose height and speed and the fall-off above it; None when the profile has no jet."""
    z = np.asarray(c["z"])
    s = np.hypot(c["u"], c["v"])
    below = np.nonzero(z <= t["jet_nose_max"])[0]
    if below.size < 2:
        return None
    k = int(below[np.argmax(s[below])])
    if k == 0:                                     # largest at the lowest level: no nose
        return None
    zmax = z[k] + t["jet_depth"]
    if zmax > z[-1]:
        raise RegimeError(f"column: the jet search (nose + 'jet_depth') reaches {zmax:.0f} m, above the column top "
                          f"{z[-1]:.0f} m")
    zz = np.concatenate([z[(z > z[k]) & (z < zmax)], [zmax]])
    smin = float(np.min(np.interp(zz, z, s)))
    drop = float(s[k] - smin)
    found = drop >= t["jet_drop"] and drop >= t["jet_drop_frac"] * s[k]
    return {"nose_height": float(z[k]), "nose_speed": float(s[k]), "drop": drop, "jet": bool(found)}


def geostrophic_shear(c, z0, depth):
    """Magnitude of the layer-mean thermal-wind shear vector (g / (f theta_v)) k x grad theta_v over z0..z0+depth,
    in m/s per km. Above z_i, where the temperature field is close to balance."""
    z = np.asarray(c["z"])
    f = 2 * OMEGA * math.sin(math.radians(c["lat"]))
    th = np.asarray(c["theta_v"])
    z0 = max(z0, z[0])
    _interp(z, th, z0 + depth, "thermal_depth above z_i")
    sx = _layer_mean(z, -G / (f * th) * np.asarray(c["dthdy"]), z0, z0 + depth)
    sy = _layer_mean(z, G / (f * th) * np.asarray(c["dthdx"]), z0, z0 + depth)
    return 1e3 * math.hypot(sx, sy)


def layer_wind(c, ztop):
    z = np.asarray(c["z"])
    _interp(z, z, ztop, "tendency_top")
    zz = np.concatenate([[0.0], z])        # the wind is 0 at the ground
    return (_layer_mean(zz, np.concatenate([[0.0], c["u"]]), 0.0, ztop),
            _layer_mean(zz, np.concatenate([[0.0], c["v"]]), 0.0, ztop))


def _n2_above_zi(c, zi, depth):
    """N^2 over max(z_i, z_1)..that + n_depth, a diagnostic only (no class uses it): NaN above the column top."""
    z0 = max(zi, c["z"][0])
    if z0 + depth > c["z"][-1]:
        return math.nan
    return buoyancy_frequency(c, z0, z0 + depth, "n_depth above z_i")


def diagnose(c, t, relief=None):
    s = c["surface"]
    qv = virtual_flux(s)
    L = obukhov_length(s, qv)
    zi = s["hpbl"]
    z = np.asarray(c["z"])
    d = {"q_v": qv, "L": L, "z_i": zi, "zi_over_L": zi / L if math.isfinite(L) else 0.0,
         "rib_low": bulk_richardson(c, t["rib_top"]),
         "N2_above_zi": _n2_above_zi(c, zi, t["n_depth"]),
         "geo_shear": geostrophic_shear(c, zi, t["thermal_depth"])}
    ctop = z <= t["cloud_top"]
    d["cloud_max"] = float(np.max(np.asarray(c["cloud_max"])[ctop])) if ctop.any() else 0.0
    d["jet"] = jet(c, t)
    H = s["relief"] if relief is None else float(relief)
    if not math.isfinite(H) or H < 0:
        raise RegimeError(f"relief must be finite and >= 0, got {H}")
    d["relief"] = H
    if H > 0:
        top = max(H, t["froude_min_depth"], z[1])
        n2 = buoyancy_frequency(c, z[0], top, "Fr (relief or 'froude_min_depth')")
        U = _layer_mean(z, np.hypot(c["u"], c["v"]), z[0], top)
        d["froude"] = U / (math.sqrt(n2) * H) if n2 > 0 else math.inf
    else:
        d["froude"] = math.inf
    return d


# ---------------------------------------------------------------------------- the class

def classify(c, t=None, before=None, after=None, relief=None):
    t = default_thresholds() if t is None else t
    check_column(c, "column")
    now = parse_time(c["valid_time"])
    hour = dt.timedelta(hours=1)
    for nb, name, delta in ((before, "before", -hour), (after, "after", hour)):
        if nb is not None:
            check_column(nb, name)
            if parse_time(nb["valid_time"]) - now != delta:
                raise RegimeError(f"{name}: valid time {nb['valid_time']} is not one hour "
                                  f"{name} {c['valid_time']}")
    d = diagnose(c, t, relief)
    eps = t["flux_eps"]
    qv, L = d["q_v"], d["L"]
    notes = []

    trans = None
    if before is None and after is None:
        notes.append("no neighbouring hours: transitions and T not tested")
    else:
        pairs = [(before, c)] if before is not None else []
        pairs += [(c, after)] if after is not None else []
        for a, b in pairs:
            qa, qb = virtual_flux(a["surface"]), virtual_flux(b["surface"])
            if max(abs(qa), abs(qb)) <= eps:      # a sign change of a flux that is near zero on both sides
                continue
            if qa > 0 >= qb:
                trans = "C7"
            elif qa <= 0 < qb:
                trans = trans or "C8"
        if before is None or after is None:
            notes.append("one neighbouring hour only: transitions tested on one side")

    j = d["jet"]
    if d["cloud_max"] > t["cloud_mixing_ratio"]:
        cls, rule = "C9", "cloud water in the box below the cloud top"
    elif trans:
        cls, rule = trans, "q_v changes sign within one hour, |q_v| > flux_eps on one side"
    elif qv > eps and -L <= t["neutral_length"]:
        cls = "C1" if -d["zi_over_L"] > t["convective_ratio"] else "C2"
        rule = "unstable, -z_i/L " + (">" if cls == "C1" else "<=") + " convective_ratio"
    elif qv < 0 and d["rib_low"] > 0 and j is not None and j["jet"]:
        cls, rule = "C5", "q_v < 0, rib_low > 0 and a low-level jet"
    elif abs(qv) <= eps:                          # L is ill-conditioned here: the layer Ri_b decides
        cls, rule = ("C6", "|q_v| <= flux_eps and rib_low > rib_crit") if d["rib_low"] > t["rib_crit"] else \
                    ("C3", "|q_v| <= flux_eps and rib_low <= rib_crit")
    elif abs(L) > t["neutral_length"]:
        cls, rule = "C3", "|L| > neutral_length"
    elif d["zi_over_L"] > t["very_stable_ratio"] or d["rib_low"] > t["rib_crit"]:
        cls, rule = "C6", "z_i/L > very_stable_ratio or rib_low > rib_crit"
    else:
        cls, rule = "C4", "stable, z_i/L <= very_stable_ratio and rib_low <= rib_crit"

    mods = []
    if d["froude"] < t["froude"]:
        mods.append("Fr")
    if d["geo_shear"] > t["geo_shear"]:
        mods.append("B")
    tend = None
    if before is not None or after is not None:
        a = before if before is not None else c
        b = after if after is not None else c
        hours = (parse_time(b["valid_time"]) - parse_time(a["valid_time"])).total_seconds() / 3600
        ua, va = layer_wind(a, t["tendency_top"])
        ub, vb = layer_wind(b, t["tendency_top"])
        tend = math.hypot(ub - ua, vb - va) / hours
        if tend > t["tendency"]:
            mods.append("T")
    d["tendency"] = tend
    return {"valid_time": c["valid_time"], "class": cls, "modifiers": mods, "rule": rule, "notes": notes,
            "diagnostics": d}


# ---------------------------------------------------------------------------- catalogue

def _classify_job(args):
    path, before, after, t = args
    load = lambda p: None if p is None else read_column(p)
    try:
        r = classify(read_column(path), t, load(before), load(after))
        r["file"] = path
        return r
    except (RegimeError, ArithmeticError, ValueError) as e:
        return {"file": path, "error": str(e) or type(e).__name__}


def read_column(path, what="column"):
    try:
        c = json.loads(_read_text(path, what))
    except ValueError as e:
        raise RegimeError(f"{what} {path}: not JSON ({e})")
    if not isinstance(c, dict):
        raise RegimeError(f"{what} {path}: not a column (a JSON object is needed)")
    check_column(c, path)
    return c


def catalogue(paths, t, jobs=1):
    cols = {}
    for p in paths:
        c = read_column(p)
        k = parse_time(c["valid_time"])
        if k in cols:
            raise RegimeError(f"{p}: valid time {c['valid_time']} also in {cols[k]}")
        cols[k] = p
    hour = dt.timedelta(hours=1)
    work = [(cols[k], cols.get(k - hour), cols.get(k + hour), t) for k in sorted(cols)]
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(jobs) as pool:
            return pool.map(_classify_job, work)
    return [_classify_job(w) for w in work]


def _fmt(x, n=3):
    if x is None:
        return ""
    if isinstance(x, float) and not math.isfinite(x):
        return "nan" if math.isnan(x) else ("inf" if x > 0 else "-inf")
    return f"{x:.{n}g}"


CSV_FIELDS = ("valid_time", "class", "modifiers", "rule", "q_v", "L", "z_i", "zi_over_L", "rib_low", "N2_above_zi",
              "jet", "jet_nose_height", "jet_nose_speed", "jet_drop", "geo_shear", "tendency", "froude", "cloud_max",
              "file", "error")


def catalogue_rows(results):
    rows = []
    for r in results:
        if "error" in r:
            rows.append({"file": r["file"], "error": r["error"]})
            continue
        d = r["diagnostics"]
        j = d["jet"] or {}
        rows.append({"valid_time": r["valid_time"], "class": r["class"], "modifiers": "+".join(r["modifiers"]),
                     "rule": r["rule"], "jet": str(int(bool(j.get("jet")))),
                     "q_v": _fmt(d["q_v"]), "L": _fmt(d["L"]), "z_i": _fmt(d["z_i"]),
                     "zi_over_L": _fmt(d["zi_over_L"]), "rib_low": _fmt(d["rib_low"]),
                     "N2_above_zi": _fmt(d["N2_above_zi"]), "jet_nose_height": _fmt(j.get("nose_height")),
                     "jet_nose_speed": _fmt(j.get("nose_speed")), "jet_drop": _fmt(j.get("drop")),
                     "geo_shear": _fmt(d["geo_shear"]), "tendency": _fmt(d["tendency"]), "froude": _fmt(d["froude"]),
                     "cloud_max": _fmt(d["cloud_max"]), "file": r["file"], "error": ""})
    return rows


def count_table(results):
    """Hours per class and month, as text."""
    months = sorted({r["valid_time"][:7] for r in results if "error" not in r})
    lines = ["class  " + " ".join(f"{m:>8s}" for m in months) + "    total"]
    for cl in CLASSES:
        n = [sum(1 for r in results if "error" not in r and r["class"] == cl and r["valid_time"][:7] == m)
             for m in months]
        lines.append(f"{cl:6s} " + " ".join(f"{x:8d}" for x in n) + f" {sum(n):8d}")
    bad = sum(1 for r in results if "error" in r)
    if bad:
        lines.append(f"refused: {bad}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------- fetch (herbie, xarray, cfgrib)

NAT_REGEX = r":(UGRD|VGRD|HGT|TMP|SPFH|PRES|CLMR|CLWMR):([1-9]|1[0-9]|20) hybrid level"
SFC_REGEX = r":(SHTFL|LHTFL|FRICV|HPBL|SFCR):|:HGT:surface|:(TMP|SPFH):2 m|:PRES:surface"


def _fetch_one(args):
    valid, lat0, lon0, half_km, grad_half_km, fxx, save, out_dir = args
    import warnings
    warnings.filterwarnings("ignore")
    from herbie import Herbie
    vt = parse_time(valid)
    init = vt - dt.timedelta(hours=fxx)
    when = init.strftime("%Y-%m-%d %H:%M")
    res, sl = {}, None
    for product, regex, prefix in (("nat", NAT_REGEX, "nat_"), ("sfc", SFC_REGEX, "sfc_")):
        dss = Herbie(when, model="hrrr", product=product, fxx=fxx, save_dir=save, verbose=False).xarray(
            regex, remove_grib=False)
        if not isinstance(dss, list):
            dss = [dss]
        for ds in dss:
            if sl is None:
                lat = ds.latitude.values
                lon = np.where(ds.longitude.values > 180, ds.longitude.values - 360, ds.longitude.values)
                d2 = (lat - lat0) ** 2 + ((lon - _wrap_lon(lon0)) * math.cos(math.radians(lat0))) ** 2
                j, i = np.unravel_index(np.argmin(d2), d2.shape)
                if math.sqrt(d2[j, i]) * math.radians(1) * R_EARTH > 3e3:
                    raise RegimeError(f"{valid}: ({lat0}, {lon0}) is outside the HRRR grid")
                # the gradient box may be rotated against the grid by up to ~25 degrees: room for its corners
                n = int(math.ceil(math.sqrt(2) * grad_half_km / 3.0)) + 3
                if j - n < 0 or i - n < 0 or j + n + 1 > lat.shape[0] or i + n + 1 > lat.shape[1]:
                    raise RegimeError(f"{valid}: the {grad_half_km} km gradient box ('--grad-half-km') passes the "
                                      "edge of the HRRR domain")
                sl = (slice(j - n, j + n + 1), slice(i - n, i + n + 1))
                res["lat"], res["lon"] = lat[sl], lon[sl]
            for name in ds.data_vars:
                a = ds[name].values
                res[prefix + name] = a[(slice(None),) + sl] if a.ndim == 3 else a[sl]
    col = column_from_crop(res, lat0, lon0, half_km, fmt_time(vt), fxx,
                           source=f"HRRR {init:%Y-%m-%d %H} UTC f{fxx:02d} (wrfnat, wrfsfc), NOAA",
                           grad_half_km=grad_half_km)
    out = os.path.join(out_dir, f"hrrr_{vt:%Y%m%d%H}_f{fxx:02d}.json")
    write_json(col, out)
    return out


def write_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")


# ---------------------------------------------------------------------------- command line

def _json_safe(o):
    if isinstance(o, float) and not math.isfinite(o):
        return "nan" if math.isnan(o) else ("inf" if o > 0 else "-inf")
    if isinstance(o, dict):
        return {k: _json_safe(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_json_safe(v) for v in o]
    return o


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download HRRR hours and write their columns")
    f.add_argument("valid", nargs="+", help="valid times, YYYY-MM-DDTHH:MMZ")
    def box(p):
        p.add_argument("--lat", type=float, required=True, help="latitude of the column, 21 to 53")
        p.add_argument("--lon", type=float, required=True, help="longitude, -135 to -60 (or 225 to 300)")
        p.add_argument("--half-km", type=float, default=12.5, help="half-width of the mean square, km (12.5)")
        p.add_argument("--grad-half-km", type=float, default=40.0,
                       help="half-width of the theta_v gradient fit, km, at least --half-km (40)")
        p.add_argument("--fxx", type=int, default=1, help="forecast hour, at least 1 (1)")
    box(f)
    f.add_argument("--save", required=True, help="Herbie download folder")
    f.add_argument("--out-dir", required=True, help="folder for the columns, hrrr_<YYYYMMDDHH>_f<fxx>.json")
    f.add_argument("--jobs", type=int, default=1, help="hours downloaded at once (1)")
    c = sub.add_parser("column", help="column from a crop saved earlier (.npz)")
    c.add_argument("npz", help="crop with nat_* [k, j, i], sfc_* [j, i], lat and lon")
    c.add_argument("--valid", required=True, help="valid time, YYYY-MM-DDTHH:MMZ")
    box(c)
    c.add_argument("--out", required=True, help="column JSON to write")
    k = sub.add_parser("classify", help="class and diagnostics of one hour")
    k.add_argument("column", help="column JSON")
    k.add_argument("--before", help="column of the hour before (transitions, tendency)")
    k.add_argument("--after", help="column of the hour after (transitions, tendency)")
    k.add_argument("--relief", type=float, help="relief H for the Froude number, m (default: the column's)")
    k.add_argument("--thresholds", help="JSON or YAML file replacing any default threshold")
    g = sub.add_parser("catalogue", help="classes of many hours, neighbours found by valid time")
    g.add_argument("columns", nargs="+", help="column JSON files, one valid time each")
    g.add_argument("--jobs", type=int, default=1, help="hours classified at once (1); same result for any value")
    g.add_argument("--thresholds", help="JSON or YAML file replacing any default threshold")
    g.add_argument("--out", help="CSV file, one row per hour")
    sub.add_parser("thresholds", help="print the default thresholds as JSON")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "thresholds":
            print(json.dumps(default_thresholds(), indent=1))
        elif a.cmd == "fetch":
            if a.jobs < 1:
                raise RegimeError("--jobs must be >= 1")
            if a.fxx < 1:
                raise RegimeError("--fxx must be >= 1: the f00 analysis carries almost no vertical motion")
            check_box(a.lat, a.lon, a.half_km, a.grad_half_km)
            for v in a.valid:
                parse_time(v, "valid")
            os.makedirs(a.out_dir, exist_ok=True)
            work = [(v, a.lat, a.lon, a.half_km, a.grad_half_km, a.fxx, a.save, a.out_dir) for v in a.valid]
            if a.jobs > 1:
                import multiprocessing as mp
                with mp.get_context("spawn").Pool(a.jobs) as pool:
                    outs = pool.map(_fetch_one, work)
            else:
                outs = [_fetch_one(w) for w in work]
            print("\n".join(outs))
        elif a.cmd == "column":
            check_box(a.lat, a.lon, a.half_km, a.grad_half_km)
            valid = fmt_time(parse_time(a.valid, "--valid"))
            if a.fxx < 1:
                raise RegimeError("--fxx must be >= 1: the f00 analysis carries almost no vertical motion")
            try:
                d = dict(np.load(a.npz))
            except (OSError, ValueError) as e:
                raise RegimeError(f"crop {a.npz}: cannot be read as .npz ({e})")
            col = column_from_crop(d, a.lat, a.lon, a.half_km, valid, a.fxx,
                                   source=f"HRRR f{a.fxx:02d} crop {os.path.basename(a.npz)}, NOAA",
                                   grad_half_km=a.grad_half_km)
            write_json(col, a.out)
        elif a.cmd == "classify":
            t = read_thresholds(a.thresholds)
            r = classify(read_column(a.column), t, a.before and read_column(a.before, "--before"),
                         a.after and read_column(a.after, "--after"), a.relief)
            print(json.dumps(_json_safe(r), indent=1))
        elif a.cmd == "catalogue":
            if a.jobs < 1:
                raise RegimeError("--jobs must be >= 1")
            t = read_thresholds(a.thresholds)
            res = catalogue(a.columns, t, a.jobs)
            rows = catalogue_rows(res)
            if a.out:
                with open(a.out, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
                    w.writeheader()
                    w.writerows(rows)
            for r in rows:
                if r["error"]:
                    print(f"{r['file']}: refused: {r['error']}")
                else:
                    print(f"{r['valid_time']}  {r['class']}  {r['modifiers'] or '-':6s} L {r['L']:>8s}  "
                          f"zi/L {r['zi_over_L']:>7s}  Rib {r['rib_low']:>7s}  "
                          f"jet {r['jet_drop'] if r['jet'] == '1' else '-':>5s}")
            print(count_table(res))
            return 1 if any("error" in r for r in res) else 0
    except RegimeError as e:
        print(f"hrrr_regime.py: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
