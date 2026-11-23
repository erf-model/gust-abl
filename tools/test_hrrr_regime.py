"""Unit tests for hrrr_regime.py (numpy only): python3 -m unittest -v test_hrrr_regime"""
import contextlib
import copy
import io
import json
import math
import os
import tempfile
import unittest

import numpy as np

import hrrr_regime as hr

Z = [10.0, 40, 90, 175, 300, 450, 650, 850, 1100, 1400, 1700, 2050, 2450, 2900, 3400, 4000, 4700, 5500]


def column(valid="2023-05-10T05:00Z", shf=0.0, lhf=0.0, ustar=0.4, hpbl=500.0, lapse=0.0, cap=None, wind=None,
           cloud=None, grad=(0.0, 0.0), lat=36.0, relief=0.0, t2m=290.0, q2m=0.0, psfc=1e5):
    """A column on Z: theta_v = 300 + lapse z (K/m) with an optional (height, jump) inversion, a log wind of 10 m/s
    at 1 km unless wind (speed per level) is given, wind from the west."""
    z = np.array(Z)
    th = 300.0 + lapse * z
    if cap is not None:
        th = th + cap[1] * (z > cap[0])
    s = np.array(wind, float) if wind is not None else 10.0 * np.log(z / 0.1) / math.log(1000 / 0.1)
    return {"valid_time": valid, "fxx": 1, "lat": lat, "lon": -97.4, "half_km": 12.5, "cells": 81, "source": "test",
            "z": list(z), "u": list(s), "v": [0.0] * len(Z), "theta_v": list(th),
            "cloud_max": list(np.zeros(len(Z)) if cloud is None else cloud),
            "dthdx": [grad[0]] * len(Z), "dthdy": [grad[1]] * len(Z),
            "surface": {"shf": shf, "lhf": lhf, "ustar": ustar, "z0": 0.1, "hpbl": hpbl, "t2m": t2m, "q2m": q2m,
                        "psfc": psfc, "relief": relief}}


def jet_wind(nose=450.0, peak=15.0, aloft=10.0):
    """Speed rising to peak at the nose and falling to aloft 500 m above it."""
    z = np.array(Z)
    low = peak * np.log(z / 0.1) / math.log(nose / 0.1)
    high = peak - (peak - aloft) * np.clip((z - nose) / 500.0, 0, 1)
    return list(np.where(z <= nose, low, high))


class SurfaceFlux(unittest.TestCase):
    def test_obukhov_closed_form(self):
        # sensible heat only, q = 0, p = 1000 hPa so theta = T: rho = p / (Rd T), q_v = H / (rho cp)
        s = column(shf=100.0, ustar=0.3, t2m=300.0)["surface"]
        rho = 1e5 / (287.04 * 300.0)
        qv = 100.0 / (rho * 1005.0)
        self.assertAlmostEqual(hr.virtual_flux(s), qv, places=12)
        L = -0.3 ** 3 * 300.0 / (0.41 * 9.81 * qv)
        self.assertAlmostEqual(hr.obukhov_length(s, hr.virtual_flux(s)), L, places=9)
        self.assertAlmostEqual(L, -23.504, places=3)        # number worked by hand

    def test_latent_heat_counts(self):
        # no sensible heat, 200 W/m2 of evaporation: q_v = 0.61 theta E / (rho Lv) > 0 (a kinematic flux gives 0)
        s = column(shf=0.0, lhf=200.0, t2m=300.0, q2m=0.01)["surface"]
        tv = 300.0 * 1.0061
        rho = 1e5 / (287.04 * tv)
        self.assertAlmostEqual(hr.virtual_flux(s), 0.61 * 300.0 * 200.0 / (rho * 2.5e6), places=12)
        self.assertGreater(hr.virtual_flux(s), 0.01)

    def test_sign(self):
        # upward heat flux (HRRR sign) is unstable: L < 0; downward is stable: L > 0
        self.assertLess(hr.diagnose(column(shf=150.0), hr.default_thresholds())["L"], 0)
        self.assertGreater(hr.diagnose(column(shf=-20.0), hr.default_thresholds())["L"], 0)


class Profiles(unittest.TestCase):
    def test_bulk_richardson_closed_form(self):
        # theta_v linear 0.01 K/m, u log: Ri_b = g/thm dtheta dz / du^2 between z[0] and 100 m
        c = column(lapse=0.01)
        t = hr.default_thresholds()
        u = lambda zz: 10 * math.log(zz / 0.1) / math.log(1e4)
        u100 = np.interp(100.0, Z, c["u"])
        dth = 0.01 * 90.0
        ri = 9.81 / (300.0 + 0.01 * 55.0) * dth * 90.0 / (u100 - u(10.0)) ** 2
        self.assertAlmostEqual(hr.bulk_richardson(c, t["rib_top"]), ri, places=10)

    def test_buoyancy_frequency(self):
        c = column(lapse=0.003)
        n2 = hr.buoyancy_frequency(c, 500.0, 1000.0, "test")
        self.assertAlmostEqual(n2, 9.81 / (300 + 0.003 * 750) * 0.003, places=12)

    def test_jet_found(self):
        t = hr.default_thresholds()
        j = hr.jet(column(wind=jet_wind()), t)
        self.assertTrue(j["jet"])
        self.assertEqual(j["nose_height"], 450.0)
        self.assertAlmostEqual(j["nose_speed"], 15.0)
        self.assertAlmostEqual(j["drop"], 5.0)

    def test_no_jet_in_log_profile(self):
        j = hr.jet(column(), hr.default_thresholds())       # speed still rising above the highest nose allowed
        self.assertLessEqual(j["drop"], 0.0)
        self.assertFalse(j["jet"])

    def test_drop_measured_above_the_nose(self):
        # a nose of 6 m/s over a 2 m/s surface and 5 m/s aloft: rise of 4, drop of 1, no jet
        z = np.array(Z)
        w = np.where(z <= 450, 2 + 4 * (z - 10) / 440, 6 - np.clip((z - 450) / 500, 0, 1))
        j = hr.jet(column(wind=list(w)), hr.default_thresholds())
        self.assertAlmostEqual(j["drop"], 1.0)
        self.assertFalse(j["jet"])

    def test_jet_needs_both_drops(self):
        t = dict(hr.default_thresholds(), jet_drop_frac=0.5)
        self.assertFalse(hr.jet(column(wind=jet_wind()), t)["jet"])     # 5 of 15 m/s is 33 % < 50 %

    def test_thermal_wind_closed_form(self):
        # |dVg/dz| = g / (f theta) |grad theta|, 1 K per 100 km at 36 N; theta_v constant
        c = column(grad=(1e-5, 0.0), hpbl=500.0)
        f = 2 * 7.292e-5 * math.sin(math.radians(36.0))
        self.assertAlmostEqual(hr.geostrophic_shear(c, 500.0, 1000.0), 1e3 * 9.81 / (f * 300.0) * 1e-5, places=9)

    def test_thermal_wind_is_a_vector_mean(self):
        # the gradient turns round over the layer: the mean shear vector cancels, the mean magnitude does not
        c = column(hpbl=500.0)
        z = np.array(Z)
        c["dthdx"] = list(np.where(z < 1000, 1e-5, -1e-5))
        self.assertLess(hr.geostrophic_shear(c, 500.0, 1000.0), 0.3)
        self.assertGreater(hr.geostrophic_shear(column(grad=(1e-5, 0)), 500.0, 1000.0), 2.0)

    def test_froude(self):
        # U / (N H): theta_v 0.01 K/m, wind 5 m/s, relief 200 m
        c = column(lapse=0.01, wind=[5.0] * len(Z), relief=200.0)
        d = hr.diagnose(c, hr.default_thresholds())
        n2 = 9.81 / (300 + 0.01 * 105) * 0.01
        self.assertAlmostEqual(d["froude"], 5.0 / (math.sqrt(n2) * 200.0), places=9)
        self.assertEqual(hr.diagnose(column(), hr.default_thresholds())["froude"], math.inf)

    def test_rotation(self):
        # angle = sin(38.5 deg) (lon + 97.5): none at -97.5, 17.12 deg at -70
        u, v = hr.rotate_to_earth(1.0, 0.0, -97.5)
        self.assertAlmostEqual(float(u), 1.0)
        self.assertAlmostEqual(float(v), 0.0)
        u, v = hr.rotate_to_earth(1.0, 0.0, -70.0)
        a = math.radians(math.sin(math.radians(38.5)) * 27.5)
        self.assertAlmostEqual(float(u), math.cos(a))
        self.assertAlmostEqual(float(v), -math.sin(a))
        self.assertAlmostEqual(math.degrees(a), 17.12, places=2)
        # 40 N, 120 W: angle -14.007 deg; checked against a Lambert conformal projection of the meridian
        u, v = hr.rotate_to_earth(1.0, 0.0, -120.0)
        self.assertAlmostEqual(float(u), 0.9703, places=4)
        self.assertAlmostEqual(float(v), 0.2420, places=4)
        u, v = hr.rotate_to_earth(0.0, 1.0, -120.0)
        self.assertAlmostEqual(float(u), -0.2420, places=4)
        self.assertAlmostEqual(float(v), 0.9703, places=4)

    def test_nan_printed_as_nan(self):
        self.assertEqual(hr._fmt(float("nan")), "nan")
        self.assertEqual(hr._json_safe({"a": [float("nan"), float("-inf")]}), {"a": ["nan", "-inf"]})


def crop(n=37, dx=3e3, nk=16, slope=0.0, gx=0.0, lapse=0.004, rot=0.0):
    """A crop around 36 N, 97.4 W on a uniform grid turned rot degrees from east: hybrid levels of height
    z_k + slope x, theta = 300 + lapse z + gx x on constant height, so the gradient at constant height is (gx, 0)
    whatever the level slope."""
    lat0, lon0 = 36.0, -97.4
    i = np.arange(n) - n // 2
    xi, yj = np.meshgrid(i * dx, i * dx)
    r = math.radians(rot)
    x, y = xi * math.cos(r) - yj * math.sin(r), xi * math.sin(r) + yj * math.cos(r)
    lat = lat0 + np.degrees(y / hr.R_EARTH)
    lon = lon0 + np.degrees(x / (hr.R_EARTH * math.cos(math.radians(lat0))))
    zk = 20.0 * 1.25 ** np.arange(nk)
    gh = zk[:, None, None] + slope * x[None]
    theta = 300 + lapse * gh + gx * x[None]
    p = 1e5 * np.exp(-gh / 8000.0)
    T = theta / (1e5 / p) ** hr.KAPPA_P
    zero = np.zeros_like(gh)
    a = np.radians(math.sin(math.radians(38.5)) * (lon + 97.5))
    # earth wind 8 m/s from the west, written grid-relative
    return {"lat": lat, "lon": lon, "nat_pres": p, "nat_t": T, "nat_q": zero, "nat_u": 8 * np.cos(a) + zero,
            "nat_v": 8 * np.sin(a) + zero, "nat_gh": gh, "nat_clwmr": zero,
            "sfc_orog": np.zeros((n, n)), "sfc_ishf": np.full((n, n), 50.0), "sfc_slhtf": np.zeros((n, n)),
            "sfc_fricv": np.full((n, n), 0.4), "sfc_fsr": np.full((n, n), 0.1), "sfc_blh": np.full((n, n), 800.0),
            "sfc_t2m": np.full((n, n), 295.0), "sfc_sh2": np.zeros((n, n)), "sfc_sp": np.full((n, n), 1e5)}


class Column(unittest.TestCase):
    def test_gradient_at_constant_height(self):
        # levels slope 2 m per km: the level gradient of theta is gx + lapse * slope, the column keeps gx only
        c = hr.column_from_crop(crop(slope=2e-3, gx=1e-5), 36.0, -97.4, 12.5, "2023-05-10T05:00Z", 1)
        k = 6
        self.assertAlmostEqual(c["dthdx"][k], 1e-5, delta=2e-7)
        self.assertAlmostEqual(c["dthdy"][k], 0.0, delta=2e-7)

    def test_earth_relative_wind(self):
        c = hr.column_from_crop(crop(), 36.0, -97.4, 12.5, "2023-05-10T05:00Z", 1)
        self.assertAlmostEqual(c["u"][3], 8.0, places=9)
        self.assertAlmostEqual(c["v"][3], 0.0, places=9)

    def test_box(self):
        c = hr.column_from_crop(crop(), 36.0, -97.4, 12.5, "2023-05-10T05:00Z", 1)
        self.assertEqual(c["cells"], 81)           # 9 x 9 cells of 3 km within 12.5 km
        with self.assertRaisesRegex(hr.RegimeError, "reaches the edge of the crop"):
            hr.column_from_crop(crop(), 36.0, -97.4, 12.5, "2023-05-10T05:00Z", 1, grad_half_km=60.0)
        for bad, msg in ((dict(grad_half_km=float("nan")), "--grad-half-km"), (dict(grad_half_km=10.0), "--grad"),
                         (dict(lat0=60.0), "--lat"), (dict(lon0=10.0), "--lon"), (dict(half_km=0.0), "--half-km")):
            kw = dict(lat0=36.0, lon0=-97.4, half_km=12.5, grad_half_km=40.0)
            kw.update(bad)
            with self.assertRaisesRegex(hr.RegimeError, msg):
                hr.column_from_crop(crop(), kw["lat0"], kw["lon0"], kw["half_km"], "2023-05-10T05:00Z", 1,
                                    grad_half_km=kw["grad_half_km"])

    def test_rotated_crop(self):
        # a 51 km crop turned 15.6 degrees (HRRR near 122 W) holds the 40 km box: its corners reach 49.5 km
        c = hr.column_from_crop(crop(n=35, rot=15.6, gx=1e-5), 36.0, -97.4, 12.5, "2023-05-10T05:00Z", 1)
        self.assertAlmostEqual(c["dthdx"][6], 1e-5, delta=2e-7)
        c2 = hr.column_from_crop(crop(), 36.0, 262.6, 12.5, "2023-05-10T05:00Z", 1)     # longitude 0..360
        self.assertEqual(c2["cells"], 81)


class Classes(unittest.TestCase):
    T = hr.default_thresholds()

    def cls(self, c, **kw):
        return hr.classify(c, self.T, **kw)

    def test_c1_convective(self):
        r = self.cls(column(shf=250.0, ustar=0.4, hpbl=1500.0, lapse=0.003))
        self.assertEqual(r["class"], "C1")
        self.assertGreater(-r["diagnostics"]["zi_over_L"], 10)

    def test_c2_weakly_unstable(self):
        r = self.cls(column(shf=30.0, ustar=0.4, hpbl=600.0, lapse=0.003))
        self.assertEqual(r["class"], "C2")

    def test_c3_neutral(self):
        r = self.cls(column(shf=-10.0, ustar=0.8, cap=(800, 3), lapse=0.003))
        self.assertEqual((r["class"], r["rule"]), ("C3", "|L| > neutral_length"))
        self.assertEqual(self.cls(column(shf=0.0, ustar=0.6, lapse=0.003))["class"], "C3")

    def test_c4_weakly_stable(self):
        r = self.cls(column(shf=-15.0, ustar=0.35, hpbl=200.0, lapse=0.003))
        d = r["diagnostics"]
        self.assertTrue(0 < d["zi_over_L"] <= 1 and d["rib_low"] <= 0.25, d)
        self.assertEqual(r["class"], "C4")

    def test_c5_jet(self):
        # cooling and a jet: C5 even though |L| is large
        r = self.cls(column(shf=-20.0, ustar=0.6, hpbl=300.0, lapse=0.003, wind=jet_wind()))
        self.assertGreater(r["diagnostics"]["L"], 500)
        self.assertEqual(r["class"], "C5")

    def test_c6_very_stable(self):
        r = self.cls(column(shf=-10.0, ustar=0.12, hpbl=100.0, lapse=0.02))
        self.assertEqual(r["class"], "C6")

    def test_c5_weak_cooling(self):
        # a jet over a weakly cooled surface (|q_v| < flux_eps, L of kilometres, Ri_b > 0) is still C5
        r = self.cls(column(shf=-2.0, ustar=0.5, hpbl=300.0, lapse=0.003, wind=jet_wind()))
        self.assertLess(abs(r["diagnostics"]["q_v"]), self.T["flux_eps"])
        self.assertGreater(r["diagnostics"]["L"], 500)
        self.assertEqual(r["class"], "C5")
        r = self.cls(column(shf=2.0, ustar=0.5, hpbl=300.0, lapse=0.0, wind=jet_wind()))     # heating: not C5
        self.assertNotEqual(r["class"], "C5")

    def test_weak_flux_uses_richardson(self):
        # |q_v| <= flux_eps: L is not used. 1 W/m2 under u* 0.05 has -z_i/L of about 150 but is not C1
        r = self.cls(column(shf=1.0, ustar=0.05, hpbl=1500.0, lapse=0.003))
        self.assertGreater(-r["diagnostics"]["zi_over_L"], 10)
        self.assertEqual(r["class"], "C3")
        c = column(shf=-3.0, ustar=0.6, hpbl=200.0, lapse=0.003)
        z = np.array(Z)
        c["theta_v"] = list(300 + 0.003 * z + 3.0 * (z >= 90))
        r = self.cls(c)
        self.assertGreater(r["diagnostics"]["L"], 500)
        self.assertEqual((r["class"], r["rule"]), ("C6", "|q_v| <= flux_eps and rib_low > rib_crit"))

    def test_low_pbl_height(self):
        # z_i below the lowest level: the hour is classified; N^2 starts at the lowest level
        d = self.cls(column(shf=-15.0, ustar=0.35, hpbl=5.0, lapse=0.003))["diagnostics"]
        self.assertAlmostEqual(d["N2_above_zi"], 9.81 / (300 + 0.003 * 260) * 0.003, places=9)

    def test_c6_by_richardson(self):
        c = column(shf=-15.0, ustar=0.35, hpbl=200.0, lapse=0.003)
        z = np.array(Z)
        c["theta_v"] = list(300 + 0.003 * z + 3.0 * (z >= 90))       # 3 K inversion at the ground
        r = self.cls(c)
        self.assertLessEqual(r["diagnostics"]["zi_over_L"], 1)
        self.assertGreater(r["diagnostics"]["rib_low"], 0.25)
        self.assertEqual(r["class"], "C6")

    def test_c7_c8_transitions(self):
        day = column("2023-05-11T00:00Z", shf=40.0, lapse=0.003)
        night = column("2023-05-11T01:00Z", shf=-20.0, lapse=0.003)
        late = column("2023-05-11T02:00Z", shf=-25.0, lapse=0.003)
        self.assertEqual(self.cls(night, before=day, after=late)["class"], "C7")
        self.assertEqual(self.cls(day, after=night)["class"], "C7")
        self.assertNotIn(self.cls(late, before=night)["class"], ("C7", "C8"))
        dawn = column("2023-05-11T12:00Z", shf=-10.0, lapse=0.003)
        morning = column("2023-05-11T13:00Z", shf=60.0, lapse=0.003)
        self.assertEqual(self.cls(morning, before=dawn)["class"], "C8")
        self.assertEqual(self.cls(dawn, after=morning)["class"], "C8")

    def test_transition_needs_flux_threshold(self):
        # +3 to -3 W/m2 stays inside flux_eps (about 0.0025 K m/s): no transition
        a = column("2023-05-11T00:00Z", shf=3.0)
        b = column("2023-05-11T01:00Z", shf=-3.0)
        self.assertNotIn(self.cls(b, before=a)["class"], ("C7", "C8"))

    def test_transition_needs_sign_change(self):
        # 6.4 to 5.9 W/m2 crosses flux_eps but stays positive: no transition either way
        a = column("2023-05-11T00:00Z", shf=6.4, lapse=0.003)
        b = column("2023-05-11T01:00Z", shf=5.9, lapse=0.003)
        q = [hr.virtual_flux(x["surface"]) for x in (a, b)]
        self.assertTrue(q[0] > self.T["flux_eps"] > q[1] > 0, q)
        self.assertNotIn(self.cls(b, before=a)["class"], ("C7", "C8"))
        self.assertNotIn(self.cls(a, after=b)["class"], ("C7", "C8"))

    def test_c9_cloud(self):
        cl = np.zeros(len(Z))
        cl[9] = 1e-4                               # at 1400 m
        self.assertEqual(self.cls(column(shf=250.0, cloud=list(cl)))["class"], "C9")
        cl = np.zeros(len(Z))
        cl[-1] = 1e-4                              # at 5500 m, above the 3 km domain top
        self.assertNotEqual(self.cls(column(shf=250.0, hpbl=1500.0, cloud=list(cl)))["class"], "C9")

    def test_modifiers(self):
        c = column(shf=-15.0, ustar=0.35, hpbl=200.0, lapse=0.01, wind=[3.0] * len(Z), relief=300.0,
                   grad=(3e-5, 0.0))
        self.assertEqual(self.cls(c)["modifiers"], ["Fr", "B"])
        self.assertEqual(self.cls(column(shf=-15.0, ustar=0.35, hpbl=200.0, lapse=0.003))["modifiers"], [])

    def test_tendency(self):
        a = column("2023-05-11T00:00Z", wind=[5.0] * len(Z))
        b = column("2023-05-11T01:00Z", wind=[6.5] * len(Z))
        c = column("2023-05-11T02:00Z", wind=[8.0] * len(Z))
        r = self.cls(b, before=a, after=c)
        # 0-1000 m mean of a wind 0 at the ground, linear to the first level (10 m), uniform above: 0.995 of the
        # uniform value; 3 m/s over 2 h
        self.assertAlmostEqual(r["diagnostics"]["tendency"], 3.0 * 0.995 / 2, places=9)
        self.assertIn("T", r["modifiers"])
        self.assertNotIn("T", self.cls(b, before=column("2023-05-11T00:00Z", wind=[6.0] * len(Z)))["modifiers"])

    def test_neighbours_must_be_one_hour(self):
        with self.assertRaisesRegex(hr.RegimeError, "not one hour before"):
            self.cls(column("2023-05-11T01:00Z"), before=column("2023-05-10T23:00Z"))

    def test_notes(self):
        self.assertIn("no neighbouring hours: transitions and T not tested", self.cls(column())["notes"])


class Thresholds(unittest.TestCase):
    def test_errors_name_the_key(self):
        for t, msg in (({"neutral_lenght": 1.0}, "unknown key 'neutral_lenght'"),
                       ({"jet_drop_frac": 1.5}, "'jet_drop_frac' = 1.5 is outside"),
                       ({"rib_crit": 0.0}, "'rib_crit' = 0.0 is outside"),
                       ({"flux_eps": -1.0}, "'flux_eps' = -1.0 is outside"),
                       ({"froude": float("nan")}, "'froude' must be a finite number"),
                       ({"tendency": True}, "'tendency' must be a finite number"),
                       ({"geo_shear": "2"}, "'geo_shear' must be a finite number")):
            with self.assertRaisesRegex(hr.RegimeError, msg):
                hr.check_thresholds(t)
        self.assertEqual(hr.check_thresholds({"flux_eps": 0.0})["flux_eps"], 0.0)

    def test_file(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.json")
            with open(p, "w") as f:
                json.dump({"neutral_length": 300}, f)
            t = hr.read_thresholds(p)
            self.assertEqual(t["neutral_length"], 300.0)
            self.assertEqual(t["rib_crit"], 0.25)

    def test_threshold_changes_class(self):
        c = column(shf=-10.0, ustar=0.6, lapse=0.003)
        L = hr.classify(c)["diagnostics"]["L"]
        self.assertEqual(hr.classify(c)["class"], "C3")
        self.assertNotEqual(hr.classify(c, hr.check_thresholds({"neutral_length": 2 * L}))["class"], "C3")

    def test_column_checks(self):
        c = column()
        c["z"][3] = c["z"][2]
        with self.assertRaisesRegex(hr.RegimeError, "'z' must be positive and increase"):
            hr.check_column(c, "col")
        c = column()
        for us in (0.0, 1e-120):                 # 1e-120 would underflow u*^3 to L = -0
            c["surface"]["ustar"] = us
            with self.assertRaisesRegex(hr.RegimeError, "surface 'ustar' must be >= 1e-3"):
                hr.check_column(c, "col")
        for edit, msg in ((lambda c: c["surface"].update(shf="10"), "surface 'shf' missing or not a number"),
                          (lambda c: c["surface"].update(hpbl=None), "surface 'hpbl'"),
                          (lambda c: c.update(lat="36"), "'lat' must be a number"),
                          (lambda c: c["z"].__setitem__(0, "a"), "'z' must be a list of numbers"),
                          (lambda c: c["surface"].update(q2m=-1 / 0.61), "'q2m' must be in"),
                          (lambda c: c["theta_v"].__setitem__(2, 0.0), "'theta_v' must be > 0")):
            c = column()
            edit(c)
            with self.assertRaisesRegex(hr.RegimeError, msg):
                hr.check_column(c, "col")


class Catalogue(unittest.TestCase):
    def test_jobs_and_neighbours(self):
        hours = [("2023-05-10T23:00Z", 40.0), ("2023-05-11T00:00Z", 30.0), ("2023-05-11T01:00Z", -20.0),
                 ("2023-05-11T02:00Z", -25.0), ("2023-06-01T12:00Z", 200.0)]
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for v, h in hours:
                p = os.path.join(d, v[:13].replace("-", "").replace("T", "") + ".json")
                hr.write_json(column(v, shf=h, lapse=0.003, hpbl=800.0), p)
                paths.append(p)
            one = hr.catalogue(paths[::-1], hr.default_thresholds(), jobs=1)
            four = hr.catalogue(paths, hr.default_thresholds(), jobs=4)
            self.assertEqual(json.dumps(one, sort_keys=True), json.dumps(four, sort_keys=True))
            self.assertEqual([r["class"] for r in one], ["C2", "C7", "C7", "C6", "C1"])
            self.assertEqual(one[4]["notes"], ["no neighbouring hours: transitions and T not tested"])
            table = hr.count_table(one)
            self.assertEqual(table.splitlines()[0].split(), ["class", "2023-05", "2023-06", "total"])
            self.assertEqual(table.splitlines()[7].split(), ["C7", "2", "0", "2"])
            self.assertEqual(table.splitlines()[1].split(), ["C1", "0", "1", "1"])
            hr.write_json(column("2023-05-11T00:00Z"), os.path.join(d, "dup.json"))
            with self.assertRaisesRegex(hr.RegimeError, "also in"):
                hr.catalogue(paths + [os.path.join(d, "dup.json")], hr.default_thresholds())

    def test_one_bad_hour_is_refused_alone(self):
        # an arithmetic error in one hour becomes that hour's refusal, not the end of the catalogue
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "a.json")
            hr.write_json(column(), p)
            real = hr.classify
            try:
                hr.classify = lambda *a, **k: 1 / 0
                r = hr._classify_job((p, None, None, hr.default_thresholds()))
            finally:
                hr.classify = real
            self.assertEqual(r, {"file": p, "error": "division by zero"})

    def test_command_line(self):
        def run(*argv):
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = hr.main(list(argv))
            return code, out.getvalue() + err.getvalue()
        with tempfile.TemporaryDirectory() as d:
            good = os.path.join(d, "a.json")
            hr.write_json(column("2023-05-10T05:00Z", shf=-20.0, lapse=0.003), good)
            short = os.path.join(d, "b.json")            # the jet search passes the top of a 1 km column
            c = column("2023-05-10T06:00Z", shf=-20.0, lapse=0.003)
            for k in ("z", "u", "v", "theta_v", "cloud_max", "dthdx", "dthdy"):
                c[k] = c[k][:9]
            hr.write_json(c, short)
            bad = os.path.join(d, "t.json")
            with open(bad, "w") as f:
                json.dump({"neutral_lenght": 1}, f)
            code, msg = run("catalogue", good, short)
            self.assertEqual(code, 1)
            self.assertIn("b.json: refused:", msg)
            self.assertIn("refused: 1", msg)
            self.assertEqual(run("catalogue", good)[0], 0)
            self.assertEqual(run("catalogue", good, "--jobs", "0"), (2, "hrrr_regime.py: --jobs must be >= 1\n"))
            code, msg = run("catalogue", good, "--thresholds", bad)
            self.assertEqual(code, 2)
            self.assertIn("unknown key 'neutral_lenght'", msg)
            code, msg = run("classify", os.path.join(d, "missing.json"))
            self.assertEqual(code, 2)
            self.assertIn("missing.json", msg)
            code, msg = run("fetch", "2023-05-10T05:00Z", "--lat", "36", "--lon", "-97", "--grad-half-km", "nan",
                            "--save", d, "--out-dir", os.path.join(d, "never"))
            self.assertEqual(code, 2)
            self.assertIn("--grad-half-km", msg)
            self.assertFalse(os.path.exists(os.path.join(d, "never")))     # refused before any folder or download


if __name__ == "__main__":
    unittest.main()
