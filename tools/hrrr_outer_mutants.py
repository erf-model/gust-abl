#!/usr/bin/env python3
"""
The one-line mutants of tools/hrrr_outer.py and tools/hrrr_outer_case.py, each run against their self-test
(tools/test_hrrr_outer.py) in a scratch copy of the mutated file. Every mutant must make the self-test fail. The
table is in docs/testing.rst. Two lines are guards whose removal changes nothing today, so they have no mutant:
w = 0 on the ground faces (the table already gives 0 at height 0) and turb.max_geom_lscale = 100 in the nest's deck
(gust_case's model.length_scale_cap = auto writes 100 m too, and the case refuses any other value).

    python3 tools/hrrr_outer_mutants.py [O1 O4 ...]
"""
import os
import shutil
import subprocess
import sys
import tempfile

O, C = "hrrr_outer.py", "hrrr_outer_case.py"
M = [
 ("O1", "an f00 crop accepted", O, "if int(fxx) != fxx or fxx < 1:", "if int(fxx) != fxx or fxx < 0:"),
 ("O2", "the crop's winds not rotated", O,
  '        u, v = grid_to_utm(f("nat_u"), f("nat_v"), lat.ravel()[None], lon.ravel()[None], self.zone)',
  '        u, v = f("nat_u"), f("nat_v")'),
 ("O3", "HRRR's grid rotation left out", O, "    u_e, v_n = hrrr_regime.rotate_to_earth(ug, vg, lon)",
  "    u_e, v_n = ug, vg"),
 ("O4", "the UTM convergence left out", O,
  "    c = np.arctan(np.tan(np.radians(np.asarray(lon) - lon_cm)) * np.sin(np.radians(lat)))",
  "    c = 0.0 * np.asarray(lon)"),
 ("O5", "the UTM convergence turned the wrong way", O,
  "return u_e * np.cos(c) - v_n * np.sin(c), u_e * np.sin(c) + v_n * np.cos(c)",
  "return u_e * np.cos(c) + v_n * np.sin(c), -u_e * np.sin(c) + v_n * np.cos(c)"),
 ("O6", "a cloudy hour accepted", O, "        if c.max() > limit:", "        if False:"),
 ("O7", "to-chk without the cloud check", O, "    field.check_cloud(half, hi[2], cloud_limit)\n", ""),
 ("O8", "cloud above the top counted", O, "low = self.h_levels[:, m] <= top", "low = self.h_levels[:, m] <= 1e9"),
 ("O9", "theta, not theta_v", O, "thv = T * (1e5 / p) ** KAPPA_P * (1.0 + 0.61 * q)", "thv = T * (1e5 / p) ** KAPPA_P"),
 ("O10", "omega taken as w (the sign)", O, 'w = -f("nat_w") * RD * tv / (p * G)', 'w = f("nat_w") * RD * tv / (p * G)'),
 ("O11", "linear, not the log law, below the first level", O,
  "a3[0, c] * np.log(zz / self.z0[c]) / math.log(h[0, c] / self.z0[c])", "a3[0, c] * self.hg[below] / h[0, c]"),
 ("O13", "heights above sea level, not above the ground", O, "    agl = zn - zn[0:1]", "    agl = zn"),
 ("O14", "u at the heights of the v faces", O, '"u": (xn, yc, _avg(_avg(agl, 0), 1))',
  '"u": (xn, yc, np.concatenate([_avg(_avg(agl, 0), 2)[:, :, :1], _avg(_avg(agl, 0), 2)], axis=2)[:, :-1])'),
 ("O15", "the score's ground from one node", O,
  "ground = 0.25 * (zn[0, :-1, :-1] + zn[0, 1:, :-1] + zn[0, :-1, 1:] + zn[0, 1:, 1:])", "ground = zn[0, :-1, :-1]"),
 ("O16", "\"not steady after\" read as steady", O, 'out["steady"] = (v[-1] == "steady") if v else None',
  'out["steady"] = ("steady after" in t) if v else None'),
 ("O17", "the iterations not read", O, "it = re.findall(r\"stage account: (\\d+) iterations\", t)", "it = []"),
 ("O18", "cloud water (CLMR) not downloaded", O, "|TKE|CLMR|CLWMR):", "|TKE|CLWMR):"),
 ("O19", "the latent part of the virtual heat flux left out", O,
  "(1.0 + 0.61 * q2) + 0.61 * theta * f(\"sfc_slhtf\") / (rho * LV))", "(1.0 + 0.61 * q2) + 0.0 * theta)"),
 ("O20", "a geostrophic fit on extrapolated pressures accepted", O, "        if g0 + z_lo < lowest:",
  "        if False:"),
 ("O21", "a top above HRRR's highest level accepted", O, "        if top > self.h_top:", "        if False:"),
 ("O22", "a top above the height table accepted", O, "0 < top <= HTOP):", "0 < top):"),
 ("O23", "the cloud check without the crop's reach", O,
  "        self.locate(np.array([-r, r, -r, r]), np.array([-r, -r, r, r]))\n", ""),
 ("O24", "w_max over the whole nest", O, '"w_max": float(np.abs(w).max())', '"w_max": float(np.abs(c["w"]).max())'),
 ("O25", "a refused command exits with 0", O, "        return 1\n", "        return 0\n"),
 ("C1", "a top below 3000 m accepted", C, "    if _num(top) and top < ho.MIN_TOP:", "    if False:"),
 ("C2", "a hierarchy of levels accepted", C, '    if m.get("hierarchy_levels", 1) != 1:', "    if False:"),
 ("C3", "a scheme key in gust_inputs accepted", C, "        if k in SCHEME_KEYS:", "        if False:"),
 ("C4", "the wind section accepted", C, '    for sec in ("wind", "atmosphere", "nests"):',
  '    for sec in ("atmosphere", "nests"):'),
 ("C5", "a nest that does not tile accepted", C, "if abs(n - round(n)) > 1e-9 * n or round(n) % 8:", "if False:"),
 ("C6", "nest faces off the parent's accepted", C, "    if (4 * m) % ratio:", "    if False:"),
 ("C7", "the parent only as wide as the nest", C, "hp = 4.0 * cp * math.ceil(2.0 * half / (4.0 * cp) - 1e-9)",
  "hp = 4.0 * cp * math.ceil(1.0 * half / (4.0 * cp) - 1e-9)"),
 ("C8", "any parent_ratio accepted", C, 'h["parent_ratio"] not in (2, 4, 8, 16):', "False:"),
 ("C13", "a float parent_ratio accepted", C, 'if type(h["parent_ratio"]) is not int or', "if"),
 ("C14", "a cloud_mixing_ratio given as text accepted", C, "        if isinstance(h[k], str):", "        if False:"),
 ("C15", "a top above the height table accepted", C, "    if _num(top) and top > ho.HTOP:", "    if False:"),
 ("C16", "a utm_zone other than the centre's accepted", C, '        if d["utm_zone"] != z0:', "        if False:"),
 ("C17", "output.prefix accepted", C, '    if (raw.get("output") or {}).get("prefix"):', "    if False:"),
 ("C9", "an unknown hrrr key accepted", C, "    if bad:", "    if False:"),
 ("C10", "the steady.* keys left in the first-guess decks", C,
  'STAGE_NAMESPACES = ("forcing.", "convection.", "simplec.", "turb.", "steady.")',
  'STAGE_NAMESPACES = ("forcing.", "convection.", "simplec.", "turb.")'),
 ("C11", "the parent's decks not stripped", C, '            if name != "nest":', '            if name == "hrrr_raw":'),
 ("C12", "hrrr.theta = frozen not written", C, '        nest_gi["simplec.theta"] = "frozen"', "        pass"),
]


def run(selected):
    here = os.path.dirname(os.path.abspath(__file__))
    survived = []
    for tag, what, f, old, new in M:
        if selected and tag not in selected:
            continue
        src = open(os.path.join(here, f)).read()
        if src.count(old) != 1:
            raise SystemExit(f"{tag}: the line to mutate is not found once in {f}: {old!r}")
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, f), "w") as fh:
                fh.write(src.replace(old, new))
            for g in (O, C, "test_hrrr_outer.py"):
                if g != f:
                    shutil.copy(os.path.join(here, g), d)
            env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": here}   # the rest of tools/ after d
            r = subprocess.run([sys.executable, "-m", "unittest", "test_hrrr_outer"], cwd=d, capture_output=True,
                               text=True, env=env)
        out = r.stdout + r.stderr
        if r.returncode != 0 and "FAILED" not in out:
            raise SystemExit(f"{tag}: the self-test did not run:\n{out[-2000:]}")
        killed = r.returncode != 0
        print(f"{tag:4s} {'killed' if killed else 'SURVIVED':8s} {what}")
        if not killed:
            survived.append(tag)
    if survived:
        raise SystemExit(f"{len(survived)} mutant(s) survived: {' '.join(survived)}")


if __name__ == "__main__":
    run(set(sys.argv[1:]))
