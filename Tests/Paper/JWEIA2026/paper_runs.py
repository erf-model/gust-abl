#!/usr/bin/env python3
"""Every run of the GUST paper (realizability sensor, JWEIA 2026), one after another, with the case decks of
Tests/Real_Terrain and their prepare / compare / sample scripts.

    python3 paper_runs.py prepare --out RUNS --wntools ~/windninja-tools --iz-data IZDIR --rushil-data RUSHILDIR
    python3 paper_runs.py run --gust build-rel/gust --out RUNS [--np 8] [GROUPS]
    python3 paper_runs.py stats --out RUNS

Every run gets its own deck copy in RUNS/runs/<tag>/inputs (keys replaced or added; None removes a key), so no key
is given twice. A run with a `done` file is not run again. Groups (default MAHSCE):
  M  main comparison: every case without and with the sensor (the three-curve figures, Tables 4-6)
  A  separation of the limit and the blending: limit only; first-order upwind without and with the limit (Table 7)
  H  empty-domain test (Section 5.1)
  S  sensitivity (Appendix B): sigma_max, Tucker wall distance, consistent inflow, C_r sweep, H_b, freeze threshold,
     the smoothed Butte terrain
  C  Durbin's C_r = 0.245 (Appendix A)
  E  grid sensitivity: coarse and fine grids, with the sensor and with the limit only (Table 8)
The Kato-Launder runs of Section 5.7 used the exploratory option turb.production of commit dfb0443, which was
removed again (cf867cf); they are not repeated here.
"""
import argparse
import glob
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RT = os.path.join(REPO, "Tests", "Real_Terrain")
IZ_PY = os.path.join(RT, "IshiharaZhou", "ishihara_zhou.py")
RU_PY = os.path.join(RT, "RUSHIL", "rushil.py")
WN_PY = os.path.join(RT, "WindNinja2019", "windninja2019.py")
Z0_PY = os.path.join(RT, "WindNinja2019", "bolund_z0_map.py")
WN = os.path.join(RT, "WindNinja2019")
ASK_TERRAIN = os.path.join(REPO, "extern_erf", "askervein.txt")


def cases(out):
    """name: (deck, keys added to the deck, H_b = 2.5 times the terrain height (m), postprocess kind)."""
    w, iz, ru = f"{out}/wn", f"{out}/iz", f"{out}/rushil"
    c = {
        "ridge_smooth": (f"{RT}/IshiharaZhou/inputs_iz_smooth", {"terrain.file": f"{iz}/ridge.xyz",
                         "bc.inflow_file": f"{iz}/inflow_smooth.txt", "turb.max_geom_lscale": "100"}, 100, ("iz", "smooth")),
        "ridge_rough": (f"{RT}/IshiharaZhou/inputs_iz_rough", {"terrain.file": f"{iz}/ridge.xyz",
                        "bc.inflow_file": f"{iz}/inflow_rough.txt", "turb.max_geom_lscale": "100"}, 100, ("iz", "rough")),
        "askervein": (f"{WN}/inputs_askervein", {"terrain.file": ASK_TERRAIN, "turb.max_geom_lscale": "100"}, 290,
                      ("wn", "askervein")),
        "bsb": (f"{WN}/inputs_bsb_srtm", {"terrain.file": f"{w}/bsb.xyz"}, 1965, ("wn", "bsb")),
        "bsb_smooth": (f"{WN}/inputs_bsb_wn_consistent", {"terrain.file": f"{w}/bsb_wn.xyz"}, 1965, ("wn", "bsb_smooth")),
    }
    for h in (3, 5, 8):
        c[f"hill_{h}"] = (f"{RT}/RUSHIL/inputs_rushil", {"terrain.file": f"{ru}/hill{h}.xyz",
                          "bc.inflow_file": f"{ru}/inflow.txt"}, 292, ("rushil", h))
    for k in (270, 239, 90):  # the water/land roughness map in place of the decks' single z0
        c[f"bolund_{k}"] = (f"{WN}/inputs_bolund_{k}", {"terrain.file": f"{w}/bolund.xyz", "turb.max_geom_lscale": "100",
                            "bc.surface_z0": None, "bc.surface_z0_file": f"{w}/bolund_z0.xyz"}, 30, ("wn", f"bolund_{k}"))
    return c


OFF = {"turb.realizability": "0", "turb.nut_relax": "1", "convection.hybrid_sigma_max": "0", "convection.hybrid_height": None}
CAPONLY = {"turb.realizability": "0.309", "turb.nut_relax": "0.5", "convection.hybrid_sigma_max": "0",
           "convection.hybrid_height": None}
UP = {"convection.limiter": "upwind"}


def sens(Hb, sig="0.1"):
    return {"turb.realizability": "0.309", "turb.nut_relax": "0.5", "convection.hybrid_sigma_max": sig,
            "convection.hybrid_height": f"{Hb:g}"}


def deck_get(text, key):
    m = re.search(r"^" + re.escape(key) + r"\s*=\s*(.*)$", text, re.M)
    return m.group(1).split("#")[0].strip() if m else None


def write_deck(src, keys, dst):
    text = open(src).read()
    for k, v in keys.items():
        pat = re.compile(r"^" + re.escape(k) + r"\s*=.*$", re.M)
        if v is None:
            text = pat.sub("", text)
        elif pat.search(text):
            text = pat.sub(f"{k} = {v}", text)
        else:
            text += f"\n{k} = {v}\n"
    open(dst, "w").write(text)
    return text


def top_height(text):
    n = int(deck_get(text, "mesh.n_cell").split()[2])
    dz0, r = float(deck_get(text, "mesh.dz0")), float(deck_get(text, "mesh.stretch_ratio"))
    return dz0 * (r ** n - 1) / (r - 1)


def consistent(text):
    """Uncapped length and a no-slip lid moving at the log law's top speed, in place of the slip top with damping."""
    us, z0 = float(deck_get(text, "bc.inflow_ustar")), float(deck_get(text, "bc.inflow_z0"))
    d = float(deck_get(text, "bc.inflow_direction") or 0)
    U = us / 0.41 * math.log((top_height(text) + z0) / z0)
    return {"turb.max_geom_lscale": "1.e6", "bc.zhi": "no_slip_wall",
            "bc.zhi_wall_velocity": f"{U * math.cos(math.radians(d)):.4f} {U * math.sin(math.radians(d)):.4f} 0.",
            "bc.rayleigh_zdamp": None, "bc.rayleigh_dampcoef": None, "bc.rayleigh_fields": None}


def refine(text, f):
    """f = 2 (fine) or 0.5 (coarse) in x, y (not the 8 spanwise cells of the 2-D cases) and z."""
    nx, ny, nz = map(int, deck_get(text, "mesh.n_cell").split())
    dz0, r = float(deck_get(text, "mesh.dz0")), float(deck_get(text, "mesh.stretch_ratio"))
    ny2 = ny if ny == 8 else int(ny * f)
    return {"mesh.n_cell": f"{int(nx * f)} {ny2} {int(nz * f)}", "mesh.dz0": f"{dz0 / f:.6g}",
            "mesh.stretch_ratio": f"{r ** (1 / f):.8f}"}


def refine_h(text, f):
    nx, ny, nz = map(int, deck_get(text, "mesh.n_cell").split())
    return {"mesh.n_cell": f"{int(nx * f)} {int(ny * f)} {nz}"}


def lastplt(d):
    p = sorted(x for x in glob.glob(f"{d}/plt[0-9]*") if "level" not in x)
    return p[-1] if p else None


class Runner:
    def __init__(self, a):
        self.a, self.out = a, os.path.abspath(a.out)
        self.cases = cases(self.out)

    def post(self, c, d, tag):
        kind, p = self.cases[c][3], lastplt(d)
        if not p:
            return "no plotfile"
        o = self.out
        if kind[0] == "iz":
            cmd = [IZ_PY, "compare", p, "--surface", kind[1], "--outdir", f"{o}/iz", "--label", tag, "--out", f"{d}/compare.csv"]
        elif kind[0] == "rushil":
            cmd = [RU_PY, "compare", p, "--hill", str(kind[1]), "--outdir", f"{o}/rushil", "--label", tag,
                   "--out", f"{d}/compare.csv"]
        else:
            site, case = {"askervein": ("askervein", "210"), "bsb": ("bsb", "222"),
                          "bsb_smooth": ("bsb", "222")}.get(kind[1]) or ("bolund", kind[1].split("_")[1])
            ter = (["--erf_terrain", ASK_TERRAIN] if site == "askervein" else
                   ["--terrain", f"{o}/wn/bsb_wn.xyz" if kind[1] == "bsb_smooth" else f"{o}/wn/{site}.xyz"])
            z0 = "0.03" if site == "askervein" else "0.01"
            cmd = [WN_PY, "sample", p, "--site", site, "--case", case, "--masts", f"{o}/wn/masts.csv", *ter,
                   "--z0", z0, "--out", f"{o}/samples/{tag}.csv"]
        r = subprocess.run([sys.executable, *cmd], capture_output=True, text=True)
        return " ".join(l for l in r.stdout.splitlines() if "reversed" in l)[:250] or r.stderr.strip()[-200:]

    def run(self, tag, c, keys, extra=None, deck=None):
        d = f"{self.out}/runs/{tag}"
        os.makedirs(d, exist_ok=True)
        os.makedirs(f"{self.out}/samples", exist_ok=True)
        src, base, _, _ = self.cases[c]
        k = dict(base)
        k.update(keys)
        k.setdefault("steady.max_iterations", str(self.a.max_iterations))
        text = write_deck(deck or src, k, f"{d}/inputs")
        if extra:
            write_deck(f"{d}/inputs", extra(text), f"{d}/inputs")
        if not os.path.exists(f"{d}/done"):
            with open(f"{d}/log.txt", "w") as lg:
                r = subprocess.run(["mpirun", "-np", str(self.a.np), os.path.abspath(self.a.gust), "inputs"], cwd=d,
                                   stdin=subprocess.DEVNULL, stdout=lg, stderr=subprocess.STDOUT)
            if r.returncode == 0:
                open(f"{d}/done", "w").close()
        st = [l for l in open(f"{d}/log.txt") if re.search(r"steady after|not steady|diverged|!!!", l)]
        print(time.strftime("%H:%M"), tag, "|", st[0].strip() if st else "no verdict", "|", self.post(c, d, tag), flush=True)

    def hhabl(self):
        for tag, deck, keys in (("H_default", "inputs_hhabl_base", {}),
                                ("H_sensor", "inputs_hhabl_base", {**sens(290), "turb.max_geom_lscale": "100"}),
                                ("H_lscale100", "inputs_hhabl_base", {"turb.max_geom_lscale": "100"}),
                                ("H_consistent", "inputs_hhabl_consistent", {"turb.max_geom_lscale": "1.e6"}),
                                ("H_consistent_sensor", "inputs_hhabl_consistent", {**sens(290), "turb.max_geom_lscale": "1.e6"})):
            d = f"{self.out}/runs/{tag}"
            os.makedirs(d, exist_ok=True)
            write_deck(f"{HERE}/{deck}", keys, f"{d}/inputs")
            if not os.path.exists(f"{d}/done"):
                with open(f"{d}/log.txt", "w") as lg:
                    r = subprocess.run(["mpirun", "-np", str(self.a.np), os.path.abspath(self.a.gust), "inputs"], cwd=d,
                                       stdin=subprocess.DEVNULL, stdout=lg, stderr=subprocess.STDOUT)
                if r.returncode == 0:
                    open(f"{d}/done", "w").close()
            print(time.strftime("%H:%M"), tag, flush=True)


def do_run(a):
    R = Runner(a)
    C, g = R.cases, a.groups
    H = {c: v[2] for c, v in C.items()}
    two_d = ["ridge_smooth", "ridge_rough", "hill_3", "hill_5", "hill_8"]
    field = ["askervein", "bolund_270", "bolund_239", "bolund_90"]
    if "M" in g:
        for c in two_d + field + ["bsb"]:
            R.run(f"M_nosensor_{c}", c, OFF)
            R.run(f"M_sensor_{c}", c, sens(H[c]))
    if "A" in g:
        for c in two_d + field + ["bsb"]:
            R.run(f"A_caponly_{c}", c, CAPONLY)
            R.run(f"A_upwind_nocap_{c}", c, {**OFF, **UP})
            R.run(f"A_upwind_cap_{c}", c, {**CAPONLY, **UP})
    if "H" in g:
        R.hhabl()
    if "S" in g:
        for c in ["ridge_smooth", "hill_3", "bolund_270"]:
            for s in ("0.05", "0.2"):
                R.run(f"S_sigma{s}_{c}", c, sens(H[c], s))
        for c in ["ridge_smooth", "ridge_rough", "hill_3", "bolund_270", "bolund_90"]:
            R.run(f"S_tucker_{c}", c, {**sens(H[c]), "turb.wall_distance": "tucker"})
        for c in field:
            R.run(f"S_consistent_nosensor_{c}", c, OFF, consistent)
            R.run(f"S_consistent_sensor_{c}", c, sens(H[c]), consistent)
        for cr in ("0.33", "0.35", "0.40"):
            for c in ["ridge_smooth", "hill_3", "hill_5", "bolund_270"]:
                R.run(f"S_cr{cr}_{c}", c, {**sens(H[c]), "turb.realizability": cr})
        for c in ["bolund_270", "bsb"]:
            for f, nm in ((1.5, "1.5x"), (4, "4x")):  # 1.5 and 4 times the terrain height instead of 2.5
                R.run(f"S_Hb{nm}_{c}", c, {**sens(H[c]), "convection.hybrid_height": f"{f * H[c] / 2.5:g}"})
        for c in ["ridge_smooth", "bsb"]:
            R.run(f"S_freeze1e-3_{c}", c, {**sens(H[c]), "convection.hybrid_freeze": "1.e-3"})
        R.run("S_sensor_bsb_smooth", "bsb_smooth", sens(H["bsb_smooth"]))
    if "C" in g:
        for c in two_d + field:
            R.run(f"C_cr245_{c}", c, {**sens(H[c]), "turb.realizability": "0.245"})
    if "E" in g:
        for c in ["ridge_smooth", "hill_3"]:
            for f, nm in ((0.5, "coarse"), (2, "fine")):
                R.run(f"E_sensor_{nm}_{c}", c, sens(H[c]), lambda text, f=f: refine(text, f))
                R.run(f"E_caponly_{nm}_{c}", c, CAPONLY, lambda text, f=f: refine(text, f))
        for f, nm in ((0.5, "coarse_h"), (2, "fine_h")):
            R.run(f"E_sensor_{nm}_bolund_270", "bolund_270", sens(H["bolund_270"]), lambda text, f=f: refine_h(text, f))
    print("PAPER RUNS DONE", flush=True)


def do_prepare(a):
    o = os.path.abspath(a.out)
    for sub in ("wn", "iz", "rushil"):
        os.makedirs(f"{o}/{sub}", exist_ok=True)
    py = [sys.executable]
    subprocess.run([*py, WN_PY, "prepare", a.wntools, f"{o}/wn"], check=True)
    subprocess.run([*py, WN_PY, "smooth", f"{a.wntools}/test_cases/big_butte.asc", f"{o}/wn/bsb_wn.xyz",
                    "--sigma", "120", "--cap", "2050", "--soft", "150"], check=True)
    subprocess.run([*py, "-I", Z0_PY, f"{o}/wn/bolund.xyz", f"{o}/wn/bolund_z0.xyz"], check=True)
    subprocess.run([*py, IZ_PY, "prepare", a.iz_data, f"{o}/iz"], check=True)
    subprocess.run([*py, RU_PY, "prepare", a.rushil_data, f"{o}/rushil"], check=True)


def do_stats(a):
    o = os.path.abspath(a.out)
    for pattern in ("M_*", "A_*", "S_*", "C_*", "E_*"):
        csvs = sorted(glob.glob(f"{o}/samples/{pattern}.csv"))
        if csvs:
            subprocess.run([sys.executable, WN_PY, "stats", *csvs, "--masts", f"{o}/wn/masts.csv"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--out", required=True)
    p.add_argument("--wntools", required=True, help="checkout of github.com/nwagenbrenner/windninja-tools")
    p.add_argument("--iz-data", required=True, help="the supplementary files of Ishihara and Zhou (2025)")
    p.add_argument("--rushil-data", required=True, help="the ERCOFTAC AC5-05 (RUSHIL) files")
    r = sub.add_parser("run")
    r.add_argument("groups", nargs="?", default="MAHSCE")
    r.add_argument("--gust", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--np", type=int, default=8)
    r.add_argument("--max-iterations", type=int, default=4000,
                   help="stalled runs stop here (the recommended setting converges in fewer than 900)")
    s = sub.add_parser("stats")
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    {"prepare": do_prepare, "run": do_run, "stats": do_stats}[a.cmd](a)


if __name__ == "__main__":
    main()
