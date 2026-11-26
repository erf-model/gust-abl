#!/usr/bin/env python3
"""The C1 and C2 options: for each convective or weakly unstable hour, the chain of case.yaml with that hour's crop
(tools/hrrr_outer.py case) and the options beside it. Every option changes only the nest, so all share the chain's
parent (a nest_<option> folder beside the chain's nest):

    scheme        the scheme as hrrr_outer.py case writes it: theta solved, neutral ground (the chain's own nest)
    solved_flux   theta solved, bc.surface_stability = heat_flux with HRRR's box-mean virtual heat flux
    frozen        simplec.theta = frozen (theta held at the parent's), neutral ground
    frozen_flux   frozen, with HRRR's heat flux at the ground
    frozen_buoyant  simplec.theta = frozen_buoyant: theta held at the parent's, its buoyancy kept in w
    frozen_buoyant_flux  frozen_buoyant with HRRR's heat flux at the ground
    frozen_tsurf, frozen_buoyant_tsurf, solved_tsurf
                  frozen, frozen_buoyant and solved theta with bc.surface_stability = surface_temperature: HRRR's skin
                  theta_v as a map (surface_theta.xyz in the option's folder), the flux from Monin-Obukhov
                  similarity cell by cell; solved_tsurf stops at SOLVED_FLUX_ITERATIONS
    <tsurf option>_kb2  the three _tsurf options with a thermal roughness bc.surface_kb_inverse = 2
    <tsurf option>_all  the three _tsurf options with kB^-1 = 2, bc.surface_stable_functions = jimenez and
                  bc.surface_wstar = beljaars with HRRR's boundary-layer height as a map (pblh.xyz in the folder)
    solved_flux_l200  solved_flux with the length-scale cap at 200 m (turb.max_geom_lscale; the scheme's is 100 m)
    frozen_buoyant_l200  frozen_buoyant with the length-scale cap at 200 m (the recommended option for C1, README)
    relaxed_flux_T  simplec.theta = relaxed_buoyant with simplec.theta_relax_time = T s (T in RELAX_TIMES): theta
                  solved and relaxed towards the parent's, its buoyancy kept in w, HRRR's heat flux at the ground

The script writes <out>/<hour>/ and <out>/run_options.sh, which runs each chain and then the nest options of its
hour, one run at a time, and scores each nest. An hour the tool refuses is reported and left out.

    python3 options.py [--data DIR] [--out DIR (options)] [--hours TAG ...] [--options NAME ...] [--overwrite]
"""
import argparse
import copy
import os
import shlex
import shutil
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools"))
sys.path.insert(0, TOOLS)
sys.path.insert(0, HERE)
import hrrr_outer as ho  # noqa: E402
import hrrr_outer_case as hc  # noqa: E402
from prepare import HOURS, LAT, LON, crop_name  # noqa: E402

# the theta-solved heated surface had no steady state in 4000 iterations (Nov 2026); 2000 are enough to show it
SOLVED_FLUX_ITERATIONS = 2000
# the relaxation times of simplec.theta = relaxed_buoyant, s (z_i / w* is about 500-600 s in these hours)
RELAX_TIMES = (300, 1000, 3000)

# bc.surface_stability = surface_temperature: the map of HRRR's skin theta_v written into the option's folder
TSURF = {"bc.surface_stability": "surface_temperature", "bc.surface_temperature_file": hc.SURFACE_THETA_FILE}
# the thermal roughness kB^-1 = ln(z0 / z0t) of the _kb2 and _all options
KB = {"bc.surface_kb_inverse": 2}
# the map of HRRR's boundary-layer height written into an _all option's folder (bc.surface_pblh_file)
PBLH_FILE = "pblh.xyz"
# the _all options: every surface feature (thermal roughness, Jimenez stable functions, Beljaars w*)
ALL = dict(KB, **{"bc.surface_stable_functions": "jimenez", "bc.surface_wstar": "beljaars",
                  "bc.surface_pblh_file": PBLH_FILE})


def nest_options(q):
    """The nest options for an hour with HRRR's box-mean kinematic virtual heat flux q (K m/s):
    name -> (keys set, keys dropped)."""
    flux = {"bc.surface_stability": "heat_flux", "bc.surface_heat_flux": round(q, 5)}
    opts = {"frozen": ({"simplec.theta": "frozen"}, ()),
            "frozen_flux": (dict(flux, **{"simplec.theta": "frozen"}), ()),
            "frozen_buoyant": ({"simplec.theta": "frozen_buoyant"}, ()),
            "frozen_buoyant_flux": (dict(flux, **{"simplec.theta": "frozen_buoyant"}), ()),
            "solved_flux": (dict(flux, **{"steady.max_iterations": SOLVED_FLUX_ITERATIONS}), ()),
            "frozen_tsurf": (dict(TSURF, **{"simplec.theta": "frozen"}), ()),
            "frozen_buoyant_tsurf": (dict(TSURF, **{"simplec.theta": "frozen_buoyant"}), ()),
            "solved_tsurf": (dict(TSURF, **{"steady.max_iterations": SOLVED_FLUX_ITERATIONS}), ()),
            "frozen_buoyant_l200": ({"simplec.theta": "frozen_buoyant", "turb.max_geom_lscale": 200.0}, ()),
            "solved_flux_l200": (dict(flux, **{"steady.max_iterations": SOLVED_FLUX_ITERATIONS,
                                                "turb.max_geom_lscale": 200.0}), ())}
    for v in ("frozen_tsurf", "frozen_buoyant_tsurf", "solved_tsurf"):
        opts[f"{v}_kb2"] = (dict(opts[v][0], **KB), ())
        opts[f"{v}_all"] = (dict(opts[v][0], **ALL), ())
    for t in RELAX_TIMES:
        relax = {"simplec.theta": "relaxed_buoyant", "simplec.theta_relax_time": t,
                 "steady.max_iterations": SOLVED_FLUX_ITERATIONS}
        opts[f"relaxed_flux_{t}"] = (dict(flux, **relax), ())
    return opts


NEST_OPTIONS = tuple(nest_options(0.1))
OPTIONS = ("scheme",) + NEST_OPTIONS


def edit_deck(path, set_keys, drop_keys=()):
    """Rewrite a deck with set_keys replaced or appended and drop_keys removed."""
    lines, seen = [], set()
    for line in open(path):
        key = line.split("=", 1)[0].strip() if "=" in line.split("#", 1)[0] else ""
        if key in drop_keys:
            continue
        if key in set_keys:
            line = f"{key} = {set_keys[key]}\n"
            seen.add(key)
        lines.append(line)
    lines += [f"{k} = {v}\n" for k, v in set_keys.items() if k not in seen]
    with open(path, "w") as fh:
        fh.write("".join(lines))


def score_cmd(crop, log="log.txt"):
    q = shlex.quote
    return (f"plt=$(ls -d plt[0-9]* | grep -E '^plt[0-9]+$' | sort | tail -1) && "
            f"$PY \"$TOOL\" score {q(crop)} \"$plt\" --lat {LAT} --lon {LON} --log {log} --json score.json "
            f"> score.txt")


def write_case(base, tag, data, out):
    raw = copy.deepcopy(base)
    name = tag
    raw["case"]["name"] = name
    raw["case"]["output_dir"] = out
    raw["terrain"]["file"] = os.path.join(data, "awaken.asc")
    raw["hrrr"]["file"] = os.path.join(data, crop_name(tag))
    yml = os.path.join(out, f"{name}.yaml")
    with open(yml, "w") as fh:
        yaml.safe_dump(raw, fh, sort_keys=False)
    return yml, raw


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=HERE, help="the folder prepare.py wrote (awaken.asc and the crops)")
    ap.add_argument("--out", default=os.path.join(HERE, "options"), help="default: options/ beside this script")
    ap.add_argument("--hours", nargs="+", default=list(HOURS), choices=list(HOURS))
    ap.add_argument("--options", nargs="+", default=list(OPTIONS), choices=list(OPTIONS))
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    data = os.path.abspath(a.data)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    base = yaml.safe_load(open(os.path.join(HERE, "case.yaml")))
    q = shlex.quote
    run = ["#!/bin/bash", "# The C1 and C2 options, written by options.py: each chain, then the nest options of its hour.",
           "#   GUST=/path/to/gust MPI=\"mpiexec -n 4\" PYTHON=python3 ./run_options.sh", "set -u -o pipefail",
           'cd "$(dirname "$0")"', ': "${GUST:?set GUST to the gust executable}"', 'MPI=${MPI-"mpiexec -n 4"}',
           "PY=${PYTHON:-python3}", f'TOOL={q(os.path.join(TOOLS, "hrrr_outer.py"))}']
    refused = []
    if not a.overwrite:
        names = list(a.hours)
        there = [n for n in names if os.path.exists(os.path.join(out, n))]
        if there:
            sys.exit(f"options.py: {', '.join(there)} already in {out}; pass --overwrite to write them again")
    for tag in a.hours:
        yml, raw = write_case(base, tag, data, out)
        try:
            info = hc.build_chain(yml, overwrite=a.overwrite)
        except ho.HrrrOuterError as err:
            print(f"{tag}: refused: {err}")
            refused.append(tag)
            continue
        chain = info["case_dir"]
        crop = raw["hrrr"]["file"]
        field = ho.HrrrField(crop, LAT, LON)
        half = raw["domain"]["half_width"]
        hf = field.surface_flux(half)
        zi = float(np.asarray(field.d["sfc_blh"], float).ravel()[field.box(half)].mean())
        print(f"{tag}: {info['hour']} -> {chain}; heat flux {hf:.4f} K m/s, boundary-layer height {zi:.0f} m")
        opts = nest_options(hf)
        run.append(f'if (cd {tag} && GUST="$GUST" MPI="$MPI" PYTHON="$PY" ./run.sh > chain.log 2>&1); then')
        for v in (o for o in a.options if o in opts):
            d = os.path.join(chain, f"nest_{v}")
            if os.path.exists(d):
                shutil.rmtree(d)
            os.makedirs(d)
            for f in ("inputs", "terrain.xyz"):
                shutil.copy(os.path.join(chain, "nest", f), d)
            edit_deck(os.path.join(d, "inputs"), *opts[v])
            if opts[v][0].get("bc.surface_stability") == "surface_temperature":
                # HRRR's columns over the nest and 6 km beyond (as hrrr_outer.py case writes it)
                n, lo, hi, mean = field.write_surface_theta(os.path.join(d, hc.SURFACE_THETA_FILE), half + 6000.0)
                print(f"  nest_{v}: skin theta_v of {n} HRRR columns, {lo:.2f} to {hi:.2f} K, mean {mean:.2f} K")
            if opts[v][0].get("bc.surface_pblh_file") == PBLH_FILE:
                n, lo, hi, mean = field.write_pblh(os.path.join(d, PBLH_FILE), half + 6000.0)
                print(f"  nest_{v}: boundary-layer height of {n} HRRR columns, {lo:.0f} to {hi:.0f} m, "
                      f"mean {mean:.0f} m")
            run.append(f"  (cd {tag}/nest_{v} && rm -rf plt[0-9]* turb[0-9]* chk[0-9]* score.* && "
                       f"$MPI \"$GUST\" inputs > log.txt 2>&1 && {score_cmd(crop)}) || echo \"{tag}/nest_{v} failed\"")
        run.append(f'else echo "{tag}: chain failed (see {tag}/chain.log); nest options not run"; fi')
    path = os.path.join(out, "run_options.sh")
    with open(path, "w") as fh:
        fh.write("\n".join(run) + "\n")
    os.chmod(path, 0o755)
    print(f"wrote {path}" + (f"; refused: {', '.join(refused)}" if refused else ""))


if __name__ == "__main__":
    try:
        main()
    except ho.HrrrOuterError as err:
        sys.exit(f"options.py: {err}")
