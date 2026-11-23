#!/usr/bin/env python3
"""The variant matrix: for each hour, the chain of case.yaml with that hour's crop (tools/hrrr_outer.py case), and
beside its nest six variants of the nest's deck that share the chain's parent:

    S  the scheme as written by hrrr_outer.py: theta solved, Coriolis on, neutral surface (the chain's own nest)
    A  theta frozen, Coriolis off
    B  theta frozen
    C  B + geostrophic forcing (HRRR's G: plane fit of p at constant height, 200-1000 m mean)
    D  geostrophic forcing, theta solved
    E  C + a relaxation zone of 8 cells inside the nest faces (bc.nest_relax_width)
    F  C + HRRR's surface heat flux (bc.surface_stability = heat_flux)

A-F are the variants of the November 2026 study, kept for comparison; B, C and F are also reachable from the hrrr
section (theta, geostrophic, surface). The script writes <out>/<hour>/ for each hour and <out>/run_matrix.sh, which
runs each chain (raw, parent, S and its score) and, when the chain succeeded, A-F with their scores, one run at a
time. An hour the tool refuses (a cloudy hour, for one) is reported and left out.

    python3 variants.py [--data DIR] [--out DIR] [--hours day20 ...] [--overwrite]
"""
import argparse
import copy
import os
import shlex
import shutil
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools"))
sys.path.insert(0, TOOLS)
sys.path.insert(0, HERE)
import hrrr_outer as ho  # noqa: E402
import hrrr_outer_case as hc  # noqa: E402
from prepare import HOURS, LAT, LON, crop_name  # noqa: E402

VARIANTS = "ABCDEF"


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


def variant_keys(v, field, half):
    frozen = {"simplec.theta": "frozen"}
    ug, vg = field.geostrophic(half, LAT)
    geo = {"forcing.geostrophic": 1, "forcing.ug": round(ug, 3), "forcing.vg": round(vg, 3)}
    q = field.surface_flux(half)
    return {"A": (frozen, ("forcing.coriolis", "forcing.latitude")),
            "B": (frozen, ()),
            "C": ({**frozen, **geo}, ()),
            "D": (geo, ()),
            "E": ({**frozen, **geo, "bc.nest_relax_width": 8}, ()),
            "F": ({**frozen, **geo, "bc.surface_stability": "heat_flux", "bc.surface_heat_flux": round(q, 5)}, ())}[v]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=HERE, help="the folder prepare.py wrote (awaken.asc and the crops)")
    ap.add_argument("--out", default=None, help="default: the data folder")
    ap.add_argument("--hours", nargs="+", default=list(HOURS), choices=list(HOURS))
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    data = os.path.abspath(a.data)
    out = os.path.abspath(a.out or data)
    base = yaml.safe_load(open(os.path.join(HERE, "case.yaml")))
    q = shlex.quote
    run = ["#!/bin/bash", "# The variant matrix, written by variants.py: each chain, then its variants A-F.",
           "#   GUST=/path/to/gust MPI=\"mpiexec -n 4\" PYTHON=python3 ./run_matrix.sh", "set -u -o pipefail",
           'cd "$(dirname "$0")"', ': "${GUST:?set GUST to the gust executable}"', 'MPI=${MPI:-"mpiexec -n 4"}',
           "PY=${PYTHON:-python3}", f'TOOL={q(os.path.join(TOOLS, "hrrr_outer.py"))}']
    refused = []
    for tag in a.hours:
        raw = copy.deepcopy(base)
        raw["case"]["name"] = tag
        raw["terrain"]["file"] = os.path.join(data, "awaken.asc")
        raw["hrrr"]["file"] = os.path.join(data, crop_name(tag))
        raw.setdefault("case", {})["output_dir"] = out
        yml = os.path.join(out, f"{tag}.yaml")
        os.makedirs(out, exist_ok=True)
        with open(yml, "w") as fh:
            yaml.safe_dump(raw, fh, sort_keys=False)
        try:
            info = hc.build_chain(yml, overwrite=a.overwrite)
        except ho.HrrrOuterError as err:
            print(f"{tag}: refused: {err}")
            refused.append(tag)
            continue
        chain = info["case_dir"]
        field = ho.HrrrField(raw["hrrr"]["file"], LAT, LON)
        half = raw["domain"]["half_width"]
        print(f"{tag}: {info['hour']} -> {chain}")
        run.append(f'if (cd {tag} && GUST="$GUST" MPI="$MPI" PYTHON="$PY" ./run.sh > chain.log 2>&1); then')
        for v in VARIANTS:
            d = os.path.join(chain, f"nest_{v}")
            if os.path.exists(d):
                shutil.rmtree(d)
            os.makedirs(d)
            for f in ("inputs", "terrain.xyz"):
                shutil.copy(os.path.join(chain, "nest", f), d)
            s, drop = variant_keys(v, field, half)
            edit_deck(os.path.join(d, "inputs"), s, drop)
            run.append(f"  (cd {tag}/nest_{v} && rm -rf plt[0-9]* turb[0-9]* chk[0-9]* score.* && "
                       f"$MPI \"$GUST\" inputs > log.txt 2>&1 && "
                       f"plt=$(ls -d plt[0-9]* | grep -E '^plt[0-9]+$' | sort | tail -1) && "
                       f"$PY \"$TOOL\" score {q(raw['hrrr']['file'])} \"$plt\" --lat {LAT} --lon {LON} "
                       f"--log log.txt --json score.json > score.txt) || echo \"{tag}/nest_{v} failed\"")
        run += [f'else echo "{tag}: chain failed (see {tag}/chain.log); variants not run"; fi']
    path = os.path.join(out, "run_matrix.sh")
    with open(path, "w") as fh:
        fh.write("\n".join(run) + "\n")
    os.chmod(path, 0o755)
    print(f"wrote {path}" + (f"; refused: {', '.join(refused)}" if refused else ""))


if __name__ == "__main__":
    try:
        main()
    except ho.HrrrOuterError as err:
        sys.exit(f"variants.py: {err}")
