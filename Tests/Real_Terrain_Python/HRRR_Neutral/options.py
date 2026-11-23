#!/usr/bin/env python3
"""The C3 options: for each neutral hour, the chain of case.yaml with that hour's crop (tools/hrrr_outer.py case)
and the options beside it.

Options that change only the nest share the chain's parent (a nest_<option> folder beside the chain's nest):

    scheme              the scheme as hrrr_outer.py case writes it: theta solved, Coriolis on, top 3 km, Rayleigh
                        layer on w and theta (the chain's own nest)
    frozen              simplec.theta = frozen (theta held at the parent's)
    no_coriolis         forcing.coriolis = 0
    frozen_no_coriolis  both
    no_rayleigh         no Rayleigh layer (bc.rayleigh_zdamp, bc.rayleigh_dampcoef, bc.rayleigh_fields removed, and
                        the inflow profile that serves only as the layer's reference)

Options that change the vertical mesh get a chain of their own (<hour>_<option>/), written and run through
wide_top.py, which widens the tool's top limits for them:

    top1500             top 1.5 km: 33 layers from 10 m
    top5000             top 5 km: 48 layers from 10 m

The script writes <out>/<hour>/, <out>/<hour>_1500/, <out>/<hour>_5000/ and <out>/run_options.sh, which runs
each chain and then the nest options of its hour, one run at a time, and scores each nest. An hour the tool refuses
is reported and left out.

    python3 options.py [--data DIR] [--out DIR (options)] [--hours TAG ...] [--options NAME ...] [--overwrite]
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
import wide_top  # noqa: E402
from prepare import HOURS, LAT, LON, crop_name  # noqa: E402

RAYLEIGH = ("bc.rayleigh_zdamp", "bc.rayleigh_dampcoef", "bc.rayleigh_fields")
# the inflow profile serves only the Rayleigh layer's reference here (four nest faces): without the layer GUST
# refuses these keys as having no effect
PROFILE = ("bc.inflow_profile", "bc.inflow_ustar", "bc.inflow_z0", "bc.inflow_direction")
NEST_OPTIONS = {
    "frozen": ({"simplec.theta": "frozen"}, ()),
    "no_coriolis": ({"forcing.coriolis": 0}, ("forcing.latitude",)),
    "frozen_no_coriolis": ({"simplec.theta": "frozen", "forcing.coriolis": 0}, ("forcing.latitude",)),
    "no_rayleigh": ({}, RAYLEIGH + PROFILE),
}
TOP_OPTIONS = {"top1500": (1500.0, 33), "top5000": (5000.0, 48)}       # top (m), layers
OPTIONS = ("scheme",) + tuple(NEST_OPTIONS) + tuple(TOP_OPTIONS)


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


def write_case(base, tag, data, out, top=None):
    raw = copy.deepcopy(base)
    name = tag if top is None else f"{tag}_{top[0]}"
    raw["case"]["name"] = name
    raw["case"]["output_dir"] = out
    raw["terrain"]["file"] = os.path.join(data, "awaken.asc")
    raw["hrrr"]["file"] = os.path.join(data, crop_name(tag))
    if top is not None:
        raw["domain"]["top"] = top[1]
        raw["mesh"]["layers"] = top[2]
        raw["mesh"].pop("stretch_ratio", None)       # found from the top, the layers and the first layer
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
    run = ["#!/bin/bash", "# The C3 options, written by options.py: each chain, then the nest options of its hour.",
           "#   GUST=/path/to/gust MPI=\"mpiexec -n 4\" PYTHON=python3 ./run_options.sh", "set -u -o pipefail",
           'cd "$(dirname "$0")"', ': "${GUST:?set GUST to the gust executable}"', 'MPI=${MPI-"mpiexec -n 4"}',
           "PY=${PYTHON:-python3}", f'TOOL={q(os.path.join(TOOLS, "hrrr_outer.py"))}']
    refused = []
    if not a.overwrite:
        names = [t for t in a.hours] + [f"{t}_{v[3:]}" for t in a.hours for v in a.options if v in TOP_OPTIONS]
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
        print(f"{tag}: {info['hour']} -> {chain}")
        run.append(f'if (cd {tag} && GUST="$GUST" MPI="$MPI" PYTHON="$PY" ./run.sh > chain.log 2>&1); then')
        for v in (o for o in a.options if o in NEST_OPTIONS):
            d = os.path.join(chain, f"nest_{v}")
            if os.path.exists(d):
                shutil.rmtree(d)
            os.makedirs(d)
            for f in ("inputs", "terrain.xyz"):
                shutil.copy(os.path.join(chain, "nest", f), d)
            edit_deck(os.path.join(d, "inputs"), *NEST_OPTIONS[v])
            run.append(f"  (cd {tag}/nest_{v} && rm -rf plt[0-9]* turb[0-9]* chk[0-9]* score.* && "
                       f"$MPI \"$GUST\" inputs > log.txt 2>&1 && {score_cmd(crop)}) || echo \"{tag}/nest_{v} failed\"")
        run.append(f'else echo "{tag}: chain failed (see {tag}/chain.log); nest options not run"; fi')
        for v in (o for o in a.options if o in TOP_OPTIONS):
            top, layers = TOP_OPTIONS[v]
            yml_t, _ = write_case(base, tag, data, out, (v[3:], top, layers))
            old = wide_top.widen()
            try:
                info_t = hc.build_chain(yml_t, overwrite=a.overwrite)
            except ho.HrrrOuterError as err:
                print(f"{tag}_{v}: refused: {err}")
                continue
            finally:
                wide_top.restore(old)
            # the chain's run.sh calls the tool for to-chk, balance and score: through wide_top.py instead
            rs = os.path.join(info_t["case_dir"], "run.sh")
            txt = open(rs).read().replace('$PY "$TOOL"', '$PY "$WIDE"')
            txt = txt.replace("TOOL=", f"WIDE={q(os.path.join(HERE, 'wide_top.py'))}; TOOL=", 1)
            with open(rs, "w") as fh:
                fh.write(txt)
            name = os.path.basename(info_t["case_dir"])
            print(f"{name}: top {top:.0f} m, {layers} layers -> {info_t['case_dir']}")
            run.append(f'(cd {name} && GUST="$GUST" MPI="$MPI" PYTHON="$PY" ./run.sh > chain.log 2>&1) || '
                       f'echo "{name} failed (see {name}/chain.log)"')
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
