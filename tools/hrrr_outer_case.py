"""The ``case`` command of hrrr_outer.py: the three decks of the HRRR chain from one gust_case YAML file.

The YAML file is a gust_case case description of the nest (domain, terrain, mesh, surface, ...) with one more
section, ``hrrr``, read here and not by make_gust_case.py:

    hrrr:
      file: awaken_2023051019_f01.npz   # the crop (hrrr_outer.py fetch); it carries the cycle and forecast hour
      parent_ratio: 4                   # bc.nest_ratio of the nest: 2, 4, 8 or 16
      raw_margin: 2000.                 # the raw HRRR domain reaches this far beyond the parent, m
      theta: solved                     # solved (gravity on) or frozen (simplec.theta = frozen)
      surface: neutral                  # neutral; heat_flux with HRRR's box-mean kinematic flux; or
                                        # surface_temperature with HRRR's skin theta_v as a map
      geostrophic: false                # true: forcing.geostrophic with HRRR's G (mean over 200-1000 m)
      cloud_mixing_ratio: 1.0e-5        # cloud water above this in the raw domain refuses the hour, kg/kg

The meshes. The nest has n = 2 half_width / cell_size cells a side, a multiple of 8. The parent's cells are
parent_ratio times larger, its half-width the smallest multiple of 4 parent cells that is at least twice the nest's;
the raw domain's cells are twice the parent's, its half-width the smallest multiple of 4 raw cells at least the
parent's plus raw_margin. Every face of a finer mesh then lies on a face of the coarser one. The blend of the
terrain into the parent's is two of the coarser cells (terrain.nest_blend_width). All three share the vertical mesh
and the terrain file; the domain top must be at least 3000 m.

The scheme keys (written to the decks; a gust_inputs key of the case that sets one of them is refused): four nest
faces with bc.nest_parent, init.first_guess = nest, Coriolis on at the centre's latitude, bc.theta_profile from
HRRR's box-mean theta_v over the nest, forcing.theta0 its ground value, bc.zhi_theta_grad its gradient over the top
100 m, the Rayleigh layer on w and theta, and turb.max_geom_lscale = 100 m. The raw and parent runs only write their
first guess (no steady stage).
"""
import copy
import math
import os
import shlex
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hrrr_outer as ho  # noqa: E402
from gust_case.builder import build  # noqa: E402
from gust_case.errors import CaseError  # noqa: E402
from gust_case.schema import load_yaml, validate  # noqa: E402

HRRR_KEYS = {
    "file": None,
    "parent_ratio": 4,
    "raw_margin": 2000.0,
    "theta": "solved",
    "surface": "neutral",
    "geostrophic": False,
    "cloud_mixing_ratio": ho.CLOUD_MIXING_RATIO,
}
LEVELS = ("hrrr_raw", "parent", "nest")
SURFACE_THETA_FILE = "surface_theta.xyz"     # hrrr.surface = surface_temperature: x y theta_v,s, nest folder
SCHEME_KEYS = ("bc.xlo", "bc.xhi", "bc.ylo", "bc.yhi", "bc.nest_parent", "bc.nest_ratio", "init.first_guess",
               "terrain.nest_blend_width", "forcing.coriolis", "forcing.latitude", "forcing.geostrophic",
               "forcing.ug", "forcing.vg", "forcing.gravity", "forcing.theta0", "bc.theta_profile_z",
               "bc.theta_profile", "bc.zhi_theta_grad", "bc.rayleigh_fields", "abl.theta", "abl.drive",
               "simplec.theta", "bc.surface_stability", "bc.surface_heat_flux", "bc.surface_temperature",
               "bc.surface_temperature_file", "steady.stage",
               "turb.max_geom_lscale", "io.write_checkpoint")


def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def read_hrrr_section(raw, base_dir):
    """The hrrr section with defaults, every key range-checked (refused by name, never clamped)."""
    sec = raw.pop("hrrr", None)
    if not isinstance(sec, dict):
        raise ho.HrrrOuterError("the case file needs an hrrr section (hrrr.file at least)")
    bad = sorted(set(sec) - set(HRRR_KEYS))
    if bad:
        raise ho.HrrrOuterError("unknown key%s %s; known: %s" % ("s" if len(bad) > 1 else "",
                                ", ".join("hrrr." + b for b in bad), ", ".join(HRRR_KEYS)))
    h = dict(HRRR_KEYS)
    h.update(sec)
    if not isinstance(h["file"], str) or not h["file"]:
        raise ho.HrrrOuterError("hrrr.file: the crop (.npz) is required")
    h["file"] = os.path.normpath(os.path.join(base_dir, h["file"]))
    if not os.path.isfile(h["file"]):
        raise ho.HrrrOuterError(f"hrrr.file: {h['file']}: no such file")
    for k in ("raw_margin", "cloud_mixing_ratio"):
        if isinstance(h[k], str):
            raise ho.HrrrOuterError(f"hrrr.{k} ({h[k]!r}) must be a number (YAML reads 1e-5 as text: write "
                                    "1.0e-5)")
    if type(h["parent_ratio"]) is not int or h["parent_ratio"] not in (2, 4, 8, 16):
        raise ho.HrrrOuterError(f"hrrr.parent_ratio ({h['parent_ratio']}) must be 2, 4, 8 or 16")
    if not (_num(h["raw_margin"]) and h["raw_margin"] > 0):
        raise ho.HrrrOuterError(f"hrrr.raw_margin ({h['raw_margin']}) must be > 0 m")
    if h["theta"] not in ("solved", "frozen"):
        raise ho.HrrrOuterError(f"hrrr.theta ({h['theta']}) must be solved or frozen")
    if h["surface"] not in ("neutral", "heat_flux", "surface_temperature"):
        raise ho.HrrrOuterError(f"hrrr.surface ({h['surface']}) must be neutral, heat_flux or surface_temperature")
    if not isinstance(h["geostrophic"], bool):
        raise ho.HrrrOuterError(f"hrrr.geostrophic ({h['geostrophic']}) must be true or false")
    if not (_num(h["cloud_mixing_ratio"]) and h["cloud_mixing_ratio"] > 0):
        raise ho.HrrrOuterError(f"hrrr.cloud_mixing_ratio ({h['cloud_mixing_ratio']}) must be > 0 kg/kg")
    return h


def check_case(raw):
    """Refuse the case keys the chain sets itself, naming them."""
    problems = []
    for sec in ("wind", "atmosphere", "nests"):
        if sec in raw:
            problems.append(f"{sec}: set by hrrr_outer.py from HRRR (remove the section)")
    d = raw.get("domain") or {}
    for k in ("west", "east", "south", "north", "center"):
        if k in d:
            problems.append(f"domain.{k}: the chain needs a square domain about domain.center_lat_lon "
                            "(domain.half_width, one number)")
    if "center_lat_lon" not in d:
        problems.append("domain.center_lat_lon is required (HRRR is placed by latitude and longitude)")
    if d.get("origin", "center") != "center":
        problems.append("domain.origin must be center")
    hw = d.get("half_width")
    if isinstance(hw, list):
        if len(hw) != 1:
            problems.append("domain.half_width: one number (a square domain)")
    top = d.get("top")
    if _num(top) and top < ho.MIN_TOP:
        problems.append(f"domain.top ({top}) must be at least {ho.MIN_TOP:.0f} m: HRRR's convergence is balanced "
                        "near 3 km, and a lower slip-wall top forces the adjustment to move every face velocity")
    if _num(top) and top > ho.HTOP:
        problems.append(f"domain.top ({top}) must be at most {ho.HTOP:.0f} m (the top of the HRRR height table)")
    if "utm_zone" in d and "center_lat_lon" in d:
        z0 = ho.geo.utm_zone(float(d["center_lat_lon"][1]))
        if d["utm_zone"] != z0:
            problems.append(f"domain.utm_zone ({d['utm_zone']}): HRRR is placed in the zone of the centre ({z0})")
    if (raw.get("output") or {}).get("prefix"):
        problems.append("output.prefix: run.sh scores the nest's pltNNNNN; leave the prefix empty")
    m = raw.get("mesh") or {}
    if m.get("hierarchy_levels", 1) != 1:
        problems.append("mesh.hierarchy_levels must be 1: nest faces refuse grid sequencing")
    if (raw.get("solver") or {}).get("grid_sequencing"):
        problems.append("solver.grid_sequencing must be false: nest faces refuse it")
    cap = (raw.get("model") or {}).get("length_scale_cap", "auto")
    if cap not in ("auto", 100, 100.0):
        problems.append(f"model.length_scale_cap ({cap}): the chain fixes 100 m")
    for k in sorted(raw.get("gust_inputs") or {}):
        if k in SCHEME_KEYS:
            problems.append(f"gust_inputs.{k}: set by the chain")
    if problems:
        raise ho.HrrrOuterError("; ".join(problems))


def plan_meshes(half, cell, ratio, margin):
    """(half-width, cell) of the nest, parent and raw domain; refuses a nest that does not tile."""
    n = 2.0 * half / cell
    if abs(n - round(n)) > 1e-9 * n or round(n) % 8:
        raise ho.HrrrOuterError(f"2 domain.half_width / mesh.cell_size ({n:g}) must be a whole multiple of 8")
    m = round(n) // 8                      # half = 4 cell m
    if (4 * m) % ratio:
        raise ho.HrrrOuterError(f"hrrr.parent_ratio ({ratio}): the nest's faces (half-width {half:g} m = "
                                f"{4 * m} cells) must lie on parent faces; 4 x {m} must be a multiple of {ratio}")
    cp = ratio * cell
    hp = 4.0 * cp * math.ceil(2.0 * half / (4.0 * cp) - 1e-9)
    cr = 2.0 * cp
    hr = 4.0 * cr * math.ceil((hp + margin) / (4.0 * cr) - 1e-9)
    return {"nest": (half, cell), "parent": (hp, cp), "hrrr_raw": (hr, cr)}


STAGE_NAMESPACES = ("forcing.", "convection.", "simplec.", "turb.", "steady.")
STAGE_KEYS = ("io.turb_plot_file", "abl.drive", "abl.theta")


def strip_stage_keys(path):
    """Remove from a first-guess deck (no steady stage) the keys GUST refuses without one (Source/App.cpp and
    Source/AblSetup.cpp: forcing., convection., simplec., turb., steady.*, io.turb_plot_file, abl.drive, abl.theta).
    Returns the removed keys."""
    keep, removed = [], []
    for line in open(path):
        key = line.split("=", 1)[0].strip() if "=" in line.split("#", 1)[0] else ""
        if key and (key.startswith(STAGE_NAMESPACES) or key in STAGE_KEYS):
            removed.append(key)
        else:
            keep.append(line)
    with open(path, "w") as f:
        f.write("".join(keep))
    return removed


def _level_case(raw, name, out_dir, half, cell, gust_inputs, wind):
    r = copy.deepcopy(raw)
    r.setdefault("case", {})
    r["case"]["name"] = name
    r["case"]["output_dir"] = out_dir
    r["domain"]["half_width"] = float(half)
    r["mesh"]["cell_size"] = float(cell)
    r["mesh"]["hierarchy_levels"] = 1
    r.setdefault("solver", {})["grid_sequencing"] = False
    r.setdefault("output", {})["checkpoint"] = True
    r["wind"] = wind
    gi = dict(r.get("gust_inputs") or {}) if name == "nest" else {}
    gi.update(gust_inputs)
    r["gust_inputs"] = gi
    return r


def build_chain(yaml_path, output_dir=None, overwrite=False, dry_run=False):
    """Write <output_dir>/<case.name>/{hrrr_raw, parent, nest}/ and run.sh; returns a dict of what was set."""
    base_dir = os.path.dirname(os.path.abspath(yaml_path))
    try:
        raw = load_yaml(yaml_path)
    except CaseError as err:
        raise ho.HrrrOuterError("; ".join(err.problems))
    h = read_hrrr_section(raw, base_dir)
    check_case(raw)
    lat0, lon0 = (float(v) for v in raw["domain"]["center_lat_lon"])
    field = ho.HrrrField(h["file"], lat0, lon0)
    try:
        case, _ = validate(copy.deepcopy(dict(raw, wind={"speed": 1.0, "height": 80.0, "direction": 0.0})))
    except CaseError as err:
        raise ho.HrrrOuterError("; ".join(err.problems))
    half = case["domain"]["half_width"]
    half = float(half[0] if isinstance(half, (list, tuple)) else half)
    top = float(case["domain"]["top"])
    if case["mesh"]["cell_size"] is None:
        raise ho.HrrrOuterError("mesh.cell_size is required")
    meshes = plan_meshes(half, float(case["mesh"]["cell_size"]), h["parent_ratio"], h["raw_margin"])
    field.check_cloud(meshes["hrrr_raw"][0], top, h["cloud_mixing_ratio"])
    # every mesh must lie inside the crop: locate() refuses a point outside it
    hr = meshes["hrrr_raw"][0]
    field.locate(np.array([-hr, hr, -hr, hr]), np.array([-hr, -hr, hr, hr]))

    # wind of the raw and parent decks' (unused) inflow profile: HRRR's box mean at 80 m over the nest
    m = field.box(max(half, 3000.0))
    k80 = int(np.searchsorted(field.hg, 80.0))
    um, vm = field.tab["u"][m, k80].mean(), field.tab["v"][m, k80].mean()
    wind = {"speed": float(max(np.hypot(um, vm), 0.5)), "height": 80.0,
            "direction": float(np.degrees(np.arctan2(-um, -vm)) % 360.0)}
    field.check_top(top)
    tz, th, grad = ho.theta_profile(field, half, top)
    case_name = (raw.get("case") or {}).get("name", "gust_case")
    out_root = output_dir or (raw.get("case") or {}).get("output_dir", ".")
    out_root = out_root if os.path.isabs(out_root) else os.path.join(base_dir if output_dir is None else os.getcwd(),
                                                                     out_root)
    chain_dir = os.path.normpath(os.path.join(out_root, case_name))
    if os.path.exists(chain_dir) and not overwrite and not dry_run:
        raise ho.HrrrOuterError(f"{chain_dir} exists; use --overwrite")

    first_only = {"steady.stage": None}
    faces = {f"bc.{f}": "nest" for f in ("xlo", "xhi", "ylo", "yhi")}
    nest_gi = dict(faces)
    nest_gi.update({
        "bc.nest_parent": "../parent/chk00000", "bc.nest_ratio": h["parent_ratio"], "init.first_guess": "nest",
        "terrain.nest_blend_width": 2.0 * meshes["parent"][1], "forcing.coriolis": 1, "forcing.latitude": lat0,
        "forcing.gravity": None, "abl.theta": None, "abl.drive": None,
        "forcing.theta0": round(th[0], 3), "bc.theta_profile_z": " ".join(f"{z:g}" for z in tz),
        "bc.theta_profile": " ".join(f"{t:.3f}" for t in th), "bc.zhi_theta_grad": round(grad, 6),
        "bc.rayleigh_fields": "w theta", "turb.max_geom_lscale": 100.0, "io.write_checkpoint": 1,
    })
    info = {"hour": field.label, "meshes": meshes, "wind": wind, "theta_z": tz, "theta": th, "top_gradient": grad}
    if h["theta"] == "frozen":
        nest_gi["simplec.theta"] = "frozen"
    if h["geostrophic"]:
        ug, vg = field.geostrophic(half, lat0)
        nest_gi.update({"forcing.geostrophic": 1, "forcing.ug": round(ug, 3), "forcing.vg": round(vg, 3)})
        info["geostrophic"] = (ug, vg)
    if h["surface"] == "heat_flux":
        q = field.surface_flux(half)
        nest_gi.update({"bc.surface_stability": "heat_flux", "bc.surface_heat_flux": round(q, 5)})
        info["heat_flux"] = q
    if h["surface"] == "surface_temperature":
        field.surface_theta()      # refuses a crop without the skin temperature before anything is written
        nest_gi.update({"bc.surface_stability": "surface_temperature",
                        "bc.surface_temperature_file": SURFACE_THETA_FILE})
    parent_gi = dict(faces)
    parent_gi.update(first_only)
    parent_gi.update({"bc.nest_parent": "../hrrr_raw/chk00000", "init.first_guess": "nest",
                      "terrain.nest_blend_width": 2.0 * meshes["hrrr_raw"][1]})
    level_gi = {"hrrr_raw": dict(first_only), "parent": parent_gi, "nest": nest_gi}

    for name in LEVELS:
        hw, cell = meshes[name]
        r = _level_case(raw, name, chain_dir, hw, cell, level_gi[name], wind)
        try:
            c, given = validate(r)
            build(c, given, base_dir, dry_run=dry_run, overwrite=overwrite, check=False)
        except CaseError as err:
            raise ho.HrrrOuterError(f"{name}: " + "; ".join(err.problems))
        if not dry_run:      # the level's own run.sh would run it alone; the chain's run.sh runs the three
            os.remove(os.path.join(chain_dir, name, "run.sh"))
            if name != "nest":
                strip_stage_keys(os.path.join(chain_dir, name, "inputs"))
    if h["surface"] == "surface_temperature" and not dry_run:
        # HRRR's columns (3 km apart) over the nest and two of them beyond, so IDW has neighbours at the faces
        info["surface_temperature"] = field.write_surface_theta(os.path.join(chain_dir, "nest", SURFACE_THETA_FILE),
                                                                half + 6000.0)
    if not dry_run:
        _write_run_script(chain_dir, h["file"], lat0, lon0, case, h["cloud_mixing_ratio"])
    info["case_dir"] = chain_dir
    return info


def _write_run_script(chain_dir, crop, lat0, lon0, case, cloud):
    run = case["run"]
    mpi = run["mpi_launcher"].replace("{ranks}", str(run["ranks"]))
    q = shlex.quote
    tool = os.path.join(HERE, "hrrr_outer.py")
    site = f"--lat {lat0!r} --lon {lon0!r}"
    lines = [
        "#!/bin/bash",
        "# The HRRR chain, written by tools/hrrr_outer.py case: raw HRRR -> adjusted parent -> nest, then the score.",
        "#   GUST=/path/to/gust MPI=\"srun -n 16\" PYTHON=python3 ./run.sh",
        "set -eo pipefail",
        'cd "$(dirname "$0")"',
        f"G=${{GUST:-{q(run['executable'])}}}; MPI=${{MPI-{q(mpi)}}}; PY=${{PYTHON:-python3}}",
        f"TOOL={q(tool)}; CROP={q(crop)}",
        'step() { (cd "$1" && rm -rf chk[0-9]* plt[0-9]* turb[0-9]* && $MPI "$G" inputs > log.txt 2>&1) || '
        '{ echo "run.sh: $1 failed (see $1/log.txt)" >&2; exit 1; }; }',
        "rm -rf score.txt score.json parent/chk00000",
        "step hrrr_raw",
        f'$PY "$TOOL" to-chk "$CROP" hrrr_raw/chk00000 {site} --cloud-mixing-ratio {cloud!r}',
        '$PY "$TOOL" balance hrrr_raw/chk00000',
        "step parent",
        "step nest",
        'plt=$(ls -d nest/plt[0-9]* | grep -E "/plt[0-9]+$" | sort | tail -1)',
        f'$PY "$TOOL" score "$CROP" "$plt" {site} --log nest/log.txt --json score.json | tee score.txt',
    ]
    path = os.path.join(chain_dir, "run.sh")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(path, 0o755)


def main_case(yaml_path, output_dir=None, overwrite=False):
    info = build_chain(yaml_path, output_dir, overwrite)
    print(f"{info['hour']}: wrote {info['case_dir']}/{{hrrr_raw, parent, nest}} and run.sh")
    for name in LEVELS:
        hw, cell = info["meshes"][name]
        print(f"  {name:8s} half-width {hw:8.0f} m, cells {cell:6.0f} m, {round(2 * hw / cell)} a side")
    print(f"  wind at 80 m (raw and parent inflow profile only) {info['wind']['speed']:.2f} m/s from "
          f"{info['wind']['direction']:.0f} degrees; theta_v {info['theta'][0]:.2f} K at the ground, "
          f"{info['top_gradient']:.4f} K/m at the top")
    if info.get("surface_temperature"):
        n, lo, hi, mean = info["surface_temperature"]
        print(f"  skin theta_v ({SURFACE_THETA_FILE}, {n} HRRR columns): {lo:.2f} to {hi:.2f} K, mean {mean:.2f} K")
