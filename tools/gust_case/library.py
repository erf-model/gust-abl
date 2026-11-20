"""The geostrophic wind of every direction of a case, from its column runs.

    tools/make_gust_case.py --collect <case folder>

A geostrophic drive (``wind.drive = geostrophic_column`` or
``geostrophic_inflow``) runs a 1-D column per direction before the 3-D run
(builder.py). The column finds the geostrophic wind G that gives
``wind.speed`` at ``wind.height`` from that direction. This module reads what
each column run left and writes, in the case folder,

* ``library.json``: per direction G (components, speed and direction), the
  turning angle between the column's wind at ``wind.height`` and G, that wind,
  the column's u*, and, for a neutral column, the neutral geostrophic drag
  law's G and angle for that u* (for comparison);
* a section "Library" at the end of ``case_summary.txt`` with the same numbers
  as a table (written again on each collect, the earlier section replaced).

``run.sh`` calls it after the last run. It needs the case tool's own Python
packages (PyYAML and numpy), as make_gust_case.py does.

Where the numbers come from:

* G: with ``wind.speed`` at ``wind.height`` (the column adjusts the forcing to
  that wind), the column log's last lines ``forcing.ug = ...`` and
  ``forcing.vg = ...`` (9 digits; the values run.sh passes to the 3-D run),
  after the line "target wind reached"; a log without that line is refused.
  With ``atmosphere.geostrophic_wind`` there is no adjustment and G is the
  column deck's ``forcing.ug`` and ``forcing.vg``.
* z0 and the latitude: the column deck's ``bc.surface_z0`` and
  ``forcing.latitude``, the values the column ran with.
* the wind at ``wind.height``, u* and 1/L: the column table
  ``column/column_guess.dat`` (z u v theta k; u*, the surface heat flux and
  1/L in its header). The wind is GUST's own (Source/ColumnGuess.cpp,
  ``wind_at``), the one the column's target was met with: linear in
  ln(z + z0) between the two cell centres around the height; below the first
  centre the log law along the first centre's wind (neutral columns only).

Directions: meteorological, the direction the wind comes FROM in degrees
clockwise from north (wind.py). The turning angle is counter-clockwise
positive: G is turned clockwise from the wind near the ground in the northern
hemisphere, so it is negative there.
"""
import json
import math
import os
import re

import yaml

from . import wind
from .errors import CaseError

LIBRARY_FILE = "library.json"
SUMMARY_FILE = "case_summary.txt"
SUMMARY_MARK = "Library (tools/make_gust_case.py --collect)"
FORMAT_VERSION = 1
GEOSTROPHIC_DRIVES = ("geostrophic_column", "geostrophic_inflow")
# The wind at wind.height must match the wind asked for to these: GUST stops adjusting at |R - 1| < 1e-4, which
# bounds the speed's relative error by 1e-4 and the direction's by 1e-4 rad (0.0057 degrees).
SPEED_TOL = 1.0e-3
DIRECTION_TOL = 0.06

_TARGET = "column guess: target wind reached"
_FORCING = re.compile(r"^\s*forcing\.(ug|vg) = (\S+)\s*$")
_HEADER = {"ustar": re.compile(r"u\* (\S+) m/s"), "heat_flux": re.compile(r"surface heat flux (\S+) K m/s"),
           "inv_L": re.compile(r"1/L (\S+) 1/m")}


def _angle(a):
    """An angle in degrees brought to (-180, 180]."""
    a = math.fmod(a, 360.0)
    if a <= -180.0:
        a += 360.0
    elif a > 180.0:
        a -= 360.0
    return a


def _finite(path, what, *values):
    if not all(math.isfinite(v) for v in values):
        raise CaseError("%s: %s is not finite (%s)" % (path, what, ", ".join(repr(v) for v in values)))


def read_column_table(path):
    """(rows, header) of a column table: rows [z, u, v, theta, k] with z increasing, the first row the ground;
    header {"ustar", "heat_flux", "inv_L"} from the comment line (None where absent)."""
    rows, header = [], dict.fromkeys(_HEADER)
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#"):
                for key, rx in _HEADER.items():
                    m = rx.search(line)
                    if m:
                        try:
                            header[key] = float(m.group(1))
                        except ValueError:
                            raise CaseError("%s: header %s = %s is not a number" % (path, key, m.group(1)))
                        _finite(path, "the header's " + key, header[key])
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                row = [float(x) for x in parts]
            except ValueError:
                raise CaseError("%s, line %d: not numbers: %s" % (path, n, line.strip()))
            if len(row) != 5:
                raise CaseError("%s, line %d: %d columns, expected 5 (z u v theta k)" % (path, n, len(row)))
            _finite(path, "line %d" % n, *row)
            rows.append(row)
    if len(rows) < 3:
        raise CaseError("%s: a column table needs at least three levels, found %d" % (path, len(rows)))
    if any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
        raise CaseError("%s: the heights z do not increase" % path)
    return rows, header


def wind_at_height(rows, height, z0, inv_L=0.0):
    """(u, v) of a column table at ``height`` above the ground, as GUST's column computes it to meet its target
    (Source/ColumnGuess.cpp, wind_at). The first row is the ground (z = 0, no wind). Between cell centres: linear in
    ln(z + z0). Below the first centre: the log law along the first centre's wind,
    W1 ln((z + z0)/z0) / ln((z1 + z0)/z0); with 1/L not 0 (a stratified column) GUST's psi_m enters there, which this
    module does not evaluate, so that case is refused."""
    z = [r[0] for r in rows]
    if not 0.0 < height <= z[-1]:
        raise CaseError("wind.height %g m is outside the column table (0 to %g m)" % (height, z[-1]))
    if height < z[1]:
        if inv_L:
            raise CaseError("wind.height %g m is below the column's first cell centre (%g m) and the column is "
                            "stratified (1/L = %g 1/m): the profile there needs psi_m, which --collect does not "
                            "evaluate" % (height, z[1], inv_L))
        f = math.log((height + z0) / z0) / math.log((z[1] + z0) / z0)
        return f * rows[1][1], f * rows[1][2]
    k = max(j for j in range(1, len(z) - 1) if z[j] <= height)
    lo, hi = math.log(z[k] + z0), math.log(z[k + 1] + z0)
    t = (math.log(height + z0) - lo) / (hi - lo)
    return ((1.0 - t) * rows[k][1] + t * rows[k + 1][1], (1.0 - t) * rows[k][2] + t * rows[k + 1][2])


def _deck_value(path, key):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("#")[0].split()
            if len(parts) >= 3 and parts[0] == key and parts[1] == "=":
                try:
                    return float(parts[2])
                except ValueError:
                    raise CaseError("%s: %s = %s is not a number" % (path, key, parts[2]))
    return None


def column_deck_values(column_dir):
    """(z0, latitude) the column ran with: its deck's bc.surface_z0 and forcing.latitude."""
    deck = os.path.join(column_dir, "inputs")
    if not os.path.isfile(deck):
        raise CaseError("%s: no column deck" % deck)
    z0, lat = _deck_value(deck, "bc.surface_z0"), _deck_value(deck, "forcing.latitude")
    if z0 is None or not z0 > 0.0 or lat is None:
        raise CaseError("%s: bc.surface_z0 (> 0) and forcing.latitude are needed" % deck)
    _finite(deck, "bc.surface_z0, forcing.latitude", z0, lat)
    return z0, lat


def column_geostrophic_wind(column_dir, target):
    """(ug, vg) of a column run. ``target``: the case asks for a wind at a height, so the column adjusted G and its
    log must say "target wind reached" (the last forcing.ug, forcing.vg lines are read); otherwise G is the column
    deck's forcing.ug, forcing.vg."""
    log = os.path.join(column_dir, "log.txt")
    if target:
        if not os.path.isfile(log):
            raise CaseError("%s: no column log (run run.sh first)" % log)
        reached, last = False, {}
        with open(log, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                reached = reached or _TARGET in line
                m = _FORCING.match(line)
                if m:
                    try:
                        last[m.group(1)] = float(m.group(2))
                    except ValueError:
                        raise CaseError("%s: forcing.%s = %s is not a number" % (log, m.group(1), m.group(2)))
        if not reached:
            raise CaseError("%s: the column did not report \"target wind reached\" (did it finish? see the log)" % log)
        if "ug" not in last or "vg" not in last:
            raise CaseError("%s: the column reached its target wind but printed no forcing.ug, forcing.vg" % log)
        ug, vg = last["ug"], last["vg"]
        _finite(log, "forcing.ug, forcing.vg", ug, vg)
        return ug, vg
    deck = os.path.join(column_dir, "inputs")
    ug, vg = _deck_value(deck, "forcing.ug"), _deck_value(deck, "forcing.vg")
    if ug is None or vg is None:
        raise CaseError("%s: no forcing.ug, forcing.vg (atmosphere.geostrophic_wind)" % deck)
    _finite(deck, "forcing.ug, forcing.vg", ug, vg)
    return ug, vg


def is_geostrophic(case_dir):
    """True when the case folder's case_resolved.yaml has a geostrophic drive (a column run per direction)."""
    path = os.path.join(case_dir, "case_resolved.yaml")
    if not os.path.isfile(path):
        return False
    with open(path, encoding="utf-8") as fh:
        case = yaml.safe_load(fh)
    return isinstance(case, dict) and (case.get("wind") or {}).get("drive") in GEOSTROPHIC_DRIVES


def _settings(case_dir):
    path = os.path.join(case_dir, "case_resolved.yaml")
    if not os.path.isfile(path):
        raise CaseError("%s: not found; --collect needs a case folder written by make_gust_case.py" % path)
    try:
        with open(path, encoding="utf-8") as fh:
            case = yaml.safe_load(fh)
        w = case["wind"]
        name = case["case"]["name"]
    except (yaml.YAMLError, TypeError, KeyError) as err:
        raise CaseError("%s: not a case description written by make_gust_case.py (%s)" % (path, err))
    if w.get("drive") not in GEOSTROPHIC_DRIVES:
        raise CaseError("wind.drive = %s has no column run; --collect needs geostrophic_column or geostrophic_inflow"
                        % w.get("drive"))
    if w.get("directions") is not None:
        directions = wind.sweep_directions(w["directions"])
        folders = ["dir_" + wind.direction_label(d) for d in directions]
    else:
        directions = [float(w["direction"])]      # the case tool requires wind.direction without a sweep
        folders = ["."]
    target = w.get("height") is not None
    return {"name": name, "drive": w["drive"], "target": target,
            "speed": float(w["speed"]) if target else None, "height": float(w["height"]) if target else None,
            "directions": directions, "folders": folders}


def entry(folder, direction, speed, height):
    """The library entry of one direction's column run (folder: the run folder holding column/). ``speed`` and
    ``height``: the wind asked for (None when atmosphere.geostrophic_wind is given)."""
    column = os.path.join(folder, "column")
    table = os.path.join(column, "column_guess.dat")
    if not os.path.isfile(table):
        raise CaseError("%s: no column table (run run.sh first)" % table)
    rows, header = read_column_table(table)
    z0, latitude = column_deck_values(column)
    target = height is not None
    ug, vg = column_geostrophic_wind(column, target)
    g = math.hypot(ug, vg)
    towards_g = math.degrees(math.atan2(vg, ug))
    inv_L = header["inv_L"] or 0.0
    at_height, turning = None, None
    if target:
        u, v = wind_at_height(rows, height, z0, inv_L)
        towards_h = math.degrees(math.atan2(v, u))
        speed_h, from_h = math.hypot(u, v), wind.meteorological_from_towards(towards_h)
        if abs(speed_h / speed - 1.0) > SPEED_TOL or abs(_angle(from_h - direction)) > DIRECTION_TOL:
            raise CaseError("%s: the column's wind at %g m is %.6g m/s from %.4f degrees, not the %g m/s from %g "
                            "degrees asked for (a table from another run?)" % (table, height, speed_h, from_h, speed,
                                                                              direction))
        turning = _angle(towards_g - towards_h)
        at_height = {"height": height, "u": u, "v": v, "speed": speed_h, "direction": from_h}
    stress = _angle(towards_g - math.degrees(math.atan2(rows[1][2], rows[1][1])))
    neutral = not inv_L and not header["heat_flux"]
    out = {
        "direction": direction,
        "folder": None,
        "geostrophic_wind": {"ug": ug, "vg": vg, "speed": g,
                             "direction": wind.meteorological_from_towards(towards_g), "adjusted_by_column": target},
        "turning_deg": turning,
        "turning_from_stress_deg": stress,
        "wind_at_height": at_height,
        "wind_at_top": {"height": rows[-1][0], "u": rows[-1][1], "v": rows[-1][2]},
        "ustar": header["ustar"],
        "inv_L": header["inv_L"],
        "surface_heat_flux": header["heat_flux"],
        "roughness": z0,
        "latitude": latitude,
        "drag_law": None,
    }
    ustar = header["ustar"]
    if neutral and ustar is not None and ustar > 0.0 and abs(wind.coriolis_parameter(latitude)) > 0.0:
        alpha = wind.drag_law_turning(ustar, z0, latitude)
        out["drag_law"] = {"speed": wind.geostrophic_speed(ustar, z0, latitude),
                           "turning_deg": -math.copysign(alpha, latitude), "A": wind.DRAG_LAW_A, "B": wind.DRAG_LAW_B}
    return out


def collect(case_dir):
    """Read every direction's column run of a case folder; returns the library as a dict (not written)."""
    st = _settings(case_dir)
    entries = []
    for direction, folder in zip(st["directions"], st["folders"]):
        e = entry(os.path.normpath(os.path.join(case_dir, folder)), direction, st["speed"], st["height"])
        e["folder"] = folder
        entries.append(e)
    return {"format_version": FORMAT_VERSION, "case": st["name"], "drive": st["drive"],
            "reference_wind": {"speed": st["speed"], "height": st["height"]},
            "conventions": "direction: meteorological, degrees clockwise from north the wind comes from; turning_deg: "
                           "from the column's wind at the reference height to G, turning_from_stress_deg: from the "
                           "wind at the first cell centre (the surface stress) to G, both counter-clockwise positive",
            "directions": entries}


def summary_text(lib):
    """The case summary's section for a library."""
    ref = lib["reference_wind"]
    asked = ("for %.6g m/s at %.6g m" % (ref["speed"], ref["height"]) if ref["height"] is not None
             else "given (atmosphere.geostrophic_wind)")
    first = lib["directions"][0]
    roughness = sorted({e["roughness"] for e in lib["directions"]})
    z0_text = ("z0 = %.6g m" % roughness[0] if len(roughness) == 1
               else "z0 per direction %.6g .. %.6g m (wind.upstream_roughness = auto)" % (roughness[0], roughness[-1]))
    head = ("from", "|G| m/s", "G from", "turning", "U_ref m/s", "U_ref from", "stress turn", "u* m/s", "drag |G|",
            "drag turn")
    fmt = "  " + " ".join(["%10s"] * len(head))
    lines = [SUMMARY_MARK,
             "  the column's geostrophic wind G %s, %s, latitude %.6g" % (asked, z0_text, first["latitude"]),
             fmt % head]

    def num(x, f):
        return "-" if x is None else f % x

    for e in lib["directions"]:
        g, h, d = e["geostrophic_wind"], e["wind_at_height"], e["drag_law"] or {}
        lines.append(fmt % ("%.6g" % e["direction"], "%.6g" % g["speed"], "%.4f" % g["direction"],
                            num(e["turning_deg"], "%.4f"), num(h and h["speed"], "%.6g"),
                            num(h and h["direction"], "%.4f"), "%.4f" % e["turning_from_stress_deg"],
                            num(e["ustar"], "%.6g"), num(d.get("speed"), "%.6g"), num(d.get("turning_deg"), "%.4f")))
    lines.append("  turning: from U_ref to G; stress turn: from the wind at the first cell centre (the surface stress) "
                 "to G; degrees, counter-clockwise positive")
    lines.append("  drag: the neutral geostrophic drag law for the column's u* (A = 1.8, B = 4.5; its angle is from "
                 "the surface stress, so compare it with the stress turn); - for a stratified column")
    return "\n".join(lines) + "\n"


def _replace(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as out:
        out.write(text)
    os.replace(tmp, path)


def write(case_dir, lib):
    """Write library.json and the case summary's section; returns the library file's path."""
    path = os.path.join(case_dir, LIBRARY_FILE)
    _replace(path, json.dumps(lib, indent=2, allow_nan=False) + "\n")
    summary = os.path.join(case_dir, SUMMARY_FILE)
    text = ""
    if os.path.isfile(summary):
        with open(summary, encoding="utf-8") as fh:
            text = fh.read()
        if text.startswith(SUMMARY_MARK + "\n"):
            text = ""
        else:
            n = text.rfind("\n" + SUMMARY_MARK + "\n")
            if n >= 0:
                text = text[:n + 1]
        if text and not text.endswith("\n\n"):
            text = text.rstrip("\n") + "\n\n"
    _replace(summary, text + summary_text(lib))
    return path
