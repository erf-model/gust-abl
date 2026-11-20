"""The case description: its keys, their defaults and ranges, and the reader.

A case is a YAML file of sections (``domain``, ``mesh``, ``wind``, ...). Every
key is declared once in :data:`SPEC` with its kind, default, unit and range;
the reader, ``make_gust_case.py --list-keys`` and docs/case_setup_keys.inc all
come from that table.

The reader follows GUST's rules for its own inputs:

* an unknown section or key is refused, with the nearest known name;
* a key given twice in the YAML file is refused;
* a value outside its range is refused, never clamped;
* every message names the key as ``section.key``.

All problems of a file are collected and reported together.
"""
import copy
import difflib
import math

import textwrap

import yaml

from .errors import CaseError

AUTO = "auto"
NONE = "none"


class Key:
    """One key of the case description.

    kind: number, integer, flag, string, path, choice, numbers, strings or
    any (checked by the builder). ``words`` are the words a number or integer
    key also takes (``auto``, ``none``). ``lo``/``hi`` bound a number, an
    integer or every entry of a list; ``lo_open``/``hi_open`` exclude the
    bound. ``length`` fixes a list's length.
    """

    def __init__(self, kind, default=None, doc="", unit="", lo=None, hi=None, lo_open=False, hi_open=False,
                 choices=None, words=(), length=None, required=False, zero_ok=False):
        self.kind = kind
        self.zero_ok = zero_ok      # 0 is accepted besides the range (GUST's "0 (off), or >= lo")
        self.default = default
        self.doc = doc
        self.unit = unit
        self.lo = lo
        self.hi = hi
        self.lo_open = lo_open
        self.hi_open = hi_open
        self.choices = tuple(choices) if choices else ()
        self.words = tuple(words)
        self.length = length
        self.required = required

    def range_text(self):
        """The accepted values, in words (for messages and the key table)."""
        if self.kind == "choice":
            return ", ".join(self.choices)
        if self.kind == "flag":
            return "true or false"
        parts = []
        if self.kind in ("number", "integer", "numbers"):
            if self.lo is not None:
                parts.append(("> " if self.lo_open else "≥ ") + _fmt(self.lo))
            if self.hi is not None:
                parts.append(("< " if self.hi_open else "≤ ") + _fmt(self.hi))
        text = {"number": "a number", "integer": "an integer", "numbers": "a list of numbers",
                "strings": "a list of words", "string": "text", "path": "a file path", "any": "see the text"}[self.kind]
        if self.kind == "numbers" and self.length:
            text = "%d numbers" % self.length
        if parts:
            text += " " + " and ".join(parts)
        if self.zero_ok:
            text = "0, or " + text
        if self.words:
            text += ", or " + " or ".join(self.words)
        return text


def _fmt(value):
    return "%g" % value


class NestList:
    """A list of sections (``nests``): each entry is validated against ``spec``."""

    def __init__(self, spec, doc=""):
        self.spec = spec
        self.doc = doc


_EXTENT = dict(unit="m", lo=0.0, lo_open=True)

SPEC = {
    "case": {
        "name": Key("string", "gust_case", "Name of the case folder, written under output_dir: letters, digits, _, - "
                    "and ., not starting with a dot."),
        "output_dir": Key("path", ".", "Folder the case folder is made in; relative to the YAML file."),
        "description": Key("string", "", "One line copied to the top of every deck."),
    },
    "domain": {
        "center": Key("numbers", None, "Centre of the domain in the terrain file's projected frame (easting, "
                      "northing). Not with center_lat_lon. Absent (and no center_lat_lon): 0 0.", unit="m", length=2),
        "center_lat_lon": Key("numbers", None, "Centre of the domain as latitude, longitude (WGS 84); converted to "
                              "UTM. Required with a terrain file in latitude and longitude.", unit="degrees", length=2),
        "utm_zone": Key("integer", None, "UTM zone of the projection. Default: the zone of center_lat_lon.", lo=1, hi=60),
        "half_width": Key("numbers", None, "Distance from the centre to the domain's faces: one number, or two "
                          "(east–west, north–south). Not with west, east, south, north.", unit="m", lo=0.0, lo_open=True),
        "west": Key("number", None, "Distance from the centre to the west face.", **_EXTENT),
        "east": Key("number", None, "Distance from the centre to the east face.", **_EXTENT),
        "south": Key("number", None, "Distance from the centre to the south face.", **_EXTENT),
        "north": Key("number", None, "Distance from the centre to the north face.", **_EXTENT),
        "top": Key("number", None, "Height of the domain above the zero level of the terrain.", required=True, **_EXTENT),
        "origin": Key("choice", "center", "Where x = y = 0 of the case is: the domain centre or its south-west corner.",
                      choices=("center", "southwest")),
    },
    "terrain": {
        "file": Key("path", None, "Elevation data. Absent: flat ground."),
        "format": Key("choice", AUTO, "geotiff (single band), esri_ascii (.asc grid) or xyz (x y z points); auto "
                      "goes by the file's extension.", choices=(AUTO, "geotiff", "esri_ascii", "xyz")),
        "coordinates": Key("choice", AUTO, "projected (metres) or geographic (longitude, latitude); auto: "
                           "geographic when the file's extent fits in ±180, ±90.",
                           choices=(AUTO, "projected", "geographic")),
        "zero_level": Key("number", "lowest", "Elevation that becomes z = 0: the lowest point kept, none (the "
                          "file's own datum) or a number.", unit="m", words=("lowest", NONE)),
        "smoothing": Key("choice", "btf", "The terrain-following coordinate: btf, stf or sullivan "
                         "(terrain.smoothing 0, 1, 2).", choices=("btf", "stf", "sullivan")),
        "point_spacing": Key("number", AUTO, "Spacing of the points written to terrain.xyz. auto: the data's own, "
                             "thinned to no finer than half a cell.", unit="m", lo=0.0, lo_open=True, words=(AUTO,)),
        "margin": Key("number", AUTO, "Width of the band of points kept outside the domain. auto: 4 cells.",
                      unit="m", lo=0.0, words=(AUTO,)),
        "idw_neighbors": Key("integer", None, "terrain.idw_neighbors (GUST's default 6).", lo=1),
        "idw_power": Key("number", None, "terrain.idw_power (GUST's default 2).", lo=0.0, lo_open=True),
        "border": {
            "flat_width": Key("number", 0.0, "Width of the flat band inside each lateral face.", unit="m", lo=0.0),
            "slope_width": Key("number", 0.0, "Width over which the ground goes from the flat band to the data. "
                               "0 with flat_width 0: the border is left as the data give it.", unit="m", lo=0.0),
            "level": Key("number", "lowest", "Elevation of the flat band: the lowest point kept, mean (the mean "
                         "elevation) or a number in the file's datum.", unit="m", words=("lowest", "mean")),
        },
    },
    "mesh": {
        "cell_size": Key("number", None, "Target horizontal cell size. The cell counts are the nearest multiples "
                         "of 8, so the actual size differs slightly.", required=True, **_EXTENT),
        "layers": Key("integer", None, "Number of layers. Required unless first_layer and stretch_ratio are given.", lo=1),
        "first_layer": Key("number", None, "Thickness of the first layer. Absent: uniform layers.", **_EXTENT),
        "stretch_ratio": Key("number", None, "Ratio between successive layers. With first_layer and no layers: the "
                             "largest ratio allowed, and the layer count follows. With both: the ratio written to "
                             "the deck (the top is then where the layers end).", lo=1.0),
        "hierarchy_levels": Key("integer", AUTO, "Levels of the grid hierarchy. auto: up to 3, as the cell counts "
                                "allow.", lo=1, hi=10, words=(AUTO,)),
        "max_grid_size": Key("integer", 64, "Largest box size in x and y.", lo=1),
    },
    "wind": {
        "speed": Key("number", None, "Wind speed at the reference height.", unit="m/s", lo=0.0, lo_open=True),
        "height": Key("number", None, "Reference height above the ground.", unit="m", lo=0.0, lo_open=True),
        "direction": Key("number", None, "Direction the wind comes from, clockwise from north (270: from the west).",
                         unit="degrees", lo=0.0, hi=360.0),
        "directions": Key("any", None, "A sweep: a list of directions, or start, stop and step. One run folder per "
                          "direction. Not with direction."),
        "drive": Key("choice", "inflow_only", "What drives the flow (abl.drive): inflow_only, couette (a moving lid "
                     "at the log law's speed), geostrophic_column or geostrophic_inflow (Coriolis and a geostrophic "
                     "wind; the inflow profile comes from a column run made first).",
                     choices=("inflow_only", "couette", "geostrophic_column", "geostrophic_inflow")),
        "profile": Key("choice", "log_law", "The inflow profile: log_law, power_law or table (a file of "
                       "z u v theta [k] rows).", choices=("log_law", "power_law", "table")),
        "power_law_exponent": Key("number", None, "Exponent of the power law (GUST's default 0.14).", lo=0.0, lo_open=True),
        "table": Key("path", None, "The profile table; copied into the case folder. With profile: table."),
        "upstream_roughness": Key("number", None, "Roughness length of the inflow profile. Default: surface.roughness.",
                                  unit="m", lo=0.0, lo_open=True),
        "side_faces": Key("choice", "slip_wall", "The two faces parallel to a wind along x or y; not with a geostrophic drive.",
                          choices=("slip_wall", "periodic")),
    },
    "atmosphere": {
        "theta": Key("choice", "neutral", "How potential temperature is treated (abl.theta): neutral (buoyancy off), "
                     "capped (an inversion or lapse rate aloft, solved), frozen (the sounding held) or surface_flux.",
                     choices=("neutral", "capped", "frozen", "surface_flux")),
        "reference_temperature": Key("number", 300.0, "Potential temperature below the inversion.", unit="K",
                                     lo=0.0, lo_open=True),
        "inversion": {
            "base": Key("number", None, "Height of the inversion's base.", unit="m", lo=0.0, lo_open=True),
            "thickness": Key("number", None, "Thickness of the inversion.", unit="m", lo=0.0, lo_open=True),
            "strength": Key("number", None, "Rise of potential temperature across the inversion.", unit="K", lo=0.0),
        },
        "lapse_rate": Key("number", None, "Gradient of potential temperature above the inversion (or above the "
                          "ground without one); also the gradient held at the top.", unit="K/m", lo=0.0),
        "sounding": {
            "heights": Key("numbers", None, "Heights of a sounding, increasing from 0. Not with inversion.", unit="m", lo=0.0),
            "theta": Key("numbers", None, "Potential temperature at those heights.", unit="K", lo=0.0, lo_open=True),
        },
        "damping": {
            "start": Key("number", AUTO, "Height where the Rayleigh layer starts. auto: 0.7 of the top. none: no layer.",
                         unit="m", lo=0.0, lo_open=True, words=(AUTO, NONE)),
            "coefficient": Key("number", 0.2, "Rayleigh coefficient.", unit="1/s", lo=0.0, lo_open=True),
            "fields": Key("strings", None, "Fields damped. Default: w (w and theta with theta: capped).",
                          choices=("u", "v", "w", "theta", "tke")),
        },
        "latitude": Key("number", None, "Latitude of the Coriolis force. Default: that of domain.center_lat_lon.",
                        unit="degrees", lo=-90.0, hi=90.0),
        "geostrophic_wind": Key("numbers", None, "Geostrophic wind (east, north components). With a geostrophic "
                                "drive: given here, or found by the column run from wind.speed at wind.height.",
                                unit="m/s", length=2),
        "surface_heat_flux": Key("number", None, "Kinematic surface heat flux; with theta: surface_flux.", unit="K m/s"),
    },
    "surface": {
        "roughness": Key("number", None, "Uniform roughness length. Required unless roughness_map is given.", unit="m",
                         lo=0.0, lo_open=True),
        "roughness_map": {
            "file": Key("path", None, "A roughness map: x y z0 points, or a land-cover raster."),
            "format": Key("choice", AUTO, "xyz (x y z0 points, in the terrain file's frame) or land_cover (a "
                          "GeoTIFF or ESRI ASCII raster of class codes). auto goes by the extension.",
                          choices=(AUTO, "xyz", "land_cover")),
            "classes": Key("any", "esa_cci", "The class table of a land-cover raster: esa_cci, or a mapping of "
                           "class code to roughness length (m)."),
            "default": Key("number", None, "Roughness length of a class the table does not hold. Absent: such a "
                           "class is refused.", unit="m", lo=0.0, lo_open=True),
            "coordinates": Key("choice", AUTO, "As terrain.coordinates.", choices=(AUTO, "projected", "geographic")),
        },
        "canopy": {
            "patches_file": Key("path", None, "Round forest patches (ERF's forest file; canopy.type = patches). "
                                "Coordinates in the case frame."),
            "height_file": Key("path", None, "x y h points of canopy height (canopy.type = map), case frame."),
            "lai_file": Key("path", None, "x y LAI points; with height_file."),
            "cd": Key("number", None, "Drag coefficient of every column; with height_file.", lo=0.0, lo_open=True),
            "cd_file": Key("path", None, "x y cd points; with height_file, not with cd."),
            "profile": Key("choice", "uniform", "Shape of the leaf area density.", choices=("uniform", "lalic_mihailovic")),
        },
        "buildings": {
            "file": Key("path", None, "Building height map (buildings.type = immersed_forcing), case frame."),
            "format": Key("choice", "erf", "erf or xyz (buildings.format).", choices=("erf", "xyz")),
        },
    },
    "model": {
        "lee_setting": Key("choice", "recommended", "recommended: turb.realizability = 0.309, turb.nut_relax = 0.5 "
                           "and the convection blend (0.1 below blend_height). none: GUST's defaults.",
                           choices=("recommended", NONE)),
        "blend_height": Key("number", AUTO, "convection.hybrid_height. auto: 2.5 times the terrain's relief, "
                            "highest less lowest point in the domain, a regional slope included (no blend over flat "
                            "ground).", unit="m", lo=0.0, lo_open=True, words=(AUTO,)),
        "length_scale_cap": Key("number", AUTO, "turb.max_geom_lscale. auto: 100 m, GUST's default (1e6 with drive: couette): "
                                "over a 40 H fetch a 30 m cap lets the near-ground wind fall "
                                "12 %, 100 m holds it within 4 % (Tests/Real_Terrain/RUSHIL/README.md).",
                                unit="m", lo=0.0, lo_open=True, words=(AUTO,)),
        "tke_floor": Key("number", 1.0e-6, "turb.tke_floor; 0 is machine epsilon.", unit="m2/s2",
                         lo=2.220446049250313e-16, zero_ok=True),
    },
    "solver": {
        "max_iterations": Key("integer", 20000, "steady.max_iterations.", lo=1),
        "residual_tolerance": Key("number", None, "steady.residual_tol (GUST's default 5e-5).", lo=0.0, lo_open=True),
        "residual_drop_orders": Key("number", None, "steady.residual_drop_orders (GUST's default 4).", lo=0.0, hi=308.0),
        "cfl_max": Key("number", 16.0, "steady.cfl_max; at least steady.cfl0 (GUST's default 1, not written).", lo=1.0),
        "sweep": Key("choice", "lusgs", "steady.sweep: lusgs (CPU), zjacobi (also GPU) or auto (GUST chooses by "
                     "its build).", choices=("lusgs", "zjacobi", AUTO)),
        "grid_sequencing": Key("flag", None, "Start on the coarse levels (steady.sequencing). Default: on with "
                               "two or more hierarchy levels."),
        "w_vnn_cap": Key("number", 16.0, "simplec.w_vnn_cap; 0 switches the cap off.", lo=16.0, zero_ok=True),
        "initial_blend": Key("numbers", (1.0, 0.5), "init.alpha_h and init.alpha_v of the mass-consistent first "
                             "guess.", length=2, lo=0.0, lo_open=True),
    },
    "output": {
        "prefix": Key("string", "", "Prefix of the plotfile and checkpoint names (prefix_plt, ...)."),
        "turbulence_plotfile": Key("flag", True, "Write the turbulence plotfile (io.turb_plot_file)."),
        "checkpoint": Key("flag", False, "Write checkpoints. Always on with nests: a patch reads its parent's."),
    },
    "run": {
        "ranks": Key("integer", 1, "MPI ranks of the run script (run.sh uses run.mpi_launcher even for 1; the "
                     "window runs gust alone).", lo=1),
        "executable": Key("string", "gust", "The GUST executable: a name on PATH or a path. run.sh runs it from "
                          "each step's folder (column/, nest_*/ are one deeper), so give an absolute path there; the "
                          "case set-up window takes a relative path from the case file's folder."),
        "mpi_launcher": Key("string", "mpiexec -n {ranks}", "The MPI command; {ranks} is replaced."),
        "resume": Key("flag", True, "Skip a wind direction whose runs reached GUST's steady verdict with the "
                      "current decks, data files and run.sh arguments and wrote a plotfile after it, so a stopped sweep "
                      "goes on where it stopped. false runs every direction (RESUME=0 ./run.sh does the same once)."),
    },
    "nests": NestList({
        "name": Key("string", None, "Name of the patch; its folder is nest_<name>.", required=True),
        "center": Key("numbers", None, "Centre of the patch, as domain.center.", unit="m", length=2),
        "center_lat_lon": Key("numbers", None, "Centre of the patch as latitude, longitude.", unit="degrees", length=2),
        "half_width": Key("numbers", None, "Distance from the centre to the patch's faces (one or two numbers). The "
                          "faces are moved to the parent's nearest cell faces.", unit="m", lo=0.0, lo_open=True,
                          required=True),
        "blend_width": Key("number", AUTO, "terrain.nest_blend_width. auto: 4 parent cells.", unit="m", lo=0.0,
                           words=(AUTO,)),
        "max_iterations": Key("integer", None, "steady.max_iterations of the patch. Default: solver.max_iterations.", lo=1),
    }, "One-way nested patches at twice the horizontal resolution, each run after its parent."),
    "gust_inputs": Key("any", None, "GUST inputs written as they are (namespace.key: value) after everything else: "
                       "a key the builder also writes is replaced, and null removes one. They go to the 3-D decks "
                       "(the run folder's and each patch's), not to column/inputs, the 1-D column run with its own "
                       "mesh and keys."),
}

SECTION_DOC = {
    "case": "Name and place of the case folder.",
    "domain": "The box: centre, extent and top.",
    "terrain": "Elevation data and how it becomes GUST's terrain.xyz.",
    "mesh": "Cell size, layers and grid hierarchy.",
    "wind": "Speed, direction and what drives the flow.",
    "atmosphere": "Potential temperature, the Rayleigh layer and rotation.",
    "surface": "Roughness, canopy and buildings.",
    "model": "Closure settings.",
    "solver": "The steady iteration.",
    "output": "Plotfiles and checkpoints.",
    "run": "The run script.",
    "nests": "Nested patches.",
    "gust_inputs": "Any other GUST input.",
}


# ---------------------------------------------------------------------------------------------------------------
# Reading the YAML file
# ---------------------------------------------------------------------------------------------------------------
class _UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that refuses a mapping key given twice (YAML keeps the last silently)."""

    def construct_mapping(self, node, deep=False):
        seen = {}
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            if key in seen:
                raise CaseError("%s: line %d: the key '%s' is given twice (first on line %d)"
                                % (getattr(self, "name", "case file"), key_node.start_mark.line + 1, key, seen[key]))
            seen[key] = key_node.start_mark.line + 1
        return super().construct_mapping(node, deep=deep)


def load_yaml(path):
    """The YAML file as nested dicts; a syntax error or a repeated key is a CaseError."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            loader = _UniqueKeyLoader(handle)
            try:
                data = loader.get_single_data()
            finally:
                loader.dispose()
    except OSError as err:
        raise CaseError("%s: cannot read the case file (%s)" % (path, err.strerror or err)) from None
    except yaml.YAMLError as err:
        raise CaseError("%s: not valid YAML: %s" % (path, " ".join(str(err).split()))) from None
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise CaseError("%s: the case file must hold sections (case:, domain:, mesh:, ...)" % path)
    return data


def apply_override(raw, assignment):
    """Apply one ``section.key=value`` override to the raw case description.

    The value is read as YAML, so ``wind.directions=[0, 90]`` and
    ``atmosphere.damping.start=none`` work. ``gust_inputs.<namespace.key>=value``
    sets one GUST input.
    """
    if "=" not in assignment:
        raise CaseError("--set %s: expected section.key=value" % assignment)
    path, text = assignment.split("=", 1)
    path = path.strip()
    try:
        value = yaml.safe_load(text)
    except yaml.YAMLError as err:
        raise CaseError("--set %s: the value is not valid YAML (%s)" % (assignment, " ".join(str(err).split()))) from None
    if path.startswith("gust_inputs."):
        parts = ["gust_inputs", path[len("gust_inputs."):]]
    else:
        parts = path.split(".")
    if not all(parts):
        raise CaseError("--set %s: expected section.key=value" % assignment)
    node = raw
    for part in parts[:-1]:
        if node.get(part) is None:
            node[part] = {}
        node = node[part]
        if not isinstance(node, dict):
            raise CaseError("--set %s: %s is not a section" % (assignment, part))
    node[parts[-1]] = value


# ---------------------------------------------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------------------------------------------
def _number(value):
    """A float from a YAML scalar, or None. Accepts the text PyYAML leaves for 1e-6."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def _in_range(key, x):
    if not math.isfinite(x):
        return False
    if key.zero_ok and x == 0.0:
        return True
    if key.lo is not None and (x < key.lo or (key.lo_open and x == key.lo)):
        return False
    if key.hi is not None and (x > key.hi or (key.hi_open and x == key.hi)):
        return False
    return True


def _check_value(name, key, value, problems):
    """The normalized value of one key; appends to ``problems`` when refused."""
    def refuse():
        problems.append("%s = %s: must be %s" % (name, _show(value), key.range_text()))
        return None

    kind = key.kind
    if kind == "any":
        return value
    if kind in ("number", "integer"):
        if isinstance(value, str) and value in key.words:
            return value
        x = _number(value)
        if x is None or not _in_range(key, x):
            return refuse()
        if kind == "integer":
            if x != int(x):
                return refuse()
            return int(x)
        return x
    if kind == "flag":
        if not isinstance(value, bool):
            return refuse()
        return value
    if kind in ("string", "path"):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            value = str(value)
        if not isinstance(value, str) or (kind == "path" and not value.strip()):
            return refuse()
        return value
    if kind == "choice":
        if value not in key.choices:
            near = difflib.get_close_matches(str(value), key.choices, n=1)
            problems.append("%s = %s: must be %s%s" % (name, _show(value), key.range_text(),
                                                         " (did you mean %s?)" % near[0] if near else ""))
            return None
        return value
    if kind == "numbers":
        items = value if isinstance(value, (list, tuple)) else [value]
        out = [_number(v) for v in items]
        if not out or any(x is None or not _in_range(key, x) for x in out):
            return refuse()
        if key.length and len(out) != key.length:
            return refuse()
        return out
    if kind == "strings":
        items = value if isinstance(value, (list, tuple)) else str(value).split()
        if not items or any(not isinstance(v, str) for v in items):
            return refuse()
        if key.choices and any(v not in key.choices for v in items):
            problems.append("%s = %s: each must be one of %s" % (name, _show(value), ", ".join(key.choices)))
            return None
        return list(items)
    raise AssertionError("unknown kind " + kind)


def _show(value):
    if isinstance(value, str):
        return "'%s'" % value
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _validate_section(prefix, spec, raw, problems, given):
    out = {}
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        problems.append("%s: must be a section of keys, not %s" % (prefix.rstrip("."), _show(raw)))
        raw = {}
    for name in raw:
        if name not in spec:
            near = difflib.get_close_matches(str(name), list(spec), n=1)
            problems.append("%s%s: unknown key%s" % (prefix, name, " (did you mean %s%s?)" % (prefix, near[0]) if near else ""))
    for name, entry in spec.items():
        full = prefix + name
        if isinstance(entry, dict):
            out[name] = _validate_section(full + ".", entry, raw.get(name), problems, given)
        elif isinstance(entry, NestList):
            items = raw.get(name)
            out[name] = []
            if items is None:
                continue
            given.add(full)
            if not isinstance(items, list):
                problems.append("%s: must be a list of patches (- name: ...)" % full)
                continue
            for i, item in enumerate(items):
                out[name].append(_validate_section("%s[%d]." % (full, i), entry.spec, item, problems, given))
        elif name in raw and raw[name] is not None:
            given.add(full)
            out[name] = _check_value(full, entry, raw[name], problems)
        else:
            if entry.required:
                problems.append("%s: required (%s)" % (full, entry.range_text()))
            out[name] = copy.deepcopy(entry.default)
            if isinstance(out[name], tuple):
                out[name] = list(out[name])
    return out


def validate(raw):
    """The normalized case description (defaults filled in) and the set of keys the user gave.

    Raises CaseError listing every unknown key and every value out of range.
    Rules that relate several keys are the builder's.
    """
    problems = []
    given = set()
    case = _validate_section("", SPEC, raw, problems, given)
    if problems:
        raise CaseError(problems)
    return case, given


# ---------------------------------------------------------------------------------------------------------------
# The key table
# ---------------------------------------------------------------------------------------------------------------
def iter_keys(spec=None, prefix=""):
    """(full name, Key) of every key, in the table's order."""
    for name, entry in (SPEC if spec is None else spec).items():
        full = prefix + name
        if isinstance(entry, dict):
            yield from iter_keys(entry, full + ".")
        elif isinstance(entry, NestList):
            yield from iter_keys(entry.spec, full + "[].")
        else:
            yield full, entry


def _default_text(key):
    if key.required:
        return "required"
    if key.default is None:
        return "none"
    if isinstance(key.default, bool):
        return "true" if key.default else "false"
    if isinstance(key.default, (list, tuple)):
        return " ".join("%g" % v for v in key.default)
    if isinstance(key.default, float):
        return "%g" % key.default
    return str(key.default) if key.default != "" else "empty"


def key_table_text():
    """The keys as plain text (``make_gust_case.py --list-keys``)."""
    lines = []
    section = None
    for name, key in iter_keys():
        top = name.split(".")[0].replace("[]", "")
        if top != section:
            section = top
            lines += ["", "%s: %s" % (section, SECTION_DOC[section])]
        unit = " [%s]" % key.unit if key.unit else ""
        lines.append("  %s%s = %s" % (name, unit, _default_text(key)))
        lines.append("      %s Accepts %s." % (key.doc, key.range_text()))
    return "\n".join(lines[1:]) + "\n"


def key_table_rst():
    """The keys as reStructuredText list-tables, one per section (docs/case_setup_keys.inc)."""
    out = [".. Generated by tools/make_gust_case.py --write-key-table; do not edit.", ""]
    section = None
    for name, key in iter_keys():
        top = name.split(".")[0].replace("[]", "")
        if top != section:
            section = top
            out += ["``%s``: %s" % (section, SECTION_DOC[section]), "",
                    ".. list-table::", "   :header-rows: 1", "   :widths: 26 12 14 48", "",
                    "   * - Key", "     - Unit", "     - Default", "     - Meaning and accepted values"]
        out += ["   * - ``%s``" % name, "     - %s" % (key.unit or "–"), "     - %s" % _default_text(key),
                "     - %s Accepts %s." % (key.doc, key.range_text())]
        nxt = None
        # a blank line closes each table: detect the last key of the section
        names = [n for n, _ in iter_keys()]
        i = names.index(name)
        if i + 1 < len(names):
            nxt = names[i + 1].split(".")[0].replace("[]", "")
        if nxt != section:
            out.append("")
    return "\n".join(out)


def _yaml_default(key):
    """(value text, note) of a key line of the YAML reference: the default as YAML, or empty with a note."""
    if key.required:
        return "", "required"
    if key.default is None:
        return "", "none"
    if isinstance(key.default, bool):
        return ("true" if key.default else "false"), ""
    if isinstance(key.default, (list, tuple)):
        return "[%s]" % ", ".join("%g" % v for v in key.default), ""
    if isinstance(key.default, float):
        text = "%g" % key.default
        if "e" in text and "." not in text:      # YAML reads 1e-06 as text; 1.0e-06 is a number
            text = text.replace("e", ".0e")
        return text, ""
    return (str(key.default) if key.default != "" else '""'), ""


YAML_REFERENCE_HEADER = """\
## GUST master YAML input reference (every key of the case description read by
## tools/make_gust_case.py, with defaults, units, ranges and meaning).
##
## Generated by tools/make_gust_case.py --write-yaml-reference; do not edit by
## hand. The CTest gust_case_selftest checks that this file is current and
## holds every key.
##
## Every line is a comment. Copy the lines you need into a case file and remove
## the leading "# " (lines that start with "##" stay comments); what is left is
## YAML, with the text below each key still a comment. A key line reads
## "key: default". An empty value with "# required" means the key has no
## default; with "# none", that the key is absent unless given. A unit follows
## in brackets. Not every key can be given with every other: the text says
## which exclude each other, and the tool refuses a case that mixes them. docs/case_setup.rst describes the sections;
## Tests/inputs_reference lists GUST's own deck inputs, which the gust_inputs
## section passes through.
"""


def key_reference_yaml(width=116):
    """Every key as a commented YAML file (Tests/inputs_yaml_reference)."""
    out = [YAML_REFERENCE_HEADER.rstrip("\n")]
    bar = "##" + "=" * (width - 2)
    open_path = []
    for name, key in iter_keys():
        parts = name.split(".")
        top = parts[0].replace("[]", "")
        in_list = parts[0].endswith("[]")
        parents = [top] + parts[1:-1]
        if not open_path or open_path[0] != top:
            out += ["", "", bar, "## %s: %s" % (top, SECTION_DOC[top]), bar, "#"]
            open_path = []
            first_in_list = in_list
            if len(parts) > 1:
                out.append("# %s:" % top)
                open_path = [top]
        else:
            first_in_list = False
            out.append("#")
        depth0 = 2 if in_list else 0          # list entries sit under "- "
        for level in range(1, len(parents)):
            if open_path[:level + 1] != parents[:level + 1]:
                out.append("# %s%s:" % (" " * (depth0 + 2 * level), parents[level]))
        open_path = parents
        indent = " " * (depth0 + 2 * (len(parts) - 1))
        if first_in_list:
            indent = indent[:-2] + "- "
        value, note = _yaml_default(key)
        notes = [n for n in (note, "[%s]" % key.unit if key.unit else "") if n]
        line = "# %s%s:%s" % (indent, parts[-1], " " + value if value else "")
        if notes:
            line += "    # " + " ".join(notes)
        out.append(line)
        pad = " " * (len(indent) + 4)
        text = "%s Accepts %s." % (key.doc, key.range_text())
        for piece in textwrap.wrap(text, width - len(pad) - 4):
            out.append("# %s# %s" % (pad, piece))
    return "\n".join(out) + "\n"
