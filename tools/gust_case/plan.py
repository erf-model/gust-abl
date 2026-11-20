"""From a validated case description to the numbers the decks need.

:func:`resolve` applies the rules that relate several keys (extent given one
way only, a moving lid needs the log law, ...) and works out the mesh, the
sounding and the patch boxes. Like the reader, it refuses instead of guessing,
names the keys, and reports every problem it finds together.
"""
import math
import os
import re

from . import atmosphere, fetch, geo, grid, wind
from .errors import CaseError
from .terrain import Frame

MARGIN_CELLS = 4          # terrain points kept outside the domain, in cells
NEST_KEEP_CELLS = 2       # parent cells between a patch and the parent's faces (bc.nest_parent)
NEST_BLEND_CELLS = 4      # default terrain.nest_blend_width, in parent cells
DAMPING_START_FRACTION = 0.7
DEFAULT_LENGTH_SCALE_CAP = 100.0
COUETTE_LENGTH_SCALE_CAP = 1.0e6
UNIFORM_TOLERANCE = 1.0e-6  # a stretch ratio this close to 1 is written as uniform layers
TOP_TOLERANCE = 0.01      # allowed relative difference between domain.top and where the layers end


class Patch:
    """A nested patch: its box in the case frame and its cell counts."""

    def __init__(self, name, x_lo, x_hi, y_lo, y_hi, nx, ny, blend_width, max_iterations):
        self.name = name
        self.x_lo, self.x_hi, self.y_lo, self.y_hi = x_lo, x_hi, y_lo, y_hi
        self.nx, self.ny = nx, ny
        self.blend_width = blend_width
        self.max_iterations = max_iterations

    def frame(self):
        f = Frame((0.0, 0.0), 1.0, 1.0, 1.0, 1.0)
        f.x_lo, f.x_hi, f.y_lo, f.y_hi = self.x_lo, self.x_hi, self.y_lo, self.y_hi
        return f


class Plan:
    """Everything resolved from the case description (plain attributes)."""


def _path(base_dir, value):
    return value if os.path.isabs(value) else os.path.normpath(os.path.join(base_dir, value))


def _resolve_domain(case, given, problems):
    d = case["domain"]
    sides = ("west", "east", "south", "north")
    have_sides = [s for s in sides if d[s] is not None]
    if d["half_width"] is not None and have_sides:
        problems.append("domain.half_width and domain.%s are both given: use one way to give the extent" % have_sides[0])
    extent = None
    if d["half_width"] is not None:
        hw = d["half_width"]
        if len(hw) > 2:
            problems.append("domain.half_width = %s: must be one number, or two (east–west, north–south)" % hw)
        else:
            hx, hy = (hw[0], hw[0]) if len(hw) == 1 else hw
            extent = (hx, hx, hy, hy)
    elif len(have_sides) == 4:
        extent = tuple(d[s] for s in sides)
    else:
        missing = [s for s in sides if d[s] is None]
        problems.append("domain: the extent is required: give domain.half_width, or all of domain.west, east, south "
                        "and north (missing: %s)" % ", ".join(missing))
    if d["center"] is not None and d["center_lat_lon"] is not None:
        problems.append("domain.center and domain.center_lat_lon are both given: use one")
    if d["utm_zone"] is not None and d["center_lat_lon"] is None:
        problems.append("domain.utm_zone has no effect without domain.center_lat_lon")
    lat_lon = zone = south = None
    center = (0.0, 0.0)
    if d["center_lat_lon"] is not None:
        lat, lon = d["center_lat_lon"]
        if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            problems.append("domain.center_lat_lon = %g %g: latitude must be in [-90, 90] and longitude in [-180, 180] "
                            "(the order is latitude, longitude)" % (lat, lon))
        else:
            east, north, zone, south = geo.latlon_to_utm(lat, lon, zone=d["utm_zone"])
            center, lat_lon = (float(east), float(north)), (lat, lon)
    elif d["center"] is not None:
        center = tuple(d["center"])
    if extent is None:
        return None
    return Frame(center, *extent, origin_at=d["origin"], lat_lon=lat_lon, zone=zone, south_hemisphere=south)


def _resolve_mesh(case, frame, problems, plan):
    m = case["mesh"]
    top = case["domain"]["top"]
    plan.nx = grid.horizontal_cells(frame.x_hi - frame.x_lo, m["cell_size"])
    plan.ny = grid.horizontal_cells(frame.y_hi - frame.y_lo, m["cell_size"])
    plan.dx = (frame.x_hi - frame.x_lo) / plan.nx
    plan.dy = (frame.y_hi - frame.y_lo) / plan.ny
    plan.max_grid_size = m["max_grid_size"]
    dz0, ratio, layers = m["first_layer"], m["stretch_ratio"], m["layers"]
    plan.write_stretch_ratio = False
    plan.top = top
    if dz0 is None:
        if ratio is not None:
            problems.append("mesh.stretch_ratio has no effect without mesh.first_layer")
        if layers is None:
            problems.append("mesh.layers: required (uniform layers); or give mesh.first_layer and mesh.stretch_ratio")
            layers = 1
        plan.layers, plan.dz0, plan.stretch_ratio = layers, None, 1.0
    elif layers is None:
        if ratio is None:
            problems.append("mesh.first_layer = %g: give mesh.layers (GUST finds the ratio that reaches domain.top) or "
                            "mesh.stretch_ratio (the layer count follows)" % dz0)
            ratio = 1.0
        if dz0 >= top:
            problems.append("mesh.first_layer = %g: must be below domain.top = %g" % (dz0, top))
            dz0 = top
        plan.layers = grid.layers_for_ratio(dz0, ratio, top)
        if dz0 * plan.layers > top * (1.0 + 1.0e-12):
            problems.append("mesh.first_layer = %g with mesh.stretch_ratio = %g: %d layers are needed to reach "
                            "domain.top = %g, and they pass it even without stretching (%g m); give mesh.layers, "
                            "a larger mesh.stretch_ratio or a mesh.first_layer that divides domain.top"
                            % (dz0, ratio, plan.layers, top, dz0 * plan.layers))
            plan.dz0, plan.stretch_ratio = dz0, 1.0
        else:
            _stretch_to_top(plan, dz0, top)
    else:
        if dz0 * layers > top * (1.0 + 1.0e-12):
            problems.append("mesh.first_layer = %g with mesh.layers = %d: the layers pass domain.top = %g even without "
                            "stretching (%g m); use fewer or thinner layers" % (dz0, layers, top, dz0 * layers))
            plan.layers, plan.dz0, plan.stretch_ratio = layers, dz0, 1.0
        elif ratio is None:
            plan.layers = layers
            _stretch_to_top(plan, dz0, top)
        else:
            plan.layers, plan.dz0, plan.stretch_ratio = layers, dz0, ratio
            plan.write_stretch_ratio = True
            plan.top = grid.stretched_top(dz0, ratio, layers)
            if abs(plan.top - top) > TOP_TOLERANCE * top:
                problems.append("mesh.first_layer = %g, mesh.stretch_ratio = %g and mesh.layers = %d end at %.6g m, not "
                                "at domain.top = %g m; give two of the three, or domain.top = %.6g"
                                % (dz0, ratio, layers, plan.top, top, plan.top))
    levels = m["hierarchy_levels"]
    if levels == "auto":
        plan.hierarchy_levels = grid.hierarchy_levels(plan.nx, plan.ny, plan.max_grid_size)
    else:
        f = 2 ** (levels - 1)
        if plan.nx % f or plan.ny % f or plan.max_grid_size % f or plan.nx // f < 2 or plan.ny // f < 2:
            problems.append("mesh.hierarchy_levels = %d: needs %d x %d cells and mesh.max_grid_size = %d divisible by "
                            "%d, with more than one coarse cell; use fewer levels or auto"
                            % (levels, plan.nx, plan.ny, plan.max_grid_size, f))
        plan.hierarchy_levels = levels
    seq = case["solver"]["grid_sequencing"]
    if seq is None:
        seq = plan.hierarchy_levels >= 2
    elif seq and plan.hierarchy_levels < 2:
        problems.append("solver.grid_sequencing = true needs mesh.hierarchy_levels of 2 or more (it is %d)"
                        % plan.hierarchy_levels)
    plan.grid_sequencing = seq


def _resolve_wind(case, given, problems, plan, base_dir):
    w = case["wind"]
    a = case["atmosphere"]
    s = case["surface"]
    drive = w["drive"]
    geostrophic = drive in ("geostrophic_column", "geostrophic_inflow")
    if w["direction"] is not None and w["directions"] is not None:
        problems.append("wind.direction and wind.directions are both given: use one")
    plan.sweep = w["directions"] is not None
    directions = []
    if w["directions"] is not None:
        spec = w["directions"]
        try:
            if isinstance(spec, dict):
                unknown = set(spec) - {"start", "stop", "step"}
                if unknown or set(spec) != {"start", "stop", "step"}:
                    raise ValueError("a sweep needs start, stop and step")
                if not float(spec["step"]) > 0.0 or float(spec["stop"]) < float(spec["start"]):
                    raise ValueError("step must be > 0 and stop not below start")
            elif not isinstance(spec, list) or not spec:
                raise ValueError("give a list, or start, stop and step")
            directions = wind.sweep_directions(spec)
            if any(not (math.isfinite(v) and 0.0 <= v <= 360.0) for v in directions):
                raise ValueError("every direction must be in [0, 360]")
            labels = [wind.direction_label(v) for v in directions]
            if len(set(labels)) != len(labels):
                raise ValueError("two directions have the same folder label dir_<deg> (labels are to 0.001 degree; "
                                 "0 and 360 are the same)")
        except (TypeError, ValueError, KeyError) as err:
            problems.append("wind.directions = %s: %s" % (spec, err))
            directions = []
    elif w["direction"] is not None:
        directions = [w["direction"]]
    else:
        problems.append("wind.direction: required (or wind.directions for a sweep)")
    plan.directions = directions

    # roughness: uniform, or a map with the inflow's own value
    rmap = s["roughness_map"]
    plan.roughness_map = rmap["file"] is not None
    if plan.roughness_map and s["roughness"] is not None:
        problems.append("surface.roughness and surface.roughness_map.file are both given: use one")
    for name in ("format", "classes", "default", "coordinates"):
        if "surface.roughness_map." + name in given and not plan.roughness_map:
            problems.append("surface.roughness_map.%s has no effect without surface.roughness_map.file" % name)
    if not plan.roughness_map and s["roughness"] is None:
        problems.append("surface.roughness: required (or surface.roughness_map.file)")
    plan.surface_roughness = s["roughness"]
    plan.upstream_auto = w["upstream_roughness"] == "auto"
    regional = s["regional_roughness_map"]
    plan.regional_map = regional["file"] is not None
    plan.fetch_length = None
    if plan.upstream_auto:
        z0_up = None              # per direction, from the regional map (fetch.py, Builder.prepare_data)
        if not plan.regional_map:
            problems.append("wind.upstream_roughness = auto needs surface.regional_roughness_map.file (the map "
                            "upwind of the domain)")
        if w["fetch_length"] is not None:
            plan.fetch_length = w["fetch_length"]
        elif w["height"] is not None:
            plan.fetch_length = fetch.FETCH_PER_HEIGHT * w["height"]
        else:
            problems.append("wind.fetch_length: required with wind.upstream_roughness = auto when wind.height is not "
                            "given (its default is %g times wind.height)" % fetch.FETCH_PER_HEIGHT)
    else:
        z0_up = w["upstream_roughness"] if w["upstream_roughness"] is not None else s["roughness"]
        if z0_up is None and plan.roughness_map:
            problems.append("wind.upstream_roughness: required with surface.roughness_map (the inflow profile's "
                            "roughness), or auto with surface.regional_roughness_map")
        for name in ("file", "format", "classes", "default", "coordinates"):
            if "surface.regional_roughness_map." + name in given:
                problems.append("surface.regional_roughness_map.%s has no effect without wind.upstream_roughness = auto"
                                % name)
        if "wind.fetch_length" in given:
            problems.append("wind.fetch_length has no effect without wind.upstream_roughness = auto")
    plan.upstream_roughness = z0_up

    # the profile
    profile = w["profile"]
    plan.profile = profile
    plan.table = None
    if profile == "table":
        if w["table"] is None:
            problems.append("wind.table: required with wind.profile = table")
        else:
            plan.table = _path(base_dir, w["table"])
            if not os.path.isfile(plan.table):
                problems.append("wind.table = %s: no such file" % plan.table)
        for name in ("speed", "height"):
            if w[name] is not None:
                problems.append("wind.%s has no effect with wind.profile = table (the table gives the wind)" % name)
        if geostrophic:
            problems.append("wind.profile = table: not with wind.drive = %s (the column run makes the table)" % drive)
    elif w["table"] is not None:
        problems.append("wind.table has no effect without wind.profile = table")
    if w["power_law_exponent"] is not None and profile != "power_law":
        problems.append("wind.power_law_exponent has no effect without wind.profile = power_law")

    plan.geostrophic_wind = a["geostrophic_wind"]
    plan.latitude = None
    if geostrophic:
        lat = a["latitude"]
        if lat is None and case["domain"]["center_lat_lon"] is not None:
            lat = case["domain"]["center_lat_lon"][0]
        if lat is None or lat == 0.0:
            problems.append("atmosphere.latitude: required and not 0 with wind.drive = %s (or give "
                            "domain.center_lat_lon)" % drive)
        plan.latitude = lat
        if profile == "power_law":
            problems.append("wind.profile = power_law: not with wind.drive = %s (the column run gives the profile)" % drive)
        if plan.geostrophic_wind is not None:
            for name in ("speed", "height"):
                if w[name] is not None:
                    problems.append("wind.%s and atmosphere.geostrophic_wind are both given: the column run takes the "
                                    "geostrophic wind, or finds it from wind.speed at wind.height" % name)
            if plan.sweep:
                problems.append("wind.directions: not with atmosphere.geostrophic_wind (one geostrophic wind has one "
                                "direction); give wind.speed and wind.height for a sweep")
        else:
            for name in ("speed", "height"):
                if w[name] is None:
                    problems.append("wind.%s: required with wind.drive = %s (the wind the column run is to reach), "
                                    "unless atmosphere.geostrophic_wind is given" % (name, drive))
        if drive == "geostrophic_column" and plan.roughness_map:
            problems.append("surface.roughness_map: not with wind.drive = geostrophic_column (GUST's column first "
                            "guess takes one roughness); use wind.drive = geostrophic_inflow")
        if "wind.side_faces" in given:
            problems.append("wind.side_faces has no effect with wind.drive = %s: the wind turns with height, so every "
                            "lateral face is an inflow or an outflow face" % drive)
    else:
        if a["latitude"] is not None:
            problems.append("atmosphere.latitude has no effect without a geostrophic wind.drive")
        if a["geostrophic_wind"] is not None:
            problems.append("atmosphere.geostrophic_wind has no effect without a geostrophic wind.drive")
        if profile != "table":
            for name in ("speed", "height"):
                if w[name] is None:
                    problems.append("wind.%s: required with wind.profile = %s" % (name, profile))
    plan.drive = drive
    plan.geostrophic = geostrophic
    plan.speed, plan.height = w["speed"], w["height"]
    plan.side_faces = w["side_faces"]
    if "wind.side_faces" in given and directions and not geostrophic and not any(
            0.0 in wind.unit_vector(wind.towards_from_meteorological(d)) for d in directions):
        problems.append("wind.side_faces has no effect: no wind.direction is along x or y, so every lateral face is "
                        "an inflow or an outflow face")
    plan.power_law_exponent = w["power_law_exponent"]

    if drive == "couette":
        if profile != "log_law":
            problems.append("wind.drive = couette needs wind.profile = log_law (the lid moves at the log law's speed)")
        if plan.roughness_map:
            problems.append("wind.drive = couette needs a uniform surface.roughness, not surface.roughness_map")
        elif plan.upstream_auto:
            problems.append("wind.upstream_roughness = auto: not with wind.drive = couette (the lid's log law takes "
                            "surface.roughness)")
        elif z0_up is not None and s["roughness"] is not None and z0_up != s["roughness"]:
            problems.append("wind.drive = couette needs wind.upstream_roughness = surface.roughness (%g and %g)"
                            % (z0_up, s["roughness"]))
        if a["theta"] != "neutral":
            problems.append("wind.drive = couette needs atmosphere.theta = neutral (it is %s)" % a["theta"])
    if plan.upstream_auto and profile != "log_law":
        problems.append("wind.upstream_roughness = auto: not with wind.profile = %s (the profile's wind is the same "
                        "for every direction; bc.inflow_z0 would set only its k)" % profile)
    plan.ustar = None
    if profile == "log_law" and not geostrophic and plan.speed and plan.height and z0_up:
        plan.ustar = wind.friction_velocity(plan.speed, plan.height, z0_up)


def _resolve_atmosphere(case, given, problems, plan):
    a = case["atmosphere"]
    theta = a["theta"]
    top = plan.top
    inv = a["inversion"]
    snd = a["sounding"]
    inv_given = [k for k in ("base", "thickness", "strength") if inv[k] is not None]
    snd_given = [k for k in ("heights", "theta") if snd[k] is not None]
    plan.theta = theta
    heights, values = [0.0], [a["reference_temperature"]]
    try:
        if snd_given:
            if len(snd_given) != 2:
                raise CaseError("atmosphere.sounding: give both heights and theta")
            if inv_given or a["lapse_rate"] is not None:
                raise CaseError("atmosphere.sounding: not with atmosphere.inversion or atmosphere.lapse_rate (the "
                                "sounding gives the whole profile)")
            if "atmosphere.reference_temperature" in given:
                raise CaseError("atmosphere.reference_temperature has no effect with atmosphere.sounding")
            atmosphere.check_sounding(snd["heights"], snd["theta"], top)
            heights, values = list(snd["heights"]), list(snd["theta"])
        elif inv_given or a["lapse_rate"] is not None:
            if inv_given and len(inv_given) != 3:
                raise CaseError("atmosphere.inversion: give base, thickness and strength (missing: %s)"
                                % ", ".join(k for k in ("base", "thickness", "strength") if inv[k] is None))
            heights, values = atmosphere.inversion_sounding(a["reference_temperature"], top, inv["base"], inv["thickness"],
                                                            inv["strength"], a["lapse_rate"])
    except CaseError as err:
        problems.extend(err.problems)
    plan.theta_heights, plan.theta_values = heights, values
    uniform = atmosphere.is_uniform(values)
    plan.top_theta_gradient = atmosphere.gradient_at_top(heights, values)
    if theta == "neutral" and not uniform:
        problems.append("atmosphere.theta = neutral: not with an inversion, a lapse rate or a sounding that varies; "
                        "use atmosphere.theta = capped or frozen")
    if theta in ("capped", "frozen") and uniform:
        problems.append("atmosphere.theta = %s needs a stratified profile: give atmosphere.inversion, "
                        "atmosphere.lapse_rate or atmosphere.sounding" % theta)
    if theta == "capped":
        if not atmosphere.never_falls(values):
            problems.append("atmosphere.theta = capped needs a sounding that never falls with height")
        if not plan.top_theta_gradient > 0.0:
            problems.append("atmosphere.theta = capped needs theta rising at the top of the domain (a lapse rate > 0, "
                            "or a sounding whose last segment rises)")
    if theta == "surface_flux":
        if a["surface_heat_flux"] is None or a["surface_heat_flux"] == 0.0:
            problems.append("atmosphere.surface_heat_flux: required and not 0 with atmosphere.theta = surface_flux")
    elif a["surface_heat_flux"] is not None:
        problems.append("atmosphere.surface_heat_flux has no effect without atmosphere.theta = surface_flux")
    plan.surface_heat_flux = a["surface_heat_flux"]

    damping = a["damping"]
    start = damping["start"]
    plan.damping_start = None
    if start == "none":
        for name in ("coefficient", "fields"):
            if "atmosphere.damping." + name in given:
                problems.append("atmosphere.damping.%s has no effect with atmosphere.damping.start = none" % name)
        if theta == "capped":
            problems.append("atmosphere.theta = capped needs the Rayleigh layer (the heat through the top needs a "
                            "sink): atmosphere.damping.start must not be none")
    else:
        plan.damping_start = DAMPING_START_FRACTION * top if start == "auto" else start
        if plan.damping_start >= top:
            problems.append("atmosphere.damping.start = %g: must be below the top of the domain (%g m)" % (start, top))
    fields = damping["fields"]
    if fields is None:
        fields = ["w", "theta"] if theta == "capped" else ["w"]
    elif theta == "capped" and not {"w", "theta"} <= set(fields):
        problems.append("atmosphere.theta = capped needs w and theta in atmosphere.damping.fields")
    if plan.geostrophic and plan.geostrophic_wind is None and {"u", "v"} & set(fields) and plan.damping_start is not None:
        problems.append("atmosphere.damping.fields: u and v are refused when the column run finds the geostrophic "
                        "wind from wind.speed (the layer would hold the wind to another profile)")
    plan.damping_fields = fields
    plan.damping_coefficient = damping["coefficient"]


def _resolve_surface(case, given, problems, plan, base_dir):
    s = case["surface"]
    c = s["canopy"]
    plan.canopy = None
    files = {}
    if c["patches_file"] is not None:
        if c["height_file"] is not None:
            problems.append("surface.canopy.patches_file and surface.canopy.height_file are both given: use one")
        for name in ("lai_file", "cd", "cd_file"):
            if c[name] is not None:
                problems.append("surface.canopy.%s has no effect with surface.canopy.patches_file (each patch line "
                                "holds its own)" % name)
        plan.canopy = "patches"
        files["canopy_patches.txt"] = ("surface.canopy.patches_file", c["patches_file"])
    elif c["height_file"] is not None:
        plan.canopy = "map"
        files["canopy_height.xyz"] = ("surface.canopy.height_file", c["height_file"])
        if c["lai_file"] is None:
            problems.append("surface.canopy.lai_file: required with surface.canopy.height_file")
        else:
            files["canopy_lai.xyz"] = ("surface.canopy.lai_file", c["lai_file"])
        if (c["cd"] is None) == (c["cd_file"] is None):
            problems.append("surface.canopy: give surface.canopy.cd or surface.canopy.cd_file, one of the two")
        elif c["cd_file"] is not None:
            files["canopy_cd.xyz"] = ("surface.canopy.cd_file", c["cd_file"])
    else:
        for name in ("lai_file", "cd", "cd_file"):
            if c[name] is not None:
                problems.append("surface.canopy.%s has no effect without surface.canopy.height_file" % name)
        if "surface.canopy.profile" in given:
            problems.append("surface.canopy.profile has no effect without surface.canopy.height_file")
    plan.canopy_cd, plan.canopy_profile = c["cd"], c["profile"]
    b = s["buildings"]
    plan.buildings = b["file"] is not None
    plan.buildings_format = b["format"]
    if plan.buildings:
        files["buildings." + ("xyz" if b["format"] == "xyz" else "txt")] = ("surface.buildings.file", b["file"])
    elif "surface.buildings.format" in given:
        problems.append("surface.buildings.format has no effect without surface.buildings.file")
    if plan.table is not None:
        files["inflow_table.dat"] = ("wind.table", plan.table)
    plan.copied_files = {}
    for target, (key, source) in files.items():
        source = _path(base_dir, source)
        if not os.path.isfile(source):
            problems.append("%s = %s: no such file" % (key, source))
        plan.copied_files[target] = source


def _resolve_nests(case, problems, plan, frame):
    plan.patches = []
    names = set()
    for i, n in enumerate(case["nests"]):
        tag = "nests[%d]" % i
        name = n["name"] or "patch%d" % i
        if not name.replace("_", "").isalnum():
            problems.append("%s.name = '%s': letters, digits and _ only (it names a folder)" % (tag, name))
        if name in names:
            problems.append("%s.name = '%s': given twice" % (tag, name))
        names.add(name)
        if (n["center"] is None) == (n["center_lat_lon"] is None):
            problems.append("%s: give %s.center or %s.center_lat_lon, one of the two" % (tag, tag, tag))
            continue
        if n["center_lat_lon"] is not None:
            if frame.lat_lon is None:
                problems.append("%s.center_lat_lon needs domain.center_lat_lon (the projection's zone)" % tag)
                continue
            cx, cy = frame.from_lat_lon(n["center_lat_lon"][0], n["center_lat_lon"][1])
            cx, cy = float(cx), float(cy)
        else:
            cx, cy = (float(v) for v in frame.to_case(n["center"][0], n["center"][1]))
        hw = n["half_width"] or [0.0]
        if len(hw) > 2:
            problems.append("%s.half_width = %s: must be one number, or two" % (tag, hw))
            continue
        hx, hy = (hw[0], hw[0]) if len(hw) == 1 else hw
        i_lo = int(round((cx - hx - frame.x_lo) / plan.dx))
        i_hi = int(round((cx + hx - frame.x_lo) / plan.dx))
        j_lo = int(round((cy - hy - frame.y_lo) / plan.dy))
        j_hi = int(round((cy + hy - frame.y_lo) / plan.dy))
        if (i_lo < NEST_KEEP_CELLS or j_lo < NEST_KEEP_CELLS or i_hi > plan.nx - NEST_KEEP_CELLS
                or j_hi > plan.ny - NEST_KEEP_CELLS):
            problems.append("%s ('%s'): the patch spans the parent's cells %d..%d in x and %d..%d in y; it must keep %d "
                            "cells from the parent's faces (%d x %d cells of %.6g x %.6g m)"
                            % (tag, name, i_lo, i_hi, j_lo, j_hi, NEST_KEEP_CELLS, plan.nx, plan.ny, plan.dx, plan.dy))
            continue
        if i_hi - i_lo < 4 or j_hi - j_lo < 4:
            problems.append("%s ('%s'): the patch covers %d x %d parent cells; it needs at least 4 x 4 (%s.half_width "
                            "against cells of %.6g x %.6g m)" % (tag, name, i_hi - i_lo, j_hi - j_lo, tag, plan.dx, plan.dy))
            continue
        x_lo, x_hi = frame.x_lo + i_lo * plan.dx, frame.x_lo + i_hi * plan.dx
        y_lo, y_hi = frame.y_lo + j_lo * plan.dy, frame.y_lo + j_hi * plan.dy
        blend = NEST_BLEND_CELLS * max(plan.dx, plan.dy) if n["blend_width"] == "auto" else n["blend_width"]
        if blend >= 0.5 * min(x_hi - x_lo, y_hi - y_lo):
            if n["blend_width"] == "auto":
                blend = 0.25 * min(x_hi - x_lo, y_hi - y_lo)
            else:
                problems.append("%s.blend_width = %g: must be below half the patch's width (%g m)"
                                % (tag, blend, min(x_hi - x_lo, y_hi - y_lo)))
        plan.patches.append(Patch(name, x_lo, x_hi, y_lo, y_hi, 2 * (i_hi - i_lo), 2 * (j_hi - j_lo), blend,
                                  n["max_iterations"] or case["solver"]["max_iterations"]))
    if plan.patches:
        if plan.buildings:
            problems.append("nests: not with surface.buildings (GUST refuses buildings with a nest face)")
        bad = sorted({"u", "v", "tke"} & set(plan.damping_fields))
        if bad and plan.damping_start is not None:
            problems.append("nests: atmosphere.damping.fields may hold only w and theta (a patch has no inflow profile "
                            "to damp %s towards)" % ", ".join(bad))
        if plan.side_faces == "periodic":
            problems.append("nests: not with wind.side_faces = periodic")


def _stretch_to_top(plan, dz0, top):
    """plan.layers layers from dz0 that end at top.

    The deck gets mesh.stretch_ratio (bisection, top met to about 1e-14), not
    mesh.dz0 alone: GUST's own Newton iteration for the ratio (ERF's
    geometric_stretch, from 1.03 with 50 steps) fails for steep stretches,
    e.g. 80 layers at a ratio of 1.1. A ratio within 1e-6 of 1 is written as
    uniform layers (no mesh.dz0): GUST refuses mesh.dz0 above top / n_cell z
    by any rounding, and its Newton stalls just below it.
    """
    r = grid.stretch_ratio_for_top(dz0, plan.layers, top)
    if r - 1.0 < UNIFORM_TOLERANCE:
        plan.dz0, plan.stretch_ratio = None, 1.0
    else:
        plan.dz0, plan.stretch_ratio = dz0, r
        plan.write_stretch_ratio = True
        plan.top = grid.stretched_top(dz0, r, plan.layers)


def resolve(case, given, base_dir):
    """The Plan of a validated case; raises CaseError listing every problem."""
    problems = []
    plan = Plan()
    plan.case, plan.given, plan.base_dir = case, given, base_dir
    name = case["case"]["name"]
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
        problems.append("case.name = '%s': letters, digits, _, - and . only, not starting with a dot (it names a folder under "
                        "case.output_dir; give the parent folder as case.output_dir)" % name)
    frame = _resolve_domain(case, given, problems)
    if frame is None:
        raise CaseError(problems)
    plan.frame = frame
    _resolve_mesh(case, frame, problems, plan)
    _resolve_wind(case, given, problems, plan, base_dir)
    _resolve_atmosphere(case, given, problems, plan)
    _resolve_surface(case, given, problems, plan, base_dir)
    _resolve_nests(case, problems, plan, frame)

    t = case["terrain"]
    plan.terrain_file = _path(base_dir, t["file"]) if t["file"] is not None else None
    if plan.terrain_file is None:
        for name in schema_terrain_keys():
            if "terrain." + name in given:
                problems.append("terrain.%s has no effect without terrain.file" % name)
        if not plan.roughness_map and plan.canopy != "map":
            for name in ("idw_neighbors", "idw_power"):
                if "terrain." + name in given:
                    problems.append("terrain.%s has no effect without terrain.file, surface.roughness_map or a "
                                    "canopy map" % name)
        if "solver.initial_blend" in given:
            problems.append("solver.initial_blend has no effect without terrain.file (the mass-consistent first "
                            "guess adjusts the wind to the terrain)")
        for i, n in enumerate(case["nests"] or []):
            if n["blend_width"] != "auto":
                problems.append("nests[%d].blend_width has no effect without terrain.file (it blends the patch's "
                                "terrain into the parent's)" % i)
    cell = min(plan.dx, plan.dy)
    plan.margin = MARGIN_CELLS * max(plan.dx, plan.dy) if t["margin"] == "auto" else t["margin"]
    plan.point_spacing = 0.5 * cell if t["point_spacing"] == "auto" else t["point_spacing"]
    b = t["border"]
    plan.border_flat, plan.border_slope, plan.border_level = b["flat_width"], b["slope_width"], b["level"]
    if "terrain.border.level" in given and b["flat_width"] == 0.0 and b["slope_width"] == 0.0:
        problems.append("terrain.border.level has no effect without terrain.border.flat_width or slope_width")
    widest = b["flat_width"] + b["slope_width"]
    if widest >= 0.5 * min(frame.x_hi - frame.x_lo, frame.y_hi - frame.y_lo) and widest > 0.0:
        problems.append("terrain.border: flat_width + slope_width = %g m reaches the middle of the domain" % widest)
    plan.rmap_file = _path(base_dir, case["surface"]["roughness_map"]["file"]) if plan.roughness_map else None
    regional = case["surface"]["regional_roughness_map"]["file"]
    plan.regional_file = _path(base_dir, regional) if regional is not None and plan.upstream_auto else None

    m = case["model"]
    plan.lee = m["lee_setting"] == "recommended"
    if not plan.lee and "model.blend_height" in given:
        problems.append("model.blend_height has no effect with model.lee_setting = none")
    cap = m["length_scale_cap"]
    if cap == "auto":
        cap = COUETTE_LENGTH_SCALE_CAP if plan.drive == "couette" else DEFAULT_LENGTH_SCALE_CAP
    plan.length_scale_cap = cap
    gi = case["gust_inputs"]
    if gi is not None:
        if not isinstance(gi, dict):
            problems.append("gust_inputs: must be a mapping of GUST keys to values (namespace.key: value)")
        else:
            for key in gi:
                if not isinstance(key, str) or key.count(".") != 1 or not all(key.split(".")):
                    problems.append("gust_inputs: '%s' is not a GUST key (namespace.key)" % key)
    if problems:
        raise CaseError(problems)
    return plan


def schema_terrain_keys():
    """The terrain keys that need a terrain file (idw_neighbors and idw_power also serve the maps)."""
    from .schema import SPEC
    names = [name for name, entry in SPEC["terrain"].items()
             if not isinstance(entry, dict) and name not in ("file", "idw_neighbors", "idw_power")]
    return names + ["border." + name for name in SPEC["terrain"]["border"]]
