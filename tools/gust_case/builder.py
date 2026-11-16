"""Write a case folder: decks, terrain and roughness files, a run script and a
summary.

Layout (one wind direction)::

    <output_dir>/<name>/
        inputs                the terrain deck
        terrain.xyz           when terrain.file is given
        roughness.xyz         when surface.roughness_map is given
        column/inputs         the column run of a geostrophic drive
        nest_<name>/inputs    one folder per nested patch (+ terrain_nest_<name>.xyz)
        run.sh                runs them in order
        case_resolved.yaml    the case description with every default filled in
        case_summary.txt      mesh, wind, files and the checker's findings

A sweep (``wind.directions``) has one ``dir_<DDD>/`` folder per direction, each
with its own ``inputs`` (and ``column/``, ``nest_*/``); the terrain files are
shared one level up.

Each deck writes the keys the case needs and leaves GUST's defaults out, so a
deck stays short enough to read. ``abl.drive`` and ``abl.theta`` name the
set-up, and GUST checks the deck against them when it starts.
"""
import os
import shutil

import yaml

from . import roughness as roughness_mod
from . import terrain as terrain_mod
from . import wind
from .deck import Deck
from .errors import CaseError
from .plan import resolve

BLEND_RELIEF_FACTOR = 2.5     # convection.hybrid_height = 2.5 terrain heights (docs/howto.rst)
REALIZABILITY = 0.309
NUT_RELAX = 0.5
BLEND_SIGMA_MAX = 0.1
COLUMN_CELLS = 4
COLUMN_ITERATIONS = 200
SMOOTHING = {"btf": 0, "stf": 1, "sullivan": 2}

COLUMN_TABLE = "column/column_guess.dat"   # the column run's table, from the terrain run's folder


class Result:
    """What :func:`build` made: the case folder, the decks and the notes."""

    def __init__(self):
        self.case_dir = None
        self.decks = []        # (path relative to the case folder, Deck, kind)
        self.files = []        # other files written, relative to the case folder
        self.notes = []
        self.findings = []     # (deck path, formatted checker report)
        self.reports = []      # (deck path, gust_input_checker Report) of every deck checked
        self.runs = []         # the run script's steps: dir, label, column, wind_from_column, patches, chk
        self.plan = None       # the resolved Plan (frame, cells, layers, patches)
        self.terrain = None    # the prepared Terrain (points in the case frame), or None without terrain.file
        self.errors = 0
        self.warnings = 0
        self.summary = ""


def _rel(path, start):
    return os.path.relpath(path, start).replace(os.sep, "/")


class Builder:
    def __init__(self, plan):
        self.plan = plan
        self.case = plan.case
        self.terrain = None
        self.roughness = None
        self.patch_terrain = {}
        self.notes = []
        self._noted = set()

    # ----------------------------------------------------------------------------------------------------------
    # data files
    # ----------------------------------------------------------------------------------------------------------
    def prepare_data(self):
        p, t = self.plan, self.case["terrain"]
        if p.terrain_file is not None:
            self.terrain = terrain_mod.prepare(p.terrain_file, t["format"], t["coordinates"], p.frame, p.margin,
                                               p.point_spacing, t["zero_level"], p.border_flat, p.border_slope,
                                               p.border_level)
            self.notes += self.terrain.notes
            self.notes.append("terrain: z from %.6g to %.6g m (relief %.6g m) after subtracting %.6g m"
                              % (self.terrain.lowest, self.terrain.highest, self.terrain.relief, self.terrain.zero_shift))
            if self.terrain.lowest < 0.0:
                self.notes.append("terrain: the ground goes below z = 0 (to %.6g m); terrain.zero_level sets the datum"
                                  % self.terrain.lowest)
            limit = p.damping_start if p.damping_start is not None else p.top
            if self.terrain.highest >= limit:
                raise CaseError("the terrain reaches %.6g m, which is not below %s = %.6g m; raise domain.top%s"
                                % (self.terrain.highest,
                                   "the Rayleigh layer (atmosphere.damping.start)" if p.damping_start is not None
                                   else "the top of the domain (domain.top)", limit,
                                   " or atmosphere.damping.start" if p.damping_start is not None else ""))
            for patch in p.patches:
                flat = self.terrain.flat_level
                self.patch_terrain[patch.name] = terrain_mod.prepare(
                    p.terrain_file, t["format"], t["coordinates"], self._patch_frame(patch),
                    p.margin, 0.25 * min(p.dx, p.dy) if t["point_spacing"] == "auto" else p.point_spacing,
                    self.terrain.zero_shift, p.border_flat, p.border_slope,
                    0.0 if flat is None else flat + self.terrain.zero_shift, border_frame=p.frame)
        if p.roughness_map:
            r = self.case["surface"]["roughness_map"]
            self.roughness = roughness_mod.prepare(p.rmap_file, r["format"], r["coordinates"], r["classes"], r["default"],
                                                   p.frame, p.margin, p.point_spacing)
            self.notes += self.roughness.notes

    def _patch_frame(self, patch):
        """The patch's box, in a frame that shares the parent's origin and projection."""
        p = self.plan
        f = terrain_mod.Frame(p.frame.center, 1.0, 1.0, 1.0, 1.0, lat_lon=p.frame.lat_lon, zone=p.frame.zone,
                              south_hemisphere=p.frame.south_hemisphere)
        f.origin = p.frame.origin
        f.x_lo, f.x_hi, f.y_lo, f.y_hi = patch.x_lo, patch.x_hi, patch.y_lo, patch.y_hi
        return f

    # ----------------------------------------------------------------------------------------------------------
    # deck sections shared by the terrain deck, the column deck and the patches
    # ----------------------------------------------------------------------------------------------------------
    def _vertical(self, deck):
        p = self.plan
        if p.dz0 is not None:          # stretched: always with the ratio (plan._stretch_to_top)
            deck.set("mesh.dz0", float(p.dz0))
            # all 17 digits: at %.10g, 80 layers at 1.1 end 7e-5 m off prob_hi z
            deck.set("mesh.stretch_ratio", "%.17g" % p.stretch_ratio, "the layers end at %.6g m" % p.top)

    def _theta(self, deck, column=False):
        p = self.plan
        deck.section("Potential temperature")
        deck.set("bc.theta_profile_z", [float(v) for v in p.theta_heights])
        deck.set("bc.theta_profile", [float(v) for v in p.theta_values])
        if p.theta == "neutral":
            deck.set("forcing.gravity", 0.0, "neutral: buoyancy off")
        else:
            if p.theta_values[0] != 300.0:
                deck.set("forcing.theta0", float(p.theta_values[0]))
        if p.theta == "capped":
            deck.set("bc.zhi_theta_grad", float(p.top_theta_gradient), "the sounding's gradient at the top")
        if p.theta == "surface_flux" and not column:
            deck.set("bc.surface_heat_flux", float(p.surface_heat_flux))

    def _damping(self, deck, fields=None):
        p = self.plan
        if p.damping_start is None:
            return
        deck.section("Rayleigh layer below the top")
        deck.set("bc.rayleigh_zdamp", float(p.damping_start))
        deck.set("bc.rayleigh_dampcoef", float(p.damping_coefficient))
        deck.set("bc.rayleigh_fields", " ".join(fields or p.damping_fields))

    def _model(self, deck, blend=True):
        p, m = self.plan, self.case["model"]
        deck.section("Closure")
        if m["tke_floor"] > 0.0:
            deck.set("turb.tke_floor", float(m["tke_floor"]))
        deck.set("turb.max_geom_lscale", float(p.length_scale_cap))
        if p.lee:
            deck.set("turb.realizability", REALIZABILITY, "recommended lee setting (docs/howto.rst)")
            deck.set("turb.nut_relax", NUT_RELAX)
            height = m["blend_height"]
            if height == "auto":
                height = BLEND_RELIEF_FACTOR * self.terrain.relief if self.terrain is not None else 0.0
            if blend and height > 0.0:
                deck.set("convection.hybrid_sigma_max", BLEND_SIGMA_MAX)
                deck.set("convection.hybrid_height", float(height),
                         "2.5 x the relief of %.6g m" % self.terrain.relief if m["blend_height"] == "auto" else None)
            elif blend and "blend-flat" not in self._noted:
                self._noted.add("blend-flat")
                self.notes.append("model.lee_setting: no convection blend over flat ground (model.blend_height = auto "
                                  "is 2.5 x the relief, which is 0); give model.blend_height to have one")

    def _surface_extras(self, deck, run_dir, case_dir):
        p = self.plan
        if p.canopy == "patches":
            deck.section("Canopy")
            deck.set("canopy.type", "patches")
            deck.set("canopy.file", _rel(os.path.join(case_dir, "canopy_patches.txt"), run_dir))
        elif p.canopy == "map":
            deck.section("Canopy")
            deck.set("canopy.type", "map")
            deck.set("canopy.height_file", _rel(os.path.join(case_dir, "canopy_height.xyz"), run_dir))
            deck.set("canopy.lai_file", _rel(os.path.join(case_dir, "canopy_lai.xyz"), run_dir))
            if p.canopy_cd is not None:
                deck.set("canopy.cd", float(p.canopy_cd))
            else:
                deck.set("canopy.cd_file", _rel(os.path.join(case_dir, "canopy_cd.xyz"), run_dir))
            if p.canopy_profile != "uniform":
                deck.set("canopy.profile", p.canopy_profile)
        if p.buildings:
            deck.section("Buildings")
            deck.set("buildings.type", "immersed_forcing")
            name = "buildings." + ("xyz" if p.buildings_format == "xyz" else "txt")
            deck.set("buildings.file", _rel(os.path.join(case_dir, name), run_dir))
            if p.buildings_format != "erf":
                deck.set("buildings.format", p.buildings_format)

    def _idw(self, deck):
        t = self.case["terrain"]
        if t["idw_neighbors"] is not None:
            deck.set("terrain.idw_neighbors", t["idw_neighbors"])
        if t["idw_power"] is not None:
            deck.set("terrain.idw_power", float(t["idw_power"]))

    def _solver(self, deck, sequencing, max_iterations):
        s = self.case["solver"]
        deck.section("Steady iteration")
        if s["w_vnn_cap"] > 0.0:
            deck.set("simplec.w_vnn_cap", float(s["w_vnn_cap"]))
        deck.set("steady.stage", "simplec")
        if s["sweep"] != "auto":
            deck.set("steady.sweep", s["sweep"])
        deck.set("steady.cfl_max", float(s["cfl_max"]))
        if sequencing:
            deck.set("steady.sequencing", 1)
        deck.set("steady.max_iterations", max_iterations)
        if s["residual_tolerance"] is not None:
            deck.set("steady.residual_tol", float(s["residual_tolerance"]))
        if s["residual_drop_orders"] is not None:
            deck.set("steady.residual_drop_orders", float(s["residual_drop_orders"]))

    def _output(self, deck, tag, checkpoint):
        o = self.case["output"]
        prefix = (o["prefix"] + "_" if o["prefix"] else "") + (tag + "_" if tag else "")
        deck.section("Output")
        if prefix:
            deck.set("io.plot_file", prefix + "plt")
            deck.set("io.chk_file", prefix + "chk")
        if o["turbulence_plotfile"]:
            deck.set("io.turb_plot_file", prefix + "turb")
        if not checkpoint:
            deck.set("io.write_checkpoint", 0)
        return prefix + "chk"

    def _title(self, what):
        c = self.case["case"]
        lines = ["GUST deck written by tools/make_gust_case.py: %s, %s" % (c["name"], what)]
        if c["description"]:
            lines.append(c["description"])
        lines.append("Change the YAML file and write the case again, rather than this file.")
        return lines

    def _apply_gust_inputs(self, deck):
        gi = self.case["gust_inputs"]
        if not gi:
            return
        deck.section("gust_inputs of the case file")
        for key, value in gi.items():
            if value is None:
                deck.remove(key)
            elif key in deck:
                deck.replace(key, value, "gust_inputs")
            else:
                deck.set(key, value)

    # ----------------------------------------------------------------------------------------------------------
    # the terrain deck
    # ----------------------------------------------------------------------------------------------------------
    def terrain_deck(self, direction, run_dir, case_dir):
        p = self.plan
        f = p.frame
        towards = wind.towards_from_meteorological(direction)
        faces = wind.lateral_faces(towards, p.side_faces)
        cx, cy = wind.unit_vector(towards)
        deck = Deck(self._title("wind from %g degrees" % direction))

        deck.section("Mesh: %d x %d x %d cells of %.6g x %.6g m" % (p.nx, p.ny, p.layers, p.dx, p.dy))
        deck.set("mesh.n_cell", [p.nx, p.ny, p.layers])
        if f.x_lo != 0.0 or f.y_lo != 0.0:
            deck.set("mesh.prob_lo", [float(f.x_lo), float(f.y_lo), 0.0])
        deck.set("mesh.prob_hi", [float(f.x_hi), float(f.y_hi), float(self.case["domain"]["top"])])
        periodic = [int(faces["xlo"] == "periodic"), int(faces["ylo"] == "periodic"), 0]
        if any(periodic):
            deck.set("mesh.is_periodic", periodic)
        deck.set("mesh.max_grid_size", p.max_grid_size)
        self._vertical(deck)
        if p.hierarchy_levels > 1:
            deck.set("mesh.hierarchy_levels", p.hierarchy_levels)

        deck.section("Terrain")
        if self.terrain is not None:
            deck.set("terrain.type", "xyz")
            deck.set("terrain.file", _rel(os.path.join(case_dir, "terrain.xyz"), run_dir))
            deck.set("terrain.smoothing", SMOOTHING[self.case["terrain"]["smoothing"]])
        else:
            deck.set("terrain.type", "flat")
        self._idw(deck)

        deck.section("Boundaries: wind from %g degrees (towards %g, counter-clockwise from x)" % (direction, towards))
        for face in ("xlo", "xhi", "ylo", "yhi"):
            deck.set("bc." + face, faces[face])
        deck.set("bc.zlo", "surface_layer")
        if p.drive == "couette":
            lid = wind.log_law_speed(p.ustar, p.top, p.upstream_roughness)
            deck.set("bc.zhi", "no_slip_wall")
            deck.set("bc.zhi_wall_velocity", [lid * cx, lid * cy, 0.0], "the log law's %.6g m/s at the top" % lid)
        else:
            deck.set("bc.zhi", "slip_wall")
        if self.roughness is not None:
            deck.set("bc.surface_z0_file", _rel(os.path.join(case_dir, "roughness.xyz"), run_dir))
        else:
            deck.set("bc.surface_z0", float(p.surface_roughness))

        deck.section("Inflow profile")
        if p.geostrophic:
            deck.set("bc.inflow_profile", "file")
            deck.set("bc.inflow_file", COLUMN_TABLE, "written by the column run")
            deck.set("bc.inflow_z0", float(p.upstream_roughness))
        elif p.profile == "table":
            deck.set("bc.inflow_profile", "file")
            deck.set("bc.inflow_file", _rel(os.path.join(case_dir, "inflow_table.dat"), run_dir))
            deck.set("bc.inflow_z0", float(p.upstream_roughness))
        elif p.profile == "power_law":
            deck.set("bc.inflow_profile", "power_law")
            deck.set("bc.inflow_speed", float(p.speed))
            deck.set("bc.inflow_zref", float(p.height))
            if p.power_law_exponent is not None:
                deck.set("bc.inflow_alpha", float(p.power_law_exponent))
            deck.set("bc.inflow_z0", float(p.upstream_roughness))
            deck.set("bc.inflow_direction", float(towards))
        else:
            deck.set("bc.inflow_profile", "log_law")
            deck.set("bc.inflow_ustar", float(p.ustar), "%.6g m/s at %.6g m over z0 = %.6g m"
                     % (p.speed, p.height, p.upstream_roughness))
            deck.set("bc.inflow_z0", float(p.upstream_roughness))
            deck.set("bc.inflow_direction", float(towards))

        self._theta(deck)
        self._damping(deck)

        deck.section("Set-up names (GUST checks the deck against them)")
        deck.set("abl.drive", p.drive)
        deck.set("abl.theta", p.theta)
        if p.geostrophic:
            deck.section("Rotation")
            deck.set("forcing.latitude", float(p.latitude))
            if p.geostrophic_wind is not None:
                deck.set("forcing.ug", float(p.geostrophic_wind[0]))
                deck.set("forcing.vg", float(p.geostrophic_wind[1]))
        self._surface_extras(deck, run_dir, case_dir)
        self._model(deck)

        deck.section("First guess")
        a_h, a_v = self.case["solver"]["initial_blend"]
        if p.drive == "geostrophic_column":
            deck.set("init.first_guess", "column")
        else:
            deck.set("init.first_guess", "profile")
        if self.terrain is not None:
            deck.set("init.alpha_h", float(a_h))
            deck.set("init.alpha_v", float(a_v))
        self._solver(deck, p.grid_sequencing, self.case["solver"]["max_iterations"])
        chk = self._output(deck, "", self.case["output"]["checkpoint"] or bool(p.patches))
        self._apply_gust_inputs(deck)
        return deck, chk

    # ----------------------------------------------------------------------------------------------------------
    # the column deck of a geostrophic drive
    # ----------------------------------------------------------------------------------------------------------
    def column_deck(self, direction):
        p = self.plan
        towards = wind.towards_from_meteorological(direction)
        cx, cy = wind.unit_vector(towards)
        deck = Deck(self._title("column run for the inflow profile, wind from %g degrees" % direction))
        deck.section("Mesh: a periodic column with the terrain deck's layers")
        deck.set("mesh.n_cell", [COLUMN_CELLS, COLUMN_CELLS, p.layers])
        deck.set("mesh.prob_hi", [COLUMN_CELLS * p.dx, COLUMN_CELLS * p.dy, float(self.case["domain"]["top"])])
        deck.set("mesh.is_periodic", [1, 1, 0])
        deck.set("mesh.max_grid_size", [COLUMN_CELLS, COLUMN_CELLS])
        self._vertical(deck)
        deck.section("Flat ground with the inflow's roughness")
        deck.set("terrain.type", "flat")
        deck.set("bc.zlo", "surface_layer")
        deck.set("bc.zhi", "slip_wall")
        deck.set("bc.surface_z0", float(p.upstream_roughness))
        self._theta(deck, column=True)
        self._damping(deck)
        deck.section("The column: GUST's implicit 1-D solve, adjusted to the wind asked for")
        deck.set("init.first_guess", "column")
        if p.geostrophic_wind is None:
            deck.set("init.column_target_u", p.speed * cx, "%.6g m/s from %g degrees" % (p.speed, direction))
            deck.set("init.column_target_v", p.speed * cy)
            deck.set("init.column_target_height", float(p.height))
        deck.set("abl.drive", "geostrophic_column")
        deck.set("abl.theta", p.theta)
        deck.set("forcing.latitude", float(p.latitude))
        if p.geostrophic_wind is not None:
            deck.set("forcing.ug", float(p.geostrophic_wind[0]))
            deck.set("forcing.vg", float(p.geostrophic_wind[1]))
        else:
            # GUST scales and turns the column's forcing by the ratio of the target to the column's wind, so the
            # column needs a wind to start from: the neutral drag law's speed, along the wind asked for.
            ustar = wind.friction_velocity(p.speed, p.height, p.upstream_roughness)
            start = wind.geostrophic_speed(ustar, p.upstream_roughness, p.latitude)
            deck.set("forcing.ug", start * cx, "where the column starts (the drag law's %.6g m/s); GUST adjusts it" % start)
            deck.set("forcing.vg", start * cy)
        deck.section("Closure")
        if self.case["model"]["tke_floor"] > 0.0:
            deck.set("turb.tke_floor", float(self.case["model"]["tke_floor"]))
        deck.set("turb.max_geom_lscale", float(p.length_scale_cap))
        deck.section("Steady iteration (the column is already the 1-D steady state)")
        deck.set("steady.stage", "simplec")
        deck.set("steady.sweep", "zjacobi")
        deck.set("steady.max_iterations", COLUMN_ITERATIONS)
        deck.section("Output")
        deck.set("io.plot_file", "column_plt")
        deck.set("io.write_checkpoint", 0)
        return deck

    # ----------------------------------------------------------------------------------------------------------
    # a nested patch
    # ----------------------------------------------------------------------------------------------------------
    def patch_deck(self, patch, direction, run_dir, case_dir):
        p = self.plan
        deck = Deck(self._title("nested patch '%s', wind from %g degrees" % (patch.name, direction)))
        deck.section("Mesh: twice the parent's horizontal resolution, the parent's layers")
        deck.set("mesh.n_cell", [patch.nx, patch.ny, p.layers])
        deck.set("mesh.prob_lo", [float(patch.x_lo), float(patch.y_lo), 0.0])
        deck.set("mesh.prob_hi", [float(patch.x_hi), float(patch.y_hi), float(self.case["domain"]["top"])])
        deck.set("mesh.max_grid_size", p.max_grid_size)
        self._vertical(deck)
        deck.section("Terrain")
        if self.terrain is not None:
            deck.set("terrain.type", "xyz")
            deck.set("terrain.file", _rel(os.path.join(case_dir, "terrain_nest_%s.xyz" % patch.name), run_dir))
            deck.set("terrain.smoothing", SMOOTHING[self.case["terrain"]["smoothing"]])
            if patch.blend_width > 0.0:
                deck.set("terrain.nest_blend_width", float(patch.blend_width))
        else:
            deck.set("terrain.type", "flat")
        self._idw(deck)
        deck.section("Boundaries: the parent's solution on the four lateral faces")
        for face in ("xlo", "xhi", "ylo", "yhi"):
            deck.set("bc." + face, "nest")
        deck.set("bc.zlo", "surface_layer")
        if p.drive == "couette":
            towards = wind.towards_from_meteorological(direction)
            cx, cy = wind.unit_vector(towards)
            lid = wind.log_law_speed(p.ustar, p.top, p.upstream_roughness)
            deck.set("bc.zhi", "no_slip_wall")
            deck.set("bc.zhi_wall_velocity", [lid * cx, lid * cy, 0.0])
        else:
            deck.set("bc.zhi", "slip_wall")
        if self.roughness is not None:
            deck.set("bc.surface_z0_file", _rel(os.path.join(case_dir, "roughness.xyz"), run_dir))
        else:
            deck.set("bc.surface_z0", float(p.surface_roughness))
        self._theta(deck)
        self._damping(deck)
        if p.geostrophic:
            deck.section("Rotation")
            deck.set("forcing.coriolis", 1)
            deck.set("forcing.geostrophic", 1)
            deck.set("forcing.latitude", float(p.latitude))
            if p.geostrophic_wind is not None:
                deck.set("forcing.ug", float(p.geostrophic_wind[0]))
                deck.set("forcing.vg", float(p.geostrophic_wind[1]))
        self._surface_extras(deck, run_dir, case_dir)
        self._model(deck)
        deck.section("First guess: the parent's field (bc.nest_parent comes from run.sh)")
        deck.set("init.first_guess", "nest")
        self._solver(deck, False, patch.max_iterations)
        self._output(deck, "", self.case["output"]["checkpoint"])
        self._apply_gust_inputs(deck)
        return deck


def _run_script(case, plan, runs):
    r = case["run"]
    launcher = r["mpi_launcher"].replace("{ranks}", str(r["ranks"]))
    lines = [
        "#!/bin/bash",
        "# Runs the case written by tools/make_gust_case.py: %s" % case["case"]["name"],
        "#   ./run.sh [key=value ...]     more GUST inputs for the terrain runs; a key the deck holds is refused by GUST",
        "#   GUST=/path/to/gust MPI=\"srun -n 64\" ./run.sh     another executable or launcher (MPI= runs without one)",
        "set -euo pipefail",
        "GUST=${GUST:-%s}" % r["executable"],
        "MPI=${MPI-%s}" % launcher,
        'cd "$(dirname "$0")"',
        "",
    ]
    for run in runs:
        lines.append("# ---- %s" % run["label"])
        lines.append("(")
        if run["dir"] != ".":
            lines.append("  cd %s" % run["dir"])
        extra = ""
        if run["column"]:
            lines.append('  ( cd column && $MPI "$GUST" inputs > log.txt )')
            lines.append("  mkdir -p column && [ -s column/column_guess.dat ] || "
                         "{ echo 'run.sh: the column run wrote no column/column_guess.dat (see column/log.txt)' >&2; exit 1; }")
            if run["wind_from_column"]:
                lines.append("  # the geostrophic wind the column run was adjusted to")
                lines.append("  ug=$(awk '$1 == \"forcing.ug\" && $2 == \"=\" {v = $3} END {print v}' column/log.txt)")
                lines.append("  vg=$(awk '$1 == \"forcing.vg\" && $2 == \"=\" {v = $3} END {print v}' column/log.txt)")
                lines.append("  [ -n \"$ug\" ] && [ -n \"$vg\" ] || "
                             "{ echo 'run.sh: no geostrophic wind in column/log.txt' >&2; exit 1; }")
                extra = ' forcing.ug="$ug" forcing.vg="$vg"'
        if run["patches"]:
            lines.append("  touch .parent_start      # the patches take a checkpoint this run writes, not an older one")
        lines.append('  $MPI "$GUST" inputs%s "$@" > log.txt' % extra)
        if run["patches"]:
            # the highest step by number (chk100000 sorts before chk99999 as text)
            lines.append("  parent=$(find . -maxdepth 1 -type d -name '%s[0-9]*' -newer .parent_start | sed 's|^\\./||' | "
                         "awk -v L=%d '{n = substr($0, L + 1) + 0; if (f == \"\" || n > m) {m = n; f = $0}} END {print f}')"
                         % (run["chk"], len(run["chk"])))
            lines.append("  [ -n \"$parent\" ] || { echo 'run.sh: the parent run wrote no checkpoint (see log.txt)' >&2; exit 1; }")
        for patch in run["patches"]:
            lines.append('  ( cd nest_%s && $MPI "$GUST" inputs bc.nest_parent="../$parent"%s "$@" > log.txt )' % (patch, extra))
        lines.append(")")
        lines.append("")
    return "\n".join(lines)


def _check(path, overrides, result, rel):
    """Run the deck checker (tools/gust_input_checker) on a deck that was written."""
    try:
        from gust_input_checker import Options, Severity, check_deck
    except ImportError:
        import sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        try:
            from gust_input_checker import Options, Severity, check_deck
        except ImportError:
            return False
        finally:
            sys.path.pop(0)
    # Each run starts in its deck's folder (run.sh). The column's table and the parent's checkpoint are
    # written by an earlier step of run.sh, so they are not looked for.
    made = [COLUMN_TABLE] + [o.split("=", 1)[1] for o in overrides if o.startswith("bc.nest_parent=")]
    report = check_deck(path, overrides=tuple(overrides),
                        options=Options(cwd=os.path.dirname(os.path.abspath(path)), made_by_run=tuple(made)))
    result.reports.append((rel, report))
    result.errors += report.count(Severity.ERROR)
    result.warnings += report.count(Severity.WARNING)
    if report.count(Severity.WARNING) or report.count(Severity.ERROR):
        result.findings.append((rel, report.format(Severity.WARNING)))
    return True


MANIFEST = ".gust_case_files"


def build(case, given, base_dir, dry_run=False, overwrite=False, check=True):
    """Write the case folder of a validated case description; returns a Result.

    ``dry_run`` resolves everything and prepares the decks, but writes nothing
    (the terrain file is still read, so its problems show).
    """
    plan = resolve(case, given, base_dir)
    builder = Builder(plan)
    result = Result()
    out_root = case["case"]["output_dir"]
    out_root = out_root if os.path.isabs(out_root) else os.path.join(base_dir, out_root)
    case_dir = os.path.normpath(os.path.join(out_root, case["case"]["name"]))
    result.case_dir = case_dir
    if os.path.exists(case_dir) and not dry_run:
        if not overwrite:
            raise CaseError("%s exists; use --overwrite to replace the files this tool writes (case.name, "
                            "case.output_dir)" % case_dir)
        if not os.path.isdir(case_dir):
            raise CaseError("%s exists and is not a folder" % case_dir)
    builder.prepare_data()

    checked = False
    runs = []
    for direction in plan.directions:
        sub = "dir_" + wind.direction_label(direction) if plan.sweep else "."
        run_dir = os.path.normpath(os.path.join(case_dir, sub))
        deck, chk = builder.terrain_deck(direction, run_dir, case_dir)
        prefix = "" if sub == "." else sub + "/"
        overrides = []
        if plan.geostrophic and plan.geostrophic_wind is None:
            overrides = ["forcing.ug=1.", "forcing.vg=0."]
        result.decks.append((prefix + "inputs", deck, "terrain", overrides))
        if plan.geostrophic:
            result.decks.append((prefix + "column/inputs", builder.column_deck(direction), "column", []))
        for patch in plan.patches:
            result.decks.append((prefix + "nest_%s/inputs" % patch.name,
                                 builder.patch_deck(patch, direction, os.path.join(run_dir, "nest_" + patch.name), case_dir),
                                 "patch", overrides + ["bc.nest_parent=../%s00000" % chk]))
        runs.append({"dir": sub, "label": "wind from %g degrees" % direction, "column": plan.geostrophic,
                     "wind_from_column": plan.geostrophic and plan.geostrophic_wind is None,
                     "patches": [patch.name for patch in plan.patches], "chk": chk})
    result.notes = builder.notes
    result.runs = runs
    result.plan, result.terrain = plan, builder.terrain
    script = _run_script(case, plan, runs)

    if not dry_run:
        os.makedirs(case_dir, exist_ok=True)
        removed = _remove_previous_write(case_dir)
        if removed:
            result.notes.append("--overwrite: removed the %d files of the previous write (%s); run outputs stay"
                                % (removed, MANIFEST))
        if builder.terrain is not None:
            builder.terrain.write(os.path.join(case_dir, "terrain.xyz"))
            result.files.append("terrain.xyz")
            for name, patch_terrain in builder.patch_terrain.items():
                patch_terrain.write(os.path.join(case_dir, "terrain_nest_%s.xyz" % name))
                result.files.append("terrain_nest_%s.xyz" % name)
        if builder.roughness is not None:
            builder.roughness.write(os.path.join(case_dir, "roughness.xyz"))
            result.files.append("roughness.xyz")
        for target, source in plan.copied_files.items():
            shutil.copyfile(source, os.path.join(case_dir, target))
            result.files.append(target)
        for rel, deck, kind, overrides in result.decks:
            path = os.path.join(case_dir, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            deck.write(path)
        with open(os.path.join(case_dir, "run.sh"), "w", encoding="utf-8") as out:
            out.write(script + "\n")
        os.chmod(os.path.join(case_dir, "run.sh"), 0o755)
        with open(os.path.join(case_dir, "case_resolved.yaml"), "w", encoding="utf-8") as out:
            out.write("# The case description with every default filled in (tools/make_gust_case.py).\n")
            yaml.safe_dump(case, out, sort_keys=False, default_flow_style=None)
        result.files += ["run.sh", "case_resolved.yaml"]
        if check:
            checked = all([_check(os.path.join(case_dir, rel), overrides, result, rel)
                           for rel, deck, kind, overrides in result.decks])
            if not checked:
                result.notes.append("the deck checker (tools/gust_input_checker) was not found: the decks are not checked")
    result.summary = _summary(case, plan, builder, result, runs, checked)
    if not dry_run:
        with open(os.path.join(case_dir, "case_summary.txt"), "w", encoding="utf-8") as out:
            out.write(result.summary)
        result.files.append("case_summary.txt")
        with open(os.path.join(case_dir, MANIFEST), "w", encoding="utf-8") as out:
            out.write("# the files make_gust_case.py wrote here; --overwrite removes them first\n")
            for rel in result.files + [rel for rel, _, _, _ in result.decks]:
                out.write(rel + "\n")
    return result


def _remove_previous_write(case_dir):
    """Remove the files listed in the previous write's manifest, and the folders they leave empty.

    Only those files: plotfiles, checkpoints and logs of earlier runs stay, and a folder that holds them stays.
    Returns the number of files removed."""
    path = os.path.join(case_dir, MANIFEST)
    if not os.path.isfile(path):
        return 0
    removed = 0
    folders = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            rel = line.strip()
            if not rel or rel.startswith("#"):
                continue
            target = os.path.normpath(os.path.join(case_dir, rel))
            if os.path.dirname(os.path.relpath(target, case_dir)).startswith(".."):
                continue                        # never outside the case folder
            if os.path.isfile(target):
                os.remove(target)
                removed += 1
            folder = os.path.dirname(target)
            while folder != case_dir and folder.startswith(case_dir + os.sep):
                folders.add(folder)
                folder = os.path.dirname(folder)
    for folder in sorted(folders, key=len, reverse=True):
        if os.path.isdir(folder) and not os.listdir(folder):
            os.rmdir(folder)
    os.remove(path)
    return removed


def _summary(case, plan, builder, result, runs, checked):
    f = plan.frame
    lines = ["Case %s" % case["case"]["name"], "folder: %s" % result.case_dir, ""]
    lines.append("Domain  x %.6g .. %.6g m, y %.6g .. %.6g m, top %.6g m" % (f.x_lo, f.x_hi, f.y_lo, f.y_hi, plan.top))
    if f.lat_lon is not None:
        lines.append("        centre %.6f N, %.6f E = UTM zone %d%s, easting %.3f m, northing %.3f m"
                     % (f.lat_lon[0], f.lat_lon[1], f.zone, "S" if f.south_hemisphere else "N", f.center[0], f.center[1]))
    cells = plan.nx * plan.ny * plan.layers
    lines.append("Mesh    %d x %d x %d = %.3g million cells of %.6g x %.6g m" % (plan.nx, plan.ny, plan.layers, cells / 1.0e6,
                                                                               plan.dx, plan.dy))
    if plan.dz0 is None:
        lines.append("        uniform layers of %.6g m" % (plan.top / plan.layers))
    else:
        last = plan.dz0 * plan.stretch_ratio ** (plan.layers - 1)
        lines.append("        layers from %.6g m at the ground to %.6g m at the top, ratio %.6g" % (plan.dz0, last, plan.stretch_ratio))
    lines.append("        %d hierarchy level%s%s" % (plan.hierarchy_levels, "" if plan.hierarchy_levels == 1 else "s",
                                                    ", grid sequencing on" if plan.grid_sequencing else ""))
    for patch in plan.patches:
        lines.append("Patch   %s: x %.6g .. %.6g m, y %.6g .. %.6g m, %d x %d x %d cells of %.6g x %.6g m"
                     % (patch.name, patch.x_lo, patch.x_hi, patch.y_lo, patch.y_hi, patch.nx, patch.ny, plan.layers,
                        plan.dx / 2.0, plan.dy / 2.0))
    lines.append("Wind    drive %s, theta %s" % (plan.drive, plan.theta))
    if plan.ustar is not None:
        lines.append("        log law: %.6g m/s at %.6g m over z0 = %.6g m gives u* = %.6g m/s"
                     % (plan.speed, plan.height, plan.upstream_roughness, plan.ustar))
    for run in runs:
        lines.append("        %s -> %s" % (run["label"], "inputs" if run["dir"] == "." else run["dir"] + "/inputs"))
    lines.append("")
    if result.notes:
        lines.append("Notes")
        lines += ["  " + note for note in result.notes]
        lines.append("")
    lines.append("Decks")
    lines += ["  " + rel for rel, _, _, _ in result.decks]
    if result.files:
        lines.append("Other files")
        lines += ["  " + name for name in result.files]
    lines.append("")
    if checked:
        if result.findings:
            lines.append("Deck checker (tools/check_gust_inputs.py): %d error%s, %d warning%s"
                         % (result.errors, "" if result.errors == 1 else "s", result.warnings,
                            "" if result.warnings == 1 else "s"))
            for rel, text in result.findings:
                lines.append("  %s" % rel)
                lines += ["    " + line for line in text.rstrip().splitlines()]
        else:
            lines.append("Deck checker (tools/check_gust_inputs.py): no errors or warnings")
        lines.append("")
    lines.append("Run it with ./run.sh in the case folder (GUST=/path/to/gust and MPI=... change the command).")
    return "\n".join(lines) + "\n"
