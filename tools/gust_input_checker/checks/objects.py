"""Canopies and buildings.

ERROR:

- canopy.* (Source/Canopy.cpp check_canopy_options, Source/App.cpp): keys
  given with canopy.type = none; the k-term coefficients with canopy.tke =
  off, or all three zero with it on; the leaf keys with canopy.heat = off;
  canopy.heat = fixed without canopy.leaf_theta; the keys of the other
  canopy type; a map without its drag coefficient or with both forms of it;
  canopy.laimax with a uniform profile; canopy.heat_frozen while theta is
  solved;
- buildings.* (Source/Buildings.cpp check_buildings_options, Source/App.cpp):
  keys given with buildings.type = none; buildings.file missing; the wall
  keys with buildings.wall_form = drag; buildings.theta_mode = fixed without
  buildings.theta_surface; buildings.limiter_zone without the wall k or
  turb.realizability; the Newton stage, the lusgs sweep, a mesh hierarchy,
  solver.poisson_backend = mlmg; a building tile with h > 0 within
  2 cells + max(bc.nest_relax_width cells, terrain.nest_blend_width) of a
  nest face (Source/NestBuildings.cpp check_nest_buildings; reported on
  that face's key);
  buildings.force_file without a stage;
- terrain.representation = immersed_forcing (Source/Inputs.cpp, Source/ImmersedTerrain.cpp
  check_immersed_terrain and check_immersed_terrain_height): terrain.type other than flat; terrain.flat_height
  missing, or outside the top of the first layer to the bottom of the top layer; the terrain keys that follow or
  smooth a ground; a hierarchy; a mesh not periodic in x and y; no surface layer, a stability,
  bc.surface_zref, bc.surface_denominator_floor, bc.surface_wall_k = source, a stratified bc.theta_profile;
  buildings, a canopy with leaf heat, nests, the column first guess, probes, the Newton and unsteady stages, the lusgs sweep and
  solver.poisson_backend = mlmg in a stage, turb.closure = none; terrain.flat_height or terrain.small_volfrac with
  a body-fitted ground.
- terrain.representation = hybrid: the same, except that terrain.smoothing and terrain.gamma_m are accepted and
  terrain.flat_height may be 0 (>= 0 and at most the bottom of the top layer; not given: the first layer's
  thickness, and terrain.hybrid_half_width not given: 2 cells);
  terrain.hybrid_half_width above the domain's larger lateral length, or given without hybrid.

The files themselves are checked by `files`.
"""
import math
import os

from ..checker import Checker, register

NEST_BUILDING_CELLS = 2     # Source/NestBuildings.H, kNestBuildingCells


def _tile_edges(x):
    """Source/Buildings.cpp tile_edges: halfway between the points, the end tiles as far outward as inward."""
    e = [0.0] * (len(x) + 1)
    for i in range(1, len(x)):
        e[i] = 0.5 * (x[i - 1] + x[i])
    e[0] = x[0] - (e[1] - x[0])
    e[-1] = x[-1] + (x[-1] - e[-2])
    return e


def read_building_tiles(path, fmt):
    """(x edges, y edges, h[i][j]) of a buildings.file (Source/Buildings.cpp parse_map / parse_xyz), or None when
    the file is not in that format (GUST then stops with its own message)."""
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
            text = fh.read()
        if fmt == "xyz":
            pts = []
            for line in text.splitlines():
                v = [float(t) for t in line.split("#", 1)[0].replace(",", " ").split()]
                if v:
                    if len(v) not in (3, 4):
                        return None
                    pts.append(v)
            xs = sorted(set(p[0] for p in pts))
            ys = sorted(set(p[1] for p in pts))
            if len(xs) < 2 or len(ys) < 2 or len(pts) != len(xs) * len(ys):
                return None
            ix = {x: i for i, x in enumerate(xs)}
            iy = {y: j for j, y in enumerate(ys)}
            h = [[0.0] * len(ys) for _ in xs]
            for p in pts:
                h[ix[p[0]]][iy[p[1]]] = p[2]
        elif fmt == "erf":
            tok = [float(t) for t in text.split()]
            if len(tok) < 2 or tok[0] != int(tok[0]) or tok[1] != int(tok[1]) or tok[0] < 2 or tok[1] < 2:
                return None
            nx, ny = int(tok[0]), int(tok[1])
            if len(tok) != 2 + nx + ny + nx * ny:
                return None
            xs, ys, flat = tok[2:2 + nx], tok[2 + nx:2 + nx + ny], tok[2 + nx + ny:]
            h = [flat[i * ny:(i + 1) * ny] for i in range(nx)]
        else:
            return None
    except (OSError, ValueError, OverflowError):
        return None
    return _tile_edges(xs), _tile_edges(ys), h


@register
class CanopyChecker(Checker):
    """The canopy (canopy.*): which keys go with patches, with a map, with the k terms and with the leaf heat."""

    name = "canopy"

    def run(self):
        c = self.ctx
        kind = c.s("canopy.type")
        if kind not in ("none", "patches", "map"):
            return
        if kind == "none":
            for key in self.deck.keys("canopy."):
                if key != "canopy.type":
                    self.no_effect(key, "with canopy.type = none")
            return
        tke, heat = c.s("canopy.tke"), c.s("canopy.heat")
        if tke == "off":
            for k in ("beta_p", "beta_d", "length_alpha"):
                self.no_effect("canopy." + k, "with canopy.tke = off")
        elif tke == "on" and c.r("canopy.beta_p") == 0 and c.r("canopy.beta_d") == 0 \
                and c.r("canopy.length_alpha") == 0:
            self.error("canopy.tke", "GUST stops: canopy.tke = on with canopy.beta_p, canopy.beta_d and "
                                     "canopy.length_alpha all 0 has no k term left.", "Set canopy.tke = off.")
        if heat == "off":
            for k in ("leaf_theta", "heat_frozen"):
                self.no_effect("canopy." + k, "with canopy.heat = off")
        elif heat == "fixed":
            leaf = c.r("canopy.leaf_theta")
            if not self.given("canopy.leaf_theta"):
                self.error("canopy.leaf_theta", "GUST stops: canopy.heat = fixed needs canopy.leaf_theta (K).",
                           "Give the leaf temperature.")
            elif leaf is not None and not 150 <= leaf <= 400:
                self.error("canopy.leaf_theta", "GUST stops: canopy.leaf_theta = %.6g: must be from 150 to 400 K."
                           % leaf, "Give the leaf temperature in kelvin.")
            if c.solve and not c.theta_frozen:
                self.no_effect("canopy.heat_frozen", "without simplec.theta = frozen (theta is solved)")
            elif c.theta_frozen and c.s("canopy.heat_frozen") == "blend" and self.given("nest.names"):
                self.error("nest.names", "GUST stops: nest.names is not supported with simplec.theta = frozen and "
                                         "canopy.heat_frozen = blend (a nest would blend the outer domain's blended "
                                         "theta again).", "Set canopy.heat_frozen = skip, or remove the nests.")
            elif c.theta_frozen and c.s("canopy.heat_frozen") == "blend" and c.first_guess == "nest":
                self.error("canopy.heat_frozen", "GUST stops: init.first_guess = nest is not supported with "
                                                 "simplec.theta = frozen and canopy.heat_frozen = blend (the parent "
                                                 "run's theta may already carry the leaves' heat).",
                           "Set canopy.heat_frozen = skip.")
        if kind == "patches":
            if not c.s("canopy.file"):
                self.error("canopy.file", "GUST stops: canopy.type = patches needs canopy.file.",
                           "Give the patch file.")
            for k in ("height_file", "lai_file", "cd", "cd_file", "profile", "laimax"):
                self.no_effect("canopy." + k, "with canopy.type = patches (the patch file gives it)")
        else:
            self.no_effect("canopy.file", "with canopy.type = map")
            for k in ("height_file", "lai_file"):
                if not c.s("canopy." + k):
                    self.error("canopy." + k, "GUST stops: canopy.type = map needs canopy.%s." % k, "Give the raster.")
            has_cd, has_file = self.given("canopy.cd"), self.given("canopy.cd_file")
            if has_cd and has_file:
                self.error("canopy.cd", "GUST stops: canopy.cd and canopy.cd_file are exclusive: give one.",
                           "Keep the uniform coefficient or the raster.")
            elif not has_cd and not has_file:
                self.error("canopy.cd", "GUST stops: canopy.type = map needs canopy.cd or canopy.cd_file.",
                           "Give the drag coefficient.")
            if c.s("canopy.profile") == "uniform":
                self.no_effect("canopy.laimax", "with canopy.profile = uniform")
        if c.solve and c.closure == "kl":
            self.info("canopy.type", "Canopy: %s, k terms %s, leaf heat %s." % (kind, tke, heat))


@register
class BuildingsChecker(Checker):
    """Buildings (buildings.*): the wall form's keys, and the solver options buildings need."""

    name = "buildings"

    def run(self):
        c = self.ctx
        kind = c.s("buildings.type")
        if kind not in ("none", "immersed_forcing"):
            return
        if kind == "none":
            for key in self.deck.keys("buildings."):
                if key != "buildings.type":
                    self.no_effect(key, "with buildings.type = none")
            return
        form = c.s("buildings.wall_form")
        stress = form == "fraction_stress"
        tke = c.s("buildings.wall_tke") if self.given("buildings.wall_tke") else ("on" if stress else "off")
        mode = c.s("buildings.theta_mode") if self.given("buildings.theta_mode") else ("ambient" if stress else "none")
        if not c.s("buildings.file"):
            self.error("buildings.file", "GUST stops: buildings.type = immersed_forcing needs buildings.file.",
                       "Give the building file.")
        if tke == "on" and not stress:
            self.error("buildings.wall_tke", "GUST stops: buildings.wall_tke = on needs buildings.wall_form = "
                                             "fraction_stress.", "Set buildings.wall_tke = off, or the wall form.")
        if form == "drag":
            for k in ("wall_speed", "z0"):
                self.no_effect("buildings." + k, "with buildings.wall_form = drag (no wall law)")
        if mode in ("ambient", "fixed") and not stress:
            self.error("buildings.theta_mode", "GUST stops: buildings.theta_mode = %s needs buildings.wall_form = "
                                               "fraction_stress." % mode, "Set buildings.theta_mode = none.")
        if mode == "fixed":
            if not c.positive("buildings.theta_surface"):
                self.error("buildings.theta_surface", "GUST stops: buildings.theta_mode = fixed needs "
                                                      "buildings.theta_surface > 0 (K).", "Give the wall temperature.")
        else:
            self.no_effect("buildings.theta_surface", "without buildings.theta_mode = fixed")
        if tke != "on":
            self.no_effect("buildings.limiter_zone", "without buildings.wall_tke = on")
            for k in ("wall_k_bound", "wall_k_speed"):
                self.no_effect("buildings." + k, "without buildings.wall_tke = on (no wall value of k)")
        if tke == "on" and c.s("buildings.wall_k_speed") == "cell_mean" and c.s("buildings.wall_speed") == "all_faces":
            self.error("buildings.wall_k_speed", "GUST stops: buildings.wall_k_speed = cell_mean needs "
                                                 "buildings.wall_speed = open_faces.",
                       "Remove buildings.wall_speed = all_faces, or use buildings.wall_k_speed = law.")
        if tke == "on" and c.solve and c.positive("buildings.limiter_zone") and self.given("buildings.limiter_zone") \
                and not (c.realizability and c.realizability > 0):
            self.no_effect("buildings.limiter_zone", "without turb.realizability > 0 (the zone lifts that limiter)")
        if self.given("buildings.force_file"):
            if not c.s("buildings.force_file"):
                self.error("buildings.force_file", "GUST stops: buildings.force_file is given without a file name.",
                           "Give the file, or remove the key.")
            elif not c.solve:
                self.no_effect("buildings.force_file", "without steady.stage = simplec")
        if c.reads_bc and c.nest_faces:
            self._nest_gap()
        if c.solve:
            self._solver(mode)

    def _nest_gap(self):
        """A building tile within the gap of a nest face (Source/NestBuildings.cpp, check_nest_buildings)."""
        c = self.ctx
        if not (c.prob_lo and c.prob_hi and c.dx and c.dy):
            return
        path, _ = c.find_file(c.s("buildings.file"))
        if path is None or os.path.isdir(path):
            return                                  # reported by `files`
        tiles = read_building_tiles(path, c.s("buildings.format") or "erf")
        if tiles is None:
            return                                  # GUST reports the format when it reads the file
        xe, ye, h = tiles
        n_relax = c.i("bc.nest_relax_width") or 0
        blend = c.r("terrain.nest_blend_width") or 0.0
        periodic = c.periodic or [0, 0, 0]
        for f in c.nest_faces:
            d = 0 if f[0] == "x" else 1
            o = 1 - d
            cell = c.dx if d == 0 else c.dy
            ocell = c.dy if d == 0 else c.dx
            gap = NEST_BUILDING_CELLS * cell + max(n_relax * cell, blend)
            lo, hi = c.prob_lo[d], c.prob_hi[d]
            olo, ohi = c.prob_lo[o], c.prob_hi[o]
            tol, otol = 1e-9 * cell, 1e-9 * ocell
            e, eo = (xe, ye) if d == 0 else (ye, xe)
            near = []
            for i in range(len(xe) - 1):
                for j in range(len(ye) - 1):
                    top = h[i][j]
                    if not top > 0:
                        continue
                    a, b = (i, j) if d == 0 else (j, i)
                    e0, e1, f0, f1 = e[a], e[a + 1], eo[b], eo[b + 1]
                    if not periodic[o] and (f1 <= olo + otol or f0 >= ohi - otol):
                        continue
                    if (f.endswith("lo") and e1 > lo + tol and e0 < lo + gap - tol) or \
                            (f.endswith("hi") and e0 < hi - tol and e1 > hi - gap + tol):
                        near.append((xe[i], xe[i + 1], ye[j], ye[j + 1], top))
            if near:
                x0, x1, y0, y1, top = near[0]
                more = " (and %d more tiles)" % (len(near) - 1) if len(near) > 1 else ""
                self.error("bc." + f, "GUST stops: buildings.file = %s: the building tile x in [%g, %g] m, y in [%g, %g] m "
                                      "(top %g m) lies within %g m of the nest face bc.%s%s (%d cells beyond the larger of "
                                      "bc.nest_relax_width = %d cells and terrain.nest_blend_width = %g m)."
                           % (c.s("buildings.file"), x0, x1, y0, y1, top, gap, f, more, NEST_BUILDING_CELLS,
                              n_relax, blend), "Move the face or the building.")

    def _solver(self, mode):
        c = self.ctx
        if c.newton:
            self.error("steady.stage", "GUST stops: buildings are not supported with steady.stage = simplec_newton.",
                       "Use steady.stage = simplec.")
        if c.sweep != "zjacobi":
            given = self.given("steady.sweep")
            self.error("steady.sweep", "GUST stops: buildings.type = immersed_forcing needs steady.sweep = zjacobi%s."
                       % ("" if given else " (the default of a CPU build is lusgs)"), "Set steady.sweep = zjacobi.")
        if c.levels is not None and c.levels != 1:
            self.error("mesh.hierarchy_levels", "GUST stops: buildings need mesh.hierarchy_levels = 1.",
                       "Remove mesh.hierarchy_levels (and steady.sequencing).")
        if c.s("solver.poisson_backend") == "mlmg":
            self.error("solver.poisson_backend", "GUST stops: buildings are not supported with solver.poisson_backend = "
                                                 "mlmg.", "Use gmres_mlmg (the default) or gmres_fft.")
        if c.theta_frozen and self.given("buildings.theta_mode") and mode != "none":
            self.no_effect("buildings.theta_mode", "with simplec.theta = frozen (theta is not solved)")


@register
class ImmersedTerrainChecker(Checker):
    """The immersed ground (terrain.representation = immersed_forcing or hybrid): what it needs and what it refuses."""

    name = "immersed_terrain"

    def run(self):
        c = self.ctx
        rep = c.s("terrain.representation")
        hybrid = rep == "hybrid"
        if not hybrid:
            self.no_effect("terrain.hybrid_half_width", "unless terrain.representation = hybrid")
        if rep == "body_fitted":
            for k in ("flat_height", "small_volfrac"):
                self.no_effect("terrain." + k, "unless terrain.representation = immersed_forcing or hybrid")
            return
        if rep not in ("immersed_forcing", "hybrid"):
            return
        with_rep = "with terrain.representation = " + rep
        # App.cpp prints this advice on a fresh start only.
        if not hybrid and not c.restart and c.terrain_type not in ("", "flat"):
            self.warn("terrain.representation",
                      "terrain.representation = immersed_forcing over terrain.type = %s: the flat mesh cuts the ground "
                      "across many layers." % c.terrain_type,
                      "terrain.representation = hybrid is recommended: the mesh follows the terrain's lower envelope "
                      "and the ground stays within a few layers of it.")

        def refuse(key, what, advice=""):
            self.error(key, "GUST stops: terrain.representation = %s %s." % (rep, what), advice)

        h = c.r("terrain.flat_height")
        if hybrid:
            w = c.r("terrain.hybrid_half_width")
            if w is not None and c.extent and w > max(c.extent[0], c.extent[1]):
                self.error("terrain.hybrid_half_width",
                           "GUST stops: terrain.hybrid_half_width must be at most the domain's larger lateral length "
                           "(%g m, got %g)." % (max(c.extent[0], c.extent[1]), w), "Use a narrower square.")
            if c.terrain_type in ("", "flat"):
                self.warn("terrain.representation",
                          "terrain.representation = hybrid on flat terrain gives the immersed_forcing mesh: the "
                          "envelope of a flat ground is the ground.", "Use terrain.representation = immersed_forcing.")
            if w is not None and w > 0 and c.n_cell and c.dx and c.dy:
                steps = [d for d, n in ((c.dx, c.n_cell[0]), (c.dy, c.n_cell[1])) if n > 1]
                if steps and all(math.floor(w / d + 1e-9) == 0 for d in steps):
                    self.warn("terrain.hybrid_half_width",
                              "terrain.hybrid_half_width = %g is below one cell (%g m): the square has no nodes "
                              "and the mesh follows the terrain itself, as with 0." % (w, min(steps)),
                              "Give 0, or a width of at least one cell (default: 2 cells).")
            # Not given: the first layer's thickness (GUST's apply_hybrid_defaults), below the top layer.
            if h is not None and c.zlevels and len(c.zlevels) > 2 and not h <= c.zlevels[-2]:
                self.error("terrain.flat_height",
                           "GUST stops: terrain.flat_height must lie at or below the bottom of the mesh's top layer "
                           "(%g m) %s (got %g)." % (c.zlevels[-2], with_rep, h), "Move the ground, or add layers.")
        elif h is None or not h > 0:
            self.error("terrain.flat_height", "GUST stops: terrain.flat_height must be finite and > 0 " + with_rep + ".",
                       "Give the ground's height above the mesh's lower face.")
        elif c.zlevels and len(c.zlevels) > 2 and not c.zlevels[1] <= h <= c.zlevels[-2]:
            self.error("terrain.flat_height",
                       "GUST stops: terrain.flat_height must lie between the top of the mesh's first layer (%g m) and "
                       "the bottom of its top layer (%g m) %s (got %g)." % (c.zlevels[1], c.zlevels[-2], with_rep, h),
                       "Move the ground, or add layers below it.")
        for k in (("hierarchy_passes", "nest_blend_width") if hybrid
                  else ("smoothing", "gamma_m", "hierarchy_passes", "nest_blend_width")):
            self.no_effect("terrain." + k, with_rep)
        if c.levels != 1:
            refuse("mesh.hierarchy_levels", "needs mesh.hierarchy_levels = 1 for now")
        if c.s("bc.outflow_transport") == "from_velocity":
            refuse("bc.outflow_transport", "refuses bc.outflow_transport = from_velocity", "Use mirror (the default).")
        if not c.surface_layer:
            refuse("bc.zlo", "needs bc.zlo = surface_layer", "Set bc.zlo = surface_layer (read with a first guess or a stage).")
        if c.reads_bc:
            if c.s("bc.surface_stability") != "neutral":
                refuse("bc.surface_stability", "needs bc.surface_stability = neutral for now")
            if c.s("bc.surface_wall_k") == "source":
                refuse("bc.surface_wall_k", "needs bc.surface_wall_k = dirichlet")
            for k in ("surface_zref", "surface_denominator_floor"):
                self.no_effect("bc." + k, with_rep)
            if c.stratified:
                refuse("bc.theta_profile", "needs a uniform bc.theta_profile for now")
        if c.buildings:
            refuse("buildings.type", "cannot be combined with buildings.type for now")
        if c.canopy and c.s("canopy.heat") == "fixed":
            refuse("canopy.heat", "cannot be combined with canopy.heat = fixed for now", "Set canopy.heat = off.")
        if self.given("nest.names"):
            refuse("nest.names", "cannot be combined with nests for now")
        if c.first_guess == "column":
            refuse("init.first_guess", "is not supported with init.first_guess = column yet", "Use profile.")
        if self.given("steady.probes_file"):
            refuse("steady.probes_file", "is not supported with steady.probes_file yet")
        if c.newton or c.unsteady:
            refuse("steady.stage", "is not supported with steady.stage = %s yet" % c.stage, "Use simplec.")
        if c.solve and c.s("solver.poisson_backend") == "mlmg":
            refuse("solver.poisson_backend", "cannot use solver.poisson_backend = mlmg", "Use gmres_mlmg or gmres_fft.")
        if c.reads_bc and c.closure == "none":
            refuse("turb.closure", "needs turb.closure = kl")
