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
  solver.poisson_backend = mlmg and a nest face with buildings;
  buildings.force_file without a stage.

The files themselves are checked by `files`.
"""
from ..checker import Checker, register


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
        elif c.solve and c.positive("buildings.limiter_zone") and self.given("buildings.limiter_zone") \
                and not (c.realizability and c.realizability > 0):
            self.no_effect("buildings.limiter_zone", "without turb.realizability > 0 (the zone lifts that limiter)")
        if self.given("buildings.force_file"):
            if not c.s("buildings.force_file"):
                self.error("buildings.force_file", "GUST stops: buildings.force_file is given without a file name.",
                           "Give the file, or remove the key.")
            elif not c.solve:
                self.no_effect("buildings.force_file", "without steady.stage = simplec")
        if c.reads_bc and c.nest_faces:
            self.error("buildings.type", "GUST stops: buildings are not supported with a nest face (bc.%s = nest)."
                       % c.nest_faces[0], "Run the buildings on a mesh without nest faces.")
        if c.solve:
            self._solver(mode)

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
