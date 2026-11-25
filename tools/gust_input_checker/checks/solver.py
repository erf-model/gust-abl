"""The pseudo-time driver, the unsteady stage and the pressure solve.

ERROR:

- steady.* (Source/Driver.cpp check_steady_options, Source/LineSweeps.cpp,
  Source/App.cpp): every key but steady.stage without a stage; the Newton and
  switch keys without steady.stage = simplec_newton; steady.cfl_max below
  steady.cfl0; steady.newton_cfl_max below steady.cfl_max;
  steady.confidence_height without steady.confidence; steady.zline_sweeps
  outside its sweep's range; the steady.line_* keys without steady.sweep =
  lusgs; steady.line_overlap with lagged ghosts; steady.line_alternate with one
  sweep;
- unsteady.* (Source/TimeDriver.cpp check_unsteady_options, Source/App.cpp):
  every key without steady.stage = unsteady; unsteady.stop_time missing; none
  or both of unsteady.dt and unsteady.cfl; unsteady.average_start not below
  the stop time; the steady stopping keys, steady.sequencing and
  io.iteration_log with the stage; intervals without their output switch;
- simplec.*, solver.* (Source/Simplec.cpp, Source/Poisson.cpp): every simplec
  key without a stage; simplec.w_vnn_cap_slope_margin without simplec.w_vnn_cap;
  simplec.theta = frozen with the Newton stage or a restart; simplec.theta_relax_time without
  simplec.theta = relaxed_buoyant, and relaxed_buoyant without it, with the Newton stage, a restart, the
  unsteady stage or nest.names; frozen_buoyant and relaxed_buoyant outside a patch of a parent run;
  solver.poisson_precond and solver.poisson_restart with a backend that does
  not use them.

WARNING:

- steady.stage = none with nothing else asked of the run (docs/running.rst:
  "With the default steady.stage = none the run skips step 4: it writes the
  mesh and initial field and stops");
- steady.residual_drop_orders beyond what round-off allows (reference: "a
  larger value runs to steady.max_iterations");
- simplec.w_vnn_cap with simplec.w_vnn_cap_slope_margin = 0 over terrain on a
  mesh whose aspect ratio is above the tested 50.78 (the run prints the same
  warning, Source/App.cpp).

ADVICE: steady.max_iterations on a restart counts the whole run (docs/howto.rst).
"""
from ..checker import Checker, register

TESTED_ASPECT = 50.78       # Source/App.cpp, w_vnn_cap_tested_aspect
THIN_ASPECT = 50.0         # docs/howto.rst: "aspect ratio far above 10"; this checker's reading of "far"
UNSTEADY_REFUSED = ("max_iterations", "residual_drop_orders", "monitor_rel_tol", "monitor_window", "confidence",
                    "confidence_height")


def _g(x):
    return "%.6g" % x


@register
class SteadyChecker(Checker):
    """The pseudo-time iteration (steady.*): the stage, the CFL ramp, the sweeps and the stopping rules."""

    name = "steady"

    def run(self):
        c = self.ctx
        if not c.solve:
            for key in self.deck.keys("steady."):
                if key != "steady.stage":
                    self.no_effect(key, "without steady.stage = simplec")
            if c.stage == "none" and not self.given("steady.stage"):
                self.warn("steady.stage", "steady.stage is not given: the default none writes the mesh and the first "
                                          "guess and stops, without solving.",
                          "Set steady.stage = simplec to solve; give steady.stage = none to say the stop is meant.")
            return
        if not c.newton:
            for key in self.deck.keys("steady."):
                tail = key.split(".", 1)[1]
                if tail.startswith("newton_") or tail.startswith("switch_"):
                    self.no_effect(key, "without steady.stage = simplec_newton")
        self._ramp()
        self._sweeps()
        orders = c.r("steady.residual_drop_orders")
        if not c.unsteady and orders is not None and orders > 15:
            self.warn("steady.residual_drop_orders", "steady.residual_drop_orders = %s: round-off limits the reachable "
                                                     "drop to about 15 orders in double precision (6 to 7 in single), "
                                                     "so the run goes on to steady.max_iterations." % _g(orders),
                      "Use a smaller drop, or rely on steady.residual_tol.")
        if c.restart and not c.unsteady and self.given("steady.max_iterations"):
            self.advise("steady.max_iterations", "steady.max_iterations counts the iterations of the whole run, the "
                                                 "iterations before the restart included.",
                        "Give the total, not the number of further iterations.")

    def _ramp(self):
        c = self.ctx
        cfl0, cmax = c.r("steady.cfl0"), c.r("steady.cfl_max")
        if cfl0 is not None and cmax is not None and cmax < cfl0:
            self.error("steady.cfl_max", "GUST stops: steady.cfl_max = %s: must be finite and >= steady.cfl0 = %s."
                       % (_g(cmax), _g(cfl0)), "Raise steady.cfl_max or lower steady.cfl0.")
        ncfl = c.r("steady.newton_cfl_max")
        if c.newton and ncfl is not None and cmax is not None and ncfl < cmax:
            self.error("steady.newton_cfl_max", "GUST stops: steady.newton_cfl_max = %s: must be finite and >= "
                                                "steady.cfl_max = %s." % (_g(ncfl), _g(cmax)),
                       "Raise steady.newton_cfl_max.")
        conf, height = c.r("steady.confidence"), c.r("steady.confidence_height")
        if not c.unsteady and conf == 0 and height is not None and height > 0:
            self.no_effect("steady.confidence_height", "without steady.confidence > 0")

    def _sweeps(self):
        c = self.ctx
        sweep, n = c.sweep, c.i("steady.zline_sweeps")
        if sweep not in ("lusgs", "zjacobi"):
            return
        low = 1 if sweep == "lusgs" else 2
        if n is not None and not low <= n <= 100:
            self.error("steady.zline_sweeps", "GUST stops: steady.zline_sweeps = %d: must be from %d to 100 with "
                                              "steady.sweep = %s." % (n, low, sweep), "Use 2.")
        if sweep != "lusgs":
            for k in ("line_ghosts", "line_overlap", "line_order", "line_alternate"):
                why = "without steady.sweep = lusgs"
                if not self.given("steady.sweep"):
                    why += " (the default sweep of a GPU build is zjacobi)"
                self.no_effect("steady." + k, why, "Give steady.sweep = lusgs with the steady.line_* keys.")
            return
        overlap = c.i("steady.line_overlap")
        if overlap and overlap > 0 and c.s("steady.line_ghosts") != "zero":
            self.error("steady.line_overlap", "GUST stops: steady.line_overlap = %d needs steady.line_ghosts = zero."
                       % overlap, "Set steady.line_ghosts = zero, or remove the overlap.")
        if overlap and overlap > 0:
            self._overlap_layout(overlap)
        if c.flag("steady.line_alternate") and n is not None and n < 2:
            self.error("steady.line_alternate", "GUST stops: steady.line_alternate = 1 needs steady.zline_sweeps >= 2.",
                       "Use steady.zline_sweeps = 2.")
        if not self.given("steady.sweep") and any(self.given("steady.line_" + k)
                                                  for k in ("ghosts", "overlap", "order", "alternate")):
            self.advise("steady.sweep", "The steady.line_* keys need the lusgs sweep, which is the default of a CPU "
                                        "build only.", "Give steady.sweep = lusgs so that the deck also runs on a GPU "
                                                       "build.")

    def _overlap_layout(self, overlap):
        """Source/LineSweeps.cpp check_line_overlap_layout, on the run's boxes and, with grid sequencing, on every
        coarse level's (Source/App.cpp). AMReX cuts n cells into ceil(n / max) boxes of n // boxes cells or one more."""
        c = self.ctx
        if not c.n_cell or not c.max_grid:
            return
        widths = []
        for d in range(2):
            if c.n_cell[d] > 1 and c.max_grid[d] >= 1:
                boxes = -(-c.n_cell[d] // c.max_grid[d])
                widths.append((c.n_cell[d] // boxes, c.n_cell[d] % boxes))
        if not widths:
            return
        levels = c.levels if (c.sequencing and not c.restart) else 1
        for level in range(levels or 1):
            r = 2 ** level
            if level and any(w % r or extra for w, extra in widths):
                return      # boxes that do not coarsen evenly: their coarse widths are not derived here
            narrowest = min(w for w, _ in widths) // r
            if overlap > narrowest:
                where = " on grid hierarchy level %d (steady.sequencing = 1)" % level if level else ""
                self.error("steady.line_overlap", "GUST stops: steady.line_overlap = %d: must be <= the narrowest box's "
                                                  "width in x and y (%d cells%s; mesh.max_grid_size)."
                           % (overlap, narrowest, where), "Lower steady.line_overlap, or use wider boxes.")
                return


@register
class UnsteadyChecker(Checker):
    """The unsteady stage (unsteady.*): the step, the stop time, the average, and the steady keys it refuses."""

    name = "unsteady"

    def run(self):
        c = self.ctx
        if not c.unsteady:
            for key in self.deck.keys("unsteady."):
                self.no_effect(key, "without steady.stage = unsteady")
            return
        stop, dt, cfl = c.r("unsteady.stop_time"), c.r("unsteady.dt"), c.r("unsteady.cfl")
        if not self.given("unsteady.stop_time"):
            self.error("unsteady.stop_time", "GUST stops: steady.stage = unsteady needs unsteady.stop_time (s).",
                       "Give the physical time to run to.")
        has_dt = self.given("unsteady.dt")
        if cfl is not None and cfl > 0 and has_dt:
            self.error("unsteady.cfl", "GUST stops: unsteady.dt and unsteady.cfl are exclusive: give one.",
                       "Keep unsteady.dt (fixed step) or unsteady.cfl (step from the CFL number).")
        if not (cfl is not None and cfl > 0) and not has_dt:
            self.error("unsteady.dt", "GUST stops: steady.stage = unsteady needs unsteady.dt (s) or unsteady.cfl > 0.",
                       "Give the time step.")
        start = c.r("unsteady.average_start")
        if stop is not None and start is not None and self.given("unsteady.stop_time") and start >= stop:
            self.error("unsteady.average_start", "GUST stops: unsteady.average_start = %s: must be below "
                                                 "unsteady.stop_time = %s." % (_g(start), _g(stop)),
                       "Start the average before the run ends.")
        if self.given("unsteady.probes_file") and not c.s("unsteady.probes_file"):
            self.error("unsteady.probes_file", "GUST stops: unsteady.probes_file is given without a file name.",
                       "Give the file, or remove the key.")
        for k in UNSTEADY_REFUSED:
            self.no_effect("steady." + k, "with steady.stage = unsteady (the run stops at unsteady.stop_time)")
        self.no_effect("steady.probes_file", "with steady.stage = unsteady (unsteady.probes_file names the probes "
                                             "of an unsteady run)", "Use unsteady.probes_file.")
        if c.flag("steady.sequencing"):
            self.error("steady.sequencing", "GUST stops: steady.sequencing = 1 is not supported with steady.stage = "
                                            "unsteady.", "Remove steady.sequencing.")
        for key, switch, what in (("unsteady.checkpoint_interval", "io.write_checkpoint", "checkpoints"),
                                  ("unsteady.plot_interval", "io.write_plotfile", "plotfiles")):
            if c.positive(key) and not c.flag(switch):
                self.error(key, "GUST stops: %s > 0 needs %s = 1 (it writes %s)." % (key, switch, what),
                           "Set %s = 1, or remove %s." % (switch, key))
        if stop is not None and stop > 0:
            if has_dt and dt:
                self.info("unsteady.dt", "%d steps of %s s to t = %s s." % (int(stop / dt + 0.5), _g(dt), _g(stop)))


@register
class PressureChecker(Checker):
    """SIMPLEC and the pressure solve (simplec.*, solver.*): keys that need an option, and the w-row cap."""

    name = "pressure"

    def run(self):
        c = self.ctx
        if c.reads_bc:
            backend = c.s("solver.poisson_backend")
            if backend in ("mlmg", "gmres_mlmg"):
                self.no_effect("solver.poisson_precond", "with solver.poisson_backend = %s (it belongs to gmres_fft)"
                               % backend)
            if backend == "mlmg":
                self.no_effect("solver.poisson_restart", "with solver.poisson_backend = mlmg (no GMRES)")
        if not c.solve:
            for key in self.deck.keys("simplec."):
                self.no_effect(key, "without steady.stage = simplec")
            return
        cap, margin = c.r("simplec.w_vnn_cap"), c.r("simplec.w_vnn_cap_slope_margin")
        if not (cap is not None and cap > 0):
            self.no_effect("simplec.w_vnn_cap_slope_margin", "without simplec.w_vnn_cap > 0")
        if c.theta_frozen:
            if c.newton:
                self.error("simplec.theta", "GUST stops: simplec.theta = frozen is not supported with steady.stage = "
                                            "simplec_newton.", "Use steady.stage = simplec.")
            if c.restart:
                self.error("simplec.theta", "GUST stops: simplec.theta = frozen is not supported with io.restart.",
                           "Start the run from a first guess.")
            if c.first_guess == "column" and c.column_evolves and c.heat_flux_surface:
                self.error("simplec.theta", "GUST stops: simplec.theta = frozen with a column first guess that sets "
                                            "theta from a surface heat flux.", "Use simplec.theta = solve.")
        theta = c.s("simplec.theta")
        # frozen_buoyant and relaxed_buoyant hold or relax theta towards a parent run's (Source/App.cpp).
        if theta in ("frozen_buoyant", "relaxed_buoyant") and c.first_guess != "nest" and not c.restart:
            self.error("simplec.theta", "GUST stops: simplec.theta = %s needs a nested patch of a parent run "
                                        "(bc.nest_parent, init.first_guess = nest)." % theta,
                       "Use simplec.theta = solve or frozen.")
        # simplec.theta = relaxed_buoyant (Source/Simplec.cpp, Source/App.cpp): the time scale is required with it and
        # refused without it; the Newton stage, a restart, the unsteady stage and nest.names are refused.
        if theta != "relaxed_buoyant":
            self.no_effect("simplec.theta_relax_time", "unless simplec.theta = relaxed_buoyant")
        else:
            if not self.given("simplec.theta_relax_time"):
                self.error("simplec.theta", "GUST stops: simplec.theta = relaxed_buoyant needs simplec.theta_relax_time.",
                           "Give simplec.theta_relax_time (s, finite and > 0).")
            if c.newton:
                self.error("simplec.theta", "GUST stops: simplec.theta = relaxed_buoyant is not supported with "
                                            "steady.stage = simplec_newton.", "Use steady.stage = simplec.")
            if c.restart:
                self.error("simplec.theta", "GUST stops: simplec.theta = relaxed_buoyant is not supported with "
                                            "io.restart.", "Start the run from a first guess.")
            if c.unsteady:
                self.error("simplec.theta", "GUST stops: simplec.theta = relaxed_buoyant is not supported with "
                                            "steady.stage = unsteady.", "Use steady.stage = simplec.")
            if self.given("nest.names"):
                self.error("simplec.theta", "GUST stops: simplec.theta = relaxed_buoyant is not supported with "
                                            "nest.names.", "Use simplec.theta = solve in a joint run.")
        aspect = c.aspect_ratio
        if cap and cap > 0 and margin == 0 and c.has_terrain and aspect is not None and aspect > TESTED_ASPECT:
            self.warn("simplec.w_vnn_cap", "simplec.w_vnn_cap = %s with simplec.w_vnn_cap_slope_margin = 0 over terrain: "
                                           "the cell aspect ratio is at least %s, above the tested %.2f. The cap may "
                                           "make the pressure correction unstable where thin cells meet steep slopes "
                                           "(the run prints the same warning)." % (_g(cap), _g(aspect), TESTED_ASPECT),
                      "Keep the default simplec.w_vnn_cap_slope_margin = 4.")
        elif not (cap and cap > 0) and aspect is not None and aspect > THIN_ASPECT:
            self.advise("simplec.w_vnn_cap", "The first cells are thin (aspect ratio %s) and simplec.w_vnn_cap is off."
                        % _g(aspect), "Consider simplec.w_vnn_cap = 16 (docs/howto.rst, thin first cells).")
