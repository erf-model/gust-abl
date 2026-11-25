"""The first guess, the forcing, the closure, the convection scheme and the named set-ups.

ERROR:

- init.* (Source/Initializer.cpp check_init_options, Source/ColumnGuess.cpp
  check_column_guess, Source/App.cpp): keys given with init.first_guess =
  none, the profile's keys with a nest, the column's keys without the column,
  the marching column's and the steady column's own keys, a target wind
  without its height, the column's needs (a stage, a horizontal forcing, a
  surface-layer ground, one roughness, a Rayleigh layer on theta when it
  solves theta), a stage without a first guess or restart, a first guess with
  a restart;
- forcing.* (Source/Forcing.cpp): every key without a stage; geostrophic
  forcing without Coriolis or at latitude 0; forcing.ug, forcing.vg without
  it; the Coriolis keys without forcing.coriolis = 1; forcing.theta0 with
  forcing.gravity = 0;
- turb.* (Source/Closure.cpp, Source/App.cpp): every key without a stage;
  turb.Rt_min not below turb.Rt_crit; the k-L keys and the options that need
  the closure with turb.closure = none; turb.frozen_switch without
  steady.stage = simplec_newton; turb.nut_relax_end with turb.nut_relax = 1 or with steady.stage = unsteady;
  turb.shear_stencil = log_law (Source/Sequencing.cpp) without bc.zlo = surface_layer, with
  turb.wall_distance = tucker or with buildings;
- convection.* (Source/Convection.cpp, Source/App.cpp): every key without a
  stage; the blend's keys without convection.hybrid_sigma_max > 0; the blend
  without convection.hybrid_height, with the upwind limiter, with the Newton
  or unsteady stage, with io.restart, or without turb.realizability > 0; the
  smooth_eps keys without convection.limiter = van_albada;
- abl.* (Source/AblSetup.cpp): the names without a stage, deck values a name
  refuses, and what each drive and each theta set-up needs.

ADVICE: the recommended lee setting on a terrain run that has none of it
(reference: turb.realizability "0.309", convection.hybrid_height
"(recommended: 2.5 terrain heights)"; docs/howto.rst).

INFO: the Coriolis parameter and the geostrophic wind.
"""
import math

from ..checker import Checker, register
from .files import inflow_table_theta

COLUMN_KEYS = ("column_tol", "column_max_iterations", "column_stability", "column_mol", "column_surface_theta",
               "column_target_u", "column_target_v", "column_target_height", "column_rib_max", "column_rib_min",
               "column_time", "column_dt", "column_wall_k", "column_z0")
LAMINAR_KEYS = ("turb.closure", "turb.nu_molecular", "turb.alpha_molecular", "turb.tke_floor")


def _g(x):
    return "%.6g" % x


@register
class InitChecker(Checker):
    """The first guess (init.*): which keys go with which guess, and what the column guess needs."""

    name = "init"

    def run(self):
        c = self.ctx
        guess = c.first_guess
        if guess not in ("none", "profile", "column", "nest"):
            return
        tu, tv, th = c.r("init.column_target_u"), c.r("init.column_target_v"), c.r("init.column_target_height")
        target = bool(tu) or bool(tv)
        if target and not (th is not None and th > 0):
            self.error("init.column_target_height", "GUST stops: init.column_target_u, init.column_target_v need "
                                                    "init.column_target_height > 0.", "Give the target's height (m).")
        if not target and th is not None and th > 0:
            self.no_effect("init.column_target_height", "without init.column_target_u or init.column_target_v")
        if not target and (self.given("init.column_target_u") or self.given("init.column_target_v")):
            key = "init.column_target_u" if self.given("init.column_target_u") else "init.column_target_v"
            self.error(key, "GUST stops: init.column_target_u, init.column_target_v are both 0: no target wind.",
                       "Leave them out.")
        if c.s("init.column_stability") == "surface_temperature" and not c.positive("init.column_surface_theta"):
            self.error("init.column_surface_theta", "GUST stops: init.column_stability = surface_temperature needs "
                                                    "init.column_surface_theta > 0 (K).", "Give the surface's theta.")
        time, dt = c.r("init.column_time"), c.r("init.column_dt")
        marching = time is not None and time > 0
        if marching and dt and dt > 0 and time / dt > 1.0e7:
            self.error("init.column_time", "GUST stops: init.column_time / init.column_dt must be at most 1e7 steps.",
                       "Use a longer step or a shorter time.")
        if guess == "none":
            for key in self.deck.keys("init."):
                if key != "init.first_guess":
                    self.no_effect(key, "without init.first_guess = profile, column or nest")
        else:
            self._with_guess(guess, marching, target)
        if c.solve and not c.initialize and not c.restart:
            self.error("init.first_guess", "GUST stops: steady.stage = %s needs a first guess: init.first_guess = "
                                           "profile, column or nest, or io.restart." % c.stage,
                       "Set init.first_guess = profile.")
        if c.initialize and c.restart:
            self.error("io.restart", "GUST stops: init.first_guess and io.restart are exclusive (a restart keeps its "
                                     "state).", "Remove init.first_guess (and the init.* keys) on a restart.")
        if guess == "column":
            self._column(marching, target)

    def _with_guess(self, guess, marching, target):
        c = self.ctx
        if guess == "nest":
            for k in ("z_agl_min", "k_decay_height", "obrien"):
                self.no_effect("init." + k, "with init.first_guess = nest (the parent gives the field)")
        if guess != "column":
            for k in COLUMN_KEYS:
                self.no_effect("init." + k, "without init.first_guess = column")
        else:
            mode, mol = c.s("init.column_stability"), c.r("init.column_mol")
            if marching:
                for k in ("column_tol", "column_max_iterations"):
                    self.no_effect("init." + k, "with init.column_time > 0 (the column marches to that time)")
                if target:
                    self.error("init.column_time", "GUST stops: init.column_target_u, init.column_target_v need the "
                                                   "steady column (init.column_time = 0).",
                               "Remove the target, or init.column_time.")
                if mode == "mol" and mol:
                    self.error("init.column_mol", "GUST stops: init.column_mol with init.column_time > 0: the "
                                                  "time-marching column takes a fixed surface.",
                               "Use init.column_stability = heat_flux or surface_temperature.")
            else:
                for k in ("column_dt", "column_wall_k"):
                    self.no_effect("init." + k, "without init.column_time > 0")
            if self.given("init.k_decay_height"):
                self.error("init.k_decay_height", "GUST stops: init.k_decay_height would decay the column's converged k "
                                                  "again: not with init.first_guess = column.", "Remove it.")
            for k, own in (("column_mol", "mol"), ("column_surface_theta", "surface_temperature")):
                if mode != own:
                    self.no_effect("init." + k, "with init.column_stability = %s" % mode)
        if c.i("init.projection") == 0:
            for k in ("alpha_h", "alpha_v"):
                self.no_effect("init." + k, "with init.projection = 0")

    def _column(self, marching, target):
        c = self.ctx
        key = "init.first_guess"
        if not c.solve:
            self.error(key, "GUST stops: init.first_guess = column needs steady.stage = simplec or simplec_newton (the "
                            "column uses the stage's closure and forcing).", "Set steady.stage.")
            return      # without the stage the forcing is not read: the other needs cannot be judged
        pg = c.rs("forcing.pressure_gradient") or [0.0, 0.0, 0.0]
        if not (c.geostrophic or (len(pg) == 3 and (pg[0] != 0 or pg[1] != 0))):
            self.error(key, "GUST stops: init.first_guess = column needs a horizontal forcing (forcing.geostrophic = 1 "
                            "or a non-zero forcing.pressure_gradient): the column has no inflow to drive it.",
                       "Give the forcing, or use init.first_guess = profile.")
        mode, mol = c.s("init.column_stability"), c.r("init.column_mol")
        if not c.surface_layer:
            self.error(key, "GUST stops: init.first_guess = column needs bc.zlo = surface_layer (the column's wall law).",
                       "Set bc.zlo = surface_layer.")
        else:
            flux_mode, mol_flux = mode == "heat_flux", mode == "mol" and bool(mol)
            if self.given("bc.surface_z0_file") and not self.given("init.column_z0"):
                self.error("init.column_z0", "GUST stops: init.first_guess = column with bc.surface_z0_file needs "
                                             "init.column_z0: the column takes one roughness, the 3-D run keeps the map.",
                           "Give init.column_z0 (m), e.g. the map's roughness upwind of the domain.")
            if self.given("init.column_z0") and not self.given("bc.surface_z0_file"):
                self.error("init.column_z0", "GUST stops: init.column_z0 has no effect without bc.surface_z0_file (the "
                                             "column takes bc.surface_z0).", "Remove init.column_z0.")
            if c.heat_flux_surface and not flux_mode and not mol_flux:
                self.error("init.column_stability", "GUST stops: init.first_guess = column with bc.surface_stability = "
                                                    "heat_flux needs init.column_stability = heat_flux, or mol with "
                                                    "init.column_mol.", "Set init.column_stability = heat_flux.")
            if c.heat_flux_surface and mol_flux and self.given("bc.surface_heat_flux"):
                self.error("bc.surface_heat_flux", "GUST stops: bc.surface_heat_flux with init.column_stability = mol: "
                                                   "the 3-D run takes the column's flux for init.column_mol.",
                           "Leave bc.surface_heat_flux out.")
            if flux_mode and not c.heat_flux_surface:
                self.error("init.column_stability", "GUST stops: init.column_stability = heat_flux needs "
                                                    "bc.surface_stability = heat_flux.",
                           "Set bc.surface_stability = heat_flux and bc.surface_heat_flux.")
            if c.temperature_surface and (mode != "mol" or mol):
                self.error("init.column_stability", "GUST stops: init.first_guess = column with bc.surface_stability "
                                                    "= surface_temperature supports a neutral column only "
                                                    "(init.column_stability = mol with init.column_mol = 0): the "
                                                    "column does not take the 3-D run's ground temperature.",
                           "Leave init.column_stability and init.column_mol out.")
            if c.s("bc.surface_wall_k") != "dirichlet":
                self.error("bc.surface_wall_k", "GUST stops: init.first_guess = column supports bc.surface_wall_k = "
                                                "dirichlet only.", "Remove bc.surface_wall_k.")
            if c.positive("bc.surface_zref"):
                self.error("bc.surface_zref", "GUST stops: init.first_guess = column takes the wall law at the first "
                                              "cell centre (bc.surface_zref is not supported).", "Remove it.")
        evolve = c.column_evolves
        if evolve and not c.positive("forcing.gravity"):
            self.error("forcing.gravity", "GUST stops: init.column_stability with a non-neutral surface needs buoyancy "
                                          "(forcing.gravity > 0).", "Remove forcing.gravity = 0.")
        for k in ("init.column_rib_max", "init.column_rib_min"):
            if not evolve:
                self.no_effect(k, ": the column's surface is neutral (init.column_stability)")
            elif marching and mode == "heat_flux":
                self.no_effect(k, "with init.column_time > 0 and init.column_stability = heat_flux")
        theta_layer = c.positive("bc.rayleigh_zdamp") and c.rayleigh_field("theta")
        if evolve and not marching and not theta_layer:
            self.error("bc.rayleigh_zdamp", "GUST stops: init.column_stability with a non-neutral surface needs a "
                                            "Rayleigh layer on theta (bc.rayleigh_zdamp > 0, theta in "
                                            "bc.rayleigh_fields), so that the column can be steady.",
                       "Add the layer, or march the column (init.column_time).")
        if target:
            tkey = "init.column_target_u" if self.given("init.column_target_u") else "init.column_target_v"
            own_column = c.inflow_profile == "column"
            if c.positive("bc.rayleigh_zdamp") and c.rayleigh_field("u"):
                self.error("bc.rayleigh_fields", "GUST stops: init.column_target_u, init.column_target_v with a Rayleigh "
                                                 "layer on u and v: the layer relaxes towards the deck's inflow "
                                                 "profile, which the target does not change.",
                           "Leave u and v out of bc.rayleigh_fields.")
            if c.positive("bc.sponge_width") and not own_column:
                self.error("bc.sponge_width", "GUST stops: init.column_target_u, init.column_target_v with a sponge "
                                              "(bc.sponge_width).",
                           "Remove the sponge, or set bc.inflow_profile = column.")
            if c.inflow_faces and not own_column:
                self.error(tkey, "GUST stops: init.column_target_u, init.column_target_v with an inflow face: the face "
                                 "keeps the deck's inflow profile.",
                           "Set bc.inflow_profile = column (the faces take the target column).")
            if evolve and not c.heat_flux_surface:
                self.error(tkey, "GUST stops: init.column_target_u, init.column_target_v with a stratified column need "
                                 "bc.surface_stability = heat_flux in the 3-D run.",
                           "Set bc.surface_stability = heat_flux.")


@register
class ForcingChecker(Checker):
    """Buoyancy, Coriolis, the geostrophic wind and the pressure gradient (forcing.*)."""

    name = "forcing"

    def run(self):
        c = self.ctx
        if not c.solve:
            for key in self.deck.keys("forcing."):
                self.no_effect(key, "without steady.stage = simplec")
            return
        cor, geo = c.flag("forcing.coriolis"), c.flag("forcing.geostrophic")
        lat = c.r("forcing.latitude")
        if geo and not cor:
            self.error("forcing.geostrophic", "GUST stops: forcing.geostrophic = 1 needs forcing.coriolis = 1 (the "
                                              "forcing is f times the geostrophic wind).", "Set forcing.coriolis = 1.")
        elif geo and lat == 0:
            self.error("forcing.latitude", "GUST stops: geostrophic forcing at forcing.latitude = 0: f = 0, so it would "
                                           "not drive the flow.", "Give the site's latitude (degrees).")
        ug, vg = c.r("forcing.ug"), c.r("forcing.vg")
        if not geo and (ug or vg):
            self.error("forcing.ug" if ug else "forcing.vg", "GUST stops: forcing.ug, forcing.vg are set but "
                                                             "forcing.geostrophic is off; they would be ignored.",
                       "Set forcing.geostrophic = 1, or remove them.")
        if c.r("forcing.gravity") == 0:
            self.no_effect("forcing.theta0", "with forcing.gravity = 0")
        if not cor:
            for k in ("latitude", "rotation_period", "coriolis_3d", "geostrophic", "ug", "vg"):
                self.no_effect("forcing." + k, "without forcing.coriolis = 1")
        elif not geo:
            for k in ("ug", "vg"):
                self.no_effect("forcing." + k, "without forcing.geostrophic = 1")
        f = c.coriolis_parameter
        if cor and f is not None:
            text = "Coriolis parameter f = %s 1/s at latitude %s" % (_g(f), _g(lat if lat is not None else 0.0))
            if geo and ug is not None and vg is not None:
                text += "; geostrophic wind %s m/s towards %s degrees from +x" \
                        % (_g(math.hypot(ug, vg)), _g(math.degrees(math.atan2(vg, ug)) % 360.0))
            self.info("forcing.coriolis", text + ".")


@register
class TurbulenceChecker(Checker):
    """The closure (turb.*): keys that need the stage or the k-L closure; the recommended lee setting."""

    name = "turbulence"

    def run(self):
        c = self.ctx
        if not c.solve:
            for key in self.deck.keys("turb."):
                self.no_effect(key, "without steady.stage = simplec")
            return
        crit, rmin = c.r("turb.Rt_crit"), c.r("turb.Rt_min")
        if crit is not None and rmin is not None and not rmin < crit:
            self.error("turb.Rt_min", "GUST stops: turb.Rt_min = %s: must be below turb.Rt_crit = %s."
                       % (_g(rmin), _g(crit)), "Keep turb.Rt_min < turb.Rt_crit.")
        if not c.newton:
            self.no_effect("turb.frozen_switch", "without steady.stage = simplec_newton")
        if c.closure == "none":
            self._laminar()
        elif c.closure == "kl":
            relax = c.r("turb.nut_relax")
            if c.r("turb.nut_relax_end") and (relax is None or relax >= 1):
                self.no_effect("turb.nut_relax_end", "with turb.nut_relax = 1 (nu_t is not relaxed)",
                               "Set turb.nut_relax below 1, or remove turb.nut_relax_end.")
            if c.unsteady and c.r("turb.nut_relax_end"):
                self.error("turb.nut_relax_end", "GUST stops: turb.nut_relax_end is not supported with steady.stage = "
                                                 "unsteady.", "Remove turb.nut_relax_end.")
            self._log_shear()
            self._lee()

    def _log_shear(self):
        c = self.ctx
        if c.s("turb.shear_stencil") != "log_law":
            return
        stops = "GUST stops: turb.shear_stencil = log_law "
        if c.buildings:
            self.error("turb.shear_stencil", stops + "does not support buildings.",
                       "Use turb.shear_stencil = central with buildings.")
        if not c.surface_layer:
            self.error("turb.shear_stencil", stops + "needs bc.zlo = surface_layer (the ground's roughness length).",
                       "Use bc.zlo = surface_layer, or turb.shear_stencil = central.")
        if c.s("turb.wall_distance") == "tucker":
            self.error("turb.shear_stencil", stops + "needs turb.wall_distance = terrain_height.",
                       "Use turb.wall_distance = terrain_height, or turb.shear_stencil = central.")

    def _laminar(self):
        c = self.ctx
        if not c.positive("turb.nu_molecular"):
            self.error("turb.nu_molecular", "GUST stops: turb.nu_molecular must be > 0 with turb.closure = none (a "
                                            "laminar run has no other viscosity).", "Give the viscosity (m^2/s).")
        for key in self.deck.keys("turb."):
            if key not in LAMINAR_KEYS and c.ref.has(key):
                self.no_effect(key, "with turb.closure = none (it belongs to the k-L closure)")
        needs = "needs turb.closure = kl"
        if c.canopy:
            self.error("canopy.type", "GUST stops: canopy.type %s (the canopy acts through the closure's k terms)."
                       % needs, "Remove the canopy, or use turb.closure = kl.")
        if c.surface_layer:
            self.error("bc.zlo", "GUST stops: bc.zlo = surface_layer %s (a laminar run has no wall law)." % needs,
                       "Use no_slip_wall or slip_wall.")
        if c.first_guess == "column":
            self.error("init.first_guess", "GUST stops: init.first_guess = column %s." % needs, "Use profile.")
        if c.newton:
            self.error("steady.stage", "GUST stops: steady.stage = simplec_newton is not supported with turb.closure = "
                                       "none.", "Use steady.stage = simplec.")
        if c.buildings:
            form = c.s("buildings.wall_form")
            if form in ("fraction_stress", "roof_stress"):
                self.error("buildings.wall_form", "GUST stops: buildings.wall_form = %s %s (a wall law)." % (form, needs),
                           "Use buildings.wall_form = drag.")
            tke = c.s("buildings.wall_tke") or ("on" if form == "fraction_stress" else "off")
            if tke == "on":
                self.error("buildings.wall_tke", "GUST stops: buildings.wall_tke = on %s." % needs,
                           "Set buildings.wall_tke = off (with buildings.wall_form = drag).")

    def _lee(self):
        c = self.ctx
        if not (c.terrain_type in ("xyz", "erf_grid", "agnesi") and not c.newton):
            return
        if c.realizability or c.hybrid or self.given("turb.nut_relax"):
            if c.realizability and not c.hybrid and not c.unsteady and not c.restart \
                    and c.s("convection.limiter") != "upwind":
                self.advise("convection.hybrid_sigma_max", "turb.realizability = %s without the upwind blend of the "
                                                           "recommended lee setting." % _g(c.realizability),
                            "Consider convection.hybrid_sigma_max = 0.1 with convection.hybrid_height of about 2.5 "
                            "terrain heights.")
            return
        if c.unsteady:
            self.advise("turb.realizability", "A terrain run with none of the recommended lee setting.",
                        "Consider turb.realizability = 0.309 and turb.nut_relax = 0.5 (steady.stage = unsteady refuses "
                        "the upwind blend and turb.nut_relax_end).")
        else:
            self.advise("turb.realizability", "A terrain run with none of the recommended lee setting.",
                        "Consider turb.realizability = 0.309, turb.nut_relax = 0.5, convection.hybrid_sigma_max = 0.1 "
                        "and convection.hybrid_height of about 2.5 terrain heights (docs/howto.rst).")


@register
class ConvectionChecker(Checker):
    """The convection scheme (convection.*): the limiter and the upwind blend."""

    name = "convection"

    def run(self):
        c = self.ctx
        if not c.solve:
            for key in self.deck.keys("convection."):
                self.no_effect(key, "without steady.stage = simplec")
            return
        key = "convection.hybrid_sigma_max"
        limiter = c.s("convection.limiter")
        if c.hybrid:
            if not self.given("convection.hybrid_height"):
                self.error("convection.hybrid_height", "GUST stops: convection.hybrid_sigma_max > 0 needs "
                                                       "convection.hybrid_height > 0 (m).",
                           "Give the height the blend fades over (about 2.5 terrain heights).")
            if limiter == "upwind":
                self.no_effect(key, "with convection.limiter = upwind")
            if c.newton:
                self.error(key, "GUST stops: convection.hybrid_sigma_max is not supported with steady.stage = "
                                "simplec_newton.", "Remove the blend, or use steady.stage = simplec.")
            if not (c.realizability and c.realizability > 0):
                self.error(key, "GUST stops: convection.hybrid_sigma_max needs turb.realizability > 0 (the blend "
                                "follows the eddy-viscosity cap).", "Set turb.realizability = 0.309.")
            if c.unsteady:
                self.error(key, "GUST stops: convection.hybrid_sigma_max is not supported with steady.stage = unsteady.",
                           "Remove the blend.")
            if c.restart:
                self.error(key, "GUST stops: convection.hybrid_sigma_max is not supported with io.restart (the blend is "
                                "not in the checkpoint).", "Remove the blend on the restart.")
        else:
            for k in ("hybrid_height", "hybrid_relax", "hybrid_freeze"):
                self.no_effect("convection." + k, "without convection.hybrid_sigma_max > 0")
        if limiter in ("van_leer", "upwind"):
            for k in ("smooth_eps_velocity", "smooth_eps_theta"):
                self.no_effect("convection." + k, "without convection.limiter = van_albada")


@register
class AblChecker(Checker):
    """The named set-ups abl.drive and abl.theta: what each needs and refuses."""

    name = "abl"

    @classmethod
    def applies(cls, ctx):
        return bool(ctx.abl_drive or ctx.abl_theta)

    def run(self):
        c = self.ctx
        for key, message in c.abl_refused:
            self.error(key, "GUST stops: %s." % message, "Remove %s, or the named set-up." % key)
        if not c.abl_valid or c.abl_refused:
            return
        if not c.solve:
            for key in ("abl.drive", "abl.theta"):
                self.no_effect(key, "without steady.stage = simplec")
            return
        self.theta_air = c.rs("bc.theta_profile") or []
        self.known = True
        if c.inflow_profile == "file":
            table = inflow_table_theta(c)
            self.known = table is not None
            self.theta_air = table or []
        if c.abl_drive:
            self._drive(c.abl_drive)
        if c.abl_drive == "couette" and c.abl_theta != "neutral":
            self.error("abl.drive", "GUST stops: abl.drive = couette needs abl.theta = neutral.",
                       "Set abl.theta = neutral.")
        if c.abl_theta:
            self._theta(c.abl_theta)

    def _need(self, key, name, value, what, advice=""):
        self.error(key, "GUST stops: %s = %s %s." % (name, value, what), advice)

    @staticmethod
    def _uniform(v):
        return all(x == v[0] for x in v)

    def _drive(self, drive):
        c = self.ctx
        need = lambda key, what, advice="": self._need(key, "abl.drive", drive, what, advice)   # noqa: E731
        pg = c.rs("forcing.pressure_gradient") or [0.0, 0.0, 0.0]
        horizontal = len(pg) == 3 and (pg[0] != 0 or pg[1] != 0)
        ug, vg = c.r("forcing.ug"), c.r("forcing.vg")
        if drive in ("geostrophic_column", "geostrophic_inflow"):
            if not (c.flag("forcing.coriolis") and c.flag("forcing.geostrophic")):
                need("forcing.geostrophic", "needs forcing.coriolis = 1 and forcing.geostrophic = 1")
            # bc.inflow_profile = column: a target wind stands for the forcing, and a restart takes it from the record.
            own_column = c.inflow_profile == "column"
            can_target = drive == "geostrophic_column" or own_column
            target = ((can_target and (c.r("init.column_target_u") or c.r("init.column_target_v")))
                      or (own_column and bool(c.restart)))
            if not ug and not vg and not target:
                need("forcing.ug", "needs a wind: forcing.ug, forcing.vg"
                     + (" or init.column_target_u, init.column_target_v" if can_target else ""),
                     "Give the geostrophic wind.")
        elif c.flag("forcing.coriolis"):
            need("forcing.coriolis", "refuses forcing.coriolis = 1", "Use abl.drive = geostrophic_column or "
                                                                     "geostrophic_inflow.")
        if drive != "pressure_driven" and any(x != 0 for x in pg):
            need("forcing.pressure_gradient", "refuses a non-zero forcing.pressure_gradient",
                 "Use abl.drive = pressure_driven.")
        if drive in ("inflow_only", "couette", "geostrophic_inflow") and not c.inflow_faces:
            need("abl.drive", "needs an inflow face (bc.xlo ... bc.yhi = inflow)", "Give an inflow face.")
        if drive == "couette":
            self._couette(need)
        elif drive == "pressure_driven":
            if not horizontal:
                need("forcing.pressure_gradient", "needs a non-zero horizontal forcing.pressure_gradient")
            if c.face["zhi"] != "slip_wall":
                need("bc.zhi", "needs bc.zhi = slip_wall (the force balances the wall stress alone)")
            damps = [c.rayleigh_damps("u"), c.rayleigh_damps("v")]
            if self._uniform(self.theta_air):
                if any(damps):
                    need("bc.rayleigh_fields", "with uniform theta refuses a Rayleigh layer on u or v",
                         "Leave u and v out of bc.rayleigh_fields.")
            elif not all(damps):
                need("bc.rayleigh_fields", "with a stratified theta needs a Rayleigh layer on u and v",
                     "Add a Rayleigh layer with u and v in bc.rayleigh_fields.")
        elif drive == "geostrophic_column":
            if c.first_guess != "column" and not c.restart:
                need("init.first_guess", "needs init.first_guess = column or io.restart")
            if c.inflow_faces and c.inflow_profile not in ("file", "column"):
                need("bc.inflow_profile", "with an inflow face needs bc.inflow_profile = file (the column's table) or "
                                          "column")

    def _couette(self, need):
        c = self.ctx
        lid = c.rs("bc.zhi_wall_velocity") or [0.0, 0.0, 0.0]
        if c.face["zhi"] != "no_slip_wall" or len(lid) != 3 or (lid[0] == 0 and lid[1] == 0):
            need("bc.zhi", "needs bc.zhi = no_slip_wall with a moving lid (bc.zhi_wall_velocity)")
        if c.inflow_profile != "log_law":
            need("bc.inflow_profile", "needs bc.inflow_profile = log_law")
        if c.s("bc.inflow_tke_length") != "unbounded":
            need("bc.inflow_tke_length", "needs bc.inflow_tke_length = unbounded")
        lmax = c.r("turb.max_geom_lscale")
        if lmax is not None and lmax < 1.0e6:
            need("turb.max_geom_lscale", "needs turb.max_geom_lscale >= 1e6 (the inflow's length is not capped)")
        z0, sz0 = c.r("bc.inflow_z0"), c.r("bc.surface_z0")
        if not c.surface_layer:
            need("bc.zlo", "needs bc.zlo = surface_layer")
        elif self.given("bc.surface_z0_file") or (z0 is not None and sz0 and abs(z0 - sz0) > 1.0e-5 * sz0):
            need("bc.inflow_z0", "needs bc.inflow_z0 equal to a uniform bc.surface_z0")
        deg = c.r("bc.inflow_direction")
        if len(lid) == 3 and deg is not None and (lid[0] != 0 or lid[1] != 0):
            turn = math.remainder(math.degrees(math.atan2(lid[1], lid[0])) - deg, 360.0)
            if abs(turn) > 0.1:
                need("bc.zhi_wall_velocity", "needs bc.zhi_wall_velocity along bc.inflow_direction (within 0.1 degrees)")

    def _theta(self, theta):
        c = self.ctx
        need = lambda key, what, advice="": self._need(key, "abl.theta", theta, what, advice)   # noqa: E731
        air, ref = self.theta_air, c.rs("bc.theta_profile") or []
        if theta != "neutral" and c.r("forcing.gravity") == 0:
            need("forcing.gravity", "refuses forcing.gravity = 0 (theta would not act on the flow)")
        if theta != "frozen" and c.theta_frozen:
            need("simplec.theta", "refuses simplec.theta = frozen", "Use abl.theta = frozen.")
        if theta != "surface_flux" and c.heat_flux_surface:
            need("bc.surface_stability", "refuses bc.surface_stability = heat_flux", "Use abl.theta = surface_flux.")
        if theta != "surface_flux" and c.temperature_surface:
            need("bc.surface_stability", "refuses bc.surface_stability = surface_temperature",
                 "Use abl.theta = surface_flux.")
        column = c.first_guess == "column"
        mode, mol = c.s("init.column_stability"), c.r("init.column_mol")
        if theta != "surface_flux" and column and (mode != "mol" or mol):
            need("init.column_stability", "refuses a stratified column first guess", "Use abl.theta = surface_flux.")
        zgrad = c.r("bc.zhi_theta_grad")
        if theta == "neutral":
            if not self._uniform(air) or not self._uniform(ref):
                need("bc.theta_profile", "needs a uniform bc.theta_profile and, with bc.inflow_profile = file, a "
                                         "uniform table")
            if zgrad:
                need("bc.zhi_theta_grad", "refuses a non-zero bc.zhi_theta_grad")
            if c.canopy and c.s("canopy.heat") == "fixed":
                need("canopy.heat", "refuses canopy.heat (the leaves would heat the air)")
        elif theta == "capped":
            if self.known and (self._uniform(air) or any(b < a for a, b in zip(air, air[1:]))):
                need("bc.theta_profile", "needs a bc.theta_profile (or bc.inflow_file table) that rises with height and "
                                         "never falls")
            if not (zgrad is not None and zgrad > 0):
                need("bc.zhi_theta_grad", "needs bc.zhi_theta_grad > 0")
            if not (c.rayleigh_damps("w") and c.rayleigh_damps("theta")):
                need("bc.rayleigh_zdamp", "needs a Rayleigh layer on w and theta (bc.rayleigh_zdamp, "
                                          "bc.rayleigh_dampcoef, bc.rayleigh_fields)")
        elif theta == "frozen":
            if self.known and self._uniform(air):
                need("bc.theta_profile", "needs a stratified bc.theta_profile or bc.inflow_file table",
                     "With a uniform theta use abl.theta = neutral.")
            if c.abl_drive == "geostrophic_column" and c.inflow_profile not in ("file", "column"):
                need("bc.inflow_profile", "with abl.drive = geostrophic_column needs bc.inflow_profile = file or column")
        else:
            if not (c.heat_flux_surface or c.temperature_surface):
                need("bc.surface_stability", "needs bc.zlo = surface_layer and bc.surface_stability = heat_flux or "
                                             "surface_temperature")
            column_flux = column and mode == "mol" and bool(mol)
            flux = c.r("bc.surface_heat_flux")
            if not c.temperature_surface and not column_flux and not (self.given("bc.surface_heat_flux") and flux):
                need("bc.surface_heat_flux", "needs a non-zero bc.surface_heat_flux or init.column_mol")
            marching = column and c.positive("init.column_time")
            if not c.outflow_faces and not c.rayleigh_damps("theta") and not marching:
                need("bc.rayleigh_zdamp", "needs a heat sink: an outflow face, a Rayleigh layer on theta or "
                                          "init.column_time > 0")
