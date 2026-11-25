"""The boundary conditions: faces, inflow profile, sounding, damping layers, the ground and nested patches.

ERROR (aborts of Source/BoundaryInputs.cpp check_boundary_inputs,
Source/InflowProfile.cpp check_inflow_profile_options, Source/Forcing.cpp
check_reference_profile, Source/SurfaceLayer.cpp, Source/ColumnGuess.cpp
check_stable_functions and Source/NestFaces.cpp check_nest_run):

- a face that is periodic when mesh.is_periodic is 0 or the reverse; a face
  type on a face that does not take it; a geostrophic top without Coriolis
  and geostrophic forcing; a wall velocity on a face that is not a no-slip
  wall, or with a normal component;
- the inflow keys of another bc.inflow_profile, bc.inflow_direction with a
  file profile, bc.inflow_file missing, and every inflow key in a run that
  uses no profile;
- sounding lists of different length, heights that do not increase, theta
  not positive; bc.zhi_theta_grad with a no-slip top;
- bc.rayleigh_zdamp at or above the top level; the Rayleigh and sponge keys
  without their layer, a layer without its coefficient, a sponge without an
  inflow or outflow face; bc.outflow_transport without an outflow face;
- the bc.surface_* keys without bc.zlo = surface_layer; bc.surface_z0 with
  bc.surface_z0_file; a heat flux on a neutral surface; the QNSE
  coefficients;
- the rules of a nested patch (parent, first guess, sequencing, walls,
  damping layers, relaxation zone, ground blend).

WARNING, each with its sentence:

- bc.inflow_z0 different from bc.surface_z0 (reference, bc.inflow_z0:
  "Normally equal to bc.surface_z0");
- an inflow face the wind of bc.inflow_direction leaves through, or an
  outflow face it enters through (docs/howto.rst: "Make the faces the wind
  enters through inflow and the faces it leaves through outflow";
  docs/boundary_conditions.rst: a zero-gradient outflow is not meant for
  backflow);
- the wind parallel to an inflow or outflow pair of faces (docs/howto.rst,
  the direction table: for a wind along x, "ylo/yhi slip or periodic");
- a stratified sounding over terrain without a Rayleigh layer
  (docs/howto.rst: "In stable air, hills launch gravity waves. Absorb them
  with a Rayleigh layer near the top").

INFO: the wind direction in both conventions and the log-law speed at 10 m
and at the top.
"""
import math

from ..checker import Checker, register
from ..context import FACES, LATERAL, SURFACE_KEYS

PARALLEL_DEG = 1.0      # a wind within this angle of a face is parallel to it


def _g(x):
    return "%.6g" % x


@register
class FacesChecker(Checker):
    """The type of each face against the periodicity and the face's place; wall velocities."""

    name = "faces"

    @classmethod
    def applies(cls, ctx):
        return ctx.reads_bc

    def run(self):
        c = self.ctx
        for o, f in enumerate(FACES):
            d = o // 2
            t, key = c.face[f], "bc." + f
            ok = True
            if c.periodic and (c.periodic[d] == 1) != (t == "periodic"):
                ok = False
                if t == "periodic":
                    self.error(key, "GUST stops: %s = periodic with mesh.is_periodic = 0 in %s (a face is periodic "
                                    "exactly when mesh.is_periodic is 1 in that direction)." % (key, "xyz"[d]),
                               "Give %s a type such as %s, or set mesh.is_periodic to 1 in %s."
                               % (key, "inflow, outflow or slip_wall" if d < 2 else "slip_wall", "xyz"[d]))
                else:
                    self.error(key, "GUST stops: %s = %s with mesh.is_periodic = 1 in %s (a face is periodic exactly "
                                    "when mesh.is_periodic is 1 in that direction)." % (key, t, "xyz"[d]),
                               "Remove %s, or set mesh.is_periodic to 0 in %s." % (key, "xyz"[d]))
            elif d < 2 and t in ("surface_layer", "geostrophic"):
                self.error(key, "GUST stops: %s = %s: only the ground (zlo) or the top (zhi) take this type." % (key, t),
                           "Use inflow, outflow, slip_wall, no_slip_wall, nest or periodic.")
            elif d == 2 and t == "nest":
                self.error(key, "GUST stops: %s = nest: only the lateral faces take this type." % key,
                           "Use the parent's %s." % key)
            elif f == "zlo" and t not in ("surface_layer", "slip_wall", "no_slip_wall"):
                self.error(key, "GUST stops: %s = %s: the ground must be surface_layer, slip_wall or no_slip_wall."
                           % (key, t), "Use one of the three.")
            elif f == "zhi" and t not in ("slip_wall", "no_slip_wall", "geostrophic"):
                self.error(key, "GUST stops: %s = %s: the top must be slip_wall, no_slip_wall or geostrophic."
                           % (key, t), "Use one of the three.")
            if ok and t == "geostrophic" and not c.geostrophic:
                self.error(key, "GUST stops: %s = geostrophic needs Coriolis with geostrophic forcing "
                                "(forcing.coriolis = 1, forcing.geostrophic = 1)." % key,
                           "Set the two forcing keys, or use slip_wall.")
            if not ok:
                continue
            wkey = key + "_wall_velocity"
            if self.given(wkey):
                wv = c.rs(wkey)
                if t != "no_slip_wall":
                    self.no_effect(wkey, "without %s = no_slip_wall" % key)
                elif wv and len(wv) == 3 and wv[d] != 0:
                    self.error(wkey, "GUST stops: %s = %s: the normal component must be 0 (the wall moves in its "
                                     "plane)." % (wkey, " ".join(c.tokens(wkey))), "Set the %s component to 0." % "uvw"[d])
        if not c.outflow_faces:
            self.no_effect("bc.outflow_transport", "without an outflow face")


@register
class InflowChecker(Checker):
    """The inflow profile's keys, and the wind direction against the inflow and outflow faces."""

    name = "inflow"

    @classmethod
    def applies(cls, ctx):
        return ctx.reads_bc

    def run(self):
        c = self.ctx
        kind = c.inflow_profile
        if kind not in ("log_law", "power_law", "file", "column"):
            return
        if kind == "column":
            self._column()
            return
        if not c.needs_profile:
            for k in ("inflow_profile", "inflow_ustar", "inflow_z0", "inflow_speed", "inflow_zref", "inflow_alpha",
                      "inflow_direction", "inflow_file", "inflow_tke_length"):
                self.no_effect("bc." + k, "without an inflow face, a damping layer or init.first_guess = profile or "
                                          "column")
            return
        for k, only in (("inflow_ustar", "log_law"), ("inflow_speed", "power_law"), ("inflow_zref", "power_law"),
                        ("inflow_alpha", "power_law"), ("inflow_file", "file")):
            if kind != only:
                self.no_effect("bc." + k, "with bc.inflow_profile = %s" % kind,
                               "Remove bc.%s, or set bc.inflow_profile = %s." % (k, only))
        if kind == "file":
            self.no_effect("bc.inflow_direction", "with bc.inflow_profile = file (the file gives u and v)")
            if not c.s("bc.inflow_file"):
                self.error("bc.inflow_file", "GUST stops: bc.inflow_file must be given for bc.inflow_profile = file.",
                           "Name the table (z u v theta, optionally k).")
            return
        z0 = c.r("bc.inflow_z0")
        sz0 = c.r("bc.surface_z0")
        # The analytic profile is imposed on inflow faces only (a nest face takes the parent's field, a periodic
        # run has no face): without one, the default bc.inflow_z0 says nothing about the deck.
        imposed = self.given("bc.inflow_z0") or bool(c.inflow_faces)
        if (imposed and c.surface_layer and z0 is not None and sz0 is not None and not self.given("bc.surface_z0_file")
                and z0 > 0 and sz0 > 0 and abs(z0 - sz0) > 1e-12 * max(z0, sz0)):
            self.warn("bc.inflow_z0", "bc.inflow_z0 = %s differs from bc.surface_z0 = %s: the inflow profile is not in "
                                      "equilibrium with the ground it flows over." % (_g(z0), _g(sz0)),
                      "The two are normally equal; keep the difference only where the upstream roughness differs.")
        self._direction()
        if kind == "log_law":
            self._log_law(z0)

    def _column(self):
        """bc.inflow_profile = column: the run's own column is the faces' profile (Source/InflowColumn.cpp)."""
        c = self.ctx
        pre = "GUST stops: bc.inflow_profile = column "
        for k in ("inflow_ustar", "inflow_speed", "inflow_zref", "inflow_alpha", "inflow_file"):
            self.no_effect("bc." + k, "with bc.inflow_profile = column (the faces take the column's table)")
        if self.given("bc.inflow_direction"):
            self.error("bc.inflow_direction", pre + "refuses bc.inflow_direction (the column gives u and v).",
                       "Remove bc.inflow_direction.")
        if self.given("bc.inflow_z0"):
            self.error("bc.inflow_z0", pre + "refuses bc.inflow_z0: the faces take the column's roughness "
                                             "(init.column_z0 with bc.surface_z0_file, else bc.surface_z0).",
                       "Remove bc.inflow_z0.")
        if self.given("bc.inflow_tke_length"):
            self.error("bc.inflow_tke_length", pre + "refuses bc.inflow_tke_length (the column's table has k).",
                       "Remove bc.inflow_tke_length.")
        if c.first_guess != "column" and not c.restart:
            self.error("bc.inflow_profile", pre + "needs init.first_guess = column (the run solves its column at "
                                                  "start-up) or io.restart (the checkpoint holds the column).",
                       "Set init.first_guess = column.")
        if not c.inflow_faces:
            self.error("bc.inflow_profile", pre + "needs an inflow face (bc.xlo, bc.xhi, bc.ylo or bc.yhi = inflow).",
                       "Without one, init.first_guess = column alone gives the first guess.")
        if c.positive("bc.rayleigh_zdamp") and any(c.rayleigh_field(f) for f in ("u", "v", "tke")):
            self.error("bc.rayleigh_fields", pre + "with a Rayleigh layer on u, v or tke: the column's layer would "
                                                   "relax towards the profile it is solving for.",
                       "Leave u, v and tke out of bc.rayleigh_fields.")
        if self.given("nest.names"):
            self.error("nest.names", pre + "is not supported with nest.names.", "Use bc.inflow_profile = file.")

    def _direction(self):
        c = self.ctx
        deg = c.r("bc.inflow_direction")
        if deg is None or not math.isfinite(deg):
            return
        cx, cy = math.cos(math.radians(deg)), math.sin(math.radians(deg))
        self.info("bc.inflow_direction", "The wind blows towards %s degrees counter-clockwise from +x: a wind from %s "
                                         "degrees (meteorological, clockwise from north, with +y north)."
                  % (_g(deg % 360.0), _g((270.0 - deg) % 360.0)))
        tol = math.sin(math.radians(PARALLEL_DEG))
        for d, comp in ((0, cx), (1, cy)):
            lo, hi = LATERAL[2 * d], LATERAL[2 * d + 1]
            open_faces = [f for f in (lo, hi) if c.face[f] in ("inflow", "outflow")]
            if not open_faces:
                continue
            if abs(comp) < tol:
                f = open_faces[0]
                self.warn("bc." + f, "bc.%s = %s with a wind along %s (bc.inflow_direction = %s): no flow crosses the "
                                     "%s faces." % (f, c.face[f], "y" if d == 0 else "x", _g(deg), "xy"[d]),
                          "Use slip_wall or periodic on the faces the wind runs along.")
                continue
            enters, leaves = (lo, hi) if comp > 0 else (hi, lo)
            if c.face[leaves] == "inflow":
                self.warn("bc." + leaves, "bc.%s = inflow, but the wind of bc.inflow_direction = %s leaves through it."
                          % (leaves, _g(deg)), "Make it outflow (and bc.%s inflow), or change bc.inflow_direction."
                          % enters)
            if c.face[enters] == "outflow":
                self.warn("bc." + enters, "bc.%s = outflow, but the wind of bc.inflow_direction = %s enters through it "
                                          "(a zero-gradient outflow is not meant for backflow)." % (enters, _g(deg)),
                          "Make it inflow (and bc.%s outflow), or change bc.inflow_direction." % leaves)

    def _log_law(self, z0):
        c = self.ctx
        ustar, kappa = c.r("bc.inflow_ustar"), (c.r("turb.kappa") if c.solve else 0.41)
        if not (ustar and z0 and kappa and ustar > 0 and z0 > 0 and kappa > 0):
            return
        text = "Log law with u* = %s m/s and z0 = %s m: %s m/s at 10 m" \
               % (_g(ustar), _g(z0), _g(ustar / kappa * math.log((10.0 + z0) / z0)))
        if c.ztop:
            text += ", %s m/s at the top (%s m)" % (_g(ustar / kappa * math.log((c.ztop + z0) / z0)), _g(c.ztop))
        self.info("bc.inflow_ustar", text + ".")


@register
class SoundingChecker(Checker):
    """The potential-temperature sounding (bc.theta_profile_z, bc.theta_profile) and the top's gradient."""

    name = "sounding"

    @classmethod
    def applies(cls, ctx):
        return ctx.reads_bc

    def run(self):
        c = self.ctx
        z, th = c.rs("bc.theta_profile_z"), c.rs("bc.theta_profile")
        if self.given("bc.theta_profile_z") or self.given("bc.theta_profile"):
            nz, nt = len(c.tokens("bc.theta_profile_z")), len(c.tokens("bc.theta_profile"))
            if nz != nt:
                key = "bc.theta_profile" if self.given("bc.theta_profile") else "bc.theta_profile_z"
                self.error(key, "GUST stops: bc.theta_profile_z has %d value%s and bc.theta_profile has %d: they need "
                                "the same number of values." % (nz, "" if nz == 1 else "s", nt),
                           "Give one potential temperature per height.")
            elif z and any(b <= a for a, b in zip(z, z[1:])):
                self.error("bc.theta_profile_z", "GUST stops: bc.theta_profile_z = %s: the heights must increase "
                                                 "strictly." % " ".join(c.tokens("bc.theta_profile_z")),
                           "List the heights from the ground up.")
        if c.face["zhi"] == "no_slip_wall":
            self.no_effect("bc.zhi_theta_grad", "with bc.zhi = no_slip_wall")
        if (c.solve and c.stratified and th and c.has_terrain and not c.positive("bc.rayleigh_zdamp")
                and c.positive("forcing.gravity") and not c.theta_frozen):
            self.warn("bc.rayleigh_zdamp", "A stratified sounding over terrain with no Rayleigh layer: the gravity waves "
                                           "the terrain launches reflect at the top.",
                      "Add bc.rayleigh_zdamp and bc.rayleigh_dampcoef (docs/howto.rst, the Rayleigh layer).")


@register
class DampingChecker(Checker):
    """The Rayleigh layer under the top and the lateral sponge."""

    name = "damping"

    @classmethod
    def applies(cls, ctx):
        return ctx.reads_bc

    def run(self):
        c = self.ctx
        zd, coef = c.r("bc.rayleigh_zdamp"), c.r("bc.rayleigh_dampcoef")
        on = zd is not None and zd > 0
        if on and c.ztop is not None and not zd < c.ztop:
            how = ""
            if c.zlevels_source == "ratio" and c.prob_hi:
                how = " (with mesh.stretch_ratio the layers end there, not at mesh.prob_hi z = %s)" % _g(c.prob_hi[2])
            self.error("bc.rayleigh_zdamp", "GUST stops: bc.rayleigh_zdamp = %s: must be in [0, the top level %s m)%s."
                       % (_g(zd), _g(c.ztop), how), "Give a depth below the top level.")
        if not on:
            for k in ("rayleigh_dampcoef", "rayleigh_fields", "rayleigh_height"):
                self.no_effect("bc." + k, "without bc.rayleigh_zdamp > 0")
        elif coef is not None and not coef > 0:
            self.error("bc.rayleigh_zdamp", "GUST stops: bc.rayleigh_zdamp > 0 needs bc.rayleigh_dampcoef > 0.",
                       "Give the layer's rate, bc.rayleigh_dampcoef (1/s).")
        elif c.ztop:
            self.info("bc.rayleigh_zdamp", "Rayleigh layer: the top %s m (%.0f %% of %s m), on %s."
                      % (_g(zd), 100.0 * zd / c.ztop, _g(c.ztop), " ".join(c.tokens("bc.rayleigh_fields"))))
        width, strength = c.r("bc.sponge_width"), c.r("bc.sponge_strength")
        sponge = width is not None and width > 0
        if sponge and not (c.inflow_faces or c.outflow_faces):
            self.error("bc.sponge_width", "GUST stops: bc.sponge_width > 0 needs an inflow or outflow face to damp "
                                          "towards.", "Remove the sponge, or give an inflow or outflow face.")
        if not sponge:
            self.no_effect("bc.sponge_strength", "without bc.sponge_width > 0")
        elif strength is not None and not strength > 0:
            self.error("bc.sponge_width", "GUST stops: bc.sponge_width > 0 needs bc.sponge_strength > 0.",
                       "Give the sponge's rate, bc.sponge_strength (1/s).")


@register
class SurfaceChecker(Checker):
    """The ground's wall law: the bc.surface_* keys."""

    name = "surface"

    @classmethod
    def applies(cls, ctx):
        return ctx.reads_bc

    def run(self):
        c = self.ctx
        if not c.surface_layer:
            for k in SURFACE_KEYS:
                self.no_effect("bc." + k, "without bc.zlo = surface_layer")
            return
        if self.given("bc.surface_z0") and self.given("bc.surface_z0_file"):
            self.error("bc.surface_z0_file", "GUST stops: bc.surface_z0 and bc.surface_z0_file cannot both be given "
                                             "(one roughness: a value or a map).", "Keep one of the two.")
        flux = c.r("bc.surface_heat_flux")
        stability = c.s("bc.surface_stability")
        if flux is not None and flux != 0 and stability in ("neutral", "surface_temperature"):
            self.error("bc.surface_heat_flux", "GUST stops: bc.surface_heat_flux is set but bc.surface_stability is "
                                               "%s; it would be ignored." % stability,
                       "Set bc.surface_stability = heat_flux, or remove the flux.")
        # SurfaceLayer.cpp, read_surface_layer_options and check_surface_layer_options: the ground temperature.
        value, ts_file = self.given("bc.surface_temperature"), self.given("bc.surface_temperature_file")
        if stability == "surface_temperature":
            if value and ts_file:
                self.error("bc.surface_temperature_file", "GUST stops: bc.surface_temperature and "
                                                          "bc.surface_temperature_file cannot both be given (one ground "
                                                          "temperature: a value or a map).", "Keep one of the two.")
            elif not value and not ts_file:
                self.error("bc.surface_stability", "GUST stops: bc.surface_stability = surface_temperature needs "
                                                   "exactly one of bc.surface_temperature (K) and "
                                                   "bc.surface_temperature_file.",
                           "Give the ground temperature, or a map of it.")
            if c.s("bc.surface_stable_functions") == "qnse":
                self.error("bc.surface_stable_functions", "GUST stops: bc.surface_stable_functions = qnse is not "
                                                          "available with bc.surface_stability = surface_temperature "
                                                          "(Businger-Dyer only).", "Remove the key.")
        else:
            for k in ("bc.surface_temperature", "bc.surface_temperature_file"):
                if self.given(k):
                    self.error(k, "GUST stops: %s is set but bc.surface_stability is %s; it would be ignored."
                               % (k, stability), "Set bc.surface_stability = surface_temperature, or remove the key.")
        q = {k: c.r("bc.surface_qnse_" + k) for k in ("am", "bm", "ah", "bh")}
        if None in q.values():
            return
        if c.s("bc.surface_stable_functions") == "qnse":
            for a, b in (("am", "bm"), ("ah", "bh")):
                if not q[b] > 0:
                    self.error("bc.surface_qnse_" + b, "GUST stops: bc.surface_qnse_%s = %s: must be > 0." % (b, _g(q[b])),
                               "Give a positive coefficient.")
                if not q[a] >= q[b]:
                    self.error("bc.surface_qnse_" + a, "GUST stops: bc.surface_qnse_%s = %s: must be >= "
                                                       "bc.surface_qnse_%s = %s." % (a, _g(q[a]), b, _g(q[b])),
                               "Keep a >= b, so that the stability function is at least 1.")
            column = c.first_guess == "column" and c.column_evolves
            if stability == "neutral" and not column:
                self.error("bc.surface_stable_functions", "GUST stops: bc.surface_stable_functions = qnse has no effect: "
                                                          "the surface is neutral and there is no stratified column.",
                           "Remove the key, or set bc.surface_stability = heat_flux.")
        else:
            for k in ("am", "bm", "ah", "bh"):
                dflt = float(c.ref.default("bc.surface_qnse_" + k)[0])
                if q[k] != dflt:
                    self.error("bc.surface_qnse_" + k, "GUST stops: bc.surface_qnse_%s is set but "
                                                       "bc.surface_stable_functions is businger_dyer; it would be "
                                                       "ignored." % k,
                               "Set bc.surface_stable_functions = qnse, or remove the coefficient.")


@register
class NestChecker(Checker):
    """A nested patch: the parent checkpoint, the nest faces, the relaxation zone and the ground blend."""

    name = "nest"

    def run(self):
        c = self.ctx
        faces = c.nest_faces if c.reads_bc else []
        first = "bc." + faces[0] if faces else ""
        guess = c.first_guess == "nest"
        parent = c.reads_bc and bool(c.s("bc.nest_parent"))
        if guess and not parent:
            self.error("bc.nest_parent", "GUST stops: init.first_guess = nest needs bc.nest_parent (the parent run's "
                                         "checkpoint).", "Name the parent's checkpoint directory.")
        elif faces and not parent:
            self.error("bc.nest_parent", "GUST stops: %s = nest needs bc.nest_parent (the parent run's checkpoint)."
                       % first, "Name the parent's checkpoint directory.")
        if parent and not guess and not faces:
            self.no_effect("bc.nest_parent", "without init.first_guess = nest or a nest face")
        if faces and not guess and not c.restart:
            self.error(first, "GUST stops: %s = nest needs init.first_guess = nest (the parent's field as the first "
                              "guess) or io.restart." % first, "Set init.first_guess = nest.")
        if (guess or faces) and c.sequencing:
            self.error("steady.sequencing", "GUST stops: a nested patch (init.first_guess = nest, a nest face) is not "
                                            "supported with steady.sequencing = 1.", "Remove steady.sequencing.")
        blend = c.r("terrain.nest_blend_width")
        if not faces:
            self.no_effect("terrain.nest_blend_width", "without a nest face")
        elif blend and blend > 0 and c.extent:
            for d in (0, 1):
                if (LATERAL[2 * d] in faces or LATERAL[2 * d + 1] in faces) and not blend < 0.5 * c.extent[d]:
                    self.error("terrain.nest_blend_width", "GUST stops: terrain.nest_blend_width = %s: must be below "
                                                           "half the patch's width in %s (%s m)."
                               % (_g(blend), "xy"[d], _g(0.5 * c.extent[d])), "Use a narrower blend.")
        if not c.reads_bc:
            return
        ratio = c.i("bc.nest_ratio")
        # GUST checks the range first, then the parent.
        if ratio is not None and ratio not in (2, 4, 8, 16):
            self.error("bc.nest_ratio", "GUST stops: bc.nest_ratio = %d: must be 2, 4, 8 or 16 (the patch's "
                                        "refinement of the parent in x and y)." % ratio, "Use 2, 4, 8 or 16.")
        elif not parent:
            self.no_effect("bc.nest_ratio", "without bc.nest_parent")
        width = c.i("bc.nest_relax_width")
        if not faces:
            self.no_effect("bc.nest_relax_width", "without a nest face")
        elif width and width > 0 and c.n_cell:
            for d in (0, 1):
                if (LATERAL[2 * d] in faces or LATERAL[2 * d + 1] in faces) and not 2 * width < c.n_cell[d]:
                    self.error("bc.nest_relax_width", "GUST stops: bc.nest_relax_width = %d: must be below half the "
                                                      "patch's cells in %s (%d cells)." % (width, "xy"[d], c.n_cell[d]),
                               "Use a narrower zone.")
        if not (width and width > 0):
            self.no_effect("bc.nest_relax_factor", "without bc.nest_relax_width > 0")
        if not faces:
            return
        if c.inflow_faces and not c.outflow_faces:
            self.error(first, "GUST stops: a nest face with an inflow face needs an outflow face (only an outflow face "
                              "can carry the difference between the two fluxes).", "Make one face outflow.")
        for f in LATERAL:
            if c.face[f] in ("slip_wall", "no_slip_wall"):
                self.error("bc." + f, "GUST stops: bc.%s = %s is not supported with a nest face: the parent's flow "
                                      "crosses that face." % (f, c.face[f]), "Use nest, outflow or (both faces) periodic.")
        if not c.inflow_faces:
            if c.positive("bc.rayleigh_zdamp") and any(c.rayleigh_field(n) for n in ("u", "v", "tke")):
                self.error("bc.rayleigh_fields", "GUST stops: with a nest face and no inflow face the Rayleigh layer may "
                                                 "damp w and theta only (u, v and tke would relax towards the "
                                                 "bc.inflow_* profile, not the parent).",
                           "Set bc.rayleigh_fields = w theta.")
            if c.positive("bc.sponge_width"):
                self.error("bc.sponge_width", "GUST stops: bc.sponge_width > 0 is not supported with a nest face and no "
                                              "inflow face.", "Remove the sponge.")
