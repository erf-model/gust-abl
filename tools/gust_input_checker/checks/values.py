"""Number, type and range of the values of each key.

One row per key of Tests/inputs_reference: the kind of value, how many values
the key takes, and the range GUST accepts. A key outside its row makes GUST
(or AMReX, for a value that is not a number) stop while it reads the inputs.

The ranges are those of the range checks in Source/ (Inputs.cpp validate,
BoundaryInputs.cpp check_boundary_inputs, Initializer.cpp check_init_options,
Driver.cpp check_steady_options, Sequencing.cpp, Forcing.cpp, Closure.cpp,
Convection.cpp, HybridConvection.cpp, Simplec.cpp, Poisson.cpp,
TimeDriver.cpp, SurfaceLayer.cpp, Canopy.cpp, Buildings.cpp, AblSetup.cpp,
Newton.cpp, Continuation.cpp), which docs/inputs.rst and the reference repeat
key by key. Where a bound depends on the precision of the build (a multiple of
machine epsilon), the row keeps the part that holds in both precisions. Rules
that involve other keys are in the area checks.

The row also covers the keys the reference marks ``required`` with no
condition (mesh.n_cell, mesh.prob_hi).
"""
import math

from ..checker import Checker, register
from ..deck import Unparsed, to_int, to_real
from ..reference import template

FACE_TYPES = "periodic,inflow,outflow,slip_wall,no_slip_wall,surface_layer,geostrophic,nest"

# key  kind  count  [range]
#   kind: int, real, flag (0 or 1), str, or a comma-separated list of the allowed words
#   count: n (exactly), a-b, n+ (at least)
#   range: alternatives joined by '|'; an alternative is a number or comparisons joined by '&'
SPEC_TEXT = """
mesh.n_cell int 3 >=1
mesh.prob_lo real 3
mesh.prob_hi real 3
mesh.is_periodic flag 3
mesh.max_grid_size int 1-2 >=1
mesh.dz0 real 1 >=0
mesh.stretch_ratio real 1 >=1
mesh.zlevels_file str 1
mesh.hierarchy_levels int 1 >=1&<=10
terrain.type flat,agnesi,xyz,erf_grid 1
terrain.file str 1
terrain.smoothing int 1 >=0&<=2
terrain.gamma_m real 1 >0&<1
terrain.hill_height real 1 >0
terrain.hill_half_width real 1 >0
terrain.hill_shape ridge_x,ridge_y,round 1
terrain.idw_neighbors int 1 >=1
terrain.idw_power real 1 >0
terrain.hierarchy_passes int 1 >=0&<=64
terrain.nest_blend_width real 1 >=0
terrain.representation body_fitted,immersed_forcing 1
terrain.flat_height real 1 >0
terrain.small_volfrac real 1 >=1e-4&<0.5
init.first_guess none,profile,column,nest 1
init.projection flag 1
init.alpha_h real 1 >0
init.alpha_v real 1 >0
init.z_agl_min real 1 >=0
init.k_decay_height real 1 >=0
init.obrien flag 1
init.column_tol real 1 >0
init.column_max_iterations int 1 >=1
init.column_stability mol,heat_flux,surface_temperature 1
init.column_mol real 1 0|>=1|<=-1
init.column_surface_theta real 1 >=0
init.column_target_u real 1
init.column_target_v real 1
init.column_target_height real 1 >=0
init.column_rib_max real 1 >0
init.column_rib_min real 1 <0
init.column_time real 1 >=0
init.column_dt real 1 >0
init.column_wall_k fixed,solved 1
bc.xlo FACES 1
bc.xhi FACES 1
bc.ylo FACES 1
bc.yhi FACES 1
bc.zlo FACES 1
bc.zhi FACES 1
bc.xlo_wall_velocity real 3
bc.xhi_wall_velocity real 3
bc.ylo_wall_velocity real 3
bc.yhi_wall_velocity real 3
bc.zlo_wall_velocity real 3
bc.zhi_wall_velocity real 3
bc.inflow_profile log_law,power_law,file 1
bc.inflow_ustar real 1 >0
bc.inflow_z0 real 1 >0
bc.inflow_speed real 1 >=0
bc.inflow_zref real 1 >0
bc.inflow_alpha real 1 >0
bc.inflow_direction real 1
bc.inflow_file str 1
bc.inflow_tke_length closure,unbounded 1
bc.theta_profile_z real 1+
bc.theta_profile real 1+ >0
bc.zhi_theta_grad real 1
bc.rayleigh_zdamp real 1 >=0
bc.rayleigh_dampcoef real 1 >=0
bc.rayleigh_fields u,v,w,theta,tke 1+
bc.rayleigh_height absolute,above_ground 1
bc.sponge_width real 1 >=0
bc.sponge_strength real 1 >=0
bc.outflow_transport mirror,from_velocity 1
bc.surface_z0 real 1 >0
bc.surface_z0_file str 1
bc.surface_zref real 1
bc.surface_stability neutral,heat_flux 1
bc.surface_heat_flux real 1
bc.surface_wall_k dirichlet,source 1
bc.surface_denominator_floor log_law,one 1
bc.surface_stable_functions businger_dyer,qnse 1
bc.surface_qnse_am real 1
bc.surface_qnse_bm real 1
bc.surface_qnse_ah real 1
bc.surface_qnse_bh real 1
bc.nest_parent str 1
bc.nest_relax_width int 1 >=0
bc.nest_relax_factor real 1 >0
bc.nest_ratio int 1
canopy.type none,patches,map 1
canopy.file str 1
canopy.height_file str 1
canopy.lai_file str 1
canopy.cd real 1 >0
canopy.cd_file str 1
canopy.profile uniform,lalic_mihailovic 1
canopy.laimax real 1
canopy.drag_diagonal derivative,picard,none 1
canopy.tke on,off 1
canopy.beta_p real 1 >=0&<=1
canopy.beta_d real 1 >=0
canopy.length_alpha real 1 >=0&<=1
canopy.heat off,fixed 1
canopy.leaf_theta real 1 >=150&<=400
canopy.heat_frozen skip,blend 1
buildings.type none,immersed_forcing 1
buildings.file str 1
buildings.format erf,xyz 1
buildings.force_file str 1
buildings.small_volfrac real 1 >=0&<0.5
buildings.z0 real 1 >0
buildings.wall_form fraction_stress,roof_stress,drag 1
buildings.wall_speed open_faces,all_faces 1
buildings.wall_tke on,off 1
buildings.theta_mode none,ambient,fixed 1
buildings.theta_surface real 1 >0
buildings.limiter_zone real 1 0|>=0.5&<=20
abl.drive inflow_only,couette,pressure_driven,geostrophic_column,geostrophic_inflow 1
abl.theta neutral,capped,frozen,surface_flux 1
forcing.gravity real 1 >=0
forcing.theta0 real 1 >0
forcing.coriolis flag 1
forcing.latitude real 1 >=-90&<=90
forcing.rotation_period real 1 >0
forcing.coriolis_3d flag 1
forcing.geostrophic flag 1
forcing.ug real 1
forcing.vg real 1
forcing.pressure_gradient real 3
turb.closure kl,none 1
turb.Cmu0 real 1 >0
turb.Cb real 1 >0
turb.Rt_crit real 1 <=0
turb.Rt_min real 1 >-3.598
turb.max_geom_lscale real 1 >0
turb.tke_floor real 1 0|>=2.220446049250313e-16
turb.sigma_k real 1 >0
turb.nu_molecular real 1 >=0
turb.alpha_molecular real 1 >=0
turb.kappa real 1 >0
turb.nut_relax real 1 >0&<=1
turb.nut_relax_end real 1 0|>0&<1
turb.frozen_switch flag 1
turb.ri_gate real 1 0|>=1e-6&<=1
turb.realizability real 1 0|>=0.01&<=10
turb.face_mean arithmetic,logarithmic 1
turb.wall_distance terrain_height,tucker 1
convection.limiter van_leer,van_albada,upwind 1
convection.smooth_eps_velocity real 1 >0
convection.smooth_eps_theta real 1 >0
convection.hybrid_sigma_max real 1 >=0&<=1
convection.hybrid_height real 1 >0
convection.hybrid_relax real 1 >0&<=1
convection.hybrid_freeze real 1 >=0&<1
steady.stage none,simplec,simplec_newton,unsteady 1
steady.max_iterations int 1 >=1
steady.residual_tol real 1 >0
steady.residual_drop_orders real 1 >=0&<=308
steady.residual_ref_iterations int 1 >=1
steady.monitor_rel_tol real 1 >0
steady.monitor_window int 1 >=2
steady.confidence real 1 >=0&<=1
steady.confidence_height real 1 >=0
steady.cfl0 real 1 >0
steady.cfl_max real 1 >0
steady.ser_exponent real 1 >=0
steady.ser_max_growth real 1 >=1
steady.dtau_floor real 1 >=0
steady.switch_window int 1 >=1
steady.switch_residual real 1 >0&<=1
steady.switch_back_growth real 1 >1
steady.switch_monitor_tol real 1 >=0&<1
steady.newton_fallback_failures int 1 >=1
steady.newton_cfl_max real 1 >0
steady.sweep zjacobi,lusgs 1
steady.zline_sweeps int 1 >=1&<=100
steady.line_ghosts lagged,zero 1
steady.line_overlap int 1 >=0&<=8
steady.line_order sequential,diagonal 1
steady.line_alternate flag 1
steady.sequencing flag 1
steady.sequencing_drop_orders real 1 >0
steady.sequencing_monitor_tol real 1 >=0&<1
steady.probes_file str 1
steady.sequencing_hold flag 1
steady.sequencing_cfl_start cfl0,previous 1
steady.machine local,kestrel-cpu,kestrel-gpu 1
steady.nodes int 1 >=1
steady.newton_gmres_restart int 1 >=5&<=200
steady.newton_fd_eps_scale real 1 >0&<=1e-2
steady.newton_krylov_rtol real 1 >0&<1
steady.newton_mem_budget_gb real 1 >0&<=80
steady.newton_gmres_max_iters int 1 >=1
steady.newton_preconditioner identity,simplec 1
steady.newton_armijo_c real 1 >0&<0.5
steady.newton_backtrack real 1 >0&<1
steady.newton_lambda_min real 1 >0&<1
steady.newton_cfl_cut real 1 >0&<1
steady.newton_max_retries int 1 >=1
steady.newton_dtau_inf_at_cap flag 1
steady.newton_ptc_growth real 1 >=1
unsteady.stop_time real 1 >0
unsteady.dt real 1 >0
unsteady.cfl real 1 >=0
unsteady.inner_drop_orders real 1 >0&<=15
unsteady.max_inner int 1 >=1
unsteady.average_start real 1 >=0
unsteady.plot_interval int 1 >=0
unsteady.checkpoint_interval int 1 >=0
unsteady.probes_file str 1
simplec.pressure_reltol real 1 >0&<1
simplec.alpha_p real 1 >0&<=1
simplec.rotation_relax real 1 >0&<=1
simplec.w_vnn_cap real 1 0|>=16
simplec.w_vnn_cap_slope_margin real 1 0|>=1
simplec.theta solve,frozen 1
solver.poisson_backend gmres_fft,mlmg,gmres_mlmg 1
solver.poisson_reltol real 1 >0&<1
solver.poisson_abstol real 1 >=0
solver.poisson_max_iter int 1 >=1
solver.poisson_restart int 1 >=1
solver.poisson_precond fft,none 1
solver.poisson_verbose int 1 >=0
io.plot_file str 1
io.chk_file str 1
io.write_plotfile flag 1
io.write_checkpoint flag 1
io.turb_plot_file str 1
io.iteration_log table,full 1
io.column_guess_file str 1
io.restart str 1
nest.names str 1+
nest.feedback none,defect,relax,inject 1
nest.feedback_factor real 1 >0
nest.feedback_flux outer,nest 1
nest.damping w_theta,outer 1
nest.<name>.lo real 2
nest.<name>.hi real 2
"""

_OPS = ((">=", lambda x, a: x >= a), ("<=", lambda x, a: x <= a), (">", lambda x, a: x > a),
        ("<", lambda x, a: x < a))


class Row:
    def __init__(self, key, kind, count, rng):
        self.key = key
        self.words = None
        if kind == "FACES":
            kind = FACE_TYPES
        if kind not in ("int", "real", "flag", "str"):
            self.words = kind.split(",")
            kind = "word"
        self.kind = kind
        if count.endswith("+"):
            self.lo, self.hi = int(count[:-1]), None
        elif "-" in count:
            self.lo, self.hi = (int(c) for c in count.split("-"))
        else:
            self.lo = self.hi = int(count)
        self.range = rng

    def count_text(self):
        if self.hi is None:
            return "at least %d value%s" % (self.lo, "" if self.lo == 1 else "s")
        if self.lo != self.hi:
            return "%d or %d values" % (self.lo, self.hi)
        return "exactly %d value%s" % (self.lo, "" if self.lo == 1 else "s")

    def count_ok(self, n):
        return n >= self.lo and (self.hi is None or n <= self.hi)

    def in_range(self, x):
        if math.isnan(x):
            return False
        if not self.range:
            return True
        for alt in self.range.split("|"):
            if alt[0] not in "<>":
                if x == float(alt):
                    return True
                continue
            ok = True
            for cond in alt.split("&"):
                for op, fn in _OPS:
                    if cond.startswith(op):
                        ok = ok and fn(x, float(cond[len(op):]))
                        break
            if ok:
                return True
        return False

    def range_text(self):
        parts = []
        for alt in self.range.split("|"):
            if alt[0] not in "<>":
                parts.append(alt)
                continue
            conds = []
            for cond in alt.split("&"):
                op = cond[:2] if cond[1] == "=" else cond[:1]
                conds.append("%s %s" % (op, cond[len(op):]))
            parts.append(" and ".join(conds))
        return " or ".join(parts)


def _rows():
    out = {}
    for line in SPEC_TEXT.strip().splitlines():
        f = line.split()
        out[f[0]] = Row(f[0], f[1], f[2], f[3] if len(f) > 3 else "")
    return out


SPEC = _rows()


@register
class ValuesChecker(Checker):
    """How many values each key takes, their type, and their range; the keys every deck needs."""

    name = "values"

    def run(self):
        for key in ("mesh.n_cell", "mesh.prob_hi"):
            if not self.deck.has(key):
                self.error(key, "GUST stops: %s is required." % key, "Add %s (three values: x, y, z)." % key)
        for key in self.deck.keys():
            row = SPEC.get(template(key))
            if row is None:
                continue
            vals = self.deck.values(key)
            shown = " ".join(vals)
            if not row.count_ok(len(vals)):
                self.error(key, "GUST stops: %s = %s: the key takes %s." % (key, shown, row.count_text()),
                           "Give %s." % row.count_text())
                continue
            for v in vals:
                self._one(key, row, v, shown)

    def _one(self, key, row, v, shown):
        if row.kind == "str":
            return
        if row.kind == "word":
            if v not in row.words:
                self.error(key, "GUST stops: %s = %s: must be %s." % (key, shown, _or(row.words)),
                           "Use one of the listed names.")
            return
        try:
            x = to_real(v) if row.kind == "real" else float(to_int(v))
        except Unparsed:
            return      # an expression AMReX's parser may accept; not evaluated here
        except ValueError:
            what = "a number" if row.kind == "real" else "an integer"
            self.error(key, "AMReX stops while it reads %s = %s: '%s' is not %s." % (key, shown, v, what),
                       "Give %s." % what)
            return
        if row.kind == "flag":
            if x not in (0.0, 1.0):
                self.error(key, "GUST stops: %s = %s: must be 0 or 1." % (key, shown), "Use 0 or 1.")
            return
        if not row.in_range(x):
            self.error(key, "GUST stops: %s = %s: must be %s." % (key, shown, row.range_text()),
                       "Tests/inputs_reference gives the range and meaning of %s." % key)


def _or(words):
    return ", ".join(words[:-1]) + " or " + words[-1] if len(words) > 1 else words[0]
