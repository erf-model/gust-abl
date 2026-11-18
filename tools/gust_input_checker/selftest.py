"""The checker's own test: small decks that trip one check each, and siblings that do not.

A case is (name, base, deck text, present, absent[, extras]). `base` names one of the decks in BASES; the
text is added to it. `present` and `absent` list findings as "S key" or "S key checker", S being E (ERROR),
W (WARNING), A (ADVICE) or I (INFO). A case passes when every `present` finding is reported, no `absent` one
is, and no finding has the key "internal". `extras` is a dict: files (name -> text, written beside the deck),
dirs (directories made beside the deck), ranks, gpu, made_by_run, overrides.

    python3 tools/check_gust_inputs.py --self-test

tools/test_gust_input_checker.py runs the same table under unittest.
"""
import os
import tempfile

from . import check_deck
from .context import Options
from .findings import Severity

PLAIN = """mesh.n_cell = 16 16 8
mesh.prob_hi = 1600 1600 800
mesh.is_periodic = 1 1 0
"""
GUESS = PLAIN + "init.first_guess = profile\n"
SOLVE = GUESS + "steady.stage = simplec\n"
INFLOW = """mesh.n_cell = 16 16 8
mesh.prob_hi = 1600 1600 800
mesh.is_periodic = 0 1 0
bc.xlo = inflow
bc.xhi = outflow
init.first_guess = profile
steady.stage = simplec
"""
HILL = INFLOW + """terrain.type = agnesi
terrain.hill_height = 100
terrain.hill_half_width = 200
"""
NEST = """mesh.n_cell = 16 16 8
mesh.prob_lo = 400 400 0
mesh.prob_hi = 1200 1200 800
mesh.is_periodic = 0 0 0
bc.xlo = nest
bc.xhi = nest
bc.ylo = nest
bc.yhi = nest
init.first_guess = nest
bc.nest_parent = parent_chk
steady.stage = simplec
"""
UNSTEADY = GUESS + """steady.stage = unsteady
unsteady.dt = 1
unsteady.stop_time = 10
"""
BASES = {"plain": PLAIN, "guess": GUESS, "solve": SOLVE, "inflow": INFLOW, "hill": HILL, "nest": NEST,
         "unsteady": UNSTEADY}

PARENT = {"dirs": ("parent_chk",), "files": {"parent_chk/Header": "gust checkpoint\n"}}
TABLE = "0 0 0 300\n100 5 0 300\n800 10 0 300\n"
NEST_A = "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
JOINT = "GUST_JOINT 1\nfeedback %s\nexchanges %d\nouter 0x0p+0\nnest %s 0x0p+0 chk_a\n"
JOINT2 = "GUST_JOINT 2\nfeedback %s\nexchanges %d\nouter 0x0p+0\nwindow %s\nnest %s 0x0p+0 chk_a\n"

CASES = [
    # --- the bases themselves ---------------------------------------------------------------------------
    ("base_plain", "plain", "", ["W steady.stage"], ["E *"]),
    ("base_guess", "guess", "steady.stage = none\n", [], ["E *", "W *"]),
    ("base_solve", "solve", "", ["I mesh.n_cell"], ["E *", "W *"]),
    ("base_inflow", "inflow", "", ["I bc.inflow_direction"], ["E *", "W *"]),
    ("base_hill", "hill", "", [], ["E *"]),
    ("base_nest", "nest", "", [], ["E *", "W *"], PARENT),
    ("base_unsteady", "unsteady", "", [], ["E *"]),
    # --- keys -------------------------------------------------------------------------------------------
    ("key_unknown", "solve", "steady.max_iteration = 10\n", ["E steady.max_iteration keys"], []),
    ("key_unknown_namespace", "solve", "stedy.max_iterations = 10\n", ["E stedy.max_iterations keys"], []),
    ("key_no_namespace", "solve", "max_iterations = 10\n", ["E max_iterations keys"], []),
    ("key_refavg", "solve", "refavg.file = a\n", ["E refavg.file keys"], []),
    ("key_amrex_passes", "solve", "amrex.fpe_trap_invalid = 1\n", [], ["E amrex.fpe_trap_invalid"]),
    ("key_repeated", "solve", "steady.max_iterations = 10\nsteady.max_iterations = 20\n",
     ["E steady.max_iterations keys"], []),
    ("key_once", "solve", "steady.max_iterations = 10\n", [], ["E steady.max_iterations"]),
    ("key_override_repeated", "solve", "steady.max_iterations = 10\n", ["E steady.max_iterations keys"], [],
     {"overrides": ("steady.max_iterations=20",)}),
    ("key_bc_not_read", "plain", "steady.stage = none\nbc.surface_z0 = 0.1\n", ["E bc.surface_z0 keys"], []),
    ("key_bc_read", "guess", "steady.stage = none\nbc.surface_z0 = 0.1\n", [], ["E bc.surface_z0"]),
    # --- values -----------------------------------------------------------------------------------------
    ("value_choice", "solve", "steady.sweep = jacobi\n", ["E steady.sweep values"], []),
    ("value_choice_ok", "solve", "steady.sweep = zjacobi\n", [], ["E steady.sweep"]),
    ("value_not_a_number", "solve", "steady.cfl_max = fast\n", ["E steady.cfl_max values"], []),
    ("value_range", "solve", "bc.surface_z0 = 0\n", ["E bc.surface_z0 values"], []),
    ("value_range_ok", "solve", "bc.surface_z0 = 0.03\n", [], ["E bc.surface_z0"]),
    ("value_count", "solve", "mesh.max_grid_size = 8 8 8\n", ["E mesh.max_grid_size values"], []),
    ("value_flag", "solve", "io.write_plotfile = 2\n", ["E io.write_plotfile values"], []),
    ("value_tke_floor_small", "solve", "turb.tke_floor = 1e-20\n", ["E turb.tke_floor values"], []),
    ("value_tke_floor_zero", "solve", "turb.tke_floor = 0\n", [], ["E turb.tke_floor"]),
    # --- mesh -------------------------------------------------------------------------------------------
    ("mesh_corner_order", "plain", "mesh.prob_lo = 2000 0 0\n", ["E mesh.prob_hi mesh"], []),
    ("mesh_z_origin", "plain", "mesh.prob_lo = 0 0 10\n", ["E mesh.prob_lo mesh"], []),
    ("mesh_corner_ok", "plain", "mesh.prob_lo = -100 -100 0\n", [], ["E mesh.prob_lo", "E mesh.prob_hi"]),
    ("mesh_not_finite", "plain", "mesh.prob_lo = -inf 0 0\n", ["E mesh.prob_lo mesh"], []),
    ("mesh_z_periodic", "guess", "mesh.is_periodic = 1 1 1\n", ["E mesh.is_periodic mesh"], [],
     {"overrides": ()}),
    ("mesh_ratio_without_dz0", "plain", "mesh.stretch_ratio = 1.1\n", ["E mesh.stretch_ratio mesh"], []),
    ("mesh_ratio_with_dz0", "plain", "mesh.stretch_ratio = 1.1\nmesh.dz0 = 10\n",
     ["I mesh.stretch_ratio"], ["E mesh.stretch_ratio"]),
    ("mesh_dz0_too_thick", "plain", "mesh.dz0 = 200\n", ["E mesh.dz0 mesh"], []),
    ("mesh_dz0_ok", "plain", "mesh.dz0 = 20\n", ["I mesh.dz0"], ["E mesh.dz0"]),
    ("mesh_strong_stretch", "plain", "mesh.dz0 = 1\n", ["A mesh.dz0 mesh"], []),
    # GUST's Newton for the ratio (from 1.03) finds 4.85 for 8 layers from 0.01 m to 800 m, and fails for 0.001 m (6.81)
    # FILE = reads the named file in place, from the run's directory, as AMReX does
    ("file_include", "plain", "FILE = mesh.inc\n", [], ["E *"],
     {"drop": "mesh.n_cell mesh.prob_hi", "files": {"mesh.inc": "mesh.n_cell = 16 16 8\nmesh.prob_hi = 1600 1600 800\n"}}),
    ("file_include_checked", "plain", "FILE = mesh.inc\n", ["E mesh.dz0 mesh"], [], {"files": {"mesh.inc": "mesh.dz0 = 200\n"}}),
    ("file_include_missing", "plain", "FILE = nowhere.inc\n", ["E deck deck"], []),
    ("file_include_itself", "plain", "FILE = inputs\n", ["E deck deck"], []),
    ("mesh_newton_converges", "plain", "mesh.dz0 = 0.01\n", ["A mesh.dz0 mesh"], ["E mesh.dz0"]),
    ("mesh_newton_fails", "plain", "mesh.dz0 = 0.001\n", ["E mesh.dz0 mesh"], []),
    ("mesh_newton_ratio_given", "plain", "mesh.dz0 = 0.001\nmesh.stretch_ratio = 6.81\n", [], ["E mesh.dz0"]),
    # --- hierarchy and sequencing -----------------------------------------------------------------------
    ("levels_not_divisible", "solve", "mesh.n_cell = 18 16 8\nmesh.hierarchy_levels = 3\n",
     ["E mesh.n_cell hierarchy"], [], {"dedupe": "mesh.n_cell"}),
    ("levels_ok", "solve", "mesh.hierarchy_levels = 3\nmesh.max_grid_size = 8\n",
     ["I mesh.hierarchy_levels"], ["E mesh.n_cell", "E mesh.max_grid_size"]),
    ("levels_box_size", "solve", "mesh.hierarchy_levels = 3\nmesh.max_grid_size = 6\n",
     ["E mesh.max_grid_size hierarchy"], []),
    ("levels_one_cell", "solve", "mesh.hierarchy_levels = 5\nmesh.max_grid_size = 16\n",
     ["E mesh.hierarchy_levels hierarchy"], []),
    ("sequencing_needs_levels", "solve", "steady.sequencing = 1\n", ["E steady.sequencing hierarchy"], []),
    ("sequencing_ok", "solve", "steady.sequencing = 1\nmesh.hierarchy_levels = 2\nmesh.max_grid_size = 8\n",
     [], ["E steady.sequencing"]),
    ("sequencing_drop_above_total", "solve",
     "steady.sequencing = 1\nmesh.hierarchy_levels = 2\nmesh.max_grid_size = 8\n"
     "steady.sequencing_drop_orders = 8\nsteady.residual_drop_orders = 4\n",
     ["E steady.sequencing_drop_orders hierarchy"], []),
    ("sequencing_key_alone", "solve", "steady.sequencing_drop_orders = 2\n",
     ["E steady.sequencing_drop_orders hierarchy"], []),
    ("levels_without_sequencing", "solve", "mesh.hierarchy_levels = 2\nmesh.max_grid_size = 8\n",
     ["A steady.sequencing hierarchy"], []),
    ("hierarchy_passes_alone", "solve", "terrain.hierarchy_passes = 2\n",
     ["E terrain.hierarchy_passes hierarchy"], []),
    # --- terrain ----------------------------------------------------------------------------------------
    ("terrain_file_missing_key", "solve", "terrain.type = xyz\n", ["E terrain.file terrain"], []),
    ("terrain_file_on_flat", "solve", "terrain.file = hill.xyz\n", ["E terrain.file terrain"], []),
    ("terrain_agnesi_needs_height", "solve", "terrain.type = agnesi\nterrain.hill_half_width = 200\n",
     ["E terrain.hill_height terrain"], ["E terrain.hill_half_width"]),
    ("terrain_hill_key_on_flat", "solve", "terrain.hill_height = 100\n", ["E terrain.hill_height terrain"], []),
    # --- layout -----------------------------------------------------------------------------------------
    ("layout_fewer_boxes_than_ranks", "solve", "", ["A mesh.max_grid_size layout"], [], {"ranks": 4}),
    ("layout_boxes_match_ranks", "solve", "mesh.max_grid_size = 8\nsteady.sweep = zjacobi\n", ["I mesh.max_grid_size"],
     ["A mesh.max_grid_size"], {"ranks": 4}),
    ("layout_gpu_default_sweep", "solve", "", ["A steady.sweep layout"], [], {"gpu": True}),
    ("layout_cpu_default_sweep", "solve", "", [], ["A steady.sweep"]),
    # --- faces ------------------------------------------------------------------------------------------
    ("face_periodic_on_open", "inflow", "bc.xlo = periodic\n", ["E bc.xlo faces"], [], {"dedupe": "bc.xlo"}),
    ("face_open_on_periodic", "solve", "bc.xlo = inflow\n", ["E bc.xlo faces"], []),
    ("face_ground_type_on_side", "inflow", "bc.xhi = surface_layer\n", ["E bc.xhi faces"], [],
     {"dedupe": "bc.xhi"}),
    ("face_nest_on_top", "solve", "bc.zhi = nest\n", ["E bc.zhi faces"], []),
    ("face_ground_inflow", "solve", "bc.zlo = inflow\n", ["E bc.zlo faces"], []),
    ("face_ground_ok", "solve", "bc.zlo = no_slip_wall\n", [], ["E bc.zlo"]),
    ("face_top_geostrophic_alone", "solve", "bc.zhi = geostrophic\n", ["E bc.zhi faces"], []),
    ("face_wall_velocity_alone", "solve", "bc.zhi_wall_velocity = 1 0 0\n", ["E bc.zhi_wall_velocity faces"], []),
    ("face_wall_velocity_normal", "solve", "bc.zhi = no_slip_wall\nbc.zhi_wall_velocity = 1 0 1\n",
     ["E bc.zhi_wall_velocity faces"], []),
    ("face_wall_velocity_ok", "solve", "bc.zhi = no_slip_wall\nbc.zhi_wall_velocity = 1 0 0\n",
     [], ["E bc.zhi_wall_velocity"]),
    ("face_outflow_transport_alone", "solve", "bc.outflow_transport = zero_gradient\n",
     ["E bc.outflow_transport faces"], []),
    # --- inflow -----------------------------------------------------------------------------------------
    ("inflow_file_missing_key", "inflow", "bc.inflow_profile = file\n", ["E bc.inflow_file inflow"], []),
    ("inflow_z0_differs", "inflow", "bc.surface_z0 = 0.5\n", ["W bc.inflow_z0 inflow"], []),
    ("inflow_z0_equal", "inflow", "bc.surface_z0 = 0.5\nbc.inflow_z0 = 0.5\n", [], ["W bc.inflow_z0"]),
    ("inflow_z0_periodic_column", "solve", "bc.surface_z0 = 0.5\n", [], ["W bc.inflow_z0"]),
    ("inflow_z0_periodic_given", "solve", "bc.surface_z0 = 0.5\nbc.inflow_z0 = 0.03\n", ["W bc.inflow_z0 inflow"], []),
    ("inflow_z0_nest", "nest", "bc.surface_z0 = 0.5\n", [], ["W bc.inflow_z0"], PARENT),
    ("inflow_wind_along_faces", "inflow", "bc.inflow_direction = 90\n", ["W bc.xlo inflow"], []),
    ("inflow_wind_leaves", "inflow", "bc.inflow_direction = 180\n", ["W bc.xlo inflow", "W bc.xhi inflow"], []),
    ("inflow_wind_enters", "inflow", "bc.inflow_direction = 20\n", [], ["W bc.xlo", "W bc.xhi"]),
    # --- sounding ---------------------------------------------------------------------------------------
    ("sounding_lengths", "solve", "bc.theta_profile_z = 0 400 800\nbc.theta_profile = 300 302\n",
     ["E bc.theta_profile sounding"], []),
    ("sounding_order", "solve", "bc.theta_profile_z = 0 800 400\nbc.theta_profile = 300 302 304\n",
     ["E bc.theta_profile_z sounding"], []),
    ("sounding_ok", "solve", "bc.theta_profile_z = 0 400 800\nbc.theta_profile = 300 302 304\n",
     [], ["E bc.theta_profile", "E bc.theta_profile_z"]),
    ("sounding_terrain_no_layer", "hill", "bc.theta_profile_z = 0 800\nbc.theta_profile = 300 304\n",
     ["W bc.rayleigh_zdamp sounding"], []),
    ("sounding_terrain_layer", "hill", "bc.theta_profile_z = 0 800\nbc.theta_profile = 300 304\n"
     "bc.rayleigh_zdamp = 200\nbc.rayleigh_dampcoef = 0.1\n", ["I bc.rayleigh_zdamp"], ["W bc.rayleigh_zdamp"]),
    # --- damping ----------------------------------------------------------------------------------------
    ("rayleigh_too_deep", "inflow", "bc.rayleigh_zdamp = 900\nbc.rayleigh_dampcoef = 0.1\n",
     ["E bc.rayleigh_zdamp damping"], []),
    ("rayleigh_no_rate", "inflow", "bc.rayleigh_zdamp = 200\nbc.rayleigh_dampcoef = 0\n",
     ["E bc.rayleigh_zdamp damping"], []),
    ("rayleigh_rate_alone", "inflow", "bc.rayleigh_dampcoef = 0.1\n", ["E bc.rayleigh_dampcoef damping"], []),
    ("sponge_no_open_face", "solve", "bc.sponge_width = 100\nbc.sponge_strength = 0.1\n",
     ["E bc.sponge_width damping"], []),
    ("sponge_ok", "inflow", "bc.sponge_width = 100\nbc.sponge_strength = 0.1\n", [], ["E bc.sponge_width"]),
    ("sponge_strength_alone", "inflow", "bc.sponge_strength = 0.1\n", ["E bc.sponge_strength damping"], []),
    # --- surface ----------------------------------------------------------------------------------------
    ("surface_key_on_wall", "solve", "bc.zlo = no_slip_wall\nbc.surface_z0 = 0.1\n",
     ["E bc.surface_z0 surface"], []),
    ("surface_two_roughnesses", "solve", "bc.surface_z0 = 0.1\nbc.surface_z0_file = z0.xyz\n",
     ["E bc.surface_z0_file surface"], []),
    ("surface_flux_neutral", "solve", "bc.surface_heat_flux = 0.02\n", ["E bc.surface_heat_flux surface"], []),
    ("surface_flux_ok", "solve", "bc.surface_heat_flux = 0.02\nbc.surface_stability = surface_flux\n",
     [], ["E bc.surface_heat_flux"]),
    ("surface_qnse_neutral", "solve", "bc.surface_stable_functions = qnse\n",
     ["E bc.surface_stable_functions surface"], []),
    # --- nest -------------------------------------------------------------------------------------------
    ("nest_parent_missing_key", "nest", "", ["E bc.nest_parent nest"], [], {"drop": "bc.nest_parent"}),
    ("nest_parent_alone", "solve", "bc.nest_parent = parent_chk\n", ["E bc.nest_parent nest"], [], PARENT),
    ("nest_face_without_guess", "nest", "init.first_guess = profile\n", ["E bc.xlo nest"], [],
     dict(PARENT, drop="init.first_guess")),
    ("nest_with_sequencing", "nest", "steady.sequencing = 1\nmesh.hierarchy_levels = 2\nmesh.max_grid_size = 8\n",
     ["E steady.sequencing nest"], [], PARENT),
    ("nest_relax_too_wide", "nest", "bc.nest_relax_width = 8\n", ["E bc.nest_relax_width nest"], [], PARENT),
    ("nest_relax_ok", "nest", "bc.nest_relax_width = 3\n", [], ["E bc.nest_relax_width"], PARENT),
    ("nest_ratio_ok", "nest", "bc.nest_ratio = 4\n", [], ["E bc.nest_ratio", "W bc.nest_ratio"], PARENT),
    ("nest_ratio_six", "nest", "bc.nest_ratio = 6\n", ["E bc.nest_ratio nest"], [], PARENT),
    ("nest_ratio_no_parent", "solve", "bc.nest_ratio = 4\n", ["E bc.nest_ratio nest"], []),
    ("nest_blend_too_wide", "nest", "terrain.nest_blend_width = 500\n",
     ["E terrain.nest_blend_width nest"], [], PARENT),
    ("nest_slip_side", "nest", "bc.ylo = slip_wall\n", ["E bc.ylo nest"], [], dict(PARENT, drop="bc.ylo")),
    # --- files ------------------------------------------------------------------------------------------
    ("file_missing", "inflow", "bc.inflow_profile = file\nbc.inflow_file = table.dat\n",
     ["E bc.inflow_file files"], []),
    ("file_present", "inflow", "bc.inflow_profile = file\nbc.inflow_file = table.dat\n",
     [], ["E bc.inflow_file", "I bc.inflow_file"], {"files": {"table.dat": TABLE}}),
    ("file_made_by_run", "inflow", "bc.inflow_profile = file\nbc.inflow_file = column/table.dat\n",
     ["I bc.inflow_file files"], ["E bc.inflow_file"], {"made_by_run": ("column/table.dat",)}),
    ("file_made_by_run_pattern", "inflow", "bc.inflow_profile = file\nbc.inflow_file = column/table.dat\n",
     ["I bc.inflow_file files"], ["E bc.inflow_file"], {"made_by_run": ("column/*",)}),
    ("file_made_by_run_other", "inflow", "bc.inflow_profile = file\nbc.inflow_file = column/table.dat\n",
     ["E bc.inflow_file files"], ["I bc.inflow_file"], {"made_by_run": ("other/table.dat",)}),
    ("file_is_directory", "inflow", "bc.inflow_profile = file\nbc.inflow_file = column\n",
     ["E bc.inflow_file files"], [], {"dirs": ("column",)}),
    ("restart_missing", "plain", "steady.stage = none\nio.restart = chk00010\n", ["E io.restart files"], []),
    ("restart_no_header", "plain", "steady.stage = none\nio.restart = chk00010\n", ["E io.restart files"], [],
     {"dirs": ("chk00010",)}),
    ("restart_present", "plain", "steady.stage = none\nio.restart = old00010\n", [], ["E io.restart"],
     {"files": {"old00010/Header": "x\n"}}),
    ("nest_parent_missing", "nest", "", ["E bc.nest_parent files"], []),
    ("nest_parent_made_by_run", "nest", "", ["I bc.nest_parent files"], ["E bc.nest_parent"],
     {"made_by_run": ("parent_chk",)}),
    ("nest_parent_made_by_run_relative", "nest", "bc.nest_parent = ../parent/chk00000\n",
     ["I bc.nest_parent files"], ["E bc.nest_parent"],
     {"drop": "bc.nest_parent", "made_by_run": ("../parent/chk00000",)}),
    # --- output -----------------------------------------------------------------------------------------
    ("output_turb_plot_no_solve", "guess", "steady.stage = none\nio.turb_plot_file = turb\n",
     ["E io.turb_plot_file output"], []),
    ("output_turb_plot_solve", "solve", "io.turb_plot_file = turb\n", [], ["E io.turb_plot_file"]),
    ("output_column_file_alone", "solve", "io.column_guess_file = col.dat\n",
     ["E io.column_guess_file output"], []),
    ("output_overwrites_restart", "plain", "steady.stage = none\nio.restart = ./chk00010/\n",
     ["E io.chk_file output"], [], {"files": {"chk00010/Header": "x\n"}}),
    ("output_other_name", "plain", "steady.stage = none\nio.restart = chk00010\nio.chk_file = new\n",
     [], ["E io.chk_file"], {"files": {"chk00010/Header": "x\n"}}),
    ("output_no_checkpoint", "plain", "steady.stage = none\nio.restart = chk00010\nio.write_checkpoint = 0\n",
     [], ["E io.chk_file"], {"files": {"chk00010/Header": "x\n"}}),
    # --- init -------------------------------------------------------------------------------------------
    ("init_solve_without_guess", "plain", "steady.stage = simplec\n", ["E init.first_guess init"], []),
    ("init_solve_with_guess", "solve", "", [], ["E init.first_guess"]),
    # --- forcing ----------------------------------------------------------------------------------------
    ("forcing_geostrophic_without_coriolis", "solve", "forcing.geostrophic = 1\n",
     ["E forcing.geostrophic forcing"], []),
    ("forcing_coriolis_info", "solve", "forcing.coriolis = 1\nforcing.latitude = 45\n",
     ["I forcing.coriolis"], ["E forcing.coriolis"]),
    # --- turbulence -------------------------------------------------------------------------------------
    ("turb_lee_setting_advice", "hill", "", ["A * turbulence"], []),
    ("turb_no_advice_flat", "inflow", "", [], ["A * turbulence"]),
    # --- steady -----------------------------------------------------------------------------------------
    ("steady_key_without_stage", "guess", "steady.stage = none\nsteady.cfl = 5\n", ["E steady.cfl steady"], []),
    ("steady_key_with_stage", "solve", "steady.cfl_max = 5\n", [], ["E steady.cfl_max"]),
    ("steady_stage_not_given", "plain", "", ["W steady.stage steady"], []),
    ("steady_newton_key_alone", "solve", "steady.newton_gmres_restart = 30\n",
     ["E steady.newton_gmres_restart steady"], []),
    ("steady_many_orders", "solve", "steady.residual_drop_orders = 16\n",
     ["W steady.residual_drop_orders steady"], []),
    ("steady_usual_orders", "solve", "steady.residual_drop_orders = 4\n", [], ["W steady.residual_drop_orders"]),
    # --- unsteady ---------------------------------------------------------------------------------------
    ("unsteady_key_without_stage", "solve", "unsteady.dt = 1\n", ["E unsteady.dt unsteady"], []),
    ("unsteady_needs_dt", "unsteady", "", ["E unsteady.dt unsteady"], [], {"drop": "unsteady.dt"}),
    # --- pressure ---------------------------------------------------------------------------------------
    ("pressure_thin_cells", "hill", "mesh.dz0 = 1\n", ["A * pressure"], []),
    ("pressure_usual_cells", "hill", "", [], ["A * pressure", "W * pressure"]),
    # --- nests (nest.*) ----------------------------------------------------------------------------------
    ("nests_ok", "inflow", "nest.names = a b\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.b.lo = 1200 0\nnest.b.hi = 1400 1600\n", ["I nest.a.lo nests"], ["E *"]),
    ("nests_misaligned", "inflow", "nest.names = a\nnest.a.lo = 450 300\nnest.a.hi = 1200 900\n",
     ["E nest.a.lo nests"], []),
    ("nests_margin", "inflow", "nest.names = a\nnest.a.lo = 100 300\nnest.a.hi = 1200 900\n",
     ["E nest.a.lo nests"], []),
    ("nests_hi_not_above", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 400 900\n",
     ["E nest.a.hi nests"], []),
    ("nests_overlap", "inflow", "nest.names = a b\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.b.lo = 1100 300\nnest.b.hi = 1400 600\n", ["E nest.b.lo nests"], ["E nest.a.lo"]),
    ("nests_missing_corner", "inflow", "nest.names = a\nnest.a.lo = 400 300\n", ["E nest.a.hi nests"], []),
    ("nests_not_named", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.c.lo = 400 300\n", ["E nest.c.lo nests"], ["E nest.a.lo"]),
    ("nests_bad_name", "inflow", "nest.names = 9x\n", ["E nest.names nests"], []),
    ("nests_name_twice", "inflow", "nest.names = a a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n",
     ["E nest.names nests"], []),
    ("nests_unknown_key", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\nnest.a.ratio = 3\n",
     ["E nest.a.ratio keys"], []),
    ("nests_three_values", "inflow", "nest.names = a\nnest.a.lo = 400 300 0\nnest.a.hi = 1200 900\n",
     ["E nest.a.lo values"], []),
    ("nests_no_stage", "guess", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n",
     ["E nest.names nests"], []),
    ("nests_no_checkpoint", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "io.write_checkpoint = 0\n", ["E io.write_checkpoint nests"], []),
    ("nests_whole_domain", "solve", "nest.names = a\nnest.a.lo = 0 0\nnest.a.hi = 1600 1600\n",
     ["E nest.a.lo nests"], []),
    ("nests_narrow_box", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 500 900\nsteady.sweep = lusgs\n"
     "steady.line_ghosts = zero\nsteady.line_overlap = 4\n", ["E steady.line_overlap nests"], []),
    ("nests_wide_box", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 600 900\nsteady.sweep = lusgs\n"
     "steady.line_ghosts = zero\nsteady.line_overlap = 4\n", [], ["E steady.line_overlap"]),
    ("nests_hill_stf", "hill", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\nsteady.stage = simplec\n",
     ["E terrain.smoothing nests"], []),
    ("nests_hill_btf", "hill", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\nsteady.stage = simplec\n"
     "terrain.smoothing = 0\n", [], ["E terrain.smoothing"]),
    ("nests_feedback_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\n", [], ["E *"]),
    ("nests_feedback_value", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = average\n", ["E nest.feedback nests"], []),
    ("nests_inject_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\n", [], ["E *"]),
    ("nests_inject_factor", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nnest.feedback_factor = 4\n", ["E nest.feedback_factor nests"], []),
    ("nests_inject_factor_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nnest.feedback_factor = 0.25\n", [], ["E *"]),
    ("nests_inject_factor_one", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nnest.feedback_factor = 1\n", [], ["E *"]),
    ("nests_inject_factor_zero", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nnest.feedback_factor = 0\n", ["E nest.feedback_factor values"],
     ["E nest.feedback_factor nests"]),
    ("nests_inject_mlmg", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nsolver.poisson_backend = mlmg\n", ["E solver.poisson_backend nests"], []),
    ("nests_relax_mlmg_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nsolver.poisson_backend = mlmg\n", [], ["E solver.poisson_backend nests"]),
    ("nests_inject_periodic", "inflow", "nest.names = a\nnest.a.lo = 400 0\nnest.a.hi = 1200 1600\n"
     "nest.feedback = inject\n", ["E nest.feedback nests"], []),
    ("nests_relax_periodic_ok", "inflow", "nest.names = a\nnest.a.lo = 400 0\nnest.a.hi = 1200 1600\n"
     "nest.feedback = relax\n", [], ["E nest.feedback"]),
    ("nests_inject_feedback_flux", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = inject\nnest.feedback_flux = nest\n", ["E nest.feedback_flux nests"], []),
    ("nests_inject_narrow", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 800 900\n"
     "nest.feedback = inject\n", ["E nest.feedback nests"], []),
    ("nests_feedback_no_nests", "inflow", "nest.feedback = defect\n", ["E nest.feedback nests"], []),
    ("nests_feedback_named_feedback", "inflow", "nest.names = feedback\n", ["E nest.names nests"], []),
    ("nests_feedback_newton", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nsteady.stage = simplec_newton\n", ["E nest.feedback nests"], []),
    ("nests_feedback_damping", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nbc.rayleigh_zdamp = 400\nbc.rayleigh_dampcoef = 0.2\nbc.rayleigh_fields = u w\n", ["E nest.feedback nests"], []),
    ("nests_feedback_damping_w", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nbc.rayleigh_zdamp = 400\nbc.rayleigh_dampcoef = 0.2\nbc.rayleigh_fields = w theta\n", [], ["E *"]),
    ("nests_feedback_narrow", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 800 900\n"
     "nest.feedback = defect\n", ["E nest.feedback nests"], []),
    ("nests_feedback_five_cells", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 900 900\n"
     "nest.feedback = defect\n", [], ["E *"]),
    ("nests_feedback_sponge", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nbc.sponge_width = 200\n", ["W nest.feedback nests"], ["E nest.feedback"]),
    ("nests_feedback_odd_boxes", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nmesh.max_grid_size = 7\n", ["E nest.feedback nests"], []),
    ("nests_relax_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nnest.feedback_factor = 4\n", [], ["E *"]),
    ("nests_relax_newton", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nsteady.stage = simplec_newton\n", ["E nest.feedback nests"], []),
    ("nests_relax_odd_boxes", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nmesh.max_grid_size = 7\n", ["E nest.feedback nests"], []),
    ("nests_relax_factor_zero", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nnest.feedback_factor = 0\n", ["E nest.feedback_factor values"], []),
    ("nests_relax_factor_with_defect", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nnest.feedback_factor = 2\n", ["E nest.feedback_factor nests"], []),
    ("nests_named_feedback_factor", "inflow", "nest.names = feedback_factor\n", ["E nest.names nests"], []),
    ("nests_named_feedback_flux", "inflow", "nest.names = feedback_flux\n", ["E nest.names nests"], []),
    # nest.damping: its values, its need of nest.names and of a layer on u, v or tke; outer lifts the joint run's
    # limit on bc.rayleigh_fields.
    ("nests_named_damping", "inflow", "nest.names = damping\n", ["E nest.names nests"], []),
    ("nests_damping_value", "inflow", NEST_A + "nest.damping = all\n", ["E nest.damping values"], []),
    ("nests_damping_no_names", "inflow", "nest.damping = outer\n", ["E nest.damping nests"], []),
    ("nests_damping_no_effect", "inflow", NEST_A + "nest.damping = outer\nbc.rayleigh_zdamp = 300\n"
     "bc.rayleigh_dampcoef = 0.1\nbc.rayleigh_fields = w\n", ["E nest.damping nests"], []),
    ("nests_damping_outer_joint", "inflow", NEST_A + "nest.feedback = defect\nnest.damping = outer\n"
     "bc.rayleigh_zdamp = 300\nbc.rayleigh_dampcoef = 0.1\nbc.rayleigh_fields = u v w\n", [],
     ["E nest.feedback", "E nest.damping"]),
    ("nests_damping_default_joint", "inflow", NEST_A + "nest.feedback = defect\n"
     "bc.rayleigh_zdamp = 300\nbc.rayleigh_dampcoef = 0.1\nbc.rayleigh_fields = u v w\n",
     ["E nest.feedback nests"], []),
    ("nests_feedback_flux_ok", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nnest.feedback_flux = nest\n", [], ["E *"]),
    ("nests_feedback_flux_value", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = defect\nnest.feedback_flux = fine\n", ["E nest.feedback_flux values"], []),
    ("nests_feedback_flux_with_relax", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = relax\nnest.feedback_flux = nest\n", ["E nest.feedback_flux nests"], []),
    # A joint restart (io.restart with nest.feedback): the checkpoint's joint record is read as the run reads it; a
    # checkpoint without one restarts the outer domain alone. (init.first_guess of the base is reported on io.restart.)
    ("nests_joint_restart_ok", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n"
     "steady.max_iterations = 10\n", [], ["E nest.feedback", "E nest.names", "E steady.max_iterations", "E io.restart nests"],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT % ("relax", 4, "a"), "chk_a/Header": "x\n"}}),
    # Written on another system: \r\n line ends and trailing spaces are not part of the nest checkpoint's path.
    ("nests_joint_restart_crlf", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n", [],
     ["E io.restart nests", "E nest.names"],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": (JOINT % ("relax", 4, "a")).replace("chk_a\n", "chk_a \t\n").replace("\n", "\r\n"),
                "chk_a/Header": "x\n"}}),
    # The record names chk_a, which is not in the directory the run starts in.
    ("nests_joint_restart_missing", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E io.restart nests"], ["E nest.names"],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT % ("relax", 4, "a")}}),
    ("nests_joint_restart_feedback", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E io.restart nests"], [],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT % ("none", 4, "a"), "chk_a/Header": "x\n"}}),
    ("nests_joint_restart_no_record", "inflow", NEST_A + "nest.feedback = relax\nio.restart = chk00010\n",
     [], ["E nest.feedback", "E nest.names"], {"files": {"chk00010/Header": "x\n"}}),
    ("nests_joint_restart_names", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E nest.names nests"], [], {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT % ("defect", 4, "b")}}),
    ("nests_joint_restart_iterations", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n"
     "steady.max_iterations = 4\n", ["E steady.max_iterations nests"], ["E nest.names"],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT % ("defect", 4, "a")}}),
    ("nests_joint_restart_version", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E io.restart nests"], [],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": "GUST_JOINT 3\nfeedback defect\nexchanges 4\n"
                                                              "outer 0x0p+0\nwindow 0\nnest a 0x0p+0 chk_a\n"}}),
    # Version 2 (the joint stop test's window, Source/JointRecord.H): a good record reads, a missing or bad window
    # line stops GUST.
    ("nests_joint_restart_v2", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     [], ["E io.restart nests", "E nest.names"],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT2 % ("defect", 4, "2 0x1p-30 0x0p+0", "a"),
                "chk_a/Header": "x\n"}}),
    ("nests_joint_restart_v2_no_window", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E io.restart nests"], [],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": (JOINT % ("defect", 4, "a"))
                                                              .replace("GUST_JOINT 1", "GUST_JOINT 2"),
                "chk_a/Header": "x\n"}}),
    ("nests_joint_restart_v2_bad_window", "inflow", NEST_A + "nest.feedback = defect\nio.restart = chk00010\n",
     ["E io.restart nests"], [],
     {"files": {"chk00010/Header": "x\n", "chk00010/Joint": JOINT2 % ("defect", 4, "2 0x1p-30", "a")}}),
    ("nests_feedback_none_odd_boxes", "inflow", "nest.names = a\nnest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
     "nest.feedback = none\nmesh.max_grid_size = 7\n", [], ["E nest.feedback"]),
    ("nests_periodic_strip", "solve", "nest.names = a\nnest.a.lo = 0 0\nnest.a.hi = 1600 1500\n",
     ["I nest.a.lo nests"], ["E *"]),
    # --- canopy and buildings ---------------------------------------------------------------------------
    ("canopy_key_without_type", "solve", "canopy.height = 15\n", ["E canopy.height canopy"], []),
    ("buildings_need_zjacobi", "solve", "buildings.type = immersed_forcing\nbuildings.file = b.txt\n",
     ["E steady.sweep buildings"], [], {"files": {"b.txt": "0 0 10 10 5\n"}}),
    ("buildings_zjacobi", "solve", "buildings.type = immersed_forcing\nbuildings.file = b.txt\n"
     "steady.sweep = zjacobi\n", [], ["E steady.sweep"], {"files": {"b.txt": "0 0 10 10 5\n"}}),
    # --- abl --------------------------------------------------------------------------------------------
    # --- turbulence --------------------------------------------------------------------------------------
    ("nut_relax_end_without_relax", "solve", "turb.nut_relax_end = 0.01\n", ["E turb.nut_relax_end turbulence"], []),
    ("nut_relax_end_range", "solve", "turb.nut_relax = 0.5\nturb.nut_relax_end = 1\n", ["E turb.nut_relax_end"], []),
    ("nut_relax_end_unsteady", "unsteady", "turb.nut_relax = 0.5\nturb.nut_relax_end = 0.01\n",
     ["E turb.nut_relax_end turbulence"], []),
    ("nut_relax_end_ok", "solve", "turb.nut_relax = 0.5\nturb.nut_relax_end = 0.01\n", [],
     ["E turb.nut_relax_end"]),
    # --- convection ---------------------------------------------------------------------------------------
    ("convection_blend_no_height", "solve", "turb.realizability = 0.309\nconvection.hybrid_sigma_max = 0.1\n",
     ["E convection.hybrid_height convection"], []),
    ("convection_blend_no_cap", "solve", "convection.hybrid_sigma_max = 0.1\nconvection.hybrid_height = 250\n",
     ["E convection.hybrid_sigma_max convection"], []),
    ("convection_blend_ok", "solve",
     "turb.realizability = 0.309\nconvection.hybrid_sigma_max = 0.1\nconvection.hybrid_height = 250\n",
     [], ["E convection.hybrid_height", "E convection.hybrid_sigma_max"]),
    ("abl_key_without_stage", "guess", "steady.stage = none\nabl.drive = geostrophic_column\n", ["E abl.drive abl"], []),
]

_SEVERITY = {"E": Severity.ERROR, "W": Severity.WARNING, "A": Severity.ADVICE, "I": Severity.INFO}


def _wanted(text):
    parts = text.split()
    return _SEVERITY[parts[0]], parts[1], parts[2] if len(parts) > 2 else None


def _matches(finding, wanted):
    severity, key, checker = wanted
    return (finding.severity == severity and (key == "*" or finding.key == key)
            and (checker is None or finding.checker == checker))


def _deck_text(base, text, extras):
    """The base deck with `text` added; a key the text gives again replaces the base's line."""
    again = set(line.split("=")[0].strip() for line in text.splitlines() if "=" in line)
    again.update(extras.get("drop", "").split())
    kept = [line for line in BASES[base].splitlines() if line.split("=")[0].strip() not in again]
    repeated = set()
    lines = []
    for line in text.splitlines():
        key = line.split("=")[0].strip()
        if key in extras.get("drop", "").split() and key not in repeated and text.count(key + " =") == 0:
            continue
        lines.append(line)
    return "\n".join(kept + lines) + "\n"


def run_case(case):
    """Run one case in a fresh directory. Returns (problems, findings): no problems means the case passed."""
    name, base, text, present, absent = case[:5]
    extras = case[5] if len(case) > 5 else {}
    with tempfile.TemporaryDirectory() as folder:
        for sub in extras.get("dirs", ()):
            os.makedirs(os.path.join(folder, sub), exist_ok=True)
        for rel, body in extras.get("files", {}).items():
            path = os.path.join(folder, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as handle:
                handle.write(body)
        deck = os.path.join(folder, "inputs")
        with open(deck, "w") as handle:
            handle.write(_deck_text(base, text, extras))
        options = Options(ranks=extras.get("ranks"), gpu=extras.get("gpu", False), cwd=folder,
                          made_by_run=extras.get("made_by_run", ()))
        report = check_deck(deck, overrides=extras.get("overrides", ()), options=options)
    problems = []
    for item in present:
        if not any(_matches(f, _wanted(item)) for f in report.findings):
            problems.append("%s: expected '%s', not reported" % (name, item))
    for item in absent:
        hits = [f for f in report.findings if _matches(f, _wanted(item))]
        if hits:
            problems.append("%s: '%s' not expected, got: %s" % (name, item, hits[0].message))
    for f in report.findings:
        if f.key == "internal":
            problems.append("%s: internal error in checker %s: %s" % (name, f.checker, f.message))
    return problems, report.findings


def checkers_covered():
    """The checker names that some case expects a finding from (by its third field)."""
    names = set()
    for case in CASES:
        for item in case[3]:
            checker = _wanted(item)[2]
            if checker:
                names.add(checker)
    return names


def main(verbose=False, out=None):
    """Run every case; print the failures and a one-line summary. Returns 0 when all pass, else 1."""
    import sys
    out = out or sys.stdout
    failed = 0
    for case in CASES:
        problems, findings = run_case(case)
        if problems:
            failed += 1
            for line in problems:
                out.write("FAIL " + line + "\n")
            for f in findings:
                if verbose or f.severity != Severity.INFO:
                    out.write("       %s %s %s | %s\n" % (f.severity.label, f.checker, f.key, f.message[:100]))
    if failed:
        out.write("self-test: %d of %d cases failed\n" % (failed, len(CASES)))
        return 1
    out.write("self-test: %d cases passed\n" % len(CASES))
    return 0
