#include "Domain.H"

#include "AblSetup.H"
#include "BoundaryConditions.H"
#include "BuildingLoads.H"
#include "Buildings.H"
#include "BuildingForcing.H"
#include "BuildingTheta.H"
#include "BuildingTke.H"
#include "Canopy.H"
#include "CanopyDrag.H"
#include "Closure.H"
#include "ColumnGuess.H"
#include "Convection.H"
#include "Driver.H"
#include "Forcing.H"
#include "StateIO.H"
#include "Initializer.H"
#include "NestBuildings.H"
#include "NestDefect.H"
#include "NestFaceFlux.H"
#include "NestRelax.H"
#include "NestInject.H"
#include "NestFaces.H"
#include "NestParent.H"
#include "Nests.H"
#include "ImmersedTerrain.H"
#include "Inputs.H"
#include "Mesh.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "Sequencing.H"
#include "TimeDriver.H"
#include "TimeTerm.H"
#include "Simplec.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TurbulenceOutput.H"
#include "UnsteadyOutput.H"
#include "WallDistance.H"
#include "ZLevels.H"

#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <cstdio>
#include <memory>
#include <limits>
#include <string>

namespace gust {

namespace {

std::string fmt (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", x);
    return b;
}

} // namespace

std::string run_domain (const DomainRole& role)
{
    // A nest takes the outer deck and replaces what Nests.H derives.
    const NestSpec* const nest = role.nest;
    const std::string role_name = nest != nullptr ? "nest." + nest->name : std::string("outer domain");
    std::string final_chk;
    // The inputs: mesh, terrain, io; with the initializer (init.first_guess =
    // profile) or the SIMPLEC stage (steady.stage = simplec) the boundary
    // conditions (bc.*) and the Poisson controls (solver.poisson_*); with the
    // stage also the closure (turb.*), the forcing (forcing.*), the
    // convection (convection.*), the stage's own (simplec.*) and the driver's
    // (steady.*).
    // The named set-up (AblSetup.H) gives its keys their defaults before any reader.
    const AblSetupOptions abl = read_abl_setup_options();
    if (!abl.errors.empty()) {
        std::string all = "gust: invalid inputs:";
        for (const auto& e : abl.errors) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    apply_abl_setup_defaults(abl);
    Inputs in = parse_inputs();
    if (nest != nullptr) { apply_nest_inputs(*nest, in); }
    // A joint restart (JointRecord.H): a nest restarts from the checkpoint the record names.
    if (nest != nullptr && role.resume != nullptr) { in.io.restart = role.resume->chk; }
    auto err = validate(in);
    InitOptions iopt = read_init_options();
    SteadyOptions sopt = read_steady_options();
    SequencingOptions seq = read_sequencing_options();
    if (nest != nullptr) {
        // From its parent on level 0; the probes are the outer domain's.
        // The outer deck's init.* keys set up the outer domain's first guess: a nest's
        // field comes from the outer domain, so they are not checked again here
        // (init.projection, alpha_h and alpha_v still apply, value-checked as before).
        iopt.first_guess = role.resume != nullptr ? InitOptions::FirstGuess::none : InitOptions::FirstGuess::nest;
        iopt.given.clear();
        iopt.obrien = false;   // the parent gives w: no O'Brien pass (nor its log line)
        seq.enabled = false;
        seq.given.clear();
        sopt.probes_file.clear();
    }
    const bool initialize = (iopt.first_guess != InitOptions::FirstGuess::none);
    const bool solve = (sopt.stage != StageKind::none);
    const bool newton = (sopt.stage == StageKind::simplec_newton);
    const bool unsteady = (sopt.stage == StageKind::unsteady);
    // steady.stage = unsteady (TimeDriver.H): the unsteady.* keys; without it every one has no effect.
    UnsteadyOptions uopt;
    if (unsteady) { uopt = read_unsteady_options(); }
    BoundaryInputs bin;
    PoissonOptions popt;
    SurfaceLayerOptions surf;
    ClosureOptions copt;
    ForcingOptions fopt;
    ConvectionInputs cin;
    SimplecOptions simplec;
    NewtonStageOptions newton_opt;
    bool surface_layer = false;
    if (initialize || solve) {
        bin = read_boundary_inputs();
        if (nest != nullptr) {
            apply_nest_boundary_inputs(role.parent_chk, in.mesh.is_periodic,
                                       read_nest_options().damping == NestDamping::outer, bin);
        }
        // The first guess uses the profile on any geometry; a restart keeps
        // the deck of the run it continues, inflow keys included.
        // A nest's first guess is its parent's field, not the profile.
        bin.profile_needed = (initialize && iopt.first_guess != InitOptions::FirstGuess::nest) || !in.io.restart.empty();
        popt = read_poisson_options();
        surface_layer = (bin.face[amrex::Orientation(2, amrex::Orientation::low)] == FaceType::surface_layer);
        // The surface layer's options are not used by the initial field; they
        // are read so that their range checks run.
        if (surface_layer) { surf = read_surface_layer_options(); }
    }
    if (solve) {
        copt = read_closure_options();
        use_closure_constants(bin, surf, copt.Cmu0, copt.kappa, copt.max_geom_lscale);
        fopt = read_forcing_options();
        // A nest: the forcing and the surface flux the outer run used (OuterOutcome).
        if (role.from_outer != nullptr) { take_outer_forcing(*role.from_outer, fopt, surf); }
        cin = read_convection_options();
        if (cin.hybrid.active() && newton) {
            err.emplace_back("convection.hybrid_sigma_max is not supported with steady.stage = simplec_newton");
        }
        if (cin.hybrid.active() && !(copt.realizability > 0)) {
            err.emplace_back("convection.hybrid_sigma_max needs turb.realizability > 0 (the blend follows the nu_t cap)");
        }
        if (cin.hybrid.active() && unsteady) {
            err.emplace_back("convection.hybrid_sigma_max is not supported with steady.stage = unsteady (the blend "
                             "changes with the inner residual, so the operator would change within and between steps)");
        }
        // isnan first: the ordered comparison traps on NaN (check_closure_options reports it).
        if (unsteady && !std::isnan(copt.nut_relax_end) && copt.nut_relax_end > 0) {
            err.emplace_back("turb.nut_relax_end is not supported with steady.stage = unsteady (the end follows the "
                             "inner residual, so nu_t would change from relaxed to unrelaxed within a physical step)");
        }
        if (cin.hybrid.active() && !in.io.restart.empty()) {
            err.emplace_back("convection.hybrid_sigma_max is not supported with io.restart (the blend is not in the "
                             "checkpoint, so a restart would not continue the run)");
        }
        simplec = read_simplec_options();
        if (newton && simplec.theta_frozen) {
            err.emplace_back("simplec.theta = frozen is not supported with steady.stage = simplec_newton (the Newton stage "
                             "solves theta)");
        }
        if (simplec.theta_frozen && !in.io.restart.empty()) {
            err.emplace_back("simplec.theta = frozen is not supported with io.restart (theta is held at the first guess; restart with "
                             "simplec.theta = solve, or start again from the first guess)");
        }
        if (newton) {
            newton_opt = read_newton_stage_options(sopt, copt);
        } else {
            // The Newton stage's keys, and the switch rules (one stage: nothing to switch to).
            for (const auto& k : sopt.given) {
                if (k == "newton_preconditioner" || k == "newton_fallback_failures" || k == "newton_cfl_max"
                    || k.rfind("switch_", 0) == 0) {
                    err.emplace_back("steady." + k + " has no effect without steady.stage = simplec_newton");
                }
            }
            for (const InputEntry& e : parmparse_entries()) {
                if (e.name.rfind("steady.newton_", 0) == 0 && e.name != "steady.newton_preconditioner"
                    && e.name != "steady.newton_fallback_failures" && e.name != "steady.newton_cfl_max") {
                    err.emplace_back(e.name + " has no effect without steady.stage = simplec_newton");
                }
            }
            if (amrex::ParmParse("turb").contains("frozen_switch")) {
                err.emplace_back("turb.frozen_switch has no effect without steady.stage = simplec_newton");
            }
        }
    } else {
        // Keys of the stage alone have no effect without it.
        for (const InputEntry& e : parmparse_entries()) {
            for (const char* ns : {"forcing.", "convection.", "simplec.", "turb."}) {
                if (e.name.rfind(ns, 0) == 0) { err.emplace_back(e.name + " has no effect without steady.stage = simplec"); }
            }
        }
        if (!in.io.turb_plot_file.empty()) { err.emplace_back("io.turb_plot_file has no effect without steady.stage = simplec"); }
        for (const auto& k : sopt.given) {
            if (k != "stage") { err.emplace_back("steady." + k + " has no effect without steady.stage = simplec"); }
        }
        // The Newton stage's keys that steady.* does not read (Newton.H, Continuation.H).
        for (const InputEntry& e : parmparse_entries()) {
            if (e.name.rfind("steady.newton_", 0) == 0 && e.name != "steady.newton_preconditioner"
                && e.name != "steady.newton_fallback_failures" && e.name != "steady.newton_cfl_max") {
                err.emplace_back(e.name + " has no effect without steady.stage = simplec_newton");
            }
        }
        for (const auto& k : seq.given) { err.emplace_back("steady." + k + " has no effect without steady.stage = simplec"); }
    }
    if (solve) {
        for (auto& e : check_sequencing_options(seq, sopt.criterion, in.mesh)) { err.push_back(std::move(e)); }
    }
    if (unsteady) {
        // The steady stop and its monitors are not the stop of an unsteady run (TimeDriver.H).
        for (const auto& k : sopt.given) {
            if (k == "max_iterations" || k == "residual_drop_orders" || k == "monitor_rel_tol" || k == "monitor_window"
                || k == "confidence" || k == "confidence_height") {
                err.emplace_back("steady." + k + " has no effect with steady.stage = unsteady (each physical step stops at "
                                 "unsteady.inner_drop_orders or unsteady.max_inner)");
            }
        }
        if (!sopt.probes_file.empty()) {
            err.emplace_back("steady.probes_file has no effect with steady.stage = unsteady (unsteady.probes_file names the "
                             "probes of an unsteady run)");
        }
        if (seq.enabled) {
            err.emplace_back("steady.sequencing = 1 is not supported with steady.stage = unsteady yet (the first step "
                             "starts from the first guess or the restart on level 0)");
        }
        if (uopt.checkpoint_interval > 0 && in.io.write_checkpoint != 1) {
            err.emplace_back("unsteady.checkpoint_interval > 0 needs io.write_checkpoint = 1");
        }
        if (uopt.plot_interval > 0 && in.io.write_plotfile != 1) {
            err.emplace_back("unsteady.plot_interval > 0 needs io.write_plotfile = 1");
        }
        if (amrex::ParmParse("io").contains("iteration_log")) {
            err.emplace_back("io.iteration_log has no effect with steady.stage = unsteady (one line per physical step)");
        }
    } else {
        for (const InputEntry& e : parmparse_entries()) {
            if (e.name.rfind("unsteady.", 0) == 0) { err.emplace_back(e.name + " has no effect without steady.stage = unsteady"); }
        }
    }
    // The forest canopy (Canopy.H): read in every run, so that a run without
    // the stage writes the drag field with its plotfile.
    const CanopyOptions canopy_opt = read_canopy_options(in.terrain.idw_neighbors, static_cast<double>(in.terrain.idw_power));
    for (auto& e : check_canopy_options(canopy_opt)) { err.push_back(std::move(e)); }
    // The buildings (Buildings.H, BuildingForcing.H). PR-39: one level and SIMPLEC with z-line sweeps (PR-43 the rest).
    BuildingsOptions buildings_opt = read_buildings_options();
    // The lee options' zone acts with the closure's cap (and the upwind blend, which needs the cap).
    buildings_opt.limiter_zone_acts = buildings_opt.active() && buildings_opt.wall_tke && buildings_opt.limiter_zone > 0
                                      && copt.realizability > 0;
    if (solve && buildings_opt.active() && buildings_opt.limiter_zone > 0 && !(copt.realizability > 0)
        && std::find(buildings_opt.given.begin(), buildings_opt.given.end(), "limiter_zone") != buildings_opt.given.end()) {
        err.emplace_back("buildings.limiter_zone has no effect without turb.realizability > 0 (the zone fades out the "
                         "eddy-viscosity cap and the upwind blend that needs it)");
    }
    for (auto& e : check_buildings_options(buildings_opt)) { err.push_back(std::move(e)); }
    if (!solve && buildings_opt.active() && !buildings_opt.force_file.empty()) {
        err.emplace_back("buildings.force_file has no effect without steady.stage (the force is the residual of a solved state)");
    }
    if (solve && buildings_opt.active()) {
        if (newton) { err.emplace_back("buildings.type = immersed_forcing is not supported with steady.stage = simplec_newton yet"); }
        if (sopt.sweep != SweepKind::zjacobi) { err.emplace_back("buildings.type = immersed_forcing needs z-line sweeps (steady.sweep) for now"); }
        if (in.mesh.hierarchy_levels > 1) { err.emplace_back("buildings.type = immersed_forcing needs mesh.hierarchy_levels = 1 for now"); }
        if (popt.backend == "mlmg") {
            err.emplace_back("solver.poisson_backend = mlmg cannot solve with buildings (no solid cells in its operator); use gmres_mlmg or gmres_fft");
        }
    }
    // The ground immersed on a flat mesh (ImmersedTerrain.H): the refusals that need other namespaces.
    if (in.terrain.representation == "immersed_forcing") {
        ImmersedTerrainFacts itf;
        itf.surface_layer = surface_layer;
        itf.neutral = surf.stability == SurfaceLayerOptions::Stability::neutral;
        itf.z0_map = !surf.z0_file.empty();
        itf.buildings = buildings_opt.active();
        itf.nests = role.nests != nullptr || nest != nullptr;
        itf.newton = newton;
        // The sweep and the Poisson backend act in a stage only.
        itf.zline_sweeps = !solve || sopt.sweep == SweepKind::zjacobi;
        itf.poisson_mlmg = solve && popt.backend == "mlmg";
        itf.laminar = copt.kind == ClosureKind::none;
        itf.unsteady = unsteady;
        itf.canopy = canopy_opt.active();
        {
            amrex::ParmParse ppbc("bc");
            std::vector<amrex::Real> th;
            ppbc.queryarr("theta_profile", th);
            for (const amrex::Real x : th) { itf.stratified = itf.stratified || x != th[0]; }
            itf.zref_given = ppbc.contains("surface_zref");
            itf.floor_given = ppbc.contains("surface_denominator_floor");
        }
        itf.wall_k_source = surf.wall_k == SurfaceLayerOptions::WallK::source;
        itf.column_guess = iopt.first_guess == InitOptions::FirstGuess::column;
        itf.probes = !sopt.probes_file.empty();
        for (auto& e : check_immersed_terrain(itf)) { err.push_back(std::move(e)); }
    }
    // turb.closure = none (laminar, Closure.H): the options that model turbulence, or need the k-L closure, are refused.
    if (solve && copt.kind == ClosureKind::none) {
        if (canopy_opt.active()) { err.emplace_back("canopy.type needs turb.closure = kl (the canopy acts through the k-L closure's k terms)"); }
        if (surface_layer) { err.emplace_back("bc.zlo = surface_layer needs turb.closure = kl (a laminar run has no wall law); use no_slip_wall or slip_wall"); }
        if (iopt.first_guess == InitOptions::FirstGuess::column) {
            err.emplace_back("init.first_guess = column needs turb.closure = kl (the column solves the k-L closure); use profile");
        }
        if (newton) { err.emplace_back("steady.stage = simplec_newton is not supported with turb.closure = none"); }
        if (buildings_opt.active() && buildings_opt.wall_form != BuildingWallForm::drag) {
            err.emplace_back("buildings.wall_form = " + buildings_opt.wall_form_name + " needs turb.closure = kl (a wall law); use drag, "
                             "the no-slip staircase wall, with turb.closure = none");
        }
        if (buildings_opt.active() && buildings_opt.wall_tke) {
            err.emplace_back("buildings.wall_tke = on needs turb.closure = kl (k is passive in a laminar run)");
        }
    }
    if (solve && buildings_opt.active() && simplec.theta_frozen && !buildings_opt.theta_mode_name.empty()
        && buildings_opt.theta_mode != BuildingThetaMode::none) {
        err.emplace_back("buildings.theta_mode = " + buildings_opt.theta_mode_name + " has no effect with simplec.theta = frozen "
                         "(theta is data, not solved; leave the key out or give none)");
    }
    if (solve && canopy_opt.active() && canopy_opt.heat.on() && !simplec.theta_frozen
        && std::find(canopy_opt.given.begin(), canopy_opt.given.end(), "heat_frozen") != canopy_opt.given.end()) {
        err.emplace_back("canopy.heat_frozen has no effect without simplec.theta = frozen (theta is solved, and the leaf "
                         "heat acts in its equation)");
    }
    for (const auto& e : unknown_inputs(parmparse_entries())) {
        // The stage's namespaces were reported above with their reason.
        bool reported = false;
        if (!solve) {
            for (const char* ns : {"'forcing.", "'convection.", "'simplec.", "'turb."}) {
                reported = reported || e.find(ns) != std::string::npos;
            }
        }
        // The Newton keys without the Newton stage: reported above.
        if (!newton) { reported = reported || e.find("'steady.newton_") != std::string::npos; }
        // The unsteady keys without the unsteady stage: reported above.
        if (!unsteady) { reported = reported || e.find("'unsteady.") != std::string::npos; }
        if (!reported) { err.push_back(e); }
    }
    if (solve && !initialize && in.io.restart.empty()) {
        err.emplace_back("steady.stage = " + std::string(newton ? "simplec_newton" : unsteady ? "unsteady" : "simplec") + " needs a first guess: init.first_guess = profile, column or nest "
                         "or io.restart (the state would start at rest with theta = 0)");
    }
    if (initialize && !in.io.restart.empty()) {
        err.emplace_back("init.first_guess and io.restart are exclusive (a restart keeps its state)");
    }
    // The nested patch (NestParent.H, NestFaces.H): the parent, its first guess and the nest faces go together.
    const bool nest_guess = iopt.first_guess == InitOptions::FirstGuess::nest;
    const bool nest_face = has_nest_face(bin);
    for (auto& e : check_nest_run(bin, in, nest_guess, solve && seq.enabled)) { err.push_back(std::move(e)); }
    // A patch's frozen theta is the parent run's, which may already carry the leaves' heat; the patch's stage would
    // then blend it again (the same reason as for nest.names, below).
    if (solve && nest == nullptr && nest_guess && simplec.theta_frozen && canopy_opt.active() && canopy_opt.heat.on()
        && canopy_opt.heat.frozen == CanopyFrozenHeat::blend) {
        err.emplace_back("init.first_guess = nest is not supported with simplec.theta = frozen and canopy.heat_frozen = "
                         "blend (the parent run's theta may already carry the leaves' heat, and the patch would blend it again; "
                         "use canopy.heat_frozen = skip)");
    }
    // Frozen theta met the unsteady stage nowhere before (it needed grid sequencing, which that stage refuses).
    if (solve && nest == nullptr && nest_guess && simplec.theta_frozen && unsteady) {
        err.emplace_back("simplec.theta = frozen is not supported with steady.stage = unsteady and init.first_guess = nest "
                         "(a held theta in an unsteady patch is not tested)");
    }
    for (const auto& e : check_column_guess(iopt, sopt, fopt, bin, surf, surface_layer)) { err.push_back(e); }
    if (nest == nullptr && amrex::ParmParse("io").contains("column_guess_file") && iopt.first_guess != InitOptions::FirstGuess::column) {
        err.emplace_back("io.column_guess_file has no effect without init.first_guess = column");
    }
    for (const auto& e : check_stable_functions(iopt, surf, surface_layer)) { err.push_back(e); }
    // A nest: the outer domain's check covers the deck's set-up.
    if (nest == nullptr) {
        AblSetupFacts facts;
        facts.bc = bin;
        facts.forcing = fopt;
        facts.init = iopt;
        facts.surface = surf;
        facts.simplec = simplec;
        if (solve && abl.active() && bin.inflow.kind == InflowProfileOptions::Kind::file && !bin.inflow.file.empty()) {
            facts.inflow_table_theta = read_inflow_table(bin.inflow.file).theta;
        }
        facts.max_geom_lscale = copt.max_geom_lscale;
        facts.surface_layer = surface_layer;
        facts.restart = !in.io.restart.empty();
        facts.solve = solve;
        facts.canopy_heat = canopy_opt.active() && canopy_opt.heat.on();
        for (auto& e : check_abl_setup(abl, facts)) { err.push_back(std::move(e)); }
    }
    if (role.nests != nullptr) {
        for (const auto& e : role.nests->errors) { err.push_back(e); }
        for (auto& e : check_nests(*role.nests, in.mesh)) { err.push_back(std::move(e)); }
        if (role.nests->damping == NestDamping::outer
            && !(bin.rayleigh_zdamp > 0 && bin.rayleigh_dampcoef > 0
                 && (bin.rayleigh_fields[0] || bin.rayleigh_fields[1] || bin.rayleigh_fields[4]))) {
            err.emplace_back("nest.damping = outer has no effect without a Rayleigh layer on u, v or tke (bc.rayleigh_zdamp > 0, "
                             "bc.rayleigh_dampcoef > 0, bc.rayleigh_fields)");
        }
        // A nest's frozen theta is the parent's, already blended towards the leaves; the nest's stage would blend it
        // again (Simplec.cpp, canopy.heat_frozen = blend takes the theta of its first iteration as the frozen field).
        if (simplec.theta_frozen && canopy_opt.active() && canopy_opt.heat.on()
            && canopy_opt.heat.frozen == CanopyFrozenHeat::blend) {
            err.emplace_back("nest.names is not supported with simplec.theta = frozen and canopy.heat_frozen = blend (a "
                             "nest would blend the outer domain's blended theta again)");
        }
        if (!solve || unsteady) {
            err.emplace_back("nest.names needs steady.stage = simplec or simplec_newton (each nest starts from the outer "
                             "domain's steady state)");
        }
        if (in.terrain.smoothing == 1 && in.terrain.type != "flat") {
            err.emplace_back("nest.names is not supported with terrain.smoothing = 1 over terrain (STF smooths the levels over each "
                             "mesh's own ground, so a nest's heights above the ground differ from the outer domain's at "
                             "its faces): use terrain.smoothing = 0 or 2");
        }
        if (in.io.write_checkpoint != 1) {
            err.emplace_back("nest.names needs io.write_checkpoint = 1 (each nest starts from the outer domain's checkpoint)");
        }
        if (buildings_opt.active()) { err.emplace_back("nest.names is not supported with buildings.type yet"); }
        if (role.nests->feedback != NestFeedback::none) {
            const std::string fbk = std::string("nest.feedback = ") + nest_feedback_name(role.nests->feedback);
            // The joint run (JointRun.H) steps every domain's SIMPLEC driver; the Newton stage's residual differs
            // from the one the source is formed with. io.restart: Run.cpp reads the checkpoint's joint record.
            if (newton) { err.push_back(fbk + " needs steady.stage = simplec (not simplec_newton)"); }
            const bool inject = role.nests->feedback == NestFeedback::inject;
            // The injection holds cells in the pressure correction (Poisson.H, held_fluid); mlmg's operator has none.
            if (inject && popt.backend == "mlmg") {
                err.push_back(fbk + " needs solver.poisson_backend = gmres_mlmg or gmres_fft (mlmg's operator has no "
                                    "held cells)");
            }
            // The restriction coarsens each nest's boxes by 2 2 1 (restrict_patch). A nest has an even number of
            // cells; AMReX's maxSize keeps the boxes even when mesh.max_grid_size is even, not when it is odd.
            if (check_nests(*role.nests, in.mesh).empty()) {
                for (const auto& nest : role.nests->nests) {
                    // The source keeps kDefectMargin outer cells from each face (defect_region, NestDefect.cpp):
                    // checked here so a narrow nest stops before the one-way runs, not after them.
                    const MeshInputs nm = nest_mesh(nest, in.mesh);
                    for (int d = 0; d < 2; ++d) {
                        if (nm.is_periodic[d] == 0 && nm.n_cell[d] / 2 < 2 * kDefectMargin + 1) {
                            err.push_back(fbk + ": nest." + nest.name + " spans "
                                          + std::to_string(nm.n_cell[d] / 2) + " outer cells in "
                                          + (d == 0 ? "x" : "y") + "; it needs at least "
                                          + std::to_string(2 * kDefectMargin + 1) + " (the feedback keeps "
                                          + std::to_string(kDefectMargin) + " cells from each face)");
                        }
                    }
                    // A held region across the whole of a periodic direction cuts the outer fluid in two, each part
                    // with its own pressure level; the masked pressure correction removes one mean only.
                    for (int d = 0; d < 2; ++d) {
                        if (inject && nm.is_periodic[d] != 0) {
                            err.push_back(fbk + ": nest." + nest.name + " spans the periodic outer direction "
                                          + (d == 0 ? "x" : "y") + "; the injection needs a nest with faces in x and y");
                        }
                    }
                    if (!make_box_array(nm).coarsenable(amrex::IntVect(2, 2, 1))) {
                        err.push_back(fbk + ": nest." + nest.name + " splits into boxes with an odd "
                                      "number of cells in x or y at mesh.max_grid_size: give an even "
                                      "mesh.max_grid_size");
                    }
                }
            }
            // The nests damp w and theta only unless nest.damping = outer (Nests.cpp). Over each footprint the
            // source carries the nest's equations, so an outer Rayleigh term on u, v or k the nests do not have would
            // act outside the footprints and not inside: a jump in the forcing at the footprint's edge above
            // bc.rayleigh_zdamp.
            if (role.nests->damping != NestDamping::outer && bin.rayleigh_zdamp > 0
                && (bin.rayleigh_fields[0] || bin.rayleigh_fields[1] || bin.rayleigh_fields[4])) {
                err.push_back(fbk + " needs bc.rayleigh_fields without u, v and tke (the nests damp w "
                              "and theta only): use bc.rayleigh_fields = w or w theta, or nest.damping = outer");
            }
            if (bin.sponge_width > 0) {
                amrex::Print() << "gust: note: " << fbk << " with bc.sponge_width > 0: the nests have no side "
                                  "sponge, so over a nest's footprint inside the sponge the joint run does not damp "
                                  "(place the nests outside bc.sponge_width)\n";
            }
        }
    }
    const std::string invalid = nest != nullptr ? "gust: invalid inputs for nest." + nest->name + " (the outer deck with "
                                                  "the nest's derived settings, Nests.H):"
                                                : std::string("gust: invalid inputs:");
    if (!err.empty()) {
        std::string all = invalid;
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }

    // The checks that need the geometry and the layout (nothing is allocated): a nest's check-only pass runs them
    // too, so that a nest stops the run before the outer domain starts.
    const amrex::Geometry geom = make_geometry(in.mesh);
    // Buildings in a nested patch keep clear of its nest faces, the relaxation zone and the ground blend (NestBuildings.H).
    if (nest_face && buildings_opt.active()) {
        const auto near = check_nest_buildings(read_building_map(buildings_opt.file, buildings_opt.format), buildings_opt.file,
                                               bin, geom, in.terrain.nest_blend_width);
        if (!near.empty()) {
            std::string all = "gust: invalid inputs:";
            for (const auto& e : near) { all += "\n  " + e; }
            amrex::Abort(all);
        }
    }
    // terrain.representation = immersed_forcing (ImmersedTerrain.H): the inputs are valid, so the ground takes the
    // buildings' wall law from here on, through options of its own; its restart record is its own (immersed_hash).
    const bool immersed = in.terrain.representation == "immersed_forcing";
    std::string immersed_hash;
    if (immersed) {
        buildings_opt = immersed_terrain_options(in.terrain, surf.z0);
        immersed_hash = immersed_terrain_hash(in.terrain, surf.z0);
    }
    // With a grid hierarchy, boxes every level coarsens exactly.
    const amrex::BoxArray ba = hierarchy_box_array(in.mesh);
    // Grid sequencing runs the coarse levels on a fresh start only.
    const bool sequencing = solve && seq.enabled && in.io.restart.empty();
    // A nest holds theta at its start: a joint run's nest the outer domain's frozen theta, interpolated; a patch
    // that starts from a parent run (init.first_guess = nest) that run's theta, interpolated, whether the parent
    // solved theta or held it.
    if (simplec.theta_frozen && nest == nullptr && !nest_guess && !(sequencing && in.mesh.hierarchy_levels > 1)) {
        amrex::Abort(invalid + "\n  simplec.theta = frozen needs grid sequencing (steady.sequencing = 1, "
                     "mesh.hierarchy_levels > 1): the frozen theta is set from the inflow profile (or the stratified "
                     "column) on the sequenced levels; or init.first_guess = nest: theta is held at the parent run's");
    }
    if (solve && sopt.sweep == SweepKind::lusgs) {
        // The overlap against the layout (LineSweeps.H), before any work;
        // with sequencing also against every coarse level's boxes.
        std::string lay = check_line_overlap_layout(ba, sopt.line_overlap);
        for (int l = 1; sequencing && lay.empty() && l < in.mesh.hierarchy_levels; ++l) {
            const int r = hierarchy_ratio(l);
            lay = check_line_overlap_layout(amrex::BoxArray(ba).coarsen(amrex::IntVect(AMREX_D_DECL(r, r, 1))), sopt.line_overlap);
            if (!lay.empty()) { lay += " (on grid hierarchy level " + std::to_string(l) + ", steady.sequencing = 1)"; }
        }
        if (!lay.empty()) { amrex::Abort(invalid + "\n  " + lay); }
    }
    if (role.check_only && (initialize || solve)) {
        // The boundary checks BoundaryConditions runs once the mesh is built, on the z-levels it will have.
        // A nest's pass has the deck's forcing; its run takes the outer run's (OuterOutcome), which these
        // checks do not test (they test the faces and the Coriolis flags, not the forcing values).
        const auto berr = check_boundary_inputs(bin, geom, fopt.coriolis, make_zlevels(in.mesh, geom).back());
        if (!berr.empty()) {
            std::string all = invalid;
            for (const auto& e : berr) { all += "\n  " + e; }
            amrex::Abort(all);
        }
    }
    if (role.after_checks) { role.after_checks(); }
    if (role.check_only) { return final_chk; }
    // The ranks of the boxes. AMReX's default map gives the heaviest boxes to the ranks with the least memory in use,
    // so a nest's map would follow what the outer run holds when the nest is built, which a joint restart changes.
    // A joint run's nest takes the same space-filling-curve map in rank order: the map then depends on the boxes and
    // the rank count alone, and a joint restart sums over the same ranks as the run it continues.
    const amrex::DistributionMapping dm = [&] {
        if (role.nest == nullptr || !role.joint) { return amrex::DistributionMapping(ba); }
        amrex::Vector<amrex::Real> cost(static_cast<std::size_t>(ba.size()));
        for (int i = 0; i < ba.size(); ++i) { cost[static_cast<std::size_t>(i)] = static_cast<amrex::Real>(ba[i].numPts()); }
        return amrex::DistributionMapping::makeSFC(cost, ba, false);
    }();

    if (abl.active()) {
        amrex::Print() << "gust: set-up: abl.drive = " << abl_drive_name(abl.drive) << ", abl.theta = "
                       << abl_theta_name(abl.theta) << "\n";
    }
    amrex::Print() << "gust: " << in.mesh.n_cell << " cells, " << ba.size()
                   << " boxes (each spans the full z column)\n";

    State state;
    state.define(ba, dm, immersed);
    Mesh mesh;
    mesh.define(ba, dm);
    const bool restart = !in.io.restart.empty();
    // The canopy's drag field on level 0's mesh (Canopy.H), once the node heights are there: built from the inputs
    // or, on a restart, after the checkpoint's canopy hash has been checked against them (before the field is
    // built). The column first guess takes the canopy's layer profile (ColumnGuess.H).
    Canopy canopy;
    Buildings buildings;
    auto build_canopy = [&] () {
        if (restart) {
            Buildings::check_restart(in.io.restart, buildings_opt.active() && !immersed ? Buildings::hash_inputs(buildings_opt) : "");
            check_immersed_terrain_restart(in.io.restart, immersed_hash);
        }
        if (immersed) {
            const auto& zl = mesh.zlevels();
            const std::string low = check_immersed_terrain_height(in.terrain, static_cast<double>(zl[1]),
                                                                  static_cast<double>(zl[zl.size() - 2]));
            if (!low.empty()) { amrex::Abort("gust: " + low); }
            // The ground's height at the nodes, from the inputs (a restart's too: its record pins them).
            fill_immersed_ground(state.z_phys_nd, geom, in.terrain);
            const std::string out = check_immersed_ground_extent(state.z_phys_nd, static_cast<double>(zl[1]),
                                                                 static_cast<double>(zl[zl.size() - 2]));
            if (!out.empty()) { amrex::Abort("gust: " + out); }
            buildings.build_from_surface(buildings_opt, "terrain.representation = immersed_forcing", geom, state.z_phys_nd, ba, dm);
            amrex::Print() << "gust: immersed terrain (terrain.representation = immersed_forcing, terrain.type = "
                           << in.terrain.type << ", terrain.flat_height = " << in.terrain.flat_height
                           << " m, terrain.small_volfrac = " << in.terrain.small_volfrac << ", ground from "
                           << state.z_phys_nd.min(state.z_phys_nd.nComp() - 1) << " m to "
                           << state.z_phys_nd.max(state.z_phys_nd.nComp() - 1) << " m): "
                           << buildings.summary(mesh.detJ, geom) << "\n"
                           << "gust: immersed terrain: the form 7 wall law with bc.surface_z0 = " << buildings_opt.z0
                           << (in.terrain.immersed_slope == "normal" ? " m at the fluid's mean distance from the ground along its normal"
                                                                      : " m at the fluid's mean height above the ground")
                           << " (terrain.immersed_slope = " << in.terrain.immersed_slope
                           << "), the wall value of k in the cut cells\n";
        } else if (buildings_opt.active()) {
            buildings.build(buildings_opt, geom, state.z_phys_nd, ba, dm);
            amrex::Print() << "gust: buildings (buildings.type = " << buildings_opt.type_name << "): "
                           << buildings.summary(mesh.detJ, geom) << "\n";
            // The choices as the run takes them, their defaults included.
            amrex::Print() << "gust: buildings: with a stage, buildings.wall_form = " << buildings_opt.wall_form_name;
            if (buildings_opt.wall_form != BuildingWallForm::drag) {
                amrex::Print() << ", buildings.wall_speed = "
                               << (buildings_opt.wall_speed == BuildingWallSpeed::open_faces ? "open_faces" : "all_faces");
            }
            amrex::Print() << ", buildings.wall_tke = " << (buildings_opt.wall_tke ? "on" : "off");
            // theta's terms act where theta is solved: with simplec.theta = frozen it is data. (Reached in
            // a patch of a parent run, init.first_guess = nest: rt_nest_faces_frozen_building. The guards here, in
            // building_heat below and in SimplecStage::compute_norms also serve buildings on sequenced levels.)
            const bool heat = buildings_opt.theta_mode != BuildingThetaMode::none && !simplec.theta_frozen;
            amrex::Print() << ", buildings.theta_mode = "
                           << (!heat ? "none" : buildings_opt.theta_mode == BuildingThetaMode::fixed ? "fixed" : "ambient");
            if (heat && buildings_opt.theta_mode == BuildingThetaMode::fixed) {
                amrex::Print() << " (buildings.theta_surface " << buildings_opt.theta_surface << " K)";
            }
            if (!heat && buildings_opt.theta_mode != BuildingThetaMode::none) { amrex::Print() << " (simplec.theta = frozen)"; }
            amrex::Print() << "\n";
        }
        if (restart) { Canopy::check_restart(in.io.restart, canopy_opt.active() ? Canopy::hash_inputs(canopy_opt) : ""); }
        if (!canopy_opt.active()) { return; }
        canopy.build(canopy_opt, geom, state.z_phys_nd, ba, dm);
        amrex::Print() << "gust: canopy (canopy.type = " << canopy_opt.type_name << "): " << canopy.summary(mesh.detJ, geom)
                       << "\ngust: canopy: the drag acts on the momentum equations when a stage runs, "
                          "canopy.drag_diagonal = " << canopy_diagonal_name(canopy_opt.drag_diagonal) << "\n";
        const CanopyTkeOptions& t = canopy_opt.tke;
        if (t.on) {
            amrex::Print() << "gust: canopy: k terms on (canopy.tke): beta_p " << t.beta_p << ", beta_d " << t.beta_d
                           << ", length cap alpha / f with canopy.length_alpha " << t.length_alpha << "\n";
        } else {
            amrex::Print() << "gust: canopy: k terms off (canopy.tke = off): the drag only\n";
        }
        const CanopyHeatOptions& hh = canopy_opt.heat;
        if (!hh.on()) {
            amrex::Print() << "gust: canopy: leaf heat off (canopy.heat = off)\n";
        } else if (!solve) {
            amrex::Print() << "gust: canopy: leaf heat (canopy.heat = fixed, canopy.leaf_theta " << hh.leaf_theta
                           << " K) acts on theta when a stage runs; the leaf area density is written as canopy_lad\n";
        } else if (!simplec.theta_frozen) {
            amrex::Print() << "gust: canopy: leaf heat on (canopy.heat = fixed): theta relaxes towards canopy.leaf_theta "
                           << hh.leaf_theta << " K at the rate 2 a / r_b\n";
        } else if (hh.frozen == CanopyFrozenHeat::skip) {
            amrex::Print() << "gust: WARNING: canopy.heat = fixed with simplec.theta = frozen and canopy.heat_frozen = skip: "
                           << "theta is not solved, so the leaf heat is left out; canopy.heat_frozen = blend sets theta in "
                           << "the canopy from the leaf exchange instead\n";
        } else {
            amrex::Print() << "gust: canopy: leaf heat with theta frozen (canopy.heat_frozen = blend): in the canopy theta = "
                           << "theta_f + w (" << hh.leaf_theta << " K - theta_f), w = lambda k / (eps + lambda k)\n";
        }
    };
    std::unique_ptr<InflowProfile> column;
    // A nested patch's parent, read once: the first guess, the nest faces' values and the ground blend.
    std::unique_ptr<NestParent> parent;
    std::unique_ptr<NestFaces> nest_faces;
    Mesh::GroundEdit ground_edit;
    Mesh::LevelsEdit levels_edit;
    const bool joint_nest = role.nest != nullptr && static_cast<bool>(role.joint);
    std::string nest_hash;
    if (nest_guess || nest_face) {
        const amrex::Real blend = in.terrain.nest_blend_width;
        parent = std::make_unique<NestParent>(bin.nest_parent, in.mesh, make_zlevels(in.mesh, geom), ba, dm,
                                              blend > 0 || joint_nest, bin.nest_ratio);
        amrex::Print() << parent->summary();
        // The hill where the parent has it, not in the middle of the patch.
        if (in.terrain.type == "agnesi") { in.terrain.hill_centre = parent->domain_centre(); }
        if (blend > 0) { ground_edit = [&, blend] (amrex::MultiFab& z) { blend_nest_ground(*parent, bin, blend, geom, z); }; }
        // A joint run's nest: its coordinate surfaces follow the parent's aloft (NestFaces.H).
        if (joint_nest) {
            levels_edit = [&] (amrex::MultiFab& z, const amrex::Vector<amrex::Real>& zl) { follow_parent_levels(*parent, geom, zl, z, role_name); };
        }
        if (nest_face) { nest_hash = NestFaces::hash_inputs(*parent, bin, blend, joint_nest); }
    }
    auto build_nest_faces = [&] () {
        if (!nest_face) { return; }
        nest_faces = std::make_unique<NestFaces>(*parent, bin, geom, mesh);
        amrex::Print() << nest_faces->summary();
    };
    if (!restart) {
        mesh.build(in, geom, state.z_phys_nd, ground_edit, levels_edit);
        build_canopy();
        build_nest_faces();
        if (initialize) {
            // The converged column as the profile (ColumnGuess.H), then the
            // mass-consistent initial field (Initializer.H) or the first
            // guess alone (init.projection = 0). With grid sequencing the
            // initializer runs on the coarsest level instead.
            if (iopt.first_guess == InitOptions::FirstGuess::column) {
                // With a canopy the column carries its layer profile (ColumnGuess.H; a canopy that differs between
                // columns stops the run there).
                const ColumnCanopy ccan = canopy.active() ? column_canopy(canopy, geom, in.terrain.type == "flat") : ColumnCanopy();
                if (canopy.active()) {
                    int layers = 0;
                    double fmax = 0;
                    for (double f : ccan.f) { layers += (f > 0) ? 1 : 0; fmax = std::max(fmax, f); }
                    amrex::Print() << "gust: column guess: with the canopy's layer profile (f > 0 in " << layers << " of "
                                   << ccan.f.size() << " layers, largest f " << fmax << " 1/m; the k terms "
                                   << (ccan.tke.on ? "on" : "off")
                                   << (canopy.heat().on() ? "; the leaf heat is not in the column" : "") << ")\n";
                }
                const ColumnGuessResult cg = run_column_guess(in, iopt, sopt, bin, surf, surface_layer, copt, fopt,
                                                            in.io.column_guess_file, ccan);
                char cl[512];
                if (cg.time > 0) {
                    std::snprintf(cl, sizeof(cl), "gust: column guess: marched to t = %.6g s in %d steps (%.3f s)%s, %d levels, "
                                  "u* %.4f m/s, u %.4f v %.4f m/s at the top\n", cg.time, cg.iterations, cg.wall_seconds,
                                  iopt.column_wall_k == InitOptions::ColumnWallK::solved ? " with the first cell's k solved" : "",
                                  static_cast<int>(cg.z.size()) - 1, static_cast<double>(cg.ustar),
                                  static_cast<double>(cg.u.back()), static_cast<double>(cg.v.back()));
                } else {
                    std::snprintf(cl, sizeof(cl), "gust: column guess: converged after %d steps in %d columns (%.3f s), %d levels, "
                                  "u* %.4f m/s, u %.4f v %.4f m/s at the top\n", cg.iterations, cg.columns, cg.wall_seconds,
                                  static_cast<int>(cg.z.size()) - 1, static_cast<double>(cg.ustar),
                                  static_cast<double>(cg.u.back()), static_cast<double>(cg.v.back()));
                }
                amrex::Print() << cl;
                // What a restart deck needs (init.* keys are refused there and the checkpoint holds the state only).
                std::string restart_lines;
                if (column_theta_evolves(iopt, surf)) {
                    const double L = cg.inv_L != 0 ? 1.0 / static_cast<double>(cg.inv_L) : std::numeric_limits<double>::infinity();
                    if (cg.time > 0) {
                        std::snprintf(cl, sizeof(cl), "gust: column guess: at t = %.6g s: surface heat flux %.4e K m/s, L %.4g m, "
                                      "surface theta %.4f K\n", cg.time, static_cast<double>(cg.heat_flux), L,
                                      static_cast<double>(cg.theta_surface));
                    } else {
                        std::snprintf(cl, sizeof(cl), "gust: column guess: surface heat flux %.4e K m/s, L %.4g m, surface theta %.4f K "
                                      "after %d surface adjustments\n", static_cast<double>(cg.heat_flux), L,
                                      static_cast<double>(cg.theta_surface), cg.flux_iterations);
                    }
                    amrex::Print() << cl;
                    const bool mol = iopt.column_stability == InitOptions::ColumnStability::mol;
                    if (surf.stability == SurfaceLayerOptions::Stability::heat_flux) {
                        if (mol) {   // the 3-D surface takes the column's flux (check_column_guess)
                            surf.heat_flux = cg.heat_flux;
                            amrex::Print() << "gust: column guess: the run takes this flux as bc.surface_heat_flux\n";
                            std::snprintf(cl, sizeof(cl), "  bc.surface_heat_flux = %.9g\n", static_cast<double>(cg.heat_flux));
                            restart_lines += cl;
                        }
                    } else if (cg.time > 0) {
                        amrex::Print() << "gust: column guess: the 3-D surface is neutral; the flux above is the time-marching "
                                          "column's at init.column_time only\n";
                    } else {
                        amrex::Print() << "gust: column guess: the 3-D surface is neutral; for the same surface set "
                                          "bc.surface_stability = heat_flux"
                                       << (mol ? " (the run then takes this flux)\n"
                                               : " and bc.surface_heat_flux to this flux with init.column_stability = heat_flux\n");
                    }
                }
                if (cg.target_iterations > 0) {
                    // The target wind: the 3-D run takes the forcing the column was adjusted to.
                    if (fopt.coriolis.enabled && fopt.coriolis.geostrophic) {
                        fopt.coriolis.ug = cg.ug;
                        fopt.coriolis.vg = cg.vg;
                        std::snprintf(cl, sizeof(cl), "  forcing.ug = %.9g\n  forcing.vg = %.9g\n", static_cast<double>(cg.ug),
                                      static_cast<double>(cg.vg));
                        restart_lines += cl;
                    }
                    fopt.pressure_gradient[0] = cg.pressure_gradient_x;
                    fopt.pressure_gradient[1] = cg.pressure_gradient_y;
                    if (cg.pressure_gradient_x != 0 || cg.pressure_gradient_y != 0) {
                        std::snprintf(cl, sizeof(cl), "  forcing.pressure_gradient = %.9g %.9g %.9g\n",
                                      static_cast<double>(cg.pressure_gradient_x), static_cast<double>(cg.pressure_gradient_y),
                                      static_cast<double>(fopt.pressure_gradient[2]));
                        restart_lines += cl;
                    }
                    std::snprintf(cl, sizeof(cl), "gust: column guess: target wind reached after %d forcing adjustments; the run "
                                  "takes geostrophic wind (%.6g, %.6g) m/s and pressure gradient (%.6g, %.6g) m/s^2\n",
                                  cg.target_iterations - 1, static_cast<double>(cg.ug), static_cast<double>(cg.vg),
                                  static_cast<double>(cg.pressure_gradient_x), static_cast<double>(cg.pressure_gradient_y));
                    amrex::Print() << cl;
                }
                if (!restart_lines.empty()) {
                    amrex::Print() << "gust: column guess: a restart from this run's checkpoints needs these inputs in place of "
                                      "the init.* keys:\n" << restart_lines;
                }
                if (cg.time == 0 && cg.tol > static_cast<double>(iopt.column_tol)) {
                    std::snprintf(cl, sizeof(cl), "gust: column guess: init.column_tol = %.3e lifted to %.3e (100 epsilon of Real)\n",
                                  static_cast<double>(iopt.column_tol), cg.tol);
                    amrex::Print() << cl;
                }
                column = std::make_unique<InflowProfile>(cg.z, cg.u, cg.v, cg.theta, cg.tke, bin.theta_ref);
            }
        }
        if (iopt.first_guess == InitOptions::FirstGuess::nest) {
            // The nested patch: the parent at bc.nest_ratio (r r 1), then the initializer's correction.
            parent->interpolate(state);
            // The nest faces' normal velocities after the net flux is removed (NestFaces.H).
            if (nest_faces) { for (int d = 0; d < 2; ++d) { nest_faces->fill(state, d, false); } }
            zero_ground_flux(state, mesh, geom);   // the patch's ground carries no flux, as in first_guess_state
            const InitResult r = correct_initial_state(iopt, bin, popt, mesh, geom, state);
            print_initializer(r, iopt, popt, "");
            if (simplec.theta_frozen && nest == nullptr) {
                // The patch of a parent run: the line the sequenced levels print (reset_frozen_theta, Sequencing.cpp).
                const amrex::Real tmin = state.theta.min(0), tmax = state.theta.max(0);
                const double tmean = static_cast<double>(state.theta.sum(0)) / static_cast<double>(ba.numPts());
                char tl[200];
                std::snprintf(tl, sizeof(tl), "gust: theta frozen at the parent run's theta (bc.nest_parent), interpolated: "
                              "min %.15e, max %.15e, mean %.6f K", static_cast<double>(tmin), static_cast<double>(tmax), tmean);
                amrex::Print() << tl << "\n";
            }
        } else if (initialize && !sequencing) {
            const InitResult r = initialize_state(iopt, bin, popt, mesh, geom, state, nullptr, fopt.coriolis, column.get(),
                                                  column_sets_theta(iopt, surf));
            print_initializer(r, iopt, popt, "");
        }
    } else {
        read_checkpoint(in.io.restart, state, in.mesh);
        NestFaces::check_restart(in.io.restart, nest_hash);
        mesh.check_restart(in, geom, state.z_phys_nd, ground_edit, levels_edit);
        amrex::Print() << "gust: restarted from " << in.io.restart
                       << " at step " << state.step << "\n";
        build_canopy();
        build_nest_faces();
    }
    // The parent's images are not needed again: NestFaces holds its own copy on the patch.
    ground_edit = nullptr;
    levels_edit = nullptr;
    if (!role.joint) { parent.reset(); }
    const amrex::MultiFab* canopy_plot = canopy.active() ? &canopy.drag() : nullptr;
    const amrex::MultiFab* lad_plot = (canopy.active() && canopy.heat().on()) ? &canopy.lad() : nullptr;
    const amrex::MultiFab* beta_plot = buildings.active() ? &buildings.beta() : nullptr;
    const amrex::MultiFab* wall_plot = buildings.active() ? &buildings.wall() : nullptr;
    // Writing a checkpoint of the same name as the one read would move it to
    // <name>.old.NNNN: refused (checked at the final step, after the solve).
    // The file system decides whether the two paths name one directory
    // (same_directory), on the I/O rank; every rank takes its answer.
    auto check_overwrite = [&] {
        if (!restart || in.io.write_checkpoint != 1) { return; }
        int same = 0;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            same = same_directory(amrex::Concatenate(in.io.chk_file, state.step, 5), in.io.restart) ? 1 : 0;
        }
        amrex::ParallelDescriptor::Bcast(&same, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        if (same != 0) {
            amrex::Abort("gust: io.chk_file = " + in.io.chk_file + " would overwrite io.restart = "
                         + in.io.restart + "; choose another io.chk_file or set io.write_checkpoint = 0");
        }
    };
    if (!solve) { check_overwrite(); }

    amrex::Print() << mesh.quality_report();
    if (solve && simplec.w_vnn_cap > 0 && simplec.w_vnn_cap_slope_margin == 0 && mesh.quality().max_ground_slope_deg > 0
        && mesh.quality().max_aspect_ratio > w_vnn_cap_tested_aspect) {
        // Without the slope margin, the cap can destabilise the pressure
        // correction where thin cells meet slopes (Predictor.H); aspect ratio
        // is the proxy.
        char tested[32];   // two decimals in either precision (50.78, not 50.77999878)
        std::snprintf(tested, sizeof(tested), "%.2f", static_cast<double>(w_vnn_cap_tested_aspect));
        amrex::Print() << "WARNING: simplec.w_vnn_cap = " << simplec.w_vnn_cap << " over terrain with max aspect ratio "
                       << mesh.quality().max_aspect_ratio << ", above the tested " << tested
                       << " (Askervein, 1 m first cell): the cap may make the pressure correction unstable where thin "
                          "cells meet steep slopes with simplec.w_vnn_cap_slope_margin = 0 (docs/solution_algorithm.rst); the "
                          "default margin of " << SimplecOptions{}.w_vnn_cap_slope_margin
                       << " converged on the tested steep cases, and without the cap (simplec.w_vnn_cap = 0) w "
                          "converges slowly on thin first cells\n";
    }

    // The x/y grid hierarchy (mesh.hierarchy_levels > 1): built from level 0's
    // node heights, reported and written with the plotfile. With
    // steady.sequencing = 1 the stage runs on every level, coarse to fine
    // (Sequencing.H); otherwise on level 0 alone.
    std::unique_ptr<MeshHierarchy> hierarchy;
    if (in.mesh.hierarchy_levels > 1) {
        hierarchy = std::make_unique<MeshHierarchy>(in, geom, state.z_phys_nd, mesh);
        amrex::Print() << hierarchy->quality_report();
    }
    auto write_hierarchy = [&] {
        if (!hierarchy) { return; }
        for (const auto& name : hierarchy->write_plotfiles(in.io.plot_file, state.step)) {
            amrex::Print() << "gust: wrote " << name << "\n";
        }
    };

    if (solve) {
        // The run's models on level 0 (Sequencing.H: LevelSolver).
        // After a column guess's adjustments: the values the nests take.
        if (role.outcome != nullptr) { keep_outer_forcing(fopt, surf, *role.outcome); }
        RunSetup setup;
        setup.bin = bin;
        setup.popt = popt;
        setup.surf = surf;
        setup.surface_layer = surface_layer;
        setup.idw_neighbors = in.terrain.idw_neighbors;
        setup.idw_power = in.terrain.idw_power;
        setup.copt = copt;
        setup.fopt = fopt;
        setup.cin = cin;
        setup.simplec = simplec;
        setup.sopt = sopt;
        setup.newton = newton_opt;
        setup.nest = nest_faces.get();
        setup.log = iteration_log(in.io.iteration_log);
        BuildingForcing building_forcing;
        BuildingTke building_tke;
        BuildingTheta building_theta;
        const bool building_heat = buildings.active() && buildings_opt.theta_mode != BuildingThetaMode::none && !simplec.theta_frozen;
        if (buildings.active()) {
            const amrex::Long sealed = count_sealed_fluid_cells(buildings, geom);
            if (sealed > 0) {
                amrex::Abort("gust: buildings.file = " + buildings_opt.file + ": " + std::to_string(sealed)
                             + " fluid cells are sealed off by buildings; the pressure correction cannot solve there. Open "
                               "the pocket, or fill it in the map");
            }
            building_forcing.build(buildings, buildings_opt, geom, mesh.detJ);
            if (buildings_opt.wall_tke) {
                building_tke.build(buildings, building_forcing, buildings_opt, geom, state.z_phys_nd, closure_constants(copt).tke_floor,
                                   copt.Cmu0);
            }
            if (building_heat) {
                building_theta.build(buildings, building_forcing, buildings_opt, geom, state.z_phys_nd, mesh.detJ, bin.theta_ref);
            }
        }
        // steady.stage = unsteady: the physical time term on level 0 (dt set by the time driver at every step;
        // with unsteady.cfl the first value is a placeholder until the first step).
        std::unique_ptr<TimeTerm> time_term;
        if (unsteady) {
            TimeTermOptions topt;
            topt.dt = uopt.cfl > 0 ? amrex::Real(1) : uopt.dt;
            time_term = std::make_unique<TimeTerm>(topt, ba, dm);
        }
        // nest.feedback = defect: the outer domain's source of its nests (NestDefect.H), zero until the joint run sets it.
        std::unique_ptr<NestDefect> defect;
        // nest.feedback = relax or relax_consistent: the outer domain's relaxation towards its nests (NestRelax.H) in its
        // place.
        std::unique_ptr<NestRelax> relax;
        // nest.feedback_flux = nest: the nests' face fluxes on the rows just outside the source (NestFaceFlux.H).
        std::unique_ptr<NestFaceFlux> face_flux;
        // nest.feedback = inject: the outer domain's injection of its nests (NestInject.H) in its place.
        std::unique_ptr<NestInject> inject;
        if (role.joint && nest == nullptr) {
            if (role.nests->feedback == NestFeedback::relax || role.nests->feedback == NestFeedback::relax_consistent) {
                relax = std::make_unique<NestRelax>(ba, dm, geom, role.nests->feedback_factor,
                                                    role.nests->feedback == NestFeedback::relax_consistent);
            } else if (role.nests->feedback == NestFeedback::inject) {
                inject = std::make_unique<NestInject>(ba, dm, geom, role.nests->feedback_factor);
            } else {
                defect = std::make_unique<NestDefect>(ba, dm);
                if (role.nests->feedback_flux == NestFeedbackFlux::nest) {
                    face_flux = std::make_unique<NestFaceFlux>(ba, dm, mesh, geom);
                }
            }
        }
        LevelSolver level0(setup, mesh, geom, state.z_phys_nd, ba, dm, "", newton, nullptr,
                           canopy.active() ? &canopy : nullptr, nullptr, buildings.active() ? &building_forcing : nullptr,
                           (buildings.active() && buildings_opt.wall_tke) ? &building_tke : nullptr,
                           building_heat ? &building_theta : nullptr, time_term.get(), defect.get(), relax.get(),
                           inject.get());
        if (buildings.active()) {
            // After the wall distance above: its Poisson solve does not see the buildings.
            mesh.set_open_faces(buildings.beta_face(0), buildings.beta_face(1), buildings.beta_face(2), buildings.beta(), geom);
        }
        NewtonStage* newton_stage = level0.newton_stage();
        const BoundaryConditions& bc = level0.bc;
        Closure& closure = level0.closure;
        const amrex::MultiFab& dist = level0.dist;
        const ResidualInputs& res_in = level0.res_in;
        const SimplecStage& simplec_stage = level0.stage();
        Driver& driver = level0.driver();
        // buildings.force_file (BuildingLoads.H): the force at the state, its ghost cells filled first.
        std::unique_ptr<BuildingLoads> loads;
        if (buildings.active() && !buildings_opt.force_file.empty()) {
            loads = std::make_unique<BuildingLoads>(buildings_opt.force_file);
            // A steady run starts its file now, so a path that cannot be written stops the run before the solve.
            if (!unsteady) { loads->open(-1.0); }
        }
        auto force_now = [&] () {
            bc.fill_ghosts(state, res_in.surface, fopt.buoyancy, false, nullptr);
            ResidualInputs fin = res_in;
            if (!simplec_stage.hybrid_sigma().empty()) { fin.hybrid_sigma = &simplec_stage.hybrid_sigma(); }
            return loads->evaluate(state, fin, mesh, geom, buildings);
        };
        if (restart) {
            // The histories: the driver's and the closure's nu_t (both
            // written with the checkpoint; the driver's must exist).
            driver.restore(read_driver_history(in.io.restart), state);
            closure.read_checkpoint(in.io.restart, ba, dm, geom);
            if (newton_stage != nullptr) { newton_stage->read_checkpoint(in.io.restart, ba, dm, geom); }
            // read_checkpoint restores the valid points only (and z_phys_nd's
            // ghost nodes); the stage's init fills the ghosts before the first
            // step (SimplecStage::init, NewtonStage::init).
            amrex::Print() << "gust: driver history restored: " << driver.history().records.size()
                           << " iterations, cfl " << driver.cfl() << ", reference hold " << driver.hold_iterations()
                           << " iterations" << (driver.history().hold == 0 ? " (released: a seeded level of grid sequencing)" : "")
                           << "\n";
            UnsteadyRecord urec;
            if (!unsteady && read_unsteady_record(in.io.restart, urec)) {
                amrex::Print() << "gust: the checkpoint is of an unsteady run at t = " << urec.time
                               << " s; the steady stage takes its state and driver history (its averages are not used)\n";
            }
        }
        amrex::Print() << "gust: SIMPLEC stage: " << sopt.zline_sweeps << " " << sweep_name(sopt.sweep) << " sweeps";
        if (sopt.sweep == SweepKind::lusgs) {
            amrex::Print() << " (ghosts " << line_ghosts_name(sopt.line_ghosts) << ", overlap " << sopt.line_overlap << ", "
                           << line_order_name(sopt.line_order) << " order"
                           << (sopt.line_alternate != 0 ? ", alternating" : "") << ")";
        }
        amrex::Print() << ", pressure " << popt.backend << " reltol "
                       << simplec.pressure_reltol << ", alpha_p " << simplec.alpha_p << ", w_vnn_cap ";
        if (simplec.w_vnn_cap > 0) {
            amrex::Print() << simplec.w_vnn_cap << " (slope margin ";
            if (simplec.w_vnn_cap_slope_margin > 0) { amrex::Print() << simplec.w_vnn_cap_slope_margin; } else { amrex::Print() << "off"; }
            amrex::Print() << ")";
        } else {
            amrex::Print() << "off";
        }
        amrex::Print() << ", cfl0 " << sopt.cfl0
                       << ", cfl_max " << sopt.cfl_max << "\n";
        if (newton) {
            amrex::Print() << "gust: Newton stage after " << sopt.switch_window << " iterations at r <= " << sopt.switch_residual
                           << (sopt.switch_monitor_tol > 0 ? " with the monitors settled" : "") << "; GMRES("
                           << newton_opt.newton.gmres_restart << "), " << newton_preconditioner_name(sopt.newton_preconditioner)
                           << " preconditioner, back to SIMPLEC after " << sopt.newton_fallback_failures << " failed steps\n";
        }
        // Grid sequencing: the coarse levels to the loose stop, then the
        // prolongation into the level-0 state (Sequencing.H).
        std::vector<LevelSummary> levels;
        if (sequencing) {
            amrex::Print() << "gust: grid sequencing over " << hierarchy->num_levels() << " levels: loose stop at "
                           << seq.drop_orders << " orders";
            if (seq.monitor_tol > 0) { amrex::Print() << " or monitors within " << seq.monitor_tol; }
            amrex::Print() << ", strict stop on level 0; "
                           << (seq.hold == 1 ? "every level holds steady.cfl0 for the reference iterations"
                                             : "the seeded levels without the reference hold")
                           << (seq.cfl_previous ? ", starting at the previous level's CFL" : "") << "\n";
            levels = sequence_coarse_levels(*hierarchy, setup, seq, iopt, column.get(), level0, state);
        } else if (solve && seq.enabled) {
            amrex::Print() << "gust: grid sequencing: the restart continues on level 0; the coarse levels are not run"
                           << " (steady.sequencing, sequencing_drop_orders, sequencing_monitor_tol, sequencing_hold and sequencing_cfl_start"
                           << " have no effect; the checkpoint keeps level 0's reference hold)\n";
        }
        SteadyVerdict v = driver.verdict();
        IterationPrinter log(setup.log);
        int done = 0;
        amrex::Real phys_time = 0;   // the plotfile's time: the physical time of an unsteady run
        // The unsteady run's outputs (UnsteadyOutput.H): averages, probes, and the record a checkpoint carries.
        std::unique_ptr<RunningAverages> averages;
        std::unique_ptr<Probes> probes;
        const std::string probes_out = probes_output_file(in.io.plot_file, unsteady ? ProbeAxis::time : ProbeAxis::iteration);
        UnsteadyRecord urec;
        // A checkpoint with everything a restart reads: the state, the histories, the models' fields and,
        // in an unsteady run, its record and the averages.
        auto write_full_checkpoint = [&] (const UnsteadyRecord* rec) {
            check_overwrite();
            const auto name = write_checkpoint(in.io.chk_file, state, in.mesh, &mesh.zlevels());
            write_driver_history(name, driver.history());
            closure.write_checkpoint(name);
            canopy.write_checkpoint(name);
            buildings.write_checkpoint(name);
            write_immersed_terrain_checkpoint(name, immersed_hash);
            NestFaces::write_checkpoint(name, nest_hash);
            if (newton_stage != nullptr) { newton_stage->write_checkpoint(name); }
            if (rec != nullptr) {
                write_unsteady_record(name, *rec);
                averages->write_checkpoint(name);
            }
            amrex::Print() << "gust: wrote " << name << "\n";
            final_chk = name;
        };
        // A plotfile of the state, with the averages in an unsteady run, at the time t.
        auto write_full_plotfile = [&] (amrex::Real t) {
            amrex::MultiFab avg_fields;
            std::vector<std::string> avg_names;
            if (averages) { averages->plot_fields(avg_fields, avg_names); }
            const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ, canopy_plot, lad_plot, beta_plot, wall_plot, t,
                                             averages ? &avg_fields : nullptr, averages ? &avg_names : nullptr);
            amrex::Print() << "gust: wrote " << name << "\n";
        };
        if (unsteady) {
            amrex::Print() << "gust: unsteady: physical steps to t = " << uopt.stop_time << " s, "
                           << (uopt.cfl > 0 ? "dt from unsteady.cfl " + std::to_string(double(uopt.cfl))
                                            : "dt " + std::to_string(double(uopt.dt)) + " s")
                           << "; each step to " << uopt.inner_drop_orders << " orders below its first residual (floor "
                           << "steady.residual_tol " << sopt.criterion.residual_tol << ") or " << uopt.max_inner
                           << " inner iterations\n";
            TimeDriver td(uopt, *time_term, driver, mesh, geom, sopt.criterion.residual_tol, [&] (State& s) {
                bc.fill_ghosts(s, res_in.surface, fopt.buoyancy, false, nullptr);
            });
            averages = std::make_unique<RunningAverages>(uopt.average_start, ba, dm);
            // A checkpoint of an unsteady run: its time, steps and averages go on.
            const bool resumed = restart && read_unsteady_record(in.io.restart, urec);
            if (resumed) {
                if (urec.average_start != double(uopt.average_start)) {
                    char msg[160];
                    std::snprintf(msg, sizeof(msg), "gust: unsteady.average_start = %.17g s differs from the checkpoint's %.17g s",
                                  double(uopt.average_start), urec.average_start);
                    amrex::Abort(std::string(msg) + " (io.restart = " + in.io.restart + "); the averages would mix two starts");
                }
                td.resume(urec.time, urec.steps, urec.dt, urec.anchor_time, urec.anchor_steps);
                averages->read_checkpoint(in.io.restart, urec.weight);
                amrex::Print() << "gust: unsteady: continues from t = " << urec.time << " s after " << urec.steps
                               << " physical steps (averaged over " << urec.weight << " s)\n";
            }
            if (loads) {
                // The first row is the state the run starts from (a restart that kept the file has it).
                if (!loads->open(resumed ? urec.time : -1.0)) { loads->write_row(td.time_exact(), force_now()); }
            }
            if (!uopt.probes_file.empty()) {
                probes = Probes::read(uopt.probes_file, state.z_phys_nd, geom, ProbeAxis::time);
                const bool kept = probes->open(probes_out, uopt.probes_file, resumed ? urec.time : -1.0);
                // The first row is the state the run starts from (a restart that kept the file has the
                // checkpoint's row).
                if (!kept) { probes->write_row(probes_out, td.time_exact(), state); }
                amrex::Print() << "gust: unsteady: " << probes->probes().size() << " probes to " << probes_out << "\n";
            }
            auto record_now = [&] {
                UnsteadyRecord r;
                r.time = td.time_exact();
                r.steps = td.total_steps();
                r.dt = uopt.cfl > 0 ? 0.0 : double(uopt.dt);
                r.average_start = double(uopt.average_start);
                r.weight = averages->weight();
                r.anchor_time = td.anchor_time();
                r.anchor_steps = td.anchor_steps();
                return r;
            };
            while (!td.finished()) {
                const double t0 = td.time_exact();
                const StepRecord& st = td.advance(state);
                done += st.inner;
                state.step += st.inner;
                amrex::Print() << "gust: " << format_step(st) << "\n";
                if (td.diverged()) { break; }
                averages->add(state, t0, td.time_exact());
                if (probes) { probes->write_row(probes_out, td.time_exact(), state); }
                if (loads) { loads->write_row(td.time_exact(), force_now()); }
                // Every so many steps; the last step's outputs are the run's final ones, below.
                const std::int64_t n = td.total_steps();
                if (!td.finished() && uopt.plot_interval > 0 && n % uopt.plot_interval == 0) { write_full_plotfile(td.time()); }
                if (!td.finished() && uopt.checkpoint_interval > 0 && n % uopt.checkpoint_interval == 0) {
                    const UnsteadyRecord r = record_now();
                    write_full_checkpoint(&r);
                }
            }
            urec = record_now();
            v = driver.verdict();
            v.reasons.clear();   // the steady criterion's unmet items are not a stop of the unsteady run
            if (td.diverged()) {
                v.status = SteadyStatus::Diverged;
                v.reasons.emplace_back("a residual became NaN or infinite in physical step " + std::to_string(td.total_steps()));
            }
            phys_time = td.time();
            const int nsteps = static_cast<int>(td.records().size());
            amrex::Print() << "gust: unsteady: " << nsteps << " physical steps to t = " << td.time() << " s, "
                           << td.inner_iterations() << " inner iterations ("
                           << (nsteps > 0 ? double(td.inner_iterations()) / nsteps : 0.0) << " per step), "
                           << td.unconverged_steps() << " steps stopped at unsteady.max_inner\n";
        }
        if (!unsteady && solve && !sopt.probes_file.empty()) {
            // steady.probes_file: a row per iteration of level 0, numbered by the run's iteration count
            // (the history's, which a restart continues). The first row is the state the iterations
            // start from; a restart that kept the file has the checkpoint's row.
            probes = Probes::read(sopt.probes_file, state.z_phys_nd, geom, ProbeAxis::iteration);
            const auto n0 = static_cast<double>(driver.history().records.size());
            if (!probes->open(probes_out, sopt.probes_file, restart ? n0 : -1.0)) { probes->write_row(probes_out, n0, state); }
            amrex::Print() << "gust: steady: " << probes->probes().size() << " probes to " << probes_out << ", a row per iteration\n";
        }
        // A joint restart (DomainRole::resume) joins the joint run at once: its checkpoint is a joint one.
        while (!unsteady && role.resume == nullptr && static_cast<int>(driver.history().records.size()) < sopt.max_iterations) {
            v = driver.step(state);
            ++done;
            ++state.step;
            if (probes) { probes->write_row(probes_out, static_cast<double>(driver.history().records.size()), state); }
            const bool last = (v.status != SteadyStatus::NotSteady
                               || static_cast<int>(driver.history().records.size()) >= sopt.max_iterations);
            log.print(driver, simplec_stage, last, "", newton_stage);
            if (v.status != SteadyStatus::NotSteady) { break; }
        }
        bool joint_pass = false;   // the outputs after the joint run (JointRun.H)
        auto finish = [&] () {
        const char* status = unsteady ? (v.status == SteadyStatus::Diverged ? "diverged" : "unsteady run ended")
                           : (v.status == SteadyStatus::Steady) ? "steady"
                           : (v.status == SteadyStatus::Diverged) ? "diverged" : "not steady";
        std::string reasons;
        for (const auto& r : v.reasons) { reasons += "\n  " + r; }
        amrex::Print() << "gust: " << (joint_pass ? "joint run: " + role_name + " " : std::string()) << status << " after "
                       << done << (joint_pass ? " joint" : "") << " iterations ("
                       << driver.history().records.size() << " over the run)" << reasons << "\n";
        for (const auto& acc : driver.history().accounts) {
            amrex::Print() << "gust: stage account: " << acc.iterations << " iterations, " << acc.wall_seconds
                           << " s, " << acc.au << " AU\n";
        }
        if (!levels.empty() && !joint_pass) { amrex::Print() << sequencing_report(levels, driver.history(), status, sopt.criterion.residual_ref_iterations); }
        if (building_heat && v.status != SteadyStatus::Diverged) {
            // The walls' heat exchange at the final state (its ghost cells and the law's rates at it).
            bc.fill_ghosts(state, res_in.surface, fopt.buoyancy, false, nullptr);
            building_forcing.update(state.u, state.v, state.w);
            amrex::Print() << "gust: buildings: the walls take " << double(building_theta.wall_heat(state.theta))
                           << " K m^3/s of heat from the air (negative: they heat it)\n";
        }
        if (loads && !unsteady && v.status != SteadyStatus::Diverged) {
            const auto F = force_now();
            loads->write_row(0.0, F);
            amrex::Print() << "gust: buildings: force on the buildings (per unit density) " << F[0] << " " << F[1] << " " << F[2]
                           << " m^4/s^2 to " << loads->file() << "\n";
        }
        check_overwrite();   // the step is final here, before any output
        if (in.io.write_plotfile == 1) {
            write_full_plotfile(phys_time);
            write_hierarchy();
        }
        // The checkpoint first, with the Newton stage as it is (a restart
        // continues it: the frozen k floor, the continuation); then a run that
        // ends in the Newton stage leaves it (the k floor thawed), so the
        // turbulence plotfile sees the run's own closure. A diverged state is
        // not evaluated again.
        if (in.io.write_checkpoint == 1) {
            if (unsteady && v.status == SteadyStatus::Diverged) {
                // Its state is not t^n's (the record's): a restart would not continue the run.
                amrex::Print() << "gust: no checkpoint: the unsteady run diverged\n";
            } else {
                write_full_checkpoint(unsteady ? &urec : nullptr);
            }
        }
        if (newton_stage != nullptr && driver.stage() == 1 && v.status != SteadyStatus::Diverged) { newton_stage->leave(state); }
        if (!in.io.turb_plot_file.empty()) {
            // The closure's fields and the wall state at the final state.
            bc.fill_ghosts(state, res_in.surface, fopt.buoyancy, false, nullptr);
            ClosureFields cf;
            closure.evaluate(state, dist, ground_roughness(res_in.surface), fopt.buoyancy, mesh, geom, cf);
            SurfaceState ss;
            if (surface_layer) { compute_surface_state(state, *res_in.surface, fopt.buoyancy, mesh, geom, ss); }
            ResidualFields R;
            ResidualInputs plot_in = res_in;   // the residual the run converged: with the stage's upwind blend
            if (!simplec_stage.hybrid_sigma().empty()) { plot_in.hybrid_sigma = &simplec_stage.hybrid_sigma(); }
            residual(state, plot_in, mesh, geom, R);
            if (simplec.theta_frozen) { R.theta.setVal(0); }   // theta is data, not an unknown (the stage's R_theta)
            const auto name = write_turbulence_plotfile(in.io.turb_plot_file, state, cf, dist, surface_layer ? &ss : nullptr,
                                                        copt, geom, &R,
                                                        res_in.surface != nullptr ? res_in.surface->z0_map : nullptr);
            amrex::Print() << "gust: wrote " << name << "\n";
        }
        if (v.status == SteadyStatus::Diverged) { amrex::Abort("gust: the run diverged" + reasons); }
        if (role.outcome != nullptr) { role.outcome->steady = (v.status == SteadyStatus::Steady); }
        };
        // A joint restart writes nothing before its joint iterations: its step is the checkpoint's.
        if (role.resume == nullptr) { finish(); }
        if (role.joint && v.status != SteadyStatus::Diverged) {
            // nest.feedback = defect, relax, relax_consistent or inject: this domain joins the joint run (Run.cpp, JointRun.H), which steps it in
            // memory; then its outputs are written again at its new step.
            JointMember m;
            m.name = role_name;
            m.state = &state;
            m.mesh = &mesh;
            m.geom = &geom;
            m.chk = role.resume != nullptr ? in.io.restart : final_chk;
            m.verdict = v;
            m.max_iterations = sopt.max_iterations;
            m.fill_ghosts = [&] (State& s) { bc.fill_ghosts(s, res_in.surface, fopt.buoyancy, false, nullptr); };
            const std::string label = nest != nullptr ? nest->name + " " : std::string("outer ");
            if (nest != nullptr && role.resume == nullptr && v.status != SteadyStatus::Steady) {
                // One value caps the one-way passes and the joint exchanges: a run stopped here and restarted with a
                // larger steady.max_iterations does not repeat a run that went straight through with that value.
                amrex::Print() << "gust: nest." << nest->name << ": its one-way run ended without reaching steady "
                                  "(steady.max_iterations = " << sopt.max_iterations << "); it joins the joint run as it is\n";
            }
            if (role.resume != nullptr) {
                driver.set_cfl_floor(role.resume->cfl_floor);
            } else {
                driver.floor_cfl_at_current();
            }
            m.cfl_floor = driver.cfl_floor();
            m.step = [&, label] () {
                driver.rearm_confidence();
                const SteadyVerdict sv = driver.step(state);
                ++state.step;
                log.print(driver, simplec_stage, sv.status != SteadyStatus::NotSteady, label, newton_stage);
                return sv;
            };
            if (defect || (relax && relax->consistent())) {
                // R_c(q*) for the defect source, or for relax_consistent's added term (NestRelax.H): the residual
                // without the feedback (one of the two pointers is null).
                if (defect) { m.defect = defect.get(); }
                m.residual_without_defect = [&] (const State& q, ResidualFields& R) {
                    ResidualInputs fin = res_in;
                    fin.defect = nullptr;
                    fin.relax = nullptr;
                    if (!simplec_stage.hybrid_sigma().empty()) { fin.hybrid_sigma = &simplec_stage.hybrid_sigma(); }
                    residual(q, fin, mesh, geom, R);
                    // theta is data with simplec.theta = frozen: no theta source, and no theta in the feedback change.
                    if (simplec.theta_frozen) { R.theta.setVal(0); }
                };
            }
            if (relax) { m.relax = relax.get(); }
            if (defect || relax || inject) {
                m.feedback_tol = sopt.criterion.monitor_rel_tol;
                m.feedback_window = sopt.criterion.monitor_window;
                if (defect) {
                    // A source change of residual_tol times the residual scale passes the test at any size.
                    amrex::Real sc[6];
                    residual_scales(res_in.jvp_scales, simplec_stage.norm_height(), sc);
                    const amrex::Real f = sopt.criterion.residual_tol / sopt.criterion.monitor_rel_tol;
                    m.feedback_floor = {f * sc[0], f * sc[1], f * sc[2], f * sc[4], f * sc[5]};
                } else {
                    // The relaxation's target and the injection's target: the jvp scales.
                    const JvpScales& js = res_in.jvp_scales;
                    m.feedback_floor = {js.velocity, js.velocity, js.velocity, js.theta, js.tke};
                }
            }
            if (face_flux) { m.face_flux = face_flux.get(); }
            m.probed_residual = [&] (const State& q, ResidualFields& R, FluxProbe* probe) {
                ResidualInputs fin = res_in;
                fin.defect = nullptr;
                fin.relax = nullptr;
                fin.flux_probe = probe;
                if (!simplec_stage.hybrid_sigma().empty()) { fin.hybrid_sigma = &simplec_stage.hybrid_sigma(); }
                residual(q, fin, mesh, geom, R);
            };
            if (inject) { m.inject = inject.get(); }
            if (nest != nullptr) {
                m.parent = parent.get();
                m.faces = nest_faces.get();
            }
            role.joint(m);
            if (m.iterations > 0) {
                v = m.verdict;
                done = m.iterations;
                joint_pass = true;
                finish();
            }
        }
        return final_chk;
    }

    // No stage: the run writes its state and stops.
    if (in.io.write_plotfile == 1) {
        const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ, canopy_plot, lad_plot, beta_plot, wall_plot);
        amrex::Print() << "gust: wrote " << name << "\n";
        write_hierarchy();
    }
    if (in.io.write_checkpoint == 1) {
        const auto name = write_checkpoint(in.io.chk_file, state, in.mesh, &mesh.zlevels());
        canopy.write_checkpoint(name);
        buildings.write_checkpoint(name);
        write_immersed_terrain_checkpoint(name, immersed_hash);
        NestFaces::write_checkpoint(name, nest_hash);
        amrex::Print() << "gust: wrote " << name << "\n";
        final_chk = name;
    }
    return final_chk;
}

} // namespace gust
