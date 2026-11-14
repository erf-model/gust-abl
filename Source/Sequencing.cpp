#include "Sequencing.H"

#include "CanopyTke.H"
#include "ColumnGuess.H"

#include "WallDistance.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace gust {

namespace {

//! The fewest iterations before the loose residual branch of a released
//! seeded level (steady.sequencing_hold = 0) may hold.
constexpr int kReleasedMinIterations = 3;

const char* const kKeys[] = {"sequencing", "sequencing_drop_orders", "sequencing_monitor_tol", "sequencing_hold",
                             "sequencing_cfl_start"};

std::string fmt (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", x);
    return b;
}

void abort_inputs (const std::vector<std::string>& err)
{
    std::string all = "gust: invalid inputs:";
    for (const auto& e : err) { all += "\n  " + e; }
    amrex::Abort(all);
}

//! Fill dist's ghost ring by zero gradient outside the domain (periodic
//! images first), for the Tucker distance (the height above the ground fills
//! its own ghosts).
void extend_by_zero_gradient (amrex::MultiFab& dist, const amrex::Geometry& geom)
{
    dist.FillBoundary(geom.periodicity());
    amrex::Box inside = geom.Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { if (geom.isPeriodic(d)) { inside.grow(d, 1024); } }
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const auto arr = dist.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (!inside.contains(amrex::IntVect(i, j, k))) {
                const int ic = amrex::max(inside.smallEnd(0), amrex::min(inside.bigEnd(0), i));
                const int jc = amrex::max(inside.smallEnd(1), amrex::min(inside.bigEnd(1), j));
                const int kc = amrex::max(inside.smallEnd(2), amrex::min(inside.bigEnd(2), k));
                arr(i, j, k) = arr(ic, jc, kc);
            }
        });
    }
}

ResidualInputs residual_inputs (const RunSetup& s, const Closure& closure, const amrex::MultiFab& dist,
                                const BoundaryConditions& bc, const SurfaceLayerOptions& surf)
{
    ResidualInputs r;
    r.closure = &closure;
    r.wall_distance = &dist;
    r.theta_ref = s.bin.theta_ref;
    r.momentum = s.cin.momentum;
    r.theta = s.cin.theta;
    r.hybrid = s.cin.hybrid;
    r.buoyancy = s.fopt.buoyancy;
    r.buoyancy_in_momentum = !s.simplec.theta_frozen;
    r.coriolis = s.fopt.coriolis;
    r.pressure_gradient = s.fopt.pressure_gradient;
    r.surface = s.surface_layer ? &surf : nullptr;
    r.bc = &bc;
    return r;
}

//! A finer level's driver before its first step: the coarsest level's
//! references, then (steady.sequencing_hold = 0) no reference hold and
//! (steady.sequencing_cfl_start = previous) the previous level's CFL.
void start_seeded (Driver& driver, const std::vector<amrex::Real>& ref, const SequencingOptions& seq, amrex::Real cfl_prev)
{
    driver.seed_reference(ref);
    if (seq.hold == 0) { driver.release_reference_hold(); }
    if (seq.cfl_previous) { driver.start_cfl(cfl_prev); }
}

} // namespace

// ---------------------------------------------------------------------------
SequencingOptions read_sequencing_options ()
{
    SequencingOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("steady");
    for (const char* key : kKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("steady.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { err.emplace_back(std::string("steady.") + key + " takes exactly one value"); }
    }
    if (!err.empty()) { abort_inputs(err); }
    int on = 0;
    if (pp.query("sequencing", on) && on != 0 && on != 1) {
        abort_inputs({"steady.sequencing = " + std::to_string(on) + ": must be 0 or 1"});
    }
    opt.enabled = (on == 1);
    pp.query("sequencing_drop_orders", opt.drop_orders);
    pp.query("sequencing_monitor_tol", opt.monitor_tol);
    if (pp.query("sequencing_hold", opt.hold) && opt.hold != 0 && opt.hold != 1) {
        abort_inputs({"steady.sequencing_hold = " + std::to_string(opt.hold) + ": must be 0 or 1"});
    }
    std::string start;
    if (pp.query("sequencing_cfl_start", start)) {
        if (start == "cfl0") { opt.cfl_previous = false; }
        else if (start == "previous") { opt.cfl_previous = true; }
        else { abort_inputs({"steady.sequencing_cfl_start = " + start + ": must be cfl0 or previous"}); }
    }
    return opt;
}

std::vector<std::string> check_sequencing_options (const SequencingOptions& opt, const SteadyParams& strict,
                                                   const MeshInputs& mesh)
{
    std::vector<std::string> err;
    // The ranges only matter with sequencing on: off, the defaults must not
    // refuse a deck main accepts (steady.residual_drop_orders may be below 2).
    // isfinite first: an ordered compare with a NaN traps under amrex.fpe_trap_invalid.
    if (opt.enabled) {
        if (!std::isfinite(opt.drop_orders) || !(opt.drop_orders > 0)
            || !(std::isfinite(strict.residual_drop_orders) && opt.drop_orders <= strict.residual_drop_orders)) {
            err.emplace_back("steady.sequencing_drop_orders = " + fmt(double(opt.drop_orders))
                             + ": must be > 0 and <= steady.residual_drop_orders = " + fmt(double(strict.residual_drop_orders))
                             + " (the loose stop of the coarse levels is not stricter than the strict one)");
        }
        if (!std::isfinite(opt.monitor_tol) || !(opt.monitor_tol >= 0) || !(opt.monitor_tol < 1)) {
            err.emplace_back("steady.sequencing_monitor_tol = " + fmt(double(opt.monitor_tol)) + ": must be 0 (off) or > 0 and < 1");
        }
    }
    const int levels = mesh.hierarchy_levels;
    if (opt.enabled && levels < 2) {
        err.emplace_back("steady.sequencing = 1 needs mesh.hierarchy_levels >= 2 (it is " + std::to_string(levels) + ")");
    }
    if (opt.enabled && levels >= 2 && levels <= kMaxHierarchyLevels) {
        // Divisibility is check_hierarchy_inputs'; here the coarsest size.
        const int r = hierarchy_ratio(levels - 1);
        const int cx = mesh.n_cell[0] / r, cy = mesh.n_cell[1] / r;
        if (cx < 2 || cy < 2) {
            err.emplace_back("steady.sequencing = 1 needs at least 2 cells in x and in y on the coarsest level: "
                             "mesh.n_cell x and y / 2^(mesh.hierarchy_levels - 1) = " + std::to_string(cx) + " "
                             + std::to_string(cy) + "; use fewer levels");
        }
    }
    // With the hold the first residual_ref_iterations iterations run at
    // cfl0 whatever the first one ran at.
    if (opt.enabled && opt.cfl_previous && opt.hold == 1) {
        err.emplace_back("steady.sequencing_cfl_start = previous needs steady.sequencing_hold = 0 "
                         "(the hold runs the first steady.residual_ref_iterations iterations at steady.cfl0)");
    }
    if (!opt.enabled) {
        for (const char* key : {"sequencing_drop_orders", "sequencing_monitor_tol", "sequencing_hold", "sequencing_cfl_start"}) {
            if (std::find(opt.given.begin(), opt.given.end(), key) != opt.given.end()) {
                err.emplace_back(std::string("steady.") + key + " has no effect unless steady.sequencing = 1");
            }
        }
    }
    return err;
}

// ---------------------------------------------------------------------------
LooseStop loose_stop (const SteadyHistory& h, const SteadyParams& p, const std::vector<MonitorSpec>& monitors,
                      const SequencingOptions& opt, int hold_iterations)
{
    LooseStop out;
    if (h.count == 0 || h.nonfinite) { return out; }
    // (a) every residual dropped drop_orders below its floored reference, or
    // below the tolerance; after the level's reference hold. Without the
    // hold (a released seeded level) the test runs against the seed raised
    // by the residuals so far, after kReleasedMinIterations iterations (or
    // residual_ref_iterations if fewer): one iteration at cfl0 barely moves
    // a prolonged field, and a stop there would hand it on unrelaxed. A
    // residual that rises after the stop is not seen; the threshold at the
    // stop is never above the held one (the reference is a maximum over
    // fewer values).
    const int first = hold_iterations > 0 ? hold_iterations : std::min(p.residual_ref_iterations, kReleasedMinIterations);
    if (h.count >= first && !h.res_last.empty()) {
        const amrex::Real drop = std::pow(amrex::Real(10), -opt.drop_orders);
        bool all = true;
        for (std::size_t q = 0; q < h.res_last.size(); ++q) {
            const amrex::Real ref = std::max(h.res_ref[q], p.residual_tol);
            const amrex::Real r = h.res_last[q];
            all = all && (r <= ref * drop || r < p.residual_tol);
        }
        if (all) {
            out.stop = true;
            out.reason = "every residual " + fmt(double(opt.drop_orders)) + " orders below its reference";
            return out;
        }
    }
    // (b) every monitor within monitor_tol of its latest value over a full window.
    const auto need = static_cast<std::size_t>(p.monitor_window) + 1;
    // Not before residual_ref_iterations, whatever the hold: the monitors right after a
    // prolongation say nothing yet.
    if (opt.monitor_tol > 0 && h.count >= p.residual_ref_iterations && !monitors.empty()
        && h.mon_hist.size() == monitors.size()) {
        bool all = true;
        for (std::size_t m = 0; m < monitors.size() && all; ++m) {
            const auto& w = h.mon_hist[m];
            if (w.size() < need) { all = false; break; }
            const amrex::Real last = w.back();
            amrex::Real change = 0;
            for (const auto x : w) { change = std::max(change, std::abs(x - last)); }
            all = change <= opt.monitor_tol * std::max(std::abs(last), monitors[m].scale);
        }
        if (all) {
            out.stop = true;
            out.reason = "every monitor within " + fmt(double(opt.monitor_tol)) + " over "
                         + std::to_string(p.monitor_window) + " iterations";
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
LevelSolver::LevelSolver (const RunSetup& s, const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                          const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const std::string& label,
                          bool newton, const amrex::MultiFab* z0_fine, const Canopy* canopy_in,
                          const amrex::MultiFab* canopy_fine_detJ, BuildingForcing* buildings, const BuildingTke* building_tke,
                          const BuildingTheta* building_theta, const TimeTerm* time)
    : bc(s.bin, geom, mesh, s.fopt.coriolis), closure(s.copt), surf(s.surf)
{
    if (time != nullptr && !label.empty()) { amrex::Abort("gust: LevelSolver: the physical time term is for level 0 only"); }
    bc.set_theta_ground_extension(s.simplec.theta_frozen);
    bc.set_nest_faces(s.nest);
    if (canopy_in != nullptr && canopy_fine_detJ != nullptr) {
        canopy_field.restrict_from(*canopy_in, *canopy_fine_detJ, geom, mesh, ba, dm);
        canopy = &canopy_field;
        amrex::Print() << "gust: " << label << "canopy (restricted from level 0): " << canopy_field.summary(mesh.detJ, geom)
                       << "\n";
    } else {
        canopy = canopy_in;
    }
    if (s.surface_layer && !s.surf.z0_file.empty()) {
        if (z0_fine == nullptr) {
            read_roughness_map(s.surf.z0_file, "bc.surface_z0_file", s.idw_neighbors, static_cast<double>(s.idw_power), geom,
                               ba, dm, z0_map_field);
        } else {
            average_roughness_map(*z0_fine, geom, ba, dm, z0_map_field);
        }
        surf.z0_map = &z0_map_field;
        amrex::Print() << "gust: " << label << "surface z0 map ("
                       << (z0_fine == nullptr ? "bc.surface_z0_file" : "geometric mean of level 0") << "): "
                       << roughness_summary(z0_map_field) << "\n";
    }
    if (s.copt.wall_distance == WallDistanceKind::tucker) {
        dist.define(ba, dm, 1, 1);
        const PoissonResult wd = tucker_wall_distance(mesh, geom, z_nd, dist, s.popt);
        amrex::Print() << "gust: " << label << "Tucker wall distance: " << wd.iterations << " Poisson iterations\n";
        extend_by_zero_gradient(dist, geom);
    } else {
        height_above_ground(z_nd, dist);
    }
    if (building_tke != nullptr) {
        // The closure's length sees the buildings (BuildingTke.H). The hybrid convection blend
        // keeps the ground's distance.
        dist_ground.define(ba, dm, 1, 1);
        amrex::MultiFab::Copy(dist_ground, dist, 0, 0, 1, 1);
        // The surface layer's z0, as the closure takes it (residual_inputs: s.surface_layer).
        const amrex::Long n = building_tke->limit_wall_distance(dist, ground_roughness(s.surface_layer ? &surf : nullptr));
        amrex::Print() << "gust: " << label << "buildings: the distance to the nearest building wall settled in "
                       << building_tke->distance_passes() << " passes; a building's wall sets the closure's length in "
                       << n << " cells\n";
        // The lee options fade out around buildings, whose edges fix their separation (BuildingTke.H).
        if (building_tke->zone_active() && s.copt.realizability > 0) {
            closure.set_cap_weight(&building_tke->limiter_zone());
            amrex::Print() << "gust: " << label << "buildings: turb.realizability" << (s.cin.hybrid.active() ? " and the upwind blend" : "")
                           << " fade out near buildings (buildings.limiter_zone): off in " << building_tke->zone_off_cells()
                           << " cells, partly in " << building_tke->zone_cells() - building_tke->zone_off_cells() << "\n";
        }
    }
    res_in = residual_inputs(s, closure, dist, bc, surf);
    res_in.buildings = buildings;
    res_in.building_tke = building_tke;
    res_in.building_theta = building_theta;
    if (building_tke != nullptr) { res_in.ground_distance = &dist_ground; }
    if (canopy != nullptr) {
        canopy_faces.build(*canopy, geom);
        res_in.canopy = &canopy_faces;
        const CanopyHeatOptions& ho = canopy->heat();
        if (ho.on() && !(s.simplec.theta_frozen && ho.frozen == CanopyFrozenHeat::skip)) {
            canopy_heat.build(*canopy);
            res_in.canopy_heat = &canopy_heat;
        }
        const CanopyTkeOptions& t = canopy->tke();
        if (t.on) {
            closure.set_canopy(&canopy->drag(), static_cast<amrex::Real>(t.beta_p), static_cast<amrex::Real>(t.beta_d),
                               static_cast<amrex::Real>(t.length_alpha));
        }
        if (label.empty()) {   // once, on level 0
            const BoundaryConditions::CanopyZoneCounts zones = bc.canopy_in_boundary_zones(canopy->drag(), z_nd);
            if (zones.inflow > 0) {
                amrex::Print() << "gust: WARNING: " << zones.inflow << " canopy cells lie in the first cell column of an "
                               << "inflow face: the inflow profile there does not feel the drag; keep the forest a few "
                               << "cells away from inflow faces (fade it out, as a terrain blend does)\n";
            }
            if (zones.damping > 0) {
                amrex::Print() << "gust: WARNING: " << zones.damping << " canopy cells lie inside a damping layer "
                               << "(bc.sponge_width, or bc.rayleigh_zdamp with u, v, theta or k damped): the layer "
                               << "relaxes them towards the inflow profile, which does not feel the drag\n";
            }
        }
        if (res_in.surface != nullptr && label.empty()) {   // once, on level 0
            const amrex::Long n = canopy_columns_on_rough_ground(canopy->drag(), ground_roughness(res_in.surface),
                                                                 kForestSizedZ0);
            if (n > 0) {
                amrex::Print() << "gust: WARNING: " << label << n << " columns with canopy have a ground z0 of "
                               << double(kForestSizedZ0) << " m or more (bc.surface_z0 or bc.surface_z0_file): the canopy "
                               << "already stands for the trees, so a forest-sized z0 counts the forest twice; give "
                               << "the soil and litter's z0 (0.01 to 0.1 m) under the canopy\n";
            }
        }
    }
    closure.check_wall_distance(dist, ground_roughness(res_in.surface));
    res_in.time = time;   // before the stage copies res_in
    auto stage = std::make_unique<SimplecStage>(s.simplec, s.sopt.zline_sweeps, res_in, &closure, s.popt, mesh, geom, z_nd,
                                                s.copt.tke_floor, ba, dm, line_sweep_options(s.sopt));
    m_stage = stage.get();
    std::vector<std::unique_ptr<Stage>> stages;
    stages.push_back(std::move(stage));
    if (newton) {
        auto ns = std::make_unique<NewtonStage>(s.sopt, s.newton, s.simplec, s.popt, res_in, &closure, mesh, geom, z_nd,
                                                s.copt.tke_floor, m_stage->norm_height(), ba, dm);
        m_newton = ns.get();
        stages.push_back(std::move(ns));
    }
    m_driver = std::make_unique<Driver>(s.sopt, std::move(stages), mesh, geom, &m_stage->nu());
}

IterationLog iteration_log (const std::string& name)
{
    return name == "full" ? IterationLog::full : IterationLog::table;
}

namespace {

//! simplec.theta = frozen: level l's theta set to the inflow profile's theta by
//! height above that level's ground (the stratified-column path of the first
//! guess, Initializer.H), in place of the sounding (coarsest level) or the
//! prolonged field, so every level holds the same prescribed stratification.
//! The profile is read at the unfloored height (init.z_agl_min = 0), as the
//! inflow face's ghost cells read it, and without O'Brien (only theta is kept).
//! Logs the level's frozen theta range and mean (valid cells).
void reset_frozen_theta (const MeshHierarchy& h, int l, const RunSetup& setup, const InitOptions& iopt,
                         const InflowProfile& inflow, State& q)
{
    if (!setup.simplec.theta_frozen) { return; }
    InitOptions topt = iopt;
    topt.z_agl_min = 0;
    topt.obrien = false;
    State fg;
    fg.define(h.boxArray(l), h.DistributionMap(l));
    amrex::MultiFab::Copy(fg.z_phys_nd, h.z_phys_nd(l), 0, 0, 1, amrex::min(fg.z_phys_nd.nGrowVect(), h.z_phys_nd(l).nGrowVect()));
    first_guess_state(topt, setup.bin, h.mesh(l), h.geom(l), fg, setup.fopt.coriolis, &inflow, true);
    amrex::MultiFab::Copy(q.theta, fg.theta, 0, 0, 1, amrex::min(q.theta.nGrowVect(), fg.theta.nGrowVect()));
    const amrex::Real tmin = q.theta.min(0), tmax = q.theta.max(0);
    const double tmean = static_cast<double>(q.theta.sum(0)) / static_cast<double>(h.boxArray(l).numPts());
    char line[200];
    std::snprintf(line, sizeof(line), "gust: level %d: theta frozen at the inflow profile's theta by height above the ground: "
                  "min %.15e, max %.15e, mean %.15e K", l, static_cast<double>(tmin), static_cast<double>(tmax), tmean);
    amrex::Print() << line << "\n";
}

//! The old "gust: it ..." line (io.iteration_log = full).
void print_full (const Driver& driver, const SimplecStage& stage, bool last, const std::string& label,
                 const NewtonStage* newton)
{
    const DriverHistory& h = driver.history();
    const IterationRecord& r = h.records.back();
    if (newton != nullptr && r.stage == 1) {
        // Every Newton step (few, and each costs many residual evaluations).
        std::string line = "gust: " + label + "it " + std::to_string(r.iteration + 1) + " newton r " + fmt(double(r.relative_residual))
                           + " res";
        for (const amrex::Real v : r.residuals) { line += " " + fmt(double(v)); }
        const auto& recs = newton->continuation().records();
        int attempts = 0, gmres = 0;
        for (std::size_t n = newton->first_record(); n < recs.size(); ++n) {   // this iteration's attempts
            ++attempts;
            gmres += recs[n].gmres_iterations;
        }
        if (newton->converged()) {
            line += " converged: no step";
        } else if (attempts > 0) {
            const ContinuationRecord& c = recs.back();
            line += " |R| " + fmt(double(c.residual)) + " -> " + fmt(double(c.residual_new)) + " cfl "
                    + (c.infinite ? std::string("inf") : fmt(double(c.cfl))) + " lambda " + fmt(double(c.lambda)) + " attempts "
                    + std::to_string(attempts) + " gmres " + std::to_string(gmres)
                    + (newton->last_status() == ContinuationStatus::failed ? " FAILED" : "");
        }
        amrex::Print() << line << "\n";
        return;
    }
    if (!(r.iteration < 10 || r.iteration % 10 == 9 || last)) { return; }
    std::string line = "gust: " + label + "it " + std::to_string(r.iteration + 1) + " cfl " + fmt(double(r.cfl)) + " r "
                     + fmt(double(r.relative_residual)) + " res";
    for (const amrex::Real v : r.residuals) { line += " " + fmt(double(v)); }
    line += " mon";
    for (const amrex::Real v : r.monitors) { line += " " + fmt(double(v)); }
    const SimplecIterationInfo& info = stage.last();
    line += " poisson " + std::to_string(info.pressure.poisson.iterations) + " div " + fmt(double(info.pressure.divergence_before))
          + " -> " + fmt(double(info.pressure.divergence_after)) + " clipped " + std::to_string(info.tke.clipped);
    if (driver.confidence_share() >= 0) {
        line += " confidence " + fmt(double(driver.confidence_share())) + " zmax " + fmt(double(driver.confidence_height()));
    }
    amrex::Print() << line << "\n";
}

//! The columns of a table row after the line solver, in the residuals' names.
const char* const kTableResiduals[] = {"u", "v", "w", "tke", "theta", "p"};

std::string table_cell (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "  %9.3e", x);
    return b;
}

//! text padded with spaces to width: on the left (right-aligned) or, with
//! left_align, on the right.
std::string table_text (const std::string& text, std::size_t width, bool left_align = false)
{
    const std::string pad(text.size() < width ? width - text.size() : 0, ' ');
    return left_align ? text + pad : pad + text;
}

} // namespace

void IterationPrinter::print (const Driver& driver, const SimplecStage& stage, bool last, const std::string& label,
                              const NewtonStage* newton)
{
    if (m_style == IterationLog::full) {
        print_full(driver, stage, last, label, newton);
        return;
    }
    const DriverHistory& h = driver.history();
    const IterationRecord& r = h.records.back();
    const bool newton_row = newton != nullptr && r.stage == 1;
    if (!newton_row && !(r.iteration < 10 || r.iteration % 10 == 9 || last)) { return; }

    if (m_rows % kHeaderRows == 0) {
        std::string head = label + table_text("Iter", 6) + "  " + table_text("CFL", 9) + "  "
                         + table_text("Line solver", 11, true);
        std::string rule = label + "------  ---------  -----------";
        for (const char* name : kTableResiduals) {
            head += "  " + table_text(std::string(name) + "-res", 9);
            rule += "  ---------";
        }
        head += "  " + table_text("r", 9);
        rule += "  ---------";
        amrex::Print() << head << "\n" << rule << "\n";
    }
    ++m_rows;

    std::string notes;   // a Newton step's outcome, after r
    std::string cfl = table_cell(double(r.cfl));
    if (newton_row) {
        const auto& recs = newton->continuation().records();
        int attempts = 0, gmres = 0;
        for (std::size_t n = newton->first_record(); n < recs.size(); ++n) {   // this iteration's attempts
            ++attempts;
            gmres += recs[n].gmres_iterations;
        }
        cfl = "  " + table_text("-", 9);   // the continuation step's CFL, "-" without a step
        if (newton->converged()) {
            notes = "  converged: no step";
        } else if (attempts > 0) {
            const ContinuationRecord& c = recs.back();
            cfl = c.infinite ? "  " + table_text("inf", 9) : table_cell(double(c.cfl));
            notes = "  |R| " + fmt(double(c.residual)) + " -> " + fmt(double(c.residual_new)) + " lambda "
                    + fmt(double(c.lambda)) + " attempts " + std::to_string(attempts) + " gmres " + std::to_string(gmres)
                    + (newton->last_status() == ContinuationStatus::failed ? " FAILED" : "");
        }
    }
    const char* solver = newton_row ? "newton" : sweep_name(stage.line_sweep().kind);
    std::string row = label + table_text(std::to_string(r.iteration + 1), 6) + cfl + "  " + table_text(solver, 11, true);
    const std::vector<std::string>& names = stage.residual_names();
    for (const char* name : kTableResiduals) {
        const auto it = std::find(names.begin(), names.end(), name);
        const auto n = static_cast<std::size_t>(it - names.begin());
        row += (it != names.end() && n < r.residuals.size()) ? table_cell(double(r.residuals[n])) : "  " + table_text("-", 9);
    }
    if (!newton_row && driver.confidence_share() >= 0) {   // steady.confidence
        notes += "  confidence " + fmt(double(driver.confidence_share())) + " zmax " + fmt(double(driver.confidence_height()));
    }
    row += table_cell(double(r.relative_residual)) + notes;
    amrex::Print() << row << "\n";
}

// ---------------------------------------------------------------------------
std::vector<LevelSummary> sequence_coarse_levels (const MeshHierarchy& h, const RunSetup& setup, const SequencingOptions& seq,
                                                  const InitOptions& iopt, const InflowProfile* column,
                                                  LevelSolver& finest, State& q)
{
    const int L = h.num_levels();
    if (L < 2) { amrex::Abort("gust: sequence_coarse_levels needs mesh.hierarchy_levels >= 2"); }
    std::vector<LevelSummary> out;
    std::unique_ptr<State> coarse;   // the state of level l + 1, once solved
    std::vector<amrex::Real> ref;    // the coarsest level's residual references
    amrex::Real cfl_prev = 0;        // the CFL the last level would have run its next iteration with
    for (int l = L - 1; l >= 1; --l) {
        const amrex::Box dom = h.geom(l).Domain();
        LevelSummary s;
        s.level = l;
        s.cells_x = dom.length(0);
        s.cells_y = dom.length(1);
        const std::string label = "level " + std::to_string(l) + ": ";
        LevelSolver ls(setup, h.mesh(l), h.geom(l), h.z_phys_nd(l), h.boxArray(l), h.DistributionMap(l), label, false,
                       finest.surf.z0_map, finest.canopy, finest.canopy != nullptr ? &h.mesh(0).detJ : nullptr);
        auto ql = std::make_unique<State>();
        ql->define(h.boxArray(l), h.DistributionMap(l));
        if (!coarse) {
            // The coarsest level: the initializer on its own mesh.
            amrex::MultiFab::Copy(ql->z_phys_nd, h.z_phys_nd(l), 0, 0, 1,
                                  amrex::min(ql->z_phys_nd.nGrowVect(), h.z_phys_nd(l).nGrowVect()));
            const InitResult r = initialize_state(iopt, setup.bin, setup.popt, h.mesh(l), h.geom(l), *ql, nullptr,
                                                  setup.fopt.coriolis, column, column_sets_theta(iopt, setup.surf));
            print_initializer(r, iopt, setup.popt, label);
            reset_frozen_theta(h, l, setup, iopt, ls.bc.profile(), *ql);
        } else {
            const auto poisson = h.make_poisson(l, setup.popt);
            s.prolong = prolong_state(h, l + 1, *coarse, l, *ql, ls.bc, ls.surface(), setup.fopt.buoyancy, *poisson,
                                      setup.copt.tke_floor);
            s.prolonged = true;
            reset_frozen_theta(h, l, setup, iopt, ls.bc.profile(), *ql);
        }
        coarse.reset();
        Driver& driver = ls.driver();
        IterationPrinter log(setup.log);   // a header per level
        if (!ref.empty()) { start_seeded(driver, ref, seq, cfl_prev); }
        const SteadyParams& strict = setup.sopt.criterion;
        while (static_cast<int>(driver.history().records.size()) < setup.sopt.max_iterations) {
            const SteadyVerdict v = driver.step(*ql);
            if (v.status == SteadyStatus::Diverged) {
                log.print(driver, ls.stage(), true, label);
                std::string reasons;
                for (const auto& r : v.reasons) { reasons += "\n  " + r; }
                amrex::Abort("gust: grid sequencing: level " + std::to_string(l) + " diverged" + reasons);
            }
            if (v.status == SteadyStatus::Steady) { s.stop = "steady"; }
            else {
                const LooseStop stop = loose_stop(driver.history().criterion, strict, ls.stage().monitor_specs(), seq,
                                                  driver.hold_iterations());
                if (stop.stop) { s.stop = stop.reason; }
            }
            const bool last = !s.stop.empty()
                              || static_cast<int>(driver.history().records.size()) >= setup.sopt.max_iterations;
            log.print(driver, ls.stage(), last, label);
            if (!s.stop.empty()) { break; }
        }
        if (s.stop.empty()) { s.stop = "steady.max_iterations"; }
        const DriverHistory& dh = driver.history();
        s.iterations = static_cast<int>(dh.records.size());
        s.wall_seconds = dh.accounts.front().wall_seconds;
        s.au = dh.accounts.front().au;
        s.relative_residual = driver.relative_residual();
        if (s.stop == "steady.max_iterations") {
            amrex::Print() << "WARNING: " << label << "stopped at steady.max_iterations = " << setup.sopt.max_iterations
                           << " before the loose stop, r " << fmt(double(s.relative_residual))
                           << "; its state is prolonged as it is\n";
        } else {
            amrex::Print() << "gust: " << label << "loose stop after " << s.iterations << " iterations (" << s.stop
                           << "), r " << fmt(double(s.relative_residual)) << "\n";
        }
        out.push_back(s);
        if (ref.empty()) { ref = dh.criterion.res_ref; }
        cfl_prev = driver.cfl();
        coarse = std::move(ql);
    }
    // Level 0: the prolongation; the run's driver does the iterations.
    LevelSummary s0;
    const amrex::Box dom = h.geom(0).Domain();
    s0.cells_x = dom.length(0);
    s0.cells_y = dom.length(1);
    {
        PoissonSolver poisson(h.mesh(0), h.geom(0), h.z_phys_nd(0), setup.popt);
        s0.prolong = prolong_state(h, 1, *coarse, 0, q, finest.bc, finest.surface(), setup.fopt.buoyancy, poisson,
                                   setup.copt.tke_floor);
    }
    reset_frozen_theta(h, 0, setup, iopt, finest.bc.profile(), q);
    s0.prolonged = true;
    start_seeded(finest.driver(), ref, seq, cfl_prev);
    out.push_back(s0);
    return out;
}

std::string sequencing_report (const std::vector<LevelSummary>& levels, const DriverHistory& finest,
                               const std::string& finest_status, int ref_iterations)
{
    std::string os = "gust: grid sequencing summary\n"
                     "  level   cells (x y)  iterations   wall (s)         AU  stop\n";
    char line[256];
    for (const auto& s : levels) {
        int its = s.iterations;
        double wall = s.wall_seconds, au = s.au;
        std::string stop = s.stop;
        if (s.level == 0) {
            its = static_cast<int>(finest.records.size());
            wall = finest.accounts.empty() ? 0.0 : finest.accounts.front().wall_seconds;
            au = finest.accounts.empty() ? 0.0 : finest.accounts.front().au;
            stop = "strict: " + finest_status;
        }
        std::snprintf(line, sizeof(line), "  %5d  %5ld %5ld  %10d  %9.3f  %9.3e  %s\n", s.level, s.cells_x, s.cells_y, its,
                      wall, au, stop.c_str());
        os += line;
        if (s.prolonged) {
            std::snprintf(line, sizeof(line), "         prolonged into level %d: divergence %.3e -> %.3e, %lld k values clipped\n",
                          s.level, static_cast<double>(s.prolong.projection.divergence_before),
                          static_cast<double>(s.prolong.projection.divergence_after),
                          static_cast<long long>(s.prolong.tke_clipped));
            os += line;
        }
    }
    // The seeded references (the coarsest level's), and the largest of level
    // 0's first ref_iterations residuals from its prolonged start, which a
    // level-0 reference of its own would have been. The second is not the
    // reference of a run without sequencing (that run starts from the
    // initializer on level 0); compare with such a run's checkpoint.
    const std::vector<amrex::Real>& seed = finest.criterion.res_ref;
    std::vector<amrex::Real> own(seed.size(), amrex::Real(0));
    for (std::size_t n = 0; n < finest.records.size() && static_cast<int>(n) < ref_iterations; ++n) {
        for (std::size_t q = 0; q < own.size() && q < finest.records[n].residuals.size(); ++q) {
            own[q] = std::max(own[q], std::abs(finest.records[n].residuals[q]));
        }
    }
    os += "  level-0 residual references (seeded from the coarsest level) / largest of level 0's first "
          + std::to_string(ref_iterations) + " residuals after the prolongation:";
    for (std::size_t q = 0; q < seed.size(); ++q) { os += " " + fmt(double(seed[q])) + " / " + fmt(double(own[q])); }
    os += "\n";
    return os;
}

void print_initializer (const InitResult& r, const InitOptions& iopt, const PoissonOptions& popt, const std::string& label)
{
    char line[256];
    if (iopt.projection) {
        std::snprintf(line, sizeof(line),
                      "initializer: %d Poisson iterations (%s, reltol %.3e), divergence %.9e before, %.9e after\n",
                      r.poisson.iterations, r.poisson_backend.c_str(), static_cast<double>(popt.reltol),
                      static_cast<double>(r.divergence_before), static_cast<double>(r.divergence_after));
    } else {
        std::snprintf(line, sizeof(line), "initializer: no projection (init.projection = 0), divergence %.9e\n",
                      static_cast<double>(r.divergence_before));
    }
    amrex::Print() << "gust: " << label << line;
    if (r.backend_fallback) {
        amrex::Print() << "gust: " << label << "initializer: init.alpha_v / init.alpha_h = " << iopt.alpha_v / iopt.alpha_h
                       << " needs the split z-face coefficient, which only solver.poisson_backend = gmres_fft has;"
                       << " the projection used gmres_fft, the SIMPLEC stage keeps " << popt.backend << "\n";
    }
    if (iopt.obrien) {
        std::snprintf(line, sizeof(line), "initializer: O'Brien removed a top residual of %.6e m/s\n",
                      static_cast<double>(r.obrien_top_residual));
        amrex::Print() << "gust: " << label << line;
    }
}

} // namespace gust
