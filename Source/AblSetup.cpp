#include "AblSetup.H"

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

namespace gust {

namespace {

using Drive = AblSetupOptions::Drive;
using Theta = AblSetupOptions::Theta;

// The length the couette set-up gives turb.max_geom_lscale: far above any
// domain, so the closure's length is the inflow's kappa (z + z0).
constexpr amrex::Real kUncappedLength = amrex::Real(1.0e6);

bool uniform (const amrex::Vector<amrex::Real>& v)
{
    return std::all_of(v.begin(), v.end(), [&] (amrex::Real x) { return x == v.front(); });
}

bool non_decreasing (const amrex::Vector<amrex::Real>& v)
{
    return std::is_sorted(v.begin(), v.end());
}

bool any_face (const BoundaryInputs& bc, FaceType type)
{
    return std::find(bc.face.begin(), bc.face.end(), type) != bc.face.end();
}

// The theta of the air: the inflow table's with a file profile (bc.theta_profile is then the reference only).
const amrex::Vector<amrex::Real>& air_theta (const AblSetupFacts& f)
{
    return f.bc.inflow.kind == InflowProfileOptions::Kind::file ? f.inflow_table_theta : f.bc.theta_ref.theta;
}

// The Rayleigh layer is on and damps the field (bc.rayleigh_fields: u v w theta tke).
bool rayleigh_damps (const BoundaryInputs& bc, int field)
{
    return bc.rayleigh_zdamp > 0 && bc.rayleigh_dampcoef > 0 && bc.rayleigh_fields[static_cast<std::size_t>(field)];
}

void check_drive (Drive drive, const AblSetupFacts& f, std::vector<std::string>& err)
{
    const std::string pre = std::string("abl.drive = ") + abl_drive_name(drive) + " ";
    const auto& cor = f.forcing.coriolis;
    const auto& pg = f.forcing.pressure_gradient;
    const bool horizontal_pg = pg[0] != 0 || pg[1] != 0;
    const bool inflow = any_face(f.bc, FaceType::inflow);
    const FaceType top = f.bc.face[amrex::Orientation(2, amrex::Orientation::high)];

    const bool rotating = drive == Drive::geostrophic_column || drive == Drive::geostrophic_inflow;
    if (rotating) {
        if (!(cor.enabled && cor.geostrophic)) {
            err.emplace_back(pre + "needs forcing.coriolis = 1 and forcing.geostrophic = 1");
        }
        // The column's target wind adjusts the forcing of a column without inflow faces (ColumnGuess.H).
        const bool target = drive == Drive::geostrophic_column
                            && (f.init.column_target_u != 0 || f.init.column_target_v != 0);
        if (cor.ug == 0 && cor.vg == 0 && !target) {
            err.emplace_back(pre + (drive == Drive::geostrophic_column
                                    ? "needs a wind: forcing.ug, forcing.vg or init.column_target_u, init.column_target_v"
                                    : "needs a wind: forcing.ug, forcing.vg"));
        }
    } else if (cor.enabled) {
        err.emplace_back(pre + "refuses forcing.coriolis = 1 (use abl.drive = geostrophic_column or geostrophic_inflow)");
    }
    if (drive != Drive::pressure_driven && (horizontal_pg || pg[2] != 0)) {
        err.emplace_back(pre + "refuses a non-zero forcing.pressure_gradient (use abl.drive = pressure_driven)");
    }
    if ((drive == Drive::inflow_only || drive == Drive::couette || drive == Drive::geostrophic_inflow) && !inflow) {
        err.emplace_back(pre + "needs an inflow face (bc.xlo ... bc.yhi = inflow)");
    }

    if (drive == Drive::couette) {
        const auto& lid = f.bc.wall_velocity[amrex::Orientation(2, amrex::Orientation::high)];
        if (top != FaceType::no_slip_wall || (lid[0] == 0 && lid[1] == 0)) {
            err.emplace_back(pre + "needs bc.zhi = no_slip_wall with a moving lid (bc.zhi_wall_velocity)");
        }
        if (f.bc.inflow.kind != InflowProfileOptions::Kind::log_law) {
            err.emplace_back(pre + "needs bc.inflow_profile = log_law");
        }
        if (f.bc.inflow.length != InflowProfileOptions::Length::unbounded) {
            err.emplace_back(pre + "needs bc.inflow_tke_length = unbounded");
        }
        if (f.max_geom_lscale < kUncappedLength) {
            err.emplace_back(pre + "needs turb.max_geom_lscale >= 1e6 (the inflow's length is not capped)");
        }
        if (!f.surface_layer) {
            err.emplace_back(pre + "needs bc.zlo = surface_layer");
        } else if (!f.surface.z0_file.empty() || std::abs(f.bc.inflow.z0 - f.surface.z0) > amrex::Real(1.0e-5) * f.surface.z0) {
            err.emplace_back(pre + "needs bc.inflow_z0 equal to a uniform bc.surface_z0 (the stress stays constant)");
        }
        // The lid moves along the inflow wind (bc.inflow_direction, degrees counter-clockwise from x).
        const double lid_deg = std::atan2(double(lid[1]), double(lid[0])) * 45.0 / std::atan(1.0);
        const double turn = std::remainder(lid_deg - double(f.bc.inflow.direction), 360.0);
        if (std::abs(turn) > 0.1) {
            err.emplace_back(pre + "needs bc.zhi_wall_velocity along bc.inflow_direction (within 0.1 degrees)");
        }
    } else if (drive == Drive::pressure_driven) {
        if (!horizontal_pg) { err.emplace_back(pre + "needs a non-zero horizontal forcing.pressure_gradient"); }
        if (top != FaceType::slip_wall) {
            err.emplace_back(pre + "needs bc.zhi = slip_wall (the force balances the wall stress alone)");
        }
        // Rayleigh damping of u, v relaxes the wind toward the inflow profile. A uniform theta needs none, and it
        // would break the balance; above an inversion the layer holds the wind, which otherwise has no steady state.
        const bool holds_wind = rayleigh_damps(f.bc, 0) && rayleigh_damps(f.bc, 1);
        const bool damps_wind = rayleigh_damps(f.bc, 0) || rayleigh_damps(f.bc, 1);
        if (uniform(air_theta(f)) && damps_wind) {
            err.emplace_back(pre + "with uniform theta refuses a Rayleigh layer on u or v (bc.rayleigh_fields): it would "
                             "take momentum from the balance of force and wall stress");
        } else if (!uniform(air_theta(f)) && !holds_wind) {
            err.emplace_back(pre + "with a stratified theta needs a Rayleigh layer on u and v (bc.rayleigh_fields): above an "
                             "inversion the wind has no steady state");
        }
    } else if (drive == Drive::geostrophic_column) {
        if (f.init.first_guess != InitOptions::FirstGuess::column && !f.restart) {
            err.emplace_back(pre + "needs init.first_guess = column or io.restart");
        }
        if (inflow && f.bc.inflow.kind != InflowProfileOptions::Kind::file) {
            err.emplace_back(pre + "with an inflow face needs bc.inflow_profile = file (the column's table, "
                             "io.column_guess_file of a periodic run)");
        }
    }
}

void check_theta (Theta theta, Drive drive, const AblSetupFacts& f, std::vector<std::string>& err)
{
    const std::string pre = std::string("abl.theta = ") + abl_theta_name(theta) + " ";
    const auto& sounding = air_theta(f);
    const bool flux_surface = f.surface.stability == SurfaceLayerOptions::Stability::heat_flux;
    const bool file_inflow = f.bc.inflow.kind == InflowProfileOptions::Kind::file;
    // No table to read: the boundary checks report the missing bc.inflow_file.
    const bool known = !(file_inflow && f.inflow_table_theta.empty());

    if (theta != Theta::neutral && f.forcing.buoyancy.gravity == 0) {
        err.emplace_back(pre + "refuses forcing.gravity = 0 (theta would not act on the flow)");
    }
    if (theta != Theta::frozen && f.simplec.theta_frozen) {
        err.emplace_back(pre + "refuses simplec.theta = frozen (use abl.theta = frozen)");
    }
    if (theta != Theta::surface_flux && flux_surface) {
        err.emplace_back(pre + "refuses bc.surface_stability = heat_flux (use abl.theta = surface_flux)");
    }

    using CS = InitOptions::ColumnStability;
    const bool column = f.init.first_guess == InitOptions::FirstGuess::column;
    if (theta != Theta::surface_flux && column && (f.init.column_stability != CS::mol || f.init.column_mol != 0)) {
        err.emplace_back(pre + "refuses a stratified column first guess (init.column_stability, init.column_mol; use "
                         "abl.theta = surface_flux)");
    }
    if (theta == Theta::neutral) {
        if (!uniform(sounding) || !uniform(f.bc.theta_ref.theta)) {
            err.emplace_back(pre + "needs a uniform bc.theta_profile and, with bc.inflow_profile = file, a uniform table");
        }
        if (f.bc.zhi_theta_grad != 0) { err.emplace_back(pre + "refuses a non-zero bc.zhi_theta_grad"); }
        if (f.canopy_heat) { err.emplace_back(pre + "refuses canopy.heat (the leaves would heat the air)"); }
    } else if (theta == Theta::capped) {
        if (known && (uniform(sounding) || !non_decreasing(sounding))) {
            err.emplace_back(pre + "needs a bc.theta_profile (or bc.inflow_file table) that rises with height and never falls");
        }
        if (!(f.bc.zhi_theta_grad > 0)) { err.emplace_back(pre + "needs bc.zhi_theta_grad > 0"); }
        if (!rayleigh_damps(f.bc, 2) || !rayleigh_damps(f.bc, 3)) {
            err.emplace_back(pre + "needs a Rayleigh layer on w and theta (bc.rayleigh_zdamp, bc.rayleigh_dampcoef, "
                             "bc.rayleigh_fields): the heat the top gradient lets in needs a sink");
        }
    } else if (theta == Theta::frozen) {
        if (!f.simplec.theta_frozen) { err.emplace_back(pre + "needs simplec.theta = frozen"); }
        if (known && uniform(sounding)) {
            err.emplace_back(pre + "needs a stratified bc.theta_profile or bc.inflow_file table (uniform theta: use "
                             "abl.theta = neutral)");
        }
        if (drive == Drive::geostrophic_column && !file_inflow) {
            err.emplace_back(pre + "with abl.drive = geostrophic_column needs bc.inflow_profile = file (the theta to hold)");
        }
    } else {
        if (!f.surface_layer || !flux_surface) {
            err.emplace_back(pre + "needs bc.zlo = surface_layer and bc.surface_stability = heat_flux");
        }
        // The column with a Monin-Obukhov length hands the 3-D run its flux.
        const bool column_flux = column && f.init.column_stability == CS::mol && f.init.column_mol != 0;
        const auto& keys = f.bc.given;
        const bool flux_given = std::find(keys.begin(), keys.end(), "surface_heat_flux") != keys.end();
        if (!column_flux && !(flux_given && f.surface.heat_flux != 0)) {
            err.emplace_back(pre + "needs a non-zero bc.surface_heat_flux or init.column_mol");
        }
        // A marching column stops at init.column_time (ColumnGuess.H): its heat budget need not close.
        const bool marching = column && f.init.column_time > 0;
        if (!any_face(f.bc, FaceType::outflow) && !rayleigh_damps(f.bc, 3) && !marching) {
            err.emplace_back(pre + "needs a heat sink: an outflow face, a Rayleigh layer on theta (bc.rayleigh_zdamp, "
                             "bc.rayleigh_dampcoef, bc.rayleigh_fields) or init.column_time > 0");
        }
    }
}

} // namespace

const char* abl_drive_name (AblSetupOptions::Drive d)
{
    switch (d) {
        case Drive::inflow_only:        return "inflow_only";
        case Drive::couette:            return "couette";
        case Drive::pressure_driven:    return "pressure_driven";
        case Drive::geostrophic_column: return "geostrophic_column";
        case Drive::geostrophic_inflow: return "geostrophic_inflow";
        case Drive::unset:              break;
    }
    return "unset";
}

const char* abl_theta_name (AblSetupOptions::Theta t)
{
    switch (t) {
        case Theta::neutral:      return "neutral";
        case Theta::capped:       return "capped";
        case Theta::frozen:       return "frozen";
        case Theta::surface_flux: return "surface_flux";
        case Theta::unset:        break;
    }
    return "unset";
}

AblSetupOptions read_abl_setup_options ()
{
    AblSetupOptions opt;
    amrex::ParmParse pp("abl");
    for (const char* key : {"drive", "theta"}) {
        if (pp.countname(key) > 1) { opt.errors.emplace_back(std::string("abl.") + key + " is given more than once"); }
        if (pp.contains(key) && pp.countval(key) != 1) { opt.errors.emplace_back(std::string("abl.") + key + " takes exactly 1 value"); }
    }
    std::string name;
    if (pp.query("drive", name)) {
        bool found = false;
        for (Drive d : {Drive::inflow_only, Drive::couette, Drive::pressure_driven, Drive::geostrophic_column,
                        Drive::geostrophic_inflow}) {
            if (name == abl_drive_name(d)) { opt.drive = d; found = true; }
        }
        if (!found) {
            opt.errors.emplace_back("abl.drive = " + name + ": must be inflow_only, couette, pressure_driven, "
                                    "geostrophic_column or geostrophic_inflow");
        }
    }
    if (pp.query("theta", name)) {
        bool found = false;
        for (Theta t : {Theta::neutral, Theta::capped, Theta::frozen, Theta::surface_flux}) {
            if (name == abl_theta_name(t)) { opt.theta = t; found = true; }
        }
        if (!found) { opt.errors.emplace_back("abl.theta = " + name + ": must be neutral, capped, frozen or surface_flux"); }
    }
    // A deck's value against the name, reported before the readers (which would answer with the keys the name adds).
    auto refuse = [&] (const char* ns, const char* key, const char* bad, const std::string& name) {
        std::string v;
        if (amrex::ParmParse(ns).query(key, v) && v == bad) {
            opt.errors.emplace_back(name + " refuses " + ns + "." + key + " = " + bad);
        }
    };
    if (opt.drive == Drive::geostrophic_column || opt.drive == Drive::geostrophic_inflow) {
        const std::string name = std::string("abl.drive = ") + abl_drive_name(opt.drive);
        refuse("forcing", "coriolis", "0", name);
        refuse("forcing", "geostrophic", "0", name);
    }
    if (opt.theta == Theta::frozen) { refuse("simplec", "theta", "solve", "abl.theta = frozen"); }
    if (opt.theta == Theta::surface_flux) { refuse("bc", "surface_stability", "neutral", "abl.theta = surface_flux"); }
    return opt;
}

void apply_abl_setup_defaults (const AblSetupOptions& opt)
{
    std::string stage = "none";
    amrex::ParmParse("steady").query("stage", stage);
    if (stage == "none" || !opt.errors.empty()) { return; }

    amrex::ParmParse forcing("forcing"), bc("bc"), turb("turb"), simplec("simplec");
    if (opt.drive == Drive::geostrophic_column || opt.drive == Drive::geostrophic_inflow) {
        if (!forcing.contains("coriolis")) { forcing.add("coriolis", 1); }
        if (!forcing.contains("geostrophic")) { forcing.add("geostrophic", 1); }
    }
    if (opt.drive == Drive::couette) {
        if (!bc.contains("inflow_tke_length")) { bc.add("inflow_tke_length", std::string("unbounded")); }
        if (!turb.contains("max_geom_lscale")) { turb.add("max_geom_lscale", kUncappedLength); }
    }
    if (opt.theta == Theta::capped && !bc.contains("rayleigh_fields")) {
        bc.addarr("rayleigh_fields", std::vector<std::string>{"w", "theta"});
    }
    if (opt.theta == Theta::frozen && !simplec.contains("theta")) { simplec.add("theta", std::string("frozen")); }
    if (opt.theta == Theta::surface_flux && !bc.contains("surface_stability")) {
        bc.add("surface_stability", std::string("heat_flux"));
    }
}

std::vector<std::string> check_abl_setup (const AblSetupOptions& opt, const AblSetupFacts& facts)
{
    // After a bad name the defaults were not set, so the other checks would report what they would have set.
    if (!opt.errors.empty()) { return opt.errors; }
    std::vector<std::string> err;
    if (!opt.active()) { return err; }
    if (!facts.solve) {
        if (opt.drive != Drive::unset) { err.emplace_back("abl.drive has no effect without steady.stage = simplec"); }
        if (opt.theta != Theta::unset) { err.emplace_back("abl.theta has no effect without steady.stage = simplec"); }
        return err;
    }
    if (opt.drive != Drive::unset) { check_drive(opt.drive, facts, err); }
    if (opt.drive == Drive::couette && opt.theta != Theta::neutral) {
        err.emplace_back("abl.drive = couette needs abl.theta = neutral (the moving lid's constant-stress layer is neutral)");
    }
    if (opt.theta != Theta::unset) { check_theta(opt.theta, opt.drive, facts, err); }
    return err;
}

} // namespace gust
