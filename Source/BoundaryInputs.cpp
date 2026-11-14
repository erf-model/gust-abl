// The boundary-condition inputs: the face type names, the bc.* reader and its
// range checks, and the Poisson boundary types of the faces (split from
// BoundaryConditions.cpp, audit A11).
#include "BoundaryConditions.H"

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace gust {

namespace {

constexpr int kNumTypes = 7;
const char* const kNames[kNumTypes] = {"periodic", "inflow", "outflow", "slip_wall", "no_slip_wall",
                                       "surface_layer", "geostrophic"};
const char* const kErfNames[kNumTypes] = {"Periodic", "Inflow", "Outflow", "SlipWall", "NoSlipWall",
                                          "surface_layer", "SlipWall"};
const char* const kFaceKeys[6] = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

} // namespace

const char* face_type_name (FaceType t) { return kNames[static_cast<int>(t)]; }
const char* erf_bc_name (FaceType t) { return kErfNames[static_cast<int>(t)]; }

FaceType face_type_from_name (const std::string& name, const std::string& key)
{
    for (int n = 0; n < kNumTypes; ++n) {
        if (name == kNames[n]) { return static_cast<FaceType>(n); }
    }
    std::string all;
    for (int n = 0; n < kNumTypes; ++n) { all += std::string(n ? ", " : "") + kNames[n]; }
    amrex::Abort("gust: " + key + " = " + name + ": unknown face type (one of " + all + ")");
    return FaceType::periodic;
}

const char* const kSurfaceKeys[] = {"surface_z0", "surface_z0_file", "surface_zref", "surface_stability",
                                    "surface_heat_flux", "surface_wall_k", "surface_denominator_floor",
                                    "surface_stable_functions", "surface_qnse_am", "surface_qnse_bm", "surface_qnse_ah",
                                    "surface_qnse_bh"};

void use_closure_constants (BoundaryInputs& bin, SurfaceLayerOptions& surf, amrex::Real Cmu0, amrex::Real kappa,
                            amrex::Real max_geom_lscale)
{
    surf.Cmu0 = Cmu0;
    surf.kappa = kappa;
    bin.inflow.Cmu0 = Cmu0;
    bin.inflow.kappa = kappa;
    bin.inflow.l_max = max_geom_lscale;
}

BoundaryInputs read_boundary_inputs ()
{
    BoundaryInputs in;
    amrex::ParmParse pp("bc");
    const char* keys[] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "inflow_profile", "inflow_ustar", "inflow_z0",
                          "inflow_speed", "inflow_zref", "inflow_alpha", "inflow_direction", "inflow_file", "inflow_tke_length",
                          "theta_profile_z", "theta_profile", "zhi_theta_grad", "rayleigh_zdamp",
                          "rayleigh_dampcoef", "rayleigh_fields", "rayleigh_height", "sponge_width", "sponge_strength",
                          "outflow_transport", "xlo_wall_velocity", "xhi_wall_velocity", "ylo_wall_velocity", "yhi_wall_velocity",
                          "zlo_wall_velocity", "zhi_wall_velocity", "nest_parent"};
    std::vector<std::string> errs;
    for (const char* key : keys) {
        if (!pp.contains(key)) { continue; }
        in.given.emplace_back(key);
        if (pp.countname(key) > 1) { errs.emplace_back(std::string("bc.") + key + " is given more than once"); }
        const std::string k(key);
        const bool array = k.rfind("theta_profile", 0) == 0;
        const bool three = k.size() > 14 && k.compare(k.size() - 14, 14, "_wall_velocity") == 0;
        if (k == "rayleigh_fields") {
            if (pp.countval(key) < 1) { errs.emplace_back("bc.rayleigh_fields takes at least one field name"); }
            continue;
        }
        if (three && pp.countval(key) != 3) { errs.emplace_back(std::string("bc.") + key + " takes exactly three values (u v w)"); }
        if (!array && !three && pp.countval(key) != 1) { errs.emplace_back(std::string("bc.") + key + " takes exactly one value"); }
    }
    for (const char* key : kSurfaceKeys) {
        if (pp.contains(key)) { in.given.emplace_back(key); }   // read by read_surface_layer_options
    }
    if (!errs.empty()) { amrex::Abort("gust: invalid inputs:\n  " + errs.front()); }
    for (int o = 0; o < 6; ++o) {
        std::string name;
        if (pp.query(kFaceKeys[o], name)) { in.face[o] = face_type_from_name(name, std::string("bc.") + kFaceKeys[o]); }
    }
    std::string kind = "log_law";
    pp.query("inflow_profile", kind);
    using K = InflowProfileOptions::Kind;
    if (kind == "log_law") { in.inflow.kind = K::log_law; }
    else if (kind == "power_law") { in.inflow.kind = K::power_law; }
    else if (kind == "file") { in.inflow.kind = K::file; }
    else { amrex::Abort("gust: bc.inflow_profile = " + kind + ": must be log_law, power_law or file"); }
    pp.query("inflow_ustar", in.inflow.ustar);
    pp.query("inflow_z0", in.inflow.z0);
    pp.query("inflow_speed", in.inflow.speed_ref);
    pp.query("inflow_zref", in.inflow.z_ref);
    pp.query("inflow_alpha", in.inflow.alpha);
    pp.query("inflow_direction", in.inflow.direction);
    pp.query("inflow_file", in.inflow.file);
    std::string length;
    if (pp.query("inflow_tke_length", length)) {
        if (length == "closure") { in.inflow.length = InflowProfileOptions::Length::closure; }
        else if (length == "unbounded") { in.inflow.length = InflowProfileOptions::Length::unbounded; }
        else { amrex::Abort("gust: bc.inflow_tke_length = " + length + ": must be closure or unbounded"); }
    }
    amrex::Vector<amrex::Real> tz, tt;
    pp.queryarr("theta_profile_z", tz);
    pp.queryarr("theta_profile", tt);
    if (!(tz.empty() && tt.empty())) {
        in.theta_ref.z = tz;
        in.theta_ref.theta = tt;
    }
    pp.query("zhi_theta_grad", in.zhi_theta_grad);
    for (int o = 0; o < 6; ++o) {
        amrex::Vector<amrex::Real> wv;
        if (pp.queryarr((std::string(kFaceKeys[o]) + "_wall_velocity").c_str(), wv) && wv.size() == 3) {
            for (int c = 0; c < 3; ++c) { in.wall_velocity[o][c] = wv[c]; }
        }
    }
    pp.query("rayleigh_zdamp", in.rayleigh_zdamp);
    pp.query("rayleigh_dampcoef", in.rayleigh_dampcoef);
    {
        std::string rh = "absolute";
        if (pp.query("rayleigh_height", rh)) {
            if (rh == "absolute") { in.rayleigh_above_ground = false; }
            else if (rh == "above_ground") { in.rayleigh_above_ground = true; }
            else { amrex::Abort("gust: bc.rayleigh_height = " + rh + ": must be absolute or above_ground"); }
        }
        std::string transport = "mirror";
        if (pp.query("outflow_transport", transport)) {
            if (transport == "mirror") { in.outflow_mirror_transport = true; }
            else if (transport == "from_velocity") { in.outflow_mirror_transport = false; }
            else { amrex::Abort("gust: bc.outflow_transport = " + transport + ": must be mirror or from_velocity"); }
        }
        std::vector<std::string> names;
        if (pp.queryarr("rayleigh_fields", names)) {
            in.rayleigh_fields = {false, false, false, false, false};
            for (const auto& n : names) {
                const int f = (n == "u") ? 0 : (n == "v") ? 1 : (n == "w") ? 2 : (n == "theta") ? 3 : (n == "tke") ? 4 : -1;
                if (f < 0) { amrex::Abort("gust: bc.rayleigh_fields = " + n + ": must be among u, v, w, theta, tke"); }
                in.rayleigh_fields[static_cast<std::size_t>(f)] = true;
            }
        }
    }
    pp.query("sponge_width", in.sponge_width);
    pp.query("sponge_strength", in.sponge_strength);
    if (pp.query("nest_parent", in.nest_parent) && in.nest_parent.empty()) {
        amrex::Abort("gust: bc.nest_parent: the parent checkpoint's directory must not be empty");
    }
    return in;
}

std::vector<std::string> check_boundary_inputs (const BoundaryInputs& in, const amrex::Geometry& geom,
                                                const CoriolisOptions& coriolis)
{
    return check_boundary_inputs(in, geom, coriolis, geom.ProbHi(2));
}

std::vector<std::string> check_boundary_inputs (const BoundaryInputs& in, const amrex::Geometry& geom,
                                                const CoriolisOptions& coriolis, amrex::Real ztop)
{
    std::vector<std::string> err;
    bool any_inflow = false, any_outflow = false;
    for (int o = 0; o < 6; ++o) {
        const amrex::Orientation ori(o % 3, o < 3 ? amrex::Orientation::low : amrex::Orientation::high);
        const int d = ori.coordDir();
        const FaceType t = in.face[o];
        const std::string key = std::string("bc.") + kFaceKeys[o];
        if (geom.isPeriodic(d) != (t == FaceType::periodic)) {
            err.emplace_back(key + " = " + face_type_name(t) + ": periodic exactly when mesh.is_periodic is 1 in that direction");
            continue;
        }
        if (d < 2) {
            if (t == FaceType::surface_layer || t == FaceType::geostrophic) {
                err.emplace_back(key + " = " + face_type_name(t) + ": only the ground (zlo) or the top (zhi) take this type");
            }
        } else if (ori.isLow()) {
            if (t != FaceType::surface_layer && t != FaceType::slip_wall && t != FaceType::no_slip_wall) {
                err.emplace_back(key + " = " + face_type_name(t) + ": the ground must be surface_layer, slip_wall or no_slip_wall");
            }
        } else if (t != FaceType::slip_wall && t != FaceType::no_slip_wall && t != FaceType::geostrophic) {
            err.emplace_back(key + " = " + face_type_name(t) + ": the top must be slip_wall, no_slip_wall or geostrophic");
        }
        if (t == FaceType::geostrophic && !(coriolis.enabled && coriolis.geostrophic)) {
            err.emplace_back(key + " = geostrophic needs Coriolis with geostrophic forcing (the wind it holds)");
        }
        any_inflow = any_inflow || t == FaceType::inflow;
        any_outflow = any_outflow || t == FaceType::outflow;
        const std::string wkey = key + "_wall_velocity";
        const bool wgiven = std::find(in.given.begin(), in.given.end(), std::string(kFaceKeys[o]) + "_wall_velocity") != in.given.end();
        if (wgiven && t != FaceType::no_slip_wall) {
            err.emplace_back(wkey + " has no effect without " + key + " = no_slip_wall");
        } else if (wgiven) {
            const auto& wv = in.wall_velocity[o];
            if (!(std::isfinite(wv[0]) && std::isfinite(wv[1]) && std::isfinite(wv[2]))) {
                err.emplace_back(wkey + ": every component must be finite (m/s)");
            } else if (wv[d] != 0) {
                err.emplace_back(wkey + " = " + std::to_string(wv[0]) + " " + std::to_string(wv[1]) + " "
                                 + std::to_string(wv[2]) + ": the normal component must be 0 (the wall moves in its plane)");
            }
        }
    }
    auto given = [&] (const char* k) { return std::find(in.given.begin(), in.given.end(), k) != in.given.end(); };
    const bool needs_profile = any_inflow || in.sponge_width > 0 || in.rayleigh_zdamp > 0 || in.profile_needed;
    if (needs_profile) {
        for (const auto& e : check_inflow_profile_options(in.inflow)) { err.push_back(e); }
        // The keys of the other kinds would be read and ignored.
        using K = InflowProfileOptions::Kind;
        const char* kind_name = (in.inflow.kind == K::log_law) ? "log_law" : (in.inflow.kind == K::power_law) ? "power_law" : "file";
        auto only_for = [&] (const char* k, K kind) {
            if (given(k) && in.inflow.kind != kind) {
                err.emplace_back(std::string("bc.") + k + " has no effect with bc.inflow_profile = " + kind_name);
            }
        };
        only_for("inflow_ustar", K::log_law);
        only_for("inflow_speed", K::power_law);
        only_for("inflow_zref", K::power_law);
        only_for("inflow_alpha", K::power_law);
        only_for("inflow_file", K::file);
        // A file profile gives u and v itself (InflowProfile.H).
        if (given("inflow_direction") && in.inflow.kind == K::file) {
            err.emplace_back("bc.inflow_direction has no effect with bc.inflow_profile = file");
        }
    } else {
        for (const char* k : {"inflow_profile", "inflow_ustar", "inflow_z0", "inflow_speed", "inflow_zref",
                              "inflow_alpha", "inflow_direction", "inflow_file", "inflow_tke_length"}) {
            if (given(k)) {
                err.emplace_back(std::string("bc.") + k
                                 + " has no effect without an inflow face, a damping layer or init.first_guess = profile or column");
            }
        }
    }
    for (const auto& e : check_reference_profile(in.theta_ref)) { err.push_back("bc.theta_profile: " + e); }
    if (in.face[2] != FaceType::surface_layer) {
        for (const char* k : kSurfaceKeys) {
            if (given(k)) { err.emplace_back(std::string("bc.") + k + " has no effect without bc.zlo = surface_layer"); }
        }
    }
    if (in.theta_ref.z.size() != in.theta_ref.theta.size()) {
        err.emplace_back("bc.theta_profile_z and bc.theta_profile need the same number of values");
    }
    if (!std::isfinite(in.zhi_theta_grad)) {
        err.emplace_back("bc.zhi_theta_grad = " + std::to_string(in.zhi_theta_grad) + ": must be finite (K/m)");
    } else if (given("zhi_theta_grad") && in.face[5] == FaceType::no_slip_wall) {
        err.emplace_back("bc.zhi_theta_grad has no effect with bc.zhi = no_slip_wall");
    }
    if (std::isnan(in.rayleigh_zdamp) || !(in.rayleigh_zdamp >= 0 && in.rayleigh_zdamp < ztop)) {
        char top[64];
        std::snprintf(top, sizeof(top), "%.1f", static_cast<double>(ztop));
        err.emplace_back("bc.rayleigh_zdamp = " + std::to_string(in.rayleigh_zdamp) + ": must be in [0, the top level "
                         + top + " m) (0: off)");
    }
    if (std::isnan(in.rayleigh_dampcoef) || !(in.rayleigh_dampcoef >= 0 && std::isfinite(in.rayleigh_dampcoef))) {
        err.emplace_back("bc.rayleigh_dampcoef = " + std::to_string(in.rayleigh_dampcoef) + ": must be finite and >= 0 (1/s)");
    } else if (given("rayleigh_dampcoef") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_dampcoef has no effect without bc.rayleigh_zdamp > 0");
    } else if (given("rayleigh_fields") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_fields has no effect without bc.rayleigh_zdamp > 0");
    } else if (given("rayleigh_height") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_height has no effect without bc.rayleigh_zdamp > 0");
    } else if (in.rayleigh_zdamp > 0 && !(in.rayleigh_dampcoef > 0)) {
        err.emplace_back("bc.rayleigh_zdamp > 0 needs bc.rayleigh_dampcoef > 0");
    }
    if (std::isnan(in.sponge_width) || !(in.sponge_width >= 0 && std::isfinite(in.sponge_width))) {
        err.emplace_back("bc.sponge_width = " + std::to_string(in.sponge_width) + ": must be finite and >= 0 (0: off)");
    } else if (in.sponge_width > 0 && !(any_inflow || any_outflow)) {
        err.emplace_back("bc.sponge_width > 0 needs an inflow or outflow face to damp towards");
    }
    if (std::isnan(in.sponge_strength) || !(in.sponge_strength >= 0 && std::isfinite(in.sponge_strength))) {
        err.emplace_back("bc.sponge_strength = " + std::to_string(in.sponge_strength) + ": must be finite and >= 0 (1/s)");
    } else if (given("sponge_strength") && !(in.sponge_width > 0)) {
        err.emplace_back("bc.sponge_strength has no effect without bc.sponge_width > 0");
    } else if (in.sponge_width > 0 && !(in.sponge_strength > 0)) {
        err.emplace_back("bc.sponge_width > 0 needs bc.sponge_strength > 0");
    }
    if (given("outflow_transport") && !any_outflow) {
        err.emplace_back("bc.outflow_transport has no effect without an outflow face");
    }
    return err;
}

PoissonBC poisson_bc_from_faces (const BoundaryInputs& in, const amrex::Geometry& geom)
{
    PoissonBC bc;
    for (int o = 0; o < 6; ++o) {
        const amrex::Orientation ori(o % 3, o < 3 ? amrex::Orientation::low : amrex::Orientation::high);
        if (geom.isPeriodic(ori.coordDir())) {
            bc.face[o] = PoissonFace::periodic;
        } else if (in.face[o] == FaceType::outflow) {
            bc.face[o] = PoissonFace::dirichlet;
        } else {
            bc.face[o] = PoissonFace::neumann;
        }
    }
    return bc;
}

} // namespace gust
