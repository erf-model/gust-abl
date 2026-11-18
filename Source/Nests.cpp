#include "Nests.H"

#include "Forcing.H"
#include "NestInject.H"
#include "NestParent.H"
#include "NestRelax.H"
#include "SurfaceLayer.H"

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

namespace gust {

namespace {

bool is_identifier (const std::string& s)
{
    if (s.empty() || std::isalpha(static_cast<unsigned char>(s[0])) == 0) { return false; }
    return std::all_of(s.begin(), s.end(), [] (char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_'; });
}

// The position of x in outer cells from x0, and whether it is a whole number
// of cells (the tolerance of check_nest_placement, NestParent.cpp).
bool whole_cells (amrex::Real x, amrex::Real x0, double dx, int& n)
{
    const double c = (static_cast<double>(x) - static_cast<double>(x0)) / dx;
    if (!std::isfinite(c) || std::abs(c) > 1.0e9) { return false; }
    n = static_cast<int>(std::lround(c));
    return std::abs(c - n) <= 1.0e-4;
}

std::string num (double x)
{
    std::ostringstream os;
    os << x;
    return os.str();
}

} // namespace

const char* nest_feedback_name (NestFeedback f)
{
    return f == NestFeedback::defect ? "defect" : f == NestFeedback::relax ? "relax" : f == NestFeedback::inject ? "inject" : "none";
}

NestOptions read_nest_options ()
{
    NestOptions opt;
    amrex::ParmParse pp("nest");
    if (pp.countname("names") > 1) { opt.errors.emplace_back("nest.names is given more than once"); }
    std::vector<std::string> names;
    pp.queryarr("names", names);
    if (pp.countname("feedback") > 1) { opt.errors.emplace_back("nest.feedback is given more than once"); }
    std::string fb = "none";
    if (pp.countval("feedback") > 1) { opt.errors.emplace_back("nest.feedback takes one value: none, defect, relax or inject"); }
    pp.query("feedback", fb);
    if (fb == "defect") {
        opt.feedback = NestFeedback::defect;
    } else if (fb == "relax") {
        opt.feedback = NestFeedback::relax;
    } else if (fb == "inject") {
        opt.feedback = NestFeedback::inject;
    } else if (fb != "none") {
        opt.errors.emplace_back("nest.feedback = '" + fb + "' is not none, defect, relax or inject");
    }
    if (pp.contains("feedback") && names.empty()) {
        opt.errors.emplace_back("nest.feedback = " + fb + " needs nest.names");
    }
    if (pp.countname("feedback_factor") > 1) { opt.errors.emplace_back("nest.feedback_factor is given more than once"); }
    if (pp.contains("feedback_factor")) {
        if (pp.countval("feedback_factor") != 1) {
            opt.errors.emplace_back("nest.feedback_factor takes one value");
        } else {
            pp.get("feedback_factor", opt.feedback_factor);
            const std::string e = opt.feedback == NestFeedback::inject ? check_inject_factor(opt.feedback_factor)
                                                                       : check_feedback_factor(opt.feedback_factor);
            if (!e.empty()) { opt.errors.push_back(e); }
        }
        if (opt.feedback != NestFeedback::relax && opt.feedback != NestFeedback::inject) {
            opt.errors.emplace_back("nest.feedback_factor has no effect without nest.feedback = relax or inject");
        }
    } else if (opt.feedback == NestFeedback::inject) {
        opt.feedback_factor = kInjectFactorDefault;
    }
    if (pp.countname("feedback_flux") > 1) { opt.errors.emplace_back("nest.feedback_flux is given more than once"); }
    if (pp.contains("feedback_flux")) {
        std::string ff;
        if (pp.countval("feedback_flux") != 1) {
            opt.errors.emplace_back("nest.feedback_flux takes one value: outer or nest");
        } else {
            pp.get("feedback_flux", ff);
            if (ff == "nest") {
                opt.feedback_flux = NestFeedbackFlux::nest;
            } else if (ff != "outer") {
                opt.errors.emplace_back("nest.feedback_flux = '" + ff + "' is not outer or nest");
            }
        }
        if (opt.feedback != NestFeedback::defect) {
            opt.errors.emplace_back("nest.feedback_flux has no effect without nest.feedback = defect");
        }
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        const std::string& n = names[i];
        if (!is_identifier(n) || n == "names" || n == "feedback" || n == "feedback_factor" || n == "feedback_flux") {
            opt.errors.emplace_back("nest.names: '" + n + "' is not a nest name (a letter, then letters, digits or '_'; "
                                    "not 'names', 'feedback', 'feedback_factor' or 'feedback_flux')");
            continue;
        }
        if (std::find(names.begin(), names.begin() + static_cast<std::ptrdiff_t>(i), n) != names.begin() + static_cast<std::ptrdiff_t>(i)) {
            opt.errors.emplace_back("nest.names: '" + n + "' is named more than once");
            continue;
        }
        NestSpec s;
        s.name = n;
        amrex::ParmParse pn("nest." + n);
        bool ok = true;
        for (const char* key : {"lo", "hi"}) {
            const std::string full = pn.prefixedName(key);
            if (pn.countname(key) > 1) { opt.errors.emplace_back(full + " is given more than once"); ok = false; }
            std::vector<amrex::Real> v;
            if (!pn.queryarr(key, v)) {
                opt.errors.emplace_back(full + " is missing (the nest's " + (key[0] == 'l' ? "low" : "high")
                                        + " corner: x y in m)");
                ok = false;
                continue;
            }
            if (v.size() != 2) {
                opt.errors.emplace_back(full + " takes two values (x y in m)");
                ok = false;
                continue;
            }
            auto& c = (key[0] == 'l') ? s.lo : s.hi;
            c = {v[0], v[1]};
        }
        if (ok) { opt.nests.push_back(s); }
    }
    return opt;
}

std::vector<std::string> check_nests (const NestOptions& opt, const MeshInputs& outer)
{
    std::vector<std::string> err;
    // The outer cells each nest spans, for the overlap test.
    std::vector<std::array<int, 4>> span;
    std::vector<std::string> spanned;
    const char* xy = "xy";
    for (const auto& s : opt.nests) {
        const std::string key = "nest." + s.name;
        std::array<int, 4> c{};
        bool ok = true;
        for (int d = 0; d < 2; ++d) {
            const std::string ax(1, xy[d]);
            const int N = outer.n_cell[d];
            if (N < 1 || !(outer.prob_hi[d] > outer.prob_lo[d])) { ok = false; continue; }   // validate() reports the mesh
            if (!std::isfinite(s.lo[d]) || !std::isfinite(s.hi[d])) {
                err.push_back(key + ".lo and " + key + ".hi " + ax + " must be finite");
                ok = false;
                continue;
            }
            const double dx = (static_cast<double>(outer.prob_hi[d]) - static_cast<double>(outer.prob_lo[d])) / N;
            int o = 0, e = 0;
            if (!whole_cells(s.lo[d], outer.prob_lo[d], dx, o) || !whole_cells(s.hi[d], outer.prob_lo[d], dx, e)) {
                err.push_back(key + ".lo and " + key + ".hi " + ax + " must lie on the outer cell faces (cells of " + num(dx)
                              + " m from " + num(outer.prob_lo[d]) + " m)");
                ok = false;
                continue;
            }
            if (e <= o) {
                err.push_back(key + ".hi " + ax + " must be above " + key + ".lo " + ax);
                ok = false;
                continue;
            }
            const int m = (outer.is_periodic[d] != 0) ? 0 : kNestMargin;
            if (o < m || e > N - m) {
                err.push_back(key + " must lie at least " + std::to_string(m) + " outer cells inside the outer " + ax
                              + " faces (outer cells " + std::to_string(m) + " to " + std::to_string(N - m) + "; the nest spans "
                              + std::to_string(o) + " to " + std::to_string(e) + ")");
                ok = false;
                continue;
            }
            c[static_cast<std::size_t>(2 * d)] = o;
            c[static_cast<std::size_t>(2 * d + 1)] = e;
        }
        if (!ok) { continue; }
        if (c[0] == 0 && c[1] == outer.n_cell[0] && c[2] == 0 && c[3] == outer.n_cell[1]) {
            err.push_back(key + " covers the whole outer domain in x and y (a nest needs at least one nest face)");
            continue;
        }
        for (std::size_t j = 0; j < span.size(); ++j) {
            const auto& b = span[j];
            if (c[0] < b[1] && b[0] < c[1] && c[2] < b[3] && b[2] < c[3]) {
                err.push_back(key + " overlaps nest." + spanned[j] + " (nests at one level may touch, not overlap)");
            }
        }
        span.push_back(c);
        spanned.push_back(s.name);
    }
    return err;
}

MeshInputs nest_mesh (const NestSpec& nest, const MeshInputs& outer)
{
    MeshInputs m = outer;
    for (int d = 0; d < 2; ++d) {
        const double dx = (static_cast<double>(outer.prob_hi[d]) - static_cast<double>(outer.prob_lo[d])) / outer.n_cell[d];
        int o = 0, e = 0;
        whole_cells(nest.lo[d], outer.prob_lo[d], dx, o);
        whole_cells(nest.hi[d], outer.prob_lo[d], dx, e);
        m.n_cell[d] = kNestRatio * (e - o);
        // On the outer cell faces exactly (a corner may be up to 1e-4 cells off them).
        m.prob_lo[d] = static_cast<amrex::Real>(static_cast<double>(outer.prob_lo[d]) + o * dx);
        m.prob_hi[d] = static_cast<amrex::Real>(static_cast<double>(outer.prob_lo[d]) + e * dx);
        // A nest on the outer high face (o > 0 is allowed across a periodic one) takes that face exactly.
        if (e == outer.n_cell[d]) { m.prob_hi[d] = outer.prob_hi[d]; }
        // Periodic only across the whole width of a periodic outer direction.
        m.is_periodic[d] = (outer.is_periodic[d] != 0 && o == 0 && e == outer.n_cell[d]) ? 1 : 0;
    }
    m.hierarchy_levels = 1;
    return m;
}

std::string nest_output_name (const std::string& base, const std::string& name)
{
    return base.empty() ? base : base + "_" + name;
}

void apply_nest_inputs (const NestSpec& nest, Inputs& in)
{
    in.mesh = nest_mesh(nest, in.mesh);
    in.io.plot_file = nest_output_name(in.io.plot_file, nest.name);
    in.io.chk_file = nest_output_name(in.io.chk_file, nest.name);
    in.io.turb_plot_file = nest_output_name(in.io.turb_plot_file, nest.name);
    in.io.restart.clear();
    in.terrain.given.clear();
}

void apply_nest_boundary_inputs (const std::string& parent_chk, const amrex::Array<int, 3>& periodic, BoundaryInputs& bin)
{
    for (int d = 0; d < 2; ++d) {
        const FaceType f = (periodic[d] != 0) ? FaceType::periodic : FaceType::nest;
        bin.face[amrex::Orientation(d, amrex::Orientation::low)] = f;
        bin.face[amrex::Orientation(d, amrex::Orientation::high)] = f;
    }
    bin.nest_parent = parent_chk;
    bin.sponge_width = 0;
    bin.sponge_strength = 0;
    // rayleigh_fields: u v w theta tke.
    bin.rayleigh_fields[0] = false;
    bin.rayleigh_fields[1] = false;
    bin.rayleigh_fields[4] = false;
    bin.nest_relax_width = 0;
    bin.nest_relax_factor = 1;
    bin.nest_ratio = kNestRatio;   // nests refine the outer cells 2 2 1 (nest_mesh)
    bin.given.clear();
}

void keep_outer_forcing (const ForcingOptions& fopt, const SurfaceLayerOptions& surf, OuterOutcome& out)
{
    out.forcing_set = true;
    out.ug = fopt.coriolis.ug;
    out.vg = fopt.coriolis.vg;
    out.pressure_gradient = {fopt.pressure_gradient[0], fopt.pressure_gradient[1]};
    out.surface_heat_flux = surf.heat_flux;
}

void take_outer_forcing (const OuterOutcome& out, ForcingOptions& fopt, SurfaceLayerOptions& surf)
{
    if (!out.forcing_set) { return; }
    fopt.coriolis.ug = out.ug;
    fopt.coriolis.vg = out.vg;
    fopt.pressure_gradient[0] = out.pressure_gradient[0];
    fopt.pressure_gradient[1] = out.pressure_gradient[1];
    surf.heat_flux = out.surface_heat_flux;
}

} // namespace gust
