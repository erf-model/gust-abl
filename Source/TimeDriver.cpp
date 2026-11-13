#include "TimeDriver.H"

#include "PseudoTime.H"

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace gust {

namespace {

std::string num (double x)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", x);
    return buf;
}

void abort_on (const std::vector<std::string>& err)
{
    if (err.empty()) { return; }
    std::string all = "gust: invalid inputs:";
    for (const auto& e : err) { all += "\n  " + e; }
    amrex::Abort(all);
}

} // namespace

const std::vector<std::string>& unsteady_keys ()
{
    static const std::vector<std::string> k{"stop_time",     "dt",
                                                "cfl",           "inner_drop_orders",
                                                "max_inner",     "average_start",
                                                "plot_interval", "checkpoint_interval",
                                                "probes_file"};
    return k;
}

std::vector<std::string> check_unsteady_options (const UnsteadyOptions& opt)
{
    std::vector<std::string> err;
    auto bad = [&] (amrex::Real x, const char* key, const char* range, bool ok) {
        // isnan first: gcc traps on an ordered compare with a NaN under the invalid trap
        if (std::isnan(x) || !ok) { err.emplace_back(std::string("unsteady.") + key + " = " + num(double(x)) + ": must be " + range); }
    };
    const auto given = [&] (const char* k) { return std::find(opt.given.begin(), opt.given.end(), k) != opt.given.end(); };
    if (!given("stop_time")) {
        err.emplace_back("unsteady.stop_time is required with steady.stage = unsteady");
    } else {
        bad(opt.stop_time, "stop_time", "finite and > 0 (s)", !std::isnan(opt.stop_time) && std::isfinite(opt.stop_time) && opt.stop_time > 0);
    }
    bad(opt.cfl, "cfl", "0 (use unsteady.dt) or finite and > 0", !std::isnan(opt.cfl) && std::isfinite(opt.cfl) && opt.cfl >= 0);
    const bool use_cfl = !std::isnan(opt.cfl) && opt.cfl > 0;
    if (use_cfl && given("dt")) {
        err.emplace_back("unsteady.dt and unsteady.cfl > 0 are exclusive (give one: a fixed step or a target CFL)");
    } else if (!use_cfl) {
        if (!given("dt")) {
            err.emplace_back("unsteady.dt is required unless unsteady.cfl > 0");
        } else {
            TimeTermOptions t;
            t.dt = opt.dt;
            for (auto& e : check_time_term_options(t)) { err.push_back(std::move(e)); }
        }
    }
    bad(opt.inner_drop_orders, "inner_drop_orders", "finite, > 0 and <= 15",
        !std::isnan(opt.inner_drop_orders) && opt.inner_drop_orders > 0 && opt.inner_drop_orders <= 15);
    if (opt.max_inner < 1) { err.emplace_back("unsteady.max_inner = " + std::to_string(opt.max_inner) + ": must be >= 1"); }
    const bool stop_ok = !std::isnan(opt.stop_time) && std::isfinite(opt.stop_time) && opt.stop_time > 0;
    bad(opt.average_start, "average_start", "finite, >= 0 and < unsteady.stop_time (s)",
        !std::isnan(opt.average_start) && std::isfinite(opt.average_start) && opt.average_start >= 0
            && (!stop_ok || opt.average_start < opt.stop_time));
    if (opt.plot_interval < 0) {
        err.emplace_back("unsteady.plot_interval = " + std::to_string(opt.plot_interval) + ": must be >= 0");
    }
    if (opt.checkpoint_interval < 0) {
        err.emplace_back("unsteady.checkpoint_interval = " + std::to_string(opt.checkpoint_interval) + ": must be >= 0");
    }
    if (given("probes_file") && opt.probes_file.empty()) { err.emplace_back("unsteady.probes_file must name a file"); }
    return err;
}

UnsteadyOptions read_unsteady_options ()
{
    UnsteadyOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("unsteady");
    for (const auto& key : unsteady_keys()) {
        if (!pp.contains(key.c_str())) { continue; }
        opt.given.push_back(key);
        if (pp.countname(key) > 1) { err.emplace_back("unsteady." + key + " is given more than once"); }
        if (pp.countval(key.c_str()) != 1) { err.emplace_back("unsteady." + key + " takes exactly one value"); }
    }
    abort_on(err);
    pp.query("stop_time", opt.stop_time);
    pp.query("dt", opt.dt);
    pp.query("cfl", opt.cfl);
    pp.query("inner_drop_orders", opt.inner_drop_orders);
    pp.query("max_inner", opt.max_inner);
    pp.query("average_start", opt.average_start);
    pp.query("plot_interval", opt.plot_interval);
    pp.query("checkpoint_interval", opt.checkpoint_interval);
    pp.query("probes_file", opt.probes_file);
    abort_on(check_unsteady_options(opt));
    return opt;
}

amrex::Real convective_dt (const State& q, const Mesh& mesh, const amrex::Geometry& geom, amrex::Real cfl)
{
    // local_dtau at CFL 1, without diffusion and without the floor, is 1 / max(rate, kMinRate) cell by cell;
    // its minimum (all-reduced) is the inverse of the largest rate.
    amrex::MultiFab inv_rate;
    local_dtau(q, nullptr, mesh, geom, amrex::Real(1), inv_rate, amrex::Real(0));
    return cfl * inv_rate.min(0);
}

std::string check_step_dt (amrex::Real dt, amrex::Real cfl)
{
    // isfinite first: no ordered compare with a NaN (the invalid trap)
    const char* src = cfl > 0 ? "unsteady.cfl" : "unsteady.dt";
    if (!std::isfinite(dt) || !(dt > 0)) {
        return std::string("the physical step from ") + src + " is " + num(double(dt)) + ": not finite and > 0 (a NaN or "
               "infinite state?)";
    }
    if (cfl > 0 && dt >= cfl / (2 * kMinRate)) {
        return "unsteady.cfl = " + num(double(cfl)) + " gives no step: the field is at rest (its largest convective rate is "
               "at the floor kMinRate); give unsteady.dt instead";
    }
    return {};
}

bool is_last_step (double left, double dt)
{
    return left <= dt * (1 + 1.0e-3);
}

TimeDriver::TimeDriver (const UnsteadyOptions& opt, TimeTerm& term, Driver& driver, const Mesh& mesh,
                        const amrex::Geometry& geom, amrex::Real residual_tol, GhostFill fill)
    : m_opt(opt), m_term(term), m_driver(driver), m_mesh(mesh), m_geom(geom), m_tol(residual_tol), m_fill(std::move(fill))
{
    abort_on(check_unsteady_options(opt));
    if (std::isnan(residual_tol) || !std::isfinite(residual_tol) || !(residual_tol > 0)) {
        amrex::Abort("gust: TimeDriver: steady.residual_tol must be finite and > 0");
    }
    if (!m_fill) { amrex::Abort("gust: TimeDriver: no ghost fill"); }
}

void TimeDriver::resume (double time, std::int64_t steps, double dt_ref, double anchor_time, std::int64_t anchor_steps)
{
    if (!m_records.empty()) { amrex::Abort("gust: TimeDriver::resume: after the first step"); }
    if (!std::isfinite(time) || !(time >= 0) || steps < 0) {
        amrex::Abort("gust: io.restart: the checkpoint's physical time " + num(time) + " s is not finite and >= 0");
    }
    if (!(time < double(m_opt.stop_time))) {
        amrex::Abort("gust: io.restart: the checkpoint is at t = " + num(time) + " s, not before unsteady.stop_time = "
                     + num(double(m_opt.stop_time)) + " s; give a later unsteady.stop_time");
    }
    if (!std::isfinite(anchor_time) || !(anchor_time >= 0) || !(anchor_time <= time) || anchor_steps < 0
        || anchor_steps > steps) {
        amrex::Abort("gust: io.restart: the checkpoint's time anchor (" + num(anchor_time) + " s, step "
                     + std::to_string(anchor_steps) + ") is not at or before its time and step");
    }
    m_time = time;
    m_steps_ref = steps;
    const bool same_dt = m_opt.cfl <= 0 && double(m_opt.dt) == dt_ref;
    m_anchor_t = same_dt ? anchor_time : time;
    m_anchor_n = same_dt ? anchor_steps : steps;
}

const StepRecord& TimeDriver::advance (State& q)
{
    if (finished()) { amrex::Abort("gust: TimeDriver::advance: the run has ended (unsteady.stop_time reached or diverged)"); }
    // 1. The ghost cells and dt.
    m_fill(q);
    amrex::Real dt = m_opt.cfl > 0 ? convective_dt(q, m_mesh, m_geom, m_opt.cfl) : m_opt.dt;
    const std::string bad_dt = check_step_dt(dt, m_opt.cfl);
    if (!bad_dt.empty()) { amrex::Abort("gust: TimeDriver::advance: " + bad_dt); }
    const double left = double(m_opt.stop_time) - m_time;
    // A step within 0.1 % of the time left takes all of it (no sliver of a last step).
    const bool last = is_last_step(left, double(dt));
    if (last) { dt = static_cast<amrex::Real>(left); }
    // 2. q^n and dt.
    m_term.set_dt(dt, m_opt.cfl > 0 ? "the step from unsteady.cfl" : "unsteady.dt");
    m_term.store(q);
    // 3. The inner iterations.
    StepRecord s;
    s.step = static_cast<int>(m_steps_ref + static_cast<std::int64_t>(m_records.size())) + 1;
    s.dt = dt;
    const amrex::Real drop = std::pow(amrex::Real(10), -m_opt.inner_drop_orders);
    std::vector<amrex::Real> threshold, ref;
    while (s.inner < m_opt.max_inner) {
        const SteadyVerdict v = m_driver.step(q);
        ++s.inner;
        if (v.status == SteadyStatus::Diverged) {
            s.diverged = true;
            break;
        }
        const std::vector<amrex::Real>& res = m_driver.history().records.back().residuals;
        if (s.inner == 1) {
            ref.resize(res.size());
            threshold.resize(res.size());
            for (std::size_t i = 0; i < res.size(); ++i) {
                ref[i] = amrex::max(std::abs(res[i]), m_tol);
                threshold[i] = amrex::max(drop * std::abs(res[i]), m_tol);
            }
        }
        bool met = true;
        s.inner_residual = 0;
        for (std::size_t i = 0; i < res.size(); ++i) {
            met = met && std::abs(res[i]) <= threshold[i];
            s.inner_residual = amrex::max(s.inner_residual, std::abs(res[i]) / ref[i]);
        }
        if (met) {
            s.converged = true;
            break;
        }
    }
    m_diverged = s.diverged;
    if (!s.diverged) {
        // t_* + (N - n_*) dt for a fixed step (no sum of rounded steps), else the sum in double. A
        // clipped last step (one that ends more than 0.1 % of dt from t_* + (N - n_*) dt) moves the
        // anchor to (stop_time, N), so a restart to a later stop goes on from there; a full last
        // step keeps it, so a chain of restarts at full steps keeps the times of one run.
        const std::int64_t N = total_steps() + 1;
        const double fixed = m_anchor_t + double(N - m_anchor_n) * double(m_opt.dt);
        m_time = last ? double(m_opt.stop_time) : (m_opt.cfl > 0 ? m_time + double(dt) : fixed);
        if (last && m_opt.cfl <= 0 && std::abs(m_time - fixed) > 1.0e-3 * double(m_opt.dt)) {
            m_anchor_t = m_time;
            m_anchor_n = N;
        }
        m_done = last;
    }
    s.time = static_cast<amrex::Real>(m_time);
    m_records.push_back(s);
    return m_records.back();
}

int TimeDriver::unconverged_steps () const
{
    int n = 0;
    for (const auto& s : m_records) { n += (!s.converged && !s.diverged) ? 1 : 0; }
    return n;
}

long TimeDriver::inner_iterations () const
{
    long n = 0;
    for (const auto& s : m_records) { n += s.inner; }
    return n;
}

std::string format_step (const StepRecord& s)
{
    std::string line = "step " + std::to_string(s.step) + ": t = " + num(double(s.time)) + " s, dt = " + num(double(s.dt))
                       + " s, " + std::to_string(s.inner) + " inner iterations, inner residual " + num(double(s.inner_residual));
    if (s.diverged) { line += " (diverged)"; }
    else if (!s.converged) { line += " (unsteady.max_inner reached)"; }
    return line;
}

} // namespace gust
