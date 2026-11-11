#include "Residual.H"

#include "BuildingForcing.H"
#include "CanopyDrag.H"
#include "CanopyHeat.H"

#include "Diffusion.H"
#include "ImplicitOperator.H"
#include "BoundaryConditions.H"
#include "MassFlux.H"
#include "Reductions.H"

#include <AMReX_iMultiFab.H>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace gust {

void ResidualFields::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    using amrex::IntVect;
    u.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    v.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    w.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    p.define(ba, dm, 1, 0);
    theta.define(ba, dm, 1, 0);
    tke.define(ba, dm, 1, 0);
}

amrex::Vector<ResidualFields::Field> ResidualFields::fields ()
{
    return {{"u", &u}, {"v", &v}, {"w", &w}, {"p", &p}, {"theta", &theta}, {"tke", &tke}};
}

amrex::Vector<const amrex::MultiFab*> ResidualFields::fields () const
{
    return {&u, &v, &w, &p, &theta, &tke};
}

void ResidualWorkspace::use_layout (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    if (m_ba == ba && m_dm == dm) { return; }
    *this = ResidualWorkspace();
    m_ba = ba;
    m_dm = dm;
}

namespace {

void check_diffusivities (const char* op, const ResidualInputs& in)
{
    std::string err;
    if (in.closure != nullptr) {
        if (in.nu != nullptr || in.gamma_theta != nullptr || in.gamma_tke != nullptr) {
            err = "with a closure, nu, gamma_theta and gamma_tke must be null (the closure computes them)";
        } else if (in.wall_distance == nullptr) {
            err = "a closure needs wall_distance";
        }
    } else if (in.nu == nullptr) { err = "nu is not set"; }
    else if (in.gamma_theta == nullptr) { err = "gamma_theta is not set"; }
    else if (in.gamma_tke == nullptr) { err = "gamma_tke is not set"; }
    if (!err.empty()) { amrex::Abort(std::string("gust: ") + op + ": ResidualInputs: " + err); }
}

//! The face average in use.
FaceMean face_mean_of (const ResidualInputs& in)
{
    return (in.closure != nullptr) ? in.closure->options().face_mean : in.face_mean;
}

//! Evaluate the closure at q into f (the surface layer's z0 or map, or 0).
void evaluate_closure (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
                       ClosureFields& f, ClosureScratch* scratch = nullptr)
{
    in.closure->evaluate(q, *in.wall_distance, ground_roughness(in.surface), in.buoyancy, mesh, geom, f, scratch);
}

#ifdef GUST_DEBUG_CHECKS
//! Abort unless the caller's wall state is the one computed from q, bit for
//! bit (ghost columns included).
void require_surface_state_of (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
                               const SurfaceState& given)
{
    SurfaceState fresh;
    compute_surface_state(q, *in.surface, in.buoyancy, mesh, geom, fresh);
    const char* names[] = {"ustar", "speed", "zref", "d0", "inv_L", "tau", "drag", "kwall", "area"};
    const amrex::MultiFab* a[] = {&fresh.ustar, &fresh.speed, &fresh.zref, &fresh.d0, &fresh.inv_L,
                                  &fresh.tau, &fresh.drag, &fresh.kwall, &fresh.area};
    const amrex::MultiFab* b[] = {&given.ustar, &given.speed, &given.zref, &given.d0, &given.inv_L,
                                  &given.tau, &given.drag, &given.kwall, &given.area};
    for (int n = 0; n < 9; ++n) {
        if (b[n]->boxArray() != a[n]->boxArray() || b[n]->DistributionMap() != a[n]->DistributionMap()) {
            amrex::Abort(std::string("gust: residual: the given surface state is not on q's layout (") + names[n] + ")");
        }
        const amrex::Long diff = count_bit_differences(*a[n], *b[n], amrex::IntVect(1, 1, 0));
        if (diff > 0) {
            amrex::Abort(std::string("gust: residual: the given surface state is not q's: ") + names[n] + " differs at "
                         + std::to_string(diff) + " points");
        }
    }
}
#endif

void check_inputs (const char* op, const ResidualInputs& in)
{
    check_diffusivities(op, in);
    std::string err;
    std::vector<std::string> all;
    for (const auto& e : check_reference_profile(in.theta_ref)) { all.push_back(e); }
    const auto& js = in.jvp_scales;
    for (amrex::Real x : {js.velocity, js.theta, js.tke, js.pressure, js.step}) {
        if (std::isnan(x) || !(x > 0 && std::isfinite(x))) {
            all.emplace_back("jvp_scales: every scale and the step must be finite and > 0");
            break;
        }
    }
    for (const auto& e : check_convection_options(in.momentum)) { all.push_back("momentum " + e); }
    for (const auto& e : check_convection_options(in.theta)) { all.push_back("theta " + e); }
    for (const auto& e : check_buoyancy_options(in.buoyancy)) { all.push_back(e); }
    for (const auto& e : check_coriolis_options(in.coriolis)) { all.push_back(e); }
    if (in.surface != nullptr) {
        for (const auto& e : check_surface_layer_options(*in.surface)) { all.push_back(e); }
    }
    if (in.bc != nullptr) {
        const bool ground_wall = (in.bc->inputs().face[2] == FaceType::surface_layer);
        if (ground_wall != (in.surface != nullptr)) {
            all.emplace_back("bc.zlo = surface_layer needs ResidualInputs::surface, and other ground types must not have it");
        }
    }
    if (err.empty() && !all.empty()) { err = all.front(); }
    if (!err.empty()) { amrex::Abort(std::string("gust: ") + op + ": ResidualInputs: " + err); }
}

// Over the owned points of q (shared faces once): sum of |(q - offset) / s|,
// sum of (v / s)^2 and the number of points, per box, written to slots
// (3 per box: box_sums[3 * box + c]). One thread per (i, j) column adds its
// k values in order and the host adds the columns in storage order, so a box
// gives the same bits whatever rank owns it (as NewtonOperator::dotProduct).
template <bool OFFSET>
void scaled_sums (const amrex::MultiFab& q, const amrex::MultiFab& v, const amrex::MultiFab* offset, amrex::Real s,
                  const amrex::Periodicity& per, amrex::Vector<double>& box_sums)
{
    const auto mask = q.OwnerMask(per);
    const double sinv = 1.0 / static_cast<double>(s);
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);
        amrex::Box cols = bx;
        cols.setBig(2, klo);
        const auto npts = static_cast<std::size_t>(cols.numPts());
        amrex::Gpu::DeviceVector<double> colsum(3 * npts);
        double* cs = colsum.data();
        const auto lo = amrex::lbound(cols);
        const auto len = amrex::length(cols);
        const auto qa = q.const_array(mfi);
        const auto va = v.const_array(mfi);
        const auto m = mask->const_array(mfi);
        const auto o = OFFSET ? offset->const_array(mfi) : qa;
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            double s0 = 0, s1 = 0, s2 = 0;
            for (int k = klo; k <= khi; ++k) {
                if (m(i, j, k) == 0) { continue; }
                const double qq = OFFSET ? static_cast<double>(qa(i, j, k) - o(i, j, k)) : static_cast<double>(qa(i, j, k));
                const double vv = static_cast<double>(va(i, j, k)) * sinv;
                s0 += std::abs(qq) * sinv;
                s1 += vv * vv;
                s2 += 1.0;
            }
            const std::size_t c = static_cast<std::size_t>((i - lo.x) + (j - lo.y) * len.x + (k0 - lo.z) * len.x * len.y);
            cs[3 * c] = s0;
            cs[3 * c + 1] = s1;
            cs[3 * c + 2] = s2;
        });
        amrex::Vector<double> h(3 * npts);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, colsum.begin(), colsum.end(), h.begin());
        amrex::Gpu::streamSynchronize();
        double t[3] = {0, 0, 0};
        for (std::size_t c = 0; c < npts; ++c) {
            for (int n = 0; n < 3; ++n) { t[n] += h[3 * c + n]; }
        }
        for (int n = 0; n < 3; ++n) { box_sums[3 * mfi.index() + n] = t[n]; }
    }
}

amrex::Vector<const amrex::MultiFab*> solution_fields (const State& q)
{
    return {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
}

} // namespace

void residual (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
               ResidualFields& R, ClosureFields* fields, const SurfaceState* surface, ResidualWorkspace* ws)
{
    check_inputs("residual", in);
    if (R.u.empty()) { R.define(q.p.boxArray(), q.p.DistributionMap()); }
    if (R.p.boxArray() != q.p.boxArray() || R.p.DistributionMap() != q.p.DistributionMap()) {
        amrex::Abort("gust: residual: R is not on q's layout");
    }
    if (surface != nullptr && in.surface == nullptr) {
        amrex::Abort("gust: residual: a surface state is given without a wall model (ResidualInputs::surface)");
    }
    const amrex::MultiFab& z = q.z_phys_nd;
    ResidualWorkspace local;
    ResidualWorkspace& w = (ws != nullptr) ? *ws : local;
    w.use_layout(q.p.boxArray(), q.p.DistributionMap());
    // The diffusivities: given, or the closure's at q.
    ClosureFields* cf = (fields != nullptr) ? fields : &w.cf;
    const amrex::MultiFab* nu = in.nu;
    const amrex::MultiFab* gamma_theta = in.gamma_theta;
    const amrex::MultiFab* gamma_tke = in.gamma_tke;
    if (in.closure != nullptr) {
        evaluate_closure(q, in, mesh, geom, *cf, &w.closure);
        nu = &cf->nu;
        gamma_theta = &cf->gamma_theta;
        gamma_tke = &cf->gamma_tke;
    }
    const FaceMean mean = face_mean_of(in);
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& F = w.F;
    transport_fluxes(q.u, q.v, q.w, z, mesh, geom, F);
    if (in.bc != nullptr) { in.bc->outflow_transport(F); }
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};

    // Momentum: C - D + G p - S.
    const amrex::Array<const amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    const amrex::Array<amrex::MultiFab*, 3> Rm{&R.u, &R.v, &R.w};
    amrex::Array<amrex::MultiFab, 3>& tmp = w.tmp;
    for (int c = 0; c < 3; ++c) {
        auto& Gc = w.Gc[static_cast<int>(velocity_loc(c))];
        convection(velocity_loc(c), *vel[c], Fp, in.momentum, mesh, geom, *Rm[c], Gc);
        if (tmp[c].empty()) { tmp[c].define(Rm[c]->boxArray(), Rm[c]->DistributionMap(), 1, 0); }
        if (in.hybrid_sigma != nullptr) {
            ConvectionOptions up = in.momentum;
            up.limiter = Limiter::upwind;
            convection(velocity_loc(c), *vel[c], Fp, up, mesh, geom, tmp[c], Gc);
            blend_towards_upwind(static_cast<int>(velocity_loc(c)), *in.hybrid_sigma, tmp[c], *Rm[c]);
        }
    }
    // Diffusion as fluxes, so the wall model can set the ground fluxes.
    auto& Gm = w.Gm;
    auto& Gt = w.Gt;
    auto& Gk = w.Gk;
    momentum_diffusive_fluxes(vel, *nu, z, mesh, geom, Gm, mean);
    scalar_diffusive_fluxes(q.theta, *gamma_theta, z, mesh, geom, Gt, mean);
    scalar_diffusive_fluxes(q.tke, *gamma_tke, z, mesh, geom, Gk, mean);
    // The wall state: the caller's (from its ghost fill of q), or computed
    // here. A state that read the top w face is computed again: the fill
    // writes that face after the wall model (SurfaceState::reads_top_face).
    const SurfaceState* ssp = (surface != nullptr && !surface->reads_top_face) ? surface : nullptr;
    if (in.surface != nullptr) {
        if (ssp == nullptr) {
            compute_surface_state(q, *in.surface, in.buoyancy, mesh, geom, w.ss);
            ssp = &w.ss;
        }
#ifdef GUST_DEBUG_CHECKS
        else { require_surface_state_of(q, in, mesh, geom, *ssp); }
#endif
        set_ground_fluxes(*ssp, *in.surface, Gm[0][2], Gm[1][2], Gt[2], Gk[2]);
    }
    // The walls' stress replaces the viscous flux through them (BuildingForcing.H).
    if (in.buildings != nullptr) {
        in.buildings->close_viscous_fluxes({&mesh.open_fraction(0), &mesh.open_fraction(1), &mesh.open_fraction(2)}, Gm);
    }
    for (int c = 0; c < 3; ++c) {
        flux_divergence(velocity_loc(c), Gm[c], mesh, geom, tmp[c]);
        amrex::MultiFab::Subtract(*Rm[c], tmp[c], 0, 0, 1, 0);
    }
    pressure_gradient(q.p, z, mesh, geom, {&tmp[0], &tmp[1], &tmp[2]});
    for (int c = 0; c < 3; ++c) { amrex::MultiFab::Add(*Rm[c], tmp[c], 0, 0, 1, 0); }
    amrex::MultiFab& theta_ref = w.theta_ref;
    if (theta_ref.empty()) { theta_ref.define(q.p.boxArray(), q.p.DistributionMap(), 1, 1); }
    fill_theta_ref(in.theta_ref, z, theta_ref);
    if (in.buoyancy_in_momentum) { subtract_buoyancy(q.theta, theta_ref, in.buoyancy, mesh, R.w); }
    subtract_coriolis(q.u, q.v, q.w, in.coriolis, mesh, R.u, R.v, R.w);
    subtract_pressure_gradient(in.pressure_gradient, mesh, R.u, R.v, R.w);

    // Continuity.
    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, R.p, w.omega);

    // Scalars: C - D.
    amrex::MultiFab& d = w.d;
    if (d.empty()) { d.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0); }
    auto& Gcc = w.Gc[static_cast<int>(Loc::cell)];
    convection(Loc::cell, q.theta, Fp, in.theta, mesh, geom, R.theta, Gcc);
    if (in.hybrid_sigma != nullptr) {
        ConvectionOptions up = in.theta;
        up.limiter = Limiter::upwind;
        convection(Loc::cell, q.theta, Fp, up, mesh, geom, d, Gcc);
        blend_towards_upwind(0, *in.hybrid_sigma, d, R.theta);
    }
    flux_divergence(Loc::cell, Gt, mesh, geom, d);
    amrex::MultiFab::Subtract(R.theta, d, 0, 0, 1, 0);
    ConvectionOptions kopt;
    kopt.limiter = Limiter::upwind;
    kopt.bounded = true;
    convection(Loc::cell, q.tke, Fp, kopt, mesh, geom, R.tke, Gcc);
    flux_divergence(Loc::cell, Gk, mesh, geom, d);
    amrex::MultiFab::Subtract(R.tke, d, 0, 0, 1, 0);
    if (in.closure != nullptr) {
        // R_k -= S_k = P_s + P_b - eps (the wall-k rows below replace the first cell's).
        amrex::MultiFab::Subtract(R.tke, cf->Ps, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(R.tke, cf->Pb, 0, 0, 1, 0);
        amrex::MultiFab::Add(R.tke, cf->eps, 0, 0, 1, 0);
        if (cf->canopy) {
            // The canopy's wake production and short-circuit sink (CanopyTke.H).
            amrex::MultiFab::Subtract(R.tke, cf->Pc, 0, 0, 1, 0);
            amrex::MultiFab::Add(R.tke, cf->Dc, 0, 0, 1, 0);
        }
    }
    if (in.surface != nullptr) { apply_wall_k(*ssp, *in.surface, q.tke, R.tke); }
    // The canopy drag before the boundary rows, which replace their rows' residuals.
    if (in.canopy != nullptr) { add_canopy_drag(*in.canopy, q.u, q.v, q.w, R.u, R.v, R.w); }
    // The leaves' heat exchange, also before the boundary rows (CanopyHeat.H).
    if (in.canopy_heat != nullptr && in.canopy == nullptr) {
        amrex::Abort("gust: residual: the leaf heat (canopy_heat) needs the canopy drag (canopy) as well");
    }
    if (in.canopy_heat != nullptr) { add_canopy_heat(*in.canopy_heat, q, R.theta); }
    if (in.buildings != nullptr) {
        in.buildings->update(q.u, q.v, q.w);
        in.buildings->add_sink(q.u, q.v, q.w, R.u, R.v, R.w);
    }
    if (in.bc != nullptr) {
        in.bc->add_damping(q, R);
        in.bc->apply_rows(q, R);
    }
    if (in.buildings != nullptr) { in.buildings->apply_rows(q.u, q.v, q.w, R.u, R.v, R.w); }

    // A face shared by two boxes (or a periodic image) gets one value: the two
    // copies can differ in the last bits where a compiler contracts
    // multiply-adds or the mesh's shared nodes differ, and a per-box update
    // would let them drift apart.
    for (amrex::MultiFab* mf : {&R.u, &R.v, &R.w}) { mf->OverrideSync(geom.periodicity()); }
}

void line_jacobian (const State& q, const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh,
                    const amrex::Geometry& geom, amrex::MultiFab& tri)
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, mesh, geom, F);
    if (in.bc != nullptr) { in.bc->outflow_transport(F); }
    SurfaceState ss;
    const bool wall = (in.surface != nullptr && eq != Equation::w);
    if (wall) { compute_surface_state(q, *in.surface, in.buoyancy, mesh, geom, ss); }
    ClosureFields cf;
    if (in.closure != nullptr) { evaluate_closure(q, in, mesh, geom, cf); }
    CanopyRows rows;
    if (in.canopy != nullptr) { canopy_rows(*in.canopy, q.u, q.v, q.w, rows); }
    if (in.canopy_heat != nullptr) { canopy_heat_rows(*in.canopy_heat, q.u, q.v, q.w, rows); }
    line_jacobian({&F[0], &F[1], &F[2]}, q.z_phys_nd, in, eq, dir, mesh, geom, tri, wall ? &ss : nullptr,
                  in.closure != nullptr ? &cf : nullptr, (in.canopy != nullptr || in.canopy_heat != nullptr) ? &rows : nullptr);
}

void implicit_rows (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& A, const SurfaceState* surface, const ClosureFields* fields, const CanopyRows* canopy)
{
    check_diffusivities("implicit_rows", in);
    const bool closure = (in.closure != nullptr);
    if (closure && (fields == nullptr || fields->empty())) {
        amrex::Abort("gust: implicit_rows: with a closure the rows need its fields at the state (Closure::evaluate)");
    }
    Loc loc = Loc::cell;
    const amrex::MultiFab* g = closure ? &fields->nu : in.nu;
    switch (eq) {
    case Equation::u: loc = Loc::xface; break;
    case Equation::v: loc = Loc::yface; break;
    case Equation::w: loc = Loc::zface; break;
    case Equation::theta: g = closure ? &fields->gamma_theta : in.gamma_theta; break;
    case Equation::tke: g = closure ? &fields->gamma_tke : in.gamma_tke; break;
    }
    implicit_coefficients(loc, F, *g, z_nd, mesh, geom, A, face_mean_of(in));
    if (closure && eq == Equation::tke) {
        // Point-implicit sources: the Patankar diagonal (eps + max(-P_b, 0))/k.
        amrex::MultiFab::Add(A, fields->dSdk, 0, coef_P, 1, 0);
    }
    const bool wall = (in.surface != nullptr && eq != Equation::w);
    if (wall) {
        if (surface == nullptr) {
            amrex::Abort("gust: implicit_rows: with a wall model the first cells' rows need the surface state");
        }
        apply_wall_to_implicit(*surface, loc, mesh, geom, A);
        if (eq == Equation::tke) { apply_wall_k_to_implicit(*surface, *in.surface, A); }
    }
    const bool momentum = (eq == Equation::u || eq == Equation::v || eq == Equation::w);
    if (in.canopy != nullptr && momentum) {
        if (canopy == nullptr) {
            amrex::Abort("gust: implicit_rows: with a canopy the momentum rows need the drag's diagonal (canopy_rows)");
        }
        const int c = (eq == Equation::u) ? 0 : (eq == Equation::v) ? 1 : 2;
        const amrex::MultiFab& ap = canopy->ap[static_cast<std::size_t>(c)];
        if (ap.boxArray() != A.boxArray() || ap.DistributionMap() != A.DistributionMap()) {
            amrex::Abort("gust: implicit_rows: the canopy rows are not on the equation's layout");
        }
        amrex::MultiFab::Add(A, ap, 0, coef_P, 1, 0);
    }
    if (in.canopy_heat != nullptr && eq == Equation::theta) {
        if (canopy == nullptr || canopy->theta.empty()) {
            amrex::Abort("gust: implicit_rows: with the leaf heat the theta rows need its diagonal (canopy_heat_rows)");
        }
        if (canopy->theta.boxArray() != A.boxArray() || canopy->theta.DistributionMap() != A.DistributionMap()) {
            amrex::Abort("gust: implicit_rows: the canopy heat rows are not on the equation's layout");
        }
        amrex::MultiFab::Add(A, canopy->theta, 0, coef_P, 1, 0);
    }
    // Last: the Dirichlet rows become identity rows.
    if (in.bc != nullptr) { in.bc->apply_to_implicit(eq, loc, z_nd, A); }
    if (in.buildings != nullptr && momentum) {
        in.buildings->add_to_implicit((eq == Equation::u) ? 0 : (eq == Equation::v) ? 1 : 2,
                                      {&mesh.open_fraction(0), &mesh.open_fraction(1), &mesh.open_fraction(2)}, A);
    }
}

void line_jacobian (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& tri, const SurfaceState* surface, const ClosureFields* fields, const CanopyRows* canopy)
{
    if (dir < 0 || dir > 2) { amrex::Abort("gust: line_jacobian: dir must be 0, 1 or 2"); }
    amrex::MultiFab A;
    implicit_rows(F, z_nd, in, eq, mesh, geom, A, surface, fields, canopy);
    const Loc loc = (eq == Equation::u) ? Loc::xface : (eq == Equation::v) ? Loc::yface : (eq == Equation::w) ? Loc::zface : Loc::cell;
    if (tri.empty()) { tri.define(A.boxArray(), A.DistributionMap(), n_line_coefs, 0); }
    check_operator_field("line_jacobian", "tri", tri, loc, mesh, 0);
    if (tri.nComp() != n_line_coefs) { amrex::Abort("gust: line_jacobian: tri needs 3 components"); }
    amrex::MultiFab::Copy(tri, A, coef_XM + 2 * dir, line_lower, 1, 0);
    amrex::MultiFab::Copy(tri, A, coef_P, line_diag, 1, 0);
    amrex::MultiFab::Copy(tri, A, coef_XP + 2 * dir, line_upper, 1, 0);
    tri.mult(amrex::Real(-1), line_lower, 1, 0);
    tri.mult(amrex::Real(-1), line_upper, 1, 0);
}

amrex::Real jvp (const State& q, const State& v, const ResidualInputs& in, const Mesh& mesh,
                 const amrex::Geometry& geom, ResidualFields& Jv, const ResidualFields* R0, bool refill,
                 ResidualWorkspace* ws)
{
    if (refill && in.bc == nullptr) {
        amrex::Abort("gust: jvp: refill needs ResidualInputs::bc (the ghost fill of q + eps v)");
    }
    if (Jv.u.empty()) { Jv.define(q.p.boxArray(), q.p.DistributionMap()); }
    const auto qf = solution_fields(q);
    const auto vf = solution_fields(v);
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        if (vf[n]->boxArray() != qf[n]->boxArray() || vf[n]->DistributionMap() != qf[n]->DistributionMap()
            || vf[n]->nGrow() < qf[n]->nGrow()) {
            amrex::Abort("gust: jvp: v must be a state on q's layout with q's ghost cells");
        }
    }
    check_inputs("jvp", in);
    ResidualWorkspace local;
    ResidualWorkspace& w = (ws != nullptr) ? *ws : local;
    w.use_layout(q.p.boxArray(), q.p.DistributionMap());
    const amrex::Periodicity per = geom.periodicity();
    amrex::MultiFab& theta_ref = w.jvp_theta_ref;
    if (theta_ref.empty()) { theta_ref.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0); }
    fill_theta_ref(in.theta_ref, q.z_phys_nd, theta_ref);
    const auto& js = in.jvp_scales;
    const amrex::Real scale[6] = {js.velocity, js.velocity, js.velocity, js.pressure, js.theta, js.tke};
    // One slot per (field, box, sum), all-reduced (one contributor per slot:
    // exact) and added in field and box order: the step, and so J v, has the
    // same bits for any DistributionMapping of the BoxArray (the GMRES of
    // the Newton stage then gives the same iterates on any rank count).
    const int nboxes = static_cast<int>(q.p.boxArray().size());
    amrex::Vector<double> slots(static_cast<std::size_t>(6 * 3 * nboxes), 0.0);
    for (int n = 0; n < 6; ++n) {
        amrex::Vector<double> box_sums(static_cast<std::size_t>(3 * nboxes), 0.0);
        if (n == 4) {
            scaled_sums<true>(*qf[n], *vf[n], &theta_ref, scale[n], per, box_sums);
        } else {
            scaled_sums<false>(*qf[n], *vf[n], nullptr, scale[n], per, box_sums);
        }
        std::copy(box_sums.begin(), box_sums.end(), slots.begin() + static_cast<std::ptrdiff_t>(3 * nboxes * n));
    }
    amrex::ParallelDescriptor::ReduceRealSum(slots.dataPtr(), static_cast<int>(slots.size()));
    double sums[3] = {0, 0, 0};
    for (int n = 0; n < 6; ++n) {
        for (int b = 0; b < nboxes; ++b) {
            for (int c = 0; c < 3; ++c) { sums[c] += slots[static_cast<std::size_t>(3 * (nboxes * n + b) + c)]; }
        }
    }
    if (sums[1] == 0) {
        for (auto& f : Jv.fields()) { f.mf->setVal(0.0); }
        return 0;
    }
    const double mean_q = sums[0] / sums[2];
    const double rms_v = std::sqrt(sums[1] / sums[2]);
    const double eps = static_cast<double>(js.step) * (1 + mean_q) / rms_v;

    if (R0 == nullptr) {
        residual(q, in, mesh, geom, w.base, nullptr, nullptr, &w);
        R0 = &w.base;
    }
    State& qp = w.qp;
    if (qp.p.empty()) { qp.define(q.p.boxArray(), q.p.DistributionMap()); }
    amrex::MultiFab::Copy(qp.z_phys_nd, q.z_phys_nd, 0, 0, 1, q.z_phys_nd.nGrow());
    const amrex::Vector<amrex::MultiFab*> pf{&qp.u, &qp.v, &qp.w, &qp.p, &qp.theta, &qp.tke};
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        amrex::MultiFab::LinComb(*pf[n], 1, *qf[n], 0, static_cast<amrex::Real>(eps), *vf[n], 0, 0, 1, qf[n]->nGrow());
    }
    // The refill's wall state goes to the residual (which computes its own when it read the top w face).
    const SurfaceState* ss = nullptr;
    if (refill) {
        in.bc->fill_ghosts(qp, in.surface, in.buoyancy, false, in.surface != nullptr ? &w.ss : nullptr);
        if (in.surface != nullptr) { ss = &w.ss; }
    }
    residual(qp, in, mesh, geom, Jv, nullptr, ss, &w);
    auto jf = Jv.fields();
    const auto rf = R0->fields();
    for (int n = 0; n < static_cast<int>(jf.size()); ++n) {
        amrex::MultiFab::Subtract(*jf[n].mf, *rf[n], 0, 0, 1, 0);
        jf[n].mf->mult(static_cast<amrex::Real>(1 / eps), 0, 1, 0);
    }
    return static_cast<amrex::Real>(eps);
}

} // namespace gust
