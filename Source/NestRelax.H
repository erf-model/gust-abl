#ifndef GUST_NEST_RELAX_H_
#define GUST_NEST_RELAX_H_

// The patches' feedback on the parent by relaxation over their footprints
// (nest.feedback = relax; plan rev. 6 tracker, PR-70's fallback PR-71). The
// joint run (JointRun.H) is the same as with nest.feedback = defect, with this
// source in place of the defect source (NestDefect.H). Each exchange:
//
//   up: q* = the parent state with I q_f on the nests' footprints
//       (restrict_patch, NestDefect.H), and on the parent's points
//
//     R_c(q_c) + c (q_c - q*) = 0,    c = w(s) / tau,    1 / tau = f U / dx,
//
// for u, v, w, theta and k, with
//   R_c    the parent's residual (Residual.H);
//   q_c    the parent's state; q* its target, fixed between exchanges;
//   f      nest.feedback_factor (dimensionless, > 0; default 1);
//   U      the mean over the footprint's cells of the horizontal speed of
//          I q_f, sqrt(u_c^2 + v_c^2) with u_c, v_c the means of the two
//          faces of the cell (one value per nest, set at each exchange);
//   dx     the parent's smaller horizontal cell size, min(dx, dy);
//   s      the point's distance from the nearest lateral face of the
//          footprint, in parent cells (i - L + 1/2 for cell i of a footprint
//          from L to H, i - L for a face point along its normal; from the
//          high face H + 1 - i - 1/2 and H + 1 - i), over the faces the nest
//          has (none across a periodic direction the nest spans);
//   w(s)   min(1, max(0, s / M)), M = kDefectMargin (2 parent cells): 0 on
//          the footprint's faces and outside it, rising linearly to 1 at M
//          cells in, 1 over the rest of the footprint (all z).
// At f U / dx the relaxation's rate matches the x-advection part U / dx of a
// parent cell's diagonal (which also has y, vertical and diffusive terms;
// near the ground diffusion dominates it): f = 1 pulls the parent towards
// the patch about as strongly as its own advection holds it, f >> 1 makes it
// follow the patch.
// The parent under the patch is a balance between its own equations and the
// relaxation: q_c - q* is O(R_c(q_c) / c), not 0 at a finite factor (with
// the defect source, q_c = I q_f on the inner box at the fixed point). On the
// Nests_Two_Way hill it falls from about f^-0.3 at small f to f^-0.6 at
// f >= 1, because q* (the nest) moves with f too.
//
// No pressure target: the defect source carries the nest's pressure into the
// parent through R_c(q*); this source does not, so the parent's pressure under
// the patch is its own, and as f grows the parent's velocity there stops
// responding to it (SIMPLEC's d = 1/dtau + c also shrinks the pressure
// correction there). The flow rate through the footprint is then set by a
// loop whose gain tends to 1 as f grows (the nest returns the flow rate it
// is given, the relaxation imposes it), so the exchanges grow with f. On
// the hill the nest moves away from the fine run as f grows from 0.25, and
// the outer domain does so above f = 1
// (docs/boundary_conditions.rst, Two-way nests by relaxation).
//
// Implicit: c is on the diagonal of each equation (implicit_rows) and does
// not depend on q_c, so the diagonal is the source's exact derivative. No
// row is fixed and the pressure correction is not masked: continuity holds
// on every parent cell, so under the patch the pressure correction acts on
// the relaxed velocity, which approaches the divergence-free I q_f only as f
// grows. The rows a boundary condition or the wall law replaces replace the
// source too: it is added right after the time term and the defect source
// (Residual.cpp), and its diagonal after the time term's (implicit_rows),
// before them.
//
// The kernels have the form of PR-66's relaxation zone (NestFaces.H: R +=
// c (q - q_parent), A_P += c); they are kept here, not shared, so the two
// features do not share code paths.

#include "NestDefect.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace gust {

class NestRelax
{
public:
    //! What set() found for one nest.
    struct Nest
    {
        std::string name;
        amrex::Real speed{0};   //!< U (m/s)
        amrex::Real rate{0};    //!< f U / dx (1/s), the rate where w = 1
    };

    //! The source on the parent's layout with factor f (nest.feedback_factor;
    //! aborts unless finite and > 0) and the parent's geometry (for dx); zero
    //! and inactive until set.
    NestRelax (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom,
               amrex::Real factor);

    //! q* = qstar (u, v, w, theta, k on the valid points; q* is read only
    //! where c > 0, inside the footprints) and the rates c of the regions'
    //! footprints (DefectRegion; a direction where inner spans the footprint
    //! has no faces). Sets the source active. Aborts if U is not finite.
    //! Records the change of the target (change(): max |q*_new - q*_old|
    //! where c > 0; size the largest |u*|, |v*|, |w*| there for each velocity
    //! component, 0 for theta, the largest k* for k; first on the first set).
    void set (const State& qstar, const std::vector<DefectRegion>& regions);

    //! R += c (q - q*) for u, v, w, theta and k (nothing while inactive).
    void add_to_residual (const State& q, ResidualFields& R) const;
    //! A(coef_P) += c for eq's rows (nothing while inactive).
    void add_to_implicit (Equation eq, amrex::MultiFab& A) const;

    [[nodiscard]] bool active () const { return m_active; }
    [[nodiscard]] amrex::Real factor () const { return m_factor; }
    //! The rate c (1/s) at the points of loc (xface, yface, zface, cell).
    [[nodiscard]] const amrex::MultiFab& rate (Loc loc) const;
    //! The target q* of eq.
    [[nodiscard]] const amrex::MultiFab& target (Equation eq) const;
    [[nodiscard]] const std::vector<Nest>& nests () const { return m_nests; }
    //! The change of the target at the last set.
    [[nodiscard]] const FeedbackChange& change () const { return m_change; }

    //! max |q - q*| and its RMS over the points with c > 0 (each point once),
    //! for u, v, w, theta and k: 10 values, the maxima first. Zeros while
    //! inactive. Every rank gets the same values.
    [[nodiscard]] std::vector<amrex::Real> differences (const State& q) const;

private:
    amrex::Geometry m_geom;
    amrex::Real m_factor{1};
    std::array<amrex::MultiFab, 5> m_target;   //!< u, v, w, theta, tke
    std::array<amrex::MultiFab, 4> m_rate;     //!< at XF, YF, ZF, CC
    std::vector<Nest> m_nests;
    FeedbackChange m_change;
    bool m_active{false};
};

//! nest.feedback_factor's check: finite and > 0. Empty: no problem.
std::string check_feedback_factor (amrex::Real factor);

} // namespace gust

#endif
