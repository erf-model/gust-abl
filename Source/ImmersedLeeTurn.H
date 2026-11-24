#ifndef GUST_IMMERSED_LEE_TURN_H_
#define GUST_IMMERSED_LEE_TURN_H_

// The lee turn of an immersed ground (terrain.lee_turn; docs/immersed_terrain.rst, "Lee turn").
//
// In the cut cells of an immersed ground (terrain.representation = immersed_forcing or hybrid) the ground
// is a staircase of partly solid cells. Where the ground falls away along the wind (the lee side), the
// cut cells pass the flow on through their open faces, but nothing turns it down with the ground: w in a
// cut cell sees the closed part of its faces and the wall law, not the slope. Over a ridge, the wind at
// 10 m above the ground in the lee then reads too slow on coarse meshes, and the deficit grows with the
// number of layers the ground crosses over one cell.
//
// The lee turn relaxes w in a lee cut cell (0 < beta < 1) towards the vertical velocity of a flow that
// follows the ground,
//   w_t = s_x u_r + s_y v_r,
// s = grad h the ground's slope (Buildings::surface(), kSurfaceSlopeX/Y) and (u_r, v_r) the horizontal
// wind at the first fully fluid cell above in the column (beta = 0, the reference cell; each component the
// open_face_mean of its two faces), as a sink on the w rows,
//   R_w += lambda (w - w_t),   a_P += lambda,
// with
//   lambda = f C_0 n^2 |U_r| / dz,   n = max(1, sqrt((s_x dx)^2 + (s_y dy)^2) / dz),
//   f = min(1, -w_t / (epsilon |U_r|)),
// dz = J dzeta the cut cell's height, |U_r| = max(sqrt(u_r^2 + v_r^2), kBuildingMinSpeed), C_0 =
// terrain.lee_turn and epsilon = terrain.lee_turn_ramp. n is the number of layers the ground crosses over one
// cell. f brings the term in from 0 where w_t changes sign (a windward cut cell, w_t >= 0, has none), so the
// rate does not switch on and off between iterations. The form of the rate, C |U| / dz with C of order one
// over a cell, is that of Kynema-SGF's drag forcing (DragForcing.cpp, Cd = 10 / dz before its 1 / dt cap,
// which a steady run does not need).
//
// A w row on a z face takes the mean of its two cells' lambda and lambda w_t, as BuildingForcing's rates do.
// The ground's own row (w at the mesh's lower face, the boundary conditions' row) and blocked faces (fixed
// rows) take nothing. The term is implicit in a_P; lambda and w_t are lagged (recomputed at each residual
// evaluation from the current wind), and the converged state is the same as with them updated in the solve.

#include "Buildings.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <cmath>

namespace gust {

//! The components of ImmersedLeeTurn::rates(): lambda (1/s) and lambda w_t (m/s^2).
inline constexpr int kLeeTurnRate = 0, kLeeTurnTarget = 1, kNumLeeTurnRates = 2;

//! The lee turn of one cut cell: lambda and lambda w_t from the ground's slope (sx, sy), the reference wind
//! (ur, vr, m/s), the cut cell's height dz and the mesh's dx, dy (m), C_0 = c0 and epsilon = eps. Both 0 on
//! the windward side (w_t >= 0). min_speed: kBuildingMinSpeed.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void lee_turn_rate (amrex::Real sx, amrex::Real sy, amrex::Real ur, amrex::Real vr, amrex::Real dz, amrex::Real dx,
                    amrex::Real dy, amrex::Real c0, amrex::Real eps, amrex::Real min_speed, amrex::Real& lam,
                    amrex::Real& lam_wt) noexcept
{
    lam = 0;
    lam_wt = 0;
    const amrex::Real wt = sx * ur + sy * vr;
    if (!(wt < amrex::Real(0))) { return; }
    const amrex::Real U = amrex::max(std::sqrt(ur * ur + vr * vr), min_speed);
    const amrex::Real f = amrex::min(amrex::Real(1), -wt / (eps * U));
    const amrex::Real step = std::sqrt((sx * dx) * (sx * dx) + (sy * dy) * (sy * dy));
    const amrex::Real n = amrex::max(amrex::Real(1), step / dz);
    lam = f * c0 * n * n * U / dz;
    lam_wt = lam * wt;
}

//! The lee turn's rates on one level, and their terms in the w residual and rows.
class ImmersedLeeTurn
{
public:
    //! Off unless c0 > 0 and the ground was built from a surface (Buildings::from_surface(): the slopes).
    //! A flat ground has zero slopes and a zero rate (the inputs refuse terrain.lee_turn > 0 there).
    //! b, detJ must outlive this object.
    void build (const Buildings& b, const amrex::Geometry& geom, const amrex::MultiFab& detJ, amrex::Real c0, amrex::Real eps);
    [[nodiscard]] bool active () const { return m_active; }
    //! The rates at the cut cells from u and v (valid faces of the reference cells) and the faces' weights
    //! in a cell's speed (BuildingForcing's, one ghost layer in x and y). Every box spans the column
    //! (Buildings::build_from_surface aborts otherwise), so the reference cell above is in the same box. 0 elsewhere, ghost cells included.
    void update (const amrex::MultiFab& u, const amrex::MultiFab& v, const std::array<amrex::MultiFab, 3>& weight);
    //! Rw += lambda w - lambda w_t on the w rows (z faces), each the mean of its two cells; not on the
    //! ground's row or on blocked faces.
    void add_sink (const amrex::MultiFab& w, amrex::MultiFab& Rw) const;
    //! a_P of the same w rows gains the same mean lambda.
    void add_to_implicit (amrex::MultiFab& A) const;
    //! kLeeTurnRate and kLeeTurnTarget at the cells, one ghost layer (empty when not active).
    [[nodiscard]] const amrex::MultiFab& rates () const { return m_rates; }

private:
    const Buildings* m_b{nullptr};
    const amrex::MultiFab* m_detJ{nullptr};
    amrex::Geometry m_geom;
    amrex::MultiFab m_rates;
    amrex::Real m_c0{0};
    amrex::Real m_eps{amrex::Real(0.05)};
    bool m_active{false};
};

}   // namespace gust

#endif
