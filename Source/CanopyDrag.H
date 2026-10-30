#ifndef GUST_CANOPY_DRAG_H_
#define GUST_CANOPY_DRAG_H_

// The canopy's drag on the momentum equations (docs/canopy.rst).
//
// Each velocity face c (c = x, y, z) takes the sink
//
//   S_c = -f_c |U| u_c,     |U| = sqrt(u^2 + v^2 + w^2) at the face,
//
// as ERF's ERF_MakeMomSources.cpp: f_c is the mean of the drag field f =
// Cd a(z) (Canopy.H) of the two cells that share the face, and the other
// two components are the four-point means of their faces around it
// (Stencil_K.H). A z-face on the domain's ground or top takes the one cell
// inside (those rows are the boundary rows of w, which the boundary
// conditions set). The residual R_c gains +f_c |U| u_c (per unit volume, as
// every residual: Residual.H).
//
// Implicit side (canopy.drag_diagonal): the u_c rows' a_P gains
//   derivative  f_c (|U| + u_c^2 / |U|) = d(f_c |U| u_c)/d u_c   (default)
//   picard      f_c |U|                    (|U| frozen)
//   none        nothing                    (ERF's explicit treatment)
// Each is >= 0, so 1/dtau + A1 keeps a positive diagonal; it reaches the
// line solves, SIMPLEC's a_P and its d = a_P - sum a_nb (implicit_rows,
// Residual.H). The fixed point does not depend on the choice. The drag does
// not enter the local pseudo-time step (PseudoTime.H).

#include "Canopy.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>

namespace gust {

//! The name of a diagonal option ("derivative", "picard", "none").
std::string canopy_diagonal_name (CanopyDiagonal d);

//! The drag field on the velocity faces: f[c] is f_c on the c-faces of the
//! level (no ghost faces), from the level's Canopy, and the diagonal option.
class CanopyFaceDrag
{
public:
    //! From canopy.drag() (one ghost column in x and y; boxes span the
    //! column) on that field's layout, with canopy.drag_diagonal(). Aborts
    //! unless canopy.active().
    void build (const Canopy& canopy, const amrex::Geometry& geom);
    [[nodiscard]] bool empty () const { return m_f[0].empty(); }
    [[nodiscard]] const amrex::MultiFab& f (int c) const { return m_f[static_cast<std::size_t>(c)]; }
    [[nodiscard]] CanopyDiagonal diagonal () const { return m_diagonal; }

private:
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> m_f;
    CanopyDiagonal m_diagonal{CanopyDiagonal::derivative};
};

//! R_c += f_c |U| u_c on every valid face of u, v, w (R on the faces, no
//! ghosts, on the drag's layout). u, v and w need one ghost cell in x, y and
//! z (filled by the caller, as for the residual; the ground and top faces of
//! w read u and v below and above the domain, rows the boundary conditions
//! then replace). Aborts on another layout or fewer ghost cells.
void add_canopy_drag (const CanopyFaceDrag& drag, const amrex::MultiFab& u, const amrex::MultiFab& v,
                      const amrex::MultiFab& w, amrex::MultiFab& Ru, amrex::MultiFab& Rv, amrex::MultiFab& Rw);

//! Component coef_P of A (component c's rows, ImplicitOperator.H) gains the
//! diagonal of drag.diagonal() at the state (u, v, w), on A's valid points.
void add_canopy_drag_to_implicit (const CanopyFaceDrag& drag, int c, const amrex::MultiFab& u,
                                  const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& A);

//! The diagonal of each component's rows at one state, the linearisation
//! point of a solve: ap[c] on the c-faces (one component, no ghosts). A
//! solver computes it once, before it updates any component (the SIMPLEC
//! predictor at q, the Newton preconditioner at its q0, not at the increment
//! it solves for), and hands it to implicit_rows (Residual.H).
struct CanopyRows
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> ap;
    //! The theta rows' diagonal of the leaf heat (canopy_heat_rows, CanopyHeat.H) on the cells; empty without it.
    //! Like ap, it is the caller's to refresh at each linearisation state (implicit_rows checks only its layout).
    amrex::MultiFab theta;
};

//! rows.ap[c] = the diagonal of drag.diagonal() at (u, v, w) (defined here if
//! empty); ghost cells of u, v, w as for add_canopy_drag.
void canopy_rows (const CanopyFaceDrag& drag, const amrex::MultiFab& u, const amrex::MultiFab& v,
                  const amrex::MultiFab& w, CanopyRows& rows);

} // namespace gust

#endif
