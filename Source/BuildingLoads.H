#ifndef GUST_BUILDING_LOADS_H_
#define GUST_BUILDING_LOADS_H_

// The force of the air on the buildings (docs/buildings.rst), summed from the momentum
// residual on the blocked faces:
//
//   F_c = - sum over the blocked c-faces f of R_c(f) V(f),
//
// R_c the residual of momentum component c per unit volume (m/s^2, Residual.H) before the
// buildings' fixed rows replace it (ResidualInputs::building_rows = false), V(f) the face's
// control volume, J_f dx dy dzeta with J_f the mean of its two cells' detJ (CVFluxes.H). On a
// fixed face u_c = 0, so R_c is the force per unit mass that holds it at rest; its sum over a
// building is the force on the fluid from the building, and F the force on the building. The
// pressure gradients telescope through the building: the pressures inside it cancel, and the
// sum keeps p on its outer faces. The convective and viscous fluxes between two fixed faces
// cancel the same way. F is kinematic (per unit density, m^4/s^2). Not in F: the faces on a
// domain face that is not periodic, and w at the ground (the boundary conditions' rows), so a
// building on the ground has no vertical force from its base; a staircase building's side
// stresses are the viscous fluxes across its blocked faces.
//
// buildings.force_file: a CSV, one row "t,Fx,Fy,Fz" per physical step of an unsteady run
// (and its first state), one row at the end of a steady run (t = 0).

#include "Buildings.H"
#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>

namespace gust {

//! F_c = - sum over the blocked c-faces of R_c V (above), with the domain faces that are not
//! periodic and the ground's w left out. R: the residual before the buildings' fixed rows.
//! A face shared by two boxes counts once. Collective.
std::array<amrex::Real, 3> building_force (const Buildings& b, const Mesh& mesh, const amrex::Geometry& geom,
                                           const ResidualFields& R);

//! The force file of a run (buildings.force_file).
class BuildingLoads
{
public:
    explicit BuildingLoads (std::string file) : m_file(std::move(file)) {}
    [[nodiscard]] const std::string& file () const { return m_file; }

    //! Start the file (the I/O rank writes; collective). keep_until >= 0 (a restart): keep the
    //! rows with t <= keep_until of an existing file; < 0: replace any file. Returns true when it
    //! kept an existing file, false when it started a new one (the caller then writes the first
    //! row). Aborts (FileOpenFailed) when the file cannot be written.
    bool open (double keep_until) const;
    //! The force at state q, whose ghost cells the caller has filled: the residual with the
    //! buildings' rows left out (in.building_rows = false, in.time = null: the time term is 0 on
    //! a fixed face). Collective.
    std::array<amrex::Real, 3> evaluate (const State& q, ResidualInputs in, const Mesh& mesh, const amrex::Geometry& geom,
                                         const Buildings& b);
    //! Append "t,Fx,Fy,Fz" (I/O rank).
    void write_row (double t, const std::array<amrex::Real, 3>& F) const;

private:
    bool open_io (double keep_until) const;   //!< open() on the I/O rank
    std::string m_file;
    ResidualFields m_R;
    ResidualWorkspace m_ws;   //!< the residual's scratch, kept between the steps of an unsteady run
};

} // namespace gust

#endif
