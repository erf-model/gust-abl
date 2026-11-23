#ifndef GUST_NEST_FACE_FLUX_H_
#define GUST_NEST_FACE_FLUX_H_

// nest.feedback_flux = nest: the conservative variant of the defect source
// (NestDefect.H; the conservative local defect
// correction of Anthonissen, Bennink, Mattheij and ten Thije Boonkkamp,
// Numerical Algorithms 32, 2003).
//
// With nest.feedback = defect alone (feedback_flux = outer) the parent's
// equations on chi_s are R_c(q_c) = R_c(q*). The row just outside chi_s (the
// ring) keeps its own equation; its face on chi_s's edge carries the parent
// scheme's flux of q*, G_c(q*), not the nests' flux through that face. On a
// lee slope the two differ, and the parent's cells next to the source see a
// flux the nest does not.
//
// Here the ring's equation takes the nests' own flux through that face:
//
//   R_c(q_c) - S = 0   on the ring, with
//   S = s (1/J) dxi^-1 [ (G_c - G_f) + w ((G_c^up - G_c) - (G_f^up - G_f)) ],
//
// with
//   G_c     the parent's flux (convective minus diffusive, as the residual
//           takes them: after the wall model and the building closures) of
//           q* through the face;
//   G_f     the nest's flux of q_f summed over the nest faces the parent face
//           covers (the mean of the nest's flux densities; NestFaceFlux.cpp
//           gives the weights on the staggered faces);
//   G^up    the same with the upwind limiter, w the ring point's upwind blend
//           weight (convection.hybrid_sigma_max; both with w = 0 without the
//           blend). The parent's weight is used for both: the ring's own
//           blend then sees the nest's flux in place of its own;
//   s       +1 on the ring below chi_s along the face's normal (the face is
//           the ring point's high face), -1 above;
//   J       the ring point's control-volume Jacobian, dxi^-1 the parent's
//           inverse cell size along the normal.
// The ring's flux through the face is then G_f + [G_c(q_c) - G_c(q*)]: the
// nest's flux plus the outer flux of what q_c and q* still differ by on the
// ring and its neighbours (q* is the restricted nest state; the ring keeps
// its own equation). Only the convective and diffusive fluxes are corrected:
// the pressure gradient, continuity and the source terms are the outer
// domain's. Where a parent face is nodal in the normal direction (u in x, v
// in y) G_f is interpolated from the two middle nest faces, not summed, and
// the blend term matches the nest's only where the nest's own blend weight
// next to the face is the ring's.
// The correction is a constant of the parent's iteration (as S on chi_s); it
// is formed with S at every exchange from the parent's residual at q* and
// one more residual of each nest at its current state (the cost of the
// option: one nest residual per exchange).
//
// Only the lateral faces of chi_s (no z faces: chi_s spans all z) in the
// directions a nest has faces (not a periodic direction it spans). k has no
// upwind blend (bounded upwind). Nests at ratio r: r must be even (named
// nests are at 2).

#include "NestDefect.H"
#include "Residual.H"

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace gust {

struct JointMember;
class Mesh;
struct State;

//! One nest for NestFaceFlux::compute.
struct FaceFluxNest
{
    const DefectRegion* region{nullptr};   //!< its region on the parent
    const amrex::Geometry* geom{nullptr};  //!< the nest's geometry (its periodicity)
    int ratio{2};
    amrex::IntVect offset;
    //! The nest's residual at its state with this probe (ResidualInputs::flux_probe).
    std::function<void(FluxProbe&)> residual;
};

class NestFaceFlux
{
public:
    //! The correction on the parent's layout (ba, dm), zero until compute.
    NestFaceFlux (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const Mesh& mesh,
                  const amrex::Geometry& geom);

    //! The correction of every nest's ring (the S above). coarse: the parent's residual at q* with this probe.
    void compute (const std::function<void(FluxProbe&)>& coarse, const std::vector<FaceFluxNest>& nests);

    //! The joint run's exchange (JointRun.cpp): the parent's residual at q* (qr) into R, the defect source on chi_s
    //! (NestDefect::set) and this correction on the rings (NestDefect::add).
    void set_defect (const std::vector<JointMember*>& members, const State& qr, const std::vector<DefectRegion>& regions,
                     ResidualFields& R);

    [[nodiscard]] const ResidualFields& correction () const { return m_C; }
    //! The largest |correction| of u, v, w, theta and k over the parent (after compute).
    [[nodiscard]] std::vector<amrex::Real> max_norms () const;
    //! "; ring max |C| u .. v .. w .. theta .. k ..": the joint run's log line ends with it.
    [[nodiscard]] std::string log_suffix () const;

private:
    const Mesh* m_mesh;
    const amrex::Geometry* m_geom;
    ResidualFields m_C;
    // The parent's and the nests' ring sums (convective, upwind, diffusive) on the parent's layout.
    std::array<ResidualFields, 3> m_A, m_F;
    std::vector<ResidualFields> m_Rf;   //!< the nests' residuals (set_defect)
    bool m_logged{false};
};

} // namespace gust

#endif
