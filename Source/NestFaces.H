#ifndef GUST_NEST_FACES_H_
#define GUST_NEST_FACES_H_

// The lateral faces of a nested patch that take their values from the parent
// run (bc.xlo, bc.xhi, bc.ylo, bc.yhi = nest; NestParent.H reads the parent).
// One way: the parent is finished and frozen, so everything here is computed
// once, at set-up.
//
// A nest face (plan sub-phase 9.2):
//  * the normal velocity on the face is the parent's (a Dirichlet row,
//    R = q - V, identity on the implicit side);
//  * the two ghost layers outside it hold the parent's u, v, w, theta and k,
//    interpolated at the ghost points themselves (NestParent::interpolate),
//    not mirrored: a ghost value does not depend on the patch's solution, so
//    the first interior row's coupling to it is dropped from the implicit
//    side (as for the inflow profile);
//  * p' has zero gradient through the face (Neumann, as an inflow face).
// A face holds inflow and outflow parts as the parent's flow gives them;
// nothing here tells them apart on a hard face.
//
// Net flux. With Neumann pressure faces all round, the pressure correction
// can only remove a divergence whose sum over the cells is zero, that is, a
// boundary flux with no net. The parent's interpolated face velocities have a
// small net (the parent conserves mass on its own cells, not on the patch's
// faces after interpolation and with the patch's own ground). With
//   W_f = a_f / dx_d      the weight of face point f (a_f: Mesh::flux_area,
//                         dx_d the cell size along the face's normal; W_f
//                         times the cell's index volume is the face's area)
//   s_f = -1 on a low face, +1 on a high face (outward)
// the net and the correction are
//   net = sum_f s_f W_f V_f,   delta = net / sum_f W_f,   V_f <- V_f - s_f delta:
// every nest face point's outward velocity drops by the same delta, so each
// face's flux changes in proportion to its area. Applied when no face is an
// outflow (an outflow face's p' = 0 carries the net instead; delta = 0 then).
// The sums are in double and in a fixed order (Reductions.H): the same bits
// for any ranks and DistributionMapping of a fixed BoxArray. The ghost layers
// of the normal velocity outside the face drop by delta with it.
//
// Relaxation zone (bc.nest_relax_width = N > 0 cells; Davies 1976, Q. J. R.
// Meteorol. Soc. 102, 405-418). Inside each nest face, where the parent's
// flow enters,
//   R += c (q - q_parent),   c = w(s) / tau,   w = 1 - s / N,
//   1 / tau = bc.nest_relax_factor * U_in / (N dx_d),
// for u, v, w, theta and k; s is the point's distance from the face in cells
// (m + 1/2 for the cell m from the face, m for the normal velocity's point m),
// so w falls linearly from 1 at the face to 0 at the zone's inner end, and
// U_in = max(0, the parent's inward normal velocity on the face) in the
// point's own row (the mean of the face points around a point staggered
// along the face): c = 0 where the parent's flow leaves. The face point
// itself stays the Dirichlet row. c is on the implicit diagonal; it does not
// depend on the solution, so the diagonal is the source's exact derivative.
// In a corner the two faces' rates add. Width 0: no zone, no kernel runs.
//
// Ground blend (terrain.nest_blend_width = W > 0 m): blend_nest_ground.

#include "BoundaryConditions.H"
#include "ImplicitOperator.H"
#include "Mesh.H"
#include "NestParent.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace gust {

//! The checks of a run with a nested patch that need more than the bc.*
//! keys (check_boundary_inputs has the rest), one message per problem naming
//! the key: bc.nest_parent given exactly when it is read; a nest face needs
//! the parent's field as the first guess (or a restart); no grid sequencing;
//! terrain.nest_blend_width finite, >= 0, below half the patch's width across
//! its nest faces, and only with a nest face.
std::vector<std::string> check_nest_run (const BoundaryInputs& bin, const Inputs& in, bool first_guess_nest,
                                         bool sequencing);

//! terrain.nest_blend_width = W: the patch's ground h (the k = 0 slab of
//! z_phys_nd, ghost nodes included, as fill_terrain_surface leaves it)
//! becomes (1 - s) h_parent + s h, s = t^2 (3 - 2 t), t = min(d / W, 1), d the
//! distance to the nearest nest face and h_parent the parent's ground at the
//! node (NestParent::interpolate_ground): the parent's ground on the face
//! (s = 0 exactly), the patch's own from W inwards (s = 1 exactly), with a
//! continuous slope between. Ghost nodes outside a non-periodic face take the
//! nearest domain node's value. Prints one line with the largest change of
//! the ground (m). Nothing with W = 0.
void blend_nest_ground (const NestParent& parent, const BoundaryInputs& bin, amrex::Real width,
                        const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd);

//! NestFaces::summary warns when the net flux removed is above this share of the flux entering.
inline constexpr double kNestNetWarning = 1.0e-3;

class NestFaces
{
public:
    //! The net flux through the nest faces and what was done with it.
    struct Flux
    {
        double net{0};            //!< sum s W V before the correction (m/s times weight)
        double inflow{0};         //!< sum W max(0, -s V): the flux entering
        double weight{0};         //!< sum W
        amrex::Real delta{0};     //!< the velocity removed (m/s), net / weight formed in double; 0 when not applied
        bool applied{false};
    };

    //! The parent's fields on the patch (mesh's layout), the flux correction
    //! and the relaxation rates. Aborts if no lateral face of bin is nest.
    NestFaces (const NestParent& parent, const BoundaryInputs& bin, const amrex::Geometry& geom, const Mesh& mesh);

    //! Whether face (d, lo) is a nest face.
    [[nodiscard]] bool is_nest (int d, bool lo) const { return m_nest[static_cast<std::size_t>(2 * d + (lo ? 0 : 1))]; }

    //! The ghost cells of u, v, w, theta and k outside the nest faces normal
    //! to d (two layers; over the ghosts of the directions before d, as the
    //! other fills) and the normal velocity's face points: the parent's
    //! values, or 0 for a homogeneous fill. p is not touched.
    void fill (State& q, int d, bool homogeneous) const;

    //! The Dirichlet rows: R = q - V on the nest faces' normal-velocity
    //! points (V = 0 when homogeneous).
    void apply_rows (const State& q, ResidualFields& R, bool homogeneous) const;

    [[nodiscard]] bool relax_on () const { return m_relax_width > 0; }
    //! The zone's source c (q - q_parent) added to R (nothing when off).
    void add_relaxation (const State& q, ResidualFields& R) const;
    //! Its diagonal c added to A (nothing when off).
    void add_relaxation_to_implicit (Equation eq, Loc loc, amrex::MultiFab& A) const;
    //! The rate c (1/s) at the points of loc (xface, yface, zface, cell); empty when off.
    [[nodiscard]] const amrex::MultiFab& rate (Loc loc) const;

    //! The parent's fields on the patch (after the flux correction).
    [[nodiscard]] const State& parent_fields () const { return m_par; }
    [[nodiscard]] const Flux& flux () const { return m_flux; }
    //! Lines for the log: the faces, the net flux and the zone.
    [[nodiscard]] std::string summary () const;

    //! The record a restart checks: the four lateral faces' types, the
    //! parent's step, cells and content (NestParent::content_digest), the
    //! patch's place in it, the zone's settings and terrain.nest_blend_width.
    [[nodiscard]] static std::string hash_inputs (const NestParent& parent, const BoundaryInputs& bin, amrex::Real blend_width);
    static void write_checkpoint (const std::string& dir, const std::string& hash);
    //! Aborts, naming io.restart, if the checkpoint dir was written with other
    //! nest settings (hash "": this run has no nest face).
    static void check_restart (const std::string& dir, const std::string& hash);

private:
    std::array<bool, 4> m_nest{};   //!< xlo, xhi, ylo, yhi
    amrex::Geometry m_geom;
    State m_par;                    //!< u, v, w, theta, tke of the parent (p and the heights released)
    Flux m_flux;
    int m_relax_width{0};
    amrex::Real m_relax_factor{1};
    std::array<amrex::MultiFab, 4> m_rate;   //!< at XF, YF, ZF, CC
};

} // namespace gust

#endif
