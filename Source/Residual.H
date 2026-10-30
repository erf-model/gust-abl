#ifndef GUST_RESIDUAL_H_
#define GUST_RESIDUAL_H_

// The steady residual and its linearisations:
//
//   R_c     = C_c(u) - D_c(u) + (G p)_c - S_c      momentum, c = x, y, z
//   R_p     = D(U)                                 continuity (the divergence)
//   R_theta = C(theta) - D_gamma(theta)
//   R_k     = C_bounded,upwind(k) - D_gamma_k(k)
//
// with C, D the convection and diffusion operators (Convection.H,
// Diffusion.H), G and D(U) the pressure gradient and divergence
// (MassFlux.H), and S the sources (Forcing.H): buoyancy on w, Coriolis and
// geostrophic forcing on u and v. The transport fluxes of C come from the
// state's velocity (transport_fluxes). p is the kinematic perturbation
// pressure. Every residual is per unit volume (m/s^2 for momentum, 1/s for
// continuity, K/s for theta, m^2/s^3 for k). Diffusion is assembled from its
// face fluxes, so that with a wall model (ResidualInputs::surface,
// SurfaceLayer.H) the ground faces carry the wall stress, the surface heat
// flux and no k flux, and the first cells' k rows are the wall-k rows; the
// surface layer is computed from q, or taken from the caller's ghost fill of
// q (the same state, unless it read the top w face: Residual below). The closure sources of k
// (production, buoyancy, dissipation) and the boundary conditions come from
// Closure.H and BoundaryConditions.H; the residual reads the state's ghost
// cells (two layers), which the caller fills, and leaves the boundary rows
// (the ground and top faces of w, faces on non-periodic sides) to the
// boundary-condition code.
//
// Contract for the later terms: everything that belongs to R enters inside
// residual() (through ResidualInputs), so that jvp, which differentiates
// residual(), sees it. With fixed diffusivities (inputs) jvp is the
// Jacobian with nu_t frozen (a Picard linearisation of diffusion); with a
// closure (ResidualInputs::closure) residual() computes nu_t from q, so jvp
// differentiates nu_t too.
//
// line_jacobian gives, for one field and one direction, the tridiagonal
// coefficients of the first-order implicit-side operator A1 along lines in
// that direction (ImplicitOperator.H), built from the state's transport.
// A1 and R are per unit volume, so the pseudo-time term the solver adds to
// the diagonal is 1/dtau (V/dtau is the same term integrated over a cell of
// volume V).
// jvp gives J v by a forward difference of the full residual.

#include "Closure.H"
#include "Convection.H"
#include "Forcing.H"
#include "Mesh.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>

namespace gust {

class CanopyFaceDrag;
class CanopyHeat;
struct CanopyRows;

class BoundaryConditions;

//! The residual of each equation, on the state's locations, without ghosts.
struct ResidualFields
{
    amrex::MultiFab u, v, w, p, theta, tke;
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
    struct Field { const char* name; amrex::MultiFab* mf; };
    amrex::Vector<Field> fields ();
    [[nodiscard]] amrex::Vector<const amrex::MultiFab*> fields () const;
};

//! Typical magnitudes of the fields, for jvp's step (> 0, finite), and the
//! step's relative size.
struct JvpScales
{
    amrex::Real velocity{10};      //!< m/s
    amrex::Real theta{10};         //!< K, for theta - theta_ref
    amrex::Real tke{amrex::Real(0.1)};  //!< m^2/s^2
    amrex::Real pressure{100};     //!< m^2/s^2 (velocity squared)
    //! The relative step h of jvp's eps (> 0, finite): sqrt(machine epsilon)
    //! of amrex::Real by default, the step that balances the truncation and
    //! rounding errors of a forward difference. The Newton core sets it
    //! from steady.newton_fd_eps_scale (Newton.H).
    amrex::Real step{std::sqrt(std::numeric_limits<amrex::Real>::epsilon())};
};

//! What the residual needs besides the state. Diffusivities are cell-centred
//! with one ghost cell filled. Either the three diffusivity fields are given
//! (held fixed), or a closure computes them from the state
//! (Closure.H): then nu, gamma_theta and gamma_tke must be null, wall_distance
//! must be set, and the residual of k gains the closure's sources.
struct ResidualInputs
{
    const amrex::MultiFab* nu{nullptr};           //!< momentum viscosity nu + nu_t (m^2/s)
    const amrex::MultiFab* gamma_theta{nullptr};  //!< theta diffusivity, kappa + nu_t / Pr_t (m^2/s)
    const amrex::MultiFab* gamma_tke{nullptr};    //!< k diffusivity, nu + nu_t / sigma_k (m^2/s)
    //! The face average of the diffusivities (Diffusion_K.H); with a closure,
    //! the closure's option is used instead.
    FaceMean face_mean{FaceMean::arithmetic};
    //! The k–L closure (null: the diffusivities above). Its roughness is the
    //! surface layer's z0 (0 without one, as ERF).
    const Closure* closure{nullptr};
    //! The wall distance the closure reads (cell-centred, one ghost cell):
    //! height_above_ground() or the Tucker distance.
    const amrex::MultiFab* wall_distance{nullptr};
    ReferenceProfile theta_ref;                   //!< hydrostatic reference theta_ref(z), a 1-D profile
    ConvectionOptions momentum;                   //!< momentum convection (k is always bounded upwind)
    ConvectionOptions theta;                      //!< theta convection
    //! The upwind blend of the momentum and theta convection (HybridConvection.H; null: none): a cell field
    //! with one ghost cell, owned by the SIMPLEC stage.
    const amrex::MultiFab* hybrid_sigma{nullptr};
    HybridOptions hybrid;                         //!< convection.hybrid_* (the stage reads them)
    BuoyancyOptions buoyancy;
    //! false with simplec.theta = frozen: theta's buoyancy acts in the closure
    //! (N^2) only, not as a force in the w equation (the prescribed
    //! stratification of van der Laan et al. 2024, whose momentum equations
    //! carry no buoyancy).
    bool buoyancy_in_momentum{true};
    CoriolisOptions coriolis;
    //! The constant kinematic pressure-gradient force (m/s^2) on u, v, w (Forcing.H).
    amrex::Array<amrex::Real, 3> pressure_gradient{0, 0, 0};
    JvpScales jvp_scales;
    //! The wall model (SurfaceLayer.H) on the ground; null: none, and the
    //! ground fluxes come from the ghost cells like any interior face.
    const SurfaceLayerOptions* surface{nullptr};
    //! The domain boundary conditions (BoundaryConditions.H);
    //! null: none. With them the residual sets the Dirichlet rows of the
    //! normal velocities and adds the damping layers, and line_jacobian
    //! gives those rows the identity and folds the ghost couplings. The
    //! caller still fills the state's ghosts (BoundaryConditions::fill_ghosts)
    //! before calling: the ghost fill is a function of q, so a perturbation
    //! vector gets the homogeneous fill. With bc set, surface must match
    //! bc.zlo = surface_layer (set) or not (null).
    const BoundaryConditions* bc{nullptr};
    //! The canopy drag on the velocity faces (CanopyDrag.H); null: none. The
    //! residual adds f_c |U| u_c to the momentum rows, and implicit_rows the
    //! diagonal of its option at the solver's linearisation state (CanopyRows).
    const CanopyFaceDrag* canopy{nullptr};
    //! The leaves' heat exchange (CanopyHeat.H); null: none. The residual adds
    //! lambda (theta - theta_l) to R_theta, and implicit_rows the theta rows'
    //! lambda at the solver's linearisation state (CanopyRows::theta). Needs
    //! canopy as well.
    const CanopyHeat* canopy_heat{nullptr};
};

//! Scratch fields of residual() and jvp, kept between calls so that a
//! caller that evaluates often (the stages, the Newton core and its
//! preconditioner) allocates them once per layout. Defined on first use; a
//! call on another BoxArray or DistributionMapping redefines them. Every
//! field is written before it is read in each call, so a workspace gives the
//! same bits as fresh fields. One workspace serves one evaluation at a time.
struct ResidualWorkspace
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;            //!< transport fluxes
    amrex::Array<amrex::MultiFab, 3> tmp;                       //!< momentum terms per component
    amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, 3> Gm;   //!< momentum diffusive fluxes
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> Gt, Gk;       //!< theta and k diffusive fluxes
    //! Convective fluxes, one set per location, indexed by the Loc value (cell 0, x 1, y 2, z 4; the rest stay empty).
    amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, 8> Gc;
    amrex::MultiFab theta_ref, d, omega;
    SurfaceState ss;              //!< the wall state, when the caller gives none
    ClosureFields cf;             //!< the closure's fields, when the caller gives none
    ClosureScratch closure;
    State qp;                     //!< jvp: q + eps v
    ResidualFields base;          //!< jvp: R(q), when the caller gives none
    amrex::MultiFab jvp_theta_ref;
    //! Clear every field when (ba, dm) differ from the layout of the last call.
    void use_layout (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
private:
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
};

//! R(q). q's ghost cells must be filled (two layers: MUSCL). Aborts, naming
//! the input, on missing fields or invalid options. With a closure, fields
//! (when given) receives what it evaluated at q.
//! surface (with a wall model; null: computed here): the wall state of q,
//! as BoundaryConditions::fill_ghosts returns it when it fills q's ghosts.
//! The wall model reads no value the fill writes, so that state is the one
//! residual() would compute; debug builds (GUST_DEBUG_CHECKS) compute it
//! again and abort on any difference. The exception, a state that read the
//! top w face (SurfaceState::reads_top_face, a z_ref in the top interval),
//! is not used: residual() computes its own.
//! ws (null: fields local to the call): the scratch fields.
void residual (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
               ResidualFields& R, ClosureFields* fields = nullptr, const SurfaceState* surface = nullptr,
               ResidualWorkspace* ws = nullptr);

//! The fields with an implicit-side operator.
enum class Equation { u, v, w, theta, tke };

//! Components of line_jacobian's output.
enum LineCoef : int { line_lower = 0, line_diag, line_upper, n_line_coefs };

//! Tridiagonal coefficients of A1 for one equation along lines in direction
//! dir (0, 1, 2): lower = -a_{dir-}, diag = a_P (every neighbour's weight),
//! upper = -a_{dir+}; the neighbours across the other directions are lagged.
//! tri: n_line_coefs components at the field's location (defined here if
//! empty). The caller adds 1/dtau and source Jacobians to the diagonal.
//! Uses the diffusivities of in, not its reference profile or options.
void line_jacobian (const State& q, const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh,
                    const amrex::Geometry& geom, amrex::MultiFab& tri);

//! The folded rows of A1 for one equation: all n_implicit_coefs coefficients
//! (ImplicitOperator.H) with the wall model's rows (k's wall-k rows included), the closure's
//! point-implicit k diagonal, the damping diagonal, the Dirichlet identity
//! rows and the ghost couplings of the boundary rows applied (the same rows
//! line_jacobian splits by direction). A: defined here if empty. surface,
//! fields and canopy as for line_jacobian below.
void implicit_rows (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& A, const SurfaceState* surface = nullptr, const ClosureFields* fields = nullptr,
                    const CanopyRows* canopy = nullptr);

//! The same from transport fluxes F already computed (transport_fluxes), for
//! line solvers that sweep every direction with one state.
//! With a wall model (and eq not w), surface must be the state's surface
//! layer (compute_surface_state): the first cells' rows lose the ground
//! neighbour, u and v gain the wall stress's Picard coefficient, and k's
//! rows are the wall-k rows. With a closure, fields must be its evaluation
//! at the state (Closure::evaluate): the diffusivities come from it and k's
//! diagonal gains the point-implicit source term dSdk. With a canopy
//! (in.canopy) the momentum rows need canopy, the drag's diagonal at the
//! solver's linearisation state (canopy_rows, CanopyDrag.H): their a_P
//! gains it. With the leaf heat (in.canopy_heat) the theta rows need canopy
//! too, with its theta diagonal (canopy_heat_rows, CanopyHeat.H).
void line_jacobian (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& tri, const SurfaceState* surface = nullptr, const ClosureFields* fields = nullptr,
                    const CanopyRows* canopy = nullptr);

//! Jv = (R(q + eps v) - R(q)) / eps, with
//!   eps = h (1 + mean |q~|) / rms(v~),   h = in.jvp_scales.step (sqrt(machine epsilon) by default)
//! over every owned point of every field (shared faces counted once), on
//! scaled variables q~ = q / s, v~ = v / s with s from in.jvp_scales and theta
//! taken as theta - theta_ref. A perturbation along a typical v then moves
//! each value by about h times its own scale, whatever
//! the field and the mesh size (a plain |q| would be dominated by theta's
//! 300 K and grow with the number of cells). Jv = 0 for v = 0. v is a state on
//! q's layout with its ghost cells filled; its z_phys_nd is not used. R0, when
//! given, must hold R(q) and saves one evaluation. Returns eps.
//! refill (needs in.bc): the ghost cells of q + eps v come from the full
//! fill of that state (BoundaryConditions::fill_ghosts) instead of q's
//! ghosts plus eps times v's. The fill is nonlinear in q where the wall law
//! sets the ground gradient, so only the refill differentiates it; the
//! Newton core (Newton.H) uses it. Without refill, v's ghosts (the
//! homogeneous fill) give the Jacobian with that gradient frozen. With
//! refill, the wall state of that fill goes to the residual.
//! ws (null: fields local to the call): the scratch fields, q + eps v
//! included.
amrex::Real jvp (const State& q, const State& v, const ResidualInputs& in, const Mesh& mesh,
                 const amrex::Geometry& geom, ResidualFields& Jv, const ResidualFields* R0 = nullptr,
                 bool refill = false, ResidualWorkspace* ws = nullptr);

} // namespace gust

#endif
