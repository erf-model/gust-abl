#ifndef GUST_COLUMNGUESS_H_
#define GUST_COLUMNGUESS_H_

// The converged-column first guess (init.first_guess =
// column): the deck's horizontally uniform column is solved as a 1-D
// problem on the deck's vertical grid and its converged u, v, theta and k
// profiles are returned as a table by height above the ground. The
// application turns the table into an InflowProfile and hands it to the
// initializer in place of the inflow profile, so the terrain case starts
// from a field that is in Ekman balance with the closure (the inflow
// profile is not: u along the geostrophic wind, v = 0, and the stage then
// spends the inertial adjustment in pseudo-time). The mass-consistent
// projection (init.projection) is a separate switch and applies to the
// broadcast field as to the profile.
//
// The 1-D column (run on every rank, host code in double precision, the
// same result on every rank): the cells of the deck's vertical grid (its
// nodal heights from Mesh::build on a flat column); u + i v, theta and k at
// the cell centres; backward Euler steps of growing size with
//   d(u + iv)/dt = d/dz(nu d(u + iv)/dz) - i f ((u + iv) - (u_g + i v_g)) + F_x + i F_y
//                  - c_r(z) ((u + iv) - (u_p + i v_p)),
//   dk/dt        = d/dz(gamma_k dk/dz) + P_s + P_b - eps - c_r(z) (k - k_p),
//   dtheta/dt    = d/dz(gamma_theta dtheta/dz) - c_r(z) (theta - theta_p)   (with a non-neutral surface),
// theta held at the sounding with the neutral surface (init.column_stability = mol, init.column_mol = 0, the
// default; a column with an inversion warms without end, docs/results.rst), so N^2 is then the sounding's;
// the diffusion, the Coriolis term (implicit through the complex unknown),
// the wall drag, the surface heat flux and the Rayleigh layer implicit, the k sinks point-implicit
// (Patankar), the closure's point kernels (Closure_K.H: the same algebra as
// the 3-D closure, with N^2 and S^2 from the cell's two faces) lagged one
// step. Ground: the neutral wall law (u* = kappa |U_1| / ln((d_1 + z0)/z0),
// the stress u*^2 along U_1 as an implicit drag, k_1 = u*^2/Cmu0^2 as the
// wall-k Dirichlet row, no heat flux), or with a non-neutral surface MOST for a surface temperature theta_s
// (ERF's Businger-Dyer psi_m, psi_h, or QNSE's stable side, bc.surface_stable_functions; zeta = d_1 / L as the
// 3-D wall model's; u* = kappa U_1 / (ln - psi_m), q = -kappa u* (theta_1 - theta_s) / (ln - psi_h) implicit in
// theta_1, k_1 = (u*^3 + kappa B d_1)^(2/3)/Cmu0^2 with B = g q / theta0 for an upward flux only, as the surface
// layer's). theta_s is init.column_surface_theta (surface_temperature), or found by the secant method so that
// 1/L is init.column_mol's (mol) or the flux is bc.surface_heat_flux (heat_flux). The heat budget: at steady state
// the surface flux is returned by the Rayleigh layer on theta, the column's only heat source or sink (required);
// the most stable column it allows (the smallest L, the largest downward flux) depends on that layer and the wind,
// and a request beyond it stops the run, naming it.
// Target wind (init.column_target_u, _v at init.column_target_height; amrTerrain's met-mast forcing): the
// geostrophic wind and the Rayleigh layer's reference are multiplied by R = U_target / U(z_t) (U(z_t) from MOST's
// profile below the first centre, linear in ln(z + z0) above), the pressure gradient by R |R|, and the column
// solved again, until |R - 1| < 1e-4; the application then runs with that
// forcing. Top: slip (zero flux; theta's flux
// gamma_theta bc.zhi_theta_grad), geostrophic (u, v held at u_g, v_g) or no
// slip (u = v = 0). The Rayleigh layer of the deck (bc.rayleigh_*) relaxes
// the fields bc.rayleigh_fields selects: u and v (together) towards the
// inflow profile, k towards its k, theta towards the sounding (with a non-neutral surface).
// The column takes nu_t relaxed by a half between steps (turb.nut_relax is
// the stage's) and the arithmetic face mean of the diffusivities
// (turb.face_mean is the stage's). Converged when the largest change of u, v
// (scaled by the geostrophic or first-cell speed) and k (its maximum) over
// one step is below
// init.column_tol (lifted to 100 epsilon of Real, since the closure's point
// kernels run in Real), or aborts after init.column_max_iterations. The column
// needs the surface layer with the Dirichlet wall k (neutral, or the column's stability), a
// horizontal forcing (geostrophic, or a pressure gradient) and the stage;
// check_column_guess names the key that is missing. Seconds for 64 levels.
//
// The time-marching column (init.column_time > 0, the form of ERF's near-1-D column): the
// same equations and closure in physical time, init.column_dt per backward Euler step (the last step ends at
// init.column_time), nu_t not relaxed between steps, no convergence test and no Rayleigh layer on theta needed (a
// stratified column need not be steady: the table is its state at init.column_time). Its surface is fixed: the neutral
// wall law (init.column_mol = 0), init.column_surface_theta (surface_temperature) or bc.surface_heat_flux itself
// (heat_flux: an explicit flux q into the first cell, MOST from column_surface_most_flux, no search for theta_s), so
// with a slip top, bc.zhi_theta_grad = 0 and no Rayleigh layer on theta the column's heat content
// sum(theta dz) grows by q t to round-off. No target wind (init.column_target_*) and no init.column_mol.
// init.column_wall_k = solved (time-marching only) solves the first cell's k with the others' equation and no k flux
// through the ground, in place of the wall-k Dirichlet row (fixed, the default).
//
// The forest canopy (Canopy.H), when the run has one, is in the column as the layer profile f_k of the level's drag
// field (ColumnCanopy, from Canopy::column_profile): the momentum equation gains the drag -f_k |W| W, linearised
// along the wind about the step's start, -f |W_n| W_n - 2 f |W_n| (W - W_n) (the derivative diagonal of the 3-D
// rows, CanopyDrag.H: at the column's large steps a lagged |W| would hand the balance speed back and forth, the
// map W -> G / (f |W|), while the derivative form is Newton's method on it), and with canopy.tke the k
// equation gains ERF's canopy pair (CanopyTke.H), P_c = beta_p f |W|^3 as a source and D_c = beta_d f |W| k on the
// Patankar diagonal, and the closure's length cap L <= alpha / f (Closure_K.H), the same terms as the 3-D stage's
// on a flat column, so the first guess and the inflow column carry the forest. The column is horizontally uniform,
// so the canopy must be: column_canopy stops the run, naming init.first_guess, when f differs between columns.

#include "BoundaryConditions.H"
#include "Canopy.H"
#include "Closure.H"
#include "Driver.H"
#include "Forcing.H"
#include "Initializer.H"
#include "Inputs.H"
#include "SurfaceLayer.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

namespace gust {

struct ColumnGuessResult
{
    bool steady{false};
    int iterations{0};   //!< the backward Euler steps of all the columns solved
    int columns{0};      //!< the columns solved (the surface search and the target wind solve several)
    double wall_seconds{0};
    double tol{0};   //!< the tolerance used: init.column_tol lifted to 100 epsilon of Real
    amrex::Real ustar{0};   //!< the wall law's friction velocity at the end (m/s)
    amrex::Real heat_flux{0};   //!< the surface heat flux at the end (K m/s, upward > 0; 0 when theta is held)
    amrex::Real inv_L{0};       //!< 1/L at the end (1/m; 0 neutral)
    amrex::Real theta_surface{0};   //!< the surface temperature at the end (K; theta at the first cell when held)
    //! The forcing the column ran with (adjusted by the target wind when init.column_target_* are given).
    amrex::Real ug{0}, vg{0}, pressure_gradient_x{0}, pressure_gradient_y{0};
    int target_iterations{0};   //!< the target-wind checks (the forcing adjustments plus one; 0 without a target)
    int flux_iterations{0};     //!< the surface-temperature adjustments of the flux modes (0 otherwise)
    double time{0};             //!< the time-marching column's end time (s; init.column_time), 0 for the steady column
    //! The table: the ground row (z = 0: u = v = 0, the first cell's theta and
    //! k) then the cell-centre heights above the ground.
    amrex::Vector<amrex::Real> z, u, v, theta, tke;
};

//! Whether the column evolves theta (any surface but the neutral one: init.column_mol = 0 holds theta at the
//! sounding, the form before init.column_stability).
bool column_theta_evolves (const InitOptions& iopt, const SurfaceLayerOptions& surf);

//! Businger-Dyer's stable side has a MOST root only for a bulk Richardson number below 1/beta = 0.2
//! (zeta = Rib ln / (1 - 5 Rib)); the column's surface holds it at this value, zeta = 3.8 ln((z1 + z0)/z0).
constexpr double kRibMaxBusingerDyer = 0.19;

//! MOST at the column's surface for a surface temperature theta_s, from the first cell's speed U1 and theta th1 at
//! height z1, as ERF's surface_temp (ERF_MOSTStress.H): the bulk Richardson number
//! Rib = z1 g (theta_1 - theta_s) / (U1^2 theta0), held in [rib_min, rib_max] (init.column_rib_min, _max; and at
//! most kRibMaxBusingerDyer with Businger-Dyer), gives zeta = z1/L from zeta = Rib (ln - psi_m)^2 / (ln - psi_h)
//! (ln = ln((z1 + z0)/z0), zeta as the 3-D wall model's d/L, clamped to [-100, 100]), solved by ERF's fixed point
//! under-relaxed by a half from the zeta of invL on entry (the previous step's) until the step is within
//! 1e-10 max(1, |zeta|), or (single precision only) until it cycles: within 64 amrex::Real roundings of |zeta|,
//! not shrinking for 50 passes and changing sign in at least 25 of them (ZetaStop, ColumnGuess.cpp). Then
//! u* = kappa U1 / (ln - psi_m) (above 1e-3 m/s), theta* = kappa (theta_1 - theta_s) / (ln - psi_h), q = -u* theta*,
//! 1/L = zeta/z1, and the transfer coefficient ch = kappa u* / (ln - psi_h) (q = -ch (theta_1 - theta_s)). The psi
//! functions are ERF's Businger-Dyer forms, or QNSE's on the stable side (qnse, StabilityFunctions.H); the denominators
//! are kept above the wall law's floor m (min(ln, 1), or 1 with bc.surface_denominator_floor = one). False when the
//! fixed point does not converge in 20000 passes or zeta reaches the clamp.
bool column_surface_most (double kappa, double z0, double z1, double U1, double th1, double theta_s, double grav,
                          double theta0, double den_floor, bool qnse, const StableCoefficients& sc, double rib_min,
                          double rib_max, double& ustar, double& thstar, double& qs, double& invL, double& ch);

//! MOST at the column's surface for a given kinematic heat flux qs (K m/s, upward > 0), the time-marching column's
//! heat_flux surface (ERF's surface_flux, ERF_MOSTStress.H): zeta = z1/L = -kappa g qs z1 / (theta0 u*^3) with
//! u* = kappa U1 / (ln - psi_m(zeta)) (above 1e-3 m/s; the denominator above the wall law's floor, as
//! column_surface_most), solved by the fixed point under-relaxed by a half from the zeta of invL on entry, to the
//! stop of column_surface_most (zeta clamped to [-100, 100]), and once more from zeta = 0 when that fails (a seed
//! above the upper stable root). Then theta* = -qs/u* and 1/L = zeta/z1. False when the fixed point does not
//! converge in 20000 passes or zeta reaches the clamp (a downward flux too large for the wind has no root).
bool column_surface_most_flux (double kappa, double z0, double z1, double U1, double qs, double grav, double theta0,
                               double den_floor, bool qnse, const StableCoefficients& sc, double& ustar, double& thstar,
                               double& invL);

//! Whether the 3-D first guess takes theta from the column: a stratified column (column_theta_evolves) whose surface
//! the 3-D run shares (bc.surface_stability = heat_flux). Otherwise theta is the sounding.
bool column_sets_theta (const InitOptions& iopt, const SurfaceLayerOptions& surf);

//! MOST's gradients at the first centre, a distance d above the ground, as the 3-D ground fill's (BoundaryConditions.cpp):
//! dU/dz = u* phi_m(zeta) / (kappa (d + z0)) and dtheta/dz = theta* phi_h(zeta) / (kappa (d + z0)), zeta = d / L clamped
//! to [-100, 100], with the stable side sc (StabilityFunctions.H most_phi; the default Businger-Dyer). The column's
//! closure takes them on the wall face.
void column_wall_gradients (double ustar, double thstar, double invL, double d, double z0, double kappa, double& dUdz,
                            double& dthdz, const StableCoefficients& sc = {});

//! bc.surface_stable_functions = qnse with nothing for it to act on (a neutral surface and no stratified
//! column): one message naming the key, else empty.
std::vector<std::string> check_stable_functions (const InitOptions& iopt, const SurfaceLayerOptions& surf, bool surface_layer);

//! The canopy as the column sees it (above): f per layer of the column's grid (the deck's layers), and the k terms'
//! options. Empty f: no canopy.
struct ColumnCanopy
{
    std::vector<double> f;
    CanopyTkeOptions tke;
};

//! The column's canopy from the level's field (Canopy::column_profile; collective). Aborts, naming
//! init.first_guess, unless the terrain is flat (flat_terrain: terrain.type = flat; on a terrain-following mesh
//! even a canopy of uniform height gives each column its own layer profile, and the column sits on the flat levels)
//! and the canopy is horizontally uniform: a relative spread of f over the columns above kColumnCanopySpread in any
//! layer is refused, naming canopy.type too.
ColumnCanopy column_canopy (const Canopy& canopy, const amrex::Geometry& geom, bool flat_terrain);

//! The relative spread of f over the columns the column guess accepts as uniform (a map's inverse-distance
//! weights of a uniform table reproduce it to rounding).
constexpr double kColumnCanopySpread = 1.0e-6;

//! The nodal heights of the deck's vertical grid over flat ground (the column's levels, nz + 1 values): the mesh of a
//! 4 x 4 column, read from its box and shared with every rank; ztop receives the top level.
std::vector<double> column_flat_levels (const Inputs& in, double& ztop);

//! Why the column cannot be formed from these inputs (empty: it can), each
//! message naming init.first_guess.
std::vector<std::string> check_column_guess (const InitOptions& iopt, const SteadyOptions& sopt, const ForcingOptions& fopt,
                                             const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer);

//! Solve the column (above) and return its profiles; aborts if it does not
//! converge within init.column_max_iterations. table_file, when not empty,
//! receives the table (z u v theta k, one row per line) on the I/O rank.
//! canopy: the column's canopy (column_canopy), or none (an empty f).
ColumnGuessResult run_column_guess (const Inputs& in, const InitOptions& iopt, const SteadyOptions& sopt,
                                    const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer,
                                    const ClosureOptions& copt, const ForcingOptions& fopt,
                                    const std::string& table_file = std::string(), const ColumnCanopy& canopy = ColumnCanopy());

} // namespace gust

#endif
