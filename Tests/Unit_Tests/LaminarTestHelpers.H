#ifndef GUST_LAMINAR_TEST_HELPERS_H_
#define GUST_LAMINAR_TEST_HELPERS_H_

// Laminar test cases on a flat mesh: the mesh, a state on it, and fixed
// uniform diffusivities with the boundary conditions and residual inputs of a
// case without a closure or wall model (test_Simplec.cpp, test_TimeDriver.cpp).

#include "BoundaryConditions.H"
#include "Poisson.H"
#include "Residual.H"
#include "State.H"
#include "TransportTestHelpers.H"

#include <AMReX_MultiFab.H>

#include <functional>
#include <memory>

namespace gust::lamtest {

using gust::FaceType;
using gust::ptest::PCase;

// A flat mesh with uniform layers: nx x ny x nz cells on Lx x Ly x H; with
// ground, the layers follow that ground height h(x, y) instead.
inline PCase flat_case (int nx, int ny, int nz, double Lx, double Ly, double H, bool periodic_x, bool periodic_y, int max_grid = 64,
                        const std::function<double(double, double)>& ground = {})
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(Lx), static_cast<amrex::Real>(Ly), static_cast<amrex::Real>(H)};
    c.in.mesh.is_periodic = {periodic_x ? 1 : 0, periodic_y ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    if (ground) {
        c.make(ground);
    } else {
        c.make([] (double, double) { return 0.0; });
    }
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

inline void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// Fixed uniform diffusivities, the boundary conditions and the residual
// inputs of a laminar case (no closure, no wall model).
struct Laminar
{
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::ResidualInputs in;
    gust::PoissonOptions popt;
    void make (const PCase& c, double visc, FaceType ground, FaceType top, FaceType xlo = FaceType::periodic,
               FaceType xhi = FaceType::periodic)
    {
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(static_cast<amrex::Real>(visc));
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = 0;
        bin.face = {xlo, c.in.mesh.is_periodic[1] ? FaceType::periodic : FaceType::slip_wall, ground,
                    xhi, c.in.mesh.is_periodic[1] ? FaceType::periodic : FaceType::slip_wall, top};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        popt.reltol = amrex::Real(0.1);
    }
};

} // namespace gust::lamtest

#endif
