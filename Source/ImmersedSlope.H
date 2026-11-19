#ifndef GUST_IMMERSED_SLOPE_H_
#define GUST_IMMERSED_SLOPE_H_

// terrain.immersed_slope = normal: the wall law of an immersed ground (BuildingForcing.H) with the
// slope's corrections (docs/immersed_terrain.rst, Sloping ground). With the ground's mean slope
// (h_x, h_y) over a cut cell's footprint (Buildings::surface()),
//     cos a = 1 / sqrt(1 + h_x^2 + h_y^2),   n = (-h_x, -h_y, 1) cos a,
//   - the distances are those along the normal, d cos a (d the fluid's mean height above the ground);
//   - the speeds are those parallel to the ground, |u - (u . n) n|, w included;
//   - the wall's area is the ground's, the footprint's part times 1 / cos a;
//   - the stress acts on w too, in the cells whose ground slopes.
// The default, vertical, is ERF's form: the heights, the horizontal speed and the footprint's area
// (roof_rates in BuildingForcing.cpp). On a ground without slope the two are bit-identical.

#include "Buildings.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>

namespace gust {

//! The rates of the ground's wall law with the slope's corrections, into rate's components 0, 1, 2
//! (u, v, w), kBuildingUstar2 and the heat components, on the valid cells. b must be built by
//! build_from_surface. weight: the faces' weights in a cell's speed (one ghost layer).
void immersed_slope_rates (const Buildings& b, const std::array<amrex::MultiFab, 3>& weight, const amrex::MultiFab& detJ,
                           const amrex::Geometry& geom, const amrex::MultiFab& u, const amrex::MultiFab& v,
                           const amrex::MultiFab& w, amrex::Real z0, amrex::Real kappa, amrex::Real umin, amrex::MultiFab& rate);

}   // namespace gust

#endif
