#ifndef GUST_NEST_BUILDINGS_H_
#define GUST_NEST_BUILDINGS_H_

// Buildings inside a nested patch (NestFaces.H, Buildings.H). One way: the
// patch takes its nest faces from the parent and its buildings from
// buildings.file on its own cells, as a whole-domain run on the same cells
// would. Nothing at a nest face knows about buildings, so a building must keep
// clear of every nest face, of the relaxation zone and of the ground blend:
//
//   gap_d = kNestBuildingCells dx_d + max(bc.nest_relax_width dx_d, terrain.nest_blend_width)
//
// for the nest faces normal to d (dx_d the patch's cell size along d). A
// building tile (BuildingHeightMap: the rectangle a point owns, with a height
// above 0) that reaches into the strip of width gap_d inside a nest face
// stops the run, naming the tile and the face. kNestBuildingCells = 2: a
// building's wall cells and their reference cells (Buildings.H: d2 is the next
// cell out) then lie inside the patch and outside the zone and the blend, and
// no fixed row of a building meets the nest face's Dirichlet row or the two
// ghost layers the parent fills.
//
// Along a nest face the tile counts where it overlaps the patch: in a
// non-periodic direction the tile's extent against the domain's, in a periodic
// one always (an image of the tile lies in the domain; Buildings::build refuses
// tiles more than a period away). The check is on the tiles, not on the cells'
// fractions, so it does not depend on buildings.small_volfrac, the box layout or
// the rank count; the map is the same on every rank (read_building_map).
//
// A tile edge on the strip's inner edge is clear, to the larger of 1e-9 of a
// cell and 64 Real round-offs of the coordinates (the allowance of
// Buildings::build; in single precision the geometry's own round-off).

#include "BoundaryConditions.H"
#include "Buildings.H"

#include <AMReX_Geometry.H>

#include <string>
#include <vector>

namespace gust {

//! The patch cells between a nest face and the nearest building, at least.
inline constexpr int kNestBuildingCells = 2;

//! One message per building tile within the gap of a nest face of bin (none
//! without a nest face): "buildings.file = <file>: the building tile x in
//! [a, b] m, y in [c, d] m (top h m) lies within g m of the nest face bc.xlo
//! (...)". At most max_messages tiles are named, then a count of the rest.
std::vector<std::string> check_nest_buildings (const BuildingHeightMap& map, const std::string& file,
                                               const BoundaryInputs& bin, const amrex::Geometry& geom,
                                               amrex::Real blend_width, int max_messages = 8);

} // namespace gust

#endif
