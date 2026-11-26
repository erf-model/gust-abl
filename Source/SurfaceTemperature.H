#ifndef GUST_SURFACE_TEMPERATURE_H_
#define GUST_SURFACE_TEMPERATURE_H_

// The ground's potential temperature theta_s for bc.surface_stability =
// surface_temperature (SurfaceLayer.H): one value (bc.surface_temperature) or a
// map theta_s(x, y) read from a file (bc.surface_temperature_file).
//
// theta_s is in the state's theta: the virtual potential temperature in a
// moist run (tools/hrrr_outer.py writes theta_v,s = T_skin (p0/p_s)^(R/cp)
// (1 + 0.61 q_2m)).
//
// The map lives on the roughness map's layout (SurfaceRoughness.H): one value
// per ground column on the ground slab of the level's cells, one ghost column
// in x and y. Level 0 reads "x y theta_s" points (terrain.file's xyz format)
// and takes, at each ground cell centre, the inverse-distance weighted mean of
// theta_s over terrain.idw_neighbors points, weight d^-terrain.idw_power (an
// arithmetic mean: theta_s enters the flux linearly). Every level reads the
// file at its own cell centres (the coarse levels only seed level 0). The ghost
// columns take the value at the centre of the column they image (periodic
// images, else the nearest domain column).
//
// The map of the boundary layer depth z_i (bc.surface_pblh_file, read_pblh_map) is read the same way.

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>

namespace gust {

//! The range every theta_s must lie in (K): the value, the file's points.
constexpr double kSurfaceTemperatureMin = 150;
constexpr double kSurfaceTemperatureMax = 350;
//! The largest bc.surface_kb_inverse = ln(z0 / z0t): z0t = z0 e^-20, far below any measured kB^-1 (about 1 to 13).
constexpr double kKbInverseMax = 20;

//! "[150, 350] K", for the messages that name the range.
inline std::string surface_temperature_range ()
{
    return "[" + std::to_string(static_cast<int>(kSurfaceTemperatureMin)) + ", "
           + std::to_string(static_cast<int>(kSurfaceTemperatureMax)) + "] K";
}

//! Aborts (naming op and what) unless map is one component on the ground slab of like's cells, on like's
//! DistributionMapping, with a ghost column in x and y.
void check_surface_temperature_map (const char* op, const amrex::MultiFab& map, const amrex::MultiFab& like,
                                    const char* what = "surface temperature map");

//! Level 0's map from file (key names the input in messages) on the ground slab of cell_ba, with
//! terrain.idw_neighbors k and terrain.idw_power. Aborts, naming key, on a file without points, a value that is
//! not finite or outside [kSurfaceTemperatureMin, kSurfaceTemperatureMax], or fewer points than k.
void read_surface_temperature_map (const std::string& file, const std::string& key, int k, double power,
                                   const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                                   const amrex::DistributionMapping& dm, amrex::MultiFab& theta_s);

//! Level 0's map of the boundary layer depth z_i (m) for bc.surface_wstar = beljaars (SurfaceWstar_K.H) from file
//! (bc.surface_pblh_file, "x y z_i" points), read as the temperature map; values in [kPblhMin, kPblhMax].
void read_pblh_map (const std::string& file, const std::string& key, int k, double power, const amrex::Geometry& geom,
                    const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm, amrex::MultiFab& zi);

//! "min a, max b, mean c <unit>" over the valid columns (every rank).
std::string surface_temperature_summary (const amrex::MultiFab& theta_s, const char* unit = "K");

} // namespace gust

#endif
