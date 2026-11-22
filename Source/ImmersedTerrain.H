#ifndef GUST_IMMERSED_TERRAIN_H_
#define GUST_IMMERSED_TERRAIN_H_

// The ground as an immersed surface on a flat mesh (terrain.representation = immersed_forcing;
// docs/immersed_terrain.rst). The mesh does not follow the ground: the cells below
// terrain.flat_height are solid, the cut cell carries ERF's form-7 wall law (BuildingForcing.H with
// an immersed ground) and the first full cell above it repeats the flat ground's shear
// (ImmersedStrain.H). The run is neutral, laterally periodic or with inflow and outflow faces.

#include "Buildings.H"
#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! What check_immersed_terrain needs from the other namespaces' inputs.
struct ImmersedTerrainFacts
{
    bool surface_layer{false};       //!< bc.zlo = surface_layer
    bool neutral{true};              //!< bc.surface_stability = neutral
    bool z0_map{false};              //!< bc.surface_z0_file given
    bool buildings{false};           //!< buildings.type given
    bool nests{false};               //!< nest.names given, or this run is a nest
    bool newton{false};              //!< steady.stage = simplec_newton
    bool zline_sweeps{true};         //!< steady.sweep = zjacobi
    bool poisson_mlmg{false};        //!< solver.poisson_backend = mlmg
    bool laminar{false};             //!< turb.closure = none
    bool unsteady{false};            //!< steady.stage = unsteady
    bool canopy{false};              //!< canopy.type given
    bool stratified{false};          //!< bc.theta_profile holds more than one value of theta
    bool column_guess{false};        //!< init.first_guess = column
    bool probes{false};              //!< steady.probes_file given
    bool zref_given{false};          //!< bc.surface_zref given
    bool floor_given{false};         //!< bc.surface_denominator_floor given
    bool wall_k_source{false};       //!< bc.surface_wall_k = source
    bool outflow_from_velocity{false};  //!< bc.outflow_transport = from_velocity
};

//! The refusals of terrain.representation = immersed_forcing that need other namespaces
//! (Inputs.cpp's validate has the terrain.* ones), with or without a stage: a run without one writes
//! the ground to its plotfile and checkpoint. Empty when the run is supported.
std::vector<std::string> check_immersed_terrain (const ImmersedTerrainFacts& f);

//! The refusal that needs the mesh's levels: the ground must lie at or above the top of the first
//! layer (first_layer_top, m), so that the cut cell is not the body-fitted ground's first cell, whose
//! k row the ground's own wall value shares (BuildingTke.H), and at or below the bottom of the top
//! layer (top_layer_bottom, m), which must hold fluid. Empty when it does.
std::string check_immersed_terrain_height (const TerrainInputs& t, double first_layer_top, double top_layer_bottom);

//! The ground's height into z_nd's last component (of two), at every node of every column, ghost
//! nodes included: terrain.flat_height plus the terrain's height (fill_terrain_surface for
//! terrain.type: 0 for flat). Component 0, the mesh's own heights, is not touched.
void fill_immersed_ground (amrex::MultiFab& z_nd, const amrex::Geometry& geom, const TerrainInputs& t);

//! The same refusal for the whole surface (z_nd's last component, valid nodes; a collective call).
//! Empty when the ground lies between the two levels.
std::string check_immersed_ground_extent (const amrex::MultiFab& z_nd, double first_layer_top, double top_layer_bottom);

//! The wall-law options of the immersed ground: buildings.type = immersed_forcing with the form 7
//! law on every wall (fraction_stress, open_faces), the wall value of k (wall_tke), no heat
//! (a neutral run), no lee zone, t.small_volfrac, the ground's roughness z0 (bc.surface_z0) and
//! immersed_ground set. t.representation must be immersed_forcing.
BuildingsOptions immersed_terrain_options (const TerrainInputs& t, double z0);

//! The ground as a height map in memory: one tile of height t.flat_height whose edges are the
//! domain's x and y faces (with a terrain file: the terrain's heights).
BuildingHeightMap immersed_terrain_map (const TerrainInputs& t, const amrex::Geometry& geom);

//! The SHA-256 of what fixes the immersed ground: flat_height, small_volfrac, z0, the terrain's type and what
//! the type reads (the hill's size and shape; the file's contents and the IDW choices). Reads terrain.file on the
//! calling rank.
std::string immersed_terrain_hash (const TerrainInputs& t, double z0);

//! <dir>/ImmersedTerrain with hash (I/O rank); nothing with an empty hash (a body-fitted run).
void write_immersed_terrain_checkpoint (const std::string& dir, const std::string& hash);

//! Refuses (abort naming io.restart) a restart whose checkpoint dir holds another hash, one when
//! hash is empty (the checkpoint's ground is immersed, the run's body-fitted), or none when hash
//! is not empty (the reverse).
void check_immersed_terrain_restart (const std::string& dir, const std::string& hash);

//! The file name in a checkpoint directory.
inline std::string immersed_terrain_checkpoint_file (const std::string& dir) { return dir + "/ImmersedTerrain"; }

}   // namespace gust

#endif
