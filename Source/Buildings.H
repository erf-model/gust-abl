#ifndef GUST_BUILDINGS_H_
#define GUST_BUILDINGS_H_

// The geometry of buildings on the terrain, from the buildings.* inputs
// (docs/buildings.rst): solid fractions and the wall cells of the form 7
// immersed wall law (ERF erf.if_wall_form = fraction_stress, Kynema-SGF
// wall_form 7). The ground keeps its body-fitted wall law and is not blanked.
// BuildingForcing.H has how it acts on momentum.
//
// buildings.file, heights above the local ground, buildings.format:
//   erf  ERF's text format (erf.buildings_file_name): nx, ny, the nx x values,
//        the ny y values, then the heights with x the slow index, h[ix ny + iy];
//   xyz  "x y height" or "x y height zstart" lines, comma or white-space
//        separated, '#' comments, holding every point of a rectangular grid once,
//        in any order (no interpolation: it would smear the walls). The building
//        is solid from zstart (default 0) to height, so zstart > 0 leaves a void
//        deck under it; height is the top, not the thickness.
// On a slope both heights are measured from the ground of the column, so
// a building follows the terrain. Each point owns a tile reaching halfway to
// its neighbours (the end tiles as far outward as inward), the height constant
// over it and 0 outside the tiles. A wall thus lies on a tile edge and is
// vertical, and a cell's fraction is an exact sum over the tiles. ERF interpolates
// the same points to the cell corners, a ramp over an edge cell; on points at every
// cell node the columns hold the same volume, in other layers.
//
// Fractions (solid part, after buildings.small_volfrac at both ends):
//   beta   each cell's, over its footprint and its heights above the ground of
//          its column (four-corner means, as in detJ);
//   x, y   a vertical face's, along it, its heights from its two-node means; on
//          a tile edge the taller tile counts;
//   z      a horizontal face's, the part under a building or at its roof.
// A face touching a solid cell then has fraction 1, and one between two fluid
// cells 0, so faces and cells agree on a slope, where the means differ.
//
// Wall cells: from the blocked fractions of opposite faces, whose difference is
// the wall area in the cell facing that way (Buildings_K.H); in a fluid cell any
// blocked face is a wall. Along z this gives ERF ib_wall::is_wall_cell's partial
// cell under fluid and fluid cell on solid. A roof and a side make an edge cell,
// two sides a corner. A building thinner than a cell has no net wall there; the log
// warns. Distances along the normal, Delta_n = dx, dy or the cell's height, with
// beta = 0 in a fluid cell:
//   d1 = (1 - beta) Delta_n / 2 (the fluid centroid),
//   d2 = (1 - beta) Delta_n + Delta_ref / 2 (the reference cell, the next one out;
//        (3/2 - beta) Delta_n on equal layers), 0 when that cell is solid or
//        outside the domain.
//
// Ghost cells of beta and wall: one layer, periodic images in a periodic
// direction, else the nearest domain column, so a domain face is no wall;
// below and above the domain beta is 0.
//
// Restart: the checkpoint holds the SHA-256 of the inputs; a restart with
// others, or with buildings where the checkpoint has none or the reverse, is
// refused.

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace gust {

//! The components of Buildings::wall().
enum BuildingWallComp { kWallMask = 0, kWallD1x, kWallD2x, kWallD1y, kWallD2y, kWallD1z, kWallD2z, kNumWallComps };

//! The bits of the mask: a roof below the cell, a side at its low or high x or y,
//! the underside of a lifted building above the cell.
enum BuildingWallBit { kWallRoof = 1, kWallXlo = 2, kWallXhi = 4, kWallYlo = 8, kWallYhi = 16, kWallUnder = 32 };

enum class BuildingsType { none, immersed_forcing };
//! buildings.wall_form: the form 7 wall law on every wall, on roofs and undersides alone,
//! or the fixed rows alone.
enum class BuildingWallForm { fraction_stress, roof_stress, drag };

//! The buildings.* inputs (read_buildings_options).
struct BuildingsOptions
{
    BuildingsType type{BuildingsType::none};
    std::string type_name{"none"};   //!< buildings.type as given
    std::string file;                //!< buildings.file
    std::string format{"erf"};       //!< buildings.format: erf or xyz
    double small_volfrac{0.005};     //!< buildings.small_volfrac, in [0, 0.5)
    double z0{0.1};                  //!< buildings.z0 (m), > 0: the roughness of the wall law
    std::string wall_form_name{"fraction_stress"};   //!< buildings.wall_form as given
    BuildingWallForm wall_form{BuildingWallForm::fraction_stress};
    //! buildings.wall_tke as given (empty: not given) and its value: on by default with
    //! fraction_stress, off with the other forms (BuildingTke.H).
    std::string wall_tke_name;
    bool wall_tke{true};
    std::vector<std::string> given;  //!< buildings.* keys in the inputs
    std::vector<std::string> parse_errors;   //!< value counts and repeats
    [[nodiscard]] bool active () const { return type != BuildingsType::none; }
};

//! Read buildings.*. No range checks.
BuildingsOptions read_buildings_options ();

//! Range checks: one message per violation, each naming the input.
std::vector<std::string> check_buildings_options (const BuildingsOptions& opt);

//! A height map as tiles: edges xe (nx + 1), ye (ny + 1), tops h[i ny + j] and
//! starts zs[i ny + j] of the solid.
struct BuildingHeightMap
{
    int nx{0}, ny{0};
    std::vector<double> xe, ye, h, zs;
};

//! The map of a file in format erf or xyz, read on the I/O rank and broadcast.
//! Aborts, naming buildings.file, on a value that is not a number, nx or ny
//! below 2, coordinates not strictly increasing, a height not finite and >= 0;
//! erf: a number of values other than 2 + nx + ny + nx ny; xyz: a line not of 3
//! or 4 numbers, a zstart < 0 or not below its height (where height > 0), or points
//! that are not each point of a rectangular grid exactly once.
BuildingHeightMap read_building_map (const std::string& file, const std::string& format = "erf");

//! A level's building geometry.
class Buildings
{
public:
    //! From the inputs (opt.active()), on the cells of z_nd's mesh. Aborts on
    //! tiles wider than a periodic domain or more than one period from it, on boxes not
    //! spanning the column, on a map with no building in the
    //! domain, and on a building reaching the top layer.
    void build (const BuildingsOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm);

    [[nodiscard]] bool active () const { return m_active; }
    //! The cells' solid fraction, one ghost layer.
    [[nodiscard]] const amrex::MultiFab& beta () const { return m_beta; }
    //! The solid fraction of the faces normal to d (0, 1, 2), no ghost faces.
    [[nodiscard]] const amrex::MultiFab& beta_face (int d) const { return m_beta_face[static_cast<std::size_t>(d)]; }
    //! The wall cells (BuildingWallComp), one ghost layer. d1 and d2 are 0
    //! except along a normal with a wall; d2 is 0 without a usable reference cell.
    [[nodiscard]] const amrex::MultiFab& wall () const { return m_wall; }

    //! The log lines: cells by kind, solid volume, and notes for a building on a
    //! non-periodic domain face and for partial cells with no wall. detJ on the mesh of beta().
    [[nodiscard]] std::string summary (const amrex::MultiFab& detJ, const amrex::Geometry& geom) const;

    //! The SHA-256 of the inputs (options and file contents).
    [[nodiscard]] const std::string& input_hash () const { return m_hash; }

    //! <dir>/Buildings with input_hash(), when active (I/O rank).
    void write_checkpoint (const std::string& dir) const;
    //! Refuses (abort naming io.restart) a restart whose checkpoint dir holds a
    //! hash other than hash (hash_inputs of the run's options; empty without
    //! buildings), a hash when hash is empty, or none when it is not.
    static void check_restart (const std::string& dir, const std::string& hash);
    static std::string hash_inputs (const BuildingsOptions& opt);

private:
    bool m_active{false};
    amrex::MultiFab m_beta;
    std::array<amrex::MultiFab, 3> m_beta_face;
    amrex::MultiFab m_wall;
    std::string m_hash;
};

//! The buildings file name in a checkpoint directory.
inline std::string buildings_checkpoint_file (const std::string& dir) { return dir + "/Buildings"; }

} // namespace gust

#endif
