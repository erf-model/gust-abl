#ifndef GUST_CANOPY_H_
#define GUST_CANOPY_H_

// The forest canopy's drag field f = Cd a(z) (1/m), cell-centred, from the
// canopy.* inputs (docs/canopy.rst).
//
// Sources (canopy.type):
//   patches  ERF's text file of round patches, one per line,
//            "type x y height diameter cd lai laimax" (ERF erf.forest_file):
//            a column is in a patch when its centre lies within diameter / 2
//            of (x, y), the periodic image nearest in a periodic direction;
//            overlapping patches add. type 1 is the uniform profile, 2 the
//            Lalic–Mihailovic profile with that line's laimax.
//   map      the canopy height, the leaf area index and (canopy.cd_file) the
//            drag coefficient as "x y value" point files, the format of
//            terrain.file = xyz and bc.surface_z0_file, each interpolated to
//            the column centres by inverse-distance weighting with
//            terrain.idw_neighbors and terrain.idw_power (TerrainSurface.H);
//            canopy.cd or canopy.cd_file, canopy.profile and canopy.laimax
//            apply to every column.
//
// Each cell takes the mean of Cd a(z) over its height above the column's
// ground (Canopy_K.H), the ground and the cell's faces being the means of
// the four corner nodes, as in the cell volume factor detJ. The column sum of
// f detJ dzeta is then Cd LAI to rounding, part-filled top cells included,
// for a canopy below the top of the domain.
//
// The field has one ghost column in x and y (no ghost layer in z, as the
// boxes span the column): periodic images in a periodic direction, else the
// nearest domain column (image_ghost_columns, SurfaceRoughness.H).
//
// Coarse levels of the grid hierarchy (mesh.hierarchy_levels > 1) take
//   f_c = (sum over the r x r fine cells of f detJ) / (r^2 detJ_c)
// in each layer, so the integral of f over each coarse cell's volume is the
// integral over the fine cells it covers (the coarse ground is smoothed, so
// the coarse cell's own volume can differ from theirs).
//
// Restart: the checkpoint holds the SHA-256 of the canopy inputs (the
// options and the contents of the files read); a restart with other inputs,
// or with a canopy where the checkpoint has none or the reverse, is refused.

#include "Mesh.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

enum class CanopyType { none, patches, map };
enum class CanopyProfile { uniform, lalic_mihailovic };
//! canopy.drag_diagonal: the momentum rows' implicit drag term (CanopyDrag.H).
enum class CanopyDiagonal { derivative, picard, none };

//! The canopy.* inputs (read_canopy_options) and the IDW keys they share.
struct CanopyOptions
{
    CanopyType type{CanopyType::none};
    std::string type_name{"none"};       //!< canopy.type as given
    std::string profile_name{"uniform"}; //!< canopy.profile as given
    std::string file;          //!< canopy.file (patches)
    std::string height_file;   //!< canopy.height_file (map)
    std::string lai_file;      //!< canopy.lai_file (map)
    std::string cd_file;       //!< canopy.cd_file (map; empty: canopy.cd)
    double cd{0};              //!< canopy.cd (map without canopy.cd_file)
    CanopyProfile profile{CanopyProfile::uniform};   //!< canopy.profile (map)
    double laimax{0.8};        //!< canopy.laimax (map with lalic_mihailovic)
    std::string drag_diagonal_name{"derivative"};   //!< canopy.drag_diagonal as given
    CanopyDiagonal drag_diagonal{CanopyDiagonal::derivative};   //!< canopy.drag_diagonal
    int idw_neighbors{6};      //!< terrain.idw_neighbors (map)
    double idw_power{2};       //!< terrain.idw_power (map)
    std::vector<std::string> given;   //!< canopy.* keys in the inputs
    std::vector<std::string> parse_errors;   //!< value counts and repeats
    [[nodiscard]] bool active () const { return type != CanopyType::none; }
};

//! Read canopy.* (and take the IDW keys given). No range checks.
CanopyOptions read_canopy_options (int idw_neighbors, double idw_power);

//! Range checks: one message per violation, each naming the input. Files are
//! checked for being readable here; their contents when the field is built.
std::vector<std::string> check_canopy_options (const CanopyOptions& opt);

//! One round patch of a canopy.file line.
struct CanopyPatch
{
    int type{1};
    double x{0}, y{0}, height{0}, diameter{0}, cd{0}, lai{0}, laimax{0.8};
};

//! The patches of an ERF forest file (key names it in messages), read on the
//! I/O rank and broadcast. Aborts, naming the file and line, on a line that
//! is not eight numbers, a type other than 1 or 2, a height or diameter not
//! finite and > 0, a cd or lai not finite and >= 0, a laimax outside [0, 1)
//! on a type 2 line (not finite on a type 1 line), or a file without patches.
std::vector<CanopyPatch> read_canopy_patches (const std::string& file, const std::string& key);

//! The level's drag field: allocated, filled and ghost columns imaged.
class Canopy
{
public:
    //! Level 0 from the inputs (opt.active()), on the cells of z_nd's mesh
    //! (cell_ba, dm; z_nd with its ghost nodes).
    void build (const CanopyOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm);

    //! A coarse level from level 0 (fine, with its mesh's detJ): see above.
    //! mesh is the coarse level's (its detJ on cell_ba, dm). Aborts unless
    //! the coarse domain coarsens the fine one by one ratio in x and y.
    void restrict_from (const Canopy& fine, const amrex::MultiFab& fine_detJ, const amrex::Geometry& geom,
                        const Mesh& mesh, const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm);

    [[nodiscard]] bool active () const { return m_active; }
    //! f (1/m), one ghost column in x and y; empty unless active().
    [[nodiscard]] const amrex::MultiFab& drag () const { return m_drag; }
    //! canopy.drag_diagonal (a coarse level takes level 0's).
    [[nodiscard]] CanopyDiagonal drag_diagonal () const { return m_diagonal; }

    //! The log line: the largest f, the columns with canopy, the mean column
    //! integral of f over them (their mean Cd LAI) and its largest value.
    [[nodiscard]] std::string summary (const amrex::MultiFab& detJ, const amrex::Geometry& geom) const;

    //! The SHA-256 of the inputs (options and file contents) the field came
    //! from; empty unless active().
    [[nodiscard]] const std::string& input_hash () const { return m_hash; }

    //! <dir>/Canopy with input_hash(), when active (I/O rank).
    void write_checkpoint (const std::string& dir) const;
    //! Refuses (abort naming io.restart) a restart whose checkpoint dir holds
    //! a canopy hash other than hash (hash_inputs of the run's options; empty
    //! without a canopy), a hash when hash is empty, or none when it is not.
    //! Static, so the run checks before it builds the field.
    static void check_restart (const std::string& dir, const std::string& hash);

    //! The SHA-256 of opt (its keys and the contents of its files).
    static std::string hash_inputs (const CanopyOptions& opt);

private:
    bool m_active{false};
    amrex::MultiFab m_drag;
    std::string m_hash;
    CanopyDiagonal m_diagonal{CanopyDiagonal::derivative};
};

//! The canopy file name in a checkpoint directory.
inline std::string canopy_checkpoint_file (const std::string& dir) { return dir + "/Canopy"; }

} // namespace gust

#endif
