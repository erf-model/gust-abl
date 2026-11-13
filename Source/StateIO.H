#ifndef GUST_STATEIO_H_
#define GUST_STATEIO_H_

// Named StateIO.H, not IO.H: Windows file systems ignore case, so an IO.H on
// the include path would stand in for the system <io.h>.

#include "State.H"

#include <array>
#include <string>
#include <vector>

namespace gust {

//! Plotfile variable names, in order. Velocities are face averages to cell
//! centres; z is the mean of the cell's eight node heights.
const amrex::Vector<std::string>& plot_var_names ();

//! Write <prefix>NNNNN (step from state). Besides the cell-centred variables,
//! the nodal heights go to <plotfile>/Level_0/z_phys_nd as a VisMF MultiFab,
//! ghost nodes included. With detJ (the mesh's cell volume factor), it is
//! added as an eighth variable "detJ". With canopy_drag (the canopy's drag
//! field f, Canopy.H) it is added after them as "canopy_drag", and with
//! canopy_lad (the leaf area density, with canopy.heat on) as "canopy_lad".
//! With buildings (Buildings.H), their solid fraction follows as "building_beta"
//! and component 0 of building_wall, the wall mask, as "building_wall".
//! time: the plotfile's time (the physical time of an unsteady run, else 0).
//! extra, with extra_names (one per component): cell-centred fields added last
//! (an unsteady run's averages, UnsteadyOutput.H).
//! Returns the directory name.
std::string write_plotfile (const std::string& prefix, const State& state,
                            const amrex::Geometry& geom, const amrex::MultiFab* detJ = nullptr,
                            const amrex::MultiFab* canopy_drag = nullptr, const amrex::MultiFab* canopy_lad = nullptr,
                            const amrex::MultiFab* building_beta = nullptr, const amrex::MultiFab* building_wall = nullptr,
                            amrex::Real time = 0, const amrex::MultiFab* extra = nullptr,
                            const std::vector<std::string>* extra_names = nullptr);

//! Write <prefix>NNNNN: a text Header (format version, bytes per real, step,
//! mesh in hex floats, the reference z-levels in hex floats, BoxArray) and
//! every State field as a VisMF MultiFab with ghosts. zlevels: the run's
//! reference z-levels (Mesh::zlevels(), nz + 1 values); when null they are
//! made from mesh (make_zlevels). A nested patch (NestParent.H) checks its
//! own z-levels against these.
//! Returns the directory name.
std::string write_checkpoint (const std::string& prefix, State& state,
                              const MeshInputs& mesh, const amrex::Vector<amrex::Real>* zlevels = nullptr);

//! The checkpoint format version write_checkpoint writes (2: with the
//! z-levels). read_checkpoint also reads version 1, which has none.
inline constexpr int kCheckpointFormatVersion = 2;

//! A checkpoint's Header, as read_checkpoint_header returns it.
struct CheckpointHeader
{
    int version{0};
    int step{0};
    amrex::IntVect n_cell{0, 0, 0};
    std::array<amrex::Real, 3> prob_lo{}, prob_hi{};
    amrex::Array<int, 3> is_periodic{0, 0, 0};
    //! The reference z-levels (nz + 1 values); empty in a version-1 checkpoint.
    amrex::Vector<amrex::Real> zlevels;
    amrex::BoxArray ba;   //!< cell-centred, as written
};

//! Read and check <dir>/Header (read on the I/O rank and broadcast). Aborts
//! with "gust: <key> = <dir>: ..." if it is missing, corrupt, of an unknown
//! format version, or was written with another precision.
CheckpointHeader read_checkpoint_header (const std::string& dir, const std::string& key);

//! Read a checkpoint into a State already defined on the current layout.
//! Aborts, naming io.restart, if the Header is corrupt, was written with
//! another precision, or its mesh differs from mesh. Valid values are
//! restored exactly; of the ghost cells, only z_phys_nd's are restored. The
//! checkpoint may come from a different box layout or rank count.
void read_checkpoint (const std::string& dir, State& state, const MeshInputs& mesh);

//! True when the checkpoint directory a write would create, target, is the
//! existing directory restart, however either path is spelled (relative or
//! absolute, "./", trailing or doubled separators, symbolic links, letter
//! case on a case-insensitive file system): std::filesystem::equivalent on
//! the two paths. False when either does not exist or cannot be examined.
//! Local (no communication); call it on the I/O rank and broadcast.
bool same_directory (const std::string& target, const std::string& restart);

} // namespace gust

#endif
