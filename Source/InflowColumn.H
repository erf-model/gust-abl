#ifndef GUST_INFLOWCOLUMN_H_
#define GUST_INFLOWCOLUMN_H_

// The run's own column as its inflow profile (bc.inflow_profile = column).
//
// With init.first_guess = column the run solves its horizontally uniform column at start-up (ColumnGuess.H),
// including the target-wind search (init.column_target_u, _v, _height) that adjusts the geostrophic wind and the
// pressure gradient. With bc.inflow_profile = column the converged column is then also the profile of the inflow
// faces and of the damping layers (InflowProfile.H, the table kind), and the run keeps the adjusted forcing. This
// is what two runs did before: a periodic column run that writes its table (io.column_guess_file) and logs the
// forcing, then the 3-D run with bc.inflow_profile = file and that forcing copied into forcing.ug, forcing.vg.
//
// The faces take the column's roughness (column_roughness: init.column_z0 under a map, else bc.surface_z0), so the
// column and the faces cannot differ; bc.inflow_z0 is refused. The column's table always has k.
//
// The column starts from the analytic inflow defaults (bc.inflow_profile = log_law with its default u* and z0,
// along x), as a periodic column run with no bc.inflow_* keys does, so the single run's column is the two-run
// column to the last bit on the same levels.
//
// Restart: the checkpoint holds the table, the faces' z0 and the forcing the column ran with (a text file
// InflowColumn in the checkpoint directory, every value in %.17g, which reads back to the last bit in double and
// single precision). A restart with bc.inflow_profile = column takes the table and z0 from it and, when the column
// had a target wind, the adjusted forcing (forcing.ug, forcing.vg, forcing.pressure_gradient x and y); init.* keys
// are refused on a restart, so the column is not solved again.

#include "BoundaryConditions.H"
#include "ColumnGuess.H"
#include "Initializer.H"
#include "InflowProfile.H"
#include "SurfaceLayer.H"

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace gust {

//! What a restart needs to continue a run with bc.inflow_profile = column.
struct InflowColumnRecord
{
    InflowTable table;        //!< z u v theta k (ncol = 5)
    amrex::Real z0{0};        //!< the faces' roughness: the column's (m)
    bool target{false};       //!< the forcing below was adjusted to a target wind
    amrex::Real ug{0}, vg{0}, pressure_gradient_x{0}, pressure_gradient_y{0};   //!< the forcing the column ran with
};

//! The roughness the column runs on: init.column_z0 under a roughness map (bc.surface_z0_file), else bc.surface_z0.
amrex::Real column_roughness (const InitOptions& iopt, const SurfaceLayerOptions& surf);

//! Why bc.inflow_profile = column cannot be used with these inputs (empty: it can, or the kind is another), each
//! message naming the key. restart: io.restart is given; nests: the run has nest.names.
std::vector<std::string> check_inflow_column (const BoundaryInputs& bin, const InitOptions& iopt, bool restart, bool nests);

//! The record of a converged column (the table, the faces' z0 and the forcing it ran with).
InflowColumnRecord make_inflow_column_record (const ColumnGuessResult& cg, amrex::Real z0);

//! The faces' profile options from a record: the column kind with its table and z0 (the other options kept).
void use_inflow_column (const InflowColumnRecord& rec, InflowProfileOptions& opt);

//! The record as text, and back (where names the source in the messages; aborts on a bad or truncated record).
std::string format_inflow_column_record (const InflowColumnRecord& rec);
InflowColumnRecord parse_inflow_column_record (const std::string& text, const std::string& where);

//! The record's file in a checkpoint directory.
std::string inflow_column_checkpoint_file (const std::string& dir);
//! Write the record into the checkpoint directory dir (the I/O rank; collective: a barrier).
void write_inflow_column_checkpoint (const std::string& dir, const InflowColumnRecord& rec);
//! The record of io.restart = dir, read on the I/O rank and broadcast (collective). active: the restart deck has
//! bc.inflow_profile = column. Aborts, naming io.restart, when the checkpoint has no record and active is true, or
//! has one and active is false. Returns the record (empty when not active).
InflowColumnRecord read_inflow_column_checkpoint (const std::string& dir, bool active);
//! The record of io.restart = dir if it has one (empty otherwise), for a restart with no stage, which reads no bc.*
//! keys and carries the record into its own checkpoint (collective).
InflowColumnRecord read_inflow_column_checkpoint_if_present (const std::string& dir);

} // namespace gust

#endif
