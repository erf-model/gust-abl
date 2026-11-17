#ifndef GUST_JOINT_RECORD_H_
#define GUST_JOINT_RECORD_H_

// The joint run's record (JointRun.H), written into the outer domain's last
// checkpoint after the joint run: how many exchanges were taken, and for the
// outer domain and every nest the CFL floor the joint run holds its driver to
// (Driver::floor_cfl_at_current) and, for a nest, the checkpoint it wrote at
// the same exchange. A restart from that checkpoint (io.restart with
// nest.feedback = defect or relax) reads it, restarts every nest from its
// checkpoint and goes on with the exchanges (Run.cpp).
//
// <dir>/Joint, text:
//   GUST_JOINT 1
//   feedback <defect|relax>
//   exchanges <n>
//   outer <cfl floor, hex float>
//   nest <name> <cfl floor, hex float> <checkpoint, the rest of the line>
//   ... one nest line per nest, in nest.names order

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace gust {

struct JointRecord
{
    struct Nest
    {
        std::string name;           //!< as in nest.names
        amrex::Real cfl_floor{0};
        std::string chk;            //!< the nest's checkpoint at the same exchange
    };
    std::string feedback;           //!< "defect" or "relax"
    int exchanges{0};               //!< exchanges taken over the run, restarts included
    amrex::Real outer_cfl_floor{0};
    std::vector<Nest> nests;
};

//! <dir>/Joint (the I/O rank writes; every rank waits). Aborts, naming the file, if it cannot be written.
void write_joint_record (const std::string& dir, const JointRecord& r);

//! Reads <dir>/Joint on the I/O rank and gives every rank its contents. Returns false when the
//! checkpoint has no Joint file (a one-way checkpoint); aborts, naming the file, when it cannot be read.
bool read_joint_record (const std::string& dir, JointRecord& r);

} // namespace gust

#endif
