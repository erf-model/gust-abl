#ifndef GUST_JOINT_RECORD_H_
#define GUST_JOINT_RECORD_H_

// The joint run's record (JointRun.H), written into the outer domain's last
// checkpoint after the joint run: how many exchanges were taken, and for the
// outer domain and every nest the CFL floor the joint run holds its driver to
// (Driver::floor_cfl_at_current) and, for a nest, the checkpoint it wrote at
// the same exchange. A restart from that checkpoint (io.restart with
// nest.feedback = defect, relax, relax_consistent or inject) reads it, restarts every nest from its
// checkpoint and goes on with the exchanges (Run.cpp).
//
// <dir>/Joint, text:
//   GUST_JOINT 3
//   feedback <defect|relax|relax_consistent|inject>
//   exchanges <n>
//   outer <cfl floor, hex float>
//   window <m> <m hex floats>
//   factor <nest.feedback_factor, hex float; 0 with defect>
//   nest <name> <cfl floor, hex float> <checkpoint, the rest of the line>
//   ... one nest line per nest, in nest.names order
//
// The window line (version 2) holds the joint stop test's last m <= steady.monitor_window changes of the feedback,
// oldest first (FeedbackWindow, NestDefect.H), and beside the record the checkpoint holds the feedback of the last
// exchange, Level_0/joint_feedback_{u,v,w,theta,tke} (S with defect, the target q* with relax and relax_consistent,
// the injection's target with inject). A restart restores
// both, so its first exchange measures the change and tests the window as the run without a restart does. The factor
// line (version 3) holds nest.feedback_factor: the window's changes and the feedback were made with it, so a restart
// with another factor (or another nest.feedback) starts with an empty window (joint_restore_refusal). A version 2
// record (no factor line) restores both with the deck's factor; a version 1 record (no window line, no field) still
// reads, and its restart starts with an empty window.

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
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
    std::string feedback;           //!< "defect", "relax", "relax_consistent" or "inject"
    int exchanges{0};               //!< exchanges taken over the run, restarts included
    amrex::Real outer_cfl_floor{0};
    std::vector<Nest> nests;
    //! Version 2: the window line was read (and the feedback field is in the checkpoint); false for version 1.
    bool has_feedback{false};
    std::vector<amrex::Real> window;   //!< the stop test's changes, oldest first
    //! Version 3: the factor line was read; false for versions 1 and 2.
    bool has_factor{false};
    amrex::Real factor{0};             //!< nest.feedback_factor (0 with defect)
};

//! <dir>/Joint (the I/O rank writes; every rank waits). Aborts, naming the file, if it cannot be written.
void write_joint_record (const std::string& dir, const JointRecord& r);

//! Reads <dir>/Joint on the I/O rank and gives every rank its contents. Returns false when the
//! checkpoint has no Joint file (a one-way checkpoint); aborts, naming the file, when it cannot be read.
bool read_joint_record (const std::string& dir, JointRecord& r);

//! Empty when a restart with nest.feedback = feedback and nest.feedback_factor = factor (0 with defect) restores the
//! record's window and last feedback; otherwise why it does not (the start of a log line).
std::string joint_restore_refusal (const JointRecord& prior, const std::string& feedback, amrex::Real factor);

//! The feedback of the last exchange into <dir>/Level_0/joint_feedback_{u,v,w,theta,tke} (valid points).
void write_joint_feedback (const std::string& dir, const std::array<const amrex::MultiFab*, 5>& f);

//! Reads the fields write_joint_feedback wrote onto the layouts of like (out defined here, no ghosts). Aborts,
//! naming io.restart and the field, when one is missing or has another index type or component count.
void read_joint_feedback (const std::string& dir, const std::array<const amrex::MultiFab*, 5>& like,
                          std::array<amrex::MultiFab, 5>& out);

} // namespace gust

#endif
