#ifndef GUST_REDUCTIONS_H_
#define GUST_REDUCTIONS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! Sum of component comp (cell-centred data only) over the valid region, identical to the last bit for
//! any number of MPI ranks and any DistributionMapping of a fixed BoxArray.
//! A different BoxArray (e.g. after a restart onto another layout) can change
//! the last bits; anything compared bitwise across runs must keep the layout.
//!
//! Each column (i, j) of a box is summed in k order by one GPU thread; the
//! column sums are added in storage order (i fastest, then j) on the host;
//! each box sum lands in its own slot of a vector of length nboxes, which is
//! all-reduced (every slot has one nonzero contributor, so adding the zeros
//! is exact) and then added in box-index order. The result does depend on
//! the BoxArray: a different box layout adds in a different order.
amrex::Real deterministic_sum (const amrex::MultiFab& mf, int comp);

//! Sum over the valid region of a(acomp) * b(bcomp), any index type, in
//! double. With an owner mask (MultiFab::OwnerMask) a point shared by two
//! boxes, or by a box and its periodic image, counts once; without one every
//! valid point of every box counts. Each column (i, j) of a box is added in k
//! order by one GPU thread, the columns in storage order on the host, the
//! boxes in box-index order after the all-reduce: the same bits for any
//! DistributionMapping of a fixed BoxArray (as deterministic_sum).
double deterministic_dot (const amrex::MultiFab& a, int acomp, const amrex::MultiFab& b, int bcomp,
                          const amrex::iMultiFab* owner);

//! One box's share of deterministic_dot: the sum over the box of mfi's
//! index (mfi.validbox() converted to a's index type) of a(acomp) * b(bcomp),
//! owned points only with an owner mask. One thread per (i, j) column adds
//! its k values in order, the host adds the columns in storage order: the
//! same bits for the same box on any rank. Local: no communication.
double box_dot (const amrex::MultiFab& a, int acomp, const amrex::MultiFab& b, int bcomp, const amrex::iMultiFab* owner,
                const amrex::MFIter& mfi);

//! The number of points of component 0 of the fields, over the valid region
//! and ngrow ghost layers, that are NaN, infinite, or not > 0 (allow_zero:
//! < 0). One fused pass per field and one all-reduce for the list, so every
//! rank gets the same count (call it on every rank). The validity check of
//! a field produced every iteration (the pseudo-time step, the SIMPLEC
//! denominators): it replaces a min and a max, two all-reduces, per field.
amrex::Long count_not_positive (const std::vector<const amrex::MultiFab*>& fields, bool allow_zero, int ngrow = 0);

//! Abort with "gust: <op>: <what> must be finite and > 0 everywhere" (allow_zero:
//! ">= 0") unless count_not_positive finds no point. Collective: every rank
//! aborts together.
void require_positive (const char* op, const std::string& what, const std::vector<const amrex::MultiFab*>& fields,
                       bool allow_zero, int ngrow = 0);

// GUST_DEBUG_CHECKS: set in Debug builds and builds with AMReX assertions
// (CMake GUST_TESTING, the CI Debug job). It guards the checks of an
// operator's input that its callers in the solver already made where the
// field is produced (the pseudo-time step in face_dtau, the wall distance at
// setup), so a release build does not repeat their all-reduces every
// iteration.
#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
#define GUST_DEBUG_CHECKS 1
#endif

//! The sum of component comp over each z level (the domain's k range, cell
//! data in z; any type in x and y) of the points the owner mask keeps: one
//! thread per (j, k) row adds its i values in order, the host adds the rows
//! in j order per box, the boxes in box-index order. sums[k - klo].
std::vector<double> deterministic_level_sums (const amrex::MultiFab& mf, int comp, const amrex::iMultiFab& owner,
                                              const amrex::Box& domain);

//! The owner mask of mf's layout (MultiFab::OwnerMask), built once per
//! layout and periodicity: the norms and the Coriolis mean need it every
//! iteration, on a layout that does not change. The last eight layouts are
//! kept (two levels of four index types), and freed at amrex::Finalize. The
//! reference is valid until the next call.
const amrex::iMultiFab& owner_mask (const amrex::MultiFab& mf, const amrex::Periodicity& period);

//! Maximum of component comp over the valid region, for any data type. NaN
//! if any valid value is NaN; zero is returned as +0. With those two cases
//! fixed, the maximum does not depend on the order of evaluation.
amrex::Real global_max (const amrex::MultiFab& mf, int comp);

} // namespace gust

#endif
