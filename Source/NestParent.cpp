#include "NestParent.H"

#include "Sha256.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace gust {

namespace {

// Position of x in units of dx from x0, and whether it is a whole number of
// cells to the alignment tolerance.
bool whole_cells (amrex::Real x, amrex::Real x0, double dx, int& n)
{
    const double c = (static_cast<double>(x) - static_cast<double>(x0)) / dx;
    if (!std::isfinite(c) || std::abs(c) > 1.0e9) { return false; }
    n = static_cast<int>(std::lround(c));
    return std::abs(c - n) <= 1.0e-4;
}

// The ghost rows below the ground and above the top (every box spans the
// column): the nearest valid row's value, over the horizontally grown box.
void zero_gradient_z (amrex::MultiFab& mf)
{
    const int ngz = mf.nGrowVect()[2];
    const amrex::IntVect ngh(mf.nGrowVect()[0], mf.nGrowVect()[1], 0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = amrex::grow(mfi.validbox(), ngh);
        const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
        const auto a = mf.array(mfi);
        amrex::Box below = vb, above = vb;
        below.setSmall(2, klo - ngz).setBig(2, klo - 1);
        above.setSmall(2, khi + 1).setBig(2, khi + ngz);
        amrex::ParallelFor(below, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = a(i, j, klo); });
        amrex::ParallelFor(above, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = a(i, j, khi); });
    }
}

} // namespace

std::vector<std::string> check_nest_placement (const CheckpointHeader& h, const MeshInputs& patch,
                                               const amrex::Vector<amrex::Real>& patch_zlevels, amrex::IntVect* offset)
{
    std::vector<std::string> err;
    const std::string pre = "bc.nest_parent: ";
    if (h.version < 2 || h.zlevels.empty()) {
        err.push_back(pre + "the parent checkpoint (format version " + std::to_string(h.version)
                      + ") has no z-levels; write it again with this version (checkpoint format version 2 or later)");
        return err;
    }
    // z: the same levels.
    if (h.n_cell[2] != patch.n_cell[2]) {
        err.push_back(pre + "the patch has " + std::to_string(patch.n_cell[2]) + " layers (mesh.n_cell z), the parent "
                      + std::to_string(h.n_cell[2]) + ": a nest keeps the parent's z-levels");
    } else if (h.prob_lo[2] != patch.prob_lo[2] || h.prob_hi[2] != patch.prob_hi[2]) {
        err.push_back(pre + "mesh.prob_lo z and mesh.prob_hi z differ from the parent's: a nest keeps the parent's z-levels");
    } else if (patch_zlevels.size() != h.zlevels.size()) {
        err.push_back(pre + "the patch has " + std::to_string(patch_zlevels.size()) + " z-levels, the parent "
                      + std::to_string(h.zlevels.size()));
    } else {
        for (int k = 0; k < static_cast<int>(h.zlevels.size()); ++k) {
            if (patch_zlevels[k] != h.zlevels[k]) {
                std::ostringstream os;
                os << pre << "z-level " << k << " is " << patch_zlevels[k] << " m in the patch and " << h.zlevels[k]
                   << " m in the parent (mesh.dz0, mesh.stretch_ratio, mesh.zlevels_file): a nest keeps the parent's z-levels";
                err.push_back(os.str());
                break;
            }
        }
    }
    if (patch.is_periodic[2] != 0 || h.is_periodic[2] != 0) {
        err.push_back(pre + "z is never periodic");
    }
    // x and y: ratio 2, aligned, inside.
    amrex::IntVect off(0);
    const char* xy = "xy";
    for (int d = 0; d < 2; ++d) {
        const std::string ax(1, xy[d]);
        const int N = h.n_cell[d];
        if (N < 1 || !(h.prob_hi[d] > h.prob_lo[d])) {
            err.push_back(pre + "the parent's mesh in " + ax + " is empty");
            continue;
        }
        // In double: the ends are exact (hex floats), a float dx is not.
        const double dx = (static_cast<double>(h.prob_hi[d]) - static_cast<double>(h.prob_lo[d])) / N;
        int o = 0, e = 0;
        const bool lo_ok = whole_cells(patch.prob_lo[d], h.prob_lo[d], dx, o);
        const bool hi_ok = whole_cells(patch.prob_hi[d], h.prob_lo[d], dx, e);
        if (!lo_ok || !hi_ok) {
            std::ostringstream os;
            os << pre << "mesh.prob_lo and mesh.prob_hi " << ax << " must lie on the parent's cell faces (cells of "
               << dx << " m from " << h.prob_lo[d] << " m)";
            err.push_back(os.str());
            continue;
        }
        const int n = e - o;
        if (n < 1) {
            err.push_back(pre + "the patch is empty in " + ax);
            continue;
        }
        if (patch.n_cell[d] != kNestRatio * n) {
            err.push_back(pre + "mesh.n_cell " + ax + " = " + std::to_string(patch.n_cell[d]) + ": the patch spans "
                          + std::to_string(n) + " parent cells, so it needs " + std::to_string(kNestRatio * n)
                          + " (refinement 2 2 1)");
            continue;
        }
        const bool pper = h.is_periodic[d] != 0;
        if (patch.is_periodic[d] != 0) {
            if (!pper || o != 0 || n != N) {
                err.push_back(pre + "mesh.is_periodic " + ax + " = 1 needs a parent periodic in " + ax
                              + " and a patch spanning its whole width");
                continue;
            }
        } else if (pper) {
            if (o < 0 || o + n > N) {
                err.push_back(pre + "the patch must lie inside the parent's domain in " + ax);
                continue;
            }
        } else if (o < kNestMargin || o + n > N - kNestMargin) {
            err.push_back(pre + "the patch must lie at least " + std::to_string(kNestMargin) + " parent cells inside the parent's "
                          + ax + " faces (parent cells " + std::to_string(kNestMargin) + " to " + std::to_string(N - kNestMargin)
                          + "; the patch spans " + std::to_string(o) + " to " + std::to_string(o + n) + ")");
            continue;
        }
        off[d] = o;
    }
    if (err.empty() && offset != nullptr) { *offset = off; }
    return err;
}

NestParent::NestParent (const std::string& dir, const MeshInputs& patch, const amrex::Vector<amrex::Real>& patch_zlevels,
                        const amrex::BoxArray& patch_ba, const amrex::DistributionMapping& patch_dm, bool with_ground)
    : m_dir(dir), m_patch_cells(patch.n_cell), m_patch_ba(patch_ba)
{
    int has_header = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        has_header = std::filesystem::is_regular_file(dir + "/Header", ec) ? 1 : 0;
    }
    amrex::ParallelDescriptor::Bcast(&has_header, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (has_header == 0) {
        amrex::Abort("gust: bc.nest_parent = " + dir + ": no checkpoint Header there");
    }
    const CheckpointHeader h = read_checkpoint_header(dir, "bc.nest_parent");
    const auto err = check_nest_placement(h, patch, patch_zlevels, &m_offset);
    if (!err.empty()) {
        std::string all = "gust: bc.nest_parent = " + dir + ":";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    m_step = h.step;
    m_parent_cells = h.n_cell;
    {
        // The digest of the parent's content (NestParent.H), read by the I/O rank.
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::string text;
            std::vector<std::string> files{dir + "/Header"};
            for (const char* f : {"u", "v", "w", "theta", "tke", "z_phys_nd"}) {
                files.push_back(amrex::MultiFabFileFullPrefix(0, dir, "Level_", f) + "_H");
            }
            for (const auto& name : files) {
                std::ifstream is(name, std::ios::binary);
                std::ostringstream os;
                os << is.rdbuf();
                text += name.substr(dir.size()) + "\n" + os.str() + "\n";
            }
            m_digest = sha256_hex(text);
        }
        long n = static_cast<long>(m_digest.size());
        const int root = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::ParallelDescriptor::Bcast(&n, 1, root);
        m_digest.resize(static_cast<std::size_t>(n));
        if (n > 0) { amrex::ParallelDescriptor::Bcast(m_digest.data(), static_cast<std::size_t>(n), root); }
    }

    MeshInputs pm;
    pm.n_cell = h.n_cell;
    pm.prob_lo = {h.prob_lo[0], h.prob_lo[1], h.prob_lo[2]};
    pm.prob_hi = {h.prob_hi[0], h.prob_hi[1], h.prob_hi[2]};
    pm.is_periodic = h.is_periodic;
    const amrex::Geometry pgeom = make_geometry(pm);
    const amrex::Box& pdom = pgeom.Domain();
    const amrex::Box fdom(amrex::IntVect(0), patch.n_cell - 1);

    // The parent cells under each patch box, in the parent's index space.
    const amrex::IntVect ratio(kNestRatio, kNestRatio, 1);
    amrex::BoxArray img = amrex::coarsen(patch_ba, ratio);
    img.shift(amrex::IntVect(m_offset[0], m_offset[1], 0));
    const amrex::DistributionMapping pdm(h.ba);

    struct Item { const char* name; amrex::MultiFab* mf; amrex::IndexType type; };
    const Item items[] = {
        {"u", &m_u, amrex::IndexType(amrex::IntVect(1, 0, 0))},
        {"v", &m_v, amrex::IndexType(amrex::IntVect(0, 1, 0))},
        {"w", &m_w, amrex::IndexType(amrex::IntVect(0, 0, 1))},
        {"p", &m_p, amrex::IndexType::TheCellType()},
        {"theta", &m_theta, amrex::IndexType::TheCellType()},
        {"tke", &m_tke, amrex::IndexType::TheCellType()},
    };
    for (const auto& it : items) {
        amrex::MultiFab src(amrex::convert(h.ba, it.type), pdm, 1, ngrow_state);
        amrex::VisMF::Read(src, amrex::MultiFabFileFullPrefix(0, dir, "Level_", it.name));
        it.mf->define(amrex::convert(img, it.type), patch_dm, 1, amrex::IntVect(kNestImageGhost, kNestImageGhost, 0));
        // Valid parent data only: the parent's own ghost cells are not read.
        it.mf->ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(kNestImageGhost, kNestImageGhost, 0), pgeom.periodicity());
    }

    if (with_ground) {
        // The node heights, as read_checkpoint reads them (their one ghost
        // layer is in the file); valid parent nodes only are copied.
        const amrex::IntVect nodes(1, 1, 1);
        amrex::MultiFab src(amrex::convert(h.ba, nodes), pdm, 1, ngrow_nodal);
        amrex::VisMF::Read(src, amrex::MultiFabFileFullPrefix(0, dir, "Level_", "z_phys_nd"));
        m_z.define(amrex::convert(img, nodes), patch_dm, 1, amrex::IntVect(kNestImageGhost, kNestImageGhost, 0));
        m_z.ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(kNestImageGhost, kNestImageGhost, 0), pgeom.periodicity());
    }

    const int r = kNestRatio;
    m_cx = transfer::make_axis(fdom.smallEnd(0), pdom, pgeom.isPeriodic(0), 0, r, false, m_offset[0]);
    m_cy = transfer::make_axis(fdom.smallEnd(1), pdom, pgeom.isPeriodic(1), 1, r, false, m_offset[1]);
    m_nx = transfer::make_axis(fdom.smallEnd(0), pdom, pgeom.isPeriodic(0), 0, r, true, m_offset[0]);
    m_ny = transfer::make_axis(fdom.smallEnd(1), pdom, pgeom.isPeriodic(1), 1, r, true, m_offset[1]);
}

void NestParent::interpolate (State& patch) const
{
    if (patch.p.boxArray() != m_patch_ba || patch.p.DistributionMap() != m_p.DistributionMap()) {
        amrex::Abort("gust: NestParent::interpolate: the patch state is not on the layout the parent was read for");
    }
    const amrex::IntVect ng(ngrow_state, ngrow_state, 0);
    transfer::interpolate<stencil::CC>(m_theta, patch.theta, m_cx, m_cy, ng);
    transfer::interpolate<stencil::CC>(m_tke, patch.tke, m_cx, m_cy, ng);
    transfer::interpolate<stencil::CC>(m_p, patch.p, m_cx, m_cy, ng);
    transfer::interpolate<stencil::XF>(m_u, patch.u, m_nx, m_cy, ng);
    transfer::interpolate<stencil::YF>(m_v, patch.v, m_cx, m_ny, ng);
    transfer::interpolate<stencil::ZF>(m_w, patch.w, m_cx, m_cy, ng);
    for (amrex::MultiFab* mf : {&patch.theta, &patch.tke, &patch.p, &patch.u, &patch.v, &patch.w}) {
        zero_gradient_z(*mf);
    }
}

void NestParent::interpolate_ground (amrex::MultiFab& hp) const
{
    if (m_z.empty()) { amrex::Abort("gust: NestParent::interpolate_ground: the parent was read without its node heights"); }
    if (hp.boxArray() != amrex::convert(m_patch_ba, amrex::IntVect(1, 1, 1)) || hp.DistributionMap() != m_z.DistributionMap()) {
        amrex::Abort("gust: NestParent::interpolate_ground: hp is not nodal on the layout the parent was read for");
    }
    const amrex::IntVect ng(std::min(hp.nGrowVect()[0], ngrow_state), std::min(hp.nGrowVect()[1], ngrow_state), 0);
    transfer::interpolate<stencil::XF | stencil::YF>(m_z, hp, m_nx, m_ny, ng);
}

std::string NestParent::summary () const
{
    std::ostringstream os;
    os << "gust: nest: parent " << m_dir << " (step " << m_step << ", " << m_parent_cells[0] << " x " << m_parent_cells[1]
       << " x " << m_parent_cells[2] << " cells); the patch's " << m_patch_cells[0] << " x " << m_patch_cells[1]
       << " cells cover parent cells " << m_offset[0] << " to " << m_offset[0] + m_patch_cells[0] / kNestRatio - 1 << " in x and "
       << m_offset[1] << " to " << m_offset[1] + m_patch_cells[1] / kNestRatio - 1 << " in y, refinement 2 2 1\n";
    return os.str();
}

} // namespace gust
