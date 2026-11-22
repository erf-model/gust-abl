#include "HybridTerrain.H"
#include "State.H"
#include "ZLevels.H"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <deque>
#include <limits>
#include <sstream>

namespace gust {

namespace {
// The tolerance of a width's node count, floor(width / d + tol): 1e-9, or more where amrex::Real is float.
double node_count_tol () { return std::max(1e-9, 64.0 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon())); }
}

namespace {

//! The node of a window index: wrapped on a periodic line of n nodes, else the index itself (the window
//! is cut at the ends by the caller).
int wrap (int i, int n) { return ((i % n) + n) % n; }

//! One pass along a line of n values (stride s from offset o in a): the running minimum (mean = false)
//! or mean (mean = true) over the window [i - p, i + p], wrapped when periodic, else cut at the ends.
void line_pass (const std::vector<double>& a, std::vector<double>& out, int n, int o, int s, int p, bool periodic, bool mean)
{
    for (int i = 0; i < n; ++i) {
        int lo = i - p, hi = i + p;
        if (!periodic) { lo = std::max(lo, 0); hi = std::min(hi, n - 1); }
        const long count = static_cast<long>(hi) - lo + 1;
        if (!mean) {
            double m = std::numeric_limits<double>::infinity();
            // A wrapped window as long as the line holds every node.
            const int last = periodic && count >= n ? lo + n - 1 : hi;
            for (int k = lo; k <= last; ++k) { m = std::min(m, a[static_cast<std::size_t>(o + s * (periodic ? wrap(k, n) : k))]); }
            out[static_cast<std::size_t>(o + s * i)] = m;
        } else {
            // A wrapped window longer than the line counts its nodes once per lap.
            double sum = 0;
            const long laps = periodic ? count / n : 0;
            if (laps > 0) {
                double lap = 0;
                for (int k = 0; k < n; ++k) { lap += a[static_cast<std::size_t>(o + s * k)]; }
                sum = static_cast<double>(laps) * lap;
            }
            const long rest = count - laps * n;
            for (long k = 0; k < rest; ++k) {
                const int kk = static_cast<int>(lo + k);
                sum += a[static_cast<std::size_t>(o + s * (periodic ? wrap(kk, n) : kk))];
            }
            out[static_cast<std::size_t>(o + s * i)] = sum / static_cast<double>(count);
        }
    }
}

//! The separable pass over the square: along x on every row, then along y on every column.
std::vector<double> square_pass (const std::vector<double>& a, int nx, int ny, int px, int py, bool perx, bool pery, bool mean)
{
    std::vector<double> t(a.size()), out(a.size());
    for (int j = 0; j < ny; ++j) { line_pass(a, t, nx, j * nx, 1, px, perx, mean); }
    for (int i = 0; i < nx; ++i) { line_pass(t, out, ny, i, nx, py, pery, mean); }
    return out;
}

}   // namespace

std::vector<double> lower_envelope (const std::vector<double>& h, int nx, int ny, int px, int py, bool perx, bool pery)
{
    if (nx < 1 || ny < 1 || static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) != h.size() || px < 0 || py < 0) {
        amrex::Abort("gust: lower_envelope: the array does not match its sizes");
    }
    if (px == 0 && py == 0) { return h; }
    const std::vector<double> m = square_pass(h, nx, ny, px, py, perx, pery, false);
    std::vector<double> e = square_pass(m, nx, ny, px, py, perx, pery, true);
    for (std::size_t n = 0; n < e.size(); ++n) { e[n] = std::min(e[n], h[n]); }
    return e;
}

std::vector<std::string> apply_hybrid_defaults (Inputs& in)
{
    std::vector<std::string> log;
    auto& t = in.terrain;
    if (t.representation != "hybrid") { return log; }
    auto given = [&t] (const char* k) { return std::find(t.given.begin(), t.given.end(), k) != t.given.end(); };
    auto metres = [] (amrex::Real v) {
        std::ostringstream o;
        o << std::setprecision(10) << static_cast<double>(v);
        return o.str();
    };
    // The lateral spacings of the directions with more than one cell (a direction of one cell does not count).
    std::vector<amrex::Real> steps;
    for (int d = 0; d < 2; ++d) {
        if (in.mesh.n_cell[d] > 1) {
            steps.push_back((in.mesh.prob_hi[d] - in.mesh.prob_lo[d]) / static_cast<amrex::Real>(in.mesh.n_cell[d]));
        }
    }
    if (!given("flat_height")) {
        const amrex::Vector<amrex::Real> zl = make_zlevels(in.mesh, make_geometry(in.mesh));
        t.flat_height = zl[1];
        log.push_back("terrain.flat_height = " + metres(t.flat_height) + " m (default: the first layer's thickness)");
    }
    if (!given("hybrid_half_width")) {
        amrex::Real cell = 0;
        for (const amrex::Real d : steps) { cell = amrex::max(cell, d); }
        t.hybrid_half_width = 2 * cell;
        log.push_back("terrain.hybrid_half_width = " + metres(t.hybrid_half_width) + " m (default: 2 cells)");
    } else if (t.hybrid_half_width > 0 && !steps.empty()) {
        bool none = true;
        for (const amrex::Real d : steps) { none = none && std::floor(static_cast<double>(t.hybrid_half_width) / static_cast<double>(d) + node_count_tol()) < 1; }
        if (none) {
            log.push_back("warning: terrain.hybrid_half_width = " + metres(t.hybrid_half_width)
                          + " m is below one cell: the square has no nodes and the mesh follows the terrain itself,"
                            " as with 0");
        }
    }
    if (t.type == "flat") {
        log.push_back("warning: terrain.representation = hybrid on flat terrain gives the immersed_forcing mesh");
    }
    return log;
}

amrex::IntVect hybrid_envelope_nodes (amrex::Real width, const amrex::Geometry& geom)
{
    amrex::IntVect p(0);
    for (int d = 0; d < 2; ++d) {
        p[d] = static_cast<int>(std::floor(static_cast<double>(width) / static_cast<double>(geom.CellSize(d)) + node_count_tol()));
    }
    return p;
}

void make_hybrid_envelope (amrex::MultiFab& z_nd, const amrex::Geometry& geom, amrex::Real width)
{
    const amrex::Box dom = geom.Domain();
    const int k0 = dom.smallEnd(2);
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    // The distinct nodes: a periodic direction's last node is its first.
    const int nx = dom.length(0) + (perx ? 0 : 1), ny = dom.length(1) + (pery ? 0 : 1);
    const int ilo = dom.smallEnd(0), jlo = dom.smallEnd(1);
    const auto idx = [=] (int i, int j) { return static_cast<std::size_t>(j - jlo) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i - ilo); };
    // Gathered onto every rank: each valid node from its box (a node shared by two boxes holds one value
    // there, so the largest is that value), the others at the lowest value.
    std::vector<amrex::Real> g(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), std::numeric_limits<amrex::Real>::lowest());
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        amrex::Box slab = mfi.validbox();
        slab.setSmall(2, k0);
        slab.setBig(2, k0);
        amrex::FArrayBox hf(slab, 1, amrex::The_Pinned_Arena());
        hf.copy<amrex::RunOn::Device>(z_nd[mfi], slab, 0, slab, 0, 1);
        amrex::Gpu::streamSynchronize();
        const auto a = hf.const_array();
        for (int j = slab.smallEnd(1); j <= slab.bigEnd(1); ++j) {
            for (int i = slab.smallEnd(0); i <= slab.bigEnd(0); ++i) {
                if (i - ilo >= nx || j - jlo >= ny) { continue; }   // a periodic direction's last node
                auto& v = g[idx(i, j)];
                v = std::max(v, a(i, j, k0));
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(g.data(), static_cast<int>(g.size()));
    const amrex::IntVect p = hybrid_envelope_nodes(width, geom);
    const std::vector<double> e = lower_envelope(std::vector<double>(g.begin(), g.end()), nx, ny, p[0], p[1], perx, pery);
    // Back into every column of every box, ghost columns included.
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        amrex::Box slab = mfi.fabbox();
        slab.setSmall(2, k0);
        slab.setBig(2, k0);
        amrex::FArrayBox hf(slab, 1, amrex::The_Pinned_Arena());
        const auto a = hf.array();
        for (int j = slab.smallEnd(1); j <= slab.bigEnd(1); ++j) {
            for (int i = slab.smallEnd(0); i <= slab.bigEnd(0); ++i) {
                const int ii = perx ? ilo + wrap(i - ilo, nx) : std::max(ilo, std::min(ilo + nx - 1, i));
                const int jj = pery ? jlo + wrap(j - jlo, ny) : std::max(jlo, std::min(jlo + ny - 1, j));
                a(i, j, k0) = static_cast<amrex::Real>(e[idx(ii, jj)]);
            }
        }
        z_nd[mfi].copy<amrex::RunOn::Device>(hf, slab, 0, slab, 0, 1);
        amrex::Gpu::streamSynchronize();
    }
}

std::string check_hybrid_ground_extent (const amrex::MultiFab& z_nd, const amrex::Geometry& geom)
{
    const int k0 = geom.Domain().smallEnd(2), ktop = geom.Domain().bigEnd(2), gc = z_nd.nComp() - 1;
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> data(op);
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        amrex::Box slab = mfi.validbox();
        slab.setSmall(2, k0);
        slab.setBig(2, k0);
        const auto z = z_nd.const_array(mfi);
        op.eval(slab, data, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real>
        {
            // How far the ground lies below the lower face (never, by construction), and above the top layer's bottom.
            const amrex::Real g = z(i, j, k0, gc);
            return {z(i, j, k0, 0) - g, g - z(i, j, ktop, 0), amrex::Math::abs(g)};
        });
    }
    auto r = data.value(op);
    amrex::Real below = amrex::get<0>(r), above = amrex::get<1>(r), mag = amrex::get<2>(r);
    amrex::ParallelDescriptor::ReduceRealMax(below);
    amrex::ParallelDescriptor::ReduceRealMax(above);
    amrex::ParallelDescriptor::ReduceRealMax(mag);
    // Round-off: a periodic face node is the envelope of its image node but the terrain at its own x.
    const amrex::Real tol = 64 * std::numeric_limits<amrex::Real>::epsilon() * amrex::max(amrex::Real(1), mag);
    std::string err;
    if (below > tol) {
        err = "the immersed ground lies up to " + std::to_string(static_cast<double>(below))
              + " m below the mesh's lower face with terrain.representation = hybrid (the envelope lies above the terrain)";
    } else if (above > 0) {
        err = "the immersed ground lies up to " + std::to_string(static_cast<double>(above))
              + " m above the bottom of the mesh's top layer with terrain.representation = hybrid; raise the top of the domain";
    }
    return err;
}

HybridGap hybrid_gap (const amrex::MultiFab& z_nd, const amrex::Vector<amrex::Real>& zlevels)
{
    const int gc = z_nd.nComp() - 1;
    const int k0 = z_nd.boxArray().minimalBox().smallEnd(2);
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real, amrex::Real> data(op);
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        amrex::Box slab = mfi.validbox();
        slab.setSmall(2, k0);
        slab.setBig(2, k0);
        const auto z = z_nd.const_array(mfi);
        op.eval(slab, data, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept -> amrex::GpuTuple<amrex::Real, amrex::Real>
        {
            const amrex::Real d = z(i, j, k0, gc) - z(i, j, k0, 0);
            return {d, d};
        });
    }
    auto r = data.value(op);
    amrex::Real lo = amrex::get<0>(r), hi = amrex::get<1>(r);
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    HybridGap gap;
    gap.min = static_cast<double>(lo);
    gap.max = static_cast<double>(hi);
    if (zlevels.size() >= 2) {
        const double d0 = static_cast<double>(zlevels[1] - zlevels[0]);
        std::size_t k = 1;
        while (k + 1 < zlevels.size() && std::abs(static_cast<double>(zlevels[k + 1] - zlevels[k]) - d0) <= std::max(1e-6, node_count_tol()) * d0) { ++k; }
        gap.uniform_top = static_cast<double>(zlevels[k] - zlevels[0]);
    }
    return gap;
}

amrex::Long hybrid_immersed_columns (const amrex::MultiFab& beta, const amrex::BoxArray& slab_ba,
                                     const amrex::DistributionMapping& dm, int ng, amrex::MultiFab& mask)
{
    if (ng > beta.nGrow()) { amrex::Abort("gust: hybrid_immersed_columns: the mask has more ghost columns than beta"); }
    if (dm != beta.DistributionMap() || slab_ba.size() != beta.boxArray().size()) {
        amrex::Abort("gust: hybrid_immersed_columns: the slab layout must be beta's lowest layer");
    }
    mask.define(slab_ba, dm, 1, amrex::IntVect(ng, ng, 0));
    mask.setVal(0);
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
        const auto m = mask.array(mfi);
        const auto b = beta.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            m(i, j, k) = b(i, j, k) > amrex::Real(0) ? amrex::Real(1) : amrex::Real(0);
        });
        op.eval(vb, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Long>
        {
            return {b(i, j, k) > amrex::Real(0) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long n = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

void hybrid_solid_pressure (const amrex::MultiFab& beta, const amrex::MultiFab& p, amrex::MultiFab& ps)
{
    if (ps.empty() || ps.boxArray() != p.boxArray() || ps.DistributionMap() != p.DistributionMap()
        || ps.nGrowVect() != p.nGrowVect()) {
        ps.define(p.boxArray(), p.DistributionMap(), 1, p.nGrowVect());
    }
    if (beta.nGrow() < 1 || p.nGrowVect()[0] < 1 || p.nGrowVect()[1] < 1 || p.nGrowVect()[2] < 1) {
        amrex::Abort("gust: hybrid mesh: the pressure fill needs one ghost cell of beta and p");
    }
    const amrex::Box all = p.boxArray().minimalBox();
    amrex::MultiFab::Copy(ps, p, 0, 0, 1, p.nGrowVect());
    for (amrex::MFIter mfi(ps, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
        if (klo != all.smallEnd(2) || khi != all.bigEnd(2)) {
            amrex::Abort("gust: hybrid mesh: the pressure fill needs boxes that span the column");
        }
        amrex::Box col = amrex::grow(vb, amrex::IntVect(1, 1, 0));
        col.setBig(2, klo);
        const auto b = beta.const_array(mfi);
        const auto a = ps.array(mfi);
        amrex::ParallelFor(col, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            amrex::Real above = a(i, j, khi);
            for (int k = khi; k >= klo; --k) {
                if (b(i, j, k) < amrex::Real(1)) {
                    above = a(i, j, k);
                } else {
                    a(i, j, k) = above;
                }
            }
            if (b(i, j, klo) >= amrex::Real(1)) { a(i, j, klo - 1) = above; }
        });
    }
}

}   // namespace gust
