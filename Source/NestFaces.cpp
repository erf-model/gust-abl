#include "NestFaces.H"

#include "BoundaryConditions_K.H"
#include "Reductions.H"
#include "Sha256.H"
#include "Stencil_K.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

namespace gust {

using namespace bc_detail;

namespace {

constexpr const char* kNestMagic = "gust_nest_faces";
constexpr int kNestVersion = 1;
const char* const kLateral[4] = {"xlo", "xhi", "ylo", "yhi"};

std::string nest_checkpoint_file (const std::string& dir) { return dir + "/NestFaces"; }

std::array<bool, 4> nest_faces_of (const BoundaryInputs& bin)
{
    std::array<bool, 4> nf{};
    for (int d = 0; d < 2; ++d) {
        nf[static_cast<std::size_t>(2 * d)] = bin.face[amrex::Orientation(d, amrex::Orientation::low)] == FaceType::nest;
        nf[static_cast<std::size_t>(2 * d + 1)] = bin.face[amrex::Orientation(d, amrex::Orientation::high)] == FaceType::nest;
    }
    return nf;
}

void bcast (std::string& s)
{
    long n = static_cast<long>(s.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

// Into the cells along face (d, lo): comp 0 the outward flux W s V of the
// cell's face point, comp 1 its weight W, comp 2 the flux entering,
// W max(0, -s V). A corner cell adds both its faces.
void accumulate_flux (const amrex::MultiFab& vn, const amrex::MultiFab& area, int d, bool lo, const amrex::Geometry& geom,
                      amrex::MultiFab& acc)
{
    const amrex::Box& dom = geom.Domain();
    const amrex::Real dxinv = geom.InvCellSize(d);
    const amrex::Real sgn = lo ? amrex::Real(-1) : amrex::Real(1);
    const int off = lo ? 0 : 1;
    const int ox = (d == 0) ? off : 0, oy = (d == 1) ? off : 0;
    const int edge = lo ? dom.smallEnd(d) : dom.bigEnd(d);
    for (amrex::MFIter mfi(acc, false); mfi.isValid(); ++mfi) {
        amrex::Box cb = mfi.validbox();
        if (cb.smallEnd(d) > edge || cb.bigEnd(d) < edge) { continue; }
        cb.setSmall(d, edge);
        cb.setBig(d, edge);
        const auto s = acc.array(mfi);
        const auto v = vn.const_array(mfi);
        const auto a = area.const_array(mfi);
        amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real W = a(i + ox, j + oy, k) * dxinv;
            const amrex::Real V = sgn * v(i + ox, j + oy, k);
            s(i, j, k, 0) += W * V;
            s(i, j, k, 1) += W;
            s(i, j, k, 2) += W * amrex::max(amrex::Real(0), -V);
        });
    }
}

// V <- V - s delta on the plane of face (d, lo) and on the ghost planes
// outside it (the parent's normal velocity there drops with the face's, so
// the gradient across the face is the parent's), over every fab that holds
// them (the ghost points along the face too: every copy of a point keeps one
// value).
void shift_face (amrex::MultiFab& vn, int d, bool lo, const amrex::Geometry& geom, amrex::Real delta)
{
    const amrex::Box& dom = geom.Domain();
    const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
    const amrex::Real change = lo ? delta : -delta;
    for (amrex::MFIter mfi(vn, false); mfi.isValid(); ++mfi) {
        amrex::Box fb = mfi.fabbox();
        if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }
        if (lo) { fb.setBig(d, face); } else { fb.setSmall(d, face); }
        vn[mfi].plus<amrex::RunOn::Device>(change, fb, 0, 1);
    }
}

// The zone of face (d, lo) for a field of type T: the index range along d
// (empty when lo > hi). The normal velocity's own face point is left out.
void zone_range (int T, int d, bool lo, const amrex::Box& dom, int N, int& z0, int& z1)
{
    const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
    if (stencil::nodal(T, d)) {
        z0 = lo ? face + 1 : face - (N - 1);
        z1 = lo ? face + N - 1 : face - 1;
    } else {
        z0 = lo ? face : face - N;
        z1 = lo ? face + N - 1 : face - 1;
    }
}

// rate += w(s) factor U_in / (N dx_d) in the zone of face (d, lo), for the
// points of type T (NestFaces.H). vn: the parent's normal velocity of
// direction d on the patch, ghost points included.
void add_rate (amrex::MultiFab& rate, int T, const amrex::MultiFab& vn, int d, bool lo, const amrex::Geometry& geom, int N,
               amrex::Real factor)
{
    const amrex::Box& dom = geom.Domain();
    const int e = 1 - d;
    const int oe = stencil::nodal(T, e) ? 1 : 0, oz = stencil::nodal(T, 2) ? 1 : 0;
    const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
    int z0 = 0, z1 = -1;
    zone_range(T, d, lo, dom, N, z0, z1);
    if (z0 > z1) { return; }
    // The face plane of vn under every box of rate (a box of the zone need
    // not touch the face), with the row before it along the face and below:
    // the box's cells, as vn's type, at the face.
    amrex::BoxList bl(vn.ixType());
    for (int b = 0; b < static_cast<int>(rate.boxArray().size()); ++b) {
        amrex::Box pb = amrex::convert(amrex::enclosedCells(rate.boxArray()[b]), vn.ixType());
        pb.setSmall(d, face);
        pb.setBig(d, face);
        bl.push_back(pb);
    }
    amrex::IntVect png(1);
    png[d] = 0;
    amrex::MultiFab plane(amrex::BoxArray(std::move(bl)), rate.DistributionMap(), 1, png);
    plane.setVal(0.0);
    // vn's ghost points too (outside the domain along the face, below the
    // ground): every copy of a point holds one value, so the order of the
    // overlapping copies does not matter.
    plane.ParallelCopy(vn, 0, 0, 1, vn.nGrowVect(), png, geom.periodicity());

    const amrex::Real inv_tau = factor / (static_cast<amrex::Real>(N) * geom.CellSize(d));
    const amrex::Real invN = amrex::Real(1) / static_cast<amrex::Real>(N);
    const amrex::Real half = stencil::nodal(T, d) ? amrex::Real(0) : amrex::Real(0.5);
    const amrex::Real inward = lo ? amrex::Real(1) : amrex::Real(-1);
    // Strides in place of a direction branch in the kernel: n = i or j, and
    // (fi, fj) is the point (i, j) moved to the face along d.
    const int di = (d == 0) ? 1 : 0, dj = 1 - di;
    const int ie = dj * oe, je = di * oe;
    const amrex::Real quarter = amrex::Real(0.25);
    for (amrex::MFIter mfi(rate, false); mfi.isValid(); ++mfi) {
        amrex::Box zb = mfi.validbox();
        zb.setSmall(d, std::max(zb.smallEnd(d), z0));
        zb.setBig(d, std::min(zb.bigEnd(d), z1));
        if (!zb.ok()) { continue; }
        const auto r = rate.array(mfi);
        const auto v = plane.const_array(mfi);
        amrex::ParallelFor(zb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int n = di * i + dj * j;
            const amrex::Real s = std::abs((static_cast<amrex::Real>(n) + half) - static_cast<amrex::Real>(face));
            const int fi = i + di * (face - i), fj = j + dj * (face - j);
            const amrex::Real U = quarter * ((v(fi, fj, k) + v(fi - ie, fj - je, k))
                                             + (v(fi, fj, k - oz) + v(fi - ie, fj - je, k - oz)));
            r(i, j, k) += (amrex::Real(1) - s * invN) * inv_tau * amrex::max(amrex::Real(0), inward * U);
        });
    }
}

// R += c (q - q_parent) on the boxes that reach a zone.
void add_source (const amrex::MultiFab& rate, const amrex::MultiFab& q, const amrex::MultiFab& par, amrex::MultiFab& R,
                 const std::vector<char>& in_zone)
{
    for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
        if (in_zone[static_cast<std::size_t>(mfi.index())] == 0) { continue; }
        const auto r = R.array(mfi);
        const auto c = rate.const_array(mfi);
        const auto a = q.const_array(mfi);
        const auto p = par.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) += c(i, j, k) * (a(i, j, k) - p(i, j, k));
        });
    }
}

void add_diagonal (const amrex::MultiFab& rate, amrex::MultiFab& A)
{
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const auto c = rate.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k, coef_P) += c(i, j, k);
        });
    }
}

// R = q - V (V null: R = q) on the face points of face (d, lo).
void face_rows (const amrex::MultiFab& q, const amrex::MultiFab* V, amrex::MultiFab& R, int d, bool lo, const amrex::Box& dom)
{
    const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
    for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
        amrex::Box fb = mfi.validbox();
        if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }
        fb.setSmall(d, face);
        fb.setBig(d, face);
        const auto r = R.array(mfi);
        const auto a = q.const_array(mfi);
        if (V != nullptr) {
            const auto v = V->const_array(mfi);
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { r(i, j, k) = a(i, j, k) - v(i, j, k); });
        } else {
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { r(i, j, k) = a(i, j, k); });
        }
    }
}

int rate_index (Loc loc)
{
    return (loc == Loc::xface) ? 0 : (loc == Loc::yface) ? 1 : (loc == Loc::zface) ? 2 : 3;
}

} // namespace

std::vector<std::string> check_nest_run (const BoundaryInputs& bin, const Inputs& in, bool first_guess_nest, bool sequencing)
{
    std::vector<std::string> err;
    const auto nf = nest_faces_of(bin);
    const bool any = nf[0] || nf[1] || nf[2] || nf[3];
    std::string first;   // the first nest face's key
    for (int n = 3; n >= 0; --n) { if (nf[static_cast<std::size_t>(n)]) { first = std::string("bc.") + kLateral[n]; } }
    if (first_guess_nest && bin.nest_parent.empty()) {
        err.emplace_back("init.first_guess = nest needs bc.nest_parent (the parent run's checkpoint)");
    }
    if (any && bin.nest_parent.empty()) {
        err.emplace_back(first + " = nest needs bc.nest_parent (the parent run's checkpoint)");
    }
    if (!first_guess_nest && !any && !bin.nest_parent.empty()) {
        err.emplace_back("bc.nest_parent has no effect without init.first_guess = nest or a nest face (bc.xlo, bc.xhi, "
                         "bc.ylo or bc.yhi = nest)");
    }
    if (any && !first_guess_nest && in.io.restart.empty()) {
        err.emplace_back(first + " = nest needs init.first_guess = nest (the parent's field as the first guess) or io.restart");
    }
    if (first_guess_nest && sequencing) {
        err.emplace_back("init.first_guess = nest is not supported with steady.sequencing = 1 (the patch starts on "
                         "level 0 from its parent)");
    }
    if (any && sequencing) {
        err.emplace_back(first + " = nest is not supported with steady.sequencing = 1 (the parent is read on level 0 only)");
    }
    // The ground blend.
    const amrex::Real W = in.terrain.nest_blend_width;
    const bool given = std::find(in.terrain.given.begin(), in.terrain.given.end(), "nest_blend_width") != in.terrain.given.end();
    if (!std::isfinite(W) || W < 0) {
        err.emplace_back("terrain.nest_blend_width = " + std::to_string(W) + ": must be finite and >= 0 (m; 0: off)");
    } else if (given && !any) {
        err.emplace_back("terrain.nest_blend_width has no effect without a nest face (bc.xlo, bc.xhi, bc.ylo or bc.yhi = nest)");
    } else if (W > 0) {
        for (int d = 0; d < 2; ++d) {
            if (!nf[static_cast<std::size_t>(2 * d)] && !nf[static_cast<std::size_t>(2 * d + 1)]) { continue; }
            const amrex::Real L = in.mesh.prob_hi[d] - in.mesh.prob_lo[d];
            if (!(W < amrex::Real(0.5) * L)) {
                std::ostringstream os;
                os << "terrain.nest_blend_width = " << W << ": must be below half the patch's width in " << (d == 0 ? "x" : "y")
                   << " (" << amrex::Real(0.5) * L << " m)";
                err.push_back(os.str());
            }
        }
    }
    return err;
}

void blend_nest_ground (const NestParent& parent, const BoundaryInputs& bin, amrex::Real width, const amrex::Geometry& geom,
                        amrex::MultiFab& z_phys_nd)
{
    if (!(width > 0)) { return; }
    const auto nf = nest_faces_of(bin);
    // 0 for a nest face, 1 for another: added to d / W, it takes the other
    // faces out of the minimum (t <= 1) with no branch in the kernel.
    const amrex::GpuArray<amrex::Real, 4> other{amrex::Real(nf[0] ? 0 : 1), amrex::Real(nf[1] ? 0 : 1),
                                                amrex::Real(nf[2] ? 0 : 1), amrex::Real(nf[3] ? 0 : 1)};
    amrex::MultiFab hp(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, amrex::IntVect(1, 1, 0));
    parent.interpolate_ground(hp);
    const amrex::Box& dom = geom.Domain();
    const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0) + 1, jlo = dom.smallEnd(1), jhi = dom.bigEnd(1) + 1;
    const int k0 = dom.smallEnd(2);
    // The index range a ghost node is clamped to: the domain's nodes across a
    // non-periodic direction, every index across a periodic one (no clamp).
    const int big = std::numeric_limits<int>::max() / 2;
    const int ci0 = geom.isPeriodic(0) ? -big : ilo, ci1 = geom.isPeriodic(0) ? big : ihi;
    const int cj0 = geom.isPeriodic(1) ? -big : jlo, cj1 = geom.isPeriodic(1) ? big : jhi;
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();
    const amrex::Real invW = amrex::Real(1) / width;
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.fabbox();
        if (b.smallEnd(2) > k0 || b.bigEnd(2) < k0) { continue; }
        b.setSmall(2, k0);
        b.setBig(2, k0);
        b &= amrex::grow(hp[mfi].box(), 2, 1);   // hp has no ghost node in z
        const auto z = z_phys_nd.array(mfi);
        const auto h = hp.const_array(mfi);
        reduce_op.eval(b, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            // A ghost node outside a non-periodic face: the nearest domain node's value.
            const int ic = amrex::max(ci0, amrex::min(ci1, i));
            const int jc = amrex::max(cj0, amrex::min(cj1, j));
            const amrex::Real x = plo[0] + static_cast<amrex::Real>(ic - ilo) * dx[0];
            const amrex::Real y = plo[1] + static_cast<amrex::Real>(jc - jlo) * dx[1];
            const amrex::Real dist[4] = {x - plo[0], phi[0] - x, y - plo[1], phi[1] - y};
            amrex::Real t = 1;
            for (int n = 0; n < 4; ++n) { t = amrex::min(t, amrex::max(amrex::Real(0), dist[n]) * invW + other[n]); }
            const amrex::Real s = t * t * (amrex::Real(3) - amrex::Real(2) * t);
            // (1 - s) h_parent + s h: the parent's exactly at s = 0, the patch's own at s = 1.
            const amrex::Real old = z(i, j, k);
            z(i, j, k) = (amrex::Real(1) - s) * h(ic, jc, k) + s * old;
            return {std::abs(z(i, j, k) - old)};
        });
    }
    amrex::Real moved = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(moved);
    char line[200];
    std::snprintf(line, sizeof(line), "gust: nest faces: terrain.nest_blend_width = %g m: the ground is the parent's on the nest "
                  "faces; it moved by up to %.3e m\n", static_cast<double>(width), static_cast<double>(moved));
    amrex::Print() << line;
}

NestFaces::NestFaces (const NestParent& parent, const BoundaryInputs& bin, const amrex::Geometry& geom, const Mesh& mesh)
    : m_nest(nest_faces_of(bin)), m_geom(geom), m_relax_width(bin.nest_relax_width), m_relax_factor(bin.nest_relax_factor)
{
    if (!(m_nest[0] || m_nest[1] || m_nest[2] || m_nest[3])) {
        amrex::Abort("gust: NestFaces: no lateral face is a nest face (bc.xlo, bc.xhi, bc.ylo or bc.yhi = nest)");
    }
    const amrex::BoxArray& ba = mesh.detJ.boxArray();
    const amrex::DistributionMapping& dm = mesh.detJ.DistributionMap();
    m_par.define(ba, dm);
    parent.interpolate(m_par);
    m_par.p.clear();
    m_par.z_phys_nd.clear();

    // The net flux through the nest faces, in a fixed order of the cells.
    amrex::MultiFab* vn[2] = {&m_par.u, &m_par.v};
    {
        amrex::MultiFab acc(ba, dm, 3, 0);
        acc.setVal(0.0);
        for (int d = 0; d < 2; ++d) {
            for (const bool lo : {true, false}) {
                if (is_nest(d, lo)) { accumulate_flux(*vn[d], mesh.flux_area(d), d, lo, geom, acc); }
            }
        }
        // In double whatever amrex::Real is: the net is a small difference of large fluxes.
        amrex::MultiFab one(ba, dm, 1, 0);
        one.setVal(1.0);
        m_flux.net = deterministic_dot(acc, 0, one, 0, nullptr);
        m_flux.weight = deterministic_dot(acc, 1, one, 0, nullptr);
        m_flux.inflow = deterministic_dot(acc, 2, one, 0, nullptr);
    }
    bool outflow = false;
    for (const FaceType t : bin.face) { outflow = outflow || t == FaceType::outflow; }
    m_flux.applied = !outflow;
    if (m_flux.applied) {
        if (!(m_flux.weight > 0) || !std::isfinite(m_flux.net)) {
            amrex::Abort("gust: bc.nest_parent = " + parent.dir() + ": the parent's velocity on the nest faces is not finite");
        }
        m_flux.delta = static_cast<amrex::Real>(m_flux.net / m_flux.weight);
        for (int d = 0; d < 2; ++d) {
            for (const bool lo : {true, false}) {
                if (is_nest(d, lo)) { shift_face(*vn[d], d, lo, geom, m_flux.delta); }
            }
        }
    }

    if (m_relax_width > 0) {
        const amrex::IntVect types[4] = {amrex::IntVect(1, 0, 0), amrex::IntVect(0, 1, 0), amrex::IntVect(0, 0, 1), amrex::IntVect(0)};
        const int T[4] = {stencil::XF, stencil::YF, stencil::ZF, stencil::CC};
        for (int n = 0; n < 4; ++n) {
            m_rate[static_cast<std::size_t>(n)].define(amrex::convert(ba, types[n]), dm, 1, 0);
            m_rate[static_cast<std::size_t>(n)].setVal(0.0);
            for (int d = 0; d < 2; ++d) {
                for (const bool lo : {true, false}) {
                    if (is_nest(d, lo)) {
                        add_rate(m_rate[static_cast<std::size_t>(n)], T[n], *vn[d], d, lo, geom, m_relax_width, m_relax_factor);
                    }
                }
            }
        }
    }
    amrex::Gpu::streamSynchronize();
}

void NestFaces::fill (State& q, int d, bool homogeneous) const
{
    if (d > 1 || !(is_nest(d, true) || is_nest(d, false))) { return; }
    if (q.p.boxArray() != m_par.theta.boxArray() || q.p.DistributionMap() != m_par.theta.DistributionMap()) {
        amrex::Abort("gust: NestFaces::fill: the state is not on the layout the parent was read for");
    }
    const amrex::Box& dom = m_geom.Domain();
    amrex::MultiFab* dst[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    const amrex::MultiFab* src[5] = {&m_par.u, &m_par.v, &m_par.w, &m_par.theta, &m_par.tke};
    for (const bool lo : {true, false}) {
        if (!is_nest(d, lo)) { continue; }
        for (int n = 0; n < 5; ++n) {
            amrex::MultiFab& mf = *dst[n];
            const int ng = mf.nGrow();
            const bool nodal_d = mf.ixType().nodeCentered(d);
            const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
            // Layers 1 to ng, and the face point (layer 0) of the normal velocity.
            const int inner = nodal_d ? 0 : 1;
            for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                if (!touches(vb, dom, d, lo)) { continue; }
                amrex::Box gb = ghost_slab(vb, dom, m_geom, d, lo, inner, ng);
                gb.setSmall(d, lo ? edge - ng : edge + inner);
                gb.setBig(d, lo ? edge - inner : edge + ng);
                gb &= mfi.fabbox();
                if (gb.isEmpty()) { continue; }
                if (homogeneous) {
                    mf[mfi].setVal<amrex::RunOn::Device>(amrex::Real(0), gb, 0, 1);
                } else {
                    mf[mfi].copy<amrex::RunOn::Device>((*src[n])[mfi], gb, 0, gb, 0, 1);
                }
            }
        }
    }
}

void NestFaces::apply_rows (const State& q, ResidualFields& R, bool homogeneous) const
{
    const amrex::MultiFab* qf[2] = {&q.u, &q.v};
    const amrex::MultiFab* pf[2] = {&m_par.u, &m_par.v};
    amrex::MultiFab* rf[2] = {&R.u, &R.v};
    for (int d = 0; d < 2; ++d) {
        for (const bool lo : {true, false}) {
            if (is_nest(d, lo)) { face_rows(*qf[d], homogeneous ? nullptr : pf[d], *rf[d], d, lo, m_geom.Domain()); }
        }
    }
}

void NestFaces::add_relaxation (const State& q, ResidualFields& R) const
{
    if (!relax_on()) { return; }
    const amrex::Box& dom = m_geom.Domain();
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    const amrex::MultiFab* pf[5] = {&m_par.u, &m_par.v, &m_par.w, &m_par.theta, &m_par.tke};
    amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const int T[5] = {stencil::XF, stencil::YF, stencil::ZF, stencil::CC, stencil::CC};
    for (int n = 0; n < 5; ++n) {
        const amrex::MultiFab& rate = m_rate[static_cast<std::size_t>(std::min(n, 3))];
        // The boxes that reach a zone (the rate is 0 elsewhere).
        std::vector<char> in_zone(static_cast<std::size_t>(rate.boxArray().size()), 0);
        for (amrex::MFIter mfi(rate, false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            for (int d = 0; d < 2; ++d) {
                for (const bool lo : {true, false}) {
                    if (!is_nest(d, lo)) { continue; }
                    int z0 = 0, z1 = -1;
                    zone_range(T[n], d, lo, dom, m_relax_width, z0, z1);
                    if (z0 <= z1 && vb.smallEnd(d) <= z1 && vb.bigEnd(d) >= z0) { in_zone[static_cast<std::size_t>(mfi.index())] = 1; }
                }
            }
        }
        add_source(rate, *qf[n], *pf[n], *rf[n], in_zone);
    }
}

void NestFaces::add_relaxation_to_implicit (Equation /*eq*/, Loc loc, amrex::MultiFab& A) const
{
    if (!relax_on()) { return; }
    add_diagonal(m_rate[static_cast<std::size_t>(rate_index(loc))], A);
}

const amrex::MultiFab& NestFaces::rate (Loc loc) const { return m_rate[static_cast<std::size_t>(rate_index(loc))]; }

std::string NestFaces::summary () const
{
    std::ostringstream os;
    os << "gust: nest faces:";
    for (int n = 0; n < 4; ++n) { if (m_nest[static_cast<std::size_t>(n)]) { os << " " << kLateral[n]; } }
    os << " take the parent's u, v, w, theta and k\n";
    char line[384];
    char what[96];
    if (m_flux.inflow > 0) {
        std::snprintf(what, sizeof(what), "net outward flux %.3e of the flux entering", m_flux.net / m_flux.inflow);
    } else {
        // Nothing enters: the net as a mean outward velocity (there is no flux to compare it with).
        std::snprintf(what, sizeof(what), "no flux enters; mean outward velocity %.3e m/s", m_flux.net / m_flux.weight);
    }
    if (m_flux.applied) {
        std::snprintf(line, sizeof(line), "gust: nest faces: %s; removed before the pressure solve: every nest face's outward "
                      "velocity changed by %.3e m/s\n", what, -static_cast<double>(m_flux.delta));
    } else {
        std::snprintf(line, sizeof(line), "gust: nest faces: %s; not removed (an outflow face carries it)\n", what);
    }
    os << line;
    const bool large = std::abs(m_flux.net) > kNestNetWarning * m_flux.inflow;
    if (m_flux.applied && m_flux.inflow > 0 && large) {
        os << "gust: WARNING: nest faces: the net flux removed is above " << kNestNetWarning << " of the flux entering: the "
           << "parent's flow on the nest faces is changed by that much. The usual cause is a patch ground that differs from the "
           << "parent's on the nest faces (terrain.nest_blend_width blends the two)\n";
    }
    if (relax_on()) {
        os << "gust: nest faces: relaxation zone of " << m_relax_width << " cells (bc.nest_relax_factor = " << m_relax_factor
           << ") where the parent's flow enters\n";
    }
    return os.str();
}

std::string NestFaces::hash_inputs (const NestParent& parent, const BoundaryInputs& bin, amrex::Real blend_width)
{
    // Every lateral face's type: a face that is not nest decides too whether the net flux is removed.
    std::string s = "faces";
    for (int n = 0; n < 4; ++n) {
        const FaceType t = bin.face[amrex::Orientation(n / 2, n % 2 == 0 ? amrex::Orientation::low : amrex::Orientation::high)];
        s += std::string(" ") + kLateral[n] + "=" + face_type_name(t);
    }
    s += "\nparent step " + std::to_string(parent.parent_step());
    for (int d = 0; d < 3; ++d) { s += " " + std::to_string(parent.parent_cells()[d]); }
    s += "\nparent content " + parent.content_digest();
    s += "\noffset " + std::to_string(parent.offset()[0]) + " " + std::to_string(parent.offset()[1]);
    s += "\nblend " + hexfloat_canonical(static_cast<double>(blend_width));
    s += "\nrelax " + std::to_string(bin.nest_relax_width) + " " + hexfloat_canonical(static_cast<double>(bin.nest_relax_factor)) + "\n";
    return sha256_hex(s);
}

void NestFaces::write_checkpoint (const std::string& dir, const std::string& hash)
{
    if (!hash.empty() && amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = nest_checkpoint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << kNestMagic << " " << kNestVersion << "\nsha256 " << hash << "\n";
    }
    amrex::ParallelDescriptor::Barrier();
}

void NestFaces::check_restart (const std::string& dir, const std::string& hash)
{
    const bool active = !hash.empty();
    // The checkpoint's hash ("" without a NestFaces file; "?" when unreadable).
    std::string rec;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(nest_checkpoint_file(dir));
        if (is) {
            std::string magic, key;
            int version = 0;
            if (!(is >> magic >> version >> key >> rec && magic == kNestMagic && version == kNestVersion && key == "sha256")) {
                rec = "?";
            }
        }
    }
    bcast(rec);
    const std::string where = "gust: io.restart = " + dir + ": ";
    if (rec == "?") { amrex::Abort(where + "its NestFaces file cannot be read"); }
    if (active && rec.empty()) {
        amrex::Abort(where + "the checkpoint was written without nest faces; a restart cannot add them (bc.xlo, bc.xhi, "
                             "bc.ylo, bc.yhi)");
    }
    if (!active && !rec.empty()) {
        amrex::Abort(where + "the checkpoint was written with nest faces; give the run's bc.* faces and bc.nest_parent");
    }
    if (active && rec != hash) {
        amrex::Abort(where + "the checkpoint was written with other nest settings (the types of the lateral faces bc.xlo, "
                             "bc.xhi, bc.ylo, bc.yhi; the parent in bc.nest_parent: its step, cells, fields or ground, or the "
                             "patch's place in it; bc.nest_relax_width, bc.nest_relax_factor or terrain.nest_blend_width); "
                             "restart with the inputs and the parent it was written with");
    }
}

} // namespace gust
