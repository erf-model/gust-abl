#include "NestFaceFlux.H"

#include "Domain.H"
#include "Mesh.H"
#include "NestParent.H"
#include "State.H"
#include "Stencil_K.H"

#include <AMReX_Print.H>

#include <cstdio>
#include <memory>
#include <string>

namespace gust {

namespace {

constexpr int kEq = 5;
constexpr int kLoc[kEq] = {stencil::XF, stencil::YF, stencil::ZF, stencil::CC, stencil::CC};

amrex::MultiFab& field (ResidualFields& R, int eq)
{
    amrex::MultiFab* f[kEq] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    return *f[eq];
}

const amrex::MultiFab& field (const ResidualFields& R, int eq)
{
    const amrex::MultiFab* f[kEq] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    return *f[eq];
}

int part_index (FluxPart p)
{
    return p == FluxPart::convective ? 0 : p == FluxPart::convective_upwind ? 1 : 2;
}

//! Whether a region has faces along d (its inner box is shrunk there: not a periodic direction the nest spans).
bool has_faces (const DefectRegion& g, int d)
{
    return g.inner.smallEnd(d) != g.footprint.smallEnd(d);
}

//! The ring of side s (0: below chi_s along d, 1: above) for a field at loc t: chi_s's points (the inner box at t)
//! one point past chi_s along d.
amrex::Box ring_box (const DefectRegion& g, int t, int d, int s)
{
    amrex::Box b = amrex::convert(g.inner, amrex::IntVect(t & 1, (t >> 1) & 1, (t >> 2) & 1));
    const int at = s == 0 ? b.smallEnd(d) - 1 : b.bigEnd(d) + 1;
    b.setRange(d, at, 1);
    return b;
}

//! acc += sgn dxinv G(face) on the ring points in ring: the parent's own flux through chi_s's face (the face is
//! the ring point's high face for sgn = +1, its low face for sgn = -1).
template <int T, int D>
void add_coarse (const amrex::MultiFab& G, const amrex::Box& ring, int sgn, amrex::Real dxinv, amrex::MultiFab& acc)
{
    constexpr int flo = stencil::flux_lo(T, D);
    const int off = (sgn > 0 ? 1 : 0) + flo;
    const amrex::Real f = amrex::Real(sgn) * dxinv;
    for (amrex::MFIter mfi(acc, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & ring;
        if (!bx.ok()) { continue; }
        const auto a = acc.array(mfi);
        const auto g = G.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) += f * g(i + (D == 0 ? off : 0), j + (D == 1 ? off : 0), k);
        });
    }
}

//! The same with the nest's flux through the parent face: acc (on the parent cells under the nest, the nest's
//! DistributionMapping) += sgn dxinv G_f with G_f the mean of the nest's flux densities over the parent face.
//! Along the normal D: the nest face on the parent face (G nodal in D), or the mean of the two nest faces either
//! side of it (G cell-centred in D: the control volume of a field nodal in D, whose faces sit at the cell centres;
//! r even). Across (t = 1 - D): the r nest faces in the parent face (G cell-centred in t), or the trapezoid rule
//! over the r + 1 nest nodes from one parent node to the next centred on it (G nodal in t: weights 1/r, 1/(2r) at
//! the two ends). z: the same index (ratio 1). G: the nest's fluxes with r ghost points filled.
template <int T, int D>
void add_fine (const amrex::MultiFab& G, int r, const amrex::IntVect& offset, const amrex::Box& ring, int sgn,
               amrex::Real dxinv, amrex::MultiFab& acc)
{
    constexpr int TG = T ^ (1 << D);
    constexpr int t = 1 - D;
    constexpr int flo = stencil::flux_lo(T, D);
    constexpr bool normal_nodal = stencil::nodal(TG, D);
    constexpr bool cross_nodal = stencil::nodal(TG, t);
    const int off = (sgn > 0 ? 1 : 0) + flo;
    const int oD = offset[D], offT = offset[t];
    const int h = r / 2;
    const amrex::Real f = amrex::Real(sgn) * dxinv;
    const amrex::Real rinv = amrex::Real(1) / amrex::Real(r);
    for (amrex::MFIter mfi(acc, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & ring;
        if (!bx.ok()) { continue; }
        const auto a = acc.array(mfi);
        const auto g = G.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int pD = (D == 0) ? i : j;
            const int pt = (D == 0) ? j : i;
            const int nD0 = normal_nodal ? r * (pD + off - oD) : r * (pD + off - oD) + h - 1;
            const int nn = normal_nodal ? 1 : 2;
            const amrex::Real wn = normal_nodal ? amrex::Real(1) : amrex::Real(0.5);
            const int t0 = cross_nodal ? r * (pt - offT) - h : r * (pt - offT);
            const int nt = cross_nodal ? r + 1 : r;
            amrex::Real s = 0;
            for (int an = 0; an < nn; ++an) {
                for (int at = 0; at < nt; ++at) {
                    const amrex::Real wt = (cross_nodal && (at == 0 || at == r)) ? amrex::Real(0.5) * rinv : rinv;
                    const int iD = nD0 + an, it = t0 + at;
                    s += wn * wt * ((D == 0) ? g(iD, it, k) : g(it, iD, k));
                }
            }
            a(i, j, k) += f * s;
        });
    }
}

template <int T>
void add_coarse_t (int d, const amrex::MultiFab& G, const amrex::Box& ring, int sgn, amrex::Real dxinv, amrex::MultiFab& acc)
{
    if (d == 0) { add_coarse<T, 0>(G, ring, sgn, dxinv, acc); } else { add_coarse<T, 1>(G, ring, sgn, dxinv, acc); }
}

template <int T>
void add_fine_t (int d, const amrex::MultiFab& G, int r, const amrex::IntVect& offset, const amrex::Box& ring, int sgn,
                 amrex::Real dxinv, amrex::MultiFab& acc)
{
    if (d == 0) { add_fine<T, 0>(G, r, offset, ring, sgn, dxinv, acc); } else { add_fine<T, 1>(G, r, offset, ring, sgn, dxinv, acc); }
}

//! The ring sums of one domain's residual: the parent's (fine == false) on its own layout, or a nest's on the
//! parent cells under it (its DistributionMapping).
class RingProbe : public FluxProbe
{
public:
    RingProbe (std::vector<const DefectRegion*> regions, const amrex::GpuArray<amrex::Real, 3>& dxinv,
               std::array<ResidualFields, 3>& acc)
        : m_regions(std::move(regions)), m_dxinv(dxinv), m_acc(acc) {}

    //! A nest's probe: its fluxes are summed onto the parent faces (ratio r, offset), with ghost points filled
    //! with the nest's periodicity.
    void set_fine (int r, const amrex::IntVect& offset, const amrex::Geometry* geom)
    {
        m_fine = true;
        m_r = r;
        m_offset = offset;
        m_geom = geom;
    }

    void take (int eq, FluxPart part, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G,
               const amrex::MultiFab* sigma) override
    {
        if (eq < 0 || eq >= kEq) { amrex::Abort("gust: NestFaceFlux: flux of an unknown equation " + std::to_string(eq)); }
        const int pi = part_index(part);
        if (part == FluxPart::convective_upwind) {
            upwind[eq] = true;
            this->sigma = sigma;
        }
        amrex::MultiFab& acc = field(m_acc[pi], eq);
        const int t = kLoc[eq];
        for (int d = 0; d < 2; ++d) {
            // A copy with the ghost points the nest's sums read, filled across its boxes.
            amrex::MultiFab Gf;
            const amrex::MultiFab* Gd = &G[d];
            bool any = false;
            for (const DefectRegion* g : m_regions) { any = any || has_faces(*g, d); }
            if (!any) { continue; }
            if (m_fine) {
                Gf.define(G[d].boxArray(), G[d].DistributionMap(), 1, m_r);
                amrex::MultiFab::Copy(Gf, G[d], 0, 0, 1, 0);
                Gf.FillBoundary(m_geom->periodicity());
                Gd = &Gf;
            }
            for (const DefectRegion* g : m_regions) {
                if (!has_faces(*g, d)) { continue; }
                for (int s = 0; s < 2; ++s) {
                    const amrex::Box ring = ring_box(*g, t, d, s);
                    const int sgn = s == 0 ? 1 : -1;
                    switch (t) {
                    case stencil::CC:
                        if (m_fine) { add_fine_t<stencil::CC>(d, *Gd, m_r, m_offset, ring, sgn, m_dxinv[d], acc); }
                        else { add_coarse_t<stencil::CC>(d, *Gd, ring, sgn, m_dxinv[d], acc); }
                        break;
                    case stencil::XF:
                        if (m_fine) { add_fine_t<stencil::XF>(d, *Gd, m_r, m_offset, ring, sgn, m_dxinv[d], acc); }
                        else { add_coarse_t<stencil::XF>(d, *Gd, ring, sgn, m_dxinv[d], acc); }
                        break;
                    case stencil::YF:
                        if (m_fine) { add_fine_t<stencil::YF>(d, *Gd, m_r, m_offset, ring, sgn, m_dxinv[d], acc); }
                        else { add_coarse_t<stencil::YF>(d, *Gd, ring, sgn, m_dxinv[d], acc); }
                        break;
                    default:
                        if (m_fine) { add_fine_t<stencil::ZF>(d, *Gd, m_r, m_offset, ring, sgn, m_dxinv[d], acc); }
                        else { add_coarse_t<stencil::ZF>(d, *Gd, ring, sgn, m_dxinv[d], acc); }
                        break;
                    }
                }
            }
        }
    }

    std::array<bool, kEq> upwind{};
    const amrex::MultiFab* sigma{nullptr};

private:
    std::vector<const DefectRegion*> m_regions;
    amrex::GpuArray<amrex::Real, 3> m_dxinv;
    std::array<ResidualFields, 3>& m_acc;
    bool m_fine{false};
    int m_r{2};
    amrex::IntVect m_offset{0};
    const amrex::Geometry* m_geom{nullptr};
};

//! C = ((A_c - A_d) - (F_c - F_d)) / J, with A the parent's ring sums and F the nests' (no blend).
template <int T>
void assemble (const amrex::MultiFab& Ac, const amrex::MultiFab& Ad, const amrex::MultiFab& Fc, const amrex::MultiFab& Fd,
               const amrex::MultiFab& detJ, amrex::MultiFab& C)
{
    for (amrex::MFIter mfi(C, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto c = C.array(mfi);
        const auto ac = Ac.const_array(mfi);
        const auto ad = Ad.const_array(mfi);
        const auto fc = Fc.const_array(mfi);
        const auto fd = Fd.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            c(i, j, k) = ((ac(i, j, k) - ad(i, j, k)) - (fc(i, j, k) - fd(i, j, k))) / stencil::interp<stencil::CC, T>(J, i, j, k);
        });
    }
}

//! The same with the upwind blend: + w ((A_u - A_c) - (F_u - F_c)) / J, w the blend weight at the point.
template <int T>
void assemble_blend (const amrex::MultiFab& Ac, const amrex::MultiFab& Au, const amrex::MultiFab& Ad,
                     const amrex::MultiFab& Fc, const amrex::MultiFab& Fu, const amrex::MultiFab& Fd,
                     const amrex::MultiFab& sigma, const amrex::MultiFab& detJ, amrex::MultiFab& C)
{
    for (amrex::MFIter mfi(C, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto c = C.array(mfi);
        const auto ac = Ac.const_array(mfi);
        const auto au = Au.const_array(mfi);
        const auto ad = Ad.const_array(mfi);
        const auto fc = Fc.const_array(mfi);
        const auto fu = Fu.const_array(mfi);
        const auto fd = Fd.const_array(mfi);
        const auto s = sigma.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real w = stencil::interp<stencil::CC, T>(s, i, j, k);
            c(i, j, k) = ((ac(i, j, k) - ad(i, j, k)) - (fc(i, j, k) - fd(i, j, k))
                          + w * ((au(i, j, k) - ac(i, j, k)) - (fu(i, j, k) - fc(i, j, k))))
                         / stencil::interp<stencil::CC, T>(J, i, j, k);
        });
    }
}

template <int T>
void assemble_any (const std::array<ResidualFields, 3>& A, const std::array<ResidualFields, 3>& F, int eq,
                   const amrex::MultiFab* sigma, const amrex::MultiFab& detJ, amrex::MultiFab& C)
{
    if (sigma != nullptr) {
        assemble_blend<T>(field(A[0], eq), field(A[1], eq), field(A[2], eq), field(F[0], eq), field(F[1], eq),
                          field(F[2], eq), *sigma, detJ, C);
    } else {
        assemble<T>(field(A[0], eq), field(A[2], eq), field(F[0], eq), field(F[2], eq), detJ, C);
    }
}

void zero (std::array<ResidualFields, 3>& a)
{
    for (auto& r : a) {
        for (int eq = 0; eq < kEq; ++eq) { field(r, eq).setVal(0.0); }
    }
}

} // namespace

NestFaceFlux::NestFaceFlux (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const Mesh& mesh,
                            const amrex::Geometry& geom)
    : m_mesh(&mesh), m_geom(&geom)
{
    m_C.define(ba, dm);
    for (auto& r : m_A) { r.define(ba, dm); }
    for (auto& r : m_F) { r.define(ba, dm); }
    for (int eq = 0; eq < kEq; ++eq) { field(m_C, eq).setVal(0.0); }
    m_C.p.setVal(0.0);
}

void NestFaceFlux::compute (const std::function<void(FluxProbe&)>& coarse, const std::vector<FaceFluxNest>& nests)
{
    const amrex::GpuArray<amrex::Real, 3> dxinv = m_geom->InvCellSizeArray();
    std::vector<const DefectRegion*> all;
    for (const FaceFluxNest& n : nests) {
        if (n.region == nullptr || n.geom == nullptr || !n.residual) {
            amrex::Abort("gust: NestFaceFlux::compute: a nest without its region, geometry or residual");
        }
        if (n.ratio < 2 || n.ratio % 2 != 0) {
            amrex::Abort("gust: " + n.region->name + ": nest.feedback_flux = nest needs an even refinement ratio, not "
                         + std::to_string(n.ratio));
        }
        all.push_back(n.region);
    }
    // The parent's ring sums.
    zero(m_A);
    RingProbe pc(all, dxinv, m_A);
    coarse(pc);
    // Each nest's, on the parent cells under it, then onto the parent's layout.
    zero(m_F);
    for (const FaceFluxNest& n : nests) {
        std::array<ResidualFields, 3> img;
        // The nest's layout comes with its fluxes; the images are made at the first one.
        struct Lazy : FluxProbe
        {
            const FaceFluxNest* n;
            const amrex::GpuArray<amrex::Real, 3>* dxinv;
            std::array<ResidualFields, 3>* img;
            std::unique_ptr<RingProbe> probe;
            void take (int eq, FluxPart part, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G,
                       const amrex::MultiFab* sigma) override
            {
                if (!probe) {
                    amrex::BoxArray ba = amrex::convert(G[0].boxArray(), amrex::IntVect(0));
                    if (!ba.coarsenable(amrex::IntVect(n->ratio, n->ratio, 1))) {
                        amrex::Abort("gust: " + n->region->name + ": the nest's boxes do not coarsen by "
                                     + std::to_string(n->ratio));
                    }
                    ba.coarsen(amrex::IntVect(n->ratio, n->ratio, 1));
                    ba.shift(amrex::IntVect(n->offset[0], n->offset[1], 0));
                    for (auto& r : *img) { r.define(ba, G[0].DistributionMap()); }
                    zero(*img);
                    probe = std::make_unique<RingProbe>(std::vector<const DefectRegion*>{n->region}, *dxinv, *img);
                    probe->set_fine(n->ratio, n->offset, n->geom);
                }
                probe->take(eq, part, G, sigma);
            }
        } lazy;
        lazy.n = &n;
        lazy.dxinv = &dxinv;
        lazy.img = &img;
        n.residual(lazy);
        if (!lazy.probe) { amrex::Abort("gust: " + n.region->name + ": its residual passed no fluxes to the probe"); }
        for (int eq = 0; eq < kEq; ++eq) {
            // A nest without the blend (or k, which has none): its upwind sum is its convective one.
            if (!lazy.probe->upwind[eq]) {
                amrex::MultiFab::Copy(field(img[1], eq), field(img[0], eq), 0, 0, 1, 0);
            }
            for (int p = 0; p < 3; ++p) {
                amrex::MultiFab& dst = field(m_F[p], eq);
                amrex::MultiFab tmp(dst.boxArray(), dst.DistributionMap(), 1, 0);
                tmp.setVal(0.0);
                tmp.ParallelCopy(field(img[p], eq), 0, 0, 1);
                amrex::MultiFab::Add(dst, tmp, 0, 0, 1, 0);
            }
        }
    }
    // The correction, one kernel per location (with or without the parent's blend).
    for (int eq = 0; eq < kEq; ++eq) {
        const amrex::MultiFab* sigma = pc.upwind[eq] ? pc.sigma : nullptr;
        amrex::MultiFab& C = field(m_C, eq);
        switch (kLoc[eq]) {
        case stencil::CC: assemble_any<stencil::CC>(m_A, m_F, eq, sigma, m_mesh->detJ, C); break;
        case stencil::XF: assemble_any<stencil::XF>(m_A, m_F, eq, sigma, m_mesh->detJ, C); break;
        case stencil::YF: assemble_any<stencil::YF>(m_A, m_F, eq, sigma, m_mesh->detJ, C); break;
        default:          assemble_any<stencil::ZF>(m_A, m_F, eq, sigma, m_mesh->detJ, C); break;
        }
    }
}

void NestFaceFlux::set_defect (const std::vector<JointMember*>& members, const State& qr,
                               const std::vector<DefectRegion>& regions, ResidualFields& R)
{
    JointMember& P = *members[0];
    if (P.defect == nullptr || !P.probed_residual) {
        amrex::Abort("gust: NestFaceFlux::set_defect: the outer domain has no defect source or no probed residual");
    }
    if (regions.size() + 1 != members.size()) {
        amrex::Abort("gust: NestFaceFlux::set_defect: one region per nest");
    }
    m_Rf.resize(regions.size());
    std::vector<FaceFluxNest> nests;
    for (std::size_t n = 1; n < members.size(); ++n) {
        JointMember& N = *members[n];
        if (!N.probed_residual || N.parent == nullptr) {
            amrex::Abort("gust: NestFaceFlux::set_defect: " + N.name + " has no probed residual or parent images");
        }
        ResidualFields& Rf = m_Rf[n - 1];
        if (Rf.u.empty()) { Rf.define(N.state->p.boxArray(), N.state->p.DistributionMap()); }
        nests.push_back({&regions[n - 1], N.geom, N.parent->ratio(), N.parent->offset(),
                         [&N, &Rf] (FluxProbe& p) { N.probed_residual(*N.state, Rf, &p); }});
    }
    compute([&] (FluxProbe& p) { P.probed_residual(qr, R, &p); }, nests);
    P.defect->set(R, regions);
    P.defect->add(m_C);
    if (!m_logged) {
        m_logged = true;
        const auto c = max_norms();
        char line[300];
        std::snprintf(line, sizeof(line), "gust: joint run: nest.feedback_flux = nest: the outer row just outside each "
                      "source region takes the nests' face fluxes; first max |correction| u %.3e v %.3e w %.3e theta %.3e "
                      "k %.3e\n", double(c[0]), double(c[1]), double(c[2]), double(c[3]), double(c[4]));
        amrex::Print() << line;
    }
}

std::vector<amrex::Real> NestFaceFlux::max_norms () const
{
    return {m_C.u.norm0(), m_C.v.norm0(), m_C.w.norm0(), m_C.theta.norm0(), m_C.tke.norm0()};
}

std::string NestFaceFlux::log_suffix () const
{
    const auto c = max_norms();
    char line[160];
    std::snprintf(line, sizeof(line), "; ring max |C| u %.3e v %.3e w %.3e theta %.3e k %.3e", double(c[0]), double(c[1]),
                  double(c[2]), double(c[3]), double(c[4]));
    return line;
}

} // namespace gust
