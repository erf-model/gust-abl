#include "BoundaryConditions.H"
#include "BoundaryConditions_K.H"

#include "MassFlux_K.H"

#include "CVFluxes.H"
#include "Heights_K.H"
#include "ImplicitOperator.H"
#include "Stencil_K.H"


#include <AMReX_BoxList.H>
#include <AMReX_ParmParse.H>
#include "NestFaces.H"

#include <AMReX_TagParallelFor.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

namespace gust {

using namespace bc_detail;

namespace {

//! How a ghost is filled along a face normal (see the header).
enum class Rule { even, odd, value, profile, gradient };

//! The rule for field f at a face of type t normal to direction d.
Rule rule_for (Field f, FaceType t, int d, bool lo)
{
    const bool normal_vel = (f == f_u && d == 0) || (f == f_v && d == 1) || (f == f_w && d == 2);
    const bool vel = (f == f_u || f == f_v || f == f_w);
    switch (t) {
    case FaceType::inflow:
        return (f == f_p) ? Rule::even : Rule::profile;
    case FaceType::nest:
        // As the inflow: the ghosts hold given values (the parent's; NestFaces
        // fills them), and the normal velocity's face point is a Dirichlet row.
        return (f == f_p) ? Rule::even : Rule::profile;
    case FaceType::outflow:
        return (f == f_p) ? Rule::odd : Rule::even;
    case FaceType::slip_wall:
        if (normal_vel) { return Rule::value; }
        // theta's gradient bc.zhi_theta_grad is the top's; a slip-wall ground is even.
        return (f == f_theta && d == 2 && !lo) ? Rule::gradient : Rule::even;
    case FaceType::no_slip_wall:
        // Every velocity component takes the wall's value: 0, or the moving
        // wall's tangential velocity (bc.<face>_wall_velocity).
        return vel ? Rule::value : Rule::even;
    case FaceType::surface_layer:
        if (normal_vel) { return Rule::value; }
        return (vel || f == f_theta) ? Rule::gradient : Rule::even;
    case FaceType::geostrophic:
        if (normal_vel) { return Rule::value; }
        if (f == f_u || f == f_v) { return Rule::value; }
        return (f == f_theta) ? Rule::gradient : Rule::even;
    default:
        return Rule::even;
    }
}

//! The inputs, after the checks (every message, each naming its key); the
//! inflow profile is built only from checked inputs.
const BoundaryInputs& checked (const BoundaryInputs& in, const amrex::Geometry& geom, const Mesh& mesh,
                               const CoriolisOptions& coriolis)
{
    if (mesh.zlevels().empty()) { amrex::Abort("gust: BoundaryConditions: the mesh has no levels (not built)"); }
    const auto err = check_boundary_inputs(in, geom, coriolis, mesh.zlevels().back());
    if (!err.empty()) {
        std::string all = "gust: invalid boundary conditions:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return in;
}

} // namespace

BoundaryConditions::BoundaryConditions (const BoundaryInputs& in, const amrex::Geometry& geom, const Mesh& mesh,
                                        const CoriolisOptions& coriolis)
    : m_in(checked(in, geom, mesh, coriolis)), m_geom(geom), m_mesh(mesh), m_coriolis(coriolis), m_profile(m_in.inflow, m_in.theta_ref)
{
    // Known only once the table is read: with its own k column the length
    // changes nothing (the fill, the damping and the first guess take k, not
    // the profile's nu_t), so the key would be read and ignored.
    if (m_profile.has_k() && std::find(m_in.given.begin(), m_in.given.end(), "inflow_tke_length") != m_in.given.end()) {
        amrex::Abort("gust: invalid boundary conditions:\n  bc.inflow_tke_length has no effect with a k column in "
                     "bc.inflow_file = " + m_in.inflow.file);
    }
}

bool BoundaryConditions::has_inflow () const
{
    for (const FaceType t : m_in.face) { if (t == FaceType::inflow) { return true; } }
    return false;
}

namespace {

//! One field's ghosts at one face of one box: both ghost layers (and the face
//! point of a normal velocity with a value or profile), so that a direction's
//! fills of every field, both faces and all boxes of the rank are one launch
//! (fill_direction). A point reads only points inside the face (for the
//! gradient rule the column's first point inside, and its own first ghost),
//! never one that another point of the launch writes, when the domain is
//! wider than the ghost layers along d (fill_direction checks that). The order
//! within the launch then does not matter, and the values are those of the
//! fill with one launch per field, face, box and layer (before audit A9).
struct FillTag
{
    amrex::Array4<amrex::Real> a;
    NodeZ nz;
    amrex::Array4<amrex::Real const> g;   //!< the wall law's gradient slab (gradient rule, per column)
    amrex::Box bx;
    amrex::Dim3 glo, ghi;                 //!< the bounds of g
    amrex::Real value;                    //!< the rule's value; the gradient rule's constant
    int loc;                              //!< location type (stencil::CC, XF, YF, ZF)
    int f;                                //!< field (Field)
    int d;                                //!< the face's normal direction
    int rule;                             //!< Rule
    int edge;                             //!< face point (nodal in d) or edge cell index
    int mirror;                           //!< the mirror index about the face (value, even, odd)
    int lo;                               //!< 1 on the low face
    int comp;                             //!< g's component
    int klo;                              //!< the ground slab's k
    int per_column;                       //!< the gradient from g (1) or the constant (0)

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Box const& box () const noexcept { return bx; }
};

//! Fill ghost point (i, j, k) of tag t, or for the gradient rule both ghosts
//! of column (i, j); see FillTag.
template <int T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fill_point (const FillTag& t, int i, int j, int k, const InflowProfile::View& prof, bool homogeneous) noexcept
{
    const auto a = t.a;
    if (t.rule == static_cast<int>(Rule::gradient)) {
        // The gradient along z (the ground or the top) from the column's
        // first point inside; the box is the first ghost layer, the second
        // is written from the first. The slab g holds the two ghost columns
        // the box spans; a face point of a nodal direction averages the two
        // cells beside it, and at the outermost ghost faces (one beyond the
        // columns) the one it has.
        const int c1 = t.lo ? t.edge + 1 : t.edge - 1;
        const int g1 = t.lo ? t.edge - 1 : t.edge + 1, g2 = t.lo ? t.edge - 2 : t.edge + 2;
        amrex::Real G = t.value;
        if (t.per_column) {
            const int ia = amrex::max(t.glo.x, amrex::min(t.ghi.x, i)), ja = amrex::max(t.glo.y, amrex::min(t.ghi.y, j));
            const int ib = amrex::max(t.glo.x, amrex::min(t.ghi.x, stencil::nodal(T, 0) ? i - 1 : i));
            const int jb = amrex::max(t.glo.y, amrex::min(t.ghi.y, stencil::nodal(T, 1) ? j - 1 : j));
            G = amrex::Real(0.5) * (t.g(ib, jb, t.klo, t.comp) + t.g(ia, ja, t.klo, t.comp));
        }
        if (homogeneous) { G = amrex::Real(0); }
        const amrex::Real zc1 = point_z<T>(t.nz, i, j, c1), zg1 = point_z<T>(t.nz, i, j, g1), zg2 = point_z<T>(t.nz, i, j, g2);
        a(i, j, g1) = a(i, j, c1) + (zg1 - zc1) * G;
        a(i, j, g2) = a(i, j, g1) + (zg2 - zg1) * G;
        return;
    }
    const amrex::Real V = homogeneous ? amrex::Real(0) : t.value;
    const int n = (t.d == 0) ? i : (t.d == 1) ? j : k;
    const bool is_profile = (t.rule == static_cast<int>(Rule::profile));
    if (n == t.edge) {
        // The face point of a normal velocity takes the value itself.
        a(i, j, k) = (is_profile && !homogeneous) ? profile_value(prof, t.f, height_above_ground<T>(t.nz, i, j, k, t.klo)) : V;
        return;
    }
    if (is_profile) {
        a(i, j, k) = homogeneous ? amrex::Real(0) : profile_value(prof, t.f, height_above_ground<T>(t.nz, i, j, k, t.klo));
        return;
    }
    // The mirror point about the face (nodal) or the face plane (cell).
    const int ex = (t.d == 0) ? 1 : 0, ey = (t.d == 1) ? 1 : 0, ez = (t.d == 2) ? 1 : 0;
    const int m = t.mirror - n;
    const amrex::Real am = a(i + ex * (m - i), j + ey * (m - j), k + ez * (m - k));
    if (t.rule == static_cast<int>(Rule::value)) {
        a(i, j, k) = 2 * V - am;
    } else {
        const int sgn = (t.rule == static_cast<int>(Rule::odd)) ? -1 : 1;
        a(i, j, k) = sgn * am;
    }
}

//! The tags of field f (location type T) at face (d, lo), one per box
//! touching it whose slab is not empty: the rule of the face's type, its value (the geostrophic
//! wind, the moving wall's velocity) and, for the gradient rule, the wall
//! law's gradient (surface layer) or bc.zhi_theta_grad.
template <int T>
void add_field_tags (amrex::Vector<FillTag>& tags, amrex::MultiFab& q, int f, int d, bool lo, const BoundaryInputs& in,
                     const amrex::MultiFab& z_nd, const amrex::Geometry& geom, const CoriolisOptions& cor,
                     const amrex::MultiFab* ground_grad)
{
    const amrex::Orientation o(d, lo ? amrex::Orientation::low : amrex::Orientation::high);
    const FaceType ft = in.face[o];
    if (ft == FaceType::nest && f != f_p) { return; }   // NestFaces::fill
    const Rule rule = rule_for(static_cast<Field>(f), ft, d, lo);
    amrex::Real value = 0;
    if (ft == FaceType::geostrophic) { value = (f == f_u) ? cor.ug : (f == f_v) ? cor.vg : amrex::Real(0); }
    if (ft == FaceType::no_slip_wall && f <= f_w) { value = in.wall_velocity[o][f]; }
    const amrex::MultiFab* grad = nullptr;
    int comp = 0;
    if (rule == Rule::gradient) {
        if (ft == FaceType::surface_layer) {
            grad = ground_grad;
            comp = (f == f_u) ? 0 : (f == f_v) ? 1 : 2;
        } else {
            value = in.zhi_theta_grad;
        }
    }
    const amrex::Box& dom = geom.Domain();
    const bool nodal_d = stencil::nodal(T, d);
    const int ng = q.nGrow();
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (!touches(vb, dom, d, lo)) { continue; }
        FillTag t{};
        t.a = q.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        t.nz = NodeZ{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        t.loc = T;
        t.f = f;
        t.d = d;
        t.rule = static_cast<int>(rule);
        t.lo = lo ? 1 : 0;
        t.klo = dom.smallEnd(2);
        t.value = value;
        if (rule == Rule::gradient) {
            t.edge = lo ? dom.smallEnd(2) : dom.bigEnd(2);
            t.bx = ghost_slab(vb, dom, geom, 2, lo, 1, ng) & mfi.fabbox();
            if (grad != nullptr) {
                t.g = grad->const_array(mfi);
                t.glo = amrex::lbound(t.g);
                t.ghi = amrex::ubound(t.g);
                t.comp = comp;
                t.per_column = 1;
            }
        } else {
            t.edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
            t.mirror = nodal_d ? 2 * t.edge : (lo ? 2 * t.edge - 1 : 2 * t.edge + 1);
            // Layers 1 to ng, and the face point (layer 0) of a normal velocity
            // with a value or profile: the slabs differ only along d.
            const int inner = (nodal_d && (rule == Rule::value || rule == Rule::profile)) ? 0 : 1;
            amrex::Box gb = ghost_slab(vb, dom, geom, d, lo, inner, ng);
            gb.setSmall(d, lo ? t.edge - ng : t.edge + inner);
            gb.setBig(d, lo ? t.edge - inner : t.edge + ng);
            t.bx = gb & mfi.fabbox();
        }
        if (!t.bx.isEmpty()) { tags.push_back(t); }
    }
}

//! Fill the ghosts of every field at both faces normal to d, on every box of
//! the rank, in one launch. When the domain is not wider than the ghost
//! layers along d, one face's ghosts read the other face's points (on two
//! cells, u(-2) mirrors the high face point u(2)), so the fill is two
//! launches, low face first, as before audit A9.
void fill_direction (State& q, int d, const BoundaryInputs& in, const amrex::Geometry& geom, const InflowProfile::View& prof,
                     const CoriolisOptions& cor, const amrex::MultiFab* ground_grad, bool homogeneous)
{
    amrex::Vector<FillTag> tags;
    auto add = [&] (bool lo) {
        add_field_tags<stencil::XF>(tags, q.u, f_u, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
        add_field_tags<stencil::YF>(tags, q.v, f_v, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
        add_field_tags<stencil::ZF>(tags, q.w, f_w, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
        add_field_tags<stencil::CC>(tags, q.p, f_p, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
        add_field_tags<stencil::CC>(tags, q.theta, f_theta, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
        add_field_tags<stencil::CC>(tags, q.tke, f_tke, d, lo, in, q.z_phys_nd, geom, cor, ground_grad);
    };
    auto launch = [&] () {
        if (tags.empty()) { return; }
        amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, FillTag const& t) noexcept
        {
            switch (t.loc) {
            case stencil::XF: fill_point<stencil::XF>(t, i, j, k, prof, homogeneous); break;
            case stencil::YF: fill_point<stencil::YF>(t, i, j, k, prof, homogeneous); break;
            case stencil::ZF: fill_point<stencil::ZF>(t, i, j, k, prof, homogeneous); break;
            default: fill_point<stencil::CC>(t, i, j, k, prof, homogeneous); break;
            }
        });
        tags.clear();
    };
    int ng = 0;
    for (const amrex::MultiFab* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { ng = std::max(ng, mf->nGrow()); }
    if (geom.Domain().length(d) <= ng) {
        add(true);
        launch();
        add(false);
        launch();
    } else {
        // The high face's tags first: a CPU build runs the tags in order (GUST
        // has no OpenMP), so if the split above were lost, the narrow domain's
        // low face would read the high face's new values and
        // BoundaryFill.NarrowDomainFillsTheLowFaceFirst would fail. In the low
        // face's order it would pass; keep this order.
        add(false);
        add(true);
        launch();
    }
}

} // namespace

namespace {
void fill_no_flux_w (State& q, const BoundaryInputs& in, const amrex::Geometry& geom);

//! The wall law's gradient of u, v (components 0, 1) and theta (2) on the
//! valid and first ghost columns of g (the surface state's slab); Z0 the
//! roughness accessor (SurfaceRoughness.H).
template <class Z0>
void ground_gradients (const SurfaceState& ss, const SurfaceLayerOptions& surface_opt, amrex::MultiFab& ground_grad)
{
    // The wall law's gradient per column: u*/(kappa (d0 + z0)) along the
    // tangential wind for u and v (components 0, 1), -q/(kappa u* (d0 +
    // z0)) for theta (component 2, 0 when neutral).
    const amrex::Real kappa = surface_opt.kappa;
    const Roughness rough = ground_roughness(&surface_opt);
    const amrex::Real qs = (surface_opt.stability == SurfaceLayerOptions::Stability::heat_flux)
                               ? surface_opt.heat_flux : amrex::Real(0);
    const StableCoefficients sc = surface_opt.stable_coefficients();
    for (amrex::MFIter mfi(ground_grad, false); mfi.isValid(); ++mfi) {
        const auto g = ground_grad.array(mfi);
        const auto us = ss.ustar.const_array(mfi);
        const auto d0 = ss.d0.const_array(mfi);
        const auto ta = ss.tau.const_array(mfi);
        const auto il = ss.inv_L.const_array(mfi);
        const Z0 zf = Z0::of(rough, mfi);
        amrex::ParallelFor(mfi.growntilebox(amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // MOST's flux-gradient relations, phi_m and phi_h of
            // zeta = d0/L (1 when neutral), with the wall model's clamp
            // and bc.surface_stable_functions' stable side (data, not a branch).
            const amrex::Real zeta = amrex::max(amrex::min(d0(i, j, k) * il(i, j, k), amrex::Real(100)), -amrex::Real(100));
            amrex::Real phi_m, phi_h;
            most_phi(zeta, sc, phi_m, phi_h);
            const amrex::Real s = kappa * (d0(i, j, k) + zf(i, j));
            const amrex::Real tmag = std::sqrt(ta(i, j, k, 0) * ta(i, j, k, 0) + ta(i, j, k, 1) * ta(i, j, k, 1)
                                               + ta(i, j, k, 2) * ta(i, j, k, 2));
            // The unit stress; with no stress (tmag = 0) the gradient is 0.
            // No select or max: the compiler evaluates both sides of one
            // and moves it past the division, and 0/0 would raise the
            // invalid trap although discarded. The added term is below
            // the resolution of any stress (bit-identical quotient).
            const amrex::Real tsafe = tmag + std::numeric_limits<amrex::Real>::min();
            const amrex::Real dir_x = ta(i, j, k, 0) / tsafe;
            const amrex::Real dir_y = ta(i, j, k, 1) / tsafe;
            g(i, j, k, 0) = us(i, j, k) / s * dir_x * phi_m;
            g(i, j, k, 1) = us(i, j, k) / s * dir_y * phi_m;
            g(i, j, k, 2) = -qs / (us(i, j, k) * s) * phi_h;
        });
    }
}

} // namespace

void BoundaryConditions::outflow_transport (amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& F) const
{
    if (!m_in.outflow_mirror_transport) { return; }
    const amrex::Box& dom = m_geom.Domain();
    // x faces, then y faces over the x ghosts they wrote (the corners).
    for (int d = 0; d < 2; ++d) {
        for (const bool lo : {true, false}) {
            if (m_in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)] != FaceType::outflow) {
                continue;
            }
            const int sh = lo ? 1 : -1;   // from the ghost towards the inside
            const int ex = (d == 0) ? 1 : 0, ey = (d == 1) ? 1 : 0;
            for (int c = 0; c < AMREX_SPACEDIM; ++c) {
                amrex::MultiFab& mf = F[c];
                // The normal flux has its face f on the boundary: the ghost faces continue it
                // linearly through the face value, F(f + m) = 2 F(f) - F(f - m) (m counted
                // outwards). The others (cell-centred along d) mirror the cells inside: ghost
                // layer L takes the cell L - 1 inside the face, F(g) = F(g + (2 L - 1) sh).
                const bool normal = (c == d);
                const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;   // nodal index of f (normal flux only)
                for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
                    const amrex::Box vb = mfi.validbox();
                    if (!touches(vb, dom, d, lo)) { continue; }
                    const auto f = mf.array(mfi);
                    for (int layer = 1; layer <= mf.nGrow(); ++layer) {
                        amrex::Box gb = ghost_slab(vb, dom, m_geom, d, lo, layer, mf.nGrow());
                        gb &= mfi.fabbox();
                        if (gb.isEmpty()) { continue; }
                        if (normal) {
                            amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                const int ie = ex ? face : i, je = ey ? face : j;
                                f(i, j, k) = 2 * f(ie, je, k) - f(ex ? 2 * face - i : i, ey ? 2 * face - j : j, k);
                            });
                        } else {
                            const int m = sh * (2 * layer - 1);   // ghost layer L mirrors the cell L - 1 inside the face
                            amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                f(i, j, k) = f(i + ex * m, j + ey * m, k);
                            });
                        }
                    }
                }
            }
        }
    }
}

namespace {

//! set_theta_ground_extension: theta's two ground ghost layers continue the
//! first two cells linearly. A free function: nvcc refuses an extended device
//! lambda in a private member.
void extend_theta_below_ground (amrex::MultiFab& theta, const amrex::MultiFab& z_nd, const amrex::Geometry& geom)
{
    // Every column of the fab (the lateral ghost columns hold their images or
    // the face fills, done above), both ghost layers, from cells klo and klo+1:
    // theta(z) = theta_0 + (z - z_0) (theta_1 - theta_0)/(z_1 - z_0), z the
    // cell-centre heights (the mean of eight nodes; NodeZ clamps laterally).
    const int klo = geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(theta, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().smallEnd(2) != klo) { continue; }
        amrex::Box b = mfi.fabbox();
        b.setRange(2, klo, 1);
        const auto a = theta.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const amrex::Real z0 = point_z<stencil::CC>(nz, i, j, klo), z1 = point_z<stencil::CC>(nz, i, j, klo + 1);
            const amrex::Real slope = (a(i, j, klo + 1) - a(i, j, klo)) / (z1 - z0);
            for (int m = 1; m <= 2; ++m) {
                a(i, j, klo - m) = a(i, j, klo) + (point_z<stencil::CC>(nz, i, j, klo - m) - z0) * slope;
            }
        });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

void BoundaryConditions::fill_ghosts (State& q, const SurfaceLayerOptions* surface_opt, const BuoyancyOptions& buoyancy,
                                      bool homogeneous, SurfaceState* surface_out) const
{
    const bool ground_wall = (m_in.face[2] == FaceType::surface_layer);
    if (ground_wall && !homogeneous && surface_opt == nullptr) {
        amrex::Abort("gust: fill_ghosts: bc.zlo = surface_layer needs the surface-layer options");
    }
    if (m_nest == nullptr && has_nest_face(m_in)) {
        amrex::Abort("gust: fill_ghosts: a nest face (bc.xlo, bc.xhi, bc.ylo or bc.yhi = nest) needs the parent's values "
                     "(set_nest_faces)");
    }
    amrex::MultiFab* fields[n_fields] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    for (auto* mf : fields) {
        if (mf->nGrow() < 2) { amrex::Abort("gust: fill_ghosts: the state needs two ghost cells"); }
        mf->FillBoundary(m_geom.periodicity());
    }
    const InflowProfile::View prof = m_profile.view();
    SurfaceState local;
    SurfaceState* ss = (surface_out != nullptr) ? surface_out : &local;
    amrex::MultiFab ground_grad;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d == 2 && ground_wall && !homogeneous) {
            // The surface layer reads the lateral ghost columns filled above.
            compute_surface_state(q, *surface_opt, buoyancy, m_mesh, m_geom, *ss);
            // Two lateral ghost columns, as the z fills span (the state's
            // ngrow). The first is computed from the surface state's own
            // ghost column; FillBoundary gives every column that images a
            // domain column (the neighbouring boxes and the periodic shifts);
            // what is left (outside a non-periodic face beyond the first
            // ghost column, and the corners there) takes the gradient of the
            // nearest column inside the domain along the non-periodic
            // directions (the state's ghost columns there are the profile's
            // or mirrors; no valid row reads those corner ghosts).
            const amrex::Real unset = std::numeric_limits<amrex::Real>::max();
            ground_grad.define(ss->ustar.boxArray(), ss->ustar.DistributionMap(), 3, amrex::IntVect(2, 2, 0));
            ground_grad.setVal(unset);
            // The wall law's gradient per column (ground_gradients).
            if (surface_opt->z0_map != nullptr) {
                check_roughness_map("fill_ghosts", *surface_opt->z0_map, m_mesh.detJ);
                ground_gradients<MapZ0>(*ss, *surface_opt, ground_grad);
            } else {
                ground_gradients<UniformZ0>(*ss, *surface_opt, ground_grad);
            }
            ground_grad.FillBoundary(m_geom.periodicity());
            {
                const amrex::Box& dom = m_geom.Domain();
                const bool perx = m_geom.isPeriodic(0), pery = m_geom.isPeriodic(1);
                for (amrex::MFIter mfi(ground_grad, false); mfi.isValid(); ++mfi) {
                    const auto g = ground_grad.array(mfi);
                    amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        if (g(i, j, k, 0) != unset) { return; }
                        const int ic = perx ? i : amrex::max(dom.smallEnd(0), amrex::min(dom.bigEnd(0), i));
                        const int jc = pery ? j : amrex::max(dom.smallEnd(1), amrex::min(dom.bigEnd(1), j));
                        for (int n = 0; n < 3; ++n) { g(i, j, k, n) = g(ic, jc, k, n); }
                    });
                }
            }
        }
        if (m_geom.isPeriodic(d)) { continue; }
        const amrex::MultiFab* gg = ground_grad.empty() ? nullptr : &ground_grad;
        // Every field at both faces in one launch (two on a narrow domain,
        // fill_direction); the no-flux w of the ground and the top reads u and
        // v there, so it comes after. The p, theta and k fills do not read w.
        fill_direction(q, d, m_in, m_geom, prof, m_coriolis, gg, homogeneous);
        if (m_nest != nullptr) { m_nest->fill(q, d, homogeneous); }
        if (d == 2) { fill_no_flux_w(q, m_in, m_geom); }
    }
    if (m_theta_ground_extension && !homogeneous) { extend_theta_below_ground(q.theta, q.z_phys_nd, m_geom); }
    for (auto* mf : fields) { sync_ghost_copies(*mf); }
}


const BoundaryConditions::GhostSlabs& BoundaryConditions::ghost_slabs (const amrex::MultiFab& mf) const
{
    for (const GhostSlabs& g : m_ghost_slabs) {
        if (g.ba == mf.boxArray() && g.dm == mf.DistributionMap()) { return g; }
    }
    // A run has a few layouts (four index types per level); a long series of
    // others (tests) starts the cache again.
    if (m_ghost_slabs.size() >= 64) { m_ghost_slabs.clear(); }
    const amrex::Box dom = amrex::convert(m_geom.Domain(), mf.ixType());
    const amrex::IntVect ng = mf.nGrowVect();
    const amrex::BoxArray& ba = mf.boxArray();
    const amrex::DistributionMapping& dm = mf.DistributionMap();
    amrex::BoxList bl(mf.ixType());
    amrex::Vector<int> pmap, owner;
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
        const amrex::Box vb = ba[b];
        amrex::Box ext = vb;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (m_geom.isPeriodic(d)) { continue; }
            if (vb.smallEnd(d) == dom.smallEnd(d)) { ext.growLo(d, ng[d]); }
            if (vb.bigEnd(d) == dom.bigEnd(d)) { ext.growHi(d, ng[d]); }
        }
        for (const amrex::Box& s : amrex::boxDiff(ext, dom)) {
            bl.push_back(s);
            pmap.push_back(dm[b]);
            owner.push_back(b);
        }
    }
    GhostSlabs g;
    g.ba = ba;
    g.dm = dm;
    if (!bl.isEmpty()) {
        g.slabs = amrex::BoxArray(std::move(bl));
        g.slab_dm = amrex::DistributionMapping(std::move(pmap));
    }
    g.owner = std::move(owner);
    m_ghost_slabs.push_back(std::move(g));
    return m_ghost_slabs.back();
}

//! Make every copy of a ghost cell hold one value. A box fills the ghosts of
//! its whole fab, but in its outermost lateral ghost layer some of what the
//! rules read is not in the fab (the heights have one lateral ghost layer, the
//! wall law's gradient two), so the box's value there differed from the one
//! the box owning that column computes. Each cell outside the domain takes the
//! value of the box whose (lateral) column it is in: the box touching that
//! face, the slab outside it, across its valid extent along the face (the
//! periodic images of those in a periodic direction). The ghost images of
//! cells inside the domain are then the neighbours' again (the face fills
//! rewrote a nodal face point in the lateral ghost columns). No valid cell
//! changes. Some of the cells that change are read: the inflow face points
//! and the corners at a box seam where an inflow face crosses terrain.
void BoundaryConditions::sync_ghost_copies (amrex::MultiFab& mf) const
{
    const GhostSlabs& g = ghost_slabs(mf);
    if (!g.slabs.empty()) {
        amrex::MultiFab slab(g.slabs, g.slab_dm, 1, 0, amrex::MFInfo(), mf.Factory());
        for (amrex::MFIter mfi(slab, false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            slab[mfi].copy<amrex::RunOn::Device>(mf[g.owner[mfi.index()]], bx, 0, bx, 0, 1);
        }
        // A nodal face shared at a seam (or with its periodic image) is in two
        // slabs; which one ParallelCopy writes last is not fixed. Both are equal
        // because every update keeps the shared faces synchronised
        // (OverrideSync), which this relies on.
        mf.ParallelCopy(slab, 0, 0, 1, amrex::IntVect(0), mf.nGrowVect(), m_geom.periodicity());
    }
    mf.FillBoundary(m_geom.periodicity());
}

namespace {

//! Whether the normal velocity of direction d has a Dirichlet value on face
//! (d, lo), and that value (profile: NaN marks "the inflow profile").
bool dirichlet_normal (const BoundaryInputs& in, int d, bool lo, FaceType& t)
{
    t = in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
    return t == FaceType::inflow || t == FaceType::slip_wall || t == FaceType::no_slip_wall
        || t == FaceType::surface_layer || t == FaceType::geostrophic;
}

//! The ground and top faces of w take the no-flux value Omega = 0,
//! w = z_x u_f + z_y v_f (MassFlux_K.H; 0 over flat ground), with the two
//! ghost layers mirrored about it (2 V - the interior), after the generic
//! value fill of w in z. u and v are filled (their lateral ghosts) before w.
//! The same formula serves the homogeneous fill (it is linear).
// With the homogeneous fill the value is linear in the perturbation's u, v
// given their homogeneous ghosts (the wall law's gradient of the nonlinear
// fill is not differentiated here; no production caller fills homogeneously).
void fill_no_flux_w (State& q, const BoundaryInputs& in, const amrex::Geometry& geom)
{
    const amrex::Box& dom = geom.Domain();
    const auto dxinv = geom.InvCellSizeArray();
    for (const bool lo : {true, false}) {
        FaceType t;
        if (!dirichlet_normal(in, 2, lo, t)) { continue; }
        const int face = lo ? dom.smallEnd(2) : dom.bigEnd(2) + 1;
        for (amrex::MFIter mfi(q.w, false); mfi.isValid(); ++mfi) {
            amrex::Box fb = mfi.validbox();
            if (fb.smallEnd(2) > face || fb.bigEnd(2) < face) { continue; }
            fb.setSmall(2, face);
            fb.setBig(2, face);
            fb.grow(0, q.w.nGrow());
            fb.grow(1, q.w.nGrow());
            fb &= mfi.fabbox();
            // Only the faces whose four u, v points and four nodal heights
            // exist in this box's fabs (the heights have one lateral ghost
            // layer, w two: the outermost lateral ghost column is set below
            // from its periodic image, and on a non-periodic side keeps the
            // lateral rule's value, which only the Dirichlet rows read).
            const amrex::Box ub = q.u[mfi].box(), vb = q.v[mfi].box(), zb = q.z_phys_nd[mfi].box();
            fb.setSmall(0, amrex::max(fb.smallEnd(0), ub.smallEnd(0), vb.smallEnd(0), zb.smallEnd(0)));
            fb.setBig(0, amrex::min(fb.bigEnd(0), ub.bigEnd(0) - 1, vb.bigEnd(0), zb.bigEnd(0) - 1));
            fb.setSmall(1, amrex::max(fb.smallEnd(1), ub.smallEnd(1), vb.smallEnd(1), zb.smallEnd(1)));
            fb.setBig(1, amrex::min(fb.bigEnd(1), ub.bigEnd(1), vb.bigEnd(1) - 1, zb.bigEnd(1) - 1));
            if (fb.isEmpty()) { continue; }
            const auto w = q.w.array(mfi);
            const auto uu = q.u.const_array(mfi), vv = q.v.const_array(mfi);
            const auto zn = q.z_phys_nd.const_array(mfi);
            const int ng = q.w.nGrow();
            const int sgn = lo ? -1 : 1;
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real V = massflux::no_flux_w(zn, uu, vv, i, j, k, dxinv[0], dxinv[1]);
                w(i, j, k) = V;
                for (int m = 1; m <= ng; ++m) { w(i, j, k + sgn * m) = 2 * V - w(i, j, k - sgn * m); }
            });
        }
    }
    // The periodic images of the face rows (the lateral ghost columns the
    // heights do not reach).
    q.w.FillBoundary(geom.periodicity());
}

} // namespace

void BoundaryConditions::apply_rows (const State& q, ResidualFields& R, bool homogeneous) const
{
    const amrex::Box& dom = m_geom.Domain();
    const int klo = dom.smallEnd(2);
    const InflowProfile::View prof = m_profile.view();
    const amrex::MultiFab* qf[3] = {&q.u, &q.v, &q.w};
    amrex::MultiFab* rf[3] = {&R.u, &R.v, &R.w};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            FaceType t;
            if (!dirichlet_normal(m_in, d, lo, t)) { continue; }
            const bool is_profile = (t == FaceType::inflow) && !homogeneous;
            const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
            for (amrex::MFIter mfi(*rf[d], false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }
                fb.setSmall(d, face);
                fb.setBig(d, face);
                const auto r = rf[d]->array(mfi);
                const auto a = qf[d]->const_array(mfi);
                const amrex::Box zb = q.z_phys_nd[mfi].box();
                const NodeZ nz{q.z_phys_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
                const int f = d;
                if (d == 0) {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real V = is_profile ? profile_value(prof, f, height_above_ground<stencil::XF>(nz, i, j, k, klo)) : amrex::Real(0);
                        r(i, j, k) = a(i, j, k) - V;
                    });
                } else if (d == 1) {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real V = is_profile ? profile_value(prof, f, height_above_ground<stencil::YF>(nz, i, j, k, klo)) : amrex::Real(0);
                        r(i, j, k) = a(i, j, k) - V;
                    });
                } else {
                    // The ground and the top: w minus the no-flux value
                    // (Omega = 0; 0 over flat ground). Linear in (u, v, w):
                    // the homogeneous fill uses the same formula.
                    const auto uu = q.u.const_array(mfi), vv = q.v.const_array(mfi);
                    const auto zn = q.z_phys_nd.const_array(mfi);
                    const auto dxinv = m_geom.InvCellSizeArray();
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        r(i, j, k) = a(i, j, k) - massflux::no_flux_w(zn, uu, vv, i, j, k, dxinv[0], dxinv[1]);
                    });
                }
            }
        }
    }
    if (m_nest != nullptr) { m_nest->apply_rows(q, R, homogeneous); }
}

void BoundaryConditions::dirichlet_rows (Equation eq, Loc loc, const Mesh& mesh, amrex::MultiFab& mask) const
{
    if (mask.empty()) {
        const amrex::MultiFab* like = (loc == Loc::xface) ? &mesh.ax : (loc == Loc::yface) ? &mesh.ay : (loc == Loc::zface) ? &mesh.az : &mesh.detJ;
        mask.define(like->boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    }
    check_operator_field("dirichlet_rows", "mask", mask, loc, mesh, 0);
    mask.setVal(0.0);
    const amrex::Box& dom = m_geom.Domain();
    const int f = (eq == Equation::u) ? f_u : (eq == Equation::v) ? f_v : (eq == Equation::w) ? f_w
                : (eq == Equation::theta) ? f_theta : f_tke;
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            const FaceType ft = m_in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            const Rule rule = rule_for(static_cast<Field>(f), ft, d, lo);
            const bool nodal_d = stencil::nodal(t, d);
            if (!(nodal_d && (rule == Rule::value || rule == Rule::profile))) { continue; }
            const int edge = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
            for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > edge || fb.bigEnd(d) < edge) { continue; }
                fb.setSmall(d, edge);
                fb.setBig(d, edge);
                mask[mfi].template setVal<amrex::RunOn::Device>(amrex::Real(1), fb, 0, 1);
            }
        }
    }
}

void BoundaryConditions::apply_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const
{
    // The damping first: the Dirichlet face points below become identity rows.
    add_damping_to_implicit(eq, loc, z_nd, A);
    const amrex::Box& dom = m_geom.Domain();
    const int f = (eq == Equation::u) ? f_u : (eq == Equation::v) ? f_v : (eq == Equation::w) ? f_w
                : (eq == Equation::theta) ? f_theta : f_tke;
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            const FaceType ft = m_in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            const Rule rule = rule_for(static_cast<Field>(f), ft, d, lo);
            // A field nodal in d is the normal velocity, with a valid point on
            // the face: a Dirichlet value makes that row the identity; an even
            // mirror about the point (outflow) couples it to its other neighbour.
            const bool nodal_d = stencil::nodal(t, d);
            const bool dirichlet_pt = nodal_d && (rule == Rule::value || rule == Rule::profile);
            const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
            const int cm = (d == 0) ? coef_XM : (d == 1) ? coef_YM : coef_ZM;
            const int cp = (d == 0) ? coef_XP : (d == 1) ? coef_YP : coef_ZP;
            const int ghost_coef = lo ? cm : cp, other_coef = lo ? cp : cm;
            // Fold sign: +1 value/odd, -1 even, 0 profile; gradient moves the
            // coefficient to the opposite neighbour.
            amrex::Real fold = 0;
            bool to_other = false;
            if (nodal_d) { to_other = true; }
            else if (rule == Rule::even) { fold = -1; }
            else if (rule == Rule::odd || rule == Rule::value) { fold = 1; }
            else if (rule == Rule::gradient) { to_other = true; }
            for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > edge || fb.bigEnd(d) < edge) { continue; }
                fb.setSmall(d, edge);
                fb.setBig(d, edge);
                const auto a = A.array(mfi);
                if (dirichlet_pt) {
                    // Dirichlet face point: the identity row.
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) = 0; }
                        a(i, j, k, coef_P) = 1;
                    });
                } else {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real ag = a(i, j, k, ghost_coef);
                        a(i, j, k, ghost_coef) = 0;
                        if (to_other) { a(i, j, k, other_coef) += ag; }
                        else { a(i, j, k, coef_P) += fold * ag; }
                    });
                }
            }
        }
    }
}

} // namespace gust
