#include "Buildings.H"

#include "Buildings_K.H"
#include "Sha256.H"
#include "SurfaceRoughness.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

namespace gust {

namespace {

const char* const kBuildingsKeys[] = {"type", "file", "format", "small_volfrac", "z0", "wall_form", "wall_tke", "wall_speed"};
constexpr const char* kBuildingsMagic = "gust_buildings";
constexpr int kBuildingsVersion = 1;

bool has (const std::vector<std::string>& v, const std::string& s)
{
    return std::find(v.begin(), v.end(), s) != v.end();
}

std::string hexfloat (double v)
{
    char b[64];
    std::snprintf(b, sizeof(b), "%a", v);
    return b;
}

bool parse_number (const std::string& tok, double& v)
{
    if (tok.empty()) { return false; }
    char* end = nullptr;
    v = std::strtod(tok.c_str(), &end);
    return *end == '\0';
}

void bcast (std::string& s)
{
    long n = static_cast<long>(s.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

// One decision for every rank: the file opens on the I/O rank.
bool readable (const std::string& file)
{
    int ok = amrex::ParallelDescriptor::IOProcessor() ? (std::ifstream(file).good() ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&ok, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return ok != 0;
}

// The tile edges of points x: halfway between them, the end tiles as far
// outward as inward.
std::vector<double> tile_edges (const std::vector<double>& x)
{
    const std::size_t n = x.size();
    std::vector<double> e(n + 1);
    for (std::size_t i = 1; i < n; ++i) { e[i] = 0.5 * (x[i - 1] + x[i]); }
    e[0] = x[0] - (e[1] - x[0]);
    e[n] = x[n - 1] + (x[n - 1] - e[n - 1]);
    return e;
}

// The file's tokens into the map; "" or what is wrong with them.
std::string parse_map (const std::string& file, std::vector<double>& x, std::vector<double>& y, std::vector<double>& h)
{
    std::ifstream is(file);
    if (!is.good()) { return "cannot be read"; }
    std::vector<std::string> tok;
    std::string t;
    while (is >> t) { tok.push_back(t); }
    double nxd = 0, nyd = 0;
    if (tok.size() < 2 || !parse_number(tok[0], nxd) || !parse_number(tok[1], nyd)
        || nxd != std::floor(nxd) || nyd != std::floor(nyd)) {
        return "must start with the integers nx and ny";
    }
    if (!(nxd >= 2 && nyd >= 2 && nxd <= 1e6 && nyd <= 1e6)) {
        return "nx and ny must be between 2 and 1000000 (got " + tok[0] + " and " + tok[1] + ")";
    }
    const auto nx = static_cast<std::size_t>(nxd), ny = static_cast<std::size_t>(nyd);
    const std::size_t want = 2 + nx + ny + nx * ny;
    if (tok.size() != want) {
        return "holds " + std::to_string(tok.size()) + " values; nx = " + tok[0] + " and ny = " + tok[1] + " need "
               + std::to_string(want);
    }
    std::vector<double> v(want - 2);
    for (std::size_t n = 2; n < want; ++n) {
        if (!parse_number(tok[n], v[n - 2])) { return "value " + std::to_string(n + 1) + " (" + tok[n] + ") is not a number"; }
    }
    x.assign(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(nx));
    y.assign(v.begin() + static_cast<std::ptrdiff_t>(nx), v.begin() + static_cast<std::ptrdiff_t>(nx + ny));
    h.assign(v.begin() + static_cast<std::ptrdiff_t>(nx + ny), v.end());
    auto increasing = [] (const std::vector<double>& c) {
        for (std::size_t i = 0; i < c.size(); ++i) {
            if (!std::isfinite(c[i]) || (i > 0 && !(c[i] > c[i - 1]))) { return false; }
        }
        return true;
    };
    if (!increasing(x)) { return "the x values must be finite and strictly increasing"; }
    if (!increasing(y)) { return "the y values must be finite and strictly increasing"; }
    for (std::size_t n = 0; n < h.size(); ++n) {
        if (!std::isfinite(h[n]) || !(h[n] >= 0)) {
            return "height " + std::to_string(n + 1) + " (" + std::to_string(h[n]) + ") must be finite and >= 0";
        }
    }
    return "";
}

// "x y height [zstart]" lines into the grid; "" or what is wrong with them.
std::string parse_xyz (const std::string& file, std::vector<double>& x, std::vector<double>& y, std::vector<double>& h,
                       std::vector<double>& zs)
{
    std::ifstream is(file);
    if (!is.good()) { return "cannot be read"; }
    std::vector<std::array<double, 4>> pts;
    std::string line;
    for (long ln = 1; std::getline(is, line); ++ln) {
        line = line.substr(0, line.find('#'));
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ls(line);
        std::vector<double> v;
        std::string t;
        double d = 0;
        while (ls >> t) {
            if (!parse_number(t, d)) { return "line " + std::to_string(ln) + ": " + t + " is not a number"; }
            v.push_back(d);
        }
        if (v.empty()) { continue; }
        if (v.size() != 3 && v.size() != 4) {
            return "line " + std::to_string(ln) + " holds " + std::to_string(v.size()) + " values; x y height [zstart] are 3 or 4";
        }
        const double start = v.size() == 4 ? v[3] : 0.0;
        for (const double c : v) {
            if (!std::isfinite(c)) { return "line " + std::to_string(ln) + " holds a value that is not finite"; }
        }
        if (!(v[2] >= 0)) { return "line " + std::to_string(ln) + ": height must be finite and >= 0"; }
        if (start < 0 || (v[2] > 0 && start >= v[2])) {
            return "line " + std::to_string(ln) + ": zstart must be >= 0 and below the height";
        }
        pts.push_back({v[0], v[1], v[2], start});
    }
    x.clear();
    y.clear();
    for (const auto& p : pts) {
        x.push_back(p[0]);
        y.push_back(p[1]);
    }
    for (auto* c : {&x, &y}) {
        std::sort(c->begin(), c->end());
        c->erase(std::unique(c->begin(), c->end()), c->end());
    }
    if (x.size() < 2 || y.size() < 2) { return "the points must span at least 2 x and 2 y values"; }
    const std::size_t n = x.size() * y.size();
    if (pts.size() != n) {
        return "holds " + std::to_string(pts.size()) + " points; its " + std::to_string(x.size()) + " x and "
               + std::to_string(y.size()) + " y values need a full grid of " + std::to_string(n);
    }
    h.assign(n, -1);
    zs.assign(n, 0);
    for (const auto& p : pts) {
        const auto i = static_cast<std::size_t>(std::lower_bound(x.begin(), x.end(), p[0]) - x.begin());
        const auto j = static_cast<std::size_t>(std::lower_bound(y.begin(), y.end(), p[1]) - y.begin());
        if (h[i * y.size() + j] >= 0) {
            return "the point x = " + std::to_string(p[0]) + ", y = " + std::to_string(p[1]) + " appears twice";
        }
        h[i * y.size() + j] = p[2];
        zs[i * y.size() + j] = p[3];
    }
    return "";
}

// The node mean (two or four nodes) of each kind of face, the columns' ground
// being the same mean in the lowest layer k0.
struct NodeMeans
{
    amrex::Array4<amrex::Real const> z;
    int k0;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double h4 (int i, int j, int k) const noexcept
    {
        return 0.25 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i + 1, j, k))
                       + static_cast<double>(z(i, j + 1, k)) + static_cast<double>(z(i + 1, j + 1, k)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double hx (int i, int j, int k) const noexcept
    {
        return 0.5 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i, j + 1, k)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double hy (int i, int j, int k) const noexcept
    {
        return 0.5 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i + 1, j, k)));
    }
};

} // namespace

BuildingsOptions read_buildings_options ()
{
    BuildingsOptions opt;
    amrex::ParmParse pp("buildings");
    for (const char* key : kBuildingsKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) {
            opt.parse_errors.emplace_back(std::string("buildings.") + key + " is given more than once");
        } else if (pp.countval(key) != 1) {
            opt.parse_errors.emplace_back(std::string("buildings.") + key + " takes exactly one value");
        }
    }
    pp.query("type", opt.type_name);
    pp.query("file", opt.file);
    pp.query("format", opt.format);
    std::string s;
    if (pp.query("small_volfrac", s) && !parse_number(s, opt.small_volfrac)) {
        opt.parse_errors.emplace_back("buildings.small_volfrac = " + s + " is not a number");
    }
    if (pp.query("z0", s) && !parse_number(s, opt.z0)) {
        opt.parse_errors.emplace_back("buildings.z0 = " + s + " is not a number");
    }
    pp.query("wall_form", opt.wall_form_name);
    if (opt.wall_form_name == "drag") { opt.wall_form = BuildingWallForm::drag; }
    if (opt.wall_form_name == "roof_stress") { opt.wall_form = BuildingWallForm::roof_stress; }
    pp.query("wall_speed", opt.wall_speed_name);
    opt.wall_speed = opt.wall_speed_name.empty()
                         ? default_wall_speed(opt.wall_form)
                         : (opt.wall_speed_name == "all_faces" ? BuildingWallSpeed::all_faces : BuildingWallSpeed::open_faces);
    pp.query("wall_tke", opt.wall_tke_name);
    opt.wall_tke = opt.wall_tke_name.empty() ? opt.wall_form == BuildingWallForm::fraction_stress : opt.wall_tke_name == "on";
    if (opt.type_name == "immersed_forcing") { opt.type = BuildingsType::immersed_forcing; }
    return opt;
}

std::vector<std::string> check_buildings_options (const BuildingsOptions& opt)
{
    std::vector<std::string> err = opt.parse_errors;
    if (opt.type_name != "none" && opt.type_name != "immersed_forcing") {
        err.emplace_back("buildings.type must be none or immersed_forcing (got " + opt.type_name + ")");
        return err;
    }
    if (opt.type == BuildingsType::none) {
        for (const char* k : kBuildingsKeys) {
            if (std::string(k) != "type" && has(opt.given, k)) {
                err.emplace_back(std::string("buildings.") + k + " has no effect with buildings.type = none");
            }
        }
        return err;
    }
    if (!(opt.z0 > 0) || !std::isfinite(opt.z0)) {
        err.emplace_back("buildings.z0 must be finite and > 0 (got " + std::to_string(opt.z0) + ")");
    }
    if (opt.wall_form_name != "fraction_stress" && opt.wall_form_name != "roof_stress" && opt.wall_form_name != "drag") {
        err.emplace_back("buildings.wall_form must be fraction_stress, roof_stress or drag (got " + opt.wall_form_name + ")");
    }
    if (!opt.wall_tke_name.empty() && opt.wall_tke_name != "on" && opt.wall_tke_name != "off") {
        err.emplace_back("buildings.wall_tke must be on or off (got " + opt.wall_tke_name + ")");
    } else if (opt.wall_tke && opt.wall_form != BuildingWallForm::fraction_stress) {
        err.emplace_back("buildings.wall_tke = on needs buildings.wall_form = fraction_stress (the wall value of k reads the "
                         "wall law's u* on every wall; got " + opt.wall_form_name + ")");
    }
    if (!opt.wall_speed_name.empty() && opt.wall_speed_name != "open_faces" && opt.wall_speed_name != "all_faces") {
        err.emplace_back("buildings.wall_speed must be open_faces or all_faces (got " + opt.wall_speed_name + ")");
    } else if (opt.wall_form == BuildingWallForm::drag && has(opt.given, "wall_speed")) {
        err.emplace_back("buildings.wall_speed has no effect with buildings.wall_form = drag (no wall law)");
    }
    if (opt.wall_form == BuildingWallForm::drag && has(opt.given, "z0")) {
        err.emplace_back("buildings.z0 has no effect with buildings.wall_form = drag (no wall law)");
    }
    if (opt.file.empty()) {
        err.emplace_back("buildings.file is required with buildings.type = immersed_forcing");
    } else if (!readable(opt.file)) {
        err.emplace_back("buildings.file = " + opt.file + ": cannot be read");
    }
    if (opt.format != "erf" && opt.format != "xyz") {
        err.emplace_back("buildings.format must be erf or xyz (got " + opt.format + ")");
    }
    if (!(opt.small_volfrac >= 0 && opt.small_volfrac < 0.5)) {
        err.emplace_back("buildings.small_volfrac must be in [0, 0.5): at 0.5 a cell would be both fluid and solid (got "
                         + std::to_string(opt.small_volfrac) + ")");
    }
    return err;
}

BuildingHeightMap read_building_map (const std::string& file, const std::string& format)
{
    std::vector<double> x, y, h, zs;
    std::string msg;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        msg = format == "xyz" ? parse_xyz(file, x, y, h, zs) : parse_map(file, x, y, h);
        zs.resize(h.size(), 0.0);
    }
    bcast(msg);
    if (!msg.empty()) { amrex::Abort("gust: buildings.file = " + file + ": " + msg); }
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    std::array<long, 3> n{static_cast<long>(x.size()), static_cast<long>(y.size()), static_cast<long>(h.size())};
    amrex::ParallelDescriptor::Bcast(n.data(), 3, root);
    x.resize(static_cast<std::size_t>(n[0]));
    y.resize(static_cast<std::size_t>(n[1]));
    h.resize(static_cast<std::size_t>(n[2]));
    zs.resize(static_cast<std::size_t>(n[2]));
    amrex::ParallelDescriptor::Bcast(x.data(), x.size(), root);
    amrex::ParallelDescriptor::Bcast(y.data(), y.size(), root);
    amrex::ParallelDescriptor::Bcast(h.data(), h.size(), root);
    amrex::ParallelDescriptor::Bcast(zs.data(), zs.size(), root);
    BuildingHeightMap m;
    m.nx = static_cast<int>(x.size());
    m.ny = static_cast<int>(y.size());
    m.xe = tile_edges(x);
    m.ye = tile_edges(y);
    m.h = std::move(h);
    m.zs = std::move(zs);
    return m;
}

std::string Buildings::hash_inputs (const BuildingsOptions& opt)
{
    std::string h;
    if (amrex::ParallelDescriptor::IOProcessor()) { h = sha256_file(opt.file); }
    bcast(h);
    if (h.empty()) { amrex::Abort("gust: buildings.file = " + opt.file + ": cannot be read"); }
    // The wall law's name as the checkpoints know it: "fraction_stress" was the law on roofs
    // and undersides alone before the sides took it, which is roof_stress now. So a checkpoint
    // of that form continues with roof_stress, and is refused with the sides' law. The
    // closure's building terms add a line only when on, so a checkpoint from before them
    // continues with buildings.wall_tke = off. Likewise the open faces' speed: a checkpoint
    // from before it continues with buildings.wall_speed = all_faces.
    const std::string law = opt.wall_form == BuildingWallForm::roof_stress ? "fraction_stress"
                          : opt.wall_form == BuildingWallForm::fraction_stress ? "fraction_stress_sides" : opt.wall_form_name;
    return sha256_hex(std::string(kBuildingsMagic) + " " + std::to_string(kBuildingsVersion) + "\ntype " + opt.type_name
                      + "\nformat " + opt.format + "\nfile " + h + "\nsmall_volfrac " + hexfloat(opt.small_volfrac) + "\nz0 " + hexfloat(opt.z0)
                      + "\nwall_form " + law + "\n" + ((opt.wall_tke && opt.wall_form == BuildingWallForm::fraction_stress) ? "wall_tke on\n" : "")
                      + ((opt.wall_speed == BuildingWallSpeed::open_faces && opt.wall_form != BuildingWallForm::drag) ? "wall_speed open_faces\n" : ""));
}

void Buildings::build (const BuildingsOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                       const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm)
{
    if (!opt.active()) { amrex::Abort("gust: Buildings::build needs buildings.type = immersed_forcing"); }
    m_active = true;
    m_hash = hash_inputs(opt);
    const BuildingHeightMap map = read_building_map(opt.file, opt.format);
    const double lx = geom.isPeriodic(0) ? static_cast<double>(geom.ProbLength(0)) : 0.0;
    const double ly = geom.isPeriodic(1) ? static_cast<double>(geom.ProbLength(1)) : 0.0;
    for (int n = 0; n < cell_ba.size(); ++n) {
        if (cell_ba[n].length(2) != geom.Domain().length(2)) {
            amrex::Abort("gust: buildings: every box must span the column (the mesh's layout does)");
        }
    }
    // Wider than the period, two images of the tiles would cover a column.
    for (int d = 0; d < 2; ++d) {
        const std::vector<double>& e = d == 0 ? map.xe : map.ye;
        const double l = d == 0 ? lx : ly;
        if (l > 0 && e.back() - e.front() > l * (1 + 1e-12)) {
            amrex::Abort("gust: buildings.file = " + opt.file + ": the tiles span " + std::to_string(e.back() - e.front())
                         + " m in " + (d == 0 ? "x" : "y") + ", wider than the periodic domain (" + std::to_string(l)
                         + " m; a file listing both x = 0 and x = L holds the seam twice)");
        }
        // Only the nearest images are searched.
        const double plo_d = static_cast<double>(geom.ProbLo(d));
        if (l > 0 && (e.front() < plo_d - l || e.back() > plo_d + 2 * l)) {
            amrex::Abort("gust: buildings.file = " + opt.file + ": the tiles lie more than one period from the domain in "
                         + (d == 0 ? "x" : "y"));
        }
    }
    amrex::Gpu::DeviceVector<double> dxe(map.xe.size()), dye(map.ye.size()), dh(map.h.size()), dzs(map.zs.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, map.xe.begin(), map.xe.end(), dxe.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, map.ye.begin(), map.ye.end(), dye.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, map.h.begin(), map.h.end(), dh.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, map.zs.begin(), map.zs.end(), dzs.begin());
    const HeightMapView view{dxe.data(), dye.data(), dh.data(), dzs.data(), map.nx, map.ny, lx, ly};

    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const double dxd = static_cast<double>(dx[0]), dyd = static_cast<double>(dx[1]);
    const double x0 = static_cast<double>(plo[0]), y0 = static_cast<double>(plo[1]);
    // At least the round-off of a sum over tiles, so small_volfrac = 0 still gives solid cells.
    const double small = std::max(opt.small_volfrac, 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()));
    const int k0 = geom.Domain().smallEnd(2);

    // The cells. Ghost layers below and above the domain stay 0.
    m_beta.clear();
    m_beta.define(cell_ba, dm, 1, 1);
    m_beta.setVal(0);
    for (amrex::MFIter mfi(m_beta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto b = m_beta.array(mfi);
        const NodeMeans z{z_nd.const_array(mfi), k0};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double g = z.h4(i, j, k0);
            const VerticalOverlap f{z.h4(i, j, k) - g, z.h4(i, j, k + 1) - g};
            const double m = area_mean(view, x0 + i * dxd, x0 + (i + 1) * dxd, y0 + j * dyd, y0 + (j + 1) * dyd, f);
            b(i, j, k) = static_cast<amrex::Real>(clip_fraction(m, small));
        });
    }
    amrex::Gpu::streamSynchronize();
    // A ghost layer in z stays 0 at the domain's ends.
    image_ghost_columns(m_beta, geom);

    // Nothing in the domain, or a building through the top layer: no roof for the wall law.
    const int ktop = geom.Domain().bigEnd(2);
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> rop;
    amrex::ReduceData<amrex::Real, amrex::Real> rdata(rop);
    for (amrex::MFIter mfi(m_beta, false); mfi.isValid(); ++mfi) {
        const auto b = m_beta.const_array(mfi);
        rop.eval(mfi.validbox(), rdata, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real, amrex::Real>
        {
            return {b(i, j, k), k == ktop ? b(i, j, k) : amrex::Real(0)};
        });
    }
    const auto rmax = rdata.value(rop);
    amrex::Real any = amrex::get<0>(rmax), top = amrex::get<1>(rmax);
    amrex::ParallelDescriptor::ReduceRealMax(any);
    amrex::ParallelDescriptor::ReduceRealMax(top);
    if (!(any > 0)) {
        amrex::Abort("gust: buildings.file = " + opt.file + ": no building lies in the domain (above buildings.small_volfrac)");
    }
    if (top > 0) {
        amrex::Abort("gust: buildings.file = " + opt.file + ": a building reaches the top layer of the mesh; "
                     "raise the top of the domain or lower the building");
    }

    // The faces. The round-off of a height is several ulp of the top of the mesh.
    const double ztol = 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon())
                        * static_cast<double>(z_nd.max(0));
    for (int d = 0; d < 3; ++d) {
        amrex::IntVect nodal(0);
        nodal[d] = 1;
        auto& bf = m_beta_face[static_cast<std::size_t>(d)];
        bf.clear();
        bf.define(amrex::convert(cell_ba, nodal), dm, 1, 0);
        bf.setVal(0);
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        const int ei = e[0], ej = e[1], ek = e[2];
        for (amrex::MFIter mfi(bf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto a = bf.array(mfi);
            const auto b = m_beta.const_array(mfi);
            const NodeMeans z{z_nd.const_array(mfi), k0};
            const double eps = std::max(1e-9 * (d == 0 ? dxd : dyd), 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon())
                                                                   * static_cast<double>(geom.ProbLength(d == 0 ? 0 : 1)));
            auto finish = [=] AMREX_GPU_DEVICE (int i, int j, int k, double m) noexcept
            {
                const double lo = static_cast<double>(b(i - ei, j - ej, k - ek)), hi = static_cast<double>(b(i, j, k));
                if (lo == 1 || hi == 1) { return amrex::Real(1); }
                if (lo == 0 && hi == 0) { return amrex::Real(0); }
                return static_cast<amrex::Real>(clip_fraction(m, small));
            };
            if (d == 0) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const double g = z.hx(i, j, k0);
                    const VerticalOverlap f{z.hx(i, j, k) - g, z.hx(i, j, k + 1) - g};
                    const double m = line_mean(view.xe, view.nx, view.lx, view.ye, view.ny, view.ly, view.h, view.zs, view.ny, 1,
                                               x0 + i * dxd, y0 + j * dyd, y0 + (j + 1) * dyd, eps, f);
                    a(i, j, k) = finish(i, j, k, m);
                });
            } else if (d == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const double g = z.hy(i, j, k0);
                    const VerticalOverlap f{z.hy(i, j, k) - g, z.hy(i, j, k + 1) - g};
                    const double m = line_mean(view.ye, view.ny, view.ly, view.xe, view.nx, view.lx, view.h, view.zs, 1, view.ny,
                                               y0 + j * dyd, x0 + i * dxd, x0 + (i + 1) * dxd, eps, f);
                    a(i, j, k) = finish(i, j, k, m);
                });
            } else {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const HorizontalBlock f{z.h4(i, j, k) - z.h4(i, j, k0), ztol};
                    const double m = area_mean(view, x0 + i * dxd, x0 + (i + 1) * dxd, y0 + j * dyd, y0 + (j + 1) * dyd, f);
                    a(i, j, k) = finish(i, j, k, m);
                });
            }
        }
        amrex::Gpu::streamSynchronize();
        // Every copy of a face on a box edge or a periodic seam holds one value.
        bf.OverrideSync(geom.periodicity());
    }

    // The wall cells, from the fractions of the cells' faces.
    const double wall_tol = small;
    m_wall.clear();
    m_wall.define(cell_ba, dm, kNumWallComps, 1);
    m_wall.setVal(0);
    for (amrex::MFIter mfi(m_wall, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto w = m_wall.array(mfi);
        const auto br = m_beta.const_array(mfi);
        const auto fx = m_beta_face[0].const_array(mfi);
        const auto fy = m_beta_face[1].const_array(mfi);
        const auto fz = m_beta_face[2].const_array(mfi);
        const NodeMeans z{z_nd.const_array(mfi), k0};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double own = static_cast<double>(br(i, j, k));
            if (own == 1) { return; }
            const double xl = static_cast<double>(fx(i, j, k)), xh = static_cast<double>(fx(i + 1, j, k));
            const double yl = static_cast<double>(fy(i, j, k)), yh = static_cast<double>(fy(i, j + 1, k));
            const double zl = static_cast<double>(fz(i, j, k)), zh = static_cast<double>(fz(i, j, k + 1));
            int mask = 0;
            if (wall_on_low_side(zl, zh, own, wall_tol)) { mask |= kWallRoof; }
            if (wall_on_low_side(zh, zl, own, wall_tol)) { mask |= kWallUnder; }
            if (wall_on_low_side(xl, xh, own, wall_tol)) { mask |= kWallXlo; }
            if (wall_on_low_side(xh, xl, own, wall_tol)) { mask |= kWallXhi; }
            if (wall_on_low_side(yl, yh, own, wall_tol)) { mask |= kWallYlo; }
            if (wall_on_low_side(yh, yl, own, wall_tol)) { mask |= kWallYhi; }
            w(i, j, k, kWallMask) = static_cast<amrex::Real>(mask);
            // Reference cells inside the domain: the ghost layer of beta holds the sides' images.
            if ((mask & (kWallXlo | kWallXhi)) != 0) {
                const double ref = static_cast<double>((mask & kWallXlo) != 0 ? br(i + 1, j, k) : br(i - 1, j, k));
                w(i, j, k, kWallD1x) = static_cast<amrex::Real>(wall_d1(own, dxd));
                w(i, j, k, kWallD2x) = static_cast<amrex::Real>(wall_d2(own, dxd, dxd, ref, true));
            }
            if ((mask & (kWallYlo | kWallYhi)) != 0) {
                const double ref = static_cast<double>((mask & kWallYlo) != 0 ? br(i, j + 1, k) : br(i, j - 1, k));
                w(i, j, k, kWallD1y) = static_cast<amrex::Real>(wall_d1(own, dyd));
                w(i, j, k, kWallD2y) = static_cast<amrex::Real>(wall_d2(own, dyd, dyd, ref, true));
            }
            if ((mask & (kWallRoof | kWallUnder)) != 0) {
                const bool roof = (mask & kWallRoof) != 0;
                const int kr = roof ? k + 1 : k - 1;
                const bool ok = kr >= k0 && kr <= ktop;
                const double dz = z.h4(i, j, k + 1) - z.h4(i, j, k);
                const double dr = ok ? z.h4(i, j, kr + 1) - z.h4(i, j, kr) : 0.0;
                w(i, j, k, kWallD1z) = static_cast<amrex::Real>(wall_d1(own, dz));
                w(i, j, k, kWallD2z) = static_cast<amrex::Real>(wall_d2(own, dz, dr, ok ? static_cast<double>(br(i, j, kr)) : 0.0, ok));
            }
        });
    }
    amrex::Gpu::streamSynchronize();
    image_ghost_columns(m_wall, geom);
}

std::string Buildings::summary (const amrex::MultiFab& detJ, const amrex::Geometry& geom) const
{
    if (!m_active) { return ""; }
    const auto dx = geom.CellSizeArray();
    const amrex::Real cell = dx[0] * dx[1] * dx[2];
    const amrex::Box& dom = geom.Domain();
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> rop;
    amrex::ReduceData<amrex::Long, amrex::Long, amrex::Long, amrex::Long, amrex::Long, amrex::Real> rdata(rop);
    for (amrex::MFIter mfi(m_beta, false); mfi.isValid(); ++mfi) {
        const auto b = m_beta.const_array(mfi);
        const auto w = m_wall.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        rop.eval(mfi.validbox(), rdata, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 -> amrex::GpuTuple<amrex::Long, amrex::Long, amrex::Long, amrex::Long, amrex::Long, amrex::Real>
        {
            const bool edge = b(i, j, k) > 0 && ((!perx && (i == dom.smallEnd(0) || i == dom.bigEnd(0)))
                                                 || (!pery && (j == dom.smallEnd(1) || j == dom.bigEnd(1))));
            return {b(i, j, k) == amrex::Real(1) ? 1 : 0, (b(i, j, k) > 0 && b(i, j, k) < 1) ? 1 : 0,
                    w(i, j, k, kWallMask) > 0 ? 1 : 0, edge ? 1 : 0,
                    (b(i, j, k) > 0 && b(i, j, k) < 1 && w(i, j, k, kWallMask) == amrex::Real(0)) ? 1 : 0, b(i, j, k) * J(i, j, k) * cell};
        });
    }
    const auto r = rdata.value(rop);
    amrex::Long solid = amrex::get<0>(r), partial = amrex::get<1>(r), walls = amrex::get<2>(r), onface = amrex::get<3>(r);
    amrex::Long thin = amrex::get<4>(r);
    amrex::Real volume = amrex::get<5>(r);
    amrex::ParallelAllReduce::Sum(thin, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(solid, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(partial, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(walls, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(onface, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(volume, amrex::ParallelDescriptor::Communicator());
    std::ostringstream os;
    os.precision(6);
    os << solid << " solid and " << partial << " partial cells, " << walls << " wall cells, solid volume " << volume << " m^3";
    if (onface > 0) {
        os << "\ngust: WARNING: buildings: " << onface << " cells with a building lie on a non-periodic side face of the "
           << "domain; keep buildings away from inflow and outflow faces";
    }
    if (thin > 0) {
        os << "\ngust: WARNING: buildings: " << thin << " partial cells have no wall: a building thinner than a cell "
           << "has opposite sides in one cell; refine the mesh around it";
    }
    return os.str();
}

void Buildings::write_checkpoint (const std::string& dir) const
{
    if (m_active && amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = buildings_checkpoint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << kBuildingsMagic << " " << kBuildingsVersion << "\nsha256 " << m_hash << "\n";
    }
    amrex::ParallelDescriptor::Barrier();
}

void Buildings::check_restart (const std::string& dir, const std::string& hash)
{
    const bool active = !hash.empty();
    // The checkpoint's hash ("" without a file; "?" when unreadable).
    std::string rec;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(buildings_checkpoint_file(dir));
        if (is) {
            std::string magic, key;
            int version = 0;
            if (!(is >> magic >> version >> key >> rec && magic == kBuildingsMagic && version == kBuildingsVersion && key == "sha256")) {
                rec = "?";
            }
        }
    }
    bcast(rec);
    const std::string where = "gust: io.restart = " + dir + ": ";
    if (rec == "?") { amrex::Abort(where + "its Buildings file cannot be read"); }
    if (active && rec.empty()) {
        amrex::Abort(where + "the checkpoint was written without buildings; a restart cannot add them (buildings.type = none, "
                             "or start a new run)");
    }
    if (!active && !rec.empty()) {
        amrex::Abort(where + "the checkpoint was written with buildings; give the run's buildings.* inputs, or start a new run");
    }
    if (active && rec != hash) {
        amrex::Abort(where + "the buildings inputs differ from the checkpoint's (buildings.type, buildings.format, buildings.file contents, "
                             "buildings.small_volfrac, buildings.z0, buildings.wall_form, buildings.wall_tke, buildings.wall_speed; a checkpoint written with "
                             "fraction_stress before the sides had the wall law continues with buildings.wall_form = roof_stress, "
                             "one from before k had building terms with buildings.wall_tke = off, and one from before the wall law "
                             "skipped fixed faces with buildings.wall_speed = all_faces)");
    }
}

} // namespace gust
