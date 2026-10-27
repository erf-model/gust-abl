#include "Canopy.H"

#include "Canopy_K.H"
#include "Sha256.H"
#include "SurfaceRoughness.H"
#include "TerrainSurface.H"

#include <ERF_ForestUtils.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>

namespace gust {

namespace {

const char* const kCanopyKeys[] = {"type", "file", "height_file", "lai_file", "cd", "cd_file", "profile", "laimax"};
constexpr const char* kCanopyMagic = "gust_canopy";
constexpr int kCanopyVersion = 1;

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

// One decision for every rank: the file opens on the I/O rank.
bool readable (const std::string& file)
{
    int ok = amrex::ParallelDescriptor::IOProcessor() ? (std::ifstream(file).good() ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&ok, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return ok != 0;
}

void bcast (std::string& s)
{
    long n = static_cast<long>(s.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

// SHA-256 of a file's contents on the I/O rank, on every rank; aborts naming
// key when it cannot be read.
std::string file_hash (const std::string& file, const std::string& key)
{
    std::string h;
    if (amrex::ParallelDescriptor::IOProcessor()) { h = sha256_file(file); }
    bcast(h);
    if (h.empty()) { amrex::Abort("gust: " + key + " = " + file + ": cannot be read"); }
    return h;
}

// A map's values at the column centres, checked (key names the input):
// every value finite and >= 0, at least idw_neighbors points.
TerrainPoints read_map_points (const std::string& file, const std::string& key, int k)
{
    TerrainPoints pts = read_xyz_points(file, key);
    std::vector<std::string> err;
    if (pts.nonfinite > 0) {
        err.push_back(std::to_string(pts.nonfinite) + " lines hold a value that is not finite; every value must be finite and >= 0");
    }
    long bad = 0;
    for (const double v : pts.z) { if (!(v >= 0)) { ++bad; } }
    if (bad > 0) { err.push_back(std::to_string(bad) + " points have a value < 0; every value must be finite and >= 0"); }
    const auto n = static_cast<long>(pts.z.size());
    if (k < 1 || k > n) {
        err.push_back("terrain.idw_neighbors = " + std::to_string(k) + " must be between 1 and the number of points ("
                      + std::to_string(n) + ")");
    }
    if (!err.empty()) {
        std::string all = "gust: " + key + " = " + file + ":";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return pts;
}

// The height above the column's ground of the cell's lower (zb) and upper
// (zt) faces: four-corner means, as in detJ. k0 is the domain's lowest layer.
struct CellHeights
{
    amrex::Array4<amrex::Real const> z;
    int k0;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double face (int i, int j, int k) const noexcept
    {
        return 0.25 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i + 1, j, k))
                       + static_cast<double>(z(i, j + 1, k)) + static_cast<double>(z(i + 1, j + 1, k)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (int i, int j, int k, double& zb, double& zt) const noexcept
    {
        const double g = face(i, j, k0);
        zb = face(i, j, k) - g;
        zt = face(i, j, k + 1) - g;
    }
};

struct UniformFraction
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (double a, double b) const noexcept { return uniform_fraction(a, b); }
};

struct LmFraction
{
    LmProfile p;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (double a, double b) const noexcept { return lm_fraction(a, b, p); }
};

// The map's drag: cols holds each column's height, LAI and Cd (components
// 0, 1, 2) on the ground slab of f's boxes. One kernel per profile.
template <typename Fraction>
void fill_map_drag (amrex::MultiFab& f, const amrex::MultiFab& cols, const amrex::MultiFab& z_nd, int k0,
                    const Fraction& frac)
{
    for (amrex::MFIter mfi(f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto out = f.array(mfi);
        const auto col = cols.const_array(mfi);
        const int kc = cols.box(mfi.index()).smallEnd(2);
        const CellHeights ch{z_nd.const_array(mfi), k0};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            double zb = 0, zt = 0;
            ch(i, j, k, zb, zt);
            out(i, j, k) = static_cast<amrex::Real>(cell_drag(zb, zt, static_cast<double>(col(i, j, kc, 0)),
                                                              static_cast<double>(col(i, j, kc, 1)),
                                                              static_cast<double>(col(i, j, kc, 2)), frac));
        });
    }
}

// A patch on the device: its centre, radius squared, height, Cd, LAI, and
// its profile (lm.c = 0 marks uniform; type 2 always has c > 0).
struct DevicePatch
{
    double x, y, r2, h, cd, lai;
    int lm;
    LmProfile prof;
};

// The patches' drag (they add where they overlap), with the nearest periodic
// image in a periodic direction (period lx, ly; 0: not periodic). A patch
// wider than the period would reach a column through two images: refused.
// Cells whose lower face is above every patch's top skip the patch loop.
void fill_patch_drag (amrex::MultiFab& f, const amrex::Geometry& geom, const amrex::MultiFab& z_nd, int k0,
                      const std::vector<CanopyPatch>& patches)
{
    double hmax = 0;
    for (std::size_t n = 0; n < patches.size(); ++n) {
        const auto& p = patches[n];
        hmax = std::max(hmax, p.height);
        for (int d = 0; d < 2; ++d) {
            if (geom.isPeriodic(d) && p.diameter > static_cast<double>(geom.ProbLength(d))) {
                amrex::Abort("gust: canopy.file: patch " + std::to_string(n + 1) + " (x " + std::to_string(p.x) + ", y "
                             + std::to_string(p.y) + ") has a diameter of " + std::to_string(p.diameter)
                             + " m, wider than the periodic domain in " + (d == 0 ? "x" : "y") + " ("
                             + std::to_string(static_cast<double>(geom.ProbLength(d))) + " m)");
            }
        }
    }
    amrex::Gpu::DeviceVector<DevicePatch> dp(patches.size());
    {
        std::vector<DevicePatch> hp;
        hp.reserve(patches.size());
        for (const auto& p : patches) {
            const double r = 0.5 * p.diameter;
            hp.push_back({p.x, p.y, r * r, p.height, p.cd, p.lai, p.type == 2 ? 1 : 0,
                          p.type == 2 ? make_lm_profile(p.laimax) : LmProfile{}});
        }
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, hp.begin(), hp.end(), dp.begin());
    }
    const DevicePatch* pd = dp.data();
    const int np = static_cast<int>(patches.size());
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const double lx = geom.isPeriodic(0) ? static_cast<double>(geom.ProbLength(0)) : 0.0;
    const double ly = geom.isPeriodic(1) ? static_cast<double>(geom.ProbLength(1)) : 0.0;
    for (amrex::MFIter mfi(f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto out = f.array(mfi);
        const CellHeights ch{z_nd.const_array(mfi), k0};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            double zb = 0, zt = 0;
            ch(i, j, k, zb, zt);
            if (zb >= hmax) {
                out(i, j, k) = 0;
                return;
            }
            const double x = static_cast<double>(plo[0]) + (i + 0.5) * static_cast<double>(dx[0]);
            const double y = static_cast<double>(plo[1]) + (j + 0.5) * static_cast<double>(dx[1]);
            double sum = 0;
            for (int n = 0; n < np; ++n) {
                const DevicePatch& p = pd[n];
                double ex = x - p.x, ey = y - p.y;
                if (lx > 0) { ex -= lx * std::round(ex / lx); }
                if (ly > 0) { ey -= ly * std::round(ey / ly); }
                if (ex * ex + ey * ey > p.r2) { continue; }
                sum += p.lm != 0 ? cell_drag(zb, zt, p.h, p.lai, p.cd, LmFraction{p.prof})
                                 : cell_drag(zb, zt, p.h, p.lai, p.cd, UniformFraction{});
            }
            out(i, j, k) = static_cast<amrex::Real>(sum);
        });
    }
    // dp is freed at the end of this scope.
    amrex::Gpu::streamSynchronize();
}

} // namespace

CanopyOptions read_canopy_options (int idw_neighbors, double idw_power)
{
    CanopyOptions opt;
    opt.idw_neighbors = idw_neighbors;
    opt.idw_power = idw_power;
    amrex::ParmParse pp("canopy");
    for (const char* key : kCanopyKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) {
            opt.parse_errors.emplace_back(std::string("canopy.") + key + " is given more than once");
        } else if (pp.countval(key) != 1) {
            opt.parse_errors.emplace_back(std::string("canopy.") + key + " takes exactly one value");
        }
    }
    pp.query("type", opt.type_name);
    pp.query("file", opt.file);
    pp.query("height_file", opt.height_file);
    pp.query("lai_file", opt.lai_file);
    pp.query("cd_file", opt.cd_file);
    pp.query("profile", opt.profile_name);
    // The numbers as text first: a value that is not a number is reported with the others.
    auto number = [&] (const char* key, double& dst) {
        std::string s;
        if (pp.query(key, s) && !parse_number(s, dst)) {
            opt.parse_errors.emplace_back(std::string("canopy.") + key + " = " + s + " is not a number");
        }
    };
    number("cd", opt.cd);
    number("laimax", opt.laimax);
    if (opt.type_name == "patches") { opt.type = CanopyType::patches; }
    if (opt.type_name == "map") { opt.type = CanopyType::map; }
    if (opt.profile_name == "lalic_mihailovic") { opt.profile = CanopyProfile::lalic_mihailovic; }
    return opt;
}

std::vector<std::string> check_canopy_options (const CanopyOptions& opt)
{
    std::vector<std::string> err = opt.parse_errors;
    auto given = [&opt] (const char* k) { return has(opt.given, k); };
    auto no_effect = [&] (const char* k, const std::string& why) {
        if (given(k)) { err.emplace_back(std::string("canopy.") + k + " has no effect " + why); }
    };
    auto need_file = [&] (const char* k, const std::string& f) {
        if (f.empty()) {
            err.emplace_back(std::string("canopy.") + k + " is required with canopy.type = " + opt.type_name);
        } else if (!readable(f)) {
            err.emplace_back(std::string("canopy.") + k + " = " + f + ": cannot be read");
        }
    };
    if (opt.type_name != "none" && opt.type_name != "patches" && opt.type_name != "map") {
        err.emplace_back("canopy.type must be none, patches or map (got " + opt.type_name + ")");
        return err;
    }
    if (opt.type == CanopyType::none) {
        for (const char* k : kCanopyKeys) {
            if (std::string(k) != "type") { no_effect(k, "with canopy.type = none"); }
        }
        return err;
    }
    if (opt.type == CanopyType::patches) {
        need_file("file", opt.file);
        for (const char* k : {"height_file", "lai_file", "cd", "cd_file", "profile", "laimax"}) {
            no_effect(k, "with canopy.type = patches (each line of canopy.file gives its patch's height, Cd, LAI, "
                         "profile and laimax)");
        }
        return err;
    }
    // map
    no_effect("file", "with canopy.type = map (canopy.file is the patch file)");
    need_file("height_file", opt.height_file);
    need_file("lai_file", opt.lai_file);
    if (given("cd") && given("cd_file")) {
        err.emplace_back("canopy.cd and canopy.cd_file cannot both be given (one drag coefficient: a value or a map)");
    } else if (given("cd_file")) {
        need_file("cd_file", opt.cd_file);
    } else if (!given("cd")) {
        err.emplace_back("canopy.type = map needs canopy.cd or canopy.cd_file (the drag coefficient has no default, as in ERF)");
    } else if (!std::isfinite(opt.cd) || !(opt.cd > 0)) {
        err.emplace_back("canopy.cd must be finite and > 0");
    }
    if (opt.profile_name != "uniform" && opt.profile_name != "lalic_mihailovic") {
        err.emplace_back("canopy.profile must be uniform or lalic_mihailovic (got " + opt.profile_name + ")");
    } else if (opt.profile == CanopyProfile::lalic_mihailovic) {
        const std::string e = erf_forest_utils::validate_laimax(static_cast<amrex::Real>(opt.laimax), "canopy.laimax");
        if (!e.empty()) { err.push_back(e); }
    } else {
        no_effect("laimax", "with canopy.profile = uniform");
    }
    return err;
}

std::vector<CanopyPatch> read_canopy_patches (const std::string& file, const std::string& key)
{
    std::vector<double> flat;   // eight values per patch
    std::string err;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(file);
        if (!is) { err = "cannot be read"; }
        int line_no = 0;
        int shown = 0;
        for (std::string line; is && std::getline(is, line);) {
            ++line_no;
            line = line.substr(0, line.find('#'));
            std::istringstream ls(line);
            std::vector<std::string> tok;
            for (std::string t; ls >> t;) { tok.push_back(t); }
            if (tok.empty()) { continue; }
            std::string bad;
            double v[8] = {};
            if (tok.size() != 8) {
                bad = std::to_string(tok.size()) + " values, not the 8 of 'type x y height diameter cd lai laimax'";
            } else {
                for (int n = 0; n < 8 && bad.empty(); ++n) {
                    if (!parse_number(tok[n], v[n]) || !std::isfinite(v[n])) { bad = "'" + tok[n] + "' is not a finite number"; }
                }
            }
            if (bad.empty()) {
                if (v[0] != 1 && v[0] != 2) {
                    bad = "type must be 1 (uniform) or 2 (Lalic-Mihailovic)";
                } else if (!(v[3] > 0)) {
                    bad = "height must be > 0";
                } else if (!(v[4] > 0)) {
                    bad = "diameter must be > 0";
                } else if (!(v[5] >= 0)) {
                    bad = "cd must be >= 0";
                } else if (!(v[6] >= 0)) {
                    bad = "lai must be >= 0";
                } else if (v[0] == 2) {
                    bad = erf_forest_utils::validate_laimax(static_cast<amrex::Real>(v[7]), "laimax");
                }
            }
            if (!bad.empty()) {
                if (shown < 10) { err += "\n  line " + std::to_string(line_no) + ": " + bad; }
                ++shown;
                continue;
            }
            flat.insert(flat.end(), v, v + 8);
        }
        if (err.empty() && flat.empty()) { err = "holds no patch"; }
    }
    bcast(err);
    if (!err.empty()) { amrex::Abort("gust: " + key + " = " + file + ": " + err); }
    long n = static_cast<long>(flat.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    flat.resize(static_cast<std::size_t>(n));
    amrex::ParallelDescriptor::Bcast(flat.data(), static_cast<std::size_t>(n), root);
    std::vector<CanopyPatch> out;
    for (std::size_t i = 0; i + 8 <= flat.size(); i += 8) {
        out.push_back({static_cast<int>(flat[i]), flat[i + 1], flat[i + 2], flat[i + 3], flat[i + 4], flat[i + 5],
                       flat[i + 6], flat[i + 7]});
    }
    return out;
}

std::string Canopy::hash_inputs (const CanopyOptions& opt)
{
    std::string s = std::string(kCanopyMagic) + " " + std::to_string(kCanopyVersion) + "\ntype " + opt.type_name + "\n";
    if (opt.type == CanopyType::patches) {
        s += "file " + file_hash(opt.file, "canopy.file") + "\n";
    } else if (opt.type == CanopyType::map) {
        s += "height_file " + file_hash(opt.height_file, "canopy.height_file") + "\n";
        s += "lai_file " + file_hash(opt.lai_file, "canopy.lai_file") + "\n";
        s += opt.cd_file.empty() ? "cd " + hexfloat(opt.cd) + "\n" : "cd_file " + file_hash(opt.cd_file, "canopy.cd_file") + "\n";
        s += "profile " + opt.profile_name + "\n";
        if (opt.profile == CanopyProfile::lalic_mihailovic) { s += "laimax " + hexfloat(opt.laimax) + "\n"; }
        s += "idw " + std::to_string(opt.idw_neighbors) + " " + hexfloat(opt.idw_power) + "\n";
    }
    return sha256_hex(s);
}

void Canopy::build (const CanopyOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                    const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm)
{
    if (!opt.active()) { amrex::Abort("gust: Canopy::build needs canopy.type = patches or map"); }
    m_active = true;
    m_hash = hash_inputs(opt);
    m_drag.clear();
    m_drag.define(cell_ba, dm, 1, amrex::IntVect(1, 1, 0));
    m_drag.setVal(0);
    const int k0 = geom.Domain().smallEnd(2);
    if (opt.type == CanopyType::patches) {
        fill_patch_drag(m_drag, geom, z_nd, k0, read_canopy_patches(opt.file, "canopy.file"));
    } else {
        // Each column's height, LAI and Cd, by IDW on the host (as the
        // roughness map), on the ground slab of the cells.
        const TerrainPoints hp = read_map_points(opt.height_file, "canopy.height_file", opt.idw_neighbors);
        const TerrainPoints lp = read_map_points(opt.lai_file, "canopy.lai_file", opt.idw_neighbors);
        TerrainPoints cp;
        if (!opt.cd_file.empty()) { cp = read_map_points(opt.cd_file, "canopy.cd_file", opt.idw_neighbors); }
        const IdwInterpolator ih(hp, opt.idw_neighbors, opt.idw_power);
        const IdwInterpolator il(lp, opt.idw_neighbors, opt.idw_power);
        std::unique_ptr<IdwInterpolator> ic;
        if (!opt.cd_file.empty()) { ic = std::make_unique<IdwInterpolator>(cp, opt.idw_neighbors, opt.idw_power); }
        amrex::MultiFab cols(ground_slab(cell_ba), dm, 3, 0);
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        for (amrex::MFIter mfi(cols, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            amrex::FArrayBox h(b, 3, amrex::The_Pinned_Arena());
            const auto ha = h.array();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const double x = static_cast<double>(plo[0]) + (i + 0.5) * static_cast<double>(dx[0]);
                const double y = static_cast<double>(plo[1]) + (j + 0.5) * static_cast<double>(dx[1]);
                ha(i, j, k, 0) = static_cast<amrex::Real>(ih(x, y));
                ha(i, j, k, 1) = static_cast<amrex::Real>(il(x, y));
                ha(i, j, k, 2) = static_cast<amrex::Real>(ic ? (*ic)(x, y) : opt.cd);
            });
            cols[mfi].copy<amrex::RunOn::Device>(h, b, 0, b, 0, 3);
            // h is freed at the end of this iteration: the copy must be done.
            amrex::Gpu::streamSynchronize();
        }
        if (opt.profile == CanopyProfile::uniform) {
            fill_map_drag(m_drag, cols, z_nd, k0, UniformFraction{});
        } else {
            fill_map_drag(m_drag, cols, z_nd, k0, LmFraction{make_lm_profile(opt.laimax)});
        }
    }
    image_ghost_columns(m_drag, geom);
}

void Canopy::restrict_from (const Canopy& fine, const amrex::MultiFab& fine_detJ, const amrex::Geometry& geom,
                            const Mesh& mesh, const amrex::BoxArray& cell_ba, const amrex::DistributionMapping& dm)
{
    if (!fine.active()) { amrex::Abort("gust: Canopy::restrict_from: the fine level has no canopy"); }
    const amrex::MultiFab& ff = fine.m_drag;
    if (fine_detJ.boxArray() != ff.boxArray() || fine_detJ.DistributionMap() != ff.DistributionMap()
        || mesh.detJ.boxArray() != cell_ba || mesh.detJ.DistributionMap() != dm) {
        amrex::Abort("gust: Canopy::restrict_from: detJ is not on the canopy's layout");
    }
    const amrex::Box fdom = ff.boxArray().minimalBox();
    const amrex::Box& cdom = geom.Domain();
    const int rx = fdom.length(0) / cdom.length(0), ry = fdom.length(1) / cdom.length(1);
    if (rx < 1 || rx != ry || fdom.length(0) != rx * cdom.length(0) || fdom.length(1) != ry * cdom.length(1)
        || fdom.length(2) != cdom.length(2) || !ff.boxArray().coarsenable(amrex::IntVect(rx, ry, 1))) {
        amrex::Abort("gust: Canopy::restrict_from: the level does not coarsen level 0 by one ratio in x and y");
    }
    m_active = true;
    m_hash = fine.m_hash;
    // sum(f detJ) / r^2 over each coarse cell's fine cells, on the fine boxes
    // coarsened (local to each box), then copied to the level's layout (which
    // may be gathered on one box) and divided by the coarse detJ.
    const amrex::IntVect ratio(rx, ry, 1);
    amrex::MultiFab sum(amrex::coarsen(ff.boxArray(), ratio), ff.DistributionMap(), 1, 0);
    const amrex::Real inv = amrex::Real(1) / static_cast<amrex::Real>(rx * ry);
    for (amrex::MFIter mfi(sum, false); mfi.isValid(); ++mfi) {
        const auto s = sum.array(mfi);
        const auto f = ff.const_array(mfi);
        const auto J = fine_detJ.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real a = 0;
            for (int b = 0; b < ry; ++b) {
                for (int c = 0; c < rx; ++c) { a += f(rx * i + c, ry * j + b, k) * J(rx * i + c, ry * j + b, k); }
            }
            s(i, j, k) = a * inv;
        });
    }
    m_drag.clear();
    m_drag.define(cell_ba, dm, 1, amrex::IntVect(1, 1, 0));
    m_drag.setVal(0);
    m_drag.ParallelCopy(sum, 0, 0, 1);
    for (amrex::MFIter mfi(m_drag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto d = m_drag.array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i, j, k) /= J(i, j, k);
        });
    }
    image_ghost_columns(m_drag, geom);
}

std::string Canopy::summary (const amrex::MultiFab& detJ, const amrex::Geometry& geom) const
{
    if (!m_active) { return "none"; }
    // Per column: the largest f, the integral of f (sum of f detJ dzeta),
    // canopy at the ground, canopy in the top cell.
    const amrex::Box& dom = geom.Domain();
    const int k0 = dom.smallEnd(2), k1 = dom.bigEnd(2);
    const amrex::Real dz = geom.CellSize(2);
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum> ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Long, amrex::Long> data(ops);
    for (amrex::MFIter mfi(m_drag, false); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.validbox();
        b.setBig(2, b.smallEnd(2));
        const auto f = m_drag.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        ops.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Long, amrex::Long> {
            amrex::Real fmax = 0, integral = 0;
            for (int k = k0; k <= k1; ++k) {
                fmax = amrex::max(fmax, f(i, j, k));
                integral += f(i, j, k) * J(i, j, k) * dz;
            }
            const bool on = f(i, j, k0) > 0;
            return {fmax, on ? integral : amrex::Real(0), integral, on ? amrex::Long(1) : amrex::Long(0),
                    f(i, j, k1) > 0 ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    const auto r = data.value();
    amrex::Real fmax = amrex::get<0>(r), total = amrex::get<1>(r), imax = amrex::get<2>(r);
    amrex::Long ncol = amrex::get<3>(r), ntop = amrex::get<4>(r);
    const auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Max(fmax, comm);
    amrex::ParallelAllReduce::Sum(total, comm);
    amrex::ParallelAllReduce::Max(imax, comm);
    amrex::ParallelAllReduce::Sum(ncol, comm);
    amrex::ParallelAllReduce::Sum(ntop, comm);
    const long ncols = static_cast<long>(dom.length(0)) * dom.length(1);
    char b[256];
    std::snprintf(b, sizeof(b), "largest f %.4g 1/m; canopy on %ld of %ld columns; Cd LAI over them: mean %.4g, largest %.4g",
                  static_cast<double>(fmax), static_cast<long>(ncol), ncols,
                  ncol > 0 ? static_cast<double>(total) / static_cast<double>(ncol) : 0.0, static_cast<double>(imax));
    std::string out = b;
    if (ncol == 0) {
        out += "\nWARNING: the canopy covers no column of this level: check that its x and y are in the mesh's frame";
    }
    if (ntop > 0) {
        out += "\nWARNING: the canopy reaches the top cell in " + std::to_string(ntop)
               + " columns: the part of it above the domain is left out";
    }
    return out;
}

void Canopy::write_checkpoint (const std::string& dir) const
{
    if (m_active && amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = canopy_checkpoint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << kCanopyMagic << " " << kCanopyVersion << "\nsha256 " << m_hash << "\n";
    }
    amrex::ParallelDescriptor::Barrier();
}

void Canopy::check_restart (const std::string& dir, const std::string& hash)
{
    const bool active = !hash.empty();
    // The checkpoint's hash ("" without a canopy file; "?" when unreadable).
    std::string rec;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(canopy_checkpoint_file(dir));
        if (is) {
            std::string magic, key;
            int version = 0;
            rec = "?";
            if (is >> magic >> version >> key >> rec && magic == kCanopyMagic && version == kCanopyVersion && key == "sha256") {
                // rec holds the hash
            } else {
                rec = "?";
            }
        }
    }
    bcast(rec);
    const std::string where = "gust: io.restart = " + dir + ": ";
    if (rec == "?") { amrex::Abort(where + "its Canopy file cannot be read"); }
    if (active && rec.empty()) {
        amrex::Abort(where + "the checkpoint was written without a canopy; a restart cannot add one (canopy.type = none, or "
                             "start a new run)");
    }
    if (!active && !rec.empty()) {
        amrex::Abort(where + "the checkpoint was written with a canopy; give the run's canopy.* inputs");
    }
    if (active && rec != hash) {
        amrex::Abort(where + "the canopy inputs differ from the checkpoint's (SHA-256 " + hash + " here, " + rec
                     + " there); a restart continues the same canopy");
    }
}

} // namespace gust
