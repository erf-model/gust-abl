#include "ImmersedTerrain.H"

#include "Sha256.H"
#include "TerrainSurface.H"

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <fstream>

namespace gust {

namespace {

constexpr const char* kMagic = "gust_immersed_terrain";
constexpr int kVersion = 3;   // 1: a flat ground (flat_height, small_volfrac, z0); 2: roofs with the
                              // buildings' tolerance; 3: roofs with none (ground_roof)

void bcast_string (std::string& s)
{
    long n = static_cast<long>(s.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

}   // namespace

std::vector<std::string> check_immersed_terrain (const ImmersedTerrainFacts& f)
{
    std::vector<std::string> err;
    const std::string rep = "terrain.representation = " + f.representation;
    if (!f.surface_layer) { err.push_back(rep + " needs bc.zlo = surface_layer (the wall law acts at the immersed ground)"); }
    if (!f.neutral) { err.push_back(rep + " needs bc.surface_stability = neutral for now"); }
    if (f.buildings) { err.push_back(rep + " cannot be combined with buildings.type for now"); }
    if (f.nests) { err.push_back(rep + " cannot be combined with nests for now"); }
    if (f.newton) { err.push_back(rep + " is not supported with steady.stage = simplec_newton yet"); }
    if (!f.zline_sweeps) { err.push_back(rep + " needs z-line sweeps (steady.sweep = zjacobi) for now"); }
    if (f.poisson_mlmg) {
        err.push_back(rep + " cannot use solver.poisson_backend = mlmg (no solid cells in its operator); use gmres_mlmg or gmres_fft");
    }
    if (f.unsteady) { err.push_back(rep + " is not supported with steady.stage = unsteady yet"); }
    if (f.canopy) { err.push_back(rep + " cannot be combined with canopy.type for now (the canopy's heights are measured from the mesh's lower face)"); }
    if (f.stratified) { err.push_back(rep + " needs a uniform bc.theta_profile for now (theta has no terms at the immersed ground)"); }
    if (f.column_guess) { err.push_back(rep + " is not supported with init.first_guess = column yet (the column's ground is the mesh's lower face); use profile"); }
    if (f.probes) { err.push_back(rep + " is not supported with steady.probes_file yet (the probes' heights are measured from the mesh's lower face)"); }
    if (f.zref_given) { err.push_back("bc.surface_zref has no effect with " + rep + " (the wall law reads the cut cell and the cell above it)"); }
    if (f.floor_given) { err.push_back("bc.surface_denominator_floor has no effect with " + rep + " (the wall law has no floor)"); }
    if (f.wall_k_source) { err.push_back(rep + " needs bc.surface_wall_k = dirichlet (the cut cell holds the wall value of k)"); }
    if (f.outflow_from_velocity) {
        err.push_back(rep + " needs bc.outflow_transport = mirror (from_velocity, kept to continue runs of v1.0.0, slows the"
                            " cut cell at an outflow face)");
    }
    if (f.laminar) { err.push_back(rep + " needs turb.closure = kl (the immersed ground is a wall law)"); }
    return err;
}

std::string check_immersed_terrain_height (const TerrainInputs& t, double first_layer_top, double top_layer_bottom)
{
    const double h = static_cast<double>(t.flat_height);
    if (t.representation == "hybrid") {
        // The lower face may be the ground (the columns where the envelope meets the terrain).
        if (h >= 0 && h <= top_layer_bottom) { return {}; }
        return "terrain.flat_height must lie between 0 and the bottom of the mesh's top layer (" + std::to_string(top_layer_bottom)
               + " m) with terrain.representation = hybrid (got " + std::to_string(h) + ")";
    }
    if (h >= first_layer_top && h <= top_layer_bottom) { return {}; }
    return "terrain.flat_height must lie between the top of the mesh's first layer (" + std::to_string(first_layer_top)
           + " m) and the bottom of its top layer (" + std::to_string(top_layer_bottom)
           + " m) with terrain.representation = " + t.representation + " (got " + std::to_string(h) + ")";
}

BuildingsOptions immersed_terrain_options (const TerrainInputs& t, double z0, bool z0_map)
{
    if (!t.ground_immersed()) {
        amrex::Abort("gust: immersed_terrain_options needs terrain.representation = immersed_forcing or hybrid");
    }
    BuildingsOptions o;
    o.type = BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.file = "(terrain.flat_height)";
    o.small_volfrac = static_cast<double>(t.small_volfrac);
    o.z0 = z0;
    o.ground_z0_map = z0_map;
    o.wall_form = BuildingWallForm::fraction_stress;
    o.wall_form_name = "fraction_stress";
    o.wall_speed = default_wall_speed(o.wall_form);
    o.wall_tke = true;
    o.theta_mode = BuildingThetaMode::none;
    o.limiter_zone = 0;
    o.limiter_zone_acts = false;
    o.immersed_ground = true;
    o.hybrid_mesh = t.representation == "hybrid";
    return o;
}

void fill_immersed_ground (amrex::MultiFab& z_nd, const amrex::Geometry& geom, const TerrainInputs& t)
{
    if (z_nd.nComp() < 2) { amrex::Abort("gust: fill_immersed_ground needs node heights with a ground component"); }
    // The terrain's height in the ground component's k0 slab, ghost columns included (each ghost node
    // takes the height of the node it stands for, and a node shared by two boxes is one expression of
    // its index: no sync is needed); then the base level, and the same height at every k of the column.
    amrex::MultiFab g(z_nd, amrex::make_alias, z_nd.nComp() - 1, 1);
    fill_terrain_surface(g, geom, t);
    const int k0 = geom.Domain().smallEnd(2);
    const amrex::Real base = t.flat_height;
    for (amrex::MFIter mfi(g, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.fabbox();
        if (bx.smallEnd(2) > k0 || bx.bigEnd(2) < k0) {
            amrex::Abort("gust: fill_immersed_ground: every box must span the column (the mesh's layout does)");
        }
        const auto a = g.array(mfi);
        amrex::Box slab = bx;
        slab.setSmall(2, k0);
        slab.setBig(2, k0);
        amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) += base; });
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k != k0) { a(i, j, k) = a(i, j, k0); }
        });
    }
    amrex::Gpu::streamSynchronize();
}

std::string check_immersed_ground_extent (const amrex::MultiFab& z_nd, double first_layer_top, double top_layer_bottom)
{
    const int gc = z_nd.nComp() - 1;
    const double lo = static_cast<double>(z_nd.min(gc)), hi = static_cast<double>(z_nd.max(gc));
    if (lo >= first_layer_top && hi <= top_layer_bottom) { return {}; }
    return "the immersed ground (terrain.flat_height plus the terrain's height: " + std::to_string(lo) + " m to "
           + std::to_string(hi) + " m) must lie between the top of the mesh's first layer (" + std::to_string(first_layer_top)
           + " m) and the bottom of its top layer (" + std::to_string(top_layer_bottom)
           + " m) with terrain.representation = immersed_forcing; change terrain.flat_height or the mesh's levels";
}

BuildingHeightMap immersed_terrain_map (const TerrainInputs& t, const amrex::Geometry& geom)
{
    BuildingHeightMap m;
    m.nx = 1;
    m.ny = 1;
    for (int d = 0; d < 2; ++d) {
        const double lo = static_cast<double>(geom.ProbLo(d)), hi = static_cast<double>(geom.ProbHi(d));
        std::vector<double>& e = d == 0 ? m.xe : m.ye;
        e = {lo, hi};
    }
    m.h.assign(1, static_cast<double>(t.flat_height));
    m.zs.assign(1, 0.0);
    return m;
}

std::string immersed_terrain_hash (const TerrainInputs& t, double z0, const std::string& z0_file)
{
    auto hx = [] (amrex::Real v) { return hexfloat_canonical(static_cast<double>(v)); };
    // "slope vertical": the wall-law form of the immersed ground, now the only one. The text stays in the
    // hashed string so that the checkpoints written while it was an input keep their hash.
    // The roughness: the value, or with a map the file's contents and what its interpolation reads (z0 is then
    // not read). A run without a map keeps the text, and the hash, it had before maps were supported.
    const std::string map_sha = z0_file.empty() ? std::string() : sha256_file(z0_file);
    if (!z0_file.empty() && map_sha.empty()) {
        // A run with no stage does not read the map; its checkpoint would hold the hash of no file.
        amrex::Abort("gust: bc.surface_z0_file = " + z0_file + " cannot be read (the immersed ground's restart record needs it)");
    }
    const std::string rough = z0_file.empty() ? hexfloat_canonical(z0)
                                              : "file " + map_sha + " idw " + std::to_string(t.idw_neighbors) + " " + hx(t.idw_power);
    std::string s = std::string(kMagic) + " " + std::to_string(kVersion) + "\nflat_height " + hx(t.flat_height)
                    + "\nsmall_volfrac " + hx(t.small_volfrac) + "\nz0 " + rough + "\nslope "
                    + "vertical" + "\ntype " + t.type + "\n";
    // Only for a hybrid mesh: an immersed_forcing checkpoint keeps its hash.
    if (t.representation == "hybrid") { s += "hybrid " + hx(t.hybrid_half_width) + "\n"; }
    // The surface: what fill_terrain_surface reads for the type.
    if (t.type == "agnesi") {
        s += "hill " + hx(t.hill_height) + " " + hx(t.hill_half_width) + " " + t.hill_shape + "\n";
    } else if (t.type == "xyz") {
        s += "file " + sha256_file(t.file) + "\nidw " + std::to_string(t.idw_neighbors) + " " + hx(t.idw_power) + "\n";
    } else if (t.type == "erf_grid") {
        s += "file " + sha256_file(t.file) + "\n";
    }
    return sha256_hex(s);
}

void write_immersed_terrain_checkpoint (const std::string& dir, const std::string& hash)
{
    if (!hash.empty() && amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = immersed_terrain_checkpoint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << kMagic << " " << kVersion << "\nsha256 " << hash << "\n";
    }
    amrex::ParallelDescriptor::Barrier();
}

void check_immersed_terrain_restart (const std::string& dir, const std::string& hash)
{
    const bool active = !hash.empty();
    // The checkpoint's hash ("" without the file; "?" when unreadable).
    std::string rec;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(immersed_terrain_checkpoint_file(dir));
        if (is) {
            std::string magic, key;
            int version = 0;
            if (!(is >> magic >> version >> key >> rec) || magic != kMagic || key != "sha256") {
                rec = "?";
            } else if (version != kVersion) {
                rec = "v" + std::to_string(version);   // not a hash: a hash is 64 hex digits
            }
        }
    }
    bcast_string(rec);
    const std::string where = "gust: io.restart = " + dir + ": ";
    if (rec == "?") { amrex::Abort(where + "its ImmersedTerrain file cannot be read"); }
    if (!rec.empty() && rec[0] == 'v') {
        amrex::Abort(where + "its ImmersedTerrain record is version " + rec.substr(1) + ", this build's is " + std::to_string(kVersion)
                     + " (the immersed ground's geometry changed between them); start the run again");
    }
    if (active && rec.empty()) {
        amrex::Abort(where + "the checkpoint's ground is body-fitted; a restart cannot change terrain.representation");
    }
    if (!active && !rec.empty()) {
        amrex::Abort(where + "the checkpoint's ground is immersed; give its terrain.representation (immersed_forcing or hybrid) and its inputs");
    }
    if (active && rec != hash) {
        amrex::Abort(where + "the immersed ground's inputs (terrain.representation, terrain.type and its keys, "
                             "terrain.flat_height, terrain.small_volfrac, terrain.hybrid_half_width, bc.surface_z0 or bc.surface_z0_file "
                             "and its terrain.idw_neighbors and terrain.idw_power) "
                             "differ from the checkpoint's, or the checkpoint was written with a wall-law form since "
                             "removed (terrain.immersed_slope, terrain.immersed_side_walls or terrain.immersed_wall_distance) "
                             "(SHA-256 " + hash + " here, " + rec + " there)");
    }
}

}   // namespace gust
