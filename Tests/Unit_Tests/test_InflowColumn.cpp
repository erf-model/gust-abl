// bc.inflow_profile = column (InflowColumn.H): the run's own column as the inflow faces' profile. The column kind
// evaluates as the in-memory table of the column guess; the inputs that make no sense with it are refused, each
// naming its key; a target wind is accepted with an inflow face only under the column kind; the checkpoint record
// reads back to the last bit, and a restart without it (or with one it does not expect) is refused; the column
// starts from the analytic defaults whatever the deck's inflow options hold.

#include "BoundaryConditions.H"
#include "ColumnGuess.H"
#include "InflowColumn.H"
#include "InflowProfile.H"
#include "Initializer.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using gust::BoundaryInputs;
using gust::FaceType;
using gust::InflowColumnRecord;
using gust::InitOptions;
using K = gust::InflowProfileOptions::Kind;

// The message of the abort f raises ("" if none).
std::string abort_message (const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
}

bool mentions (const std::vector<std::string>& err, const std::string& what)
{
    for (const auto& e : err) { if (e.find(what) != std::string::npos) { return true; } }
    return false;
}

// A column table: values with no short decimal form, so a lossy record would show.
gust::InflowTable sample_table ()
{
    gust::InflowTable t;
    for (int i = 0; i < 6; ++i) {
        const amrex::Real z = (i == 0) ? amrex::Real(0) : amrex::Real(10) * static_cast<amrex::Real>(i * i) / amrex::Real(3);
        t.z.push_back(z);
        t.u.push_back(amrex::Real(8) * std::log1p(z) / amrex::Real(7));
        t.v.push_back(amrex::Real(-1) / amrex::Real(3) + amrex::Real(0.1) * static_cast<amrex::Real>(i));
        t.theta.push_back(amrex::Real(300) + amrex::Real(0.003) * z);
        t.k.push_back(amrex::Real(0.9) / (amrex::Real(1) + amrex::Real(0.01) * z) + amrex::Real(1e-7));
    }
    t.ncol = 5;
    return t;
}

InflowColumnRecord sample_record ()
{
    InflowColumnRecord r;
    r.table = sample_table();
    r.z0 = amrex::Real(0.398);
    r.target = true;
    r.ug = amrex::Real(15.387417118);
    r.vg = amrex::Real(6.4327676801) / amrex::Real(3);
    r.pressure_gradient_x = amrex::Real(1) / amrex::Real(7000);
    r.pressure_gradient_y = amrex::Real(-2.5e-5);
    r.z0_map = true;
    r.coriolis = true;
    r.latitude = amrex::Real(39.7) / amrex::Real(3);
    r.rotation_period = amrex::Real(86164.1);
    return r;
}

void expect_same_record (const InflowColumnRecord& a, const InflowColumnRecord& b)
{
    EXPECT_EQ(a.z0, b.z0);
    EXPECT_EQ(a.target, b.target);
    EXPECT_EQ(a.ug, b.ug);
    EXPECT_EQ(a.vg, b.vg);
    EXPECT_EQ(a.pressure_gradient_x, b.pressure_gradient_x);
    EXPECT_EQ(a.pressure_gradient_y, b.pressure_gradient_y);
    EXPECT_EQ(a.z0_map, b.z0_map);
    EXPECT_EQ(a.coriolis, b.coriolis);
    EXPECT_EQ(a.latitude, b.latitude);
    EXPECT_EQ(a.rotation_period, b.rotation_period);
    EXPECT_EQ(b.table.ncol, 5);
    ASSERT_EQ(a.table.z.size(), b.table.z.size());
    for (amrex::Long i = 0; i < a.table.z.size(); ++i) {
        EXPECT_EQ(a.table.z[i], b.table.z[i]) << i;
        EXPECT_EQ(a.table.u[i], b.table.u[i]) << i;
        EXPECT_EQ(a.table.v[i], b.table.v[i]) << i;
        EXPECT_EQ(a.table.theta[i], b.table.theta[i]) << i;
        EXPECT_EQ(a.table.k[i], b.table.k[i]) << i;
    }
}

// An open-face deck with the column kind that check_inflow_column accepts.
struct Deck
{
    BoundaryInputs bin;
    InitOptions iopt;
    Deck ()
    {
        bin.face = {FaceType::inflow, FaceType::inflow, FaceType::surface_layer, FaceType::outflow, FaceType::outflow,
                    FaceType::slip_wall};
        bin.inflow.kind = K::column;
        iopt.first_guess = InitOptions::FirstGuess::column;
    }
};

} // namespace

// The column kind is the table it holds: the same values as the column guess's in-memory table at every height
// (between rows, beyond the top, at the ground), k from the table; an empty table aborts naming the kind.
TEST(InflowColumn, ColumnKindIsTheTable)
{
    const gust::InflowTable t = sample_table();
    gust::ReferenceProfile ref;
    ref.z = {0};
    ref.theta = {300};
    gust::InflowProfileOptions opt;
    opt.kind = K::column;
    opt.table = t;
    opt.z0 = amrex::Real(0.398);
    const gust::InflowProfile col(opt, ref);
    const gust::InflowProfile mem(t.z, t.u, t.v, t.theta, t.k, ref);
    EXPECT_EQ(col.table_size(), 6);
    EXPECT_TRUE(col.has_k());
    const auto a = col.host_view();
    const auto b = mem.host_view();
    EXPECT_EQ(a.kind, 2);
    for (const amrex::Real z : {amrex::Real(0), amrex::Real(1.7), amrex::Real(13.3), amrex::Real(40), amrex::Real(500)}) {
        amrex::Real ua, va, ta, ka, na, ub, vb, tb, kb, nb;
        a.eval(z, ua, va, ta, ka, na);
        b.eval(z, ub, vb, tb, kb, nb);
        EXPECT_EQ(ua, ub) << z;
        EXPECT_EQ(va, vb) << z;
        EXPECT_EQ(ta, tb) << z;
        EXPECT_EQ(ka, kb) << z;
    }
    amrex::Real u, v, th, k, nut;
    a.eval(amrex::Real(10) / amrex::Real(3), u, v, th, k, nut);   // the second row exactly
    EXPECT_EQ(u, t.u[1]);
    EXPECT_EQ(k, t.k[1]);
    gust::InflowProfileOptions empty;
    empty.kind = K::column;
    EXPECT_NE(abort_message([&] { gust::InflowProfile p(empty, ref); }).find("bc.inflow_profile = column: the column's table is not set"),
              std::string::npos);
    gust::InflowProfileOptions bad = opt;
    bad.table.z[0] = 1;   // not from the ground: the in-memory table's check, named after the kind
    EXPECT_NE(abort_message([&] { gust::InflowProfile p(bad, ref); }).find("gust: bc.inflow_profile = column: the first height must be 0"),
              std::string::npos);
    // The column kind needs no bc.inflow_file (the file kind does).
    EXPECT_TRUE(gust::check_inflow_profile_options(opt).empty());
    gust::InflowProfileOptions file;
    file.kind = K::file;
    EXPECT_TRUE(mentions(gust::check_inflow_profile_options(file), "bc.inflow_file must be given"));
}

// check_inflow_column: the column kind needs the column guess (or a restart) and an inflow face; a Rayleigh layer on
// u, v or tke and nests are refused; other kinds are not its concern.
TEST(InflowColumn, Checks)
{
    {
        Deck d;
        EXPECT_TRUE(gust::check_inflow_column(d.bin, d.iopt, false, false).empty());
        d.bin.rayleigh_zdamp = 1000;
        d.bin.rayleigh_fields = {false, false, true, false, false};   // w only
        EXPECT_TRUE(gust::check_inflow_column(d.bin, d.iopt, false, false).empty());
    }
    {
        Deck d;
        d.iopt.first_guess = InitOptions::FirstGuess::profile;
        EXPECT_TRUE(mentions(gust::check_inflow_column(d.bin, d.iopt, false, false),
                             "bc.inflow_profile = column needs init.first_guess = column"));
        EXPECT_TRUE(gust::check_inflow_column(d.bin, d.iopt, true, false).empty());   // a restart: the checkpoint's column
    }
    {
        Deck d;
        d.bin.face[0] = FaceType::periodic;
        d.bin.face[3] = FaceType::periodic;
        d.bin.face[1] = FaceType::slip_wall;
        d.bin.face[4] = FaceType::slip_wall;
        EXPECT_TRUE(mentions(gust::check_inflow_column(d.bin, d.iopt, false, false), "bc.inflow_profile = column needs an inflow face"));
    }
    for (int f : {0, 1, 4}) {   // u, v, tke
        Deck d;
        d.bin.rayleigh_zdamp = 1000;
        d.bin.rayleigh_fields = {false, false, false, false, false};
        d.bin.rayleigh_fields[f] = true;
        EXPECT_TRUE(mentions(gust::check_inflow_column(d.bin, d.iopt, false, false), "with a Rayleigh layer on u, v or tke")) << f;
        d.bin.rayleigh_zdamp = 0;   // no layer: the fields do nothing
        EXPECT_TRUE(gust::check_inflow_column(d.bin, d.iopt, false, false).empty()) << f;
    }
    {
        Deck d;
        EXPECT_TRUE(mentions(gust::check_inflow_column(d.bin, d.iopt, false, true), "is not supported with nest.names"));
        d.bin.inflow.kind = K::log_law;   // another kind: nothing to check here
        d.iopt.first_guess = InitOptions::FirstGuess::none;
        EXPECT_TRUE(gust::check_inflow_column(d.bin, d.iopt, false, true).empty());
    }
}

// The bc keys the column kind does not read are refused, each named: bc.inflow_z0 (the faces take the column's
// roughness), bc.inflow_tke_length (the table has k), bc.inflow_direction (the table has u and v), bc.inflow_file.
TEST(InflowColumn, BoundaryKeysRefused)
{
    const amrex::Box dom(amrex::IntVect(0), amrex::IntVect(7, 7, 7));
    const amrex::RealBox rb({0., 0., 0.}, {800., 800., 800.});
    const amrex::Geometry geom(dom, rb, 0, {0, 0, 0});
    gust::CoriolisOptions cor;
    Deck d;
    d.bin.theta_ref.z = {0};
    d.bin.theta_ref.theta = {300};
    EXPECT_TRUE(gust::check_boundary_inputs(d.bin, geom, cor).empty());
    const char* const keys[][2] = {
        {"inflow_z0", "bc.inflow_z0 has no effect with bc.inflow_profile = column (the faces take the column's roughness"},
        {"inflow_tke_length", "bc.inflow_tke_length has no effect with bc.inflow_profile = column (the column's table has k)"},
        {"inflow_direction", "bc.inflow_direction has no effect with bc.inflow_profile = column"},
        {"inflow_file", "bc.inflow_file has no effect with bc.inflow_profile = column"},
        {"inflow_ustar", "bc.inflow_ustar has no effect with bc.inflow_profile = column"},
    };
    for (const auto& kv : keys) {
        BoundaryInputs b = d.bin;
        b.given = {"inflow_profile", kv[0]};
        EXPECT_TRUE(mentions(gust::check_boundary_inputs(b, geom, cor), kv[1])) << kv[0];
        b.inflow.kind = K::log_law;   // the log law reads bc.inflow_z0, bc.inflow_tke_length and bc.inflow_direction
        const bool read = std::string(kv[0]) != "inflow_file";
        EXPECT_EQ(gust::check_boundary_inputs(b, geom, cor).empty(), read) << kv[0];
    }
    // The reader takes the name.
    amrex::ParmParse pp("bc");
    pp.add("inflow_profile", std::string("column"));
    EXPECT_EQ(gust::read_boundary_inputs().inflow.kind, K::column);
    while (pp.remove("inflow_profile") > 0) {}
    pp.add("inflow_profile", std::string("columns"));
    EXPECT_NE(abort_message([] { (void) gust::read_boundary_inputs(); }).find("must be log_law, power_law, file or column"),
              std::string::npos);
    while (pp.remove("inflow_profile") > 0) {}
}

// The target wind with an inflow face or a sponge: refused with the deck's own profile (the face would keep it),
// accepted with the column kind (the faces and the sponge take the target column).
TEST(InflowColumn, TargetWithInflowFace)
{
    InitOptions o;
    o.first_guess = InitOptions::FirstGuess::column;
    o.column_target_u = 6;
    o.column_target_height = 10;
    gust::SteadyOptions so;
    so.stage = gust::StageKind::simplec;
    gust::ForcingOptions fopt;
    fopt.pressure_gradient = {amrex::Real(1e-4), 0, 0};
    gust::SurfaceLayerOptions cs;
    Deck d;
    d.bin.inflow.kind = K::log_law;
    auto err = gust::check_column_guess(o, so, fopt, d.bin, cs, true);
    EXPECT_TRUE(mentions(err, "with an inflow face: the face keeps the deck's inflow profile; set bc.inflow_profile = column"));
    d.bin.inflow.kind = K::column;
    EXPECT_TRUE(gust::check_column_guess(o, so, fopt, d.bin, cs, true).empty());
    d.bin.sponge_width = 100;
    d.bin.sponge_strength = amrex::Real(0.1);
    EXPECT_TRUE(gust::check_column_guess(o, so, fopt, d.bin, cs, true).empty());
    d.bin.inflow.kind = K::file;
    err = gust::check_column_guess(o, so, fopt, d.bin, cs, true);
    EXPECT_TRUE(mentions(err, "with a sponge (bc.sponge_width)"));
    EXPECT_TRUE(mentions(err, "with an inflow face"));
    // The Rayleigh layer on u and v stays refused with a target under the column kind too.
    d.bin.inflow.kind = K::column;
    d.bin.rayleigh_zdamp = 500;
    EXPECT_TRUE(mentions(gust::check_column_guess(o, so, fopt, d.bin, cs, true), "with a Rayleigh layer on u and v"));
}

// The record as text reads back to the last bit (%.17g), and a bad or truncated record aborts with its source named.
TEST(InflowColumn, RecordRoundTrip)
{
    const InflowColumnRecord r = sample_record();
    const std::string text = gust::format_inflow_column_record(r);
    expect_same_record(r, gust::parse_inflow_column_record(text, "test"));
    InflowColumnRecord none = r;
    none.target = false;
    EXPECT_FALSE(gust::parse_inflow_column_record(gust::format_inflow_column_record(none), "test").target);
    auto aborts = [] (const std::string& t, const std::string& what) {
        const std::string m = abort_message([&] { (void) gust::parse_inflow_column_record(t, "src"); });
        EXPECT_NE(m.find("gust: src: " + what), std::string::npos) << m;
    };
    aborts("GustInflowColumn 2\n", "not a GustInflowColumn record of version 1");
    aborts("Something 1\n", "not a GustInflowColumn record of version 1");
    aborts(text.substr(0, text.rfind("\n", text.size() / 2) + 1), "the record ends before");   // at a line end
    aborts(text + "1 2 3\n", "unexpected text after the last row");
    std::string flag = text;
    flag.replace(flag.find("forcing 1"), 9, "forcing 2");
    aborts(flag, "the forcing flag must be 0 or 1");
    std::string map = text;
    map.replace(map.find(" 1\nforcing"), 2, " 3");
    aborts(map, "the map flag must be 0 or 1");
    std::string cor = text;
    cor.replace(cor.find("coriolis 1"), 10, "coriolis 2");
    aborts(cor, "the Coriolis flag must be 0 or 1");
    InflowColumnRecord flat = r;
    flat.z0_map = false;
    flat.coriolis = false;
    const InflowColumnRecord fr = gust::parse_inflow_column_record(gust::format_inflow_column_record(flat), "test");
    EXPECT_FALSE(fr.z0_map);
    EXPECT_FALSE(fr.coriolis);
    std::string rows = text;
    rows.replace(rows.find("rows 6"), 6, "rows 1");
    aborts(rows, "the row count must be at least 2");
    std::string z0 = text;
    z0.replace(z0.find("z0 "), 3, "z0 -");
    aborts(z0, "z0 must be > 0");
    std::string word = text;
    word.replace(word.find("rows"), 4, "lines");
    aborts(word, "expected 'rows'");
    std::string nan = text;
    nan.replace(nan.find("forcing 1 ") + 10, 1, "nan");
    aborts(nan, "bad ug");
    std::string fflag = text;
    fflag.replace(fflag.find("forcing 1"), 9, "forcing 1.0");
    aborts(fflag, "bad the forcing flag '1.0'");
    std::string frows = text;
    frows.replace(frows.find("rows 6"), 6, "rows 6.5");
    aborts(frows, "bad the row count '6.5'");
    // The rows as the table kind needs them, the record named.
    auto row_aborts = [&] (int row, int col, const char* value, const std::string& what) {
        InflowColumnRecord b = sample_record();
        std::vector<amrex::Real>* cols[] = {&b.table.z, nullptr, nullptr, &b.table.theta, &b.table.k};
        (*cols[col])[row] = static_cast<amrex::Real>(std::strtod(value, nullptr));
        aborts(gust::format_inflow_column_record(b), "row " + std::to_string(row + 1) + ": " + what);
    };
    row_aborts(0, 0, "1", "the first height must be 0");
    row_aborts(3, 0, "1", "the heights must increase");
    row_aborts(2, 3, "0", "theta must be > 0");
    row_aborts(4, 4, "-1e-9", "k must be >= 0");
    if (sizeof(amrex::Real) == 4) {   // finite in double, not in float
        std::string big = text;
        big.replace(big.find("forcing 1 ") + 10, 0, "1e39 ");
        big.erase(big.find("1e39 ") + 5, big.find(' ', big.find("1e39 ") + 5) - big.find("1e39 ") - 5 + 1);
        aborts(big, "bad ug '1e39'");
    }
}

// The record in a checkpoint directory: written by the I/O rank, read on every rank; a restart with the column kind
// needs it, and one without the kind refuses it.
TEST(InflowColumn, CheckpointRecord)
{
    const std::string dir = "ut_inflow_column_chk";
    const std::string bare = "ut_inflow_column_bare";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::UtilCreateDirectory(dir, 0755);
        amrex::UtilCreateDirectory(bare, 0755);
    }
    amrex::ParallelDescriptor::Barrier();
    const InflowColumnRecord r = sample_record();
    gust::write_inflow_column_checkpoint(dir, r);
    EXPECT_EQ(gust::inflow_column_checkpoint_file(dir), dir + "/InflowColumn");
    expect_same_record(r, gust::read_inflow_column_checkpoint(dir, true));
    EXPECT_NE(abort_message([&] { (void) gust::read_inflow_column_checkpoint(dir, false); })
                  .find("io.restart = " + dir + ": the checkpoint was written with bc.inflow_profile = column"),
              std::string::npos);
    EXPECT_NE(abort_message([&] { (void) gust::read_inflow_column_checkpoint(bare, true); })
                  .find("io.restart = " + bare + ": the checkpoint was written without bc.inflow_profile = column"),
              std::string::npos);
    EXPECT_TRUE(gust::read_inflow_column_checkpoint(bare, false).table.z.empty());
    // A restart with no stage: the record if there is one, nothing otherwise.
    expect_same_record(r, gust::read_inflow_column_checkpoint_if_present(dir));
    EXPECT_TRUE(gust::read_inflow_column_checkpoint_if_present(bare).table.z.empty());
    // A record that exists but cannot be opened is a read error, not a checkpoint without the record.
    const std::string locked = "ut_inflow_column_locked";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::UtilCreateDirectory(locked, 0755);
        std::ofstream(locked + "/InflowColumn") << "x";
        chmod((locked + "/InflowColumn").c_str(), 0);
    }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor() && !std::ifstream(locked + "/InflowColumn")) {   // not as root
        EXPECT_NE(abort_message([&] { (void) gust::read_inflow_column_checkpoint_if_present(locked); })
                      .find("InflowColumn exists but could not be opened"),
                  std::string::npos);
    }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        chmod((locked + "/InflowColumn").c_str(), 0644);
        amrex::FileSystem::RemoveAll(locked);
    }
    // An empty record is not written.
    EXPECT_NE(abort_message([&] { gust::write_inflow_column_checkpoint(bare, InflowColumnRecord{}); })
                  .find("the inflow column record is not set"),
              std::string::npos);
    // use_inflow_column: the column kind with the record's table and z0, the other options kept.
    gust::InflowProfileOptions opt;
    opt.l_max = 77;
    gust::use_inflow_column(r, opt);
    EXPECT_EQ(opt.kind, K::column);
    EXPECT_EQ(opt.z0, r.z0);
    EXPECT_EQ(opt.l_max, 77);
    EXPECT_EQ(opt.table.k[3], r.table.k[3]);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::FileSystem::RemoveAll(dir);
        amrex::FileSystem::RemoveAll(bare);
    }
    amrex::ParallelDescriptor::Barrier();
}

// The faces' roughness is the column's: init.column_z0 under a map, else bc.surface_z0.
TEST(InflowColumn, Roughness)
{
    InitOptions io;
    io.column_z0 = amrex::Real(0.398);
    gust::SurfaceLayerOptions s;
    s.z0 = amrex::Real(0.05);
    EXPECT_EQ(gust::column_roughness(io, s), amrex::Real(0.05));
    s.z0_file = "z0.xyz";
    EXPECT_EQ(gust::column_roughness(io, s), amrex::Real(0.398));
}

// The column starts from the analytic defaults under the column kind (as a periodic column deck without bc.inflow_*
// keys does), whatever the deck's inflow options hold: the converged column is the same to the last bit. A column
// started from another profile converges to another column within the tolerance, so the comparison would show it.
TEST(InflowColumn, ColumnStartsFromTheDefaults)
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(4, 4, 32);
    in.mesh.prob_lo = {0, 0, 0};
    in.mesh.prob_hi = {256, 256, 640};
    in.mesh.is_periodic = {1, 1, 0};
    in.mesh.max_grid_size = {4, 4};
    in.terrain.type = "flat";
    InitOptions io;
    io.first_guess = InitOptions::FirstGuess::column;
    io.column_target_u = 5;
    io.column_target_v = 2;
    io.column_target_height = 30;
    gust::SteadyOptions so;
    so.stage = gust::StageKind::simplec;
    gust::SurfaceLayerOptions surf;
    surf.z0 = amrex::Real(0.1);
    gust::ClosureOptions copt;
    gust::ForcingOptions fopt;
    fopt.coriolis.enabled = true;
    fopt.coriolis.latitude = 45;
    fopt.coriolis.three_d = false;
    fopt.coriolis.geostrophic = true;
    fopt.coriolis.ug = 10;
    BoundaryInputs plain;
    plain.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic, FaceType::periodic,
                  FaceType::slip_wall};
    plain.theta_ref.z = {0};
    plain.theta_ref.theta = {300};
    plain.inflow.l_max = copt.max_geom_lscale;
    const auto a = gust::run_column_guess(in, io, so, plain, surf, true, copt, fopt);
    ASSERT_TRUE(a.steady);
    BoundaryInputs own = plain;
    own.inflow.kind = K::column;
    own.inflow.table = sample_table();
    own.inflow.ustar = amrex::Real(0.9);
    own.inflow.z0 = amrex::Real(0.5);
    own.inflow.direction = 30;
    const auto b = gust::run_column_guess(in, io, so, own, surf, true, copt, fopt);
    ASSERT_TRUE(b.steady);
    ASSERT_EQ(a.z.size(), b.z.size());
    for (amrex::Long i = 0; i < a.z.size(); ++i) {
        EXPECT_EQ(a.u[i], b.u[i]) << i;
        EXPECT_EQ(a.v[i], b.v[i]) << i;
        EXPECT_EQ(a.tke[i], b.tke[i]) << i;
    }
    EXPECT_EQ(a.ug, b.ug);
    EXPECT_EQ(a.vg, b.vg);
    EXPECT_EQ(a.iterations, b.iterations);   // the same path, not only the same end
    // The log law with these options (another start) gives another column (in u, v or k, in single precision too):
    // the comparisons above can fail.
    BoundaryInputs other = own;
    other.inflow.kind = K::log_law;
    other.inflow.table = gust::InflowTable{};
    const auto c = gust::run_column_guess(in, io, so, other, surf, true, copt, fopt);
    ASSERT_TRUE(c.steady);
    bool differs = c.ug != a.ug;
    for (amrex::Long i = 0; i < a.z.size(); ++i) {
        differs = differs || c.u[i] != a.u[i] || c.v[i] != a.v[i] || c.tke[i] != a.tke[i];
    }
    EXPECT_TRUE(differs);
}

// A restart deck must keep what the stored column depends on: Coriolis and its latitude, the surface roughness, and
// (without a target wind) the forcing; with a target the record's forcing replaces the deck's.
TEST(InflowColumn, RestartChecks)
{
    InflowColumnRecord r = sample_record();
    r.z0_map = false;
    gust::ForcingOptions f;
    f.coriolis.enabled = true;
    f.coriolis.geostrophic = true;
    f.coriolis.latitude = r.latitude;
    f.coriolis.rotation_period = r.rotation_period;
    f.coriolis.ug = 10;
    gust::SurfaceLayerOptions surf;
    surf.z0 = r.z0;
    auto msg = [&] (const InflowColumnRecord& rec, const gust::ForcingOptions& fopt, const gust::SurfaceLayerOptions& sopt) {
        std::string m;
        for (const auto& e : gust::check_inflow_column_restart(rec, fopt, sopt, "chk9")) { m += e + "\n"; }
        return m;
    };
    EXPECT_EQ(msg(r, f, surf), "");   // a target wind: the deck's ug, vg are replaced
    auto has = [&] (const std::string& m, const std::string& what) {
        EXPECT_NE(m.find("io.restart = chk9: " + what), std::string::npos) << m;
    };
    { auto g = f; g.coriolis.latitude += amrex::Real(1); has(msg(r, g, surf), "the checkpoint's column ran at forcing.latitude"); }
    { auto g = f; g.coriolis.rotation_period *= 2; has(msg(r, g, surf), "the checkpoint's column ran at forcing.latitude"); }
    { auto g = f; g.coriolis.enabled = false; has(msg(r, g, surf), "the checkpoint's column ran with Coriolis"); }
    { auto g = f; g.coriolis.geostrophic = false; has(msg(r, g, surf), "the checkpoint's column ran with a geostrophic wind"); }
    { auto o = surf; o.z0 = amrex::Real(0.5); has(msg(r, f, o), "the checkpoint's column ran on bc.surface_z0"); }
    { auto o = surf; o.z0_file = "z0.txt"; has(msg(r, f, o), "the checkpoint's column ran on bc.surface_z0"); }
    {   // under a map: the deck's map, whatever bc.surface_z0 says (init.column_z0 is refused on a restart)
        auto m = r; m.z0_map = true;
        auto o = surf; o.z0 = amrex::Real(0.5); o.z0_file = "z0.txt";
        EXPECT_EQ(msg(m, f, o), "");
        has(msg(m, f, surf), "the checkpoint's column ran under a roughness map");
    }
    {   // without Coriolis in the record, a deck with it is refused
        auto n = r; n.coriolis = false; n.latitude = 0; n.rotation_period = 0; n.ug = 0; n.vg = 0;
        has(msg(n, f, surf), "the checkpoint's column ran without Coriolis");
    }
    {   // no target wind: the deck's forcing must be the record's, value for value
        auto n = r; n.target = false;
        auto g = f; g.coriolis.ug = n.ug; g.coriolis.vg = n.vg;
        g.pressure_gradient[0] = n.pressure_gradient_x; g.pressure_gradient[1] = n.pressure_gradient_y;
        EXPECT_EQ(msg(n, g, surf), "");
        auto h = g; h.coriolis.ug = std::nextafter(h.coriolis.ug, amrex::Real(100));
        has(msg(n, h, surf), "the checkpoint's column ran with forcing.ug");
        h = g; h.coriolis.vg += amrex::Real(1); has(msg(n, h, surf), "the checkpoint's column ran with forcing.ug");
        h = g; h.pressure_gradient[0] = 0; has(msg(n, h, surf), "the checkpoint's column ran with forcing.ug");
        h = g; h.pressure_gradient[1] = 0; has(msg(n, h, surf), "the checkpoint's column ran with forcing.ug");
        // apply: nothing replaced without a target
        auto a = h; gust::apply_inflow_column_forcing(n, a); EXPECT_EQ(a.pressure_gradient[1], amrex::Real(0));
    }
    // With a target the record's forcing replaces the deck's.
    auto a = f;
    gust::apply_inflow_column_forcing(r, a);
    EXPECT_EQ(a.coriolis.ug, r.ug);
    EXPECT_EQ(a.coriolis.vg, r.vg);
    EXPECT_EQ(a.pressure_gradient[0], r.pressure_gradient_x);
    EXPECT_EQ(a.pressure_gradient[1], r.pressure_gradient_y);
    auto p = f; p.coriolis.geostrophic = false; p.coriolis.ug = 3;
    gust::apply_inflow_column_forcing(r, p);
    EXPECT_EQ(p.coriolis.ug, amrex::Real(3));   // no geostrophic forcing: ug left alone
    EXPECT_EQ(p.pressure_gradient[0], r.pressure_gradient_x);
}

// A wind that leaves through an inflow face is reported per face, with the rows and heights; the ground row is skipped.
TEST(InflowColumn, OutflowRowsWarned)
{
    gust::InflowTable t;
    t.ncol = 5;
    t.z = {0, 10, 100, 1000};
    t.u = {0, 5, 3, -2};    // turns westward aloft
    t.v = {0, 1, 2, 4};
    t.theta = {300, 300, 300, 300};
    t.k = {1, 1, 1, 1};
    gust::BoundaryInputs bin;
    using O = amrex::Orientation;
    bin.face[O(0, O::low)] = gust::FaceType::inflow;
    bin.face[O(1, O::low)] = gust::FaceType::inflow;
    bin.face[O(0, O::high)] = gust::FaceType::outflow;
    bin.face[O(1, O::high)] = gust::FaceType::outflow;
    auto w = gust::inflow_column_outflow_rows(t, bin);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_EQ(w[0], "bc.xlo = inflow: the column's wind leaves through this face in 1 of 3 rows (z 1000 to 1000 m)");
    bin.face[O(1, O::high)] = gust::FaceType::inflow;   // v > 0 leaves through yhi in every row
    w = gust::inflow_column_outflow_rows(t, bin);
    ASSERT_EQ(w.size(), 2u);
    EXPECT_EQ(w[1], "bc.yhi = inflow: the column's wind leaves through this face in 3 of 3 rows (z 10 to 1000 m)");
    t.u[3] = 1;
    bin.face[O(1, O::high)] = gust::FaceType::outflow;
    EXPECT_TRUE(gust::inflow_column_outflow_rows(t, bin).empty());
}
