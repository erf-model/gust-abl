// The named set-ups abl.drive and abl.theta (AblSetup.H): the names, the
// defaults they give, and the table of what each one needs and refuses. The
// runs behind the names are in Tests/Infrastructure/Abl_Setup.
#include "AblSetup.H"
#include "ColumnGuess.H"
#include "Driver.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>

namespace {

using gust::AblSetupFacts;
using gust::AblSetupOptions;
using gust::FaceType;
using Drive = AblSetupOptions::Drive;
using Theta = AblSetupOptions::Theta;

const Drive kDrives[] = {Drive::inflow_only, Drive::couette, Drive::pressure_driven, Drive::geostrophic_column,
                         Drive::geostrophic_inflow};
const Theta kThetas[] = {Theta::neutral, Theta::capped, Theta::frozen, Theta::surface_flux};

const int kXlo = amrex::Orientation(0, amrex::Orientation::low);
const int kXhi = amrex::Orientation(0, amrex::Orientation::high);
const int kZhi = amrex::Orientation(2, amrex::Orientation::high);

AblSetupOptions options (Drive d, Theta t)
{
    AblSetupOptions o;
    o.drive = d;
    o.theta = t;
    return o;
}

// A deck that is exactly the set-up; each test below changes one input of it.
AblSetupFacts deck (Drive d, Theta t)
{
    AblSetupFacts f;
    f.solve = true;
    f.surface_layer = true;
    f.bc.face[kXlo] = FaceType::inflow;
    f.bc.face[kXhi] = FaceType::outflow;
    f.bc.face[kZhi] = FaceType::slip_wall;
    if (d == Drive::couette) {
        f.bc.face[kZhi] = FaceType::no_slip_wall;
        f.bc.wall_velocity[kZhi][0] = 15;
        f.bc.inflow.length = gust::InflowProfileOptions::Length::unbounded;
        f.bc.inflow.z0 = f.surface.z0;
        f.max_geom_lscale = 1.0e6;
    } else if (d == Drive::pressure_driven) {
        f.forcing.pressure_gradient[0] = amrex::Real(1.0e-4);
    } else if (d == Drive::geostrophic_column || d == Drive::geostrophic_inflow) {
        f.forcing.coriolis.enabled = true;
        f.forcing.coriolis.geostrophic = true;
        f.forcing.coriolis.latitude = 45;
        f.forcing.coriolis.ug = 10;
        if (d == Drive::geostrophic_column) {
            f.init.first_guess = gust::InitOptions::FirstGuess::column;
            f.bc.inflow.kind = gust::InflowProfileOptions::Kind::file;
        }
    }
    if (t == Theta::capped) {
        f.bc.theta_ref.z = {0, 700, 1000};
        f.bc.theta_ref.theta = {300, 300, 309};
        f.bc.zhi_theta_grad = amrex::Real(0.003);
        f.bc.rayleigh_zdamp = 800;
        f.bc.rayleigh_dampcoef = amrex::Real(0.2);
    } else if (t == Theta::frozen) {
        f.simplec.theta_frozen = true;
        f.bc.theta_ref.z = {0, 1000};
        f.bc.theta_ref.theta = {300, 305};
        if (d == Drive::pressure_driven) {
            f.bc.rayleigh_zdamp = 800;
            f.bc.rayleigh_dampcoef = amrex::Real(0.2);
        }
    } else if (t == Theta::surface_flux) {
        f.surface.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
        f.surface.heat_flux = amrex::Real(0.05);
        f.bc.given.emplace_back("surface_heat_flux");
        if (d == Drive::geostrophic_column) {
            f.init.column_stability = gust::InitOptions::ColumnStability::heat_flux;
            f.init.column_time = 3600;
        }
    }
    f.inflow_table_theta = f.bc.theta_ref.theta;
    return f;
}

bool mentions (const std::vector<std::string>& msgs, const std::string& s)
{
    return std::any_of(msgs.begin(), msgs.end(), [&] (const std::string& m) { return m.find(s) != std::string::npos; });
}

// Adds keys the way an inputs file does, and removes them and the defaults
// the set-up added when the test ends.
struct ScopedDeck
{
    explicit ScopedDeck (const std::string& text)
    {
        const std::string fname = "ut_abl_inputs_tmp";
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream os(fname);
            os << text;
        }
        amrex::ParallelDescriptor::Barrier();
        amrex::ParmParse::addfile(fname);
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
    }
    ~ScopedDeck ()
    {
        const std::vector<std::pair<std::string, std::string>> keys{
            {"abl", "drive"}, {"abl", "theta"}, {"steady", "stage"}, {"forcing", "coriolis"}, {"forcing", "geostrophic"},
            {"bc", "inflow_tke_length"}, {"bc", "rayleigh_fields"}, {"bc", "surface_stability"},
            {"turb", "max_geom_lscale"}, {"simplec", "theta"}};
        for (const auto& [prefix, name] : keys) {
            amrex::ParmParse pp(prefix);
            while (pp.remove(name) > 0) {}
        }
    }
};

std::string value (const char* prefix, const char* key)
{
    std::string v;
    amrex::ParmParse(prefix).query(key, v);
    return v;
}

TEST (AblSetup, ReadsEveryNameOfTheInputTable)
{
    for (Drive d : kDrives) {
        ScopedDeck guard(std::string("abl.drive = ") + gust::abl_drive_name(d) + "\n");
        const auto o = gust::read_abl_setup_options();
        EXPECT_EQ(o.drive, d) << gust::abl_drive_name(d);
        EXPECT_TRUE(o.errors.empty());
    }
    for (Theta t : kThetas) {
        ScopedDeck guard(std::string("abl.theta = ") + gust::abl_theta_name(t) + "\n");
        const auto o = gust::read_abl_setup_options();
        EXPECT_EQ(o.theta, t) << gust::abl_theta_name(t);
        EXPECT_TRUE(o.errors.empty());
    }
}

TEST (AblSetup, RefusesBadNamesAndRepeatedKeys)
{
    {
        ScopedDeck guard("abl.drive = ekman\nabl.theta = warm\n");
        const auto o = gust::read_abl_setup_options();
        EXPECT_TRUE(mentions(o.errors, "abl.drive = ekman: must be inflow_only"));
        EXPECT_TRUE(mentions(o.errors, "abl.theta = warm: must be neutral"));
        EXPECT_FALSE(o.active());
    }
    {
        ScopedDeck guard("abl.drive = couette\nabl.drive = inflow_only\n");
        EXPECT_TRUE(mentions(gust::read_abl_setup_options().errors, "abl.drive is given more than once"));
    }
}

TEST (AblSetup, WithoutTheKeysNothingIsChecked)
{
    // Coriolis without a geostrophic wind and a heated top: every set-up would refuse this deck.
    AblSetupFacts f = deck(Drive::inflow_only, Theta::neutral);
    f.forcing.coriolis.enabled = true;
    f.bc.zhi_theta_grad = 1;
    EXPECT_TRUE(gust::check_abl_setup(AblSetupOptions{}, f).empty());
}

TEST (AblSetup, RefusedWithoutAStage)
{
    AblSetupFacts f = deck(Drive::inflow_only, Theta::neutral);
    f.solve = false;
    const auto err = gust::check_abl_setup(options(Drive::inflow_only, Theta::neutral), f);
    EXPECT_TRUE(mentions(err, "abl.drive has no effect without steady.stage"));
    EXPECT_TRUE(mentions(err, "abl.theta has no effect without steady.stage"));
}

// The table of docs/inputs.rst: every pair is accepted on its own deck but
// couette with stratification.
TEST (AblSetup, CompatibilityTable)
{
    for (Drive d : kDrives) {
        for (Theta t : kThetas) {
            const AblSetupFacts f = deck(d, t);
            auto err = gust::check_abl_setup(options(d, t), f);
            // The application's own check of the column must accept the same deck.
            gust::SteadyOptions sopt;
            sopt.stage = gust::StageKind::simplec;
            for (auto& e : gust::check_column_guess(f.init, sopt, f.forcing, f.bc, f.surface, f.surface_layer)) {
                err.push_back(std::move(e));
            }
            const std::string pair = std::string(gust::abl_drive_name(d)) + " + " + gust::abl_theta_name(t);
            if (d == Drive::couette && t != Theta::neutral) {
                EXPECT_TRUE(mentions(err, "abl.drive = couette needs abl.theta = neutral")) << pair;
            } else {
                EXPECT_TRUE(err.empty()) << pair << ": " << (err.empty() ? "" : err.front());
            }
        }
    }
}

TEST (AblSetup, DriveRefusals)
{
    struct Case { Drive drive; const char* message; void (*change) (AblSetupFacts&); };
    const Case cases[] = {
        {Drive::inflow_only, "refuses forcing.coriolis = 1", [] (AblSetupFacts& f) { f.forcing.coriolis.enabled = true; }},
        {Drive::inflow_only, "refuses a non-zero forcing.pressure_gradient",
         [] (AblSetupFacts& f) { f.forcing.pressure_gradient[1] = 1; }},
        {Drive::inflow_only, "refuses a non-zero forcing.pressure_gradient",
         [] (AblSetupFacts& f) { f.forcing.pressure_gradient[2] = 1; }},
        {Drive::inflow_only, "needs an inflow face", [] (AblSetupFacts& f) { f.bc.face[kXlo] = FaceType::periodic; }},
        {Drive::couette, "needs an inflow face", [] (AblSetupFacts& f) { f.bc.face[kXlo] = FaceType::periodic; }},
        {Drive::couette, "moving lid", [] (AblSetupFacts& f) { f.bc.face[kZhi] = FaceType::slip_wall; }},
        {Drive::couette, "moving lid", [] (AblSetupFacts& f) { f.bc.wall_velocity[kZhi] = {0, 0, 0}; }},
        {Drive::couette, "bc.inflow_profile = log_law",
         [] (AblSetupFacts& f) { f.bc.inflow.kind = gust::InflowProfileOptions::Kind::power_law; }},
        {Drive::couette, "bc.inflow_tke_length = unbounded",
         [] (AblSetupFacts& f) { f.bc.inflow.length = gust::InflowProfileOptions::Length::closure; }},
        {Drive::couette, "turb.max_geom_lscale", [] (AblSetupFacts& f) { f.max_geom_lscale = 30; }},
        {Drive::couette, "bc.inflow_z0", [] (AblSetupFacts& f) { f.bc.inflow.z0 = amrex::Real(0.03); }},
        {Drive::couette, "bc.inflow_z0", [] (AblSetupFacts& f) { f.surface.z0_file = "map"; }},
        {Drive::couette, "along bc.inflow_direction", [] (AblSetupFacts& f) { f.bc.inflow.direction = 30; }},
        {Drive::couette, "along bc.inflow_direction", [] (AblSetupFacts& f) { f.bc.wall_velocity[kZhi] = {0, 15, 0}; }},
        {Drive::pressure_driven, "refuses a Rayleigh layer on u or v",
         [] (AblSetupFacts& f) { f.bc.rayleigh_zdamp = 800; f.bc.rayleigh_dampcoef = 1; }},
        {Drive::pressure_driven, "needs a Rayleigh layer on u and v",
         [] (AblSetupFacts& f) { f.bc.theta_ref.z = {0, 1}; f.bc.theta_ref.theta = {300, 301}; }},
        {Drive::geostrophic_inflow, "needs a wind", [] (AblSetupFacts& f) { f.forcing.coriolis.ug = 0; f.init.column_target_u = 8; }},
        {Drive::couette, "bc.zlo = surface_layer", [] (AblSetupFacts& f) { f.surface_layer = false; }},
        {Drive::pressure_driven, "needs a non-zero horizontal", [] (AblSetupFacts& f) { f.forcing.pressure_gradient = {0, 0, 1}; }},
        {Drive::pressure_driven, "bc.zhi = slip_wall", [] (AblSetupFacts& f) { f.bc.face[kZhi] = FaceType::no_slip_wall; }},
        {Drive::pressure_driven, "refuses forcing.coriolis = 1", [] (AblSetupFacts& f) { f.forcing.coriolis.enabled = true; }},
        {Drive::geostrophic_column, "needs forcing.coriolis = 1", [] (AblSetupFacts& f) { f.forcing.coriolis.enabled = false; }},
        {Drive::geostrophic_column, "needs forcing.coriolis = 1", [] (AblSetupFacts& f) { f.forcing.coriolis.geostrophic = false; }},
        {Drive::geostrophic_column, "needs a wind", [] (AblSetupFacts& f) { f.forcing.coriolis.ug = 0; }},
        {Drive::geostrophic_column, "needs init.first_guess = column",
         [] (AblSetupFacts& f) { f.init.first_guess = gust::InitOptions::FirstGuess::profile; }},
        {Drive::geostrophic_column, "bc.inflow_profile = file",
         [] (AblSetupFacts& f) { f.bc.inflow.kind = gust::InflowProfileOptions::Kind::log_law; }},
        {Drive::geostrophic_column, "refuses a non-zero forcing.pressure_gradient",
         [] (AblSetupFacts& f) { f.forcing.pressure_gradient[0] = 1; }},
        {Drive::geostrophic_inflow, "needs an inflow face", [] (AblSetupFacts& f) { f.bc.face[kXlo] = FaceType::periodic; }},
        {Drive::geostrophic_inflow, "needs a wind", [] (AblSetupFacts& f) { f.forcing.coriolis.ug = 0; }},
    };
    for (const auto& c : cases) {
        AblSetupFacts f = deck(c.drive, Theta::neutral);
        c.change(f);
        const auto err = gust::check_abl_setup(options(c.drive, Theta::unset), f);
        EXPECT_TRUE(mentions(err, c.message))
            << gust::abl_drive_name(c.drive) << ": " << c.message;
    }
}

TEST (AblSetup, DriveAccepts)
{
    // A column target wind stands for the geostrophic wind, on a deck without inflow faces.
    AblSetupFacts f = deck(Drive::geostrophic_column, Theta::neutral);
    f.bc.face[kXlo] = FaceType::periodic;
    f.bc.face[kXhi] = FaceType::periodic;
    f.bc.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    f.forcing.coriolis.ug = 0;
    f.init.column_target_u = 8;
    EXPECT_TRUE(gust::check_abl_setup(options(Drive::geostrophic_column, Theta::neutral), f).empty());
    // A restart replaces the column first guess.
    f = deck(Drive::geostrophic_column, Theta::neutral);
    f.init.first_guess = gust::InitOptions::FirstGuess::none;
    f.restart = true;
    EXPECT_TRUE(gust::check_abl_setup(options(Drive::geostrophic_column, Theta::neutral), f).empty());
    // A periodic deck needs no table.
    f = deck(Drive::geostrophic_column, Theta::neutral);
    f.bc.face[kXlo] = FaceType::periodic;
    f.bc.face[kXhi] = FaceType::periodic;
    f.bc.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    EXPECT_TRUE(gust::check_abl_setup(options(Drive::geostrophic_column, Theta::neutral), f).empty());
}

TEST (AblSetup, ThetaRefusals)
{
    struct Case { Theta theta; const char* message; void (*change) (AblSetupFacts&); };
    const Case cases[] = {
        {Theta::neutral, "uniform bc.theta_profile", [] (AblSetupFacts& f) { f.bc.theta_ref.z = {0, 1}; f.bc.theta_ref.theta = {300, 301}; }},
        {Theta::neutral, "bc.zhi_theta_grad", [] (AblSetupFacts& f) { f.bc.zhi_theta_grad = amrex::Real(0.003); }},
        {Theta::neutral, "canopy.heat", [] (AblSetupFacts& f) { f.canopy_heat = true; }},
        {Theta::neutral, "refuses simplec.theta = frozen", [] (AblSetupFacts& f) { f.simplec.theta_frozen = true; }},
        {Theta::neutral, "refuses bc.surface_stability = heat_flux",
         [] (AblSetupFacts& f) { f.surface.stability = gust::SurfaceLayerOptions::Stability::heat_flux; }},
        {Theta::capped, "rises with height", [] (AblSetupFacts& f) { f.bc.theta_ref.theta = {300, 300, 300}; }},
        {Theta::capped, "rises with height", [] (AblSetupFacts& f) { f.bc.theta_ref.theta = {305, 300, 309}; }},
        {Theta::capped, "bc.zhi_theta_grad > 0", [] (AblSetupFacts& f) { f.bc.zhi_theta_grad = 0; }},
        {Theta::capped, "Rayleigh layer", [] (AblSetupFacts& f) { f.bc.rayleigh_zdamp = 0; }},
        {Theta::capped, "Rayleigh layer", [] (AblSetupFacts& f) { f.bc.rayleigh_dampcoef = 0; }},
        {Theta::capped, "Rayleigh layer", [] (AblSetupFacts& f) { f.bc.rayleigh_fields[2] = false; }},
        {Theta::capped, "Rayleigh layer", [] (AblSetupFacts& f) { f.bc.rayleigh_fields[3] = false; }},
        {Theta::capped, "forcing.gravity = 0", [] (AblSetupFacts& f) { f.forcing.buoyancy.gravity = 0; }},
        {Theta::frozen, "needs simplec.theta = frozen", [] (AblSetupFacts& f) { f.simplec.theta_frozen = false; }},
        {Theta::frozen, "stratified bc.theta_profile", [] (AblSetupFacts& f) { f.bc.theta_ref.theta = {300, 300}; }},
        {Theta::surface_flux, "needs bc.zlo = surface_layer", [] (AblSetupFacts& f) { f.surface_layer = false; }},
        {Theta::surface_flux, "needs bc.zlo = surface_layer",
         [] (AblSetupFacts& f) { f.surface.stability = gust::SurfaceLayerOptions::Stability::neutral; }},
        {Theta::surface_flux, "bc.surface_heat_flux or init.column_mol", [] (AblSetupFacts& f) { f.surface.heat_flux = 0; }},
        {Theta::surface_flux, "bc.surface_heat_flux or init.column_mol", [] (AblSetupFacts& f) { f.bc.given.clear(); }},
        {Theta::surface_flux, "a heat sink", [] (AblSetupFacts& f) { f.bc.face[kXhi] = FaceType::slip_wall; }},
    };
    for (const auto& c : cases) {
        AblSetupFacts f = deck(Drive::inflow_only, c.theta);
        c.change(f);
        EXPECT_TRUE(mentions(gust::check_abl_setup(options(Drive::unset, c.theta), f), c.message))
            << gust::abl_theta_name(c.theta) << ": " << c.message;
    }
}

TEST (AblSetup, ThetaWithTheColumnDrive)
{
    using D = Drive;
    // A frozen theta needs a table to hold: the column's theta is not the one held.
    AblSetupFacts f = deck(D::geostrophic_column, Theta::frozen);
    f.bc.face[kXlo] = FaceType::periodic;
    f.bc.face[kXhi] = FaceType::periodic;
    f.bc.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(D::geostrophic_column, Theta::frozen), f), "(the theta to hold)"));
    // Column_Marching's form: periodic, no Rayleigh layer, the column marches to a fixed time.
    f = deck(D::geostrophic_column, Theta::surface_flux);
    f.bc.face[kXlo] = FaceType::periodic;
    f.bc.face[kXhi] = FaceType::periodic;
    EXPECT_TRUE(gust::check_abl_setup(options(D::geostrophic_column, Theta::surface_flux), f).empty());
    // The same deck as a steady column: its heat has no way out.
    f.init.column_time = 0;
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(D::geostrophic_column, Theta::surface_flux), f), "a heat sink"));
    // A steady column with a fixed flux and a Rayleigh layer on theta is steady.
    f.bc.rayleigh_zdamp = 800;
    f.bc.rayleigh_dampcoef = amrex::Real(0.2);
    EXPECT_TRUE(gust::check_abl_setup(options(D::geostrophic_column, Theta::surface_flux), f).empty());
    // A column that takes a Monin-Obukhov length hands the run its flux.
    f.init.column_stability = gust::InitOptions::ColumnStability::mol;
    f.init.column_mol = 200;
    f.bc.given.clear();
    f.surface.heat_flux = 0;
    EXPECT_TRUE(gust::check_abl_setup(options(D::geostrophic_column, Theta::surface_flux), f).empty());
    // The other treatments of theta refuse a stratified column.
    f = deck(D::geostrophic_column, Theta::neutral);
    f.init.column_mol = 200;
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(D::geostrophic_column, Theta::neutral), f), "stratified column"));
    f = deck(D::geostrophic_column, Theta::capped);
    f.init.column_stability = gust::InitOptions::ColumnStability::surface_temperature;
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(D::geostrophic_column, Theta::capped), f), "stratified column"));
    // A Rayleigh layer on theta is the heat sink of a periodic deck.
    f = deck(D::inflow_only, Theta::surface_flux);
    f.bc.face[kXhi] = FaceType::slip_wall;
    f.bc.rayleigh_zdamp = 800;
    f.bc.rayleigh_dampcoef = amrex::Real(0.2);
    EXPECT_TRUE(gust::check_abl_setup(options(D::inflow_only, Theta::surface_flux), f).empty());
}

// With a file inflow the air's theta is the table's; bc.theta_profile is the reference only.
TEST (AblSetup, TheInflowTableIsTheTheta)
{
    using Kind = gust::InflowProfileOptions::Kind;
    AblSetupFacts f = deck(Drive::inflow_only, Theta::neutral);
    f.bc.inflow.kind = Kind::file;
    f.inflow_table_theta = {300, 305};
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(Drive::unset, Theta::neutral), f), "a uniform table"));

    f = deck(Drive::inflow_only, Theta::capped);
    f.bc.inflow.kind = Kind::file;
    f.bc.theta_ref.z = {0};
    f.bc.theta_ref.theta = {300};
    f.inflow_table_theta = {300, 300, 309};
    EXPECT_TRUE(gust::check_abl_setup(options(Drive::unset, Theta::capped), f).empty());
    f.bc.theta_ref.z = {0, 1000};
    f.bc.theta_ref.theta = {300, 309};
    f.inflow_table_theta = {300, 300, 300};
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(Drive::unset, Theta::capped), f), "rises with height"));

    f = deck(Drive::inflow_only, Theta::frozen);
    f.bc.inflow.kind = Kind::file;
    f.inflow_table_theta = {300, 300};
    EXPECT_TRUE(mentions(gust::check_abl_setup(options(Drive::unset, Theta::frozen), f), "stratified"));
    f.bc.theta_ref.theta = {300, 300};
    f.inflow_table_theta = {300, 305};
    EXPECT_TRUE(gust::check_abl_setup(options(Drive::unset, Theta::frozen), f).empty());
}

TEST (AblSetup, AfterABadNameOnlyThatIsReported)
{
    // The defaults were not set, so a good name's checks would report them.
    ScopedDeck guard("abl.drive = geostrophic_inflow\nabl.theta = warm\n");
    const auto o = gust::read_abl_setup_options();
    const auto err = gust::check_abl_setup(o, deck(Drive::inflow_only, Theta::neutral));
    ASSERT_EQ(err.size(), 1U);
    EXPECT_TRUE(mentions(err, "abl.theta = warm"));
}

TEST (AblSetup, RefusesExtraValuesAndContradictingDeckValues)
{
    {
        ScopedDeck guard("abl.drive = couette neutral\n");
        EXPECT_TRUE(mentions(gust::read_abl_setup_options().errors, "abl.drive takes exactly 1 value"));
    }
    {
        ScopedDeck guard("abl.drive = geostrophic_inflow\nforcing.coriolis = 0\nabl.theta = frozen\nsimplec.theta = solve\n");
        const auto err = gust::read_abl_setup_options().errors;
        EXPECT_TRUE(mentions(err, "abl.drive = geostrophic_inflow refuses forcing.coriolis = 0"));
        EXPECT_TRUE(mentions(err, "abl.theta = frozen refuses simplec.theta = solve"));
    }
    {
        ScopedDeck guard("abl.drive = geostrophic_column\nforcing.geostrophic = 0\n");
        EXPECT_TRUE(mentions(gust::read_abl_setup_options().errors, "refuses forcing.geostrophic = 0"));
    }
    {
        ScopedDeck guard("abl.theta = surface_flux\nbc.surface_stability = neutral\n");
        EXPECT_TRUE(mentions(gust::read_abl_setup_options().errors, "refuses bc.surface_stability = neutral"));
    }
}

TEST (AblSetup, GeostrophicDrivesSwitchOnTheForcing)
{
    for (const char* name : {"geostrophic_column", "geostrophic_inflow"}) {
        ScopedDeck guard(std::string("steady.stage = simplec\nabl.drive = ") + name + "\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        int cor = 0, geo = 0;
        amrex::ParmParse("forcing").query("coriolis", cor);
        amrex::ParmParse("forcing").query("geostrophic", geo);
        EXPECT_EQ(cor, 1) << name;
        EXPECT_EQ(geo, 1) << name;
    }
}

TEST (AblSetup, CouetteUncapsTheLength)
{
    ScopedDeck guard("steady.stage = simplec\nabl.drive = couette\n");
    gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
    EXPECT_EQ(value("bc", "inflow_tke_length"), "unbounded");
    amrex::Real l = 0;
    amrex::ParmParse("turb").query("max_geom_lscale", l);
    EXPECT_EQ(l, amrex::Real(1.0e6));
}

TEST (AblSetup, ThetaDefaults)
{
    {
        ScopedDeck guard("steady.stage = simplec\nabl.theta = capped\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        std::vector<std::string> fields;
        amrex::ParmParse("bc").queryarr("rayleigh_fields", fields);
        EXPECT_EQ(fields, (std::vector<std::string>{"w", "theta"}));
    }
    {
        ScopedDeck guard("steady.stage = simplec\nabl.theta = frozen\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        EXPECT_EQ(value("simplec", "theta"), "frozen");
    }
    {
        ScopedDeck guard("steady.stage = simplec\nabl.theta = surface_flux\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        EXPECT_EQ(value("bc", "surface_stability"), "heat_flux");
    }
}

TEST (AblSetup, TheDeckKeepsItsValue)
{
    ScopedDeck guard("steady.stage = simplec\nabl.drive = geostrophic_inflow\nforcing.coriolis = 0\nabl.theta = frozen\nsimplec.theta = solve\n");
    gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
    amrex::ParmParse forcing("forcing");
    EXPECT_EQ(forcing.countname("coriolis"), 1);
    int cor = 1;
    forcing.query("coriolis", cor);
    EXPECT_EQ(cor, 0);
    EXPECT_EQ(amrex::ParmParse("simplec").countname("theta"), 1);
    EXPECT_EQ(value("simplec", "theta"), "solve");
}

TEST (AblSetup, NoDefaultsWithoutAStageOrWithABadName)
{
    {
        ScopedDeck guard("abl.drive = geostrophic_inflow\nabl.theta = frozen\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        EXPECT_FALSE(amrex::ParmParse("forcing").contains("coriolis"));
        EXPECT_FALSE(amrex::ParmParse("simplec").contains("theta"));
    }
    {
        ScopedDeck guard("steady.stage = none\nabl.drive = geostrophic_inflow\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        EXPECT_FALSE(amrex::ParmParse("forcing").contains("coriolis"));
    }
    {
        ScopedDeck guard("steady.stage = simplec\nabl.drive = geostrophic_inflow\nabl.theta = warm\n");
        gust::apply_abl_setup_defaults(gust::read_abl_setup_options());
        EXPECT_FALSE(amrex::ParmParse("forcing").contains("coriolis"));
    }
}

} // namespace
