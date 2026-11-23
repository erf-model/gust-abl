#include "InflowColumn.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace gust {

namespace {

const char* const kMagic = "GustInflowColumn";
constexpr int kVersion = 1;

//! A string from the I/O rank on every rank.
void bcast_string (std::string& s)
{
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Long n = static_cast<amrex::Long>(s.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

std::string g17 (amrex::Real x)
{
    char b[40];
    std::snprintf(b, sizeof(b), "%.17g", static_cast<double>(x));
    return b;
}

} // namespace

amrex::Real column_roughness (const InitOptions& iopt, const SurfaceLayerOptions& surf)
{
    return surf.z0_file.empty() ? surf.z0 : iopt.column_z0;
}

std::vector<std::string> check_inflow_column (const BoundaryInputs& bin, const InitOptions& iopt, bool restart, bool nests)
{
    std::vector<std::string> err;
    if (bin.inflow.kind != InflowProfileOptions::Kind::column) { return err; }
    const std::string pre = "bc.inflow_profile = column ";
    if (!restart && iopt.first_guess != InitOptions::FirstGuess::column) {
        err.emplace_back(pre + "needs init.first_guess = column (the run solves its column at start-up) or io.restart (the "
                         "checkpoint holds the column)");
    }
    const bool inflow = std::find(bin.face.begin(), bin.face.end(), FaceType::inflow) != bin.face.end();
    if (!inflow) {
        err.emplace_back(pre + "needs an inflow face (bc.xlo, bc.xhi, bc.ylo or bc.yhi = inflow); without one, "
                         "init.first_guess = column alone gives the first guess");
    }
    // The column's own Rayleigh layer relaxes u, v and k towards the deck's inflow profile (ColumnGuess.H), which is
    // the column being solved.
    if (bin.rayleigh_zdamp > 0 && (bin.rayleigh_fields[0] || bin.rayleigh_fields[1] || bin.rayleigh_fields[4])) {
        err.emplace_back(pre + "with a Rayleigh layer on u, v or tke: the column's layer would relax towards the profile it "
                         "is solving for; leave u, v and tke out of bc.rayleigh_fields");
    }
    if (nests) {
        err.emplace_back(pre + "is not supported with nest.names (the nests' record of the outer run's forcing and faces "
                         "does not hold the column)");
    }
    return err;
}

InflowColumnRecord make_inflow_column_record (const ColumnGuessResult& cg, amrex::Real z0)
{
    InflowColumnRecord rec;
    rec.table.z = cg.z;
    rec.table.u = cg.u;
    rec.table.v = cg.v;
    rec.table.theta = cg.theta;
    rec.table.k = cg.tke;
    rec.table.ncol = 5;
    rec.z0 = z0;
    rec.target = cg.target_iterations > 0;
    rec.ug = cg.ug;
    rec.vg = cg.vg;
    rec.pressure_gradient_x = cg.pressure_gradient_x;
    rec.pressure_gradient_y = cg.pressure_gradient_y;
    return rec;
}

void use_inflow_column (const InflowColumnRecord& rec, InflowProfileOptions& opt)
{
    opt.kind = InflowProfileOptions::Kind::column;
    opt.table = rec.table;
    opt.z0 = rec.z0;
}

std::string format_inflow_column_record (const InflowColumnRecord& rec)
{
    std::string s = std::string(kMagic) + " " + std::to_string(kVersion) + "\n";
    s += "z0 " + g17(rec.z0) + "\n";
    s += "forcing " + std::string(rec.target ? "1" : "0") + " " + g17(rec.ug) + " " + g17(rec.vg) + " "
         + g17(rec.pressure_gradient_x) + " " + g17(rec.pressure_gradient_y) + "\n";
    s += "rows " + std::to_string(rec.table.z.size()) + "\n";
    for (amrex::Long i = 0; i < rec.table.z.size(); ++i) {
        s += g17(rec.table.z[i]) + " " + g17(rec.table.u[i]) + " " + g17(rec.table.v[i]) + " " + g17(rec.table.theta[i])
             + " " + g17(rec.table.k[i]) + "\n";
    }
    return s;
}

InflowColumnRecord parse_inflow_column_record (const std::string& text, const std::string& where)
{
    const std::string pre = "gust: " + where + ": ";
    std::istringstream is(text);
    InflowColumnRecord rec;
    std::string magic, key;
    int version = 0;
    if (!(is >> magic >> version) || magic != kMagic || version != kVersion) {
        amrex::Abort(pre + "not a " + std::string(kMagic) + " record of version " + std::to_string(kVersion));
    }
    // Every number through strtod (the C locale), each a whole token, as the table reader takes them.
    auto number = [&] (const char* what) {
        std::string tok;
        if (!(is >> tok)) { amrex::Abort(pre + "the record ends before " + what); }
        char* stop = nullptr;
        const double x = std::strtod(tok.c_str(), &stop);
        // Finite after the cast to Real too (1e39 overflows a float).
        if (stop == tok.c_str() || *stop != '\0' || !std::isfinite(x) || !std::isfinite(static_cast<amrex::Real>(x))) {
            amrex::Abort(pre + "bad " + what + " '" + tok + "'");
        }
        return static_cast<amrex::Real>(x);
    };
    // A whole token of decimal digits (the forcing flag, the row count).
    auto count = [&] (const char* what) {
        std::string tok;
        if (!(is >> tok)) { amrex::Abort(pre + "the record ends before " + what); }
        if (tok.empty() || tok.size() > 9 || tok.find_first_not_of("0123456789") != std::string::npos) {
            amrex::Abort(pre + "bad " + what + " '" + tok + "'");
        }
        return std::stol(tok);
    };
    auto expect = [&] (const char* k) {
        if (!(is >> key) || key != k) { amrex::Abort(pre + "expected '" + std::string(k) + "'"); }
    };
    expect("z0");
    rec.z0 = number("z0");
    expect("forcing");
    const long flag = count("the forcing flag");
    if (flag != 0 && flag != 1) { amrex::Abort(pre + "the forcing flag must be 0 or 1"); }
    rec.target = flag == 1;
    rec.ug = number("ug");
    rec.vg = number("vg");
    rec.pressure_gradient_x = number("the pressure gradient (x)");
    rec.pressure_gradient_y = number("the pressure gradient (y)");
    expect("rows");
    const long n = count("the row count");
    if (n < 2) { amrex::Abort(pre + "the row count must be at least 2"); }
    for (long i = 0; i < n; ++i) {
        rec.table.z.push_back(number("z"));
        rec.table.u.push_back(number("u"));
        rec.table.v.push_back(number("v"));
        rec.table.theta.push_back(number("theta"));
        rec.table.k.push_back(number("k"));
    }
    std::string extra;
    if (is >> extra) { amrex::Abort(pre + "unexpected text after the last row: '" + extra + "'"); }
    rec.table.ncol = 5;
    if (!(rec.z0 > 0)) { amrex::Abort(pre + "z0 must be > 0"); }
    // The table's rows as the table kind needs them, with the record named (InflowProfile checks them again).
    for (long i = 0; i < n; ++i) {
        const std::string row = "row " + std::to_string(i + 1) + ": ";
        if (i == 0 && rec.table.z[0] != 0) { amrex::Abort(pre + row + "the first height must be 0"); }
        if (i > 0 && !(rec.table.z[i] > rec.table.z[i - 1])) { amrex::Abort(pre + row + "the heights must increase"); }
        if (!(rec.table.theta[i] > 0)) { amrex::Abort(pre + row + "theta must be > 0"); }
        if (!(rec.table.k[i] >= 0)) { amrex::Abort(pre + row + "k must be >= 0"); }
    }
    return rec;
}

std::string inflow_column_checkpoint_file (const std::string& dir)
{
    return dir + "/InflowColumn";
}

void write_inflow_column_checkpoint (const std::string& dir, const InflowColumnRecord& rec)
{
    if (rec.table.ncol != 5 || rec.table.z.size() < 2) {
        amrex::Abort("gust: " + dir + ": the inflow column record is not set (bc.inflow_profile = column)");
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = inflow_column_checkpoint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << format_inflow_column_record(rec);
        if (!os.good()) { amrex::Abort("gust: " + name + ": the record could not be written"); }
    }
    amrex::ParallelDescriptor::Barrier();
}

namespace {
//! The record's text on the I/O rank; found broadcast to every rank (collective).
int read_record_text (const std::string& dir, std::string& text)
{
    int found = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(inflow_column_checkpoint_file(dir), std::ios::binary);
        if (is) {
            found = 1;
            std::ostringstream ss;
            ss << is.rdbuf();
            text = ss.str();
        }
    }
    amrex::ParallelDescriptor::Bcast(&found, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return found;
}
} // namespace

InflowColumnRecord read_inflow_column_checkpoint_if_present (const std::string& dir)
{
    std::string text;
    if (read_record_text(dir, text) == 0) { return InflowColumnRecord{}; }
    bcast_string(text);
    return parse_inflow_column_record(text, "io.restart = " + dir + " (InflowColumn)");
}

InflowColumnRecord read_inflow_column_checkpoint (const std::string& dir, bool active)
{
    std::string text;
    const int found = read_record_text(dir, text);
    const std::string where = "io.restart = " + dir;
    if (active && found == 0) {
        amrex::Abort("gust: " + where + ": the checkpoint was written without bc.inflow_profile = column (no "
                     "InflowColumn record); a restart cannot solve the column (init.* keys are refused on a restart)");
    }
    if (!active && found != 0) {
        amrex::Abort("gust: " + where + ": the checkpoint was written with bc.inflow_profile = column; keep "
                     "bc.inflow_profile = column in the restart deck");
    }
    if (!active) { return InflowColumnRecord{}; }
    bcast_string(text);
    return parse_inflow_column_record(text, where + " (InflowColumn)");
}

} // namespace gust
