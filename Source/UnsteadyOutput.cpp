#include "UnsteadyOutput.H"

#include "FileMagic.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

namespace gust {

namespace {

constexpr const char* kRecordMagic = "gust-unsteady";
constexpr int kRecordVersion = 1;

std::string hexfloat (double x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", x);
    return buf;
}

std::string exact (double x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", x);
    return buf;
}

// One add of the averages on the valid cells of one box: x from q (face values averaged to the
// centre), then West's update with a = w/W (W after the add) and the weight w.
void add_box (const amrex::Box& bx, const amrex::Array4<amrex::Real>& mean, const amrex::Array4<amrex::Real>& m2,
              const amrex::Array4<const amrex::Real>& u, const amrex::Array4<const amrex::Real>& v,
              const amrex::Array4<const amrex::Real>& w, const amrex::Array4<const amrex::Real>& th,
              const amrex::Array4<const amrex::Real>& k, amrex::Real a, amrex::Real wt)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        constexpr amrex::Real half = amrex::Real(0.5);
        const amrex::Real x[kUnsteadyFields] = {half * (u(i, j, kk) + u(i + 1, j, kk)), half * (v(i, j, kk) + v(i, j + 1, kk)),
                                                half * (w(i, j, kk) + w(i, j, kk + 1)), th(i, j, kk), k(i, j, kk)};
        for (int n = 0; n < kUnsteadyFields; ++n) {
            const amrex::Real d = x[n] - mean(i, j, kk, n);
            mean(i, j, kk, n) += a * d;
            m2(i, j, kk, n) += wt * d * (x[n] - mean(i, j, kk, n));
        }
    });
}

void variance_box (const amrex::Box& bx, const amrex::Array4<amrex::Real>& out, const amrex::Array4<const amrex::Real>& m2,
                   amrex::Real inv_w)
{
    amrex::ParallelFor(bx, kUnsteadyFields, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        out(i, j, k, kUnsteadyFields + n) = m2(i, j, k, n) * inv_w;
    });
}

// The height above the ground of cell (i, j, k)'s centre (the mean of its eight nodes) in the
// column whose ground is node level klo.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ground_of (const amrex::Array4<const amrex::Real>& z, int i, int j, int k)
{
    return amrex::Real(0.25) * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real centre_of (const amrex::Array4<const amrex::Real>& z, int i, int j, int k)
{
    return amrex::Real(0.5) * (ground_of(z, i, j, k) + ground_of(z, i, j, k + 1));
}

// Over the column (i, j, klo..khi) in the local boxes: the ground and top heights (the mean of
// the column's four bottom and top nodes) and the smallest |h_c - h|. Ranks without the column
// give +inf. Collective.
struct ColumnFacts { amrex::Real ground, top, dmin; };

ColumnFacts column_facts (const amrex::MultiFab& z_phys_nd, const amrex::Box& col, amrex::Real h)
{
    const int klo = col.smallEnd(2);
    const int khi = col.bigEnd(2);
    constexpr amrex::Real big = std::numeric_limits<amrex::Real>::max();
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpMin> op;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> data(op);
    using T = amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real>::Type;
    // The nodal MultiFab's boxes are the cells' boxes made nodal: the cell box is the valid box
    // shrunk back to cells.
    for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        const amrex::Box b = amrex::enclosedCells(mfi.validbox()) & col;
        if (!b.ok()) { continue; }
        const auto z = z_phys_nd.const_array(mfi);
        op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> T
        {
            const amrex::Real g = ground_of(z, i, j, klo);
            const amrex::Real ground = (k == klo) ? g : big;
            const amrex::Real top = (k == khi) ? ground_of(z, i, j, khi + 1) : big;
            return {ground, top, amrex::Math::abs(centre_of(z, i, j, k) - g - h)};
        });
    }
    T r = data.value(op);
    amrex::Real v[3] = {amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r)};
    amrex::ParallelAllReduce::Min(v, 3, amrex::ParallelDescriptor::Communicator());
    return {v[0], v[1], v[2]};
}

// The lowest k of the column whose |h_c - h| is at most dmin, and its h_c. Collective.
std::pair<int, amrex::Real> nearest_cell (const amrex::MultiFab& z_phys_nd, const amrex::Box& col, amrex::Real h,
                                          amrex::Real dmin)
{
    const int klo = col.smallEnd(2);
    amrex::ReduceOps<amrex::ReduceOpMin> op;
    amrex::ReduceData<int> data(op);
    using T = amrex::ReduceData<int>::Type;
    for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        const amrex::Box b = amrex::enclosedCells(mfi.validbox()) & col;
        if (!b.ok()) { continue; }
        const auto z = z_phys_nd.const_array(mfi);
        op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> T
        {
            const amrex::Real g = ground_of(z, i, j, klo);
            return {amrex::Math::abs(centre_of(z, i, j, k) - g - h) <= dmin ? k : std::numeric_limits<int>::max()};
        });
    }
    int k = amrex::get<0>(data.value(op));
    amrex::ParallelAllReduce::Min(k, amrex::ParallelDescriptor::Communicator());
    // The cell's height: a sum over its one cell (0 on the other ranks).
    amrex::ReduceOps<amrex::ReduceOpSum> sop;
    amrex::ReduceData<amrex::Real> sdata(sop);
    using S = amrex::ReduceData<amrex::Real>::Type;
    const amrex::Box cell(amrex::IntVect(col.smallEnd(0), col.smallEnd(1), k), amrex::IntVect(col.smallEnd(0), col.smallEnd(1), k));
    for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        const amrex::Box b = amrex::enclosedCells(mfi.validbox()) & cell;
        if (!b.ok()) { continue; }
        const auto z = z_phys_nd.const_array(mfi);
        sop.eval(b, sdata, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept -> S
        {
            return {centre_of(z, i, j, kk) - ground_of(z, i, j, klo)};
        });
    }
    amrex::Real hc = amrex::get<0>(sdata.value(sop));
    amrex::ParallelAllReduce::Sum(hc, amrex::ParallelDescriptor::Communicator());
    return {k, hc};
}

void probe_error (const std::string& file, const std::string& what)
{
    amrex::Abort("gust: unsteady.probes_file = " + file + ": " + what);
}

} // namespace

const std::vector<std::string>& unsteady_field_names ()
{
    static const std::vector<std::string> names{"u", "v", "w", "theta", "tke"};
    return names;
}

// ---------------------------------------------------------------------------------------------
// RunningAverages

RunningAverages::RunningAverages (amrex::Real average_start, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : m_start(average_start), m_mean(ba, dm, kUnsteadyFields, 0), m_m2(ba, dm, kUnsteadyFields, 0)
{
    m_mean.setVal(0);
    m_m2.setVal(0);
}

void RunningAverages::add (const State& q, double t0, double t1)
{
    const double w = t1 - std::max(t0, double(m_start));
    if (!(w > 0)) { return; }
    m_weight += w;
    const auto a = static_cast<amrex::Real>(w / m_weight);
    const auto wt = static_cast<amrex::Real>(w);
    for (amrex::MFIter mfi(m_mean, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        add_box(mfi.tilebox(), m_mean.array(mfi), m_m2.array(mfi), q.u.const_array(mfi), q.v.const_array(mfi),
                q.w.const_array(mfi), q.theta.const_array(mfi), q.tke.const_array(mfi), a, wt);
    }
}

void RunningAverages::plot_fields (amrex::MultiFab& out, std::vector<std::string>& names) const
{
    out.define(m_mean.boxArray(), m_mean.DistributionMap(), 2 * kUnsteadyFields, 0);
    amrex::MultiFab::Copy(out, m_mean, 0, 0, kUnsteadyFields, 0);
    if (m_weight > 0) {
        const auto inv_w = static_cast<amrex::Real>(1.0 / m_weight);
        for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            variance_box(mfi.tilebox(), out.array(mfi), m_m2.const_array(mfi), inv_w);
        }
    } else {
        out.setVal(0, kUnsteadyFields, kUnsteadyFields, 0);
    }
    names.clear();
    for (const auto& n : unsteady_field_names()) { names.push_back(n + "_mean"); }
    for (const auto& n : unsteady_field_names()) { names.push_back(n + "_var"); }
}

void RunningAverages::write_checkpoint (const std::string& dir) const
{
    amrex::VisMF::Write(m_mean, amrex::MultiFabFileFullPrefix(0, dir, "Level_", "unsteady_mean"));
    amrex::VisMF::Write(m_m2, amrex::MultiFabFileFullPrefix(0, dir, "Level_", "unsteady_m2"));
}

void RunningAverages::read_checkpoint (const std::string& dir, double weight)
{
    for (auto* pair : {&m_mean, &m_m2}) {
        const char* name = (pair == &m_mean) ? "unsteady_mean" : "unsteady_m2";
        const std::string file = amrex::MultiFabFileFullPrefix(0, dir, "Level_", name);
        if (!amrex::VisMF::Exist(file)) {
            amrex::Abort("gust: io.restart = " + dir + ": the unsteady record is there but " + name + " is missing");
        }
        amrex::MultiFab src;   // defined from the file's own layout
        amrex::VisMF::Read(src, file);
        if (src.nComp() != kUnsteadyFields) {
            amrex::Abort("gust: io.restart = " + dir + ": " + name + " has " + std::to_string(src.nComp()) + " components, not "
                         + std::to_string(kUnsteadyFields));
        }
        pair->ParallelCopy(src, 0, 0, kUnsteadyFields, 0, 0);
    }
    m_weight = weight;
}

// ---------------------------------------------------------------------------------------------
// Probes

std::vector<Probe> parse_probes (const std::string& text, const std::string& file, std::vector<std::string>& err)
{
    std::vector<Probe> out;
    std::istringstream is(text);
    std::string line;
    int lineno = 0;
    auto bad = [&] (const std::string& what) {
        err.emplace_back("unsteady.probes_file = " + file + ", line " + std::to_string(lineno) + ": " + what);
    };
    while (std::getline(is, line)) {
        ++lineno;
        const auto hash = line.find('#');
        if (hash != std::string::npos) { line.erase(hash); }
        std::istringstream ls(line);
        std::vector<std::string> tok;
        for (std::string t; ls >> t;) { tok.push_back(t); }
        if (tok.empty()) { continue; }
        if (tok.size() != 3) {
            bad("expected three numbers \"x y h\", found " + std::to_string(tok.size()) + " fields");
            continue;
        }
        double v[3];
        bool ok = true;
        for (int n = 0; n < 3; ++n) {
            char* end = nullptr;
            v[n] = std::strtod(tok[n].c_str(), &end);
            if (end == tok[n].c_str() || *end != '\0' || !std::isfinite(v[n])) {
                bad("'" + tok[n] + "' is not a finite number");
                ok = false;
            }
        }
        if (!ok) { continue; }
        if (!(v[2] > 0)) {
            bad("h = " + tok[2] + ": the height above the ground must be > 0 (m)");
            continue;
        }
        Probe p;
        p.x = static_cast<amrex::Real>(v[0]);
        p.y = static_cast<amrex::Real>(v[1]);
        p.h = static_cast<amrex::Real>(v[2]);
        out.push_back(p);
    }
    if (out.empty() && err.empty()) { err.emplace_back("unsteady.probes_file = " + file + ": no probe (\"x y h\" lines)"); }
    return out;
}

Probes::Probes (std::vector<Probe> probes, const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom,
                const std::string& probes_file)
    : m_probes(std::move(probes))
{
    const amrex::Box dom = geom.Domain();
    // The column scans take the ground from the box that holds the column's lowest cell: every box
    // must span the column (make_box_array never splits z), as the line solves need.
    const amrex::BoxArray& nba = z_phys_nd.boxArray();
    for (int b = 0; b < nba.size(); ++b) {
        const amrex::Box cb = amrex::enclosedCells(nba[b]);
        if (cb.smallEnd(2) != dom.smallEnd(2) || cb.bigEnd(2) != dom.bigEnd(2)) {
            amrex::Abort("gust: unsteady.probes_file: every box must span the whole column (mesh.max_grid_size is horizontal)");
        }
    }
    for (std::size_t n = 0; n < m_probes.size(); ++n) {
        Probe& p = m_probes[n];
        const std::string which = "probe " + std::to_string(n + 1) + " (" + exact(double(p.x)) + ", " + exact(double(p.y)) + ", "
                                  + exact(double(p.h)) + ")";
        int ij[2];
        const amrex::Real xy[2] = {p.x, p.y};
        for (int d = 0; d < 2; ++d) {
            const amrex::Real lo = geom.ProbLo(d);
            const amrex::Real hi = geom.ProbHi(d);
            if (xy[d] < lo || xy[d] > hi) {
                probe_error(probes_file, which + ": " + (d == 0 ? "x" : "y") + " is outside the domain [" + exact(double(lo)) + ", "
                                         + exact(double(hi)) + "]");
            }
            // A point on the high face takes the last cell.
            ij[d] = std::min(static_cast<int>(std::floor((xy[d] - lo) * geom.InvCellSize(d))) + dom.smallEnd(d), dom.bigEnd(d));
        }
        const amrex::Box col(amrex::IntVect(ij[0], ij[1], dom.smallEnd(2)), amrex::IntVect(ij[0], ij[1], dom.bigEnd(2)));
        const ColumnFacts f = column_facts(z_phys_nd, col, p.h);
        if (p.h > f.top - f.ground) {
            probe_error(probes_file, which + ": h is above the column's top (" + exact(double(f.top - f.ground))
                                     + " m above the ground)");
        }
        const auto [k, hc] = nearest_cell(z_phys_nd, col, p.h, f.dmin);
        if (k == std::numeric_limits<int>::max()) {
            probe_error(probes_file, which + ": no cell of the column is nearest (are the node heights finite?)");
        }
        p.cell = amrex::IntVect(ij[0], ij[1], k);
        p.h_cell = hc;
    }
}

std::vector<amrex::Real> Probes::sample (const State& q) const
{
    std::vector<amrex::Real> out(m_probes.size() * kUnsteadyFields, amrex::Real(0));
    using Sum = amrex::ReduceOpSum;
    for (std::size_t n = 0; n < m_probes.size(); ++n) {
        const amrex::Box cell(m_probes[n].cell, m_probes[n].cell);
        amrex::ReduceOps<Sum, Sum, Sum, Sum, Sum> op;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> data(op);
        using T = decltype(data)::Type;
        for (amrex::MFIter mfi(q.theta); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox() & cell;
            if (!b.ok()) { continue; }
            const auto u = q.u.const_array(mfi);
            const auto v = q.v.const_array(mfi);
            const auto w = q.w.const_array(mfi);
            const auto th = q.theta.const_array(mfi);
            const auto k = q.tke.const_array(mfi);
            op.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept -> T
            {
                constexpr amrex::Real half = amrex::Real(0.5);
                return {half * (u(i, j, kk) + u(i + 1, j, kk)), half * (v(i, j, kk) + v(i, j + 1, kk)),
                        half * (w(i, j, kk) + w(i, j, kk + 1)), th(i, j, kk), k(i, j, kk)};
            });
        }
        const T r = data.value(op);
        amrex::Real* o = out.data() + n * kUnsteadyFields;
        o[0] = amrex::get<0>(r);
        o[1] = amrex::get<1>(r);
        o[2] = amrex::get<2>(r);
        o[3] = amrex::get<3>(r);
        o[4] = amrex::get<4>(r);
    }
    // One cell is on one rank: the sum over the ranks is its value.
    amrex::ParallelAllReduce::Sum(out.data(), static_cast<int>(out.size()), amrex::ParallelDescriptor::Communicator());
    return out;
}

std::string Probes::header () const
{
    std::string h = "# GUST probes: x y h (m, h above the ground) -> the cell (i j k) and its centre's height above the ground (m)\n";
    std::string cols = "t";
    for (std::size_t n = 0; n < m_probes.size(); ++n) {
        const Probe& p = m_probes[n];
        h += "# probe " + std::to_string(n + 1) + ": " + exact(double(p.x)) + " " + exact(double(p.y)) + " " + exact(double(p.h))
             + " -> " + std::to_string(p.cell[0]) + " " + std::to_string(p.cell[1]) + " " + std::to_string(p.cell[2]) + " "
             + exact(double(p.h_cell)) + "\n";
        for (const auto& f : unsteady_field_names()) { cols += ",p" + std::to_string(n + 1) + "_" + f; }
    }
    return h + cols + "\n";
}

bool Probes::open (const std::string& out, const std::string& probes_file, double keep_until) const
{
    int kept_file = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) { kept_file = open_io(out, probes_file, keep_until) ? 1 : 0; }
    amrex::ParallelDescriptor::Bcast(&kept_file, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return kept_file != 0;
}

bool Probes::open_io (const std::string& out, const std::string& probes_file, double keep_until) const
{
    const std::string head = header();
    std::string kept;
    const bool exists = std::filesystem::exists(out);
    if (keep_until >= 0 && !exists) {
        amrex::Print() << "gust: unsteady: no probe file " << out << " to continue; a new one starts at the checkpoint's time\n";
    } else if (keep_until < 0 && exists) {
        amrex::Print() << "gust: unsteady: replaces the probe file " << out << "\n";
    }
    if (keep_until >= 0 && exists) {
        std::ifstream is(out);
        if (!is.good()) { amrex::FileOpenFailed(out); }
        std::string old_head, line;
        // The header: the '#' lines and the column line after them.
        while (std::getline(is, line)) {
            old_head += line + "\n";
            if (line.empty() || line[0] != '#') { break; }
        }
        if (old_head != head) {
            probe_error(probes_file, "the restart's probe file " + out + " was written for other probes (its header differs); "
                                     "move it away or give the probes of the run it continues");
        }
        // The rows up to the checkpoint's time; those after it are the interrupted run's and are run again.
        double last_kept = -1;
        while (std::getline(is, line)) {
            if (line.empty()) { continue; }
            const double t = std::strtod(line.c_str(), nullptr);
            if (t <= keep_until) {
                kept += line + "\n";
                last_kept = t;
            }
        }
        // The header shows the probes, not the run: the last kept row's time lets a user see that
        // the file is the checkpoint's run (its row at the checkpoint's time).
        amrex::Print() << "gust: unsteady: keeps the probe rows of " << out << " up to t = " << exact(last_kept)
                       << " s (the checkpoint is at " << exact(keep_until) << " s)\n";
    }
    std::ofstream os(out, std::ios::trunc);
    if (!os.good()) { amrex::FileOpenFailed(out); }
    os << head << kept;
    os.flush();
    if (!os.good()) { amrex::Abort("gust: could not write " + out); }
    return keep_until >= 0 && exists;
}

void Probes::write_row (const std::string& out, double t, const State& q) const
{
    const std::vector<amrex::Real> s = sample(q);
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::string row = exact(t);
    for (const amrex::Real x : s) { row += "," + exact(double(x)); }
    std::ofstream os(out, std::ios::app);
    if (!os.good()) { amrex::FileOpenFailed(out); }
    os << row << "\n";
    os.flush();
    if (!os.good()) { amrex::Abort("gust: could not write " + out); }
}

// ---------------------------------------------------------------------------------------------
// The checkpoint record

void write_unsteady_record (const std::string& dir, const UnsteadyRecord& r)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream os(dir + "/Unsteady");
    if (!os.good()) { amrex::FileOpenFailed(dir + "/Unsteady"); }
    os << kRecordMagic << " " << kRecordVersion << "\n"
       << "time " << hexfloat(r.time) << "\n"
       << "steps " << r.steps << "\n"
       << "dt " << hexfloat(r.dt) << "\n"
       << "average_start " << hexfloat(r.average_start) << "\n"
       << "weight " << hexfloat(r.weight) << "\n"
       << "anchor_time " << hexfloat(r.anchor_time) << "\n"
       << "anchor_steps " << r.anchor_steps << "\n";
    os.flush();
    if (!os.good()) { amrex::Abort("gust: could not write " + dir + "/Unsteady"); }
}

bool read_unsteady_record (const std::string& dir, UnsteadyRecord& r)
{
    const std::string file = dir + "/Unsteady";
    int exists = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) { exists = std::filesystem::exists(file) ? 1 : 0; }
    amrex::ParallelDescriptor::Bcast(&exists, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (exists == 0) { return false; }
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(file, buf);
    std::istringstream is(std::string(buf.dataPtr()));
    auto fail = [&] (const std::string& what) { amrex::Abort("gust: io.restart = " + dir + ": Unsteady: " + what); };
    std::string magic;
    int version = 0;
    is >> magic >> version;
    if (!is || !is_file_magic(magic, kRecordMagic) || version != kRecordVersion) {
        fail("not an unsteady record of format version " + std::to_string(kRecordVersion));
    }
    auto real = [&] (const char* key) {
        std::string k, tok;
        is >> k >> tok;
        char* end = nullptr;
        const double x = std::strtod(tok.c_str(), &end);
        if (!is || k != key || tok.empty() || *end != '\0' || !std::isfinite(x)) { fail(std::string("bad or missing '") + key + "'"); }
        return x;
    };
    r.time = real("time");
    {
        std::string k;
        is >> k >> r.steps;
        if (!is || k != "steps" || r.steps < 0) { fail("bad or missing 'steps'"); }
    }
    r.dt = real("dt");
    r.average_start = real("average_start");
    r.weight = real("weight");
    r.anchor_time = real("anchor_time");
    {
        std::string k;
        is >> k >> r.anchor_steps;
        if (!is || k != "anchor_steps" || r.anchor_steps < 0) { fail("bad or missing 'anchor_steps'"); }
    }
    if (!(r.time >= 0) || !(r.dt >= 0) || !(r.weight >= 0) || !(r.anchor_time >= 0)) {
        fail("negative time, dt, weight or anchor_time");
    }
    return true;
}

std::string probes_output_file (const std::string& plot_file)
{
    return plot_file + "_probes.csv";
}

} // namespace gust
