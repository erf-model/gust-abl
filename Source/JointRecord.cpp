#include "JointRecord.H"

#include "Sha256.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace gust {

namespace {

constexpr const char* kJointMagic = "GUST_JOINT";
constexpr int kJointVersion = 3;
constexpr const char* kFeedbackNames[5] = {"u", "v", "w", "theta", "tke"};

std::string feedback_file (const std::string& dir, int n)
{
    return amrex::MultiFabFileFullPrefix(0, dir, "Level_", std::string("joint_feedback_") + kFeedbackNames[n]);
}

std::string joint_file (const std::string& dir) { return dir + "/Joint"; }

//! A finite hex float > 0 or exactly 0 (a floor not set), or false.
bool parse_floor (const std::string& s, amrex::Real& out)
{
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0' || !std::isfinite(v) || v < 0) { return false; }
    out = static_cast<amrex::Real>(v);
    return true;
}

} // namespace

void write_joint_record (const std::string& dir, const JointRecord& r)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string name = joint_file(dir);
        std::ofstream os(name);
        if (!os.good()) { amrex::FileOpenFailed(name); }
        os << kJointMagic << " " << kJointVersion << "\nfeedback " << r.feedback << "\nexchanges " << r.exchanges
           << "\nouter " << hexfloat_canonical(static_cast<double>(r.outer_cfl_floor)) << "\nwindow " << r.window.size();
        for (const amrex::Real d : r.window) { os << " " << hexfloat_canonical(static_cast<double>(d)); }
        os << "\nfactor " << hexfloat_canonical(static_cast<double>(r.factor)) << "\n";
        for (const auto& n : r.nests) {
            os << "nest " << n.name << " " << hexfloat_canonical(static_cast<double>(n.cfl_floor)) << " " << n.chk << "\n";
        }
        os.close();
        if (!os.good()) { amrex::Abort("gust: cannot write " + name); }
    }
    amrex::ParallelDescriptor::Barrier();
}

bool read_joint_record (const std::string& dir, JointRecord& r)
{
    // The file's text on the I/O rank ("" without a file), then on every rank.
    std::string text;
    int found = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(joint_file(dir));
        if (is) {
            found = 1;
            std::ostringstream ss;
            ss << is.rdbuf();
            text = ss.str();
        }
    }
    amrex::ParallelDescriptor::Bcast(&found, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (found == 0) { return false; }
    amrex::Vector<char> buf(text.begin(), text.end());
    amrex::BroadcastArray(buf, amrex::ParallelDescriptor::MyProc(), amrex::ParallelDescriptor::IOProcessorNumber(),
                          amrex::ParallelDescriptor::Communicator());
    text.assign(buf.begin(), buf.end());

    const std::string bad = "gust: io.restart = " + dir + ": its Joint file cannot be read";
    std::istringstream is(text);
    std::string magic, key, floor;
    int version = 0;
    r = JointRecord{};
    if (!(is >> magic >> version) || magic != kJointMagic || version < 1 || version > kJointVersion) {
        amrex::Abort(bad + " (expected " + kJointMagic + " 1 to " + std::to_string(kJointVersion) + ")");
    }
    if (!(is >> key >> r.feedback) || key != "feedback"
        || (r.feedback != "defect" && r.feedback != "relax" && r.feedback != "relax_consistent" && r.feedback != "inject")) {
        amrex::Abort(bad + " (the feedback line)");
    }
    if (!(is >> key >> r.exchanges) || key != "exchanges" || r.exchanges < 1) { amrex::Abort(bad + " (the exchanges line)"); }
    if (!(is >> key >> floor) || key != "outer" || !parse_floor(floor, r.outer_cfl_floor)) {
        amrex::Abort(bad + " (the outer line)");
    }
    if (version >= 2) {
        // The window holds at most steady.monitor_window values; the cap only keeps a corrupt count from allocating.
        long m = -1;
        if (!(is >> key >> m) || key != "window" || m < 0 || m > 1000000) { amrex::Abort(bad + " (the window line)"); }
        r.window.resize(static_cast<std::size_t>(m));
        for (auto& d : r.window) {
            if (!(is >> floor) || !parse_floor(floor, d)) { amrex::Abort(bad + " (the window line)"); }
        }
        r.has_feedback = true;
    }
    if (version >= 3) {
        // 0 with defect (no factor), > 0 with the others.
        if (!(is >> key >> floor) || key != "factor" || !parse_floor(floor, r.factor)
            || ((r.feedback == "defect") != (r.factor == 0))) {
            amrex::Abort(bad + " (the factor line)");
        }
        r.has_factor = true;
    }
    while (is >> key) {
        JointRecord::Nest n;
        if (key != "nest" || !(is >> n.name >> floor) || !parse_floor(floor, n.cfl_floor)) {
            amrex::Abort(bad + " (a nest line)");
        }
        std::getline(is, n.chk);
        const auto first = n.chk.find_first_not_of(' ');
        n.chk = first == std::string::npos ? std::string() : n.chk.substr(first);
        // A record edited on another system may end its lines with \r or spaces.
        n.chk.erase(n.chk.find_last_not_of(" \t\r") + 1);
        if (n.chk.empty()) { amrex::Abort(bad + " (nest " + n.name + " names no checkpoint)"); }
        r.nests.push_back(std::move(n));
    }
    if (r.nests.empty()) { amrex::Abort(bad + " (no nest line)"); }
    return true;
}

std::string joint_restore_refusal (const JointRecord& prior, const std::string& feedback, amrex::Real factor)
{
    if (!prior.has_feedback) { return "its joint record (GUST_JOINT 1) holds no feedback"; }
    if (prior.feedback != feedback) { return "its joint record holds the feedback of nest.feedback = " + prior.feedback; }
    if (prior.has_factor && prior.factor != factor) {
        const auto text = [] (amrex::Real v, int digits) {
            std::ostringstream os;
            os << std::setprecision(digits) << v;
            return os.str();
        };
        // The digits a deck is written with, unless the two values then print alike.
        using Limits = std::numeric_limits<amrex::Real>;
        const int digits = text(prior.factor, Limits::digits10) == text(factor, Limits::digits10) ? Limits::max_digits10
                                                                                                  : Limits::digits10;
        return "its joint record holds the feedback of nest.feedback_factor = " + text(prior.factor, digits) + ", not "
               + text(factor, digits);
    }
    return {};
}

void write_joint_feedback (const std::string& dir, const std::array<const amrex::MultiFab*, 5>& f)
{
    for (int n = 0; n < 5; ++n) { amrex::VisMF::Write(*f[static_cast<std::size_t>(n)], feedback_file(dir, n)); }
}

void read_joint_feedback (const std::string& dir, const std::array<const amrex::MultiFab*, 5>& like,
                          std::array<amrex::MultiFab, 5>& out)
{
    for (int n = 0; n < 5; ++n) {
        const auto i = static_cast<std::size_t>(n);
        const std::string file = feedback_file(dir, n);
        const std::string what = "gust: io.restart = " + dir + ": its joint record has a window, but the feedback field "
                                 + std::string("joint_feedback_") + kFeedbackNames[n];
        if (!amrex::VisMF::Exist(file)) { amrex::Abort(what + " is missing"); }
        amrex::MultiFab src;   // defined from the file's own layout
        amrex::VisMF::Read(src, file);
        if (src.nComp() != 1 || src.ixType() != like[i]->ixType()) {
            amrex::Abort(what + " has another index type or component count");
        }
        out[i].define(like[i]->boxArray(), like[i]->DistributionMap(), 1, 0);
        out[i].setVal(0.0);
        out[i].ParallelCopy(src, 0, 0, 1, 0, 0);
    }
}

} // namespace gust
