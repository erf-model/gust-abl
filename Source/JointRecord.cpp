#include "JointRecord.H"

#include "Sha256.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace gust {

namespace {

constexpr const char* kJointMagic = "GUST_JOINT";
constexpr int kJointVersion = 1;

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
           << "\nouter " << hexfloat_canonical(static_cast<double>(r.outer_cfl_floor)) << "\n";
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
    if (!(is >> magic >> version) || magic != kJointMagic || version != kJointVersion) {
        amrex::Abort(bad + " (expected " + kJointMagic + " " + std::to_string(kJointVersion) + ")");
    }
    if (!(is >> key >> r.feedback) || key != "feedback" || (r.feedback != "defect" && r.feedback != "relax")) {
        amrex::Abort(bad + " (the feedback line)");
    }
    if (!(is >> key >> r.exchanges) || key != "exchanges" || r.exchanges < 1) { amrex::Abort(bad + " (the exchanges line)"); }
    if (!(is >> key >> floor) || key != "outer" || !parse_floor(floor, r.outer_cfl_floor)) {
        amrex::Abort(bad + " (the outer line)");
    }
    while (is >> key) {
        JointRecord::Nest n;
        if (key != "nest" || !(is >> n.name >> floor) || !parse_floor(floor, n.cfl_floor)) {
            amrex::Abort(bad + " (a nest line)");
        }
        std::getline(is, n.chk);
        const auto first = n.chk.find_first_not_of(' ');
        n.chk = first == std::string::npos ? std::string() : n.chk.substr(first);
        if (n.chk.empty()) { amrex::Abort(bad + " (nest " + n.name + " names no checkpoint)"); }
        r.nests.push_back(std::move(n));
    }
    if (r.nests.empty()) { amrex::Abort(bad + " (no nest line)"); }
    return true;
}

} // namespace gust
