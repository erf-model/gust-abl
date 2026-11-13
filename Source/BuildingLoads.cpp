#include "BuildingLoads.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>

namespace gust {

namespace {

constexpr const char* kHeader = "t,Fx,Fy,Fz";

std::string exact (double x)
{
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", x);
    return buf;
}

} // namespace

std::array<amrex::Real, 3> building_force (const Buildings& b, const Mesh& mesh, const amrex::Geometry& geom,
                                           const ResidualFields& R)
{
    const amrex::MultiFab* Rc[3] = {&R.u, &R.v, &R.w};
    const amrex::Box& dom = geom.Domain();
    const auto dx = geom.CellSizeArray();
    const amrex::Real cell = dx[0] * dx[1] * dx[2];
    std::array<amrex::Real, 3> F{0, 0, 0};
    for (int c = 0; c < 3; ++c) {
        // The c-faces that are rows of the domain's boundary conditions: the two domain faces normal to c
        // (unless periodic), and w at the ground (BuildingForcing's ground_row; inside the first range).
        const bool periodic = geom.isPeriodic(c);
        const int lo = dom.smallEnd(c), hi = dom.bigEnd(c) + 1;
        const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0, dk = (c == 2) ? 1 : 0;
        amrex::MultiFab f(Rc[c]->boxArray(), Rc[c]->DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto out = f.array(mfi);
            const auto r = Rc[c]->const_array(mfi);
            const auto beta = b.beta_face(c).const_array(mfi);
            const auto J = mesh.detJ.const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int n = (c == 0) ? i : (c == 1) ? j : k;
                const bool row = !periodic && (n == lo || n == hi);
                out(i, j, k) = (beta(i, j, k) == amrex::Real(1) && !row)
                                   ? -r(i, j, k) * amrex::Real(0.5) * (J(i - di, j - dj, k - dk) + J(i, j, k)) * cell
                                   : amrex::Real(0);
            });
        }
        // A face on a box boundary (or a periodic image) is in two boxes: count it once.
        F[static_cast<std::size_t>(c)] = f.sum_unique(0, false, geom.periodicity());
    }
    return F;
}

bool BuildingLoads::open (double keep_until) const
{
    int kept_file = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) { kept_file = open_io(keep_until) ? 1 : 0; }
    amrex::ParallelDescriptor::Bcast(&kept_file, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return kept_file != 0;
}

bool BuildingLoads::open_io (double keep_until) const
{
    std::string kept;
    const bool exists = std::filesystem::exists(m_file);
    if (keep_until >= 0 && exists) {
        std::ifstream is(m_file);
        if (!is.good()) { amrex::FileOpenFailed(m_file); }
        std::string line;
        std::getline(is, line);
        if (line != kHeader) {
            amrex::Abort("gust: buildings.force_file = " + m_file + ": the restart's file has another header (" + line + ")");
        }
        // The rows up to the checkpoint's time; those after it are run again.
        while (std::getline(is, line)) {
            if (!line.empty() && std::strtod(line.c_str(), nullptr) <= keep_until) { kept += line + "\n"; }
        }
        amrex::Print() << "gust: buildings: keeps the force rows of " << m_file << " up to t = " << exact(keep_until) << " s\n";
    } else if (keep_until >= 0) {
        amrex::Print() << "gust: buildings: the restart finds no force file " << m_file << "; starts a new one\n";
    } else if (exists) {
        amrex::Print() << "gust: buildings: replaces the force file " << m_file << "\n";
    }
    std::ofstream os(m_file, std::ios::trunc);
    if (!os.good()) { amrex::FileOpenFailed(m_file); }
    os << kHeader << "\n" << kept;
    return keep_until >= 0 && exists;
}

std::array<amrex::Real, 3> BuildingLoads::evaluate (const State& q, ResidualInputs in, const Mesh& mesh,
                                                    const amrex::Geometry& geom, const Buildings& b)
{
    in.building_rows = false;
    in.time = nullptr;
    residual(q, in, mesh, geom, m_R, nullptr, nullptr, &m_ws);
    return building_force(b, mesh, geom, m_R);
}

void BuildingLoads::write_row (double t, const std::array<amrex::Real, 3>& F) const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream os(m_file, std::ios::app);
    if (!os.good()) { amrex::FileOpenFailed(m_file); }
    os << exact(t) << "," << exact(double(F[0])) << "," << exact(double(F[1])) << "," << exact(double(F[2])) << "\n";
}

} // namespace gust
