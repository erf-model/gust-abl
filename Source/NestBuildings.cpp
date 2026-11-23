#include "NestBuildings.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace gust {

std::vector<std::string> check_nest_buildings (const BuildingHeightMap& map, const std::string& file,
                                               const BoundaryInputs& bin, const amrex::Geometry& geom,
                                               amrex::Real blend_width, int max_messages)
{
    std::vector<std::string> err;
    int count = 0;
    const char* const names[2][2] = {{"bc.xlo", "bc.xhi"}, {"bc.ylo", "bc.yhi"}};
    for (int d = 0; d < 2; ++d) {
        const int o = 1 - d;   // the direction along the face
        const double dx = static_cast<double>(geom.CellSize(d));
        // The wall cells and their reference cells (two cells) stay outside the zone and the blend.
        const double gap = kNestBuildingCells * dx
                           + std::max(static_cast<double>(bin.nest_relax_width) * dx, static_cast<double>(blend_width));
        const double lo = static_cast<double>(geom.ProbLo(d)), hi = static_cast<double>(geom.ProbHi(d));
        const double olo = static_cast<double>(geom.ProbLo(o)), ohi = static_cast<double>(geom.ProbHi(o));
        // The geometry is in Real: in single precision its round-off is far above 1e-9 of a cell
        // (the form of the allowance of Buildings::build_from_map, which takes the domain's length).
        const double eps_real = 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
        const double tol = std::max(1e-9 * dx, eps_real * std::max(std::abs(lo), std::abs(hi)));
        const double otol = std::max(1e-9 * static_cast<double>(geom.CellSize(o)),
                                      eps_real * std::max(std::abs(olo), std::abs(ohi)));
        const bool operiodic = geom.isPeriodic(o);
        const std::vector<double>& e = d == 0 ? map.xe : map.ye;
        const std::vector<double>& eo = d == 0 ? map.ye : map.xe;
        for (int side = 0; side < 2; ++side) {
            if (bin.face[amrex::Orientation(d, side == 0 ? amrex::Orientation::low : amrex::Orientation::high)] != FaceType::nest) {
                continue;
            }
            for (int i = 0; i < map.nx; ++i) {
                for (int j = 0; j < map.ny; ++j) {
                    const double h = map.h[static_cast<std::size_t>(i) * static_cast<std::size_t>(map.ny) + static_cast<std::size_t>(j)];
                    if (!(h > 0)) { continue; }
                    const int a = d == 0 ? i : j, b = d == 0 ? j : i;
                    const double e0 = e[static_cast<std::size_t>(a)], e1 = e[static_cast<std::size_t>(a) + 1];
                    const double f0 = eo[static_cast<std::size_t>(b)], f1 = eo[static_cast<std::size_t>(b) + 1];
                    if (!operiodic && (f1 <= olo + otol || f0 >= ohi - otol)) { continue; }
                    const bool near = side == 0 ? (e1 > lo + tol && e0 < lo + gap - tol)
                                                : (e0 < hi - tol && e1 > hi - gap + tol);
                    if (!near) { continue; }
                    ++count;
                    if (count > max_messages) { continue; }
                    const double x0 = d == 0 ? e0 : f0, x1 = d == 0 ? e1 : f1;
                    const double y0 = d == 0 ? f0 : e0, y1 = d == 0 ? f1 : e1;
                    char line[512];
                    std::snprintf(line, sizeof(line),
                                  "buildings.file = %s: the building tile x in [%.6g, %.6g] m, y in [%.6g, %.6g] m (top %.6g m) "
                                  "lies within %.6g m of the nest face %s (%d cells beyond the larger of bc.nest_relax_width = %d "
                                  "cells and terrain.nest_blend_width = %.6g m); move the face or the building",
                                  file.c_str(), x0, x1, y0, y1, h, gap, names[d][side], kNestBuildingCells,
                                  bin.nest_relax_width, static_cast<double>(blend_width));
                    err.emplace_back(line);
                }
            }
        }
    }
    if (count > max_messages) {
        err.push_back("buildings.file = " + file + ": " + std::to_string(count - max_messages)
                      + " more pairs of a building tile and a nest face lie within the gap (a corner tile counts once per face)");
    }
    return err;
}

} // namespace gust
