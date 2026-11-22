#!/usr/bin/env python3
"""
The one-line mutants of the wall value of k beside buildings
(buildings.wall_k_bound, buildings.wall_k_speed: Source/BuildingTke.cpp,
Source/BuildingForcing.H, Source/BuildingForcing.cpp, Source/Buildings.cpp).

Run against BuildingTke.* and BuildingLaw.* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). The table is in
docs/testing.rst.

    python3 tools/wall_k_mutants.py BUILD_DIR [K1 K4 ...]
"""
import mutants

M = [
 ("K1", "the bound dropped: the law's value alone", "Source/BuildingTke.cpp",
  'return amrex::max(c.k_wall(i, j, k), c.k_impinging(t, i, j, k));',
  'return c.k_wall(i, j, k);'),
 ("K2", "a wall at low x reads the k of the cell at lower x, not the one out along its normal", "Source/BuildingTke.cpp",
  'if ((mask & kWallXlo) != 0) { kr = amrex::max(kr, to(i, j, k, 0) * t(i + 1, j, k)); }',
  'if ((mask & kWallXlo) != 0) { kr = amrex::max(kr, to(i, j, k, 0) * t(i - 1, j, k)); }'),
 ("K3", "a roof takes the share of the walls along y", "Source/BuildingTke.cpp",
  'if ((mask & kWallRoof) != 0 && k < khi) { kr = amrex::max(kr, to(i, j, k, 2) * t(i, j, k + 1)); }',
  'if ((mask & kWallRoof) != 0 && k < khi) { kr = amrex::max(kr, to(i, j, k, 1) * t(i, j, k + 1)); }'),
 ("K4", "a fixed face: half of the open face's square in place of a third", "Source/BuildingForcing.H",
  '? s * s : s * s / amrex::Real(3);',
  '? s * s : s * s / amrex::Real(2);'),
 ("K5", "the rows of k read the law's u*^2 with cell_mean", "Source/BuildingForcing.H",
  'return m_k_cell_mean ? kBuildingUstar2K : kBuildingUstar2; }',
  'return kBuildingUstar2; }'),
 ("K6", "no line in the restart hash for wall_k_bound = impinging", "Source/Buildings.cpp",
  '? "wall_k_bound impinging\\n" : "")',
  '? "" : "")'),
 ("K15", "no line in the restart hash for wall_k_speed = cell_mean", "Source/Buildings.cpp",
  '? "wall_k_speed cell_mean\\n" : "")',
  '? "" : "")'),
 ("K7", "cell_mean accepted with wall_speed = all_faces", "Source/Buildings.cpp",
  'opt.wall_k_speed == BuildingWallKSpeed::cell_mean && opt.wall_speed != BuildingWallSpeed::open_faces',
  'false && opt.wall_speed != BuildingWallSpeed::open_faces'),
 ("K8", "the side at high y takes the open face's speed with cell_mean", "Source/BuildingForcing.cpp",
  'if ((mask & kWallYhi) != 0) { top = amrex::max(top, form7_ustar2(own, std::sqrt(u2(i, j - 1, k) + w2(i, j - 1, k)), b0, d1, d2, z0, kappa)); }',
  'if ((mask & kWallYhi) != 0) { top = amrex::max(top, form7_ustar2(std::sqrt(amrex::Real(3)) * own, std::sqrt(u2(i, j - 1, k) + w2(i, j - 1, k)), b0, d1, d2, z0, kappa)); }'),
 ("K9", "a flow that leaves the wall counts as one towards it", "Source/BuildingForcing.H",
  'const amrex::Real towards = amrex::max(n, amrex::Real(0));',
  'const amrex::Real towards = n;'),
 ("K10", "the share is of the speed, not of the kinetic energy", "Source/BuildingForcing.H",
  'return towards * towards / (n * n + t1 * t1 + t2 * t2 + kBuildingMinSpeed * kBuildingMinSpeed);',
  'return towards / std::sqrt(n * n + t1 * t1 + t2 * t2 + kBuildingMinSpeed * kBuildingMinSpeed);'),
 ("K11", "a side at high y takes towards as the positive y of a side at low y", "Source/BuildingForcing.cpp",
  'if ((mask & kWallYhi) != 0) { sy = amrex::max(sy, towards_wall_share(vc(i, j - 1, k), uc(i, j - 1, k), wc(i, j - 1, k))); }',
  'if ((mask & kWallYhi) != 0) { sy = amrex::max(sy, towards_wall_share(-vc(i, j - 1, k), uc(i, j - 1, k), wc(i, j - 1, k))); }'),
 ("K12", "the roof's share reads the wall cell, not the cell above", "Source/BuildingForcing.cpp",
  'sz = amrex::max(sz, towards_wall_share(-wc(i, j, k + 1), uc(i, j, k + 1), vc(i, j, k + 1))); }',
  'sz = amrex::max(sz, towards_wall_share(-wc(i, j, k), uc(i, j, k), vc(i, j, k))); }'),
 ("K14", "the share of a side at low x reads the wall cell's own wind", "Source/BuildingForcing.cpp",
  'sx = amrex::max(sx, towards_wall_share(-uc(i + 1, j, k), vc(i + 1, j, k), wc(i + 1, j, k))); }',
  'sx = amrex::max(sx, towards_wall_share(-uc(i, j, k), vc(i, j, k), wc(i, j, k))); }'),
 ("K13", "the shares are not filled", "Source/BuildingForcing.cpp",
  'if (m_k_towards) { towards_wall_shares(*m_b, m_geom, u, v, w, m_rate); }',
  'if (false) { towards_wall_shares(*m_b, m_geom, u, v, w, m_rate); }'),
]

if __name__ == "__main__":
    mutants.run(M, "BuildingTke.*:BuildingLaw.*", set(), __doc__)
