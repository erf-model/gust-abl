# Nest_Shapes

## Purpose
One-way nested patches around six canonical shapes, against a uniform fine run of the whole domain (tracker row
PR-76; sub-phase 9.4 in the repository's numbering, 7.4 in the tracker's; beyond plan rev. 6). Each patch is at
refinement 2 2 1 (4 4 1 in the ratio runs below) of a coarse parent run, with the same z-levels (for the solar
farm, a chain of two such nests, and 4 4 1 beside it), and its `nest` faces lie at a distance D from the shape's
centre. The reference is the whole domain on the patch's cells. Three terrain shapes test patches over steep ground
(the ridge's lee separates); a cube, a raised plate and a solar farm (rows of flat raised panels, tracker row PR-68)
test buildings inside a patch (`buildings.type` with `nest` faces, allowed when every building keeps clear of the
faces, `docs/boundary_conditions.rst`). `Tests/Canonical_RANS/Nest_Validation` is the same study on a
gentler hill; `Tests/Canonical_RANS/Nest_Faces` tests the faces themselves.

## Files
| File | Role |
|---|---|
| `inputs_terrain_domain` | The whole domain for the terrain shapes: 3200 m long, periodic in y, 800 m high, 16 layers from 10 m, log-law inflow along x (u* 0.4 m/s, z0 0.1 m), neutral, SIMPLEC to the steady verdict. |
| `inputs_terrain_nest` | A patch over a terrain shape: `nest` faces, first guess from the parent. |
| `inputs_building_domain` | The whole domain for the cube, the plate and the solar farm: flat, periodic in y, log-law inflow along x, the building with the form 7 wall law, z-line sweeps. |
| `inputs_building_nest` | A patch around a building. |
| `nest_shapes.py` | Writes each shape's files (`files`), runs a whole domain (`run`) or a patch (`nest`), compares a patch with the reference (`compare`), or does the whole study (`study`). Pure Python. |

The decks hold the recommended setting (`turb.realizability` 0.309, `turb.nut_relax` 0.5,
`convection.hybrid_sigma_max` 0.1, `turb.max_geom_lscale` 100 m). The script gives every run its shape: the terrain
or building (`buildings.file`, `buildings.format`), the domain size and cells, the roughness (`bc.surface_z0`,
`bc.inflow_z0`, `buildings.z0`) and the closure's `convection.hybrid_height` (2.5 shape heights). For the buildings
it also gives `turb.max_geom_lscale` = 40 m (the value of `Tests/Canonical_RANS/Building_Cube`), the CFL cap and the
iteration limit: `steady.cfl_max` 8 and 4000 iterations for the cube; `steady.cfl_max` 4, 20000 iterations and
`steady.confidence` = 0.95 for the plate and the solar farm, with `mesh.zlevels_file`. On the plate's 1 m cells the
cap of 8 took 10023 iterations to the steady verdict against 5055 with 4. The solar farm also gets
`simplec.pressure_reltol` = 1e-6 (default 1e-3); Expected Results gives the reason.

| Shape | Domain (m) | Shape | Parent / reference cells | D (m) | U_ref (m/s) |
|---|---|---|---|---|---|
| ridge | 3200 x 200 x 800, 16 layers | Witch of Agnesi ridge along y, h = 150 m, half-width a = 150 m (largest slope 33 degrees) | 50 / 25 m | 300, 600, 900 | 7.313 |
| hill | 3200 x 2400 x 800, 16 layers | round Witch of Agnesi hill, h = 150 m, a = 150 m | 50 / 25 m | 300, 600, 900 | 7.313 |
| escarpment | 3200 x 200 x 800, 16 layers | z = (H/2)(1 + tanh((x - 1600 m)/L)), H = 100 m, L = 100 m, from an ERF terrain file (`terrain.type = erf_grid`) | 50 / 25 m | 300, 600, 900 | 6.908 |
| cube | 480 x 240 x 120, 24 layers of 5 m | the block of `Building_Cube`, 30 x 30 x 20 m, centre (130, 120) m, z0 0.1 m | 10 / 5 m | 40, 70, 100 | 5.298 |
| plate | 64 x 4 x 20.6, 93 layers (0.025 m to 1.5 m, then x 1.15) | a plate 4 m long, 0.1 m thick, 1 m above the ground, across the periodic width, z0 0.01 m | 1 / 0.5 m | 4, 8, 12 | 4.700 |
| solar | 160 x 2 x 63.9, 219 layers (0.02 m to 3.2 m, then x 1.1) | six rows of flat panels, each 4 m long in x and 0.1 m thick (five layers), underside 1.5 m above the ground, across the periodic width (rows without ends), 10 m apart, centres at x = 55 to 105 m; z0 0.01 m | 1 m, middle patch 0.5 m (D = 64 m) / 0.25 m | 32, 40, 48 | 5.075 |

The patch is a square of half-width D around the hill and the cube (four `nest` faces), and spans the periodic width
around the ridge, the escarpment, the plate and the solar farm (two `nest` faces, periodic in y). D is measured from
the shape's centre (x = 80 m for the solar farm): the nearest building wall is D - 15 m from a face for the cube,
D - 2 m for the plate and D - 27 m for the solar farm.

The solar farm is a chain of two nests at 2 2 1: a 0.5 m middle patch of the 1 m parent (D = 64 m), and 0.25 m
patches of the middle one at each D. The targets judge a 0.25 m patch against its own parent, the middle patch; the
1 m parent is reported too. Beside the chain, the study runs 0.25 m patches straight from the 1 m parent at 4 4 1
(`bc.nest_ratio = 4`, runs `solar_r4_dD`) at the same D, for the record, with no target. The solar farm's 0.02 m
layers reach twice the panel top (3.2 m), so each panel is five layers thick. The panels are a stand-in for the
tilted panels of the PV phase (tracker row PR-52). The middle patch's faces are 16 m from the domain's inflow and
outflow faces. The solar runs set `simplec.pressure_reltol` = 1e-6; see Expected Results.

## Data and acknowledgements
None. Every shape is analytic; the escarpment's terrain file, the building maps and the plate's and solar farm's
z-levels are written by `nest_shapes.py files SHAPE`.

## How to run
In a scratch directory:
```bash
python3 nest_shapes.py study -- mpiexec -n 4 <build>/gust
python3 nest_shapes.py study --shapes cube --ranks cube -- mpiexec -n 4 <build>/gust
```
Single runs:
```bash
python3 nest_shapes.py run cube cube_ref 5 -- mpiexec -n 4 <build>/gust inputs_building_domain
python3 nest_shapes.py run cube cube_par 10 -- mpiexec -n 4 <build>/gust inputs_building_domain
python3 nest_shapes.py nest cube cube_d70 cube_par 70 -- mpiexec -n 4 <build>/gust inputs_building_nest
python3 nest_shapes.py compare cube cube_d70 cube_ref cube_par
```
A study without `--shapes` runs every shape but the solar farm. There is no time limit on a run: the solver's output
goes to `NAME.log` as it runs, and `steady.max_iterations` ends a run that does not converge. The solar farm's runs
of the table below (the parent and middle patch on three ranks, the patches on two, the reference on eight):
```bash
python3 nest_shapes.py run solar solar_par 1 -- mpiexec -n 3 <build>/gust inputs_building_domain
python3 nest_shapes.py run solar solar_ref 0.25 -- mpiexec -n 8 <build>/gust inputs_building_domain
python3 nest_shapes.py nest solar solar_mid solar_par 64 -- mpiexec -n 3 <build>/gust inputs_building_nest
python3 nest_shapes.py nest solar solar_d40 solar_mid 40 -- mpiexec -n 2 <build>/gust inputs_building_nest
python3 nest_shapes.py nest solar solar_r4_d40 solar_par 40 --nest-ratio 4 -- mpiexec -n 2 <build>/gust inputs_building_nest
python3 nest_shapes.py nest solar solar_np4 solar_mid 40 -- mpiexec -n 4 <build>/gust inputs_building_nest "mesh.max_grid_size=160 4"
python3 nest_shapes.py compare solar solar_d40 solar_ref solar_mid
```
(D = 32 and 48 m likewise; `solar_np1` and `solar_np2` with `-n 1` and `-n 2`.) `study --shapes solar --ranks solar`
runs all of them with one launch command.

## Measures
The difference between two runs on the same cells is the length of the velocity-vector difference as a fraction of the
shape's speed scale U_ref = (u*/0.4) ln(h/z0), h the shape's height (the building's top for the cube and the plate,
the panels' top, 1.6 m, for the solar farm),
and the difference of k as a fraction of U_ref^2. Both are given as RMS and largest value over four regions:

| Region | Cells |
|---|---|
| roi | around the shape: within a of the crest (ridge, hill), within 200 m of the escarpment's middle, x in [100, 160] m and y in [90, 150] m (cube), x in [19, 25] m (plate), x in [53, 107] m (solar, the array from the first row's leading edge to the last row's trailing edge); layers centred below the shape's height (ridge, hill, escarpment), 2 cube heights (cube), 3 m (plate) or 3.2 m (solar) over flat ground |
| lee | downwind: x from a to 3 a past the crest (ridge, hill; within a of the axis for the hill), 200 to 600 m past the escarpment's middle, x in (145, 165] m behind the cube within its width, x in (24, 28] m behind the plate, x in (107, 111] m behind the solar farm; the same layers |
| faces | the two rings of patch cells beside the `nest` faces, all layers |
| patch | every patch cell |

The parent is compared with the reference averaged over each parent cell, over the same regions. Its ground (and its
buildings) are on the coarser cells, so its difference includes them. Patch and reference have the same ground and
buildings on the same cells (`test_NestBuildings.cpp` checks the buildings' fractions and wall cells).

## Targets
Fixed in this file while the first cube and ridge studies were running, before any of their numbers were read
(commit "Nest_Shapes: case folder and targets").

| | Target |
|---|---|
| T1 | Ridge, hill, escarpment, cube and solar farm, at the largest D: ROI velocity RMS <= 0.5 % of U_ref and <= 1/4 of the parent's |
| T2 | Ridge, hill, escarpment, cube and solar farm: ROI velocity RMS does not grow with D (5 % allowance) |
| T3 | With `--ranks`: the middle D in four boxes on 1, 2 and 4 ranks reaches the steady verdict, and the velocity fields of 2 and 4 ranks differ from that of 1 rank by at most 1e-8 of U_ref anywhere in the patch |
| T4 | Every run reaches the steady verdict |
| — | Plate: no accuracy target (tracker row PR-76). Reported: the four regions and the iterations of each patch and of the reference |

The solar farm was added after the first five shapes (tracker row PR-68). Its T1 and T2 were applied after one trial
study had been read, on another mesh: one nest (1 m parent, 0.5 m patches and reference), 72 layers, a 4 m periodic
width, `simplec.pressure_reltol` = 1e-3, the parent on four ranks (ROI RMS 0.269, 0.121 and 0.065 % at D = 32, 40
and 48 m, parent 1.438 %; 564 iterations for the parent). Unlike the other shapes, the targets were not fixed before
its numbers were known.

The iterations of each patch against those of the reference are reported (PT-9.4.1), with no target. A target that
is missed is reported as missed, with the number.

T3 first asked the three ROI RMS values to agree to 10 %. The review of the PR found that too loose to catch an
error of the decomposition (the measured differences were round-off), and T3 was changed to the field difference
above after the cube, hill and plate rank runs had been read (largest differences 2.3e-10, 1.9e-10 and 2.1e-15 of
U_ref).

The regions compare cells with the same indices. Over terrain the parent's ground is sampled on its own cells, so a
parent layer and the averaged reference layer are not at quite the same height; the parent's difference includes
that offset, as it includes its coarser buildings.

## Expected Results
macOS arm64, double precision, the study on four ranks (ridge and escarpment on two), from 2026-11-15; the solar
farm on 2026-11-16 and 17 (parent and middle patch on three ranks, patches on two, reference on eight). Velocity RMS
of the difference from the reference, % of U_ref; iterations to the steady verdict.

| Shape | Run | Cells | Iterations | ROI | Lee | Faces |
|---|---|---|---|---|---|---|
| cube | reference, 5 m | 110 592 | 131 | | | |
| | parent, 10 m | 27 648 | 154 | 18.431 | 9.164 | |
| | D = 40 m | 6 144 | 108 | 4.599 | 5.149 | 2.521 |
| | D = 70 m | 18 816 | 110 | 1.394 | 1.431 | 1.004 |
| | **D = 100 m** | 38 400 | 110 | **0.454** | 0.505 | 0.659 |
| hill | reference, 25 m | 196 608 | 106 | | | |
| | parent, 50 m | 49 152 | 149 | 3.392 | 5.660 | |
| | D = 300 m | 9 216 | 118 | 0.571 | 3.569 | 1.405 |
| | D = 600 m | 36 864 | 165 | 0.134 | 0.614 | 0.856 |
| | **D = 900 m** | 82 944 | 165 | **0.043** | 0.101 | 0.625 |
| escarpment | reference, 25 m | 16 384 | 129 | | | |
| | parent, 50 m | 4 096 | 136 | 0.530 | 0.112 | |
| | D = 300 m | 3 072 | 112 | 0.019 | 0.053 | 0.052 |
| | D = 600 m | 6 144 | 111 | 0.005 | 0.033 | 0.035 |
| | **D = 900 m** | 9 216 | 111 | **0.004** | 0.013 | 0.036 |
| plate | reference, 0.5 m | 95 232 | 2235 | | | |
| | parent, 1 m | 23 808 | 5055 | 2.530 | 1.551 | |
| | D = 4 m | 11 904 | 1951 | 0.491 | 0.674 | 0.560 |
| | D = 8 m | 23 808 | 2412 | 0.180 | 0.255 | 0.367 |
| | D = 12 m | 35 712 | 2241 | 0.082 | 0.110 | 0.324 |
| solar | reference, 0.25 m | 1 121 280 | 2574 | | | |
| | parent, 1 m | 70 080 | 4671 | 1.845 | 1.571 | |
| | middle patch, 0.5 m, D = 64 m | 224 256 | 2914 | 0.910 | 0.724 | 0.119 |
| | D = 32 m | 448 512 | 1350 | 0.092 | 0.249 | 0.255 |
| | D = 40 m | 560 640 | 1435 | 0.036 | 0.088 | 0.200 |
| | **D = 48 m** | 672 768 | 1511 | **0.011** | 0.026 | 0.144 |
| | 4 4 1, D = 32 m | 448 512 | 1565 | 0.196 | 0.513 | 0.600 |
| | 4 4 1, D = 40 m | 560 640 | 1588 | 0.082 | 0.169 | 0.381 |
| | 4 4 1, D = 48 m | 672 768 | 1668 | 0.042 | 0.077 | 0.206 |
| ridge | reference, 25 m | 16 384 | 2333 | | | |
| | parent, 50 m | 4 096 | 677 | 8.499 | 13.895 | |
| | D = 300 m | 3 072 | 494 | 5.652 | 14.381 | 7.800 |
| | D = 600 m | 6 144 | 2309 | 8.144 | 20.171 | 8.232 |
| | D = 900 m | 9 216 | 1488 | 8.606 | 21.380 | 7.770 |

The parent's lee and faces regions differ a little between rows (the faces are the patch's); the table gives the
largest D's. The solar farm's targets judge each 0.25 m patch against the middle patch, whose own row gives its
regions (its faces at x = 16 and 144 m); around the panels it is 16 % of U_ref off the reference at most (the 1 m
parent 27 %), the 0.25 m patches 0.25 % or less. The walls of the cube are D - 15 m from the faces: within 0.5 % of
U_ref around the cube took 85 m (2.8 cube widths).

| Target | Verdict |
|---|---|
| T1 | Met for the cube (0.454 %, 0.025 of the parent's), the hill (0.043 %, 0.013) and the escarpment (0.004 %, 0.008). **Missed for the ridge** (8.606 %, 1.01 of the parent's). |
| T2 | Met for the cube, the hill and the escarpment. **Missed for the ridge** (5.652, 8.144, 8.606 %). |
| T1, solar | Met: 0.011 % at D = 48 m, 0.012 of the middle patch's (0.910 %). |
| T2, solar | Met: 0.092, 0.036 and 0.011 % at D = 32, 40 and 48 m. |
| T3 | Met for the cube (D = 70 m), the hill (D = 600 m), the plate (D = 8 m) and the solar farm (D = 40 m): 2.3e-10, 1.9e-10, 2.1e-15 and 1.4e-15 of U_ref; the same iterations on 1, 2 and 4 ranks (110, 161, 2412, 1435). |
| T4 | Met: every run above reached the steady verdict. The solar farm's reference did so in a second run; the first, on four ranks of a shared laptop, was stopped at iteration 2550 (confidence 0.947) by the 8 h limit per run that the script had then. |

The patches take about as many iterations as the reference, or fewer (PT-9.4.1): 108–110 against 131 (cube), 118–165 against
106 (hill), 111–112 against 129 (escarpment), 1951–2412 against 2235 (plate), 1350–1511 against 2574 (solar
farm, chain) and 1565–1668 (solar farm, 4 4 1).

**The ridge.** Its lee separates. Near the ground (the first layer, at y = 100 m) the reference's reversed flow runs
from x = 1687.5 to 2312.5 m and the parent's from 1725 to 2475 m. The patches end theirs at different places:
1837.5–2312.5 m at D = 900 m (its face at x = 2500 m), and 1837.5–2212.5 m for a patch of D = 1500 m (faces at 100
and 3100 m, ROI 8.553 %). Two runs on the domain of that patch, both with the reference's log-law inflow and an
outflow face and on the reference's cells, separate the cause:

| Run (x from 100 to 3100 m, 25 m cells) | First guess | Reversed flow (m) | ROI |
|---|---|---|---|
| inflow and outflow faces | log-law profile, as the reference | 1687.5–2312.5 | 0.184 % |
| inflow and outflow faces | the parent's field (`init.first_guess = nest`) | 1837.5–2187.5 | 8.559 % |

Both reach the steady verdict (2152 and 1180 iterations). The first matches the reference to the effect of a 100 m
shorter fetch. The second matches the patches. So on these cells the steady problem over the ridge has more than one
converged solution, and the parent's field leads to another one than the profile does; the nest faces are not
needed for the difference. Deeper convergence (seven orders, 5152 iterations for the reference) and an outflow face
in place of the nest face at x = 3100 m left the patch's numbers unchanged (8.553 and 8.550 %). With the ground blend of the hybrid closure off
the parent gives 8.44 % and the patches of D = 900 and 1500 m 5.20 and 5.27 %. At D ≤ 900 m the face at the
downwind end lies in or next to the parent's reversed flow, so the faces play a part there as well.

**Refinement ratio (`bc.nest_ratio`, beyond the plan).** The cube from a 20 m parent, to 5 m cells in one step
(ratio 4) or in two (ratio 2 to 10 m, then ratio 2 to 5 m). D takes two half-widths, `DX,DY`, because the cube's
centre (130, 120) m does not lie on a 20 m face in y. ROI velocity RMS against the 5 m reference, four boxes on two
ranks:

| Run | Cells | Iterations | ROI | Largest on the face cells |
|---|---|---|---|---|
| 20 m parent | 24 x 12 x 24 | -- | 42.03 % | 31.1 % |
| 5 m patch, ratio 4, D = 70,60 m | 28 x 24 x 24 | 106 | 3.969 % | 32.4 % |
| 10 m patch, ratio 2, D = 90,100 m | 18 x 20 x 24 | 153 | 19.20 % | 14.6 % |
| 5 m patch, ratio 2 from the 10 m patch, D = 70,60 m | 28 x 24 x 24 | 108 | **2.598 %** | 15.5 % |

The ratio-4 patch and the chain differ by 1.399 %. The ratio-4 patch on 1, 2 and 4 ranks gives the same field
(largest difference 0 of U_ref). Commands:
```bash
python3 nest_shapes.py run cube cube_par20 20 -- mpiexec -n 2 <build>/gust inputs_building_domain "mesh.max_grid_size=12 6"
python3 nest_shapes.py nest cube cube_r4 cube_par20 70,60 --nest-ratio 4 -- mpiexec -n 2 <build>/gust inputs_building_nest "mesh.max_grid_size=14 12"
python3 nest_shapes.py nest cube cube_c10 cube_par20 90,100 -- mpiexec -n 2 <build>/gust inputs_building_nest "mesh.max_grid_size=9 10"
python3 nest_shapes.py nest cube cube_chain cube_c10 70,60 -- mpiexec -n 2 <build>/gust inputs_building_nest "mesh.max_grid_size=14 12"
python3 nest_shapes.py compare cube cube_r4 cube_ref cube_par20
```

The solar farm's 0.25 m patches straight from the 1 m parent (4 4 1) against the chain's (2 2 1 twice), velocity RMS
of their difference, % of U_ref:

| D (m) | ROI | Lee | Faces | Largest in the patch |
|---|---|---|---|---|
| 32 | 0.287 | 0.756 | 0.736 | 1.874 |
| 40 | 0.117 | 0.256 | 0.498 | 0.979 |
| 48 | 0.052 | 0.103 | 0.316 | 0.657 |

The two shrink together as D grows. Against the reference: the chain's ROI RMS is 0.092, 0.036 and 0.011 %, the 4 4 1 patches' 0.196, 0.082 and 0.042 %; the chain is closer at every D, as for the cube above. The chain costs one more run (the middle
patch, 224 256 cells, 2914 iterations); the 4 4 1 patches took 10 to 16 % more iterations than the chain's.

**Pressure tolerance of the solar farm.** The solar runs set `simplec.pressure_reltol` = 1e-6 (default 1e-3). At
1e-3 (the same runs with that key left out of `SHAPES["solar"]`) the k residual rose by a factor of 3 or more between
two printed iterations (every 10th) 19 times in the parent after iteration 500, the last at 2190, and 10 times in the
middle patch, the last at 2860, while u, w and r went on falling. At 1e-6 there were no such jumps in any solar run.
In the parent at iteration 1200, 90 % of the k residual's sum of squares was in the first layer behind the array
(x = 110 to 160 m). The ground's k row there sets k to its wall value k_w(u*) with the weight u*/d0, d0 = 0.01 m the
first cell's centre and u* about 0.35 m/s, so about 35 1/s. k follows k_w each iteration, and the row's residual shows
the change of the first cell's velocity from one iteration to the next, which the pressure solve's error at 1e-3 made
non-monotone. At iteration 2400 of the parent the fields of 1e-3 and 1e-6 differ by at most 3.5e-8 of U_ref in
velocity and 3.5e-9 of U_ref^2 in k, and the iterations to the steady verdict hardly change (4668 and 4671 for the
parent, 2927 and 2914 for the middle patch, 1500 and 1435 for the chain's D = 40 m patch). An iteration costs more at
1e-6: on the shared laptop the 0.25 m patches took 3 to 4 times as long per iteration, the middle patch 1.3 times and
the parent the same (not a timing study). The first cell's centre, at 0.01 m, lies within the roughness sublayer of
z0 = 0.01 m; the wall law stays defined (it uses ln((d + z0)/z0)), and the case measures the nesting error, not
accuracy.

## CTest
`rt_nest_shapes_cube` (smoke) runs the parent, the reference and the D = 70 m patch of the cube on two ranks, four
boxes each, and passes when the patch is within 2 % of U_ref of the reference around the cube (RMS; measured
1.394 %) and at most a quarter of the parent's difference (18.431 %); fixtures `rt_nest_shapes_clean`,
`rt_nest_shapes_cube_parent`, `_cube_reference` and `_cube_d70`. The ratio runs above are smoke tests too:
`rt_nest_shapes_cube_ratio4` (ROI at most 5 % and 0.15 of the parent's), `rt_nest_shapes_cube_chain2` (at most 3.5 %
and 0.1 of the parent's), `rt_nest_shapes_cube_ratio4_vs_chain` (at most 2 %) and `rt_nest_shapes_cube_ratio4_ranks1`
and `_ranks4` (the ratio-4 patch on 1 and 4 ranks within 1e-8 of U_ref of its 2-rank run on every patch cell);
fixtures `_cube_parent20`, `_cube_r4`, `_cube_r4_np1`, `_cube_r4_np4`, `_cube_c10` and `_cube_chain`.
`rt_nest_shapes_study` (nightly; fixture `rt_nest_shapes_study_clean`) runs the study for the hill, the escarpment, the cube and the plate on four ranks, with
`--ranks hill,cube`; the ridge is left out because it misses T1 and T2. The solar farm is not in it either: its
runs take hours on a laptop (the 0.25 m reference has 1.1 M cells), so it is run by hand with
`study --shapes solar --ranks solar`. The buildings in a patch are also tested by
`rt_nest_faces_building_inside` and `abort_nest_faces_buildings` (`Tests/Canonical_RANS/Nest_Faces`) and by the
unit tests `NestBuildings.*`.
