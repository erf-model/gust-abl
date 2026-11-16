# Nest_Shapes

## Purpose
One-way nested patches around five canonical shapes, against a uniform fine run of the whole domain (tracker row
PR-76, sub-phase 7.4, beyond plan rev. 6). Each patch is at refinement 2 2 1 of a coarse parent run, with the same
z-levels, and its `nest` faces lie at a distance D from the shape's centre. The reference is the whole domain on the
patch's cells. Three terrain shapes test patches over steep ground with a separated lee; a cube and a raised plate
test buildings inside a patch (`buildings.type` with `nest` faces, allowed since this PR when every building keeps
clear of the faces, `docs/boundary_conditions.rst`). `Tests/Canonical_RANS/Nest_Validation` is the same study on a
gentler hill; `Tests/Canonical_RANS/Nest_Faces` tests the faces themselves.

## Files
| File | Role |
|---|---|
| `inputs_terrain_domain` | The whole domain for the terrain shapes: 3200 m long, periodic in y, 800 m high, 16 layers from 10 m, log-law inflow along x (u* 0.4 m/s, z0 0.1 m), neutral, SIMPLEC to the steady verdict. |
| `inputs_terrain_nest` | A patch over a terrain shape: `nest` faces, first guess from the parent. |
| `inputs_building_domain` | The whole domain for the cube and the plate: flat, periodic in y, log-law inflow along x, the building with the form 7 wall law, z-line sweeps. |
| `inputs_building_nest` | A patch around a building. |
| `nest_shapes.py` | Writes each shape's files (`files`), runs a whole domain (`run`) or a patch (`nest`), compares a patch with the reference (`compare`), or does the whole study (`study`). Pure Python. |

The script gives every run its shape: the terrain or building, the domain size and cells, the roughness, the closure's
`convection.hybrid_height` (2.5 shape heights) and, for the buildings, `turb.max_geom_lscale` = 40 m (the value of
`Tests/Canonical_RANS/Building_Cube`; the terrain decks use the recommended 100 m).

| Shape | Domain (m) | Shape | Parent / reference cells | D (m) | U_ref (m/s) |
|---|---|---|---|---|---|
| ridge | 3200 x 200 x 800, 16 layers | Witch of Agnesi ridge along y, h = 150 m, half-width a = 150 m (largest slope 33 degrees) | 50 / 25 m | 300, 600, 900 | 7.313 |
| hill | 3200 x 2400 x 800, 16 layers | round Witch of Agnesi hill, h = 150 m, a = 150 m | 50 / 25 m | 300, 600, 900 | 7.313 |
| escarpment | 3200 x 200 x 800, 16 layers | z = (H/2)(1 + tanh((x - 1600 m)/L)), H = 100 m, L = 100 m, from an ERF terrain file (`terrain.type = erf_grid`) | 50 / 25 m | 300, 600, 900 | 6.908 |
| cube | 480 x 240 x 120, 24 layers of 5 m | the block of `Building_Cube`, 30 x 30 x 20 m, centre (130, 120) m, z0 0.1 m | 10 / 5 m | 40, 70, 100 | 5.298 |
| plate | 64 x 4 x 20.6, 93 layers (0.025 m to 1.5 m, then x 1.15) | a plate 4 m long, 0.1 m thick, 1 m above the ground, across the periodic width, z0 0.01 m | 1 / 0.5 m | 4, 8, 12 | 4.700 |

The patch is a square of half-width D around the hill and the cube (four `nest` faces), and spans the periodic width
around the ridge, the escarpment and the plate (two `nest` faces, periodic in y). D is measured from the shape's
centre: the nearest building wall is D - 15 m from a face for the cube and D - 2 m for the plate.

## Data and acknowledgements
None. Every shape is analytic; the escarpment's terrain file, the building maps and the plate's z-levels are written by
`nest_shapes.py files SHAPE`.

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

## Measures
The difference between two runs on the same cells is the length of the velocity-vector difference as a fraction of the
shape's speed scale U_ref = (u*/0.4) ln(h/z0), h the shape's height (the building's top for the cube and the plate),
and the difference of k as a fraction of U_ref^2. Both are given as RMS and largest value over four regions:

| Region | Cells |
|---|---|
| roi | around the shape: within a of the crest (ridge, hill), within 200 m of the escarpment's middle, x in [100, 160] m and y in [90, 150] m (cube), x in [19, 25] m (plate); layers centred below the shape's height (ridge, hill, escarpment), 2 cube heights (cube) or 3 m (plate) over flat ground |
| lee | downwind: x from a to 3 a past the crest (ridge, hill; within a of the axis for the hill), 200 to 600 m past the escarpment's middle, x in (145, 165] m behind the cube within its width, x in (24, 28] m behind the plate; the same layers |
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
| T1 | Ridge, hill, escarpment and cube, at the largest D: ROI velocity RMS <= 0.5 % of U_ref and <= 1/4 of the parent's |
| T2 | Ridge, hill, escarpment and cube: ROI velocity RMS does not grow with D (5 % allowance) |
| T3 | With `--ranks`: the middle D in four boxes on 1, 2 and 4 ranks reaches the steady verdict, and the three ROI RMS values agree to 10 % |
| T4 | Every run reaches the steady verdict |
| — | Plate: no accuracy target (tracker row PR-76). Reported: the four regions and the iterations of each patch and of the reference |

The iterations of each patch against those of the reference are reported (PT-7.4.1), with no target. A target that
is missed is reported as missed, with the number.

## Expected Results
To be filled in from the runs.

## CTest
None runs the study. The buildings in a patch are tested by `rt_nest_faces_building_inside` and
`abort_nest_faces_buildings` (`Tests/Canonical_RANS/Nest_Faces`) and by the unit tests `NestBuildings.*`.
