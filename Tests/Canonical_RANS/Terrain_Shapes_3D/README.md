# Terrain_Shapes_3D

## Purpose
Which terrain shapes the steady stage settles on, and how soon. Beyond the plan. Fourteen analytic shapes, all
760 m high, stand in the same neutral log-law wind: round buttes of four slopes, a butte steep on one side only,
ridges across, oblique to and along the wind, two buttes with a gap, a butte with a ridge behind it, a narrow
butte, a crater rim, and a ridge that spans the domain. Each run records the wind at 2560 points after every
iteration (`steady.probes_file`), so the residual verdict can be set against what the wind does at the points:
`tools/probe_settle.py` gives the iteration after which every point stays within 0.5 m/s and within 0.25 m/s of
the end state.

The case asks whether a steep body alone makes the steady iteration repeat a cycle at the points. None of the
three-dimensional shapes does.

The domain is 16 km x 12.8 km x 6 km (160 x 128 x 40 cells of 100 m, y periodic, 40 layers from `mesh.dz0` by a
fixed ratio, STF). The wind is along +x: a log law over the sea (u* = 0.3031 m/s, z0 = 0.0002 m) on xlo,
outflow on xhi, a no-slip lid moving at the log-law speed of the top (12.728 m/s), Rayleigh damping of w above
4.5 km, neutral. The recommended closure setting (`turb.realizability = 0.309`, `turb.nut_relax = 0.5`, the
near-ground blend) and a steep-terrain solver setting (line LU-SGS, two sweeps, CFL cap 16, `simplec.w_vnn_cap = 16`).

A body is h = H/2 (1 - tanh(d / W)) with d the distance outside an ellipse; W sets the slope (600 m: 32 deg,
300 m: 51 deg, 150 m: 67 deg, 90 m: 75 deg on 100 m cells). `make_terrain.py` has the shapes.

## Files
| File | Role |
|---|---|
| `make_terrain.py` | Writes `<shape>.xyz` for the fourteen shapes and `probes.txt` (a lattice every 400 m at 80 m and at 20 m above the ground), and prints the layer ratios for `mesh.dz0` = 2 and 8 m. `--dx 200` for the mesh of the CTest. |
| `inputs_terrain_shapes` | The deck. `terrain.file`, `mesh.dz0`, `mesh.stretch_ratio`, `steady.max_iterations` and `steady.cfl_max` on the command line. |
| `check_terrain_shapes.py` | The CTest: the gentle and the steep round butte on 200 m cells. |
| `tools/probe_settle.py` | Reads a `steady.probes_file` output: when the points settle. Not specific to this case. |

## Data
None: the terrain is analytic.

## How to run
```bash
python3 Tests/Canonical_RANS/Terrain_Shapes_3D/make_terrain.py out
mkdir out/round_steep && cd out/round_steep && cp ../probes.txt .
mpirun -np 8 ../../build-rel/gust ../../Tests/Canonical_RANS/Terrain_Shapes_3D/inputs_terrain_shapes \
    terrain.file=../round_steep.xyz mesh.dz0=8 mesh.stretch_ratio=1.119188 steady.max_iterations=1500 \
    steady.cfl_max=16
python3 ../../tools/probe_settle.py .
```
With 2 m first layers: `mesh.dz0=2 mesh.stretch_ratio=1.168454`. `probe_settle.py A:B` measures run A against
the end state of run B.

## Expected Results
Measured on 8 ranks, 8 m first layers unless noted, at most 1500 iterations. "Verdict" is the steady verdict
(every variable 3 orders below the reference and below 5e-5); "residual" is the last column `r` of the iteration
line of a run without it. The ground step is counted in the mesh's own first layer, which STF thins to about 4 m
at the edge of a crest (1 m with 2 m layers): with BTF the steep butte has 33 where the table has 53. The
"within" columns are the iteration after which the velocity at every one of the 2560 points stays within the
tolerance of the mean of the run's last 100 iterations. "Last 200" is the largest distance of a point's velocity
in the last 200 iterations from its last value.

| Shape | Slope | Ground step (first layers) | Verdict | Residual | Within 0.5 m/s | Within 0.25 m/s | Last 200 (m/s) |
|---|---|---|---|---|---|---|---|
| `round_gentle` | 32 deg | 8.2 | steady after 137 | | 21 | 29 | |
| `round_mid` | 51 deg | 18.8 | steady after 1322 | | 126 | 207 | 0.001 |
| `round_steep` | 67 deg | 53.2 | none at 1500 | 1.0e-3 | 147 | 251 | 0.009 |
| `round_steep`, 5000 iterations | 67 deg | 53.2 | none at 5000 | 6.5e-5 | 154 | 280 | 0.000 |
| `round_steep`, 2 m layers | 67 deg | 208.7 (warning) | none at 1500 | 1.7e-2 | 84 | 135 | 0.005 |
| `round_sharp` | 75 deg | 79.5 (warning) | none at 1500 | 3.4e-3 | 85 | 117 | 0.002 |
| `steep_lee` | 32 / 67 deg | 54.3 | none at 1500 | 1.0e-3 | 250 | 651 | 0.033 |
| `steep_windward` | 67 / 32 deg | 54.3 | none at 1500 | 5.0e-4 | 98 | 135 | 0.000 |
| `ridge_cross` | 68 deg | 53.2 | none at 1500 | 7.7e-4 | 140 | 215 | 0.010 |
| `ridge_oblique` | 68 deg | 48.1 | none at 1500 | 4.4e-4 | 171 | 196 | 0.000 |
| `ridge_stream` | 68 deg | 53.2 | none at 1500 | 1.9e-4 | 27 | 31 | 0.000 |
| `two_gap` | 67 deg | 52.9 | none at 1500 | 9.9e-4 | 220 | 385 | 0.014 |
| `butte_canyon` | 68 deg | 54.0 | none at 1500 | 1.1e-3 | 565 | 758 | 0.033 |
| `narrow` | 67 deg | 54.0 | none at 1500 | 8.8e-4 | 38 | 65 | 0.002 |
| `rim` | 67 deg | 53.5 | none at 1500 | 1.0e-3 | 800 | 911 | 0.012 |

- No three-dimensional shape repeats a cycle. On every one of them the points move by at most 0.033 m/s over
  the last 200 iterations.
- The points settle long before the residual verdict. The steep round butte has no verdict at 5000 iterations
  (the residual of w is 7.5e-5 against 5e-5), yet every point is within 0.25 m/s from iteration 280, and its state at 1500 is within
  0.03 m/s of its state at 5000.
- The slope sets the cost of the verdict: 137 iterations at 32 deg, 1322 at 51 deg, more than 5000 at 67 deg.
  The points settle after 29, 207 and 251 iterations: between 51 and 67 deg the verdict moves out of reach and
  the points barely notice.
- What is behind the body matters more than its slope for the points: a second ridge in the wake
  (`butte_canyon`, 758) and a crater (`rim`, 911) take the longest; a steep lee side (651) takes longer than a
  steep windward side (135).
- Thinner layers and a sharper edge (ground step 209 and 80, both above the tested 64) raise the residual at 1500
  iterations and do not delay the points.
- The "within" columns measure a run against its own end state, so a run that still drifts can look settled. The
  "last 200" column and the 5000-iteration row are the guard.

**The ridge that spans the domain (`ridge_2d`) is not in the table.** At 67, 51 and 32 deg (W = 150, 300 and
600 m) its residual rises and the CFL number falls to about 0.02 by iteration 1500; the points still move by 0.1
to 0.55 m/s per 200 iterations. The same happens at 32 deg on a domain two cells wide, so it is not an effect of
the periodic direction. One difference from the three-dimensional shapes is that the flow has no way round a ridge across
the whole domain under a lid at 6 km; the 2-D cliff of `Rampart_2D`, which has no lee side, converges. The cause
was not looked for further.

## CTest
- `rt_terrain_shapes_settle` (label `regression`, 2 ranks, about 1 to 2 min): `round_gentle` and `round_steep`
  on 200 m cells (80 x 64 x 40, 8 m layers). The gentle butte reaches the verdict within 250 iterations (136
  measured) with no ground-step warning. The steep one prints the warning (94 first layers) and has no verdict at
  200 iterations, while every point is within 0.5 m/s of the mean of rows 151 to 200 from iteration 40 at the
  latest (26 measured) and within 0.25 m/s from iteration 100 (59 measured), and no point moves more than 0.2 m/s
  between iterations 100 and 200. The gentle butte's points are within 0.25 m/s of its last row from iteration
  40 at the latest (29 measured). The rows are tied to the run: the first guess is more than 1 m/s from row 200
  somewhere, and the gentle butte's first probe ends within 0.02 m/s of the inflow log law at its cell centre
  (9.466 m/s at 72.65 m). Mutant: with `steady.cfl_max = 2` the steep butte is within 0.5 m/s from
  iteration 48 and the gentle one within 0.25 m/s from 53, and both limits of 40 fail.
