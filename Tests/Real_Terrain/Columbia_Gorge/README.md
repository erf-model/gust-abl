# Columbia_Gorge

## Purpose
Two-way nesting (PR-70, `nest.feedback = defect`) on steep real terrain, beyond the plan. The domain is the Columbia
River Gorge around the WFIP2 site near Wasco, Oregon (Shaw et al. 2019): the river at about 50 m above sea level, the
Columbia Hills to the north rising about 850 m above it, and the plateau to the south at about 300 m. The nest covers
the river bend and the south face of the hills. Its north face cuts across the ridge and its east face crosses the
river, so the faces stand where the flow changes most. The case tests the nesting under these conditions: it is not a
comparison with the WFIP2 observations, and none are used.

Three runs, one after another, all with the recommended lee setting of `tools/make_gust_case.py`. Its
`turb.nut_relax = 0.5` does not change the steady state, but it decides whether the outer run reaches it: without
`turb.nut_relax_end` the outer domain stalls at a relative residual of 3.9e-3 in a two-iteration cycle on the south
face of the Columbia Hills (`docs/turbulence_model.rst`, Ending the relaxation of nu_t). `prepare` therefore adds
`turb.nut_relax_end = 0.01` to the decks, which ends the relaxation below a relative residual of 0.01. The key is not
part of the recommended setting: on the Ishihara and Zhou rough ridge it stops a run that converges without it.
1. the reference: the whole domain on 100 m cells;
2. one way: the outer domain on 200 m cells, then the nest on 100 m cells (`nest.feedback = none`);
3. joint: the same, then the joint run (`nest.feedback = defect`).

With `--relax F ...`, `run` then makes one joint run per factor with feedback by relaxation (PR-71,
`nest.feedback = relax`, `nest.feedback_factor = F`; `docs/boundary_conditions.rst`, Two-way nests by relaxation),
and `compare` adds them to its tables.

`compare` gives the RMS difference of each run from the reference: over the nest, over the outer cells outside the
nest's footprint, over the footprint, and downstream (east) of the nest. It then counts the cells in the lowest
200 m whose horizontal speed differs from the reference's by more than 0.25, 0.5 and 1 m/s
(`tools/speed_difference.py`), and with `--plot map.png` draws each column's largest difference on the terrain. The
nest is compared at its own cell heights: a joint nest's levels follow the parent's and differ from the
reference's.

## Files
| File | Role |
|---|---|
| `case.yaml` | The case for `tools/make_gust_case.py`: 19.2 km x 16 km centred on 45.70 N, 120.80 W, top 3000 m, 200 m cells, first layer 10 m, 37 layers; wind 8 m/s at 10 m from the west (270 degrees), z0 = 0.05 m, neutral; a 1000 m flat border and a 2000 m slope to the terrain on every side. |
| `columbia_gorge.py` | `prepare` (fetch the terrain, write the decks), `run` (the three runs; with `--relax F ...` also the relax runs, with `--relax-only` only those), `compare` (the tables below and, with `--plot`, the map). |

The deck written by `prepare`: 96 x 80 x 37 cells; inflow on xlo, outflow on xhi, slip walls on y; a Rayleigh layer
on w above 2100 m (as every nest has, so the outer domain and the nest are damped alike); three hierarchy levels with
grid sequencing. The terrain spans 0 to 906 m after the lowest point (49 m) is subtracted. The reference deck
`inputs_ref` is the same with 192 x 160 x 37 cells. The nest `rg` runs from (-4000, -2400) to (5600, 4800) m:
96 x 72 x 37 cells of 100 m.

## Data and acknowledgements
The terrain is the SRTM 1 arc-second set from the Terrain Tiles on the AWS Open Data Registry
(<https://registry.opendata.aws/terrain-tiles/>, skadi format, made by Mapzen / Tilezen joerd). Over the United
States the tiles combine 3DEP (formerly NED), SRTM and GMTED2010 terrain data courtesy of the U.S. Geological Survey;
the joerd attribution list (<https://github.com/tilezen/joerd/blob/master/docs/attribution.md>) gives every source.
The U.S. sources (3DEP, SRTM, GMTED2010) are U.S. Government works in the public domain; the tiles may be used
with the attribution in the joerd list. We thank the U.S. Geological Survey, NASA and Mapzen / Tilezen for making these data freely available.

`prepare` downloads tile N45W121 and N45W122 (no key or account is needed) with `tools/fetch_terrain.py`,
which writes the attribution beside `dem.asc` (`dem.sources.txt`). No terrain is stored in the repository.

The site is that of the Second Wind Forecast Improvement Project (WFIP2): Shaw, W.J., et al. (2019), The Second Wind
Forecast Improvement Project (WFIP2): general overview, *Bull. Amer. Meteor. Soc.* 100, 1687-1699,
doi:10.1175/BAMS-D-18-0036.1. No WFIP2 data are used.

## How to run
```bash
python3 columbia_gorge.py prepare
python3 columbia_gorge.py run --gust <build>/gust --mpi "mpiexec -n 4"
python3 columbia_gorge.py compare --plot gorge_speed_map.png
```
The relax runs at factors 0.5, 1 and 4, after the three runs above (`--relax-only` skips those):
```bash
python3 columbia_gorge.py run --gust <build>/gust --mpi "mpiexec -n 4" --relax 0.5 1 4 --relax-only
```
Everything is written to `columbia_gorge_run/` in the current directory (`--work` changes it).

## Expected Results
Measured on 4 ranks with the decks `prepare` writes (terrain 0 to 906 m, `turb.nut_relax_end = 0.01`). Every run
reaches the steady verdict: the reference after 197 iterations, the outer domain after 164 and the one-way nest after
207; the joint run after 279 exchanges. RMS against the 100 m reference (m/s, k in m^2/s^2):

| Field | nest, one way | nest, joint | outer, one way | outer, joint | footprint, one way | footprint, joint | downstream, one way | downstream, joint |
|---|---|---|---|---|---|---|---|---|
| u | 0.122 | **0.120** | 0.102 | **0.088** | 0.202 | **0.159** | 0.148 | **0.110** |
| v | 0.033 | **0.030** | 0.034 | **0.031** | 0.055 | **0.039** | 0.046 | **0.039** |
| w | 0.021 | **0.018** | 0.026 | **0.025** | **0.043** | 0.067 | 0.030 | **0.027** |
| k | **0.028** | 0.031 | 0.036 | **0.030** | 0.056 | **0.035** | 0.054 | **0.036** |

Cells in the lowest 200 m whose horizontal speed differs from the reference's by more than 0.25, 0.5 and 1 m/s
(reference mean 9.6 m/s):

| Run | Region | Cells | Largest (m/s) | > 0.25 m/s | > 0.5 m/s | > 1 m/s |
|---|---|---|---|---|---|---|
| one way | nest | 86999 | 1.37 | 11005 (12.7 %) | 3269 (3.8 %) | 3 (0.00 %) |
| one way | footprint | 21752 | 2.87 | 7337 (33.7 %) | 2726 (12.5 %) | 134 (0.6 %) |
| one way | outer | 76387 | 2.32 | 6204 (8.1 %) | 1526 (2.0 %) | 87 (0.11 %) |
| joint | nest | 87009 | 1.44 | 10704 (12.3 %) | 3128 (3.6 %) | 6 (0.01 %) |
| joint | footprint | 21752 | 2.47 | 4813 (22.1 %) | 1296 (6.0 %) | 111 (0.5 %) |
| joint | outer | 76387 | 2.32 | 5105 (6.7 %) | 1050 (1.4 %) | 79 (0.10 %) |

The joint run lowers every RMS difference except w over the footprint (0.043 to 0.067) and k in the nest (0.028 to
0.031). In the nest it lowers the counts above 0.25 and 0.5 m/s but has a larger largest difference (1.44 against
1.37 m/s, 6 cells above 1 m/s against 3). The footprint w is the cost of the nest's levels: they follow the parent's surfaces aloft (decay height
200 m here), so the outer cells under the nest see the nest's flow on a mesh that bends faster than their own. A
slower fade (a decay height of 1 km, or a linear fade as ERF's `amr_terrain_refinement = transform`) lowers the
nest's and downstream speed differences a little further but raises the footprint w RMS to 0.12-0.15; keeping the
parent's interpolated ground in the nest (as ERF's `interpolate`) did not converge. These runs were made while
choosing the rule (PR-70, on a 903 m terrain mosaic) and are not options of the code. The largest outer differences
(2.3 m/s) sit where 200 m cells cannot follow the terrain; nesting does not change them.

### Feedback by relaxation (`--relax 0.5 1 4`)
Measured on 4 ranks with the same decks (RMS of u against the reference, m/s; exchanges of the joint run):

| Run | Exchanges | nest | outer | footprint | downstream |
|---|---|---|---|---|---|
| one way | | 0.125 | 0.105 | 0.204 | 0.152 |
| joint, `defect` | 276 | **0.122** | 0.093 | **0.164** | 0.117 |
| relax, f = 0.5 | 174 | 0.126 | **0.092** | 0.166 | 0.116 |
| relax, f = 1 | 180 | 0.128 | **0.092** | 0.167 | **0.115** |
| relax, f = 4 | 439 | 0.147 | 0.099 | 0.180 | 0.133 |

Cells in the lowest 200 m whose horizontal speed differs from the reference's by more than 1 m/s: over the footprint
118 (f = 0.5), 114 (f = 1) and 104 (f = 4), against 169 with `defect` and 197 one way; over the nest 60, 62 and 74,
against 59 and 71. At f = 0.5 and 1 the outer regions match the defect run in fewer exchanges, and the nest's u is
slightly further from the reference than one way. The footprint w is 0.045 at f = 0.5 and 0.050 at f = 1 (`defect`
0.068). At f = 4 every region moves away from the reference (`docs/boundary_conditions.rst`, Two-way nests by
relaxation).

These numbers use the terrain of the PR-70 runs. On the terrain the current `prepare` writes (relief 906 m rather
than 903 m), the one-way outer run does not reach steady with `turb.nut_relax = 0.5`; this is being followed up
separately.

## CTest
None: the three runs take about 15 minutes on 4 ranks of a laptop. The CTest of two-way nesting is `rt_nests_two_way`
(`Tests/Canonical_RANS/Nests_Two_Way`).
