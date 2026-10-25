# Outflow_Steep_Terrain

## Purpose
Steep terrain at an outflow face. A 2-D domain (2400 m x 32 m, y periodic) ends on the steepest part of a rising
cos^2 ridge (the ridge of Ishihara and Zhou, 40 m high with a 100 m half-width, crest 50 m beyond the outflow face;
max slope 32 deg at the face). With the face cutting the slope, the steady iteration stops at a fixed state. With the
terrain levelled beyond x = 2380 m (8.2 m high, the slope 26 deg just inside the domain), so that only the face has
no normal slope, it converges. Found on Big Southern Butte, where the eastern outflow face crosses a ridge and a valley
(`Tests/Paper/JWEIA2026`). Beyond the plan.

## Files
| File | Role |
|---|---|
| `make_terrain.py` | Writes `outflow_slope.xyz` (the cut slope) and `outflow_slope_flat.xyz` (levelled beyond x = 2380 m). |
| `inputs_outflow_slope` | Log-law inflow (u* 0.21 m/s, z0 0.01 m), outflow, slip top, the grid of the Ishihara-Zhou case (4 m cells, 0.5 m first layer stretched by 1.081), the realizability limit without the near-ground blending. `terrain.file` on the command line. |

## How to run
```bash
python3 Tests/Canonical_RANS/Outflow_Steep_Terrain/make_terrain.py out
mpirun -np 2 build-rel/gust Tests/Canonical_RANS/Outflow_Steep_Terrain/inputs_outflow_slope \
    terrain.file=out/outflow_slope.xyz steady.max_iterations=400
mpirun -np 2 build-rel/gust Tests/Canonical_RANS/Outflow_Steep_Terrain/inputs_outflow_slope \
    terrain.file=out/outflow_slope_flat.xyz steady.max_iterations=400
```
Other settings need a copy of the deck (a key in the deck cannot be given again on the command line).

## Expected Results
| Setting | Face cuts the slope | Levelled at the face |
|---|---|---|
| Without the sensor (`turb.realizability = 0`) | fixed state, r = 0.126 | steady, 745 iterations |
| Limit only (this deck) | fixed state, r = 0.178 | steady, 161 iterations |
| Limit and near-ground blending (`turb.nut_relax = 0.5`, `convection.hybrid_sigma_max = 0.1`, `convection.hybrid_height = 100`) | fixed state, r = 0.177 | steady, 162 iterations |

The stall therefore depends on the slope at the outflow face, not on the model setting, and the blending does not
remove it. On Big Southern Butte the stall sits in the last two columns at the eastern outflow face; giving the
terrain no normal slope over the last 680 m makes the limit alone converge with the same station errors as the run
with the blending. The likely cause is the outflow condition, which copies values along terrain-following grid
lines; on a sloping face that is not a zero gradient normal to the boundary. Not fixed here. Common practice is to
level the terrain towards the domain boundaries (Sorensen et al. 2012).

## CTest
- `outflow_steep_terrain_prepare`: writes the terrain files (fixture).
- `rt_outflow_steep_terrain_flat`: terrain levelled at the outflow face, steady within 400 iterations (161).
- `rt_outflow_steep_terrain_cut`: the cut slope is not steady after 400 iterations and its relative residual is
  still above 0.1 there (the known limitation; when the
  outflow condition is fixed, this test fails and should be changed). The cut-slope terrain given to the first test
  makes it fail, which is its mutant.

All label `regression`, 2 ranks, about 3 min together.
