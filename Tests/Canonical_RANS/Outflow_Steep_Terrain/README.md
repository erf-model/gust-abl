# Outflow_Steep_Terrain

## Purpose
Steep terrain at an outflow face. A 2-D domain (2400 m x 32 m, y periodic) ends on the steepest part of a rising
cos^2 ridge (the ridge of Ishihara and Zhou, 40 m high with a 100 m half-width, crest 50 m beyond the outflow face;
max slope 32 deg). With the face cutting the slope, the steady iteration stops at a fixed state; with no slope normal
to the face over the last 100 m it converges. Found on Big Southern Butte, where the eastern outflow face crosses a
ridge and a valley (`Tests/Paper/JWEIA2026`). Beyond the plan.

## Files
| File | Role |
|---|---|
| `make_terrain.py` | Writes `outflow_slope.xyz` (the cut slope) and `outflow_slope_flat.xyz` (no normal slope beyond x = 2300 m). |
| `inputs_outflow_slope` | Log-law inflow (u* 0.21 m/s, z0 0.01 m), outflow, slip top, the grid of the Ishihara-Zhou case (4 m cells, 0.5 m first layer stretched by 1.081), the realizability limit without the near-ground blending. `terrain.file` on the command line. |

## How to run
```bash
python3 Tests/Canonical_RANS/Outflow_Steep_Terrain/make_terrain.py out
mpirun -np 2 build-rel/gust Tests/Canonical_RANS/Outflow_Steep_Terrain/inputs_outflow_slope terrain.file=out/outflow_slope.xyz
mpirun -np 2 build-rel/gust Tests/Canonical_RANS/Outflow_Steep_Terrain/inputs_outflow_slope terrain.file=out/outflow_slope_flat.xyz
```

## Expected Results
- No normal slope at the outflow: steady after 206 iterations.
- Cut slope: the relative residual stays at 0.178 from iteration 250 on (fixed state, not an oscillation).
- With the tunnel inflow table of the Ishihara-Zhou case instead of the log law (not in the repository), the same
  result: without the sensor r = 0.126, with the limit r = 0.180, with the limit and the near-ground blending
  r = 0.174 (all fixed), and 214 iterations with no normal slope. The blending therefore does not remove this stall.
- On Big Southern Butte the stall sits in the last two columns at the eastern outflow face; giving the terrain no
  normal slope over the last 680 m makes the limit alone converge with the same station errors as the run with the
  blending.

The likely cause is the outflow condition, which copies values along terrain-following grid lines; on a sloping face
that is not a zero gradient normal to the boundary. Not fixed here. Common practice, which this case follows, is to
level the terrain towards the domain boundaries (Sorensen et al. 2012).

## CTest
- `outflow_steep_terrain_prepare`: writes the terrain files (fixture).
- `rt_outflow_steep_terrain_flat`: no normal slope at the outflow, steady within 400 iterations.
- `rt_outflow_steep_terrain_cut`: the cut slope is not steady after 400 iterations (the known limitation; when the
  outflow condition is fixed, this test fails and should be changed). The cut-slope terrain given to the first test
  makes it fail, which is its mutant.

All label `regression`, 2 ranks, about 3 min together.
