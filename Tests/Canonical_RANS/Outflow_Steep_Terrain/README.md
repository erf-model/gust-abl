# Outflow_Steep_Terrain

## Purpose
Steep terrain at an outflow face. A 2-D domain (2400 m x 32 m, y periodic) ends on the steepest part of a rising
cos^2 ridge (the ridge of Ishihara and Zhou, 40 m high with a 100 m half-width, crest 50 m beyond the outflow face;
max slope 32 deg at the face). Up to v1.0.0, with the face cutting the slope, the steady iteration alternated
between two states in the last column. With the terrain levelled beyond x = 2380 m (8.2 m high, the slope 26 deg
just inside the domain), so that only the face has no normal slope, it converged. The ghost cells outside the outflow
face now carry the transport fluxes of the cells inside, mirrored through the face (`bc.outflow_transport = mirror`,
the default), and both terrains converge. Found on Big Southern Butte, where the eastern outflow face crosses a ridge and
a valley (`Tests/Paper/JWEIA2026`). Beyond the plan.

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
| **Limit only (this deck), `bc.outflow_transport = mirror` (default)** | **steady, 148 iterations** | **steady, 162 iterations** |
| Limit only, `bc.outflow_transport = from_velocity` (up to v1.0.0) | two alternating states, r = 0.176 (0.178 before the layout-free norms of #60) | steady, 161 iterations |
| Without the sensor (`turb.realizability = 0`), `from_velocity` | two alternating states, r = 0.126 | steady, 745 iterations |
| Limit and near-ground blending (`turb.nut_relax = 0.5`, `convection.hybrid_sigma_max = 0.1`, `convection.hybrid_height = 100`), `from_velocity` | two alternating states, r = 0.177 | steady, 162 iterations |

The last two rows were measured with v1.0.0, before #60.

**Why `from_velocity` stalls.** The ghost node columns outside the face take the face node's height, so the ghost cell is
flat across the face, and with w copied its contravariant velocity Omega = w - z_x u - z_y v is w itself (+0.89 m/s
next to the ground) while the edge cell's is -0.16 m/s. The control volume of the face-point u reaches half a cell into the ghost
cell and sees a vertical flux the domain does not have, which the implicit diagonal does not hold; the face-point
rows overshoot and the state alternates between two values (u at the first cell 1.149 and 1.271 m/s, period 2). The
residual norms are equal on both phases, so the log looks frozen. From the alternating state, `steady.cfl_max = 1`
converges in 37 iterations: a steady solution exists. `bc.outflow_transport = mirror` makes each ghost cell's
transport fluxes the mirror image of the cell inside it (Omega and the tangential flux mirrored, the normal flux
continued linearly through the face), so the face-point control volume has the edge cell's net flux; the velocity
ghosts stay copies (`docs/boundary_conditions.rst`).

On Big Southern Butte without levelled edges (`A_caponly_bsb`), `from_velocity` was not steady after 4000
iterations (r = 3.5e-3, the residual in the last column at the eastern face near the ground); with `mirror` it is
steady after 793 iterations, as with levelled edges, and no residual is left at the outflow columns. Levelling the terrain towards the boundaries (Sorensen et al. 2012) is
no longer needed for convergence.

## CTest
- `outflow_steep_terrain_prepare`: writes the terrain files (fixture).
- `rt_outflow_steep_terrain_flat`: terrain levelled at the outflow face, steady within 400 iterations (162).
- `rt_outflow_steep_terrain_cut`: the cut slope, steady within 400 iterations (148). Fails with the v1.0.0
  condition (its mutant: `bc.outflow_transport = from_velocity`, the next test).
- `rt_outflow_steep_terrain_from_velocity`: the cut slope with `bc.outflow_transport = from_velocity` is not steady after 400 iterations and
  its relative residual is still above 0.1 (the old behaviour, kept reachable).

All label `regression`, 2 ranks, about 3 min together.
