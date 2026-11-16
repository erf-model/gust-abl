# Rampart_2D

## Purpose
How abrupt ground slows the steady iteration, on a 2-D analytic cliff. Beyond the plan. The mesh report prints the
largest step of the ground between two neighbouring columns in units of their thinner first layer
("max ground step / first layer thickness") and warns above 64, the largest value of the runs of this case that
reached the steady verdict within a few thousand iterations. The case shows what that number does and does not predict.

A cliff h(x) = H/2 (1 + tanh((x - x_c) / W)), H = 900 m, stands in a 16 km x 400 m domain
(80 x 2 x 60 cells of 200 m, y periodic, 60 layers from `mesh.dz0` by a fixed ratio to a top at 10 km, BTF). The wind
climbs the cliff: a log law over the sea (u* = 0.3031 m/s, z0 = 0.0002 m) on xlo, outflow on the plateau,
a no-slip lid moving at the log-law speed of the top, Rayleigh damping of w above 8 km, neutral. The recommended
closure setting (`turb.realizability = 0.309`, `turb.nut_relax = 0.5`, the near-ground blend) and a solver
setting for steep terrain (line LU-SGS, two sweeps, CFL cap 16, `simplec.w_vnn_cap = 16`).

## Files
| File | Role |
|---|---|
| `make_terrain.py` | Writes `rampart_steep.xyz` (W = 170 m: 476 m between two nodes, 69 deg) and `rampart_gentle.xyz` (W = 750 m: 119 m, 31 deg), and prints the layer ratios for `mesh.dz0` = 2 and 8 m. |
| `inputs_rampart` | The deck. `terrain.file`, `mesh.dz0`, `mesh.stretch_ratio` and `steady.max_iterations` on the command line. |
| `check_rampart.py` | The CTest: three runs, the closed-form ground step and aspect ratio, the warning, the iteration counts. |

## Data
None: the terrain is analytic.

## How to run
```bash
python3 Tests/Canonical_RANS/Rampart_2D/make_terrain.py out
cd out
mpirun -np 2 ../build-rel/gust ../Tests/Canonical_RANS/Rampart_2D/inputs_rampart terrain.file=rampart_steep.xyz \
    mesh.dz0=2 mesh.stretch_ratio=1.111114 steady.max_iterations=20000
```
The other rows of the table: `rampart_gentle.xyz`, and `mesh.dz0=8 mesh.stretch_ratio=1.079941`.

## Expected Results
Measured on 2 ranks; iterations to the steady verdict (3 orders below the reference and 5e-5).

| Cliff | First layer | Max slope | Max aspect ratio | Ground step (first layers) | Steady after |
|---|---|---|---|---|---|
| gentle (W = 750 m) | 8 m | 30.8 deg | 27.5 | 15.7 | 268 |
| gentle | 2 m | 30.8 deg | 109.9 | 62.9 | 260 |
| steep (W = 170 m) | 8 m | 67.2 deg | 27.5 | 63.9 | 2570 |
| steep | 2 m | 67.2 deg | 109.9 | 255.5 (warning) | 15383 |

- The aspect ratio alone does not set the cost: on the gentle cliff four times thinner layers change nothing (268
  and 260); on the steep cliff they cost six times the iterations (2570 and 15383).
- The ground step alone does not either: at the same value (63 and 64) the steep cliff takes ten times the
  iterations of the gentle one. The number marks meshes outside the tested range; it is not a predictor of the
  iteration count.
- Every run converges, the residuals falling steadily. The case reproduces the slow convergence of abrupt ground;
  it has no lee side and no repeating cycle.
- The report's slope is that of the mesh (the steepest pair of nodes), below the 69 deg of the analytic cliff.

**The ground step and the answer.** The ground step is the ratio of Mahrer's (1984) consistency condition for
terrain-following coordinates, which asks that the ground change by no more than a layer between neighbouring
columns (a ratio of 1 at most). The gentle cliff with 8 m layers on 200, 100, 50 and 25 m cells (80, 160, 320 and
640 x 2 x 60 on 4 ranks with `mesh.max_grid_size = 80 2`; ground step 15.7, 7.9, 3.9 and 2.0; steady after 231, 252,
304 and 933 iterations, where the table's 268 is of boxes of 40 on 2 ranks) tests what a
ratio above 1 costs in this neutral case. At 24 points (six stations from the foot to the plateau, 4, 20, 80 and
300 m above the ground) the largest difference of u to the 25 m mesh is 0.14 m/s at 200 m, 0.03 m/s at 100 m and
0.007 m/s at 50 m (speeds of 4 to 14 m/s), and it falls by a factor of 4 to 5 per halving. A ground step of 16 is
then a discretisation error of about 1 % here, not a failure. The sweep is neutral and has no lee side; a
stratified run, where the condition comes from the pressure-gradient error, was not made. To repeat it:
`make_terrain.py out 25` and a copy of the deck with `mesh.n_cell = N 2 60` and `mesh.max_grid_size = 80 2`.

The closed form of the ground step: the largest of |h(x + dx) - h(x)| / (dz0 (1 - h / z_top)) over the node pairs,
with h the higher of the two (the BTF first layer at a node is dz0 (1 - h / z_top)).

## CTest
- `rt_rampart_ground_step` (label `regression`, 2 ranks, about 30 s): gentle with 2 m layers and steep with 8 m to
  the verdict, steep with 2 m for one iteration. The ground step and the aspect ratio within 0.1 % of the closed
  forms; the warning on the steep cliff with 2 m layers only; gentle steady within 400 iterations, steep with 8 m
  within 4000 and more than four times the gentle count.
