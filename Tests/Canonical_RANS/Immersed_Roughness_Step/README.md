# Immersed_Roughness_Step

## Purpose
A roughness map (`bc.surface_z0_file`, [docs](../../../docs/immersed_terrain.rst)) under an immersed ground
(`terrain.representation = immersed_forcing`) and under the ground of a hybrid mesh
(`terrain.representation = hybrid`), against the same flow on a body-fitted mesh with the same map.

The channel is that of `../Immersed_Open_Faces`: 2560 m long and 1000 m high above a flat ground, neutral
(uniform θ, no gravity), `bc.xlo = inflow` (log law, u* = 0.35 m/s, z0 = 0.03 m), `bc.xhi = outflow`, periodic in
y, a slip top, 64 x 4 columns and layers of 15.625 m. The roughness length is 0.03 m up to x = 1120 m and 0.3 m
beyond: the flow slows down over the rougher ground and a new internal layer grows from the step. The step lies
inside a box (the boxes are 16 columns long and 2 wide), so a box holds both values.

Three grounds are compared with the body-fitted one:

- immersed, three more layers below the ground, the ground at solid fraction F of the third layer;
- hybrid with `terrain.flat_height = 0`: the terrain is flat, so every column's lower face is the ground and the
  surface layer acts on it (the run equals the body-fitted one);
- hybrid with `terrain.flat_height = H dz`: every column's ground cuts the first layer.

In the immersed columns the map enters the wall law of the cut cell, the wall value of k, the cut cell's strain
and, with `turb.shear_stencil = log_law`, the log-law gradient above it. A column's z0 is read by all four; a
run that read a uniform value, or another column's, differs from the body-fitted run downstream of the step.

## Files
| File | Role |
|---|---|
| `z0_step.xyz` | The map: two points, z0 = 0.03 m at x = 640 m and 0.3 m at x = 1600 m. With `terrain.idw_neighbors = 1` each column takes its nearest point's value: the step is at x = 1120 m. |
| `z0_step_other.xyz` | The same with 0.2 m downstream: a changed map for the restart test. |
| `inputs_body_fitted_step` | The body-fitted channel. |
| `inputs_immersed_step` | The immersed channel. `terrain.flat_height` is set by the script. |
| `inputs_hybrid_step` | The hybrid channel. `terrain.flat_height` is set by the script. |
| `check_roughness_step.py` | Runs the three decks with the map and the body-fitted deck with a uniform z0 = 0.03 m, and compares every column. With `--limits` it is the test. |
| `inputs_map_checkpoint`, `inputs_map_restart` | Two iterations that write a checkpoint with the map, and a restart from it on another box layout. They leave out the convection blend, which a restart does not support. |

## Data and acknowledgements
None: the case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "z0_map|rt_immersed_roughness_step"
```
By hand:
```bash
python3 check_roughness_step.py <build>/gust <folder> [--limits DU DK DKW] [--mpi "mpiexec -n"] [--fresh] \
        [key=value ...] [--immersed F ...] [--hybrid H ...]
```
`key=value` arguments replace or add deck keys in every run (each run folder holds its deck as `deck`; `key=` with
no value removes a key); `terrain.*` keys go to the immersed and hybrid runs only. The default grounds are
F = 0, 0.3, 0.7 and H = 0, 0.5. With `--mpi` the runs take 2, 4, 1, 2, ... ranks in turn. A run already in the
folder is read, not repeated, unless `--fresh` is given.

## Expected Results
Every run reaches the steady verdict. The columns are compared with the body-fitted run at equal height above the
ground (the body-fitted column interpolated in ln(z + z0), z0 the column's), in percent of the body-fitted value:

- `du`: u in every column, the fluid cells above the wall cell;
- `dk far`: tke the same way, in the columns upstream of the step and more than 400 m downstream of it;
- `dk near`: tke in the first 400 m downstream of the step. There the new internal layer grows from less than one
  layer to a few layers, and tke changes across it within one or two cells: a cut cell's centre and a value between
  two cell centres of the body-fitted run sample it at different places. The number is printed and has no limit;
- `dk wall`: tke in each column's wall cell (the cut cell, or the first cell of a ground on a face) against the
  body-fitted first cell, the columns at the step included.

With `turb.shear_stencil=log_law terrain.small_volfrac=0.1` (the test's setting; the decks hold the realizability
cap 0.309, `turb.nut_relax = 0.5` and the convection blend). The body-fitted run takes 121 iterations with the map
and 119 with the uniform z0; the map changes u of the body-fitted run by up to 24.3 %.

| Ground | iterations | du (%) | dk far (%) | dk near (%) | dk wall (%) |
|---|---|---|---|---|---|
| immersed, F = 0 | 120 | 0.07 | 0.62 | 0.47 | 0.21 |
| immersed, F = 0.3 | 121 | 0.96 | 8.3 | 27.1 | 4.6 |
| immersed, F = 0.7 | 126 | 2.81 | 14.0 | 51.6 | 10.2 |
| hybrid, H = 0 | 121 | 0.00 | 0.00 | 0.00 | 0.00 |
| hybrid, H = 0.5 | 125 | 2.04 | 12.5 | 40.8 | 5.9 |

The differences of the cut grounds are those of the cut cell on a uniform ground (`../Immersed_Open_Faces`: 1.4 %
and 8.7 % at F = 0.5 in the outflow column), larger here where the flow is still adjusting to the new roughness.
A wall law that reads the first column's z0 in every box (a one-line change) gives `dk wall` = 64 % and
`dk far` up to 25 %.

A ground whose cut cell `terrain.small_volfrac` merges into the ground (F = 0.9 with 0.1) has its wall cell one
layer higher; the script's `dk wall` does not follow that and such an F is not in the default list.

## CTest
| Test | Label | What it checks |
|---|---|---|
| `rt_immersed_roughness_step` | regression | The comparison above with `--limits 4 18 13` (du, dk far, dk wall), the runs on 2, 4 and 1 ranks in turn (8 boxes); the map's effect on u must be at least 20 %. Needs numpy. |
| `immersed_z0_map_clean`, `hybrid_z0_map_clean` | regression, smoke | Empty the run folders of the tests below. |
| `immersed_z0_map_run`, `hybrid_z0_map_run` | regression, smoke | Two iterations on two ranks with the map; the log names the map as the wall law's z0; a checkpoint is written. |
| `immersed_z0_map_restart`, `hybrid_z0_map_restart` | regression, smoke | A restart on another box layout and three ranks with the same map is accepted. |
| `abort_immersed_z0_map_restart_other_map`, `..._uniform`, `..._neighbors` and the three `abort_hybrid_z0_map_restart_*` | regression, smoke | A restart with another map file, with a uniform `bc.surface_z0`, or with another `terrain.idw_neighbors` is refused. |
| `abort_immersed_z0_map_missing`, `abort_hybrid_z0_map_missing` | regression, smoke | A map file that cannot be read stops the run naming `bc.surface_z0_file`. |
