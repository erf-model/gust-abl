# Immersed_Canopy

## Purpose
A forest canopy (`canopy.type`, [docs](../../../docs/canopy.rst)) over an immersed ground
(`terrain.representation = immersed_forcing`) and on a hybrid mesh (`terrain.representation = hybrid`), against the
same flow on a body-fitted mesh.

The folder holds no deck of its own. It runs the decks of three canopy cases as they are, and again with the
ground offset from the mesh's lower face:

| Case | Deck | Flow |
|---|---|---|
| `homogeneous` | `../Homogeneous_Canopy/inputs_homogeneous_canopy` | A uniform 20 m canopy, LAI 5, flat ground, 0.25 m first layer |
| `edge` | `../Forest_Edge/inputs_forest_edge` | Forest, clearing, forest (18 m, LAI 2), flat ground, 1 m first layer |
| `hill` | `../Forested_Hill/inputs_deep_hill` | A 20 m canopy, LAI 10, on cosine hills 4 m high and 160 m apart, 1 m first layer |

Each case is run on six grounds:

- body-fitted: the deck;
- immersed, F = 0, 0.3, 0.7: the deck's layers with two more layers of the first layer's thickness below them
  (`mesh.zlevels_file`) and `terrain.flat_height = (2 + F) dz0`. The ground cuts the fraction F off the third
  layer. On the hill the mesh is flat and the hill is immersed in it;
- hybrid, H = 0, 0.5: the deck's mesh and `terrain.flat_height = H dz0`. On flat ground with H = 0 every column's
  lower face is the ground and the run equals the body-fitted one. On the hill the mesh follows the terrain's
  lower envelope and the ground lies up to 0.62 m (H = 0) or 1.12 m (H = 0.5) above its lower face.

The canopy's height above the ground is measured from the immersed ground. A cell that the ground cuts holds the
leaf area of its part above the ground; a cell below the ground holds none. A canopy measured from the mesh's
lower face would stand 2 to 2.7 layers too low in the immersed runs and leave a part of its leaf area in solid
cells.

## Files
| File | Role |
|---|---|
| `check_immersed_canopy.py` | Writes each case's tables (canopy height, LAI, terrain), runs the six grounds and compares every column with the body-fitted run. With `--limits` it is the test. |

## Data and acknowledgements
None: the cases are analytic. The geometries are those of the three case folders, whose READMEs name the papers
they follow (Finnigan 2007; Dupont and Brunet 2008; Finnigan and Belcher 2004).

## How to run
```bash
ctest --test-dir <build> -R rt_immersed_canopy
```
By hand:
```bash
python3 check_immersed_canopy.py <build>/gust <folder> [--cases NAME ...] [--limits DU DK DKTOP] \
        [--mpi "mpiexec -n"] [--fresh] [key=value ...] [--immersed F ...] [--hybrid H ...]
```
Every run takes `init.first_guess = profile` (the column first guess is not available over an immersed ground),
`steady.sweep = zjacobi`, `steady.zline_sweeps = 16`, `steady.cfl_max = 4` and
`solver.poisson_backend = gmres_fft`. `key=value` arguments replace or add deck keys in every run (each run folder
holds its deck as `deck`; `key=` with no value removes a key); `terrain.*` keys go to the immersed and hybrid runs
only. With `--mpi` the runs take 2, 4, 1, 2, ... ranks in turn. A run already in the folder is read, not repeated,
unless `--fresh` is given.

## Expected Results
Every run reaches the steady verdict. The columns are compared with the body-fitted run at equal height above the
ground (the body-fitted column interpolated linearly in height), over the fluid cells above the wall cell, in
percent of the body-fitted run's largest value of the field:

- `du`: u at every height;
- `dk in`, `dk above`: tke below and above the canopy top, away from it;
- `dk top`: tke within 1.5 layers of the canopy top. There tke rises across one layer (on the hill from 0.20 to
  0.64 m²/s² between the centres at 18.6 m and 20.9 m), and a cell centre between two body-fitted centres samples
  the rise at another place;
- `drag`: the integral of the canopy's drag field over a column's cells with a part above the ground, against
  the body-fitted column's, in percent of C_d LAI. Limit 1e-6 %. (Cells that `terrain.small_volfrac` turns solid
  count here; their share, at most that fraction of one layer, acts on no flow.)

With `turb.realizability=0.309 turb.nut_relax=0.5 turb.shear_stencil=log_law terrain.small_volfrac=0.1
steady.max_iterations=30000` (the tests' setting), measured on macOS arm64 (AppleClang, double precision):

| Case | Ground | iterations | du (%) | dk in (%) | dk above (%) | dk top (%) | drag (%) |
|---|---|---|---|---|---|---|---|
| homogeneous | body-fitted | 293 | | | | | |
| | immersed, F = 0 | 293 | 0.00 | 0.00 | 0.00 | 0.00 | 8e-11 |
| | immersed, F = 0.3 | 292 | 0.14 | 0.31 | 0.10 | 0.28 | 8e-11 |
| | immersed, F = 0.7 | 292 | 0.23 | 0.50 | 0.15 | 0.46 | 8e-11 |
| | hybrid, H = 0 | 293 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| | hybrid, H = 0.5 | 292 | 0.20 | 0.43 | 0.14 | 0.39 | 0 |
| edge | body-fitted | 7580 | | | | | |
| | immersed, F = 0 | 7577 | 0.02 | 0.05 | 0.00 | 0.01 | 1e-10 |
| | immersed, F = 0.3 | 7498 | 0.74 | 3.8 | 0.84 | 2.5 | 1e-10 |
| | immersed, F = 0.7 | 7480 | 1.83 | 8.6 | 1.98 | 2.7 | 1e-10 |
| | hybrid, H = 0 | 7580 | 0.00 | 0.00 | 0.00 | 0.00 | 0 |
| | hybrid, H = 0.5 | 7485 | 1.15 | 6.3 | 1.46 | 2.8 | 1e-14 |
| hill | body-fitted | 10332 | | | | | |
| | immersed, F = 0 | 10373 | 3.1 | 3.2 | 2.9 | 18.0 | 1e-10 |
| | immersed, F = 0.3 | 10390 | 2.9 | 3.3 | 3.5 | 17.8 | 1e-10 |
| | immersed, F = 0.7 | 10421 | 3.2 | 3.2 | 4.6 | 16.0 | 1e-10 |
| | hybrid, H = 0 | 10385 | 2.8 | 1.4 | 4.3 | 17.1 | 6e-14 |
| | hybrid, H = 0.5 | 10538 | 3.0 | 1.9 | 2.5 | 18.6 | 7e-14 |

The drag integral of the immersed runs differs by the twelve digits of the level file.

**Flat ground.** The differences grow with the cut: they are those of the cut cell
(`../Immersed_Roughness_Step`, `../Immersed_Open_Faces`). The largest tke difference of the edge case is in the
first cell above the cut cell, 20 to 40 m behind the forest's leading edge (0.8 m above the ground at F = 0.7),
where tke falls from the canopy flow to the wall value within one layer.

**Hill.** The differences do not depend on the cut: the three meshes place their layers differently against the
canopy top, which follows the ground. On the body-fitted mesh the top lies at the same place in a layer in every
column; on the flat mesh and on the hybrid mesh it moves through the 2.3 m layers along the hill. Away from the
top the fields agree within 5 %.

A canopy measured from the mesh's lower face (a one-line change) gives, for the homogeneous case, a drag integral
2.5 % low over the fluid cells, 5.0 % in u and 19 % in tke at F = 0.7.

**Limits of the case.**
- The decks keep their own closure settings except the keys above; they do not use the convection blend of the
  recommended setting, which acts where the eddy-viscosity cap holds over terrain.
- The homogeneous case starts from a profile here and stops at the verdict after 293 iterations; its own test
  starts from the column first guess. The wind at a quarter of the canopy height is 0.158 m/s, the balance
  √(F/f) of that case's README.
- The hill has 2.3 m layers at the canopy top; a finer mesh was not run.
- `canopy.heat = fixed` is refused over an immersed ground and on a hybrid mesh.
- No restart of a run with a canopy over an immersed ground is tested here; the canopy's and the ground's
  checkpoint records are each tested on their own.
- The limits of the cut-cell rule (slopes, the drag on w above a cut cell) are in the
  [docs](../../../docs/canopy.rst).

## CTest
| Test | Label | What it checks |
|---|---|---|
| `rt_immersed_canopy_homogeneous` | regression | The homogeneous case: `--limits 0.5 1 1`. About 30 s on 4 cores. The limits are 1.5 to 2.2 times the measured differences. |
| `rt_immersed_canopy_edge` | nightly | The edge case: `--limits 3 13 5`. About 15 minutes. |
| `rt_immersed_canopy_hill` | nightly | The hill case: `--limits 5 7 28`. About 20 minutes. |
| `abort_immersed_canopy_heat` | regression, smoke | `canopy.heat = fixed` over an immersed ground is refused. |

The unit test `Canopy.StandsOnTheImmersedGround` checks the drag field cell by cell and the column integral for
flat grounds at four cut fractions and for a sloping ground (a map, and a round patch against the map), and the
column count of the rough-ground warning, on 1, 2 and 4 ranks.
