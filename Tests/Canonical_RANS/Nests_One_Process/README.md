# Nests_One_Process

## Purpose
Nests run in the same process as the outer domain. The inputs give only each nest's name and corners:
```
nest.names   = hill edge
nest.hill.lo = 400 300
nest.hill.hi = 1200 900
nest.edge.lo = 1200 300
nest.edge.hi = 1500 600
```
Every other setting is the outer deck's, or derived from it (`Source/Nests.H`, `docs/inputs.rst`): refinement 2 2 1, the
four lateral faces `nest` with the outer run's final checkpoint as parent, `init.first_guess = nest`, a Rayleigh layer
on w and theta only, and output names with `_<name>` appended. The outer domain runs to steady and writes its
checkpoint; then each nest runs from it, in the order of `nest.names`.

The case checks four things:
1. the nests change nothing in the outer domain: its plotfile is bit for bit that of the same deck without nests;
2. a nest is the run its derived settings describe: its plotfile is bit for bit that of a separate run of
   `inputs_nest` (the same nest written as a deck of its own) from the outer run's final checkpoint;
3. a nest stands on the outer ground: the hill stays at the outer domain's centre (800 m, 600 m), so the nest's ground
   equals the outer ground on the nodes they share;
4. a nest refines the outer solution: its u, averaged over the 2 x 2 nest cells of each outer cell, is close to the
   outer u there.

The outer domain is the hill parent of `Nest_Faces` with a Rayleigh layer on every field above 400 m, so that the
nests must drop the layer on u, v and k.

## Files
| File | Role |
|---|---|
| `inputs_outer` | The outer domain: 32 x 24 x 16 cells of 50 m over 1600 m x 1200 m x 800 m (first layer 10 m), a 60 m round Witch of Agnesi hill (half-width 250 m) in the middle, periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m, 20 degrees) on xlo, outflow on xhi, neutral, Rayleigh layer above 400 m. No nest keys: the check script appends them. |
| `inputs_nest` | One nest as a deck of its own: the settings a nest derives from `inputs_outer`. The mesh extent, `bc.nest_parent` and the output names come from the command line. |
| `check_nests_one_process.py` | Runs the outer deck with and without nests and each nest separately, and compares the plotfiles (below). |

The two nests: `hill`, outer cells 8 to 24 in x and 6 to 18 in y (32 x 24 x 16 cells of 25 m), and `edge`, outer
cells 24 to 30 in x and 6 to 12 in y (12 x 12 x 16 cells). They touch along x = 1200 m; `edge` is 2 outer cells
from the outflow face, the smallest margin allowed.

## Data and acknowledgements
None. The case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "nests_"
```
By hand, in a scratch directory:
```bash
python3 check_nests_one_process.py <build>/amrex_fcompare -- mpiexec -n 2 <build>/gust inputs_outer
```
The script makes `joint/` (the outer deck with the nests), `outer/` (the outer deck alone) and `sep_<name>/` (each
nest from `inputs_nest`).

## Expected Results
`check_nests_one_process.py` requires:
1. each run reaches the steady verdict: once in `outer/`, three times in `joint/` (the outer domain and each nest);
2. the joint log names each nest after the outer run, in the order of `nest.names`, with its corners and the outer
   run's final checkpoint;
3. `amrex_fcompare` with zero tolerances finds no difference between the joint run's outer plotfile and that of
   `outer/`, and between each nest's final plotfile and that of `sep_<name>/`;
4. each nest writes `outer_plt_<name>` and `outer_chk_<name>`;
5. each nest's ground (`z_phys_nd` at k = 0, written next to the plotfile) equals the outer ground on the shared
   nodes to 1e-6 m;
6. max |<u>_nest - u_outer| over the nest's outer cells is below the bound in `DU_TOL`.

Measured (macOS arm64, AppleClang, double precision, 2 ranks):

| Domain | Cells | Iterations to steady | Joint vs separate | Ground vs outer (shared nodes) | max \|<u>_nest - u_outer\| (bound) |
|---|---|---|---|---|---|
| Outer | 32 x 24 x 16 | 84 | bit for bit | | |
| `hill` | 32 x 24 x 16 | 96 | bit for bit | 0 m (221 nodes) | 0.0959 m/s (0.12) |
| `edge` | 12 x 12 x 16 | 73 | bit for bit | 0 m (49 nodes) | 0.0115 m/s (0.02) |

With the hill centred on each nest instead (the behaviour of a one-way patch before this case), the `edge` nest's
ground is 50.3 m off the outer ground and the check fails.

The faults of the nest inputs stop the run before the outer domain starts, each naming the key: a corner off the
outer cell faces, a nest less than 2 outer cells from a non-periodic outer face, two nests overlapping, a missing
corner, an unknown `nest.<name>.*` key, `io.write_checkpoint = 0`, `steady.stage = none`, `simplec.theta = frozen`,
and a nest box narrower than `steady.line_overlap`. An outer domain that does not reach steady stops the run after the
outer domain, and no nest runs. A nest that does not reach steady does not stop the other nests; the run then
stops with an error that names it.

## CTest
| Test | Ranks | Checks |
|---|---|---|
| `rt_nests_one_process` | 2 | The checks above. |
| `abort_nests_misaligned` | 1 | `nest.a.lo = 410 300`: off the 50 m outer cell faces. |
| `abort_nests_margin` | 1 | `nest.a.lo = 50 300`: 1 outer cell from xlo. |
| `abort_nests_overlap` | 1 | Two nests sharing outer cells. |
| `abort_nests_missing_corner` | 1 | `nest.a.hi` not given. |
| `abort_nests_unknown_key` | 1 | `nest.a.ratio`: not an input. |
| `abort_nests_no_checkpoint` | 1 | `io.write_checkpoint = 0`: the nests need the outer run's checkpoint. |
| `abort_nests_stage` | 1 | `steady.stage = none`. |
| `abort_nests_frozen` | 1 | `simplec.theta = frozen`. |
| `abort_nests_layout_before_outer` | 1 | `nest.a` 1 outer cell wide (2 nest cells) with `steady.sweep = lusgs`, `steady.line_overlap = 4`: stops before the outer run. |
| `abort_nests_outer_not_steady` | 1 | `steady.max_iterations = 3`: the outer domain does not reach steady; no nest runs. |
| `abort_nests_nest_not_steady` | 1 | `steady.max_iterations = 90`: the outer domain (84) and `edge` (73) reach steady, `hill` (96) does not; both nests run, then the run stops naming `nest.hill`. |
| `abort_nests_stf` | 1 | Another deck (`Steep_Ridge_Thin_Cell`, a ridge with `terrain.smoothing` at its default 1): STF is refused with nests. |
| `rt_nests_outer_init_keys` | 1 | One nest with `init.z_agl_min = 5` (an outer first-guess key): the nest runs. |

The settings are not the recommended ones (`turb.realizability = 0.309`, `turb.nut_relax = 0.5`, the hybrid blend):
the case reuses the `Nest_Faces` hill parent and tests equality between runs, not the flow.
