# Nests_Two_Way

## Purpose
Two-way nesting (PR-70, `nest.feedback = defect`; and PR-71, `nest.feedback = relax`, below). The outer domain and its nest first run one way, as in
`Nests_One_Process` except that the nest's levels already follow the joint rule (they bend to the outer domain's
surfaces aloft), so "one way" in the tables below is that first pass. Then they run jointly: every iteration the nest's faces are refreshed from the outer domain, and
the outer domain carries the defect source of the nest over the nest's footprint
(`docs/boundary_conditions.rst`, Two-way nests). The case checks that the joint run moves both domains towards the
uniform fine run:
1. the joint nest is closer to the uniform 25 m run than the one-way nest;
2. the outer domain outside the footprint, where the nest is not, is closer too: the nest's flow reaches the rest of
   the outer domain.

The hill is steeper than that of `Nests_One_Process` (100 m high, half-width 150 m, against 60 m and 250 m), so that
the 50 m outer cells miss a part of the flow the 25 m nest resolves.

## Files
| File | Role |
|---|---|
| `inputs_outer` | The outer domain: 32 x 24 x 16 cells of 50 m over 1600 m x 1200 m x 800 m (first layer 10 m), the round Witch of Agnesi hill in the middle, periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m, 20 degrees) on xlo, outflow on xhi, neutral, Rayleigh layer above 400 m on w and theta only (as every nest). No nest keys: the check script appends them. |
| `inputs_abort` | `inputs_outer` without `mesh.max_grid_size` and `steady.stage`, for the abort tests. |
| `check_nests_two_way.py` | Runs the uniform 25 m reference (the deck with 64 x 48 x 16 cells) and the joint run, and compares them (below). |
| `check_joint_restart.py` | Runs a joint run straight, the same run stopped after `--stop` exchanges, and its restart from the joint record, and compares them byte for byte (Joint restart, below). |

The nest `hill` covers outer cells 8 to 24 in x and 6 to 18 in y (`nest.hill.lo = 400 300`,
`nest.hill.hi = 1200 900`): 32 x 24 x 16 cells of 25 m. The second joint run has two nests that touch at
x = 1000 m: `crest` (`300 300` to `1000 900`, the windward side and the crest) and `lee` (`1000 250` to `1500 950`).

## Data
None. The case is analytic.

## How to run
```bash
ctest --test-dir <build> -R rt_nests_two_way
```
By hand, in a scratch directory:
```bash
python3 check_nests_two_way.py -- mpiexec -n 2 <build>/gust inputs_outer
```
The script makes `ref/` (the 25 m run), `joint/` (the outer domain and the nest `hill`: one way, then jointly) and
`joint2/` (the same with the nests `crest` and `lee`).

## Expected Results
`check_nests_two_way.py` requires:
1. both runs exit 0; the reference reaches the steady verdict; the joint log has the joint run's start line, the
   largest |S| per exchange and `joint run: steady after N exchanges`;
2. two plotfiles each for the outer domain and the nest: the one-way ones, then the joint ones;
3. against the reference, the RMS of u:
   - over the nest (on the 25 m cells): joint at most 0.5 of one way;
   - over the outer cells outside the footprint (the reference averaged over 2 x 2 cells): joint at most 0.6 of one
     way;
   - with the two nests: over each nest, joint at most 0.7 of one way; outside both footprints, at most 0.9.

Measured on 2 ranks (RMS of each field against the 25 m run, m/s; k in m^2/s^2):

| | nest, one way | nest, joint | outer outside, one way | outer outside, joint | footprint, one way | footprint, joint |
|---|---|---|---|---|---|---|
| u | 8.23e-3 | **2.03e-3** | 1.44e-2 | **4.89e-3** | 2.89e-2 | **1.32e-2** |
| w | 3.21e-3 | **7.64e-4** | | | 8.38e-3 | **7.90e-3** |

The joint run reaches steady after 208 exchanges.

With the two nests (u unless stated):

| | nest, one way | nest, joint | footprint, one way | footprint, joint |
|---|---|---|---|---|
| crest | 6.11e-3 | **2.56e-3** | 2.66e-2 | **1.45e-2** |
| lee | 3.00e-2 | **1.82e-2** | 3.00e-2 | **1.79e-2** |
| lee, w | **6.42e-3** | 1.27e-2 | **6.72e-3** | 1.19e-2 |

Outside both footprints u moves from 8.53e-3 to **6.83e-3**; the run is steady after 254 exchanges. The lee nest's
w is further from the reference after the joint run: over the footprint the outer domain carries the w of the
restricted contravariant flux, not a cell mean of the nest's w, and in the hill's wake the two differ most.

Where the nest's downstream face stands matters. With the crest nest alone, ending at x = 800, 900, 1000, 1100 or
1200 m (crest at 800 m, half-width 150 m), the joint u error over the nest is 1.54, 3.31, 0.40, 0.27 and 0.22 times
the one-way one, and outside the footprint 1.59, 2.74, 0.66, 0.40 and 0.34 times. A face on the lee slope, where the
flow changes most, takes the 50 m solution there as its boundary values, and the joint run then moves both domains
away from the fine run. The two-nest layout therefore ends the crest nest 1.33 half-widths behind the crest.

## Feedback by relaxation (`--relax`)
`check_nests_two_way.py --relax` (CTest `rt_nests_relax`, PR-71) runs the reference, the defect run with the nest
`hill`, and the joint runs with `nest.feedback = relax` at `nest.feedback_factor` = 1, 4 and 16 (`relax1/`,
`relax4/`, `relax16/`; `docs/boundary_conditions.rst`, Two-way nests by relaxation). It requires:
1. every run exits 0 (the two-nest `joint2/` run is skipped), and each relax log's start line names its factor;
2. the outer domain minus the restricted nest where the relaxation acts (the log's last RMS of u) falls from factor 1
   to 4 to 16, and at 16 is at most 0.35 of its value at 1;
3. at factor 1, against the 25 m run, the RMS of u over the nest at most 0.5 of one way and over the outer cells
   outside the footprint at most 0.6.

Measured on 2 ranks (RMS of u, m/s; "outer - nest" where the relaxation acts):

| run | exchanges | outer - nest | nest | outer outside |
|---|---|---|---|---|
| one way | | | 8.23e-3 | 1.44e-2 |
| defect | 208 | | **2.03e-3** | **4.89e-3** |
| relax, f = 0.25 | 150 | 1.33e-2 | 2.32e-3 | 5.90e-3 |
| relax, f = 0.5 | 154 | 1.07e-2 | 2.34e-3 | 5.10e-3 |
| relax, f = 1 | 160 | 8.04e-3 | 3.38e-3 | 5.03e-3 |
| relax, f = 2 | 159 | 5.52e-3 | 5.60e-3 | 6.14e-3 |
| relax, f = 4 | 177 | 3.49e-3 | 9.94e-3 | 9.54e-3 |
| relax, f = 16 | 508 | 1.50e-3 | 3.27e-2 | 2.94e-2 |

(0.25, 0.5 and 2 are not in the test.) The difference between the two domains falls with the factor, but above
f = 1 both move away from the fine run: the relaxation's target carries no nest pressure gradient, so as it tightens
the outer velocity under the nest stops responding to the outer pressure, and the nest, which takes its faces from
it, follows. The exchanges grow with f for the same reason: the flow rate through the footprint is set by a loop
(the nest returns the flow rate it is given, the relaxation imposes it) whose gain tends to 1. The bounds on the
errors therefore apply at f = 1 only.

## Joint restart (`check_joint_restart.py`)
A joint run that stops after `steady.max_iterations` exchanges writes a joint record (`Joint`) in the outer
domain's last checkpoint and names that checkpoint in its error message (`docs/boundary_conditions.rst`, Two-way
nests, Restart). The script runs, each in its own directory:
1. `straight/`: the deck with the nests `crest` and `lee` and `nest.feedback = defect` (with `--relax`: the nest
   `hill` and `nest.feedback = relax`), to the steady verdict;
2. `stopped/`: the same with `steady.max_iterations = M` (`--stop M`), which stops after M exchanges;
3. `stopped/`: the restart, with the deck less its `init.*` lines and `io.restart` set to the named checkpoint.

It requires:
1. the straight run and the restart exit 0, the stopped run does not; the record has the feedback, M exchanges and
   one line per nest in `nest.names` order, each naming a checkpoint that exists;
2. the restart logs `restarted after exchange M`, and is steady after the straight run's exchanges with every
   domain's iterations over the run the same;
3. the last plotfile of every domain and the last checkpoint's state, node heights, closure history, nest faces and
   joint record are byte for byte the straight run's (the driver record holds wall times and is not compared).

M must be above the one-way runs' iterations (`steady.max_iterations` caps them too: 99 for the outer domain, up to
121 for a nest here) and below the straight run's exchanges.

| | ranks | M | exchanges | result |
|---|---|---|---|---|
| crest and lee, defect | 3 | 150 | 254 | byte for byte |
| hill, relax | 2 | 140 | 160 | byte for byte |

## CTest
`rt_nests_two_way` (2 ranks, labels `regression;smoke`). `rt_nests_joint_restart` (3 ranks, labels
`regression;smoke`) and `rt_nests_joint_restart_relax` (2 ranks, label `regression`) run `check_joint_restart.py`
(above). The abort tests `abort_nests_joint_restart_names`, `_iterations`, `_version` and `_missing` check that a
joint record naming other nests, a `steady.max_iterations` at or below its exchanges, a record of another version,
and a record naming a nest checkpoint that is not there stop the run and name `io.restart`. The abort tests `abort_nests_feedback_value`,
`abort_nests_feedback_no_nests`, `abort_nests_feedback_newton`,
`abort_nests_feedback_damping`, `abort_nests_feedback_odd_boxes` and `abort_nests_feedback_narrow` check that the inputs `nest.feedback` refuses stop
the run and name the key. `inputs_abort` is `inputs_outer` without `mesh.max_grid_size` and `steady.stage`, which
two of them set on the command line (a key given twice stops the run).

`rt_nests_relax` (2 ranks, label `regression`; above): about 2 minutes on an Apple M-series laptop, against under
1 minute for `rt_nests_two_way`, most of it the f = 16 run (508 exchanges). `--relax` skips the two-nest `joint2/`
run. `abort_nests_feedback_factor`, `_factor_nan`, `_factor_unused` and `_relax_damping` check that
`nest.feedback_factor = 0`, `nan`, a factor with `nest.feedback = defect`, and `relax` with a Rayleigh layer on u, v
and k stop the run and name the key.
