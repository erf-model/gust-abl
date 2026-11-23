# Column_Inflow

## Purpose
`bc.inflow_profile = column` in one run, compared with the two runs it replaces. Before this option, a 3-D run with
inflow faces and a target wind needed two runs:

1. A periodic column run with `init.first_guess = column` and a target wind (`init.column_target_u`, `_v`,
   `_height`). It wrote its converged column (`io.column_guess_file`) and logged the adjusted geostrophic wind.
2. The 3-D run, with `bc.inflow_profile = file` on that table and the logged forcing copied by hand into
   `forcing.ug` and `forcing.vg`.

With `bc.inflow_profile = column`, the 3-D run solves the same column at start-up. It then uses that column as the
inflow faces' profile, and keeps the adjusted forcing. The case checks four things:

- **Same column.** The single run writes the two-run table to the last bit (every value in `%.17g`), and its log
  prints the same forcing.
- **Checkpoint record.** The checkpoint holds that column, its roughness and its forcing.
- **Restart.** A restart from the checkpoint continues the run to the same answer, bit for bit.
- **Restart with no stage.** A restart that reads no `bc.*` keys still writes the record into its own checkpoint.

The 3-D deck is a round Agnesi hill (100 m high, 400 m half-width) on a 3.2 km square domain, 16 x 16 x 32 cells
on four boxes. Inflow is on x-lo and y-lo, outflow on x-hi and y-hi, and the top is slip. The surface roughness is
0.05 m. A Rayleigh layer damps w above 1200 m. There is no gravity, and Coriolis is on at 45 deg with
`abl.drive = geostrophic_inflow`. The target is 8 m/s at 100 m from the south-west (5.657, 5.657 m/s). The column
deck is the same column on a 4 x 4 x 32 periodic mesh with the same levels.

What "same" means here: the table and the forcing are bit-equal between the two ways. A true two-run 3-D result
copies the forcing from the log, which prints 9 significant digits, so its 3-D fields differ from the single run's
by about 1e-11 (relative). The single run keeps the forcing in full precision.

## Files
| File | Role |
|---|---|
| `inputs_column` | First of the two runs: the periodic column with the target wind; writes `two_run_inflow.txt` and logs the forcing. |
| `inputs_inflow_column` | The single run: the hill with `bc.inflow_profile = column`; writes `one_run_inflow.txt`, logs the forcing, and writes checkpoints. `steady.max_iterations` is given on the command line. |
| `inputs_inflow_column_restart` | The single run's deck for a restart: no `init.*` keys and no `io.column_guess_file` (a restart refuses them), with its own output names. |
| `inputs_inflow_column_nostage` | A restart with no stage: mesh, terrain and output names only. It writes a plotfile and a checkpoint of the state it read. |
| `check_inflow_column.py` | `same-column` compares the two tables value for value and the forcing keys both logs print. `record` checks a checkpoint's `InflowColumn` record against the table, the log and z0. `drop-record` copies a checkpoint without its record. |

## Data and acknowledgements
No external data.

## How to run
From a build directory (`gust` built), in an empty folder:
```bash
mpiexec -n 2 ../gust ../Tests/Canonical_RANS/Column_Inflow/inputs_column > two.log 2>&1
```
```bash
mpiexec -n 2 ../gust ../Tests/Canonical_RANS/Column_Inflow/inputs_inflow_column steady.max_iterations=4 > one.log 2>&1
```
```bash
python3 ../Tests/Canonical_RANS/Column_Inflow/check_inflow_column.py same-column two_run_inflow.txt two.log one_run_inflow.txt one.log
```
```bash
python3 ../Tests/Canonical_RANS/Column_Inflow/check_inflow_column.py record chk00004 one_run_inflow.txt one.log 0.05
```
The column run ends "not steady after 1 iterations" and the single run "not steady after 4 iterations"; both are
expected. Each run takes about a second.

## Expected Results
The two tables (33 rows) are equal value for value. Both logs print the same forcing: `forcing.ug = 9.26889497`,
`forcing.vg = 5.13447251`. The single run's log names the faces' column:
`gust: inflow faces: the converged column (bc.inflow_profile = column), 33 rows, z0 0.05 m, with the target wind's forcing`.

A restart from the checkpoint of iteration 2, run to iteration 4, gives a plotfile equal to the uninterrupted run's
at iteration 4 (fcompare with zero tolerances). The restart must use the same rank count. GUST's results change with
the rank count by about 1e-5 (relative) after four iterations, with any inflow profile. That dependence was there
before this option.

The same comparison on the Perdigão deck of `Tests/Real_Terrain/Sites_128m` (z0 from `init.column_z0` under its
roughness map) also gave the two-run table to the last bit, with forcing 15.3874171 / 6.43276768 (run locally,
2026-11-23; not a CTest).

Limits, shared with `bc.inflow_profile = file`:

- Below the first cell centre, the faces interpolate linearly from the table's z = 0 row; they do not follow the
  column's log law there.
- Nothing checks that the column's wind, turned by Coriolis, enters through the faces set to `inflow`.
- The column's surface heat flux is not in the record. A stratified restart takes the heat flux from its deck.

## CTest
All tests have the labels `regression;smoke` and run on 1 or 2 ranks, about 5 s in all.

| Test | Check |
|---|---|
| `rt_column_inflow_clean` | Makes the run folders. |
| `rt_column_inflow_two_run` | The periodic column run. |
| `rt_column_inflow_one_run` | The single run, 4 iterations; its log names the column (33 rows, z0 0.05 m). |
| `rt_column_inflow_same_column` | `same-column`: the tables are equal to the last bit, and so is the forcing. |
| `rt_column_inflow_record` | `record`: the checkpoint `chk00004` holds the table, z0 and forcing. |
| `rt_column_inflow_first`, `rt_column_inflow_restart` | 2 iterations, then a restart to 4 ("4 over the run"). |
| `rt_column_inflow_restart_matches` | fcompare with zero tolerances: the restart's plotfile at 4 equals the single run's. |
| `rt_column_inflow_drop_record`, `abort_column_inflow_no_record` | A checkpoint without the record is refused, naming `io.restart`. |
| `rt_column_inflow_nostage`, `rt_column_inflow_nostage_record` | A restart with no stage from `chk00004`; its checkpoint `nostage_chk00004` holds the same record. |
| `abort_column_inflow_no_face` | `bc.inflow_profile = column` on the periodic deck is refused ("needs an inflow face"). |

Mutants: each of these makes a test fail.

- The faces are not given the column (`use_inflow_column` skipped after the column guess): `rt_column_inflow_one_run`
  and `rt_column_inflow_first` abort.
- The restart does not take the record's geostrophic wind: `rt_column_inflow_restart_matches` fails.
- The record is written in 9 digits instead of 17: `rt_column_inflow_record`, `rt_column_inflow_restart_matches` and
  `rt_column_inflow_nostage_record` fail.
- The no-stage restart does not write the record: `rt_column_inflow_nostage_record` fails.
- The faces take z0 from `bc.inflow_z0`'s default instead of the column's: `rt_column_inflow_one_run` and
  `rt_column_inflow_restart` fail.

The unit tests `InflowColumn.*` (`Tests/Unit_Tests/test_InflowColumn.cpp`) cover the checks, the record's round
trip and parser, the checkpoint record, the roughness, and the column's start from the analytic defaults, each with
its own mutant.
