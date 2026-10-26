# Kestrel_Dryrun

## Purpose
The Kestrel scripts (`tools/kestrel/`) cannot run here. This case checks that they parse, and runs them with stub
`module`, `sbatch` and `srun` commands to exercise the submit and archive paths.

## Files
| File | Role |
|---|---|
| `kestrel_dryrun.sh` | `kestrel_dryrun.sh SRC REFAVG PLOTFILE WORK`: builds a fake ERF tree and finished run in `WORK`, then submits and archives it. |
| `nightly_gpu_dryrun.sh` | `nightly_gpu_dryrun.sh SRC WORK`: runs `tools/kestrel/nightly_gpu.sh` with stub `git`, `module`, `sbatch`, `srun`, `cmake` and `ctest`. |

## CTest
`kestrel_parse_*` (one per script, `bash -n`), `kestrel_dryrun` (needs `smoke_run` for its plotfile),
`kestrel_nightly_dryrun` (label `unit`, no build products).

## Expected Results
- The Askervein job script has the expected run length, plot period, rank count and no 2-D plotfiles.
- The archive averages only the 3-D plotfiles in the final window and keeps the inputs but not the run's outputs.
- A second case with the default prefix `plt` keeps its `plt2d_` plotfiles out.
- The nightly GPU job refuses to run outside an allocation; a passing run clones erf-model/gust-abl at `main`,
  builds with `build_kestrel.sh gpu`, configures one rank per GPU (`--gpus-per-task=1`), runs `ctest -L smoke` with a
  JUnit file, writes `pass` and prunes old run logs; a smoke run of `gust` without the GPU runtime flags stops the job
  before the tests (`fail: arena flags`, naming the test); failing tests give `fail: tests` and exit 1; `submit` needs
  `GUST_ACCOUNT` and its job script requests the GPUs; a bad `GUST_NIGHTLY_GPUS` is refused.

## How to run
```bash
ctest --test-dir <build> -R "^kestrel_"
```
