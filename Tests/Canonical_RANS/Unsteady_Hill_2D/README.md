# Unsteady_Hill_2D

## Purpose
RT-6.2.1 and RT-6.3.1 (plan rev. 6, sub-phases 6.2 and 6.3, beyond the plan): the 2-D hill of `Neutral_Hill_2D` (ERF's Neutral_Hill_2D: the
Witch of Agnesi ridge, 100 m high, half-width 500 m, 128 x 1 x 64 cells, the surface layer, Coriolis and the
geostrophic wind) run as unsteady RANS by dual time stepping (`steady.stage = unsteady`, docs: Solution algorithm,
Time driver). Ten backward-Euler steps of 60 s from the log-law first guess to t = 600 s; each step converges two
orders below its starting residual by SIMPLEC pseudo-time iterations (at most 60). The test checks that the time
driver gives the same steps and the same answer on 1, 2 and 4 ranks, and that the unsteady outputs (running averages
and variances, probes, plotfiles and checkpoints every N steps) are right and a restart continues the run bit for bit
(docs: Solution algorithm, Unsteady outputs).

## Files
| File | Role |
|---|---|
| `inputs_unsteady_hill2d` | The deck: `inputs_steady_hill2d_restart` with `steady.stage = unsteady`, the inner-stop keys; the script gives `unsteady.stop_time`, `unsteady.dt` and the output keys. |
| `probes_hill2d.txt` | Four probes (`x y h`, h above the ground): the hill top at 50 m, upstream at 10 m, the lee at 30 m and the domain's high corner at 500 m. |
| `probes_outside_hill2d.txt` | One probe outside the domain, for the CTest `abort_unsteady_probe_outside`. |
| `probes_restart_hill2d.txt` | One probe (the hill top), for the CTest `abort_unsteady_restart_probes`. |
| `check_unsteady_hill2d.py` | Runs the deck on 2, 1 and 4 ranks, one short step and a restart; checks the step log, the answers and the outputs. |

## Data and acknowledgements
No external data. The case is ERF's Neutral_Hill_2D (`extern_erf/`, ERF developers).

## How to run
```bash
ctest --test-dir <build> -R rt_unsteady_hill2d
```

## Expected Results
- 2 ranks: ten steps of 60 s, the last ending at t = 600 s, none stopped at `unsteady.max_inner`, the summary line, and
  the time 600 s in the plotfile's header.
- 1 and 4 ranks: the same inner iteration count in every step; u, v, w within 1e-9 of the largest |u| and p, θ, k
  within 1e-9 of their own largest values, against 2 ranks.
- u at the hill top, 50 m above the ground, differs from the first guess by more than 0.1 m/s at t = 600 s.
- Outputs (RT-6.3.1): the probe file has 11 rows (t = 0 and one per step); plotfiles at t = 300 s and 600 s carry
  `u_mean ... tke_mean` and `u_var ... tke_var`, equal at each probe's cell to the two-pass weighted mean and variance
  of that probe's rows, with weights 30 s for step 2 and 60 s after it (averages from t = 90 s, mid-step 2). The
  unit test `UnsteadyOutput.AveragesAreTheWeightedSeries` checks every cell.
- Restart: from the 2-rank checkpoint of step 5 the run continues at t = 300 s with steps 6 to 10, the same inner
  counts, a byte-identical final plotfile and probe file, and a final checkpoint identical except the driver history's
  wall-clock columns.

The case checks the time driver's bookkeeping and its rank independence. It does not compare with measurements or a
closed form; the first-order accuracy in dt is checked against Stokes' first problem in the unit test
`TimeDriver.StokesFirstProblemFirstOrder`.

## CTest
`rt_unsteady_hill2d` (label `regression`, MPI builds only). The checkpoint fixture `rt_unsteady_hill2d_chk` and the
restart refusals `abort_unsteady_restart_average_start`, `abort_unsteady_restart_past_stop` and
`abort_unsteady_restart_probes` use the same deck.
