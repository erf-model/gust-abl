# Unsteady_Hill_2D

## Purpose
RT-6.2.1 (plan rev. 6, sub-phase 6.2, beyond the plan): the 2-D hill of `Neutral_Hill_2D` (ERF's Neutral_Hill_2D: the
Witch of Agnesi ridge, 100 m high, half-width 500 m, 128 x 1 x 64 cells, the surface layer, Coriolis and the
geostrophic wind) run as unsteady RANS by dual time stepping (`steady.stage = unsteady`, docs: Solution algorithm,
Time driver). Ten backward-Euler steps of 60 s from the log-law first guess to t = 600 s; each step converges two
orders below its starting residual by SIMPLEC pseudo-time iterations (at most 60). The test checks that the time
driver gives the same steps and the same answer on 1, 2 and 4 ranks.

## Files
| File | Role |
|---|---|
| `inputs_unsteady_hill2d` | The deck: `inputs_steady_hill2d_restart` with `steady.stage = unsteady`, the inner-stop keys and no checkpoint; the script gives `unsteady.stop_time` and `unsteady.dt`. |
| `check_unsteady_hill2d.py` | Runs the deck on 2, 1 and 4 ranks and one short step; checks the step log and the answers. |

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

The case checks the time driver's bookkeeping and its rank independence. It does not compare with measurements or a
closed form; the first-order accuracy in dt is checked against Stokes' first problem in the unit test
`TimeDriver.StokesFirstProblemFirstOrder`.

## CTest
`rt_unsteady_hill2d` (label `regression`, MPI builds only).
