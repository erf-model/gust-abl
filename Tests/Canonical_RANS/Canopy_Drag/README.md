# Canopy_Drag

## Purpose
The forest canopy's drag on the momentum equations (`canopy.drag_diagonal`, `docs/canopy.rst`). Neutral flow driven
by a geostrophic wind of 10 m/s over flat ground, through two ERF-format forest patches in a periodic domain: a
20 m patch with the uniform profile and LAI 3, and a 25 m patch with the Lalic–Mihailovic profile and LAI 6 (the dense
end, C_d LAI = 1.2). The check runs the case with each implicit diagonal of the drag and once with the drag
explicit, as in ERF:

- the `derivative` and `picard` diagonals reach the same steady state (the diagonal changes the pseudo-time path,
  not the fixed point);
- the drag slows the wind inside the dense patch;
- with the drag explicit (`none`) the iteration departs within 20 iterations: in the dense patch the drag rate
  f |U| exceeds 1/Δτ (the limit of an explicit drag on the component along the wind; 2/Δτ across it) already at
  the starting CFL of 1;
- grid sequencing over two levels (the coarse level takes the restricted drag field) reaches the same steady state;
- a restart continues bit for bit (the canopy and its face drag are rebuilt from the inputs), and a restart with
  another `canopy.drag_diagonal` is accepted.

A second check, `check_canopy_tke.py`, runs the case with the canopy terms of the k–L closure on (`canopy.tke`, the
default) and off: the length cap L <= alpha/f holds in every canopy cell, k at the crown of the dense patch is near
and the k terms cut the crown's eddy viscosity. `rt_canopy_drag` runs the drag checks with `canopy.tke = off` (the
drag alone, as before the k terms).

## Files
| File | Role |
|---|---|
| `inputs_canopy_drag` | 16 x 16 x 40 cells on 1280 m x 1280 m, a 2 m first layer stretched by 1.08, four boxes, periodic sides; converged well past the default verdict (`steady.residual_tol = 1e-8`, 7 orders), so two runs can be compared at their fixed point. `canopy.file` and `init.first_guess` are given on the command line. |
| `canopy_patches.txt` | The two patches, ERF's `type x y height diameter cd lai laimax`. |
| `check_canopy_drag.py` | Runs the drag runs below in their own folders and checks them; the exit code is the verdict. |
| `check_canopy_tke.py` | Runs `on/` and `off/` (`canopy.tke`) to the verdict with the turbulence plotfile and checks the length cap, the crown balance and the crown's eddy viscosity. |
| `inputs_canopy_column` | The deck of `inputs_canopy_drag` without its canopy line, for `check_canopy_column.py`. |
| `check_canopy_column.py` | Writes a uniform map canopy (20 m, LAI 3, Cd 0.25 on every column) and runs `column/` (`init.first_guess = column`, with the table `column.txt`), `profile/` and `patches/` (the two-patch canopy with the column guess, refused); checks that the column's first residuals are below 2e-3 of the profile's (1e-9 measured), that the converged 3-D column is within 2e-3 of the table (1.3e-3), that the column run takes at most 0.7 of the profile run's iterations (1082 against 2046), and the refusal (UT-8.5.1, docs/initial_conditions.rst). |

The runs: `derivative/`, `picard/` and `sequenced/` to steady; `none/` for 20 iterations; `first/` for 40
iterations, then continued from its checkpoint to 80; `other/` continued from the same checkpoint with `picard`;
`straight/` for 80 iterations.

## Data and acknowledgements
No external data: the patches are made up for the test. The patch format is ERF's `erf.forest_file`
(`Source/SourceTerms/ERF_ForestDrag.cpp`, erf-model/ERF), and the sink is ERF's (`ERF_MakeMomSources.cpp`).

## How to run
```bash
ctest --test-dir <build> -R rt_canopy_drag --output-on-failure
```
or by hand, from an empty folder:
```bash
python3 <repo>/Tests/Canonical_RANS/Canopy_Drag/check_canopy_drag.py --tol 1e-6 --max-ratio 0.6 -- \
    mpiexec -n 2 <build>/gust <repo>/Tests/Canonical_RANS/Canopy_Drag/inputs_canopy_drag \
    canopy.file=<repo>/Tests/Canonical_RANS/Canopy_Drag/canopy_patches.txt
```

## Expected Results
Measured on macOS arm64 (AppleClang, two ranks, double precision), `check_canopy_drag.py` with `canopy.tke = off`
(the CTest's setting; with the k terms on the same checks pass with 5586 and 5386 iterations and a wind ratio of
0.56):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict, `derivative` / `picard` | 5710 / 5160 | both steady |
| `derivative` vs `picard`, max difference / max value | u 2.1e-10, v 2.6e-9, w 8.2e-10, p 1.4e-8, k 8.4e-10 | 1e-6 |
| Grid sequencing (two levels) vs `derivative` | u 9.2e-11, v 5.7e-11, w 7.1e-11, p 1.2e-8, k 1.5e-10 | 1e-6 |
| u at 5.3 m in patch 2 / on open ground | 1.14 / 2.16 m/s (0.53) | at most 0.6 |
| r after 20 iterations, `none` / `derivative` | 3.1e3 / 0.16 | `none` above 1 and above `derivative` |
| Restart: continued vs straight, every field | 0 | bit for bit |
| Restart with `canopy.drag_diagonal = picard` | accepted | |

`check_canopy_tke.py` (`--nut-ratio 0.1`):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict, k terms on / off | 5586 / 5710 | both steady |
| Canopy cells with L above alpha/f (on) / at it / above it in off | 0 / 176 / 176 | 0 / > 0 / > 0 |
| Crown eddy viscosity, on / off | 0.033 / 3.99 m^2/s (0.008) | at most 0.1 |

At the default verdict (r dropped 4 orders) with both patches at one LAI, the iterations are (LAI 1, 3, 6, 10):
`derivative` 2312, 1737, 1820, 1831; `picard` 2299, 1655, 1855, 2017. Without a canopy the same case takes 3362.

## CTest
`rt_canopy_drag`, `rt_canopy_tke` and `rt_canopy_column` (label `regression`, double precision only; about 4, 3 and
3 minutes on two ranks).
