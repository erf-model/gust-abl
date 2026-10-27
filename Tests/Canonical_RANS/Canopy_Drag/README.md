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

## Files
| File | Role |
|---|---|
| `inputs_canopy_drag` | 16 x 16 x 40 cells on 1280 m x 1280 m, a 2 m first layer stretched by 1.08, four boxes, periodic sides; converged well past the default verdict (`steady.residual_tol = 1e-8`, 7 orders), so two runs can be compared at their fixed point. `canopy.file` and `init.first_guess` are given on the command line. |
| `canopy_patches.txt` | The two patches, ERF's `type x y height diameter cd lai laimax`. |
| `check_canopy_drag.py` | Runs the five runs below in their own folders and checks them; the exit code is the verdict. |

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
Measured on macOS arm64 (AppleClang, two ranks, double precision):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict, `derivative` / `picard` | 5710 / 5160 | both steady |
| `derivative` vs `picard`, max difference / max value | u 2.1e-10, v 2.6e-9, w 8.2e-10, p 1.4e-8, k 8.4e-10 | 1e-6 |
| Grid sequencing (two levels) vs `derivative` | u 9.2e-11, v 5.7e-11, w 7.1e-11, p 1.2e-8, k 1.5e-10 | 1e-6 |
| u at 5.3 m in patch 2 / on open ground | 1.14 / 2.16 m/s (0.53) | at most 0.6 |
| r after 20 iterations, `none` / `derivative` | 3.1e3 / 0.16 | `none` above 1 and above `derivative` |
| Restart: continued vs straight, every field | 0 | bit for bit |
| Restart with `canopy.drag_diagonal = picard` | accepted | |

At the default verdict (r dropped 4 orders) with both patches at one LAI, the iterations are (LAI 1, 3, 6, 10):
`derivative` 2312, 1737, 1820, 1831; `picard` 2299, 1655, 1855, 2017. Without a canopy the same case takes 3362.

## CTest
`rt_canopy_drag` (label `regression`, double precision only; about 4 minutes on two ranks).
