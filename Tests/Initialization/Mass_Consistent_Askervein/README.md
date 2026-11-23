# Mass_Consistent_Askervein

## Purpose
The mass-consistent initial field on the Askervein mesh (30 x 30 x 18): a log-law first guess (u* = 0.6 m/s,
z0 = 0.03 m, along x), inflow at xlo, outflow at xhi, alpha_h = 1, alpha_v = 0.5. The divergence drops from
4.05e-1 to 6.1e-12 1/s.

## Files
| File | Role |
|---|---|
| `inputs_init_askervein` | The deck. Terrain file on the command line: `terrain.file=extern_erf/askervein.txt`. |
| `Gold/` | `init_askervein.json` (iterations, divergence, extrema), `init_askervein_plt00000/` (the state). |

Check script: `Supporting_Files/check_init.py askervein` (shared with `Bump_Speedup_Trend`).

## Data and acknowledgements
The terrain is ERF's Askervein terrain of `Exec/CanonicalTests/Real_Terrain/Askervein` (`extern_erf/askervein.txt`,
recorded in `extern_erf/EXTRACTED.md`; BSD 3-clause with the LBNL clause, `extern_erf/LICENSE.md`). Thanks to the ERF
developers. The file is redistributed in this repository; nothing has to be downloaded.

## CTest
`rt_init_askervein`, `rt_init_askervein_matches` (double precision only, two ranks).

## Expected Results
- Iterations at most 10 % above the gold's 13; divergence after at most 10 x `solver.poisson_reltol` times before;
  divergence before at the gold's printed 10 digits.
- Maximum |w| and the largest horizontal speed at the gold's values (1e-6); theta = 300 K everywhere; p' = 0;
  k > 0.
- The state matches `Gold/init_askervein_plt00000` at relative 1e-8. The gold's k is for the 100 m default of
  `turb.max_geom_lscale` (rewritten in November 2026; a run without the SIMPLEC stage refuses `turb.*`). The other
  fields did not change with the default.

lambda is not in the plotfile: the state plotfile is the regression artefact, and the unit tests pin lambda.

## How to run
```bash
ctest --test-dir <build> -R "^rt_init_askervein"
```
