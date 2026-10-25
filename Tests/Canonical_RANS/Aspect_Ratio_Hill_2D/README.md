# Aspect_Ratio_Hill_2D

## Purpose
How thin first cells change the convergence of the steady iteration. A 2-D Agnesi ridge (100 m high, 250 m
half-width, max ground slope 14.6 deg) with dx = dy = 100 m and the first layer dz1 = dx/ratio for dx/dz1 = 10 to
5000; the layers are stretched by 1.1 to the same top of about 1.5 km. Log-law inflow (u* 0.4 m/s, z0 0.01 m), outflow,
slip top with Rayleigh damping of w above 1000 m, no gravity, the basic terrain-following mapping, and the recommended
setting of the realizability sensor. The quantity is the number of iterations until the relative residual falls below
1e-4. Beyond the plan.

## Files
| File | Role |
|---|---|
| `inputs_aspect_ratio_hill` | The deck (the dx/dz1 = 10 member). |
| `aspect_ratio_sweep.py` | Copies the deck with `mesh.dz0` = dx/ratio and enough layers for the same top, runs each ratio, and prints the iterations to r < 1e-4 and the verdict. `--check` exits non-zero unless every ratio up to 200 reaches 1e-4 within `--max-iterations`. |

## How to run
```bash
python3 Tests/Canonical_RANS/Aspect_Ratio_Hill_2D/aspect_ratio_sweep.py --gust build-rel/gust --out ar-runs --np 2 \
    --ratios 10 25 50 100 200 500 1000 2000 5000 --max-iterations 8000
```
The roughness runs below copy a member's `inputs` and set `bc.surface_z0` and `bc.inflow_z0` to the same value.

## Expected Results
Iterations to r < 1e-4 with z0 = 0.01 m (2 ranks, 2026-10-25):

| dx/dz1 | 10 | 25 | 50 | 100 | 200 | 500 | 1000 | 2000 | 5000 |
|---|---|---|---|---|---|---|---|---|---|
| dz1 (m) | 10 | 4 | 2 | 1 | 0.5 | 0.2 | 0.1 | 0.05 | 0.02 |
| Iterations | 234 | 183 | 189 | 181 | 284 | 855 | 2509 | 5480 | not converged |

The largest cell aspect ratio in the mesh is about 6 % above the nominal ratio (terrain-following compression).

With `steady.cfl_max = 4` instead of 16, dx/dz1 = 10, 100 and 200 converge in 82, 116 and 127 iterations: on this
case the default cap is not the fastest choice.

Thickness of the first layer relative to z0 (the same deck with other z0):

| dx/dz1 | 100 | 100 | 1000 | 1000 | 2000 | 2000 | 2000 | 5000 | 5000 |
|---|---|---|---|---|---|---|---|---|---|
| dz1/z0 | 2 | 100 | 2 | 10 | 3 | 4 | 5 | 2 | 20 |
| Iterations | 983 | 181 | 6540 | 2509 | > 8000 (r 1.2e-4, decreasing) | 6580 | 5480 | not converged (fixed state) | about 10 500 (r 3.0e-4 at 8000) |

The case dx/dz1 = 5000 with dz1 = 2 z0 also stops at a fixed state (r 0.377) on a domain twice as long and wide.
Its cause is not known yet (the residual sits in k and u near the ground over the hill; no reversal, no spanwise
drift). The paper (`Tests/Paper/JWEIA2026`, Appendix D) recommends dx/dz1 up to about 500 for routine use and dz1 of at
least 5 z0 at higher ratios.

With the Newton-Krylov stage (`steady.stage = simplec_newton`, switch at r = 1e-2, limit only, `turb.nut_relax = 1`)
at dx/dz1 = 500, every step reached 2000 GMRES(30) iterations at a Krylov relative residual of 0.95-0.995; the stage
does not support `convection.hybrid_sigma_max`.

## CTest
`rt_aspect_ratio_hill` (label `regression`, 2 ranks, about 3 min): dx/dz1 = 10, 100 and 200 each reach r < 1e-4
within 400 iterations. Mutant: without `simplec.w_vnn_cap` (the w-row cap for thin cells), dx/dz1 = 10 still converges
(266 iterations) but 100 and 200 do not reach 1e-4 within 400, and the check fails.
