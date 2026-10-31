# Homogeneous_Canopy

## Purpose
A horizontally uniform forest in a neutral flow driven by a constant pressure gradient F over flat ground: the
simplest flow with a published canopy relation, and no fitted inputs. A 20 m canopy with C_d 0.2 at LAI 2, 5 and 10.
The check tests the drag and the canopy k terms against:

- the column's momentum budget, τ_s + ∫ f U² dz = F H;
- the balance deep in a dense canopy, where the drag carries the pressure gradient: U = √(F/f), a closed form;
- the exponential upper canopy of Inoue (1963), U = U_h exp(−n (1 − z/h)), with n = f h / (2β²), β = u*(h)/U_h and
  u*(h)² = F (H − h);
- β in [0.25, 0.40], around the values of about 0.3 in Poggi et al. (2008);
- the same answer on 1, 2 and 4 ranks at LAI 10.

## Files
| File | Role |
|---|---|
| `inputs_homogeneous_canopy` | The deck: 8 x 8 x 143 cells, a 0.25 m first layer stretched by 1.02 (top at 199.7 m), `abl.drive = pressure_driven`, the column first guess, `steady.dtau_floor = 1`. |
| `check_homogeneous_canopy.py` | Writes the canopy tables, runs the five runs, and checks them; the exit code is the verdict. |
| `compare_erf.py` | ERF on the same column against GUST's runs (not a CTest; it needs an ERF build with erf-model/ERF#4340). |

## Data and acknowledgements
No external data. β is compared with D. Poggi et al. (2008), Q. J. R. Meteorol. Soc. 134, 1095–1112, https://doi.org/10.1002/qj.276. The relations are from E. Inoue (1963), J. Meteor. Soc. Japan 41(6), 317–326, and I. N. Harman and
J. J. Finnigan (2007), Boundary-Layer Meteorol. 123, 339–363, https://doi.org/10.1007/s10546-006-9145-6.

## How to run
```bash
ctest --test-dir <build> -R rt_homogeneous_canopy --output-on-failure
```

## Expected Results
Measured on macOS arm64 (AppleClang, double precision), about 45 s:

| Check | LAI 2 | LAI 5 | LAI 10 | Tolerance |
|---|---|---|---|---|
| Budget, ratio − 1 | +3.7e-5 | +1.5e-5 | −8.0e-5 | 0.01 |
| U(h/4) / √(F/f) − 1 | not reached | +1.2e-3 | −5.7e-5 | 0.02 (LAI 5, 10) |
| n fitted / f h/(2β²), on U/U_h ≥ 0.35 | 2.40 / 2.58 | 6.22 / 5.70 | 12.42 / 9.96 | 15 %, fit rms 0.05 (LAI 2, 5) |
| β | 0.278 | 0.296 | 0.317 | [0.25, 0.40] |
| Layers in the fit | 16 | 6 | 3 | |

Iterations from the column first guess: 144, 210 and 595. On 1 and 4 ranks at LAI 10 the wind differs by 1.8e-15 m/s
from 2 ranks (tolerance 1e-8). From a profile first guess the same mesh does not reach the verdict in 10000 iterations, with the same
wind: the whole column adjusts slowly.

**Limits.** At LAI 10 the fitted rate is 25 % above the relation and is printed, not checked: the e-folding depth is
about 2.5 layers. On a finer mesh (0.05 m first layer, stretched by 1.0135, lid at 117 m, not a CTest) the ratio n fitted / relation is 0.82, 0.97 and 1.07
at LAI 2, 5 and 10, so the LAI 10 excess falls with the mesh, and the LAI 2 ratio moves below the tolerance. This exclusion and the 15 % tolerance were set after the first runs. The theory takes β as a
constant; here it rises with LAI by 14 % from LAI 2 to 10. Sources differ on neutral β (about 0.3 to 0.39), and Harman and
Finnigan's own value has not been read from the paper.

**What the check detects** (each makes it fail): `canopy.cd` 0.4 in the deck, `canopy.tke = off`,
`canopy.length_alpha = 0`. It does not detect `canopy.beta_d = 2` or a ground z0 of 0.1 m: these move the numbers by less than
their tolerances.

**ERF on the same column** (`compare_erf.py`; the benchmark ERF build, 2 x 2 x 200 cells of 1 m,
12000 s from rest). GUST minus ERF:

| LAI | u(h/4) | u(h/2) | u(h) | k(h) | u(1.5 h) | u(5 h) | ERF drift in its last 6000 s at 5 h |
|---|---|---|---|---|---|---|---|
| 2 | −0.2 % | −0.3 % | −0.1 % | +0.7 % | −0.7 % | −0.4 % | 0.035 m/s |
| 5 | −0.0 % | −0.9 % | +3.2 % | +6.1 % | +2.2 % | +1.0 % | 0.029 m/s |

The closures share their canopy terms (extern_erf), so this tests the two solvers and meshes, not the closure. ERF's
upper column was still accelerating when it stopped; its layers are 1 m at the canopy top where GUST's are 0.65 m, and
the LAI 5 difference at the top has not been separated into these two.

## CTest
`rt_homogeneous_canopy` (label `regression`, MPI builds).
