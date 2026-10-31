# Forested_Hill

## Purpose
A uniform canopy on a periodic train of cosine hills in a neutral flow driven by a constant pressure gradient. Two
geometries, classified on the regime map of Poggi et al. (2008), fig. 1:

- **Deep canopy, long hill** (regime I, the range of Finnigan and Belcher 2004): H = 4 m, L = 40 m, a 20 m canopy, L_c = 10 m
  (H_c/L_c 2, L/L_c 4, envelope 3.2). The lower canopy should follow the balance of eq. 3 of Poggi et al.,
  −dp/dx = (u|u| − U_b²)/L_c, with L_c = 1/(C_d a). Here the background is the constant forcing, U_b² = F L_c, so
  u|u| = L_c (F − dp/dx).
- **The flume's geometry** of Poggi et al. (2008) scaled by 100: H = 8 m, L = 80 m, a 10 m canopy, at H_c/L_c = 0.13 (their
  regime V), 0.26 and 0.52. The paper measures reversed flow low in the canopy on the lee slope at 0.13. At 0.52 the case is long,
  deep and at fixed pressure, regime I: the predicted lee reversal, not the flume's regime.

## Files
| File | Role |
|---|---|
| `inputs_deep_hill` | The deep canopy: 40 x 4 x 39 cells, a 1 m first layer stretched by 1.07 (top at about 186 m), LAI 10. |
| `inputs_flume_hill` | The flume: 80 x 4 x 53 cells, four boxes, a 0.25 m first layer stretched by 1.05 (top at 61.4 m). |
| `check_forested_hill.py` | Writes the tables, runs the six runs, and checks them; the exit code is the verdict. |

C_d 0.2 throughout; the LAI gives H_c/L_c. The check passes what follows the executable to every run after the deck.
F only sets the speed scale (u* = 0.5 m/s): no published value sets it. The paper states no flow depth: the 61 m lid is a
choice, from its Q = 120 l/s, 0.9 m width and Fr_b < 0.1 (depth above about 0.57 m in the flume), and is a rigid lid.

## Data and acknowledgements
No data files: the geometry and the relations come from the papers, whose measurements appear only as figures, so no
measured profile is compared here. D. Poggi, G. G. Katul, J. J. Finnigan and S. E. Belcher (2008), Q. J. R. Meteorol. Soc.
134, 1095–1112, https://doi.org/10.1002/qj.276; J. J. Finnigan and S. E. Belcher (2004), Q. J. R. Meteorol. Soc. 130,
1–29. We thank the authors for reporting the flume's geometry and regime map.

## How to run
```bash
ctest --test-dir <build> -R rt_forested_hill --output-on-failure
```

## Expected Results
Measured on macOS arm64 (AppleClang, double precision), about 2 minutes:

| Check | Measured | Tolerance |
|---|---|---|
| Deep hill: rms(u − sign(R)√\|R\|) / rms(u) at h/4, R = L_c (F − dp/dx) | 0.195 | 0.25 |
| H_c/L_c = 0.52: reversed columns at h/4 | 25 of 80, x/L from 2.92 to 4.12 (crest at 2, valley at 4); u down to −0.34 m/s | at least 10 %, within x/L 2 to 4.25 |
| H_c/L_c = 0.52: 1 and 4 ranks against 2, largest u difference / max u | 3.2e-7, 7.0e-7 | 1e-5 |
| Printed only: H_c/L_c = 0.13, 0.26, u(h/4) | 0.66 to 1.37; 0.17 to 0.91 m/s, **no reversal** | none |
| Printed only: deep hill, u(h/4) | −0.08 to 0.21 m/s, 13 of 40 columns reversed, x/L 3.75 to 4.95 | none |

Iterations to the verdict: deep 8744, H_c/L_c = 0.13 1523, 0.52 1417.

**Limits.**
- **The flume's measured reversal is not reproduced.** At its H_c/L_c = 0.13 this case does not reverse, nor at twice the drag; at
  four times it does, in another regime. The threshold between 2 and 4 and the cause are not tested. Candidates are the
  flume's top-heavy rods, the lid height and the closure.
- The balance leaves 20 % of the speed's rms on the deep hill: it neglects advection and stress. The tolerance 0.25 was set after
  the first runs. It catches C_d 0.3 and 0.15 (0.275, 0.271) and does not catch C_d 0.25 (0.203).
- On the deep hill the reversal lies from just before the valley up the next windward slope, downstream of the lee slope
  where linear theory with pressure in phase with the hill places it. The run's own pressure is used in the check, so
  the check tests the local balance, not the pressure's predicted phase. This phase shift is not measured here.
- Before the verdict the answers on different rank counts differ (8e-4 after 200 iterations on an earlier deck): the pressure
  solve stops at its tolerance. They agree at the verdict.

**What the check detects** (each makes it fail): `canopy.cd` 0.3 and 0.15 (the balance), `canopy.tke = off` (0.715, and no
reversal).

## CTest
`rt_forested_hill` (label `regression`, MPI builds).
