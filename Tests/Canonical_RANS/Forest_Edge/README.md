# Forest_Edge

## Purpose
A forest-clearing-forest pattern, the geometry of Dupont and Brunet (2008), sec. 2, in a neutral flow driven by a constant
pressure gradient: 690 m periodic in x, an 18 m forest of uniform leaf density (306 m on the 6 m columns), and a 21.2 h
clearing. The LES has a geostrophic drive and a Rayleigh layer where this has a pressure gradient and a slip lid. Behind
the leading edge the flow in the canopy slows and lifts out through the canopy top. Dupont and Brunet end the adjustment
where the mean upward motion at the canopy top falls below w/u_ref = 0.01, with u_ref = u(x = −8.5 h, z = 2 h). Their LES
gives about 9 h at LAI 2; its mean profiles were validated against the wind-tunnel data of Raupach et al. (1987). Denser canopies adjust faster (Belcher
et al. 2003; Dupont and Brunet 2008).

## Files
| File | Role |
|---|---|
| `inputs_forest_edge` | The deck: 115 x 4 x 49 cells (6 m columns, four uneven boxes), a 1 m first layer stretched by 1.05 (top at 198.4 m), `steady.dtau_floor = 1`. |
| `check_forest_edge.py` | Writes the canopy tables on the columns' centres, runs LAI 2 and 5 and LAI 2 on 1 and 4 ranks to the verdict, and checks them; the exit code is the verdict. |
| `mesh_study.py` | The 3 m mesh, with its deck and table (not a CTest). |

## Data and acknowledgements
No data files; the measurements appear in the papers as figures and are not compared here. The geometry and the
adjustment length are from S. Dupont and Y. Brunet (2008), Boundary-Layer Meteorol. 126, 51–71,
https://doi.org/10.1007/s10546-007-9216-3; the scales from S. E. Belcher, N. Jerram and J. C. R. Hunt (2003), J. Fluid Mech.
488, 369–398, https://doi.org/10.1017/S0022112003005019. We thank the authors, and M. R. Raupach and colleagues, whose
1987 wind-tunnel data the LES was validated against.

## How to run
```bash
ctest --test-dir <build> -R rt_forest_edge --output-on-failure
```

## Expected Results
Measured on macOS arm64 (AppleClang, double precision), about 12 minutes on three test slots:

| Check | LAI 2 | LAI 5 | Tolerance |
|---|---|---|---|
| Budget / (F H L) − 1 | −1.9e-3 | −4.9e-3 | 0.01 |
| Volume flux spread between columns | 7.0e-9 | 8.1e-9 | 1e-6 |
| Adjustment length | 8.42 h | 7.14 h | LAI 2 in [7.5, 10.5] h; LAI 5 shorter |
| Printed: peak w/u_ref | 0.070 at 0.5 h | 0.116 at 0.2 h | none |
| Printed: x0 = L_c ln(U_h h/(u* L_c)), Belcher et al. eq. 3.3 | 1.23 h | 1.22 h | none |

x0 is the origin of the internal boundary layer, a different length and an order of magnitude: it is printed for reference.
Dupont and Brunet's denser case 2 gives 8 h, with another leaf profile; 7.14 h here is not checked against it.
On 1 and 4 ranks at LAI 2 at the verdict the velocities differ by 1.4e-7 and 1.5e-8 of the largest u, within 1e-6.
A 3 m mesh (`mesh_study.py`, 17898 iterations to the verdict) also gives 8.42 h.

**Limits.**
- Dupont and Brunet set C_d LAI equal to that of Raupach et al. (1987) without giving the value here. **C_d 0.2 is an
  assumption.**
- The range [7.5, 10.5] h was set after the first runs and mutants: their 9 h, and 8 to 10 h for their cases 1, 2 and 4 (other
  leaf profiles; their case 3, with a sparse trunk space, had not adjusted at the end of the forest, 20 h), with half a canopy height either side.
- **The low-k zone behind the edge is not reproduced.** In the LES, k at the canopy top falls behind the edge and recovers
  further on. Here it rises at once to its forest value (k/u_ref² 0.03 in the clearing, 0.06 at 0.8 h, 0.068 from 3 h).
- The enhanced gust zone and the skewness of the LES have no counterpart in a steady RANS model.

**What the check detects** (each makes it fail): `canopy.cd` 0.4 (budget −48 %; its adjustment 7.44 h was on the first tables, and is
close to the range's edge, so the budget is what catches it), `canopy.cd` 0.02 (budget, no adjustment found), `canopy.tke = off`
(adjustment 7.12 h; budget and flux pass). The tables' first form, off the columns' centres, gave a smeared edge with a flux spread of
8.6e-5.

## CTest
`rt_forest_edge` (label `regression`, MPI builds).
