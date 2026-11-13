# Cylinder_Laminar

## Purpose
RT-6.4.1 to RT-6.4.3 (plan rev. 6, sub-phase 6.4): laminar flow past a circular cylinder, run as steady RANS and as
unsteady RANS (dual time stepping, `steady.stage = unsteady`), with no turbulence model (`turb.closure = none`). The
cylinder is a staircase of solid cells (`buildings.wall_form = drag`: no-slip walls on the blocked faces), and its
drag and lift come from `buildings.force_file` (docs: Buildings, Force on the buildings).

- **Re = 40:** the flow is steady. Both approaches give the same answer: two attached eddies behind the cylinder,
  compared with the drag and wake length of Dennis & Chang (1970) and Fornberg (1980).
- **Re = 100:** the true flow sheds vortices (the Kármán street). The steady solver converges to the symmetric
  steady solution (Dennis & Chang 1970, Fornberg 1980: lower drag, a long wake). The unsteady run sheds; its time
  mean is compared with Park, Kwon & Choi (1998), and its shedding frequency with Williamson (1996).

The case shows what the steady solver gives for a flow that is not steady, and that the unsteady driver recovers the
shedding.

## Files
| File | Role |
|---|---|
| `inputs_cylinder` | The deck: x from −6 to 14 m, z from 0 to 12 m, y periodic and two cells thick (0.1 m); a uniform inflow at 1 m/s, an outflow face, slip ground and top. The command line gives `mesh.n_cell`, `buildings.file`, `bc.inflow_file`, `turb.nu_molecular = 1/Re`, the stage and the unsteady keys. |
| `make_cylinder.py` | Writes the staircase map `cylinder_<n>.txt` (`buildings.format = xyz`, `x y height zstart`: a cell is solid when its centre is inside the circle) for n cells per diameter, and the inflow file. |
| `cylinder_10.txt`, `cylinder_20.txt` | The maps for D/d = 10 and 20 (cells of 0.1 and 0.05 m). |
| `inflow_uniform.txt` | The inflow: u = 1 m/s, θ = 300 K, k = 0 at every height. |
| `probes_cylinder.txt` | Probes 2 D and 5 D behind the axis, and 2 D behind, 0.5 D above it (`unsteady.probes_file`). |
| `inputs_abort_laminar` | The deck for the CTest `abort_laminar_cross_checks`: the surface layer and the buildings' wall law and wall k, which a laminar run refuses. |
| `inputs_abort_laminar_other` | The deck for the CTest `abort_laminar_other_options`: a canopy, the column first guess and the Newton correction, also refused. |
| `compare_cylinder.py` | Cd, Cl (or the mean Cd, the lift amplitude and St) from the force file, Lw/D from the last plotfile; prints the published values; `--check-*` bands for CTest. |

## Geometry and numbers
The cylinder has D = 1 m and lies along y, so the flow plane is x–z. Its axis is at x = 0, z = 6 − d (d the cell
size): on a cell face, so the staircase is symmetric about the axis and spans the full diameter, one cell closer to
the ground than to the top. That asymmetry, d in 12 m, is the only seed of the shedding. (An axis on a cell centre
gives an odd number of cells across and a frontal height of D − d: at D/d = 10 it lowered Cd from 1.70 to 1.59 and
Lw/D from 1.96 to 1.54.)

- U = 1 m/s, so Re = 1/ν: `turb.nu_molecular = 0.025` for Re 40, `0.01` for Re 100.
- Cd = 2 Fx / (U² D Ly) and Cl = 2 Fz / (U² D Ly), with Ly = 0.1 m the span; St = f D / U.
- Lw is measured along the axis row from the rear of the circle (x = D/2) to where u (u_mean when unsteady) turns
  positive.
- Blockage D/H = 1/12 between slip walls; the published values are for unbounded flow. The domain reaches 6 D
  upstream and 14 D downstream.

## Data and acknowledgements
No data are downloaded; the case compares with published numbers only, quoted in `compare_cylinder.py`:

- S. C. R. Dennis and G.-Z. Chang (1970), Numerical solutions for steady flow past a circular cylinder at Reynolds
  numbers up to 100, J. Fluid Mech. 42(3), 471–489, https://doi.org/10.1017/S0022112070001428 (table 2: Cd and the
  wake length, given there in radii).
- B. Fornberg (1980), A numerical study of steady viscous flow past a circular cylinder, J. Fluid Mech. 98(4),
  819–855, https://doi.org/10.1017/S0022112080000419 (Cd at Re 20, 40 and 100).
- J. Park, K. Kwon and H. Choi (1998), Numerical solutions of flow past a circular cylinder at Reynolds numbers up to
  160, KSME Int. J. 12(6), 1200–1205, https://doi.org/10.1007/BF02942594 (table 1: mean Cd and the lift amplitude).
- C. H. K. Williamson (1996), Vortex dynamics in the cylinder wake, Annu. Rev. Fluid Mech. 28, 477–539,
  https://doi.org/10.1146/annurev.fl.28.010196.002401 (St = −3.3265/Re + 0.1816 + 1.6×10⁻⁴ Re for laminar parallel
  shedding).
- Y.-C. Li, S. Peng and T. Kouser (2022), Effect of wall slip on laminar flow past a circular cylinder,
  arXiv:2205.15886, https://doi.org/10.48550/arXiv.2205.15886 (table 1: Fornberg's and Park et al.'s wake lengths,
  read from their figures).

We thank these authors.

## How to run
```bash
python3 make_cylinder.py --cells-per-d 20 --out cylinder_20.txt     # maps for other resolutions
# Re = 40, steady, D/d = 10
mpiexec -n 2 gust inputs_cylinder "mesh.n_cell=200 2 120" buildings.file=cylinder_10.txt \
    bc.inflow_file=inflow_uniform.txt turb.nu_molecular=0.025 steady.stage=simplec steady.max_iterations=3000
python3 compare_cylinder.py --run . --re 40 --cells-per-d 10
# Re = 40, unsteady (the same answer): 120 s in steps of 0.05 s, averages from 80 s
mpiexec -n 2 gust inputs_cylinder "mesh.n_cell=200 2 120" buildings.file=cylinder_10.txt \
    bc.inflow_file=inflow_uniform.txt turb.nu_molecular=0.025 steady.stage=unsteady unsteady.stop_time=120 \
    unsteady.dt=0.05 unsteady.inner_drop_orders=2 unsteady.max_inner=40 unsteady.average_start=80
python3 compare_cylinder.py --run . --re 40 --cells-per-d 10 --t-avg 80
# Re = 100, unsteady, D/d = 10: 160 s in steps of 0.05 s, averages from 100 s
mpiexec -n 2 gust inputs_cylinder "mesh.n_cell=200 2 120" buildings.file=cylinder_10.txt \
    bc.inflow_file=inflow_uniform.txt turb.nu_molecular=0.01 steady.stage=unsteady unsteady.stop_time=160 \
    unsteady.dt=0.05 unsteady.inner_drop_orders=2 unsteady.max_inner=40 unsteady.average_start=100 \
    unsteady.probes_file=probes_cylinder.txt
python3 compare_cylinder.py --run . --re 100 --cells-per-d 10 --t-avg 100
```
D/d = 20 uses `"mesh.n_cell=400 2 240"` and `cylinder_20.txt`. A Re 100 steady run uses `steady.stage=simplec` with
`turb.nu_molecular=0.01`; it stops at `steady.max_iterations`.

## Expected Results
The case is a consistency check of the two drivers, not an accuracy study: the published values are quoted for
context. All runs below are D/d = 10, H = 12, 2 ranks.

**Re = 40: steady and unsteady agree.**

| Run | Cd | Cl | Lw/D | Iterations |
|---|---|---|---|---|
| **steady (`simplec`)** | **1.6985** | **0.001** | **1.963** | **733** |
| unsteady, dt 0.05 s to 120 s, mean over 80–120 s | 1.6985 | lift amplitude 0.000 | 1.963 (u_mean) | 4430 (1.8 per step) |
| Published (unbounded) | 1.50–1.52 | 0 | 2.24–2.35 | |

The unsteady run settles to the steady answer: the same Cd to four digits and the same wake, with no lift
oscillation.

**Re = 100: the unsteady run sheds; the steady run does not settle.**

| Run | Cd | Lift | St | Lw/D | Iterations |
|---|---|---|---|---|---|
| steady (`simplec`) | 1.43 at the stop | Cl 0.13 at the stop | – | – | 6000 (stop; residual ~1e-2) |
| **unsteady, dt 0.05 s** | **1.428 (mean)** | **amplitude 0.209** | **0.164** | **1.86 (mean)** | 14449 (4.5 per step) |
| unsteady, dt 0.025 s | 1.445 (mean) | amplitude 0.240 | 0.167 | 1.72 (mean) | 23688 (3.7 per step) |
| Published | 1.33 (mean; steady branch 1.06) | amplitude 0.33 | 0.164 | 1.41 (mean; steady branch 6.6) | |

- Unsteady runs: 160 s, averages over 100–160 s (9 shedding periods). The cylinder's axis one cell below
  mid-height is the only seed. The lift stays within ±0.003 for the first 30 s, grows by a factor of 2–3 every
  10 s, and reaches a constant amplitude by about 100 s (dt 0.05; about 90 s with dt 0.025).
- Halving dt changes St by 2 % and raises the lift amplitude by 15 %: BDF1 damps the oscillation at this step.
- The steady solver does not reach the symmetric steady branch of Dennis & Chang (1970) and Fornberg (1980):
  its residual stalls near 1e-2 with an asymmetric wake, the signature of a flow with no stable steady state.

**Sensitivity of the Re 40 answer (context for the offset from the published values).**

| Change from the bold row | Cd | Lw/D |
|---|---|---|
| D/d = 20 (1313 iterations) | 1.693 | 2.05 |
| H = 24 (blockage 1/24; 847 iterations) | 1.646 | 2.05 |
| `convection.limiter = van_albada` | 1.710 | 1.95 |
| `convection.limiter = upwind` (first order everywhere) | 1.881 | 1.75 |

Cd stays about 10 % above the published values. Resolution moves it by 0.3 % and halving the blockage by 3 %. The
limiter rows show the size of the upwinding error: first order everywhere raises Cd by 11 % and shortens the wake
by 11 %, while the two smooth limiters differ by 0.7 %. The default MUSCL falls back to upwind only where van Leer
clips (the stagnation point, the separation points, the wake's minimum of u). Three errors of the discretization
contribute, and the first two change in opposite directions as the cells shrink, so the small change with
resolution does not show that each is small:
- upwinding at the clipped points: numerical viscosity lowers the effective Re, which raises Cd and shortens the
  wake;
- the staircase holds the tangential velocity at 0 half a cell inside the wall (docs: Buildings, Force on the
  buildings), a slip that lowers the friction drag;
- the domain reaches only 6 D upstream.

The bold Re 40 steady row is the CTest `rt_cylinder_re40`. Its bands reject the `upwind` run (Cd and Lw/D out).

## CTest
- `rt_cylinder_re40` (label `regression`, after `rt_cylinder_re40_run`): Re 40, D/d = 10, steady on 2 ranks; Cd,
  Cl and Lw/D within the bands pinned to this deck.
- `abort_laminar_nu_zero`, `abort_laminar_kl_key`, `abort_laminar_cross_checks`, `abort_laminar_other_options`,
  `abort_force_file_without_stage`
  (labels `regression;smoke`): the refusals of `turb.closure = none` and `buildings.force_file`.
- The unit tests `Closure.LaminarClosureGivesTheMolecularDiffusivities`,
  `BuildingLoads.PressureGradientTelescopesThroughABlock`, `BuildingLoads.ShearOverASolidLayerGivesThePlateDrag`
  and `BuildingLoads.FaceVolumeIsTheMeanOfItsTwoCells` check the closure and the force against closed forms; the
  `force_file_cube*` CTests (Canonical_RANS/Building_Cube) check the force file in unsteady runs and restarts.

The unsteady runs (15 min to 2.5 h on 2 ranks) and the sensitivity runs are not CTests.
