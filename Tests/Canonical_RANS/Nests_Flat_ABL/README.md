# Nests_Flat_ABL

## Purpose
Nests on a flat, horizontally uniform atmospheric boundary layer, neutral, buoyant (unstable) and stable, one way
(`nest.feedback = none`) and joint (`nest.feedback = defect`). Beyond the plan.

The outer domain is doubly periodic over flat ground and starts from its own steady column
(`init.first_guess = column`). Its steady 3-D state is then the column, the same at every (x, y). A nest at 2 2 1
inside it has the outer cells' heights, so it must return the outer domain's level means, one way and after the
joint run. The check compares the nest with the outer domain, and the outer domain with the column. Published
numbers for these set-ups are context only: the column is not compared with LES.

## Files
| File | Role |
|---|---|
| `inputs_neutral_s` | The neutral variant of case S of Moeng and Sullivan (1994), as rerun by Churchfield et al. (2010) (MS94's case S has a small surface heat flux, about 0.05 K m/s): geostrophic wind 15 m/s, f = 1.028e-4 1/s (latitude 45°), z0 = 0.16 m, theta 300 K to 468 m, +8 K over the next 83 m, then 0.003 K/m; no surface heat flux. 3000 m x 3000 m x 1000 m, 16 x 16 x 64 cells. Nest from 750 m to 2250 m in x and y. theta frozen at the inflow table's, which is the sounding: a neutral column keeps the sounding, so this is also the column's. The Rayleigh layer from 450 m damps u, v and w towards `column_neutral_s.dat`, in the nest too (`nest.damping = outer`); `theta` in `bc.rayleigh_fields` acts in the 1-D column only. |
| `inputs_buoyant_b` | Case B: geostrophic wind 10 m/s, surface heat flux 0.24 K m/s, theta 300 K to 937 m, +8 K over the next 167 m, then 0.003 K/m. 5000 m x 5000 m x 2000 m, 16 x 16 x 64 cells; nest from 1250 m to 3750 m. theta frozen at the stratified column's; layer from 700 m on u, v and w towards `column_buoyant_b.dat` (theta in the 1-D column only), `nest.damping = outer`. |
| `inputs_stable_gabls1` | The set-up of GABLS1 (Beare et al. 2006): geostrophic wind 8 m/s, f = 1.39e-4 1/s, z0 = 0.1 m, theta 265 K to 100 m, then 0.01 K/m. GABLS1 cools the surface by 0.25 K/h and has no steady state; here the surface carries a fixed flux of -0.0025 K m/s, which the layer on theta above 200 m supplies. 400 m cube, 16 x 16 x 64 cells; nest from 100 m to 300 m. theta solved; the layer damps w and theta only (no table). |
| `column_neutral_s.dat`, `column_buoyant_b.dat` | The decks' `bc.inflow_file` tables (z u v theta k): the steady column's u, v and k, with theta replaced by the deck's sounding (below). |
| `make_column_table.py` | Writes those tables from a deck (below). |
| `check_nests_flat_abl.py` | Runs a deck one way or one way then joint, and checks the plotfiles (Expected Results). |

All three decks run grid sequencing on the outer domain (`mesh.hierarchy_levels = 2`, `steady.sequencing = 1`) and
keep the closure defaults: they are a consistency check of the nests, not a run of the recommended setting
(`turb.realizability = 0.309`, `turb.nut_relax = 0.5` and the blend).

### Why theta is frozen in S and B
With theta solved, entrainment at the inversion warms the mixed layer without end: the 1-D column and the 3-D run
have no steady state. `simplec.theta = frozen` holds theta at the steady column's (`init.first_guess = column` with
a non-neutral `init.column_stability`); the closure keeps N^2.

### Why u and v are damped in S and B
Without a layer on u and v below the inversion, rolls grow there in the 3-D run. The layer damps u, v and w towards
the column's table (theta is frozen), so it removes departures from the column and leaves the column itself unchanged. The
nests must damp the same fields (`nest.damping = outer`): with the default (w and theta only) the forcing jumps at
each footprint's edge, the neutral joint run is not steady after 10000 exchanges, and the buoyant nest has rolls of
4.5e-3 m/s.

### The tables
The column takes its start theta and its theta target from the inflow table. A table with the column's own theta
would change the column's start (the buoyant column moved 0.26 K off the table). The tables therefore keep the
column's u, v and k and carry the deck's sounding theta:
```bash
python3 make_column_table.py inputs_buoyant_b --column-deck column_deck
mpiexec -n 2 <build>/gust column_deck            # writes column_guess.dat
python3 make_column_table.py inputs_buoyant_b --table column_guess.dat --out column_buoyant_b.dat
```
`--column-deck` writes the deck without its nest, its table and its layer on u and v, stopped after one iteration.
S: u* 0.5716 m/s, 132 column steps; B: u* 0.5339 m/s, 364 column steps.

### Why GABLS has no table
With u and v also damped towards a column table, the column's surface flux, as a function of the surface theta, has
a maximum of about -2.63e-3 K m/s near 265.75 K. The table's own solution (265.6 K, -2.50e-3 K m/s) lies past that
maximum, where the column's search for the surface theta does not look, and the run aborts. Without the u/v layer
the flux is monotone over 265.2 to 267 K. The stable deck therefore damps w and theta only and solves theta. Its
column start is close to the steady state, so the stop takes no drop of the residuals
(`steady.residual_drop_orders = 0`); every residual must still fall below `steady.residual_tol`.

## Data and acknowledgements
No data are downloaded or redistributed. The decks use published set-ups:
- Moeng, C.-H. and Sullivan, P. P. (1994). A comparison of shear- and buoyancy-driven planetary boundary layer
  flows. *J. Atmos. Sci.* 51, 999–1022. https://doi.org/10.1175/1520-0469(1994)051<0999:ACOSAB>2.0.CO;2 (cases S and B).
- Churchfield, M. J., Vijayakumar, G., Brasseur, J. G. and Moriarty, P. J. (2010). Wind energy-related atmospheric
  boundary layer large-eddy simulation using OpenFOAM. NREL/CP-500-48905. https://www.nrel.gov/docs/fy10osti/48905.pdf
  (the reruns of S and B these decks follow).
- Beare, R. J. et al. (2006). An intercomparison of large-eddy simulations of the stable boundary layer.
  *Boundary-Layer Meteorol.* 118, 247–272. https://doi.org/10.1007/s10546-004-2820-6 (GABLS1).

We thank the authors for the published set-ups.

## How to run
```bash
ctest --test-dir <build> -R rt_nests_flat_abl
ctest --test-dir <build> -L nightly -R rt_nests_flat_abl_stable_joint
```
By hand, in a scratch directory:
```bash
python3 check_nests_flat_abl.py -- mpiexec -n 2 <build>/gust inputs_buoyant_b
python3 check_nests_flat_abl.py --one-way -- mpiexec -n 2 <build>/gust inputs_stable_gabls1
```
The script copies the deck's table into the current directory and appends `nest.feedback = defect` (or `none` with
`--one-way`).

## Expected Results
`check_nests_flat_abl.py` requires:
1. the run exits 0; the log has "steady after" for the outer domain and the nest, and (joint) "joint run: steady";
2. one plotfile each of the outer domain and the nest (one way), or two (one way, then joint);
3. each domain horizontally uniform: u, v and k within 1e-6 of the level mean (S, B); for GABLS within 1e-4 on the
   levels whose mean k is above 1e-4 of the largest (above them the departure is printed, not bounded: see below);
4. the outer level means within 0.3 m/s (u, v), 10 % of the largest k and 0.02 K of the column's
   (`column_guess.dat`);
5. the nest's level means within 1e-6 m/s (S, B) or 2e-4 m/s (GABLS), 1e-4 of the largest k and 1e-5 K of the
   outer domain's.

Measured on 2 ranks with the CTest flags (Apple M-series laptop):

| Case | Outer iterations | Joint exchanges | Nest vs outer, one way (m/s) | Nest vs outer, joint (m/s) | Outer vs column, u/v (m/s) | Outer vs column, k | Largest departure from uniform (bounded levels) | Time |
|---|---|---|---|---|---|---|---|---|
| S | 701 | 514 | 1.1e-7 | **8.5e-10** | 0.149 | 2.7 % | 2.0e-7 | 35 s |
| B | 175 | 94 | 2e-8 | **6e-10** | 0.005 | 1.4 % | 3.7e-8 | 9 s |
| GABLS, one way | 567 | | 3.2e-5 | | 0.061 | 3.5 % | 2.4e-5 (4.3e-5 above) | 11 s |
| GABLS, joint | 567 | 8516 | 3.2e-5 | **6.5e-5** | 0.091 | 3.7 % | 2.1e-5 (4.0e-4 above) | 7 min |

theta: frozen in S and B (the outer domain is the column's to 2e-12 K); in GABLS the outer domain is 2.0e-3 K (one
way) and 0.012 K (joint) from the column, and the nest within 1.6e-6 K of the outer domain.

### The joint stable run
The outer domain's stop (`steady.residual_tol = 5e-5`, no drop) comes early for this deck. Run alone to 20000
iterations (the residuals reach about 8e-8), its mean moves 0.066 m/s in u near 316 m, 0.032 m/s in v and 0.013 K
in theta near the ground. The joint run iterates on and moves the outer mean by the same amounts: that is the
change of the outer domain against the column between the one-way and joint rows (0.061 to 0.091 m/s, 2.0e-3 to
0.012 K). The slow mode is an inertial oscillation of the outer mean in pseudo-time (2π/f is about 7200
iterations); it is why the joint run takes 8516 exchanges. A tighter stop would make the outer domain's own run
take most of those iterations, so the deck keeps the default and the joint test is labelled nightly.

u and v are not damped aloft in GABLS. Above the turbulent layer (mean k below 1e-4 of the largest) the joint run
settles with a roll in u at 320 to 340 m, uniform in x and varying in y, of 4.0e-4 m/s in both domains; the outer
domain run alone stays uniform to 2e-6 there. With k about 0, an x-uniform u(y) in geostrophic balance is a
neutral mode: nothing damps it, and its size depends on the iteration history. The check therefore prints the
departure above the turbulent layer and bounds it only below (2.1e-5 after the joint run).

### The column and the 3-D run
The 3-D outer domain is not a fixed point of the 1-D column for these capped, stratified columns: the first 3-D
residuals are 0.06 to 0.3 of the first-iteration scale, and the steady outer domain is 0.005 to 0.15 m/s from the
column. The bound on the outer domain against the column is loose for that reason. It catches a 3-D state far from
the column, such as theta frozen at the sounding in place of the column (5.1 K and 1.46 m/s off in B).

## CTest
| Test | Deck | Mode | Label |
|---|---|---|---|
| `rt_nests_flat_abl_neutral` | `inputs_neutral_s` | one way, then joint | regression |
| `rt_nests_flat_abl_buoyant` | `inputs_buoyant_b` | one way, then joint | regression |
| `rt_nests_flat_abl_stable` | `inputs_stable_gabls1` | one way | regression |
| `rt_nests_flat_abl_stable_joint` | `inputs_stable_gabls1` | one way, then joint | nightly |
