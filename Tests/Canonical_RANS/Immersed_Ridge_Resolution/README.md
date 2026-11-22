# Immersed_Ridge_Resolution

## Purpose
How fine the columns must be for an immersed ground (`terrain.representation = immersed_forcing`,
[docs](../../../docs/immersed_terrain.rst)) to give the crest speed-up of the same ground on a terrain-following
(TF) mesh. The ground is a 2-D Witch of Agnesi ridge, 200 m high with half-width 500 m (slopes up to 0.26, 15°),
steeper than the ridge of `../Immersed_Ridge_2D` (100 m, 7°). The case sweeps the column width dx at fixed 2 m
layers, and the layer thickness dz at dx = 64 m.

Neutral (uniform θ, no gravity), geostrophic wind 10 m/s with Coriolis, periodic in x (8192 m) with one column in
y, z0 = 0.1 m, slip top. The immersed mesh puts its bottom 2 m below the terrain's zero level
(`terrain.flat_height = 2`). Both decks use the recommended turbulence setting (docs/turbulence_model.rst),
`turb.shear_stencil = log_law` and `solver.poisson_backend = gmres_fft`. The immersed deck sets
`terrain.small_volfrac = 0.1`, as the measured cases of
`../../Real_Terrain/Immersed_Terrain_Cases` do. The stop rule is four orders of residual drop.

## Files
| File | Role |
|---|---|
| `inputs_ridge_immersed` | The immersed ground, dx = 64 m, 149 layers (dz = 2 m). |
| `inputs_ridge_body_fitted` | Its TF partner. |
| `make_levels.py` | Writes the layer file: dz up to 224 m, then stretched by about 1.1 per layer to the top (1002 m immersed, 1000 m TF). |
| `compare_dx.py` | Reads the runs of one or more dx and prints the crest speed-up, the lee minimum, k, the immersed minus TF gap per height with the observed order in dx, and the immersed / TF far-field and crest speeds; `--plot` draws dS along x. |

## Data
None: the case is analytic.

## How to run
One folder per run, each with its deck copied as `deck` and its layers as `zlevels.txt`. The decks are set for
dx = 64 m; for dx = 16, 32 and 64 m (512, 256 and 128 columns) at dz = 2 m:
```bash
case=<src>/Tests/Canonical_RANS/Immersed_Ridge_Resolution
for nx in 512 256 128; do
  dx=$((8192 / nx))
  for kind in body_fitted immersed; do
    top=1000; [ $kind = immersed ] && top=1002
    mkdir -p dx$dx/$kind && cd dx$dx/$kind
    python3 $case/make_levels.py 2 $top > zlevels.txt
    sed -e "s/^mesh.n_cell .*/mesh.n_cell = $nx 1 149/" \
        -e "s/^mesh.max_grid_size .*/mesh.max_grid_size = $((nx / 4)) 1/" $case/inputs_ridge_$kind > deck
    mpirun -np 4 <build>/gust deck > log
    cd ../..
  done
done
python3 $case/compare_dx.py 16=dx16/body_fitted,dx16/immersed 32=dx32/body_fitted,dx32/immersed \
  64=dx64/body_fitted,dx64/immersed --plot dS.png
```
Other layer thicknesses (dz must divide 224 m): `make_levels.py DZ TOP` prints the layer count, which goes in
`mesh.n_cell`. The immersed mesh's first layer must lie below the ground, so the immersed deck also needs
`terrain.flat_height` = dz + 2 and a top of 1000 m plus that (`mesh.prob_hi` and TOP); the dz = 8, 16 and 32 m runs
of Expected Results used flat heights of 10, 18 and 34 m. The dx = 64 m pair takes about 4000 iterations each;
dx = 16 m takes 9000 to 14000 (Expected Results).

## Expected Results
dS is the speed-up over the run's own far field at d m above the local ground (`compare_dx.py`). All runs below
reach the steady verdict.

**Column width, dz = 2 m** (iterations TF / immersed: dx 16: 14352 / 9255, dx 32: 7081 / 5083, dx 64: 4005 /
3653). At 2 m layers a wall-cell tolerance of 0.1 would give the same wall cells as the test without one (the
blocked part changes by more than 0.1 per layer, docs, *Wall cells*).

| dx (m) | Ground rise per column (m) | Crest dS10, TF / immersed | Gap at 10 m | Gap at 20 m | Gap at 40 m | Crest k at 10 m, TF / immersed (m²/s²) |
|---|---|---|---|---|---|---|
| 64 | 17 | 1.133 / 0.652 | −0.48 | −0.27 | −0.10 | 0.97 / 1.87 |
| 32 | 8 | 1.148 / 0.857 | −0.29 | −0.17 | −0.05 | 0.98 / 1.55 |
| **16** | **4** | **1.159 / 1.032** | **−0.12** | **−0.07** | **−0.02** | 0.99 / 1.22 |

The ground rise per column is the largest slope times dx. TF changes by 0.03 over the sweep. The immersed peak
stays in the crest column and rises toward TF: each halving of dx shrinks the gap by a factor of 1.6 to 2.4
(observed order 0.7 to 1.2 at 10 m, 1.0 to 1.4 at 40 m), about first order and not yet in the asymptotic range. At
80 m the gap is 0.014 or less. The lee minimum at 10 m is deeper with the immersed ground at every dx (−0.44 against
−0.36 at dx = 16 m) and nearer the crest; on the steep slopes at coarse dx the 10 m values come from cells in or
near the cut band and are less certain than the crest values.

Each run is divided by its own far field. The immersed far field is itself slower than TF's, by 3.0, 2.6 and 1.6 %
at dx = 64, 32 and 16 m, at every height (the immersed ridge takes more momentum out of the periodic domain). The
crest speed at 10 m, immersed over TF, is 0.751, 0.842 and 0.928; against TF's far field the 10 m gaps would be
−0.53, −0.34 and −0.16.

**Layer thickness, dx = 64 m.** At 40 m the immersed crest dS is 0.71 for dz = 2, 8, 16 and 32 m. At 10 m (crest
dS10, TF / immersed): dz = 2 m: 1.133 / 0.652; dz = 8 m: 1.151 / 0.525; dz = 16 m: 1.095 / 0.512; dz = 32 m:
0.935 / 0.580 (at dz = 32 m both runs' first cell centres are above 10 m, and these are the first centres' values).
The 10 m gap ranges from −0.36 to −0.63 with no monotone trend in dz, while the ratio dx/dz goes from 32 to 2: the
column width, not the cell aspect ratio, sets the gap.

**Why.** The geometry is exact for the piecewise-linear ground: an x face is blocked only below the ground at its
own node. The horizontal velocity lives on the column faces only, so where the ground rises by Δh across a column,
every layer of that column below the downstream node height has its downstream face shut, and the air must turn
up inside one column. In an inflow/outflow variant (dx = 64 m, dz = 1 m, Coriolis off) the pressure in the first
open cell above that band on the upper windward slope (x = −224 m) is 0.40 m²/s² above the value 60 m higher,
against 0.07 for TF, about 4 % of the near-ground dynamic pressure: a local pressure rise, which the crest sees as
an adverse pressure gradient. Those open x faces keep the gap (crest dS10 0.95 TF, 0.42 immersed), so it is not
the periodic domain returning the wake. The extra k follows the shear at the top of the slow layer; capping the
shear production removed most of the extra k but not the speed gap. The inflow/outflow and capped-production runs
are not in this folder.

**The guide for users.** Choose dx so the ground rises by about 5 m or less across one column on the slopes that
matter (dx ≲ 5 m / largest slope). At that point this ridge keeps a crest gap of about −0.12 at 10 m (about 10 %
of the speed-up). The layer thickness can then be a few metres. The lee converges more slowly than the crest. The
guide comes from one ridge: the rise per column and the columns per half-width change together in this sweep, so
it does not separate the two.

## CTest
None: the finest pair needs tens of minutes. The case is a record of the resolution study behind
`docs/immersed_terrain.rst` (section *Grid resolution*).
