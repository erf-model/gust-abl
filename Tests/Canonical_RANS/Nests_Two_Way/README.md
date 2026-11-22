# Nests_Two_Way

## Purpose
Two-way nesting (PR-70, `nest.feedback = defect`; PR-71, `nest.feedback = relax`; `nest.feedback = relax_consistent` and `nest.feedback = inject`, below). The outer domain and its nest first run one way, as in
`Nests_One_Process` except that the nest's levels already follow the joint rule (they bend to the outer domain's
surfaces aloft), so "one way" in the tables below is that first pass. Then they run jointly: every iteration the nest's faces are refreshed from the outer domain, and
the outer domain carries the defect source of the nest over the nest's footprint
(`docs/boundary_conditions.rst`, Two-way nests). The case checks that the joint run moves both domains towards the
uniform fine run:
1. the joint nest is closer to the uniform 25 m run than the one-way nest;
2. the outer domain outside the footprint, where the nest is not, is closer too: the nest's flow reaches the rest of
   the outer domain.

The hill is steeper than that of `Nests_One_Process` (100 m high, half-width 150 m, against 60 m and 250 m), so that
the 50 m outer cells miss a part of the flow the 25 m nest resolves.

## Files
| File | Role |
|---|---|
| `inputs_outer` | The outer domain: 32 x 24 x 16 cells of 50 m over 1600 m x 1200 m x 800 m (first layer 10 m), the round Witch of Agnesi hill in the middle, periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m, 20 degrees) on xlo, outflow on xhi, neutral, Rayleigh layer above 400 m on w and theta only (as every nest). No nest keys: the check script appends them. |
| `inputs_abort` | `inputs_outer` without `mesh.max_grid_size` and `steady.stage`, for the abort tests. |
| `check_nests_two_way.py` | Runs the uniform 25 m reference (the deck with 64 x 48 x 16 cells) and the joint run, and compares them (below). |
| `check_joint_restart.py` | Runs a joint run straight, the same run stopped after `--stop` exchanges, and its restart from the joint record, and compares them byte for byte (Joint restart, below). |
| `relax_factor_sweep.py` | Runs `nest.feedback = relax` at six factors on four nests and tabulates the errors and exchanges (Choice of the relaxation factor, below). |

The nest `hill` covers outer cells 8 to 24 in x and 6 to 18 in y (`nest.hill.lo = 400 300`,
`nest.hill.hi = 1200 900`): 32 x 24 x 16 cells of 25 m. The second joint run has two nests that touch at
x = 1000 m: `crest` (`300 300` to `1000 900`, the windward side and the crest) and `lee` (`1000 250` to `1500 950`).

## Data and acknowledgements
The hill and the two-nest layout are analytic. `relax_factor_sweep.py` also runs the Columbia Gorge and Askervein
cases; their terrain and observations are credited, with thanks to their providers, under "Data and
acknowledgements" in `Tests/Real_Terrain/Columbia_Gorge/README.md` and `Tests/Real_Terrain/Askervein_Nest/README.md`.
Nothing more is redistributed here.

## How to run
```bash
ctest --test-dir <build> -R rt_nests_two_way
```
By hand, in a scratch directory:
```bash
python3 check_nests_two_way.py -- mpiexec -n 2 <build>/gust inputs_outer
```
The script makes `ref/` (the 25 m run), `joint/` (the outer domain and the nest `hill`: one way, then jointly) and
`joint2/` (the same with the nests `crest` and `lee`).

## Expected Results
`check_nests_two_way.py` requires:
1. both runs exit 0; the reference reaches the steady verdict; the joint log has the joint run's start line, the
   largest |S| per exchange and `joint run: steady after N exchanges`;
2. two plotfiles each for the outer domain and the nest: the one-way ones, then the joint ones;
3. against the reference, the RMS of u:
   - over the nest (on the 25 m cells): joint at most 0.5 of one way;
   - over the outer cells outside the footprint (the reference averaged over 2 x 2 cells): joint at most 0.6 of one
     way;
   - with the two nests: over each nest, joint at most 0.7 of one way; outside both footprints, at most 0.9.

Measured on 2 ranks (RMS of each field against the 25 m run, m/s; k in m^2/s^2):

| | nest, one way | nest, joint | outer outside, one way | outer outside, joint | footprint, one way | footprint, joint |
|---|---|---|---|---|---|---|
| u | 8.23e-3 | **2.03e-3** | 1.44e-2 | **4.89e-3** | 2.89e-2 | **1.32e-2** |
| w | 3.21e-3 | **7.64e-4** | | | 8.38e-3 | **7.90e-3** |

The joint run reaches steady after 208 exchanges.

With the two nests (u unless stated):

| | nest, one way | nest, joint | footprint, one way | footprint, joint |
|---|---|---|---|---|
| crest | 6.11e-3 | **2.56e-3** | 2.66e-2 | **1.45e-2** |
| lee | 3.00e-2 | **1.82e-2** | 3.00e-2 | **1.79e-2** |
| lee, w | **6.42e-3** | 1.27e-2 | **6.72e-3** | 1.19e-2 |

Outside both footprints u moves from 8.53e-3 to **6.83e-3**; the run is steady after 254 exchanges. The lee nest's
w is further from the reference after the joint run: over the footprint the outer domain carries the w of the
restricted contravariant flux, not a cell mean of the nest's w, and in the hill's wake the two differ most.

Where the nest's downstream face stands matters. With the crest nest alone, ending at x = 800, 900, 1000, 1100 or
1200 m (crest at 800 m, half-width 150 m), the joint u error over the nest is 1.54, 3.31, 0.40, 0.27 and 0.22 times
the one-way one, and outside the footprint 1.59, 2.74, 0.66, 0.40 and 0.34 times. A face on the lee slope, where the
flow changes most, takes the 50 m solution there as its boundary values, and the joint run then moves both domains
away from the fine run. The two-nest layout therefore ends the crest nest 1.33 half-widths behind the crest.

## Feedback by relaxation (`--relax`)
`check_nests_two_way.py --relax` (CTest `rt_nests_relax`, PR-71) runs the reference, the defect run with the nest
`hill`, and the joint runs with `nest.feedback = relax` at `nest.feedback_factor` = 1, 4 and 16 (`relax1/`,
`relax4/`, `relax16/`; `docs/boundary_conditions.rst`, Two-way nests by relaxation). It requires:
1. every run exits 0 (the two-nest `joint2/` run is skipped), and each relax log's start line names its factor;
2. the outer domain minus the restricted nest where the relaxation acts (the log's last RMS of u) falls from factor 1
   to 4 to 16, and at 16 is at most 0.35 of its value at 1;
3. at factor 1, against the 25 m run, the RMS of u over the nest at most 0.5 of one way and over the outer cells
   outside the footprint at most 0.6;
4. old against new default (tracker N-114 item 13): a relax run without `nest.feedback_factor` (`relax_default/`)
   names 0.5 on its start line, and its RMS of u over the nest is at most 0.8 of the factor-1 run's (measured 0.694:
   2.35e-3 against 3.38e-3, 158 exchanges) and over the outer cells outside the footprint at most 1.05 (1.00).

Measured on 2 ranks (RMS of u, m/s; "outer - nest" where the relaxation acts):

| run | exchanges | outer - nest | nest | outer outside |
|---|---|---|---|---|
| one way | | | 8.23e-3 | 1.44e-2 |
| defect | 208 | | **2.03e-3** | **4.89e-3** |
| relax, f = 0.25 | 150 | 1.33e-2 | 2.32e-3 | 5.90e-3 |
| relax, f = 0.5 (default) | 158 | 1.07e-2 | 2.35e-3 | 5.02e-3 |
| relax, f = 1 | 160 | 8.04e-3 | 3.38e-3 | 5.03e-3 |
| relax, f = 2 | 159 | 5.52e-3 | 5.60e-3 | 6.14e-3 |
| relax, f = 4 | 177 | 3.49e-3 | 9.94e-3 | 9.54e-3 |
| relax, f = 16 | 508 | 1.50e-3 | 3.27e-2 | 2.94e-2 |

(0.25 and 2 are not in the test and come from PR-71, #116; the other rows are from this branch, where f = 1, 4 and
16 reproduce PR-71 and f = 0.5 is the `relax_default/` run.) The difference between the two domains falls with the factor, but above
f = 1 both move away from the fine run: the relaxation's target carries no nest pressure gradient, so as it tightens
the outer velocity under the nest stops responding to the outer pressure, and the nest, which takes its faces from
it, follows. The exchanges grow with f for the same reason: the flow rate through the footprint is set by a loop
(the nest returns the flow rate it is given, the relaxation imposes it) whose gain tends to 1. The bounds on the
errors therefore apply at f = 1 only.

## Pressure-consistent relaxation (`--relax-consistent`)
`check_nests_two_way.py --relax-consistent` (CTest `rt_nests_relax_consistent`, tracker N-114 item 12) runs the
reference, the defect run with the nest `hill`, `nest.feedback = relax` at factor 1 (`relax1/`) and
`nest.feedback = relax_consistent` at factors 1 and 4 (`rc1/`, `rc4/`; `docs/boundary_conditions.rst`,
Pressure-consistent relaxation). It requires:
1. every run exits 0, and each relax_consistent log's start line ends `with the weighted residual of the target`;
2. at factor 1 the outer domain minus the restricted nest (the log's last RMS of u) at most 0.25 of relax's, and
   smaller at factor 4 than at 1;
3. at factor 1, against the 25 m run, the RMS of u over the nest at most 0.6 of the one-way one and over the outer
   cells outside the footprint at most 0.45.

Measured on 2 ranks (RMS of u, m/s):

| run | exchanges | outer - nest | nest | outer outside |
|---|---|---|---|---|
| one way | | | 8.23e-3 | 1.44e-2 |
| defect | 208 | | **2.03e-3** | **4.89e-3** |
| relax, f = 1 | 160 | 8.04e-3 | 3.38e-3 | 5.03e-3 |
| relax_consistent, f = 1 | 169 | 1.32e-3 | 4.24e-3 | 5.18e-3 |
| relax_consistent, f = 4 | 174 | 1.22e-3 | 1.07e-2 | 1.01e-2 |

The outer domain follows the nest about six times more closely than with relax, as the added term intends. Both
domains are slightly further from the fine run than with relax at f = 1, and further at f = 4: the coupling loop
described above is unchanged. The bounds on the errors apply at f = 1 only.

## Choice of the relaxation factor (`relax_factor_sweep.py`)
Tracker row N-114, item 13. In PR-71 (#116) the hill alone suggested a factor below 1. `relax_factor_sweep.py`
runs `nest.feedback = relax` at `nest.feedback_factor` = 0.1, 0.25, 0.5, 0.75, 1 and 2 on four nests, each beside
its one-way run (`nest.feedback = none`), its defect run and its uniform fine run:

| case | layout | ranks |
|---|---|---|
| hill | this case, the nest `hill` | 2 |
| lee | this case, the nest `crest` alone ending at x = 900 m, on the lee slope (`--feedback-flux`, above) | 2 |
| gorge | `Tests/Real_Terrain/Columbia_Gorge`, the decks `prepare` writes | 4 |
| askervein | `Tests/Real_Terrain/Askervein_Nest`, the `two_way` layout (20 m nest of half-width 1000 m) | 4 |

Every run uses the recommended setting (`docs/turbulence_model.rst`): `turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, `turb.max_geom_lscale = 100`, `convection.hybrid_sigma_max = 0.1` and
`convection.hybrid_height` = 2.5 terrain heights. The Gorge and Askervein decks set it already; the script adds it
to the hill and lee runs, so their numbers differ from the tables above, which use the deck as it is. The nest is
compared with the fine run at its own cell heights (`tools/speed_difference.py`), the outer domain outside the
footprint with the fine run's 2 x 2 means.

```bash
python3 relax_factor_sweep.py run --gust <build>/gust --work <dir>
python3 relax_factor_sweep.py report --work <dir> --plot relax_factor_sweep.png
```
`run` takes `--cases` and `--factors` to run a part, and `--keep` to skip runs that have their verdict. `report`
prints the table below and writes it to `<dir>/relax_factor_sweep.md` and `.csv`.

Measured on macOS arm64 (AppleClang, double precision). Every run reached the steady verdict. RMS against the fine
run (m/s; k in m^2/s^2); in brackets, u over the one-way run's u:

| case | run | exchanges | nest u | nest v | nest w | nest k | outer u | outer v | outer w | outer k |
|---|---|---|---|---|---|---|---|---|---|---|
| hill | one way |  | 7.46e-03 | 5.13e-03 | 2.97e-03 | 9.94e-04 | 0.0131 | 5.31e-03 | 2.45e-03 | 4.09e-03 |
| hill | defect | 195 | 2.24e-03 (0.30) | 1.60e-03 | 7.68e-04 | 6.47e-04 | 4.11e-03 (0.31) | 2.44e-03 | 1.56e-03 | 2.17e-03 |
| hill | relax, f = 0.1 | 139 | 3.47e-03 (0.47) | 2.72e-03 | 1.27e-03 | 7.12e-04 | 7.18e-03 (0.55) | 3.53e-03 | 2.00e-03 | 3.15e-03 |
| hill | relax, f = 0.25 | 138 | 2.19e-03 (0.29) | 1.83e-03 | 7.89e-04 | 6.35e-04 | 4.86e-03 (0.37) | 2.80e-03 | 1.83e-03 | 2.58e-03 |
| hill | **relax, f = 0.5** | 146 | 2.58e-03 (0.35) | 1.75e-03 | 8.53e-04 | 6.46e-04 | 4.20e-03 (0.32) | 2.50e-03 | 1.70e-03 | 2.26e-03 |
| hill | relax, f = 0.75 | 147 | 3.18e-03 (0.43) | 1.96e-03 | 9.42e-04 | 6.76e-04 | 4.20e-03 (0.32) | 2.46e-03 | 1.61e-03 | 2.13e-03 |
| hill | relax, f = 1 | 151 | 3.79e-03 (0.51) | 2.25e-03 | 1.04e-03 | 7.12e-04 | 4.40e-03 (0.34) | 2.52e-03 | 1.57e-03 | 2.07e-03 |
| hill | relax, f = 2 | 160 | 6.14e-03 (0.82) | 3.53e-03 | 1.49e-03 | 8.91e-04 | 5.86e-03 (0.45) | 3.09e-03 | 1.61e-03 | 2.00e-03 |
| lee | one way |  | 6.66e-03 | 3.91e-03 | 3.67e-03 | 1.16e-03 | 0.0152 | 6.64e-03 | 3.54e-03 | 6.03e-03 |
| lee | defect | 204 | 0.023 (3.46) | 0.0144 | 0.0155 | 5.44e-03 | 0.0448 (2.95) | 0.0115 | 6.52e-03 | 0.0101 |
| lee | relax, f = 0.1 | 131 | 6.12e-03 (0.92) | 3.59e-03 | 3.55e-03 | 1.12e-03 | 0.0127 (0.84) | 6.17e-03 | 3.23e-03 | 6.02e-03 |
| lee | relax, f = 0.25 | 137 | 5.71e-03 (0.86) | 3.35e-03 | 3.62e-03 | 1.09e-03 | 0.0105 (0.69) | 5.82e-03 | 3.12e-03 | 6.05e-03 |
| lee | **relax, f = 0.5** | 145 | 5.85e-03 (0.88) | 3.49e-03 | 4.15e-03 | 1.21e-03 | 9.83e-03 (0.65) | 5.64e-03 | 3.26e-03 | 6.13e-03 |
| lee | relax, f = 0.75 | 148 | 6.37e-03 (0.96) | 3.83e-03 | 4.75e-03 | 1.38e-03 | 0.0108 (0.71) | 5.62e-03 | 3.53e-03 | 6.20e-03 |
| lee | relax, f = 1 | 148 | 6.91e-03 (1.04) | 4.16e-03 | 5.23e-03 | 1.52e-03 | 0.0119 (0.79) | 5.64e-03 | 3.84e-03 | 6.25e-03 |
| lee | relax, f = 2 | 151 | 8.25e-03 (1.24) | 4.95e-03 | 6.22e-03 | 1.79e-03 | 0.0148 (0.98) | 5.85e-03 | 4.88e-03 | 6.30e-03 |
| gorge | one way |  | 0.122 | 0.0331 | 0.0212 | 0.0278 | 0.102 | 0.0336 | 0.0256 | 0.0362 |
| gorge | defect | 279 | 0.12 (0.98) | 0.0301 | 0.018 | 0.0308 | 0.0877 (0.86) | 0.031 | 0.0248 | 0.0297 |
| gorge | relax, f = 0.1 | 212 | 0.121 (0.99) | 0.0291 | 0.0186 | 0.031 | 0.0909 (0.90) | 0.0323 | 0.0244 | 0.0317 |
| gorge | relax, f = 0.25 | 163 | 0.121 (0.99) | 0.029 | 0.0185 | 0.0311 | 0.0893 (0.88) | 0.0318 | 0.0243 | 0.0306 |
| gorge | **relax, f = 0.5** | 176 | 0.121 (0.99) | 0.029 | 0.0186 | 0.0314 | 0.0883 (0.87) | 0.0314 | 0.0244 | 0.0301 |
| gorge | relax, f = 0.75 | 175 | 0.122 (1.00) | 0.0293 | 0.0189 | 0.0316 | 0.0881 (0.87) | 0.0313 | 0.0245 | 0.03 |
| gorge | relax, f = 1 | 184 | 0.123 (1.01) | 0.0296 | 0.0191 | 0.0317 | 0.0882 (0.87) | 0.0313 | 0.0247 | 0.0299 |
| gorge | relax, f = 2 | 242 | 0.127 (1.04) | 0.0314 | 0.0204 | 0.0323 | 0.09 (0.89) | 0.0316 | 0.0256 | 0.0302 |
| askervein | one way |  | 3.37e-03 | 5.68e-03 | 2.87e-03 | 1.56e-03 | 4.62e-03 | 9.37e-03 | 2.57e-03 | 3.73e-03 |
| askervein | defect | 156 | 3.37e-03 (1.00) | 6.03e-03 | 1.48e-03 | 2.76e-03 | 3.21e-03 (0.69) | 6.93e-03 | 2.51e-03 | 3.12e-03 |
| askervein | relax, f = 0.1 | 87 | 2.81e-03 (0.83) | 4.74e-03 | 1.28e-03 | 2.61e-03 | 3.16e-03 (0.68) | 6.72e-03 | 2.44e-03 | 3.14e-03 |
| askervein | relax, f = 0.25 | 113 | 2.85e-03 (0.85) | 4.77e-03 | 1.26e-03 | 2.61e-03 | 3.17e-03 (0.69) | 6.71e-03 | 2.44e-03 | 3.10e-03 |
| askervein | **relax, f = 0.5** | 143 | 3.04e-03 (0.90) | 5.15e-03 | 1.34e-03 | 2.63e-03 | 3.32e-03 (0.72) | 6.91e-03 | 2.46e-03 | 3.10e-03 |
| askervein | relax, f = 0.75 | 168 | 3.27e-03 (0.97) | 5.69e-03 | 1.44e-03 | 2.65e-03 | 3.52e-03 (0.76) | 7.20e-03 | 2.48e-03 | 3.13e-03 |
| askervein | relax, f = 1 | 191 | 3.53e-03 (1.05) | 6.33e-03 | 1.54e-03 | 2.67e-03 | 3.75e-03 (0.81) | 7.56e-03 | 2.50e-03 | 3.16e-03 |
| askervein | relax, f = 2 | 272 | 4.69e-03 (1.39) | 9.30e-03 | 1.96e-03 | 2.79e-03 | 4.96e-03 (1.07) | 9.52e-03 | 2.60e-03 | 3.34e-03 |

The Gorge runs at f = 1 and with `defect`, and the Askervein runs at f = 1 and with `defect`, reproduce the
exchanges and errors of their case READMEs.

**What the sweep shows.**
- f = 0.5 against f = 1: fewer exchanges on all four cases (146 against 151, 145 against 148, 176 against 184, 143
  against 191). u over the nest is lower on all four (by 32, 15, 1 and 14 %); u over the outer domain outside the
  footprint is lower on three and equal on the Gorge (0.0883 against 0.0882).
- Not every field improves. On the hill, the outer domain's w and k outside the footprint are 8 and 9 % higher at
  f = 0.5 than at 1 (1.70e-3 against 1.57e-3, 2.26e-3 against 2.07e-3). On the Gorge, k over the outer domain is
  under 1 % higher.
- f = 0.25 is closer still on the lee layout and on Askervein, but on the hill the outer u error is 0.37 of one way
  (0.32 at f = 0.5). At f = 0.1 the hill's errors grow again (0.47 and 0.55 of one way) and the Gorge takes 212
  exchanges.
- Above f = 1 every case moves away from the fine run, and the Gorge and Askervein need 242 and 272 exchanges.
- Relax at f = 0.5 against `defect`, u: closer on the lee layout (where `defect` is 3 times one way); on Askervein
  closer over the nest (3.04e-3 against 3.37e-3) and further outside it (3.32e-3 against 3.21e-3); further on the
  hill; within 2 % on the Gorge.

On these four cases the table supports 0.5 over 1 as the default of `nest.feedback_factor` with relax, and 0.5 is
the default since this item; `nest.feedback_factor = 1` keeps the old answer, and `rt_nests_relax` runs both
(Feedback by relaxation, above). The relaxation's target has no nest pressure gradient (above); a target that
carries one (tracker N-114 item 12) removes the cause of the errors above f = 1 and could move the best factor up,
so that change should re-run this sweep.

Limits of the sweep:
- Each factor ran once; the exchanges vary by about 5 % from one factor to the next without a trend (the Gorge:
  163, 176, 175 at 0.25, 0.5, 0.75), so smaller differences in exchanges are not gains. 0.35 was not run.
- The outer domain is compared with the reference's 2 x 2 means level by level, without interpolating to the outer
  cell heights; where the two meshes follow the terrain differently this adds an error that is the same at every
  factor, and draws every ratio towards 1 (most on the Gorge).
- Not covered: two touching nests with relax, stratified decks (theta is relaxed at the same rate), canopy and
  Coriolis decks.

## Conservative faces (`--feedback-flux`)
`check_nests_two_way.py --feedback-flux` (CTest `rt_nests_feedback_flux`, tracker N-114 item 1) runs the reference
and three joint runs with `nest.feedback = defect`: the nest `crest` ending at x = 900 m, on the lee slope, with
`nest.feedback_flux = outer` (`flux_outer/`, the defect source alone) and `= nest` (`flux_nest/`), and the nest
`hill` with `= nest` (`flux_hill/`; `docs/boundary_conditions.rst`, Conservative faces). It requires:
1. every run exits 0; the line `nest.feedback_flux = nest: ... first max |correction|` is in the two `nest`
   logs, and their joint lines end with the ring's largest correction (`; ring max |C| u ...`); the `outer` log
   has neither;
2. on the crest nest, against the 25 m run, the RMS of u with `outer` at least 2.0 times the one-way one over the
   nest; with `nest` at most 1.6 times over the nest and 1.25 times over the outer cells outside it, and at most 0.6
   of the `outer` run's on both;
3. on the hill, with `nest`, at most 1.6 and 1.25 times the one-way one.

Measured on 2 ranks, the deck as it is (RMS against the 25 m run, m/s; in brackets the joint over the one-way
value). The row "hill, outer" is the joint run of `rt_nests_two_way` (above), not a run of this script:

| nest | `feedback_flux` | exchanges | nest u | nest w | nest k | outer u outside |
|---|---|---|---|---|---|---|
| crest to 900 m | outer | 206 | 2.37e-2 (3.31) | 1.62e-2 (4.02) | 5.41e-3 (4.47) | 4.55e-2 (2.74) |
| crest to 900 m | **nest** | 207 | **8.74e-3 (1.22)** | 6.95e-3 (1.73) | 2.15e-3 (1.78) | **1.60e-2 (0.96)** |
| hill | outer | 208 | 2.03e-3 (0.25) | 7.64e-4 (0.24) | 5.61e-4 (0.49) | 4.89e-3 (0.34) |
| hill | nest | 209 | 2.10e-3 (0.26) | 6.32e-4 (0.20) | 5.60e-4 (0.49) | 4.76e-3 (0.33) |

On the crest nest the faces cut the flow where it turns down the lee slope, and the defect source alone leaves the
joint run 3.3 times worse in u than one way. The conservative faces bring u to 1.2 times one way on the nest and
below it outside; w and k stay 1.7 to 1.8 times worse. On the hill, whose faces lie where the flow is level, u is
3 % worse and w 17 % better. The other layouts of `docs/boundary_conditions.rst` are not in the test.

## Feedback by injection (`--inject`)
`check_nests_two_way.py --inject` (CTest `rt_nests_inject`, tracker row N-114, item 2) runs the reference, the
defect run with the nest `hill`, and the joint run with `nest.feedback = inject` at the default
`nest.feedback_factor` of 0.5 (`inject/`; `docs/boundary_conditions.rst`, Two-way nests by injection). Every
exchange the outer state, from 2 outer cells inside the nest's faces, moves that fraction of the way to the
restricted nest and is held there through the outer iteration. It requires:
1. every run exits 0 (the two-nest `joint2/` run is skipped), and the inject log's start line names its factor;
2. against the 25 m run, the RMS of u over the nest at most 0.55 of one way and over the outer cells outside the
   footprint at most 0.42.

Measured on 2 ranks (RMS of u, m/s, against the 25 m run):

| run | exchanges | nest | outer outside |
|---|---|---|---|
| one way | | 8.23e-3 | 1.44e-2 |
| defect | 208 | **2.03e-3** | **4.89e-3** |
| relax, f = 1 | 160 | 3.38e-3 | 5.03e-3 |
| inject, f = 0.25 | 602 | 4.22e-3 | 5.41e-3 |
| inject, f = 0.5 | 344 | 4.22e-3 | 5.41e-3 |
| inject, f = 0.6 | 300 | 4.22e-3 | 5.41e-3 |
| inject, f = 0.75 | not steady after 3000 | | |
| inject, f = 1 | not steady after 600 | | |

(Only 0.5 is in the test.) Up to f = 0.6 injection reaches the same steady state: the factor sets the number of
exchanges, not the answer. At 0.75 and 1 the joint run oscillates and does not reach steady. On this hill its errors are larger than those of the defect run and of relax at f = 1,
and it takes more exchanges.

## Joint restart (`check_joint_restart.py`)
A joint run that stops after `steady.max_iterations` exchanges writes a joint record (`Joint`) in the outer
domain's last checkpoint and names that checkpoint in its error message (`docs/boundary_conditions.rst`, Two-way
nests, Restart). The script runs, each in its own directory:
1. `straight/`: the deck with the nests `crest` and `lee` and `nest.feedback = defect` (with `--relax`, `--relax-consistent` or
   `--inject`: the nest `hill` and `nest.feedback = relax`, `relax_consistent` or `inject`), to the steady verdict;
2. `stopped/`: the same with `steady.max_iterations = M` (`--stop M`), which stops after M exchanges;
3. `stopped/`: the restart, with the deck less its `init.*` lines and `io.restart` set to the named checkpoint;
4. with `--relax`, before 3, two restarts of one exchange with `nest.feedback_factor = 0.75`, each from a copy of
   `stopped/`: `other/` from the record as written, `older/` from the record rewritten as version 2 (no factor).

It requires:
1. the straight run and the restart exit 0, the stopped run does not; the record has the feedback, M exchanges, the
   factor (0 with defect) and one line per nest in `nest.names` order, each naming a checkpoint that exists;
2. the restart logs `restarted after exchange M`, and is steady after the straight run's exchanges with every
   domain's iterations over the run the same;
3. the last plotfile of every domain and the last checkpoint's state, node heights, closure history, nest faces and
   joint record are byte for byte the straight run's (the driver record holds wall times and is not compared);
4. with `--relax`, `other/` logs 0.5 and 0.75 and that the window fills again, `older/` logs that the record holds
   no factor and restores the window, and both write a version 3 record with M + 1 exchanges and the factor 0.75.

M must be above the one-way runs' iterations (`steady.max_iterations` caps them too: 99 for the outer domain, up to
121 for a nest here) and below the straight run's exchanges.

| | ranks | M | exchanges | result |
|---|---|---|---|---|
| crest and lee, defect | 3 | 150 | 254 | byte for byte |
| hill, relax (default factor 0.5) | 2 | 140 | 158 | byte for byte |
| hill, inject | 2 | 150 | 344 | byte for byte |

## CTest
`rt_nests_two_way` (2 ranks, labels `regression;smoke`). `rt_nests_joint_restart` (3 ranks, labels
`regression;smoke`) `rt_nests_joint_restart_relax` and `rt_nests_joint_restart_inject` (2 ranks, label `regression`) run `check_joint_restart.py`
(above). The abort tests `abort_nests_joint_restart_names`, `_iterations`, `_version` and `_missing` check that a
joint record naming other nests, a `steady.max_iterations` at or below its exchanges, a record of another version,
and a record naming a nest checkpoint that is not there stop the run and name `io.restart`. The abort tests `abort_nests_feedback_value`,
`abort_nests_feedback_no_nests`, `abort_nests_feedback_newton`,
`abort_nests_feedback_damping`, `abort_nests_feedback_odd_boxes` and `abort_nests_feedback_narrow` check that the inputs `nest.feedback` refuses stop
the run and name the key. `inputs_abort` is `inputs_outer` without `mesh.max_grid_size` and `steady.stage`, which
two of them set on the command line (a key given twice stops the run).

`rt_nests_relax` (2 ranks, label `regression`; above): about 50 s alone on an Apple M-series laptop (about 2 minutes
under `ctest -j 4`), most of it the f = 16 run (508 exchanges). `--relax` skips the two-nest `joint2/`
run. `abort_nests_feedback_factor`, `_factor_nan`, `_factor_unused` and `_relax_damping` check that
`nest.feedback_factor = 0`, `nan`, a factor with `nest.feedback = defect`, and `relax` with a Rayleigh layer on u, v
and k stop the run and name the key.

`rt_nests_feedback_flux` (2 ranks, label `regression`; above): about 1.5 minutes on an Apple M-series laptop.
`abort_nests_feedback_flux_value` and `_unused` check that `nest.feedback_flux = fine` and `nest` without
`nest.feedback = defect` stop the run and name the key.
`rt_nests_inject` (2 ranks, label `regression`; above): about 30 s on an Apple M-series laptop.
`abort_nests_inject_factor` checks that `nest.feedback_factor = 2` with `nest.feedback = inject` stops the run and
names the key; `abort_nests_feedback_factor_unused` covers a factor with `defect` for both relax and inject.
