# Nest_Placement_Runs

## Purpose
Joint runs (`nest.feedback = defect`) that test the step rule of `tools/nest_placement.py` and its verdicts for an
oblique wind (tracker N-114, item 11; beyond the plan). The tool's angle and misfit rules were fitted on a channel
site with the wind normal to the faces (`docs/boundary_conditions.rst`, Where to put a joint nest). The step rule
(escarpments, plateau edges, steep valley sides) and winds at an angle to the faces had no run behind them. Here
one escarpment, two nests and three winds give six pairs of a tool verdict and a joint run. Each joint run is also
repeated with the conservative nest faces (`nest.feedback_flux = nest`, item 1). The case checks:
1. that a nest the tool flags makes the joint run worse than the one-way run, and a nest it clears does not;
2. whether the conservative faces remove that harm.

## Files
| File | Role |
|---|---|
| `inputs_outer` | The outer domain: 60 x 48 x 16 cells of 50 m over 3000 m x 2400 m x 800 m (first layer 10 m), log-law inflow (u* 0.4 m/s, z0 0.1 m), neutral, Rayleigh layer above 400 m on w and theta, the recommended setting with `convection.hybrid_height` = 150 m (2.5 escarpment heights). No terrain file, wind, side faces or nest keys: the script gives them. |
| `nest_placement_runs.py` | Writes the terrain, runs the tool and the runs, and prints the table below (also written to `nest_placement_runs.csv`). |

The escarpment is z = (H/2)(1 + tanh(s/L)) with H = L = 60 m. Here s is the distance from the line through
(1500, 1200) m at 35° to the x axis. The script writes it to `escarpment.xyz` on points 12.5 m apart. The 50 m outer
cells miss the ground by at most 1.8 m. That is above the tool's misfit gate at this cell size (1.25 m). The nests:

| Nest | lo (m) | hi (m) | Faces that cut the step | Angle of the step to their normal |
|---|---|---|---|---|
| `xcut` | 1250 800 | 1750 1600 | west, east | 35° |
| `ycut` | 750 1000 | 2250 1400 | south, north | 55° |

The winds are given as `nest_placement.py` reads them (degrees from north, the direction the wind comes from):
- 270°: GUST `bc.inflow_direction = 0`, slip walls on the y faces;
- 180°: `= 90`, slip walls on the x faces;
- 225°: `= 45`, inflow on the west and south faces, outflow on the east and north faces.

## Data
None. The case is analytic.

## How to run
From an empty directory, with `<repo>` the GUST checkout and `<build>` its build directory:
```bash
python3 <repo>/Tests/Canonical_RANS/Nest_Placement_Runs/nest_placement_runs.py -- mpiexec -n 2 <build>/gust <repo>/Tests/Canonical_RANS/Nest_Placement_Runs/inputs_outer
```
`--only 270 180` runs some winds only. For each wind the script runs:
- the uniform 25 m reference (`ref_<wind>/`);
- for each nest, the joint run (`<nest>_<wind>/`);
- for each nest, the joint run with `nest.feedback_flux = nest` (`<nest>_<wind>_flux/`).

That makes 15 runs. A folder whose log already has the steady verdict is reused. The script has no verdict of its
own. It exits 0 when every run exits 0 and is steady.

## Expected Results
Measured on 2 ranks (macOS arm64, AppleClang). Each entry is the RMS against the 25 m run, joint over one way:
- the nest is compared at its own cell heights;
- the outer domain's cells outside the footprint are compared against the reference's 2 x 2 means;
- each entry reads "defect alone / with the conservative faces";
- bold marks a defect run more than 5 % worse than one way.

The last column counts the outer cells of the lowest 200 m whose horizontal speed differs from the reference's by
more than 0.25 m/s: one way → defect / conservative faces. No nest cell differs by more than 0.25 m/s in any run.

| wind from | nest | tool | exchanges | nest u | nest v | nest w | nest k | outer u | outer v | outer w | outer k | outer cells > 0.25 m/s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 270° | `xcut` | flagged (east, west) | 170 / 168 | 1.04 / 0.96 | 1.03 / 0.90 | 0.99 / 0.96 | 1.03 / 0.98 | **1.10** / 1.01 | **1.10** / 1.00 | **1.71** / 1.00 | 1.04 / 1.00 | 16 → 16 / 16 |
| 270° | `ycut` | clear | 180 / 183 | 1.02 / 1.03 | 0.99 / 0.97 | 0.96 / 0.93 | 1.01 / 1.02 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 16 → 16 / 16 |
| 180° | `xcut` | clear | 198 / 200 | 0.85 / 0.84 | 0.89 / 0.90 | 0.77 / 0.76 | 0.93 / 0.91 | 0.95 / 0.96 | 1.00 / 1.00 | 1.01 / 1.01 | 0.99 / 0.99 | 5 → 5 / 5 |
| 180° | `ycut` | flagged (north, south) | 208 / 211 | **1.21** / 1.02 | **1.35** / 0.98 | **1.29** / 0.98 | **1.11** / 1.00 | **1.14** / 0.97 | **1.41** / 1.08 | **1.75** / 1.02 | **1.16** / 1.03 | 5 → 11 / 5 |
| 225° | `xcut` | flagged (east, west) | 173 / 165 | 0.88 / 0.90 | 0.94 / 0.95 | 0.94 / 0.93 | 0.96 / 0.96 | 0.95 / 0.95 | 0.94 / 0.94 | 0.98 / 0.97 | 0.95 / 0.95 | 0 → 0 / 0 |
| 225° | `ycut` | flagged (north, south) | 169 / 171 | 0.95 / 0.92 | 1.03 / 0.94 | 0.91 / 0.87 | 0.99 / 0.97 | 0.95 / 0.94 | 0.94 / 0.93 | 0.97 / 0.95 | 0.94 / 0.94 | 0 → 0 / 0 |

One-way RMS of the streamwise component, for scale: u over the nest 2.0e-2 (`xcut`) and 5.1e-3 m/s (`ycut`) at
270°; v over the nest 6.4e-3 and 1.9e-2 m/s at 180°.

**Winds along an axis.** The tool's verdict matches the runs. The two nests it flags make the joint run worse than
the one-way run. The harm is larger where the step meets the faces at 55° (wind from 180°) than at 35° (wind from
270°):
- at 55°: up to 1.35 in the nest and 1.75 outside it, and 5 → 11 outer cells above 0.25 m/s;
- at 35°: at most 1.10 in the outer u and 1.71 in the outer w.

The two nests it clears stay within 1.02 of one way or improve. With the conservative faces the two flagged nests
stay within 1.02 of one way in the nest and 1.08 outside it. Their outer cells above 0.25 m/s stay at 16 and 5.

**Oblique wind.** For the wind from 225° the tool flags both nests: every face is inflow or outflow, and both nests
cut the step at more than 20°. Neither joint run is worse than one way (at most 1.03, mostly 0.88 to 0.98). This wind
blows within 10° of the step's contour lines, so little flow crosses the step at the faces. The tool measures the
angle against the face normal, not the wind, so it does not see this. The flag is then on the safe side.

**Thresholds.** The tool's thresholds are not changed. The runs confirm the flags for the winds across the step, and
the one case where a flag was not needed errs on the safe side. The effects are small here because the misfit
(1.8 m) is just above the gate (1.25 m). With the conservative faces the exchanges change by at most 8 (165 to 211
in all runs). Wall times are not reported: the runs shared the machine.

## CTest
None. The case is a study of 15 runs with no verdict of its own; the table above is its record. The tool itself is
tested by `nest_placement_selftest`.
