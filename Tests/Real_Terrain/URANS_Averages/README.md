# URANS_Averages

## Purpose
The time average of an unsteady RANS run against the steady answer of the same deck,
on two cases that have a reference:

- RUSHIL Hill 5 (`../RUSHIL`), against the wind tunnel, with the default closure (steady) and with a setting
  whose steady run stops at the iteration limit;
- Askervein with ERF's setup (`../Askervein_210`), against ERF's unsteady RANS averaged over 2400–4800 s on
  lines A and AA.

The folder holds no deck of its own: the unsteady deck is written from the steady deck of the case.

## Files
| File | Role |
|---|---|
| `urans_averages.py` | `deck` (the unsteady deck of a steady deck; the keys an unsteady run refuses are dropped and printed), `cost` (iterations of a steady log, with those on the coarser levels of a sequenced run; steps and inner iterations of an unsteady log), `fields` (an unsteady plotfile against a steady plotfile on the same mesh: largest and RMS difference of the running means `u_mean`, `v_mean` from `u`, `v`; largest difference of the last step; largest standard deviation of u, from `u_var`). |
| `check_mean_option.py` | The CTest's check of the time-mean option: the column reader (`tools/askervein_lines.py`, `Columns(mean=True)`) and the three command lines that forward to it. |

The comparisons with the references use the case's own tools with `--mean`, which reads `u_mean`, `v_mean` in
place of `u`, `v`: `../RUSHIL/rushil.py compare|plot --mean` and `tools/askervein_lines.py sample|check --mean`.
Their output says that it holds the time mean. A plotfile written before `unsteady.average_start` (means of zero)
is refused.

## Data and acknowledgements
Nothing is downloaded for this folder beyond what the two cases need, and no data are added to the repository.

- Hill 5: the RUSHIL measurements of Khurshudyan, L.H., Snyder, W.H., Nekrasov, I.V. (1981), "Flow and dispersion
  of pollutants over two-dimensional hills", EPA-600/4-81-067, as distributed by the ERCOFTAC Knowledge Base
  (Classic Collection case 69 / Application Challenge AC5-05). The user downloads the files; they are not in the
  repository (`../RUSHIL/README.md`); the files carry no licence statement. We thank ERCOFTAC and the QNET-CFD Knowledge Base for keeping the data
  available, and the US EPA Fluid Modeling Facility for the measurements.
- Askervein: the reference is ERF's own run (`../Askervein_210/Gold/erf_urans_lines_10m.csv`, in the repository;
  `../Askervein_210/README.md`). The mast positions of lines A and AA (Taylor, P.A., Teunissen, H.W. (1987), "The
  Askervein Hill project: overview and background data", *Boundary-Layer Meteorology* 39, 15) are as tabulated in
  N. Wagenbrenner's windninja-tools, used with her permission. We thank the ERF developers for the code and its
  Askervein case, and N. Wagenbrenner for the tabulated positions.

## How to run
A key given in the deck and again on the command line aborts, so the two runs of a pair go in two directories
(`s`, `u`) and keep the deck's plotfile names.

Hill 5, from the directory `rushil` that `rushil.py prepare` makes (`../RUSHIL/README.md`). The default closure
is the RUSHIL deck with `turb.realizability = 0`, `turb.nut_relax = 1` and no blending (the deck's
`turb.max_geom_lscale = 100` is kept). The stalled setting is `turb.realizability = 0.33`, the deck's
`turb.nut_relax = 0.5` and no blending, with the steady run's limit set to 4000 iterations (the deck's is 20000):
replace the first two `sed` expressions by `-e 's/^turb.realizability.*/turb.realizability = 0.33/'` and
`-e 's/^steady.max_iterations.*/steady.max_iterations = 4000/'`.

```bash
T=../Tests/Real_Terrain
G="terrain.file=../hill5.xyz bc.inflow_file=../inflow.txt"
mkdir s u
sed -e 's/^turb.realizability.*/turb.realizability = 0/' -e 's/^turb.nut_relax.*/turb.nut_relax = 1/' \
    -e 's/^convection.hybrid_sigma_max.*/convection.hybrid_sigma_max = 0/' -e '/^convection.hybrid_height/d' \
    $T/RUSHIL/inputs_rushil > s/inputs
python3 $T/URANS_Averages/urans_averages.py deck s/inputs u/inputs
(cd s && mpirun -np 8 ../../build-rel/gust inputs $G > log.txt)
(cd u && mpirun -np 8 ../../build-rel/gust inputs $G unsteady.stop_time=4800 unsteady.dt=2.5 \
    unsteady.average_start=2400 unsteady.inner_drop_orders=2 > log.txt)
python3 $T/RUSHIL/rushil.py compare s/plt --hill 5 --outdir .
python3 $T/RUSHIL/rushil.py compare u/plt --hill 5 --outdir . --mean
python3 $T/URANS_Averages/urans_averages.py fields s/plt u/plt
python3 $T/URANS_Averages/urans_averages.py cost s/log.txt u/log.txt
```

Askervein, from a scratch directory `run` in the repository root:

```bash
A=../Tests/Real_Terrain/Askervein_210
G="terrain.file=../../extern_erf/askervein.txt bc.inflow_file=../$A/askervein_210_inflow.txt"
mkdir s u
cp $A/inputs_steady_askervein_210 s/inputs
python3 ../Tests/Real_Terrain/URANS_Averages/urans_averages.py deck s/inputs u/inputs
(cd s && mpirun -np 8 ../../build-rel/gust inputs $G > log.txt)
(cd u && mpirun -np 8 ../../build-rel/gust inputs $G unsteady.stop_time=4800 unsteady.dt=4 \
    unsteady.average_start=2400 unsteady.inner_drop_orders=2 > log.txt)
python3 ../tools/askervein_lines.py check u/steady_askervein_210_plt $A/Gold/erf_urans_lines_10m.csv \
    --terrain ../extern_erf/askervein.txt --lines $A/askervein_lines.csv --out lines_mean.csv \
    --lee_d 300 --lee_rel_tol 0.10 --mean
python3 ../Tests/Real_Terrain/URANS_Averages/urans_averages.py fields s/steady_askervein_210_plt u/steady_askervein_210_plt
python3 ../Tests/Real_Terrain/URANS_Averages/urans_averages.py cost s/log.txt u/log.txt
```

## Expected Results
Measured on 8 ranks of one laptop (Apple M2 Pro). Costs are iterations; wall times are not reported, and the cost
of one inner iteration of the unsteady run was not compared with that of one steady iteration.

Every unsteady run starts from the deck's first guess (the inflow profile), not from the steady field, without
grid sequencing, with BDF1 in time. A step stops two orders below its first residual or at
`steady.residual_tol`, whichever comes first.

**Hill 5** (936 × 8 × 64 cells; unsteady: 1920 steps of 2.5 s to 4800 s, means over 2400–4800 s, 960 samples).

| Setting | Run | Verdict or steps | Iterations | Reversed flow at 0.1 H | U rms vs tunnel |
|---|---|---|---|---|---|
| Default closure | steady | steady | 248 (and 42 + 8 on the two coarser levels) | none | 0.225 m/s |
| Default closure | unsteady, time mean | 1920 steps | 5354 (2.79 per step) | none | 0.225 m/s |
| `turb.realizability = 0.33`, `turb.nut_relax = 0.5`, no blending | steady | not steady at the limit of 4000 | 4000 (and 45 + 9) | none | 0.214 m/s |
| The same | unsteady, time mean | 1920 steps | 5438 (2.83 per step) | none | 0.214 m/s |

- The unsteady run reaches the steady field. At 4800 s it differs from the steady field by at most 1.7e-4 m/s
  (0.004 % of the largest speed, 3.93 m/s) with the default closure and by 1.3e-3 m/s (0.03 %) with the stalled
  setting. From about 4200 s (default) and 4300 s (stalled) every step is one iteration: the state is within
  `steady.residual_tol` in every field.
- The time mean is further from the steady field than the last step is: at most 8.1e-3 m/s (0.21 %; RMS
  1.7e-3 m/s) with the default closure and 1.1e-2 m/s (0.28 %; RMS 2.1e-3 m/s) with the stalled setting, both in
  the column next to the outflow face. The mean includes the end of the start-up: 2400 s is one flow-through time
  at the domain top but 0.6 of one at 57 m above the ground, where the largest difference is. A later
  `unsteady.average_start` would bring the mean to the steady field; that run was not made.
- The standard deviation of u over 2400–4800 s is zero upstream and grows to 0.009 m/s (default closure) and
  0.013 m/s (stalled setting) in the outflow column. It is this decaying start-up, not an oscillation: the last
  step lies on one side of the mean almost everywhere.
- Stalled setting: the steady run's residuals repeat with a period of 30 iterations from about iteration 330 (the
  k residual between 5.1e-5 and 8.0e-5; `steady.residual_tol` is 5e-5). The field it stops on is the unsteady
  run's final field to 0.03 %. Both have a layer of reversed flow in the ground cells from 0.8a to 1.0a behind the
  summit (at most −0.024 m/s), below 0.1 H.
- The unsteady run to 4800 s takes 22 times the iterations of the steady run that reaches the verdict (5354
  against 248). The ratio depends on the stop time.
- Limit of these runs: BDF1 with a step of 2.5 s (Δt·U/Δx about 1) damps an oscillation of period 195 s (a
  Strouhal number of 0.2 on the hill height) by 22 % per period, and shorter periods more. A weakly unstable mode
  could be suppressed at this step. No run with a smaller step was made.

**Askervein, ERF's setup** (300 × 300 × 18 cells; unsteady: 1200 steps of 4 s to 4800 s, means over 2400–4800 s,
the window of ERF's averages).

| Run | Verdict or steps | Iterations | Lines A and AA against ERF, to 300 m past the hilltop | Lee masts (d > 300 m) |
|---|---|---|---|---|
| steady | steady | 135 | at most 4.31 % (AA, d = 284 m) | 6.86 to 9.85 % slower |
| unsteady, time mean | 1200 steps | 2248 (1.87 per step) | at most 4.31 % | 6.86 to 9.85 % slower |

- The time mean and the steady field differ by at most 9.8e-4 m/s (0.004 % of the largest speed, 22.3 m/s; RMS
  6.6e-5 m/s), in the outflow corner. At the 30 masts the two give the same speed to the three decimals printed.
- From step 210 (840 s) every step is one iteration. The standard deviation of u over 2400–4800 s is at most
  1e-11 m/s: the unsteady run is steady well before the averaging starts.
- Against ERF's unsteady RANS the 5 % target is met at the 26 masts up to 300 m past the
  hilltop and missed at the four lee masts, by the same amount in both runs: line A d = 393 m −9.85 %; line AA
  d = 382, 477, 575 m −6.91, −8.19, −6.86 %. The time average does not change the lee difference. It is the
  difference recorded for the steady run in `../Askervein_210/README.md`, with the reasons given there; the check
  above passes with that folder's 10 % lee tolerance.
- The unsteady run takes 17 times the iterations of the steady run (2248 against 135).

**Laminar cylinder** (`../../Canonical_RANS/Cylinder_Laminar/README.md`, 2 ranks; measured with that case). At
Re = 40 the steady run takes 733 iterations and the unsteady run 4430 (2400 steps, 1.8 per step) for the same
drag coefficient to four digits. At Re = 100 the steady run stops at its limit of 6000 with no verdict and the
unsteady run takes 14 449 (3200 steps, 4.5 per step) and sheds.

## CTest
`rt_urans_mean_option` (with `rt_urans_mean_option_clean` and `rt_urans_mean_option_run`; `regression;smoke`,
double precision): five steps of 60 s of `../../Canonical_RANS/Unsteady_Hill_2D` with two cells in y, on 2 ranks,
with `unsteady.average_start = 0` and a plotfile every step.

- Reader: the mean read with the option equals the arithmetic mean of the five steps' fields in three columns
  (measured 1.8e-15 m/s; limit 1e-10 of the largest |u|) and differs from the last step by more than 1e-3 of the
  largest |u| (measured 0.207 m/s of 8.03 m/s). A mean that is not finite fails.
- Command lines: `askervein_lines.py sample --mean` at two masts gives the speed of those mean columns (4e-7 m/s,
  the CSV's six decimals) and differs from plain `sample` (0.084 m/s); `check --mean` passes against those speeds
  and plain `check` fails; `rushil.py compare --mean` prints the reference speed of the mean columns (4.049 m/s;
  4.001 m/s without the option).

Mutants that fail: the option ignored in the reader (`u`, `v` read); `u_mean` read for `v_mean`; `mean` not
passed on in `sample`, in `check`, in `rushil.py compare`.

The Hill 5 and Askervein runs above are not CTests: Hill 5 needs the downloaded tunnel data, and the Askervein
unsteady run is 1200 steps on the full mesh.
