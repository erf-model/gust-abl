# Nest_Faces

## Purpose
A nested patch whose four lateral faces take their values from the parent run (`bc.xlo`, `bc.xhi`, `bc.ylo`,
`bc.yhi = nest`; sub-phase 9.2, beyond plan rev. 6). The parent runs to steady on 50 m cells; the patch refines its
centre at 2 2 1 (25 m cells), starts from the parent's checkpoint and runs to steady with the parent's normal velocity
on its faces and the parent's u, v, w, theta and k in its ghost cells. The wind is from 20 degrees off the x axis, so
every face of the patch has the parent's flow entering on one part or leaving on another. Two grounds: flat (the patch
must keep the parent's flow) and a 60 m round Witch of Agnesi hill (the finer cells change the flow over the hill).
The relaxation zone (`bc.nest_relax_width`), the ground blend (`terrain.nest_blend_width`) and the restart record are
exercised too.

## Files
| File | Role |
|---|---|
| `inputs_flat_parent`, `inputs_hill_parent` | The parents: 32 x 24 x 16 cells, 1600 m x 1200 m x 800 m, periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m, 20 degrees) on xlo, outflow on xhi, neutral. |
| `inputs_flat_patch`, `inputs_hill_patch` | The patches: 32 x 24 x 16 cells from (400, 300) m to (1200, 900) m (parent cells 8 to 23 in x, 6 to 17 in y), four nest faces, `init.first_guess = nest`, boxes of at most 13 x 9 cells. |
| `inputs_flat_patch_open` | The flat patch with `bc.xlo = nest`, `bc.xhi = outflow` and periodic y: parent cells 8 to 23 in x and the parent's whole width in y (32 x 48 x 16 cells). |
| `inputs_hill_parent_pair` | The hill parent without `bc.inflow_ustar` or an iteration count: two parents on the same mesh at the same step with different flows (u* 0.4 and 0.6 m/s, 30 iterations), for the restart tests. |
| `inputs_hill_patch_restart` | The hill patch without a first guess, box size, parent or iteration count: the restart tests give them on the command line. |
| `check_nest_faces.py` | Runs a parent or a patch and checks it (below). |

The decks do not name their outputs; `check_nest_faces.py` appends `io.plot_file` and `io.chk_file`, and copies the
last plotfile and checkpoint to `<prefix>_plt_final` and `<prefix>_chk_final` (the names the patch decks use for
`bc.nest_parent`).

## Data
None. The cases are analytic.

## How to run
```bash
ctest --test-dir <build> -R "nest_faces"
```
`check_nest_faces.py open ...` takes the arguments of `patch` for a patch with an outflow face (check 1 then requires
the net flux to be reported as not removed).

By hand, in a scratch directory:
```bash
python3 check_nest_faces.py parent hill_parent -- mpiexec -n 2 <build>/gust inputs_hill_parent
python3 check_nest_faces.py patch hill_parent_plt_final hill_np2 8 6 0.15 0.06 1e-8 1e-3 -- mpiexec -n 2 <build>/gust inputs_hill_patch
```

## Expected Results
`check_nest_faces.py patch` requires:
1. the log names the nest faces and reports the net flux through them as removed before the pressure solve, with
   |net flux| / (flux entering) in the range given;
2. the divergence after the first pressure correction is at most 1e-8 x the divergence before. With zero-gradient p'
   on all four faces this needs a boundary flux with no net;
3. the steady verdict;
4. the mean of u, v, w and k over the 2 x 2 patch cells of every parent cell is within the tolerances of the parent
   cell's value.

Measured (macOS arm64, AppleClang, double precision):

| Case | Iterations to steady | Net flux / flux entering | Divergence before, after | max abs(2 x 2 mean - parent): u, v, w (m/s), k (m2/s2) |
|---|---|---|---|---|
| Flat parent | 79 | | | |
| Flat patch, 3 ranks | 64 | 2.0e-13 | | 2.1e-4, 9.5e-5, 9.7e-5, 3.5e-4 |
| Flat patch, zone of 4 cells | 63 | | | within the same tolerances |
| Flat patch, nest xlo, outflow xhi, periodic y, 3 ranks | 80 | -1.000, not removed | 1.8e-4, 9.0e-15 | 1.6e-3, 5.7e-4, 1.0e-3, 3.8e-4 |
| Hill parent | 84 | | | |
| Hill patch, 1 / 2 / 4 ranks | 103 / 102 / 102 | 6.3e-6 | 4.5e-2, 1.7e-12 | 6.2e-2, 2.7e-2, 1.6e-2, 2.4e-2 |
| Hill patch, zone of 4 cells and blend of 150 m | 104 | | | within the same tolerances |

- Flat ground: the patch differs from its parent by the steady tolerance of the two runs (tolerances 1e-3 m/s and
  2e-3 m2/s2). The parent's flow is uniform along x and y, so the interpolated flux has no net beyond round-off.
- Nest with an outflow face: all the flux enters through the one nest face and the outflow face carries it, so nothing
  is removed. The largest differences are in the cells beside the outflow face, which the parent does not have there
  (tolerances 5e-3 m/s and 2e-3 m2/s2).
- Hill: the interpolated face velocities have a net of 6.3e-6 of the flux entering; 1.5e-5 m/s is removed from the
  outward velocity of every nest face point. The differences from the parent are what the finer cells change over
  the hill (tolerances 0.15 m/s and 0.06 m2/s2).
- Ranks: the hill patch on 2 and 4 ranks agrees with 1 rank to 1e-6 (relative or absolute; measured 5.0e-7 absolute,
  in p). It is not bit for bit: the line sweeps are local to a rank's boxes, and the iteration count differs by one.
- Restart: with a zone of 2 cells and a blend of 100 m, 10 iterations, then a restart to 20, equal 20 iterations
  straight, bit for bit, on the same boxes and ranks. The restarted run's log names the nest faces, the zone and the
  blended ground (moved by up to 0.12 m). A restart with another `bc.nest_relax_width` (3) or
  `terrain.nest_blend_width` (150), or from a parent of the same mesh, place and step with another flow, stops,
  naming `io.restart`.
- Refused with a nest face: lateral walls, buildings, and (with no inflow face) a Rayleigh layer on u, v or k and a
  side sponge.

## CTest
Double precision only.
| CTest | What it checks |
|---|---|
| `rt_nest_faces_clean`, `rt_nest_faces_flat_parent`, `rt_nest_faces_hill_parent` | Fixtures: remove old outputs; run the parents to steady on 2 ranks. |
| `rt_nest_faces_flat` | The flat patch on 3 ranks: checks 1 to 4, net flux at most 1e-10. |
| `rt_nest_faces_flat_zone` | The same with `bc.nest_relax_width=4` on 2 ranks. |
| `rt_nest_faces_flat_open` | One nest face, an outflow face and periodic y on 3 ranks: the net flux is reported as not removed; checks 2 to 4. |
| `rt_nest_faces_hill_np1`, `_np2`, `_np4` | The hill patch on 1, 2 and 4 ranks: checks 1 to 4, net flux between 1e-8 and 1e-3. |
| `rt_nest_faces_hill_np2_matches`, `_np4_matches` | `fcompare` of the final plotfiles against 1 rank, tolerance 1e-6. |
| `rt_nest_faces_hill_zone` | The hill patch with `bc.nest_relax_width=4` and `terrain.nest_blend_width=150`. |
| `rt_nest_faces_pair_a`, `_pair_b` | Fixtures: the hill parent for 30 iterations with u* 0.4 and 0.6 m/s. |
| `rt_nest_faces_restart_first`, `_continued`, `_straight`, `_matches` | Zone of 2 cells and blend of 100 m; 10 + 10 iterations against 20 straight: `fcompare` with zero tolerance; the restarted log names the zone and the blend. |
| `abort_nest_faces_restart_zone`, `_restart_blend`, `_restart_parent` | A restart with a changed zone or blend, or from the other parent of the pair, stops, naming `io.restart`. |
| `abort_nest_faces_no_parent`, `_no_guess`, `_sequencing`, `_widths`, `_blend_width`, `_keys_unused`, `_blend_unused`, `_buildings`, `_damping` | Input faults stop the run, naming the key. |
