# Infrastructure

## Purpose
Cases that check the application around the solver: a run and its restart, input decks that must stop, reference
archives, case set-up from YAML, wind libraries and the Kestrel scripts.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Smoke_Restart` | `smoke_*`, `abort_restart_*` | A 2-rank run, then a 3-rank restart on another box layout; the two plotfiles identical. Two decks that must stop on a bad restart. |
| `Input_Validation` | `abort_*`, `warn_w_vnn_cap_*` | Decks that must stop, with the input named in the message. |
| `Abl_Setup` | `rt_abl_*`, `abort_abl_*` | The named set-ups `abl.drive` and `abl.theta`: a named run equals the written-out deck; contradictions stop. |
| `Building_Field` | `buildings_*`, `abort_buildings_*` | The building geometry from a height map: column volumes, wall mask, restart and refusals. |
| `Canopy_Field` | `canopy_*`, `abort_canopy_*` | The forest canopy's drag field: column integrals against a closed form, restart and refusals. |
| `Reference_Archive` | `refavg_*`, `abort_refavg_*` | `gust_refavg` end to end on copies of the smoke plotfile. |
| `Case_Setup` | `rt_case_setup_*` | YAML case descriptions written by `tools/make_gust_case.py`, run through `run.sh` and checked against closed forms. |
| `Library_Hill` | `rt_library_hill` | A hill sweep run through `run.sh`; wind libraries built from it are blended between directions and checked against the run in between. |
| `Library_Ensemble` | `rt_library_ensemble` | Ensemble statistics and WindNinja-format output of a wind library, against the same statistics from the runs' plotfiles. |
| `Kestrel_Dryrun` | `kestrel_*` | A dry run of the Kestrel submit and archive scripts with stub `module`/`sbatch`/`srun`. |
