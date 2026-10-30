# Abl_Setup

## Purpose
Runs behind the named set-ups `abl.drive` and `abl.theta` (`Source/AblSetup.H`, documented in `docs/howto.rst`).
A name gives a few keys a default and checks the deck, and adds no physics. So a named run must match the same deck
with the keys written out, bit for bit, and a deck that contradicts its name must stop.

## Files
| File | Description |
|---|---|
| `inputs_abl_column` | A flat periodic column with a geostrophic wind; every key except `forcing.coriolis` and `forcing.geostrophic`. |

The other runs use `Tests/Canonical_RANS/Neutral_Column/inputs_steady_column` and its check script.

## Data
None.

## How to run
```bash
ctest --test-dir <build> -R "abl_"
```

## Expected Results
- `rt_abl_named` and `rt_abl_explicit` give identical plotfiles (`rt_abl_named_matches`, zero tolerance).
- `rt_abl_pressure_driven`: the neutral column under `abl.drive = pressure_driven` reaches the closed form
  $u_*^2 = F H$ (`check_steady_column.py`).
- The `abort_abl_*` decks stop, naming the input.

## CTest
`rt_abl_named`, `rt_abl_explicit`, `rt_abl_named_matches`, `rt_abl_pressure_driven`, `abort_abl_names`,
`abort_abl_contradiction`, `abort_abl_deck_value`, `abort_abl_no_inflow`, `abort_abl_couette_stratified`, `abort_abl_capped_no_sink`,
`abort_abl_no_stage`. The tables of what each name needs are in the unit test `AblSetup.*`.
