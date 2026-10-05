# 3_PREPARE_R_TABLES

Stage 3 of the BCI stem-reconstruction pipeline. This stage consumes the
reconstructed stem identities from `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/` and generates the
ForestGEO-format census tables plus QC exports used for downstream analysis.

## Scripts (run in order)

- `1_prepare_posteriors_BCI.R` — Consolidates `_paths.feather` posterior files
  from a completed stage 2 run into `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`,
  applying to every sample the measurement-discontinuity joins of the merge
  step (`DATA/PROCESSED/measurement_rejoin_pairs.csv`).
- `2_create_R_tables_BCI.R` — Builds the final census tables and species table.
  It assigns each stem a corrected status (`Rstatus`) in every census,
  gives every tree one location, imputes the dates that are needed, and
  exports ForestGEO-format `.Rdata` (and supporting `.csv`) tables.
- `rstatus_functions.R` — The `Rstatus` and `dbh` rules as small functions,
  sourced by `2_create_R_tables_BCI.R` and by the tests.

## Rstatus and dbh rules

`Rstatus` is the corrected status; `DFstatus` keeps the raw field status.

| Rstatus | Meaning |
| --- | --- |
| `P` | before the stem's first record |
| `A` | alive: from the first record to the last evidence of life |
| `G` | stem dead, its tree alive (some stem of the tree is `A` then or later) |
| `D` | stem dead and its whole tree dead |

- Evidence of life: an `alive` record, or a DBH on a `broken below`,
  `missing` or status-less record.
- Dead record: `dead` / `stem dead` (with or without DBH), and `broken below`
  without DBH.
- Missed censuses and false deaths between a stem's first record and its last
  evidence of life are `A` (alive later ⇒ never dead). A dead first record of
  a stem alive later is `A` too.
- After the last evidence of life a stem is dead from the first census
  without one. A stem never alive is `P` until its first record, then `G`/`D`.
  A stem never recorded (only `missing` or empty rows) is `P` in every census.
- `dbh` is exported exactly as recorded: no DBH is removed or imputed. A `G`
  or `D` cell carries a DBH only when it was recorded on a `dead` / `stem dead`
  record of a stem never alive later (listed in
  `CHECKS/dbh_on_dead_records.csv`); a `P` cell never has one; an `A` without
  DBH is a missed measurement.

## Other exported columns

- `stemID`: the reconstructed stem number within its tree (1, 2, ...), so
  `treeID` + `stemID` identify a stem. Internally the script keys stems on
  `TreeID_ReconstructedStemID`; a stem without a DP identity (records with no
  status and no DBH) gets `stemID` NA.
- `ExactDate`: every recorded date is kept. A missing date is imputed (modal
  field date of the tree, else the quadrat, else the census) only on `A` rows,
  measured rows and the first `G`/`D` row of a stem; `P` rows and later `G`/`D`
  rows have a date only when one was recorded. Stage-4 scripts fill the dates
  they need themselves.
- `gx`, `gy`, `quadrat`: one location per tree, the same for every stem and
  census: the position most censuses agree on (one vote per census; ties go to
  the most recent census), with the quadrat derived from it.

## Tests

`tests/test_Rstatus.R` (testthat) checks the rules with hand-written
fixtures, every combination of raw states for small trees against a simple
reference implementation (`tests/reference_rstatus.R`), the validity of every
`Rstatus` sequence (`tests/rstatus_validity.R`), and the exported R tables
cell by cell (`Rstatus`, `dbh`, `DFstatus`, `ExactDate`, location and
`stemID`). Run from the project root:

```r
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R")
```

`RSTATUS_TEST_LEVEL=quick` runs the smaller enumerations only (about a
minute; `full`, the default, takes about 15 minutes). `RSTATUS_RTABLES_DIR`
and `RSTATUS_STAGE2_FILE` point the real-data tests at other tables (default:
`DATA/RTABLES` and `DATA/PROCESSED`).

## Outputs

- `BCI_stem_reconstruction/DATA/RTABLES/<site>.stemN.Rdata`
- `BCI_stem_reconstruction/DATA/RTABLES/<site>.spptable.rdata`
- `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`
- QC exports in `BCI_stem_reconstruction/DATA/CHECKS/`, among them
  `dbh_on_dead_records.csv`, `location_conflicts.csv`,
  `subset_stems_never_alive.csv` and `subset_stems_never_recorded.csv`

## Notes

- `1_prepare_posteriors_BCI.R` must be run before `2_create_R_tables_BCI.R`.
