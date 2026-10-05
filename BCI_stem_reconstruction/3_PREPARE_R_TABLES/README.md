# 3_PREPARE_R_TABLES

Stage 3 of the BCI stem-reconstruction pipeline. This stage consumes the
reconstructed stem identities from `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/` and generates the
ForestGEO-format census tables plus QC exports used for downstream analysis.

## Scripts

- `1_prepare_posteriors_BCI.R` — Consolidates `_paths.feather` posterior files
  from a completed stage 2 run into `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`,
  applying to every sample the measurement-discontinuity joins of the merge
  step (`DATA/PROCESSED/measurement_rejoin_pairs.csv`).
- `2_create_R_tables_BCI.R` — Builds the final census tables from
  `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds`.
  It assigns each stem a corrected status (`Rstatus`) in every census,
  gives every tree one location, fills the `date` column, and
  exports ForestGEO-format `.Rdata` (and supporting `.csv`) tables.
- `rstatus_functions.R` — The `Rstatus` and `dbh` rules as small functions,
  sourced by `2_create_R_tables_BCI.R` and by the tests.

Both scripts read the output of the stage-2 merge
(`2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R`), so run it first. They
do not depend on each other: `2_create_R_tables_BCI.R` does not read the
posterior file, which is used by
`4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R`.

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

- `stemID`: the reconstructed stem number within its tree, so `treeID` +
  `stemID` identify a stem. The engine numbers the stems of a tree 1, 2, ...
  in order of first appearance; in a tree changed by the measurement rejoin of
  the stage-2 merge, the joined stem keeps one of its two numbers, so the
  tree skips one number per join. Internally the script keys stems on
  `TreeID_ReconstructedStemID`; a stem without a DP identity (records with no
  status and no DBH) gets `stemID` NA.
- `ExactDate`: the field date exactly as recorded (`NA` where there is no
  record); never imputed, like `dbh` and `DFstatus`.
- `date`: the ForestGEO R-table date, in days since 1960-01-01 (the unit of
  the ForestGEO database's `Date` field and of Condit's tables;
  `as.Date(date, origin = "1960-01-01")` gives the calendar date). Every row
  has one:
  1. the recorded date (`ExactDate`) where there is one;
  2. otherwise the most common recorded date of the same tree in that census;
  3. if the tree has none, the most common recorded date of the same quadrat
     in that census;
  4. if the quadrat has none, the most common recorded date of the census.

  Only recorded dates vote, and ties go to the earliest date. The most common
  date (the mode) is used rather than a mean or median because it is always a
  day on which the field crew was recording. Condit's Dryad tables fill `date`
  with the mean recorded date of the quadrat instead: the two agree within
  about a day in 1985–2015, but in 1982 the records of one quadrat span months,
  so the tree's own date is closer. Census intervals (growth, recruitment,
  mortality) are computed from `date`.
- `gx`, `gy`, `quadrat`: one location per tree, the same for every stem and
  census. Each census casts one vote, the most common (`PX`, `PY`) pair among
  the tree's stems in that census (ties go to the smallest x, then y); the
  tree takes the pair with the most votes (ties go to the most recent census).
  `quadrat` is derived from that pair (20 m quadrats named `XXYY`); a tree
  without coordinates keeps `NA` coordinates and its most common raw quadrat.

## Tests

`tests/test_Rstatus.R` (testthat) checks the rules with hand-written
fixtures, every combination of raw states for small trees against a simple
reference implementation (`tests/reference_rstatus.R`), the validity of every
`Rstatus` sequence (`tests/rstatus_validity.R`), and the exported R tables
cell by cell (`Rstatus`, `dbh`, `DFstatus`, `ExactDate`, `date`, location and
`stemID`). Run from the project root:

```r
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R")
```

`RSTATUS_TEST_LEVEL=quick` runs the smaller enumerations only (a few minutes
with the real-data tests; `full`, the default, takes about 10 minutes on 16
cores). `RSTATUS_RTABLES_DIR`
and `RSTATUS_STAGE2_FILE` point the real-data tests at other tables (default:
`DATA/RTABLES` and `DATA/PROCESSED`); if the tables or the stage-2 file are
missing, the real-data tests are skipped. `RSTATUS_CORES` sets the cores of
the reference implementation (default: all but two).

## Outputs

- `BCI_stem_reconstruction/DATA/RTABLES/<site>.stemN.Rdata` (and `<site>.stemN.csv`)
- `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`
- QC exports in `BCI_stem_reconstruction/DATA/CHECKS/`, among them
  `dbh_on_dead_records.csv`, `location_conflicts.csv`,
  `subset_stems_never_alive.csv`, `subset_stems_never_recorded.csv` and
  `duplicate_measurements.csv` (likely duplicates: two stems of one tree
  measured on the same date at heights at least 0.5 m apart, taper-corrected
  DBHs within 20 %; the trunk measured at the old and the new height when its
  point of measurement was raised)

The species table is not written by this stage (Section 15 of
`2_create_R_tables_BCI.R` is commented out). The stage-4 scripts read
`DATA/RTABLES/bci.spptable.rdata`, so copy the ForestGEO-format table there
(for example `data_paper_and_repo_publication/RTABLES/bci.spptable.rdata`).
