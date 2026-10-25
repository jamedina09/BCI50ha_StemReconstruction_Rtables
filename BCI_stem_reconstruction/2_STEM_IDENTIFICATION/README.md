# 2_STEM_IDENTIFICATION

Stage 2 of the pipeline. It reconstructs the stems of the multiple-stem trees
with the `dp_global` engine and merges the result, with the single-stem trees,
into the stem table read by stage 3 (`3_PREPARE_R_TABLES/`).

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [What it does](#what-it-does)
  - [`1_main_cpp_chunk_bci.R`](#1_main_cpp_chunk_bcir)
  - [`2_merge_chunks_to_datatable.R`](#2_merge_chunks_to_datatabler)
- [Checks](#checks)
- [Settings](#settings)
- [Tests](#tests)
- [Data flow](#data-flow)
- [Requirements](#requirements)

## Contents

| File | What it does |
| --- | --- |
| `1_main_cpp_chunk_bci.R` | Engine driver: sets the pins, estimates the biological parameters and runs the engine on the multiple-stem tags in parallel chunks |
| `2_merge_chunks_to_datatable.R` | Merges the chunk outputs, adds the single-stem tags, joins stems split by a moved point of measurement and checks the field tags |
| `run_chunk_bci.md` | Run command, settings and resume instructions of the driver |
| `run_records/` | The record of every engine run, written by the driver: settings with their defaults, rerun command, copy of the driver, provenance, and `runs_index.csv` (one line per run); described in `run_records/README.md` |
| `tests/test_measurement_rejoin.R` | Tests of the measurement rejoin |

## How to run

From the project root:

1. Run the driver with the command in `run_chunk_bci.md`. It writes to a new
   run folder under `BASE_OUT_DIR` (outside the repository).
2. Set `home_dir` (the `BASE_OUT_DIR` of the run) and `run_code` (the run
   folder name) at the top of `2_merge_chunks_to_datatable.R`, then run it:

```sh
Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R
```

## Inputs

| File | Content |
| --- | --- |
| `DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds` | Census table of stage 1 |

## Outputs

| File | Content |
| --- | --- |
| `<BASE_OUT_DIR>/<run folder>/` | One Feather file per chunk (`stem_reconstruction_dp_global_rcpp_chunk_NNN.feather`), its `_done.txt` marker, `posteriors/` (one path file per tag), `pin_stop_at_fast_link.csv`, `run_started.txt`, `run_finished.txt`, `run_log.txt`, `run_parameters_full.txt` (command line, input file with its modification time and size, the key settings and the parameter structure; read by the merge checks) and the files of the run record |
| `run_records/<run folder>/` | Run record: `run_settings.txt` (method summary and every setting with its default and whether it was changed), `rerun_command.sh`, a copy of the driver, `run_summary.txt` (time, git commit, R and package versions, machine, input file with its MD5, chunks, status) and `uncommitted_changes.patch` when the code had uncommitted changes |
| `run_records/runs_index.csv` | One line per run: status, times, git commit, input MD5, test run or not, changed and key settings |
| `DATA/<run_code>/merged_output.parquet`, `.rds` | Engine output of every chunk, merged |
| `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` | Final stem table: every record with its reconstructed stem (`ReconstructedStemID`) and how it was assigned (`ReconstructionMethod`); input of stage 3 |
| `DATA/PROCESSED/measurement_rejoin_audit.csv` | Stems joined by the measurement rejoin, with the sign of a moved point of measurement found on each tree (`evidence`) |
| `DATA/PROCESSED/measurement_rejoin_pairs.csv` | Observation pairs of those joins, applied to the posterior samples in stage 3 |
| `DATA/PROCESSED/engine_run.txt` | The engine run behind the final table (`run_code`, `home_dir`, finish and merge times); stages 3 and 4 check that they use the same run |
| `DATA/CHECKS/check_merge_bci.txt` | Every check of the merge (✓ passed, ❌ failed) |

## What it does

### `1_main_cpp_chunk_bci.R`

1. Loads the census table (DBH converted from mm to cm), assigns the growth
   forms (tree, shrub, palm, strangler, fern) and keeps the multiple-stem tags
   for the engine (single-stem tags are added by the merge). Every growth form
   and species code named in `DP_FALLBACK_GROWTH_FORMS`,
   `NON_TAPER_CORRECTED_GROWTH_FORMS` and `PROB_SPECIES` must occur in the
   data; a value that does not (e.g. `palms` for the label `palm`) stops the
   run. Censuses 1 and 2 must show the 5 mm DBH classes (`DBH_ROUND_CENSUSES`).
2. Estimates the biological parameters (growth, mortality, recruitment) from
   the stems of all tags in censuses 7–9 (2010–2023), whose identity is known
   (`dp_global/README.md`, "Parameter estimation"): one set per species with
   enough data (trees, shrubs, palms), and pooled sets for the other species
   of each growth form (all trees, all shrubs, palms with tree ferns, strangler
   figs). A species gets its own set when its diameters fill the size range,
   either from the smallest to the largest or over their middle 95%
   (`COVERAGE_RULE = "union"`), so one out-of-range diameter does not send a
   well-sampled species to a pooled set. The growth SD is fitted as a constant
   when it falls with size, and the recruitment rate is new stems per
   established tree per year (`RECRUIT_RATE_UNIT = "tree"`).
3. Sets the pins (`TrueStemID`), the records whose stem the database fixes.
   The database `StemID` is used only here and for unmeasured records left
   without a stem after the engine (a death, stump or missing record joins the
   stem of the previous record with its `StemID`; a stem never measured keeps
   its `StemID`); the engine never reads it.
   - Step 1: every record with a `StemTag` is pinned to its `StemID`, and so
     is every other record from census 7 (2010) on, except a record without a
     DBH whose status is dead, stem dead or broken below. A field tag that the
     database recorded under two `StemID`s takes the first one.
   - Step 2: unmeasured records after a `StemID`'s last DBH take the pin of a
     pinned record of that phase.
   - Step 3a: a broken-below or resprout-coded record with a DBH is pinned to
     its `StemID` when that `StemID` has another record
     (`PIN_BB_REQUIRE_OTHER_RECORD`). Before 2010 the database gives a broken
     record a new `StemID`; when it has no other record, the pin carries no
     identity but would forbid linking the record to the pinned stem it grows
     into.
   - Step 3a.5: a `StemID` with no pin, alive measured records and a dead or
     broken-below record without a DBH gets its alive records pinned.
   - Step 3b: the pin of a `StemID`'s measured records fills its other records
     within each run of consecutive censuses; a census without a record of the
     `StemID` stops it (`PIN_STOP_AT_GAP`): the database sometimes attached a
     2010 stem to an old `StemID` with no record in between.
   - Step 3c (`PIN_STOP_AT_FAST_LINK`): a pinned link of a `StemID` between two
     consecutive censuses before 2010 that grows faster than 99.9% of the
     tagged links of its species is contradicted by the measurements when
     another record of the tree in the earlier census (unpinned, its `StemID`
     without a measured record in the later census) grows into the later record
     within [max shrink, 99% of the tagged links]. The records of that `StemID`
     up to the earlier census lose their pin and the engine places them. This
     is the case where the database gave the `StemID` of a small 2005 record to
     the 2010 stem that continues another 2005 record. The links freed are
     written to `pin_stop_at_fast_link.csv`.
4. Runs the engine (`dp_global`) on each (tag, species) group, in parallel
   chunks: the exact DP, or the probabilistic matcher for palm clumps
   (*Oenocarpus mapora*, *Bactris major*), strangler figs and tags too large
   for the DP (`PROB_SPECIES`, `DP_MAX_STATES`). Both engines use the same
   biological model, read the DBH of 1982 and 1985 below 55 mm as a 5 mm class
   (`DBH_ROUND_CENSUSES`), and draw 200 posterior samples per tag. The matcher
   lets every stem continue, die or be recruited in every census pair
   (`PROB_BIRTH_DEATH`) and exports its most representative sample
   (`PROB_EXPORT_RULE = "mea"`; `"fewest_me"` first takes the samples with the
   fewest links that shrink beyond measurement error). The DP posterior counts
   each reconstruction once (`DP_POSTERIOR_LABEL_FREE`): tracks with no stem
   in 2023 and no pinned measurement are interchangeable labels.
5. Applies the BCI break convention of the censuses up to 2010: one
   broken-below record per tree, written on the main stem's record with the
   DBH of the largest living stem (or no DBH). The main stem of the census
   before ends there; the broken record is matched like any other measurement
   but cannot continue that main stem (`BB_DBH_RULE = "main_ends"`). The
   tree's other stems are not cut at the break (`DP_ALLOW_SEGMENT_SPLIT =
   FALSE`).
6. Keeps every pin after the engine (`HARD_PINS`), so each field-tagged stem
   stays one stem. Field tags are the `StemTag`s from 2010
   (`FIELD_TAG_FIRST_CENSUS = 7`), the first census that tagged every stem; a
   broken-below record without a DBH that carries a field tag is that stem with
   no measurable DBH at 1.3 m, so it stays with its tag (e.g. a 2010 stump whose
   tag is measured in 2015). Exception, so that a stump stays the end record of
   the stem that broke: a stump, death or missing record without a DBH that
   comes before any life of its own `StemID` goes to the stem that broke when
   it carries no field tag, and so do all records of a `StemID` that never
   lives. (Up to 2010 the database stores a break under the new `StemID` of the
   resprout measured later, and shows that `StemID`'s later tag on its earlier
   records; those tags are not field tags.)
7. Applies the post-engine helpers to each chunk (described in
   `dp_global/README.md`): pin-track rejoin, terminal backfills, terminal
   records to the stem that ended, broken-below invariants, renumbering of
   each tag's stems to 1..N, and the posterior path files.
8. Writes one Feather file per chunk and the posterior path files to
   `BASE_OUT_DIR/<run folder>/` (the folder name starts with the run's time
   stamp). Chunks with a `_done.txt` marker are skipped when a run is resumed.
9. Keeps the run record in `run_records/<run folder>/` and in the run folder:
   written at the start, completed at the finish (or at an error), with the
   run's line in `run_records/runs_index.csv` (`run_records/README.md`).

A blind-test mode checks the engines against the tagged censuses
(`BLIND_MIN_CENSUS`, `BIO_ESTIMATION_START_CENSUS`; diagnostics only): the
engine sees only the censuses from `BLIND_MIN_CENSUS` on and every pin before
the anchor is removed. All settings are listed in `run_chunk_bci.md`.

### `2_merge_chunks_to_datatable.R`

1. Reads the chunk files of `home_dir/run_code` (the run must have
   `run_finished.txt`), converts them in batches into temporary Parquet parts
   (to cap memory), merges them and writes `merged_output.parquet` and `.rds`
   to `DATA/<run_code>/`.
2. Compares the merged table with the stage-1 table, restores the recorded
   DBH (mm) and adds the single-stem tags: their one stem gets
   `ReconstructedStemID = 1` (`ReconstructionMethod =
   "single_stem_tag_no_reconstructed"`). Records of tags without any DBH take
   their `StemID`.
3. Joins a stem that the engine ended only because the point of measurement
   of its trunk moved (buttresses, deformities, strangler figs), so that the
   next measurement fell outside the growth bounds. The later stem must be an
   impossible recruit (above the recruit size limit), the size change must be
   the one a moved point of measurement gives, the tree must show a sign of it
   (a height of measurement away from 1.3 m, a buttress code B, a strangler
   growth form, or in 1982–1985 a species that usually has buttresses), and
   the pair must be clearly the best candidate (same height of measurement,
   then the closest size). The later stem may also start two censuses later
   when the trunk was not measured in between and its recorded height of
   measurement rose. Recording errors and trunks only missed in one census are
   not joined, and database `StemID`s are not used (`apply_measurement_rejoin()`;
   rules and their justification in `dp_global/R/measurement_rejoin.R`). The
   joined stem keeps one of its two IDs, so the tree's IDs skip one number per
   join (`ReconstructionMethod = "measurement_rejoin"`).
4. Restores the recorded `ExactDate`. The driver fills a missing date with
   the mean date of the tag in its census (else of the census), because the
   engine and the rejoin need a date to compute growth; the dates are restored
   after the rejoin, so the rejoin sees the same dates as the engine. Stage 3
   gives every row a date in its own column `date`.
5. Checks that every pin of a record relabelled by the broken-below rules is a
   database `StemID` of its tree (a pin that is not takes the record's own
   `StemID`).
6. Checks the field tags and writes the final table, the rejoin audit and the
   rejoin pairs.

## Checks

A failed check prints `❌ CHECK FAILED`, raises a warning, writes to the log
and stops.

The driver stops when the stage-1 species table is newer than the census table
(stage 1 run out of order: rerun `2_prepare_viewfulltable.R`), when a
setting names a growth form or species code that is not in the data, and when
a resumed run has another input file (MD5), chunk size or setting outside run
control than at its first start (read from its run record).

The merge writes every check to `DATA/CHECKS/check_merge_bci.txt`. Before it
reads anything it checks that the run is the right one and complete:

- the stage-1 table exists, and the run folder `home_dir/run_code` exists and
  has `run_finished.txt`;
- every chunk 1..N has its Feather file and its `_done.txt` marker, and no
  chunk failed without being rerun;
- the run is a full run: no `TAG_FILTER_FILE`, `TAG_FILTER_N`,
  `BLIND_MIN_CENSUS`, `BIO_ESTIMATION_START_CENSUS`, `DP_CHUNK_START` or
  `DP_CHUNK_END`;
- the run read the current stage-1 table (same modification time and size as
  recorded in `run_parameters_full.txt`; for a run without that record, the run
  started after the table was written). After rerunning stage 1, the engine
  must be run again;
- the rejoin bounds of the merge equal the growth bounds and recruit limit of
  the run.

Then it checks that:

- the engine output holds every multiple-stem record of the stage-1 table once,
  and every stage-1 column of the final table equals the stage-1 table
  (`growth_form` aside, the engine's label);
- stage 2 changed no recorded DBH or date;
- no stem has two records in one census and no stem holds two pins after the
  measurement rejoin;
- every record with a field tag (2010 on) is pinned, every field tag is one
  reconstructed stem, and no reconstructed stem holds two field tags.

The log also counts the records before 2010 whose copied `StemTag` lost its
pin in Step 3c.

## Settings

The driver settings (identity rules, engine limits, posterior sampling,
chunking, cores, output folder) are listed with their defaults in
`run_chunk_bci.md`. The merge has `home_dir` (the `BASE_OUT_DIR` of the run)
and `run_code` (the run folder), and the rejoin bounds `REJOIN_MAX_SHRINK`,
`REJOIN_MAX_GROWTH` (cm/yr) and `REJOIN_RECRUIT_MAX_MM`, which must equal the
driver's `MAX_SHRINK_FIXED`, `MAX_GROWTH_FIXED` and `RECRUIT_MAX_FIXED`.

## Tests

```r
testthat::test_file("BCI_stem_reconstruction/2_STEM_IDENTIFICATION/tests/test_measurement_rejoin.R")
```

The tests check each rule of the measurement rejoin on constructed trees.

## Data flow

```text
DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
        │
        ▼
1_main_cpp_chunk_bci.R
        ├──▶ BASE_OUT_DIR/<run folder>/  (chunk Feather files, posteriors/, run record)
        └──▶ run_records/<run folder>/, run_records/runs_index.csv

2_merge_chunks_to_datatable.R
        ├──▶ DATA/<run_code>/merged_output.{parquet,rds}
        └──▶ DATA/PROCESSED/
               complete_dataset_final_with_reconstructed_stemids.rds
               measurement_rejoin_audit.csv, measurement_rejoin_pairs.csv
```

The comparisons of the outputs with the published stem tables of Condit et al.
(2019) are in `COMPARE_DATASETS/` at the project root, one folder per stage
(not under version control); the stem trajectories, link by link, are in
`COMPARE_DATASETS/2_STEM_IDENTIFICATION/`.

## Requirements

R packages: `data.table`, `parallel`, `arrow` (Feather and Parquet), and those
loaded by `dp_global/R/dp_global_main.R` (`Rcpp`, `igraph`, `MASS`,
`lpSolve`); `testthat` for the tests. A C++ compiler for the engine.
