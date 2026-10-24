# 2_STEM_IDENTIFICATION

Stage 2 of the pipeline. It reconstructs the stems of the multiple-stem trees
with the `dp_global` engine and merges the result, with the single-stem trees,
into the stem table read by stage 3 (`3_PREPARE_R_TABLES/`).

## Contents

| File | What it does |
| --- | --- |
| `1_main_cpp_chunk_bci.R` | Engine driver: sets the pins, estimates the biological parameters and runs the engine on the multiple-stem tags in parallel chunks |
| `2_merge_chunks_to_datatable.R` | Merges the chunk outputs, adds the single-stem tags, joins stems split by a moved point of measurement and checks the field tags |
| `run_chunk_bci.md` | Run command, settings and resume instructions of the driver |
| `tests/test_measurement_rejoin.R` | Tests of the measurement rejoin |
| `comparissons/` | Comparisons with the published stem tables of Condit et al. (2019); not under version control |

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
| `<BASE_OUT_DIR>/<run folder>/` | One Feather file per chunk (`stem_reconstruction_dp_global_rcpp_chunk_NNN.feather`), its `_done.txt` marker, `posteriors/` (one path file per tag), `pin_stop_at_fast_link.csv`, `run_log.txt` and `run_parameters_full.txt` |
| `DATA/<run_code>/merged_output.parquet`, `.rds` | Engine output of every chunk, merged |
| `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` | Final stem table: every record with its reconstructed stem (`ReconstructedStemID`) and how it was assigned (`ReconstructionMethod`); input of stage 3 |
| `DATA/PROCESSED/measurement_rejoin_audit.csv` | Stems joined by the measurement rejoin, with the sign of a moved point of measurement found on each tree (`evidence`) |
| `DATA/PROCESSED/measurement_rejoin_pairs.csv` | Observation pairs of those joins, applied to the posterior samples in stage 3 |

## What it does

### `1_main_cpp_chunk_bci.R`

1. Loads the census table and keeps the multiple-stem tags. Every growth
   form and species code named in `DP_FALLBACK_GROWTH_FORMS`,
   `NON_TAPER_CORRECTED_GROWTH_FORMS` and `PROB_SPECIES` must occur in the
   data.
2. Estimates the biological parameters (growth, mortality, recruitment) from
   the stems of censuses 7–9 (2010–2023), whose identity is known: one set per
   species with enough data whose diameters cover the size range, and pooled
   sets for all trees, all shrubs, palms with tree ferns, and strangler figs.
3. Sets the pins (`TrueStemID`), the records whose stem the database fixes.
   The database `StemID` is used only here and for unmeasured records left
   without a stem after the engine; the engine never reads it.
   - Step 1: every record with a `StemTag` is pinned to its `StemID`, and so
     is every other record from census 7 (2010) on, except a record without a
     DBH whose status is dead, stem dead or broken below. A field tag that the
     database recorded under two `StemID`s takes the first one.
   - Step 2: unmeasured records after a `StemID`'s last DBH take the pin of a
     pinned record of that phase.
   - Step 3a: a broken-below or resprout-coded record with a DBH is pinned to
     its `StemID` when that `StemID` has another record.
   - Step 3a.5: a `StemID` with no pin, alive measured records and a dead or
     broken-below record without a DBH gets its alive records pinned.
   - Step 3b: the pin of a `StemID`'s measured records fills its other records
     within each run of consecutive censuses; a census without a record of the
     `StemID` stops it.
   - Step 3c: a pinned link of a `StemID` before 2010 that grows faster than
     99.9% of the tagged links of its species loses its pin when another
     record of the tree fits the later record (details in `run_chunk_bci.md`).
4. Runs the engine (`dp_global`) on each tag: the exact DP, or the
   probabilistic matcher for palm clumps (*Oenocarpus mapora*, *Bactris
   major*), strangler figs and tags too large for the DP. Both engines use
   the same biological model, read the DBH of 1982 and 1985 below 55 mm as a
   5 mm class, and draw 200 posterior samples per tag.
5. Applies the BCI break convention of the censuses up to 2010: one
   broken-below record per tree, written on the main stem's record with the
   DBH of the largest living stem (or no DBH). The main stem of the census
   before ends there; the broken record is matched like any other measurement
   but cannot continue that main stem. The tree's other stems are not cut.
6. Keeps every pin after the engine, so each field-tagged stem stays one stem,
   its broken-below records without a DBH included. A stump, death or missing
   record without a DBH and without a field tag that comes before any life of
   its own `StemID` goes to the stem that broke.
7. Applies the post-engine helpers to each chunk (described in
   `dp_global/README.md`): pin-track rejoin, terminal backfills, terminal
   records to the stem that ended, broken-below invariants, renumbering of
   each tag's stems to 1..N, and the posterior path files.
8. Writes one Feather file per chunk. Chunks with a `_done.txt` marker are
   skipped when a run is resumed.

### `2_merge_chunks_to_datatable.R`

1. Reads the chunk files of `home_dir/run_code` (the run must have
   `run_finished.txt`) and writes them, merged, to `DATA/<run_code>/`.
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
   the pair must be clearly the best candidate (`apply_measurement_rejoin()`;
   rules in `dp_global/R/measurement_rejoin.R`). The joined stem keeps one of
   its two IDs (`ReconstructionMethod = "measurement_rejoin"`).
4. Restores the recorded `ExactDate`. The driver fills a missing date with
   the mean date of the tag in its census, because the engine needs a date to
   compute growth.
5. Checks that every pin of a record relabelled by the broken-below rules is a
   database `StemID` of its tree (a pin that is not takes the record's own
   `StemID`).
6. Checks the field tags and writes the final table, the rejoin audit and the
   rejoin pairs.

## Checks

A failed check prints `❌ CHECK FAILED`, raises a warning and stops. The merge
checks that:

- the engine output holds every multiple-stem record of the stage-1 table,
  and every stage-1 column of the final table equals the stage-1 table
  (`growth_form` aside, the engine's label);
- stage 2 changed no recorded DBH or date;
- no stem has two records in one census and no stem holds two pins after the
  measurement rejoin;
- every record with a field tag (2010 on) is pinned, every field tag is one
  reconstructed stem, and no reconstructed stem holds two field tags.

The log also counts the records before 2010 whose copied `StemTag` lost its
pin in Step 3c.

## Tests

```r
testthat::test_file("BCI_stem_reconstruction/2_STEM_IDENTIFICATION/tests/test_measurement_rejoin.R")
```

The tests check each rule of the measurement rejoin on constructed trees.

## Requirements

R packages: `data.table`, `parallel`, `arrow` (Feather and Parquet), and those
loaded by `dp_global/R/dp_global_main.R` (`Rcpp`, `igraph`, `MASS`,
`lpSolve`); `testthat` for the tests. A C++ compiler for the engine.
