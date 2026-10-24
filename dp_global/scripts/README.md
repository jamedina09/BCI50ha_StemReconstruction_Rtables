# dp_global/scripts

Driver scripts of the `dp_global` engine, and a basal-area uncertainty script
for their outputs. The BCI production run uses its own driver,
`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`.

## Contents

| File | What it does |
| --- | --- |
| `main_cpp.R` | Driver for one tag or a small data set; can also run sensitivity sweeps and a realism report |
| `main_cpp_chunk.R` | Driver for large data sets: processes the tags in chunks, writes each chunk to disk and resumes interrupted runs |
| `main_cpp_bci.R` | Debug driver for single tags of the BCI sample file |
| `basal_area_uncertainty.R` | Basal-area stocks and changes per tag, with the uncertainty from the posterior path files of a run |

## How to run

From the project root:

```sh
Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WHICH_TAG=20
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --PLOT_PDF_ONE_TAG_ONLY=FALSE
Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=<file>
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
Rscript dp_global/scripts/basal_area_uncertainty.R --RUN_DIR=dp_global/output/<run folder>
```

Any top-level variable of a driver can be set as `--NAME=value` (case does not
matter; `-` and `_` are equivalent; write text values such as `WHICH_TAG`
with `_`). An unknown name gives a warning and is ignored. `--help` lists the
main settings. `source("dp_global/scripts/main_cpp.R")` defines the settings
and functions without starting a run (`run_main()` starts it).

## Inputs

| Script | Input |
| --- | --- |
| `main_cpp.R`, `main_cpp_chunk.R` | `INPUT_FILE`, a CSV with at least `Tag`, `CensusID`, `DBH` (cm), `ExactDate` and `TrueStemID`; `OriginalStemID` (or `StemID`), `Status`, `ListOfTSM`, `growth_form` and `HOM` are used when present |
| `main_cpp_bci.R` | `data_simulation/sample_data_BCI/multistem_tags.rds`: BCI multiple-stem tags with `Species`, `Tag`, `StemTag`, `OriginalStemID`, `CensusID`, `ExactDate`, `DBH` (cm), `ListOfTSM`, `HOM`, `Status` and the `Bio_*` columns |
| `basal_area_uncertainty.R` | `<RUN_DIR>/stem_reconstruction_dp_global_rcpp.csv` and the CSV path files `<RUN_DIR>/posteriors/*_paths.csv` |

## Outputs

Each driver run writes to a new folder under `dp_global/output/` (name built
from the time stamp and the settings; `OUT_DIR_OVERRIDE` sets the folder):

| File | Content |
| --- | --- |
| `run_started.txt`, `run_finished.txt` | Start and end time stamps |
| `run_parameters_full.txt` | Settings, command line and parameter sets |
| `run_log.txt` | Log of the run |
| `stem_reconstruction_dp_global_rcpp.csv` | Reconstruction (`WRITE_DP_CSV`; the chunked driver appends each chunk) |
| `stem_reconstruction_dp_global_rcpp.rds`, `.feather` | The same in binary form (`main_cpp.R`; `WRITE_DP_RDS`, `WRITE_DP_FEATHER`) |
| `stem_reconstruction_dp_global_rcpp_chunk_NNN.rds`, `.feather` | One file per chunk (`main_cpp_chunk.R`), with a `_done.txt` marker, or `_failed.txt` when the chunk raised an error |
| `stem_reconstruction_dp_global_rcpp.pdf`, `_chunk_NNN.pdf` | Reconstructed trajectories per tag (`WRITE_DP_PDF`) |
| `posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<csv\|rds\|feather>` | Posterior path file of each tag (format in `dp_global/README.md`) |
| `bio_pars_report.pdf` | Plots of the estimated parameters (`main_cpp.R`) |
| realism and sensitivity files | With `RUN_REALISM_REPORT`, `SENSITIVITY_MODE` or `RUN_K_SWEEP_DEMO` (`main_cpp.R`) |

`basal_area_uncertainty.R` writes to `<RUN_DIR>/`:

| File | Content |
| --- | --- |
| `basal_area_tag_census.csv` | Per tag and census: total basal area (m²), number of stems, year |
| `basal_area_tag_change.csv` | Per tag and census pair: interval, growth, loss, gain and total change of the exported reconstruction, numbers of survivors, deaths and recruits, and the posterior mean, SD and 95% interval of each component |
| `basal_area_figures.pdf` | Summary pages, per-tag panels (trajectory, stem count, change components, stem demography), uncertainty histograms and posterior densities |

## What it does

### Drivers

1. Load the input, set the parameter set of each row (`species`; one set for
   all rows with `FORCE_ONE_SPECIES_PARAMETERS = TRUE`, label `"all"`) and
   estimate the biological parameters per set (`estimate_bio_pars()`). The
   BCI debug driver estimates none: its input carries the `Bio_*` columns.
2. Run the engine on each (`Tag`, `species`) group through
   `run_dp_one_group()`: the DP, or the probabilistic matcher on fallback.
   A DP error sends the group to the matcher (`DP_FallbackReason =
   "error:<message>"`).
3. Apply the pin sweep and the post-engine helpers, then write the posterior
   path files (order and rules in `dp_global/README.md`, *Post-engine
   helpers*).
4. Write the outputs.

`main_cpp.R` runs one tag (`RUN_ALL_TAGS = FALSE`, `WHICH_TAG`) or all tags in
parallel, and holds the full result in memory. `main_cpp_chunk.R` always runs
every group: it splits the groups into chunks of `DP_CHUNK_SIZE`, runs each
chunk in parallel, applies the helpers to the chunk, writes it and frees it
before the next one. A chunk with a `_done.txt` marker is skipped when
`DP_CHUNK_RESUME = TRUE`. `merge_chunks_to_csv(<run folder>)` (defined in
`main_cpp_chunk.R`) streams the chunk files into one CSV.

`main_cpp_bci.R` sources `main_cpp.R` for its functions and sets the pins
before the engine from the BCI database conventions:

- Step 1: rows with a `StemTag`, and rows from census 7 on (except rows
  without a DBH whose status is dead, stem dead or broken below), are pinned
  to their `OriginalStemID`.
- Step 2: unmeasured rows after an `OriginalStemID`'s last DBH take the pin of
  a pinned row of that phase.
- Step 3a: rows with a DBH that are broken below or carry a resprout code are
  pinned to their `OriginalStemID`.
- Step 3a.5: an `OriginalStemID` with no pin, alive measured rows and a dead
  or broken-below row without a DBH gets its alive rows pinned.
- Step 3b: the pin of an `OriginalStemID`'s measured rows fills its other
  rows.

It has none of the later pin rules of the production driver (Step 3a only with
another record, no pin across a census gap, Step 3c, hard pins, field tags).

### `basal_area_uncertainty.R`

1. Reads the reconstruction and the posterior path files of the run; without
   path files it gives the exported reconstruction only. DBH is in cm.
2. Sums the basal area of the measured stems per tag and census,
   π/4 × (DBH/100)² in m², which does not depend on stem identity.
3. Splits each change between censuses into growth of surviving stems, loss by
   mortality and gain by recruitment, which do, for the exported
   reconstruction and for every posterior path. The paths hold the
   observations up to the anchor; the observations after it (and before the
   paths' first census) are spliced into every path on the path stem of their
   exported stem-mate at the window edge, so every census pair is decomposed
   from one consistent set of stems (a tag with a measured observation missing
   inside the window uses the exported reconstruction).
4. Summarises the posterior per tag and census pair, weighting each path by
   its share of the tag's samples, `path_count / sum(path_count)`.

## Settings

| Setting | `main_cpp.R` | `main_cpp_chunk.R` | What it does |
| --- | --- | --- | --- |
| `INPUT_FILE` | `data_simulation/data/simulated_data_1.csv` | same | Input CSV |
| `RUN_ALL_TAGS`, `WHICH_TAG` | `FALSE`, `"20"` | not used (all groups) | Tags to run |
| `FORCE_ONE_SPECIES_PARAMETERS` | `TRUE` | `TRUE` | One parameter set for all rows |
| `USE_MEASUREMENT_ERROR` | `TRUE` | `TRUE` | Measurement-error mixture in the estimation and the DP |
| `MAX_GROWTH_FIXED`, `MAX_SHRINK_FIXED` | `7.5`, `-0.5` | `5`, `-0.5` | Hard growth bounds (cm/yr) |
| `K_SHRINK_FIXED`, `K_GROWTH_FIXED` | `0`, `0` | `0`, `0` | Soft penalties |
| `RECRUIT_MAX_FIXED` | `MAX_GROWTH_FIXED × 5 + 0.9999` | same | Largest DBH of a new stem (cm) |
| `PRUNE_BOUND_FACTOR` | `1.25` | `1.25` | Pruning bounds = factor × the hard bounds and the recruit limit |
| `ANCHOR_START_CENSUS` | `7` | `7` | Anchor census |
| `DP_MODE` | `"marginals+bins"` | same | `"map"` and `"marginals"` run the same DP; `"+bins"` adds `DP_PosteriorBin`; `"none"` skips the DP (`main_cpp.R` only) |
| `DP_MAX_STATES` | `40000` | `40000` | DP state limit per census |
| `DP_MAX_TRACKS`, `DP_SLACK_TRACKS` | auto, `1` | same | Track limit and slack tracks |
| `ALLOW_PROVISIONAL_DP_ANCHOR` | `TRUE` | `TRUE` | Provisional anchor when the anchor census has no pins |
| `PIN_TRUESTEMID` | `TRUE` | `TRUE` | Pins in the DP enumeration and the matcher |
| `DP_FALLBACK_GROWTH_FORMS`, `PROB_SPECIES` | none | none | Growth forms and species sent to the matcher |
| `NON_TAPER_CORRECTED_GROWTH_FORMS` | palm, strangler_fig, tree_fern | same | Growth forms with their own pruning bounds and HOM widening |
| `HOM_TOLERANCE_SCALE` | `2.0` | `2.0` | HOM widening of the pruning bounds (cm per m) |
| `PROB_N_SAMPLES`, `PROB_LOOKAHEAD_WEIGHT` | `200`, `1` | same | Matcher samples and lookahead weight |
| `PROB_N_SIGMA_ME` | `3` | engine default `3` | Matcher shrinkage-run cut (`Inf` = off) |
| `USE_BIO_HARD_SHRINK_IN_PROB`, `USE_BIO_HARD_GROWTH_IN_PROB` | `TRUE`, `TRUE` | same | Hard `Bio_*` bounds in the matcher |
| `POSTERIOR_SAMPLES`, `POSTERIOR_SAMPLES_FORMAT` | `200`, `"csv"` | same | Posterior samples per tag and file format |
| `POSTERIOR_SAMPLE_SEED` | `123` | `123` | Seed of the samples |
| `WRITE_DP_CSV`, `WRITE_DP_RDS`, `WRITE_DP_FEATHER`, `WRITE_DP_PDF` | `TRUE`, `TRUE`, `FALSE`, `TRUE` | same | Output files |
| `DP_CHUNK_SIZE` | — | `7` | Groups per chunk |
| `DP_CHUNK_RESUME`, `DP_CHUNK_OVERWRITE`, `DP_CHUNK_START`, `DP_CHUNK_END` | — | `TRUE`, `FALSE`, none, none | Resume and chunk range |
| `MANUAL_CORES`, `MANUAL_CORES_VALUE` | `TRUE`, `1` | same | Workers (`FALSE`: all cores but one; forked, Unix-alikes only) |
| `OUT_DIR_OVERRIDE`, `CONFIG_NAME`, `BATCH_TS` | none | none | Output folder and its name |
| `SENSITIVITY_MODE`, `RUN_REALISM_REPORT`, `RUN_K_SWEEP_DEMO` | `"none"`, `FALSE`, `FALSE` | — | Sensitivity sweeps, realism report, k-sweep plot |

`main_cpp_bci.R` sets `WHICH_TAG = "080297"`, `FORCE_ONE_SPECIES_PARAMETERS =
FALSE`, `MAX_GROWTH_FIXED = 5`, `DP_MAX_STATES = 1039`,
`PRUNE_BOUND_FACTOR = 5`, `USE_MEASUREMENT_ERROR = FALSE`,
`POSTERIOR_SAMPLES = 0` and `POSTERIOR_SAMPLE_SEED = 123` (the seed of the
matcher's draws, so runs are reproducible).

## Requirements

R packages `data.table`, `here`, `Rcpp`, `igraph`, `MASS`, `lpSolve`,
`parallel`; optional `arrow` (Feather), `ggplot2` and `cowplot` (PDF plots),
`reshape2` and `scales` (parameter report, k-sweep plot). `ggplot2`,
`patchwork` and `scales` for `basal_area_uncertainty.R`.
