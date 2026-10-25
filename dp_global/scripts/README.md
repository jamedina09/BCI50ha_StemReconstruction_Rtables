# dp_global/scripts

Driver scripts of the `dp_global` engine, and a basal-area uncertainty script
for their outputs. The drivers load the data, estimate the biological
parameters, run the engine tag by tag, apply the pin sweep and the post-engine
helpers, and write the reconstruction and the posterior path files. The engine
itself is described in `dp_global/README.md`. The BCI production run uses its
own driver, `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [What it does](#what-it-does)
  - [`main_cpp.R`](#main_cppr)
  - [`main_cpp_chunk.R`](#main_cpp_chunkr)
  - [Anchor scoping and pins](#anchor-scoping-and-pins)
  - [Pin sweep, `SweepAuditOverride` and duplicate-aware pinning](#pin-sweep-sweepauditoverride-and-duplicate-aware-pinning)
  - [Post-engine helper chain](#post-engine-helper-chain)
  - [`main_cpp_bci.R`](#main_cpp_bcir)
  - [`basal_area_uncertainty.R`](#basal_area_uncertaintyr)
- [Settings](#settings)
- [Checks](#checks)
- [Requirements](#requirements)

## Contents

| File | What it does |
| --- | --- |
| `main_cpp.R` | Driver for one tag or a small data set; holds the full result in memory; can also run sensitivity sweeps, a realism report and a k-sweep plot |
| `main_cpp_chunk.R` | Driver for large data sets: processes the (Tag, species) groups in chunks, writes each chunk to disk and resumes interrupted runs |
| `main_cpp_bci.R` | Debug driver for single tags of the BCI sample file; sources `main_cpp.R` and sets BCI pins before the engine |
| `basal_area_uncertainty.R` | Basal-area stocks and changes per tag, with the uncertainty from the posterior path files of a run |

## How to run

From the project root (`dp_global/R/dp_global_main.R` finds the modules and
the C++ file from the working directory):

```sh
# One tag
Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WHICH_TAG=20
# All tags of a small data set (every tag in the PDF)
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --PLOT_PDF_ONE_TAG_ONLY=FALSE
# Large data set, in chunks
Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=<file> --DP_CHUNK_SIZE=7 --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=8
# Resume a chunked run in its folder
Rscript dp_global/scripts/main_cpp_chunk.R --OUT_DIR_OVERRIDE=dp_global/output/<run_dir> --DP_CHUNK_RESUME=TRUE
# BCI debug driver
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=187064 --DP_VERBOSE=FALSE
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=000184 --DP_MAX_STATES=2   # tight state limit: the matcher
# Basal area of a finished run
Rscript dp_global/scripts/basal_area_uncertainty.R --RUN_DIR=dp_global/output/<run_dir>
```

**Command line.** Any top-level variable of a driver can be set as
`--NAME=value`: case does not matter and `-` and `_` are equivalent (write text
values such as `WHICH_TAG` with `_`); `true` / `false` become logical, whole
numbers integer, decimals numeric; a flag without a value is `TRUE`. An
unknown name gives a warning and is ignored. `--help` prints the keys of
`CLI_REFERENCE` with their current values (orchestrators can build flags from
these canonical names). `main_cpp_chunk.R` has a reduced `CLI_REFERENCE` (no
`WHICH_TAG`, `PROB_N_SIGMA_ME` or sensitivity and realism keys; it adds the
`DP_CHUNK_*` and `WRITE_DP_PDF_PER_CHUNK` keys) but accepts any variable it
defines. `main_cpp_bci.R` inherits the command line of `main_cpp.R` and
applies the overrides again after its own defaults.

**Interactive use.** `source("dp_global/scripts/main_cpp.R")` defines the
settings and functions without starting a run; `run_main()` starts it. The
same holds for `main_cpp_chunk.R` (`run_main_chunked()`).

## Inputs

| Script | Input |
| --- | --- |
| `main_cpp.R`, `main_cpp_chunk.R` | `INPUT_FILE`, a CSV read with `data.table::fread()`, with at least `Tag`, `CensusID`, `DBH` (cm), `ExactDate` and `TrueStemID`; `OriginalStemID` (or `StemID`), `Status`, `ListOfTSM`, `growth_form` and `HOM` are used when present. Default: `data_simulation/data/simulated_data_1.csv` |
| `main_cpp_bci.R` | `data_simulation/sample_data_BCI/multistem_tags.rds` (`readRDS()`): BCI multiple-stem tags with `Species`, `Tag`, `StemTag`, `OriginalStemID`, `CensusID`, `ExactDate`, `DBH` (cm), `ListOfTSM`, `HOM`, `Status` and the `Bio_*` columns (no `growth_form` column) |
| `basal_area_uncertainty.R` | `<RUN_DIR>/stem_reconstruction_dp_global_rcpp.csv` (required: `Tag`, `CensusID`, `DBH` in cm, `ExactDate`, `ReconstructedStemID`, `obs_row_id`) and the CSV path files `<RUN_DIR>/posteriors/tag_<Tag>_posterior_samples_*_paths.csv` |

## Outputs

Each driver run writes to a new folder under `dp_global/output/`. The name is
built by `build_out_dir_name()` from `BATCH_TS`, `CONFIG_NAME`, the tag
selection (`allT`, `T0` or `T<tag>`), the DP mode and the key settings (the
driver prints it at start); `OUT_DIR_OVERRIDE` writes into a given folder
instead. `PROJECT_ROOT` changes where `dp_global/output/` is created.

| File | Content |
| --- | --- |
| `run_started.txt`, `run_finished.txt` | Start and end time stamps, for job watchers |
| `run_parameters_full.txt` | Run-level settings, command line, overrides and the parameter structure |
| `run_log.txt` | Lines of `log_msg()`; every write through `maybe_write()` is logged as success or failure (e.g. `Wrote RDS chunk 2: <path>`) |
| `stem_reconstruction_dp_global_rcpp.csv` | Reconstruction (`WRITE_DP_CSV`; the chunked driver appends one chunk at a time and writes the header once) |
| `stem_reconstruction_dp_global_rcpp.rds`, `.feather` | The same in binary form (`main_cpp.R` only; `WRITE_DP_RDS`, `WRITE_DP_FEATHER`) |
| `stem_reconstruction_dp_global_rcpp.pdf` | Trajectories per tag (`main_cpp.R`, `WRITE_DP_PDF`) |
| `stem_reconstruction_dp_global_rcpp_chunk_NNN.rds`, `.feather`, `.pdf` | One file per chunk (`main_cpp_chunk.R`; `WRITE_DP_RDS`, `WRITE_DP_FEATHER`, `WRITE_DP_PDF` with `WRITE_DP_PDF_PER_CHUNK`) |
| `stem_reconstruction_dp_global_rcpp_chunk_NNN_done.txt` / `_failed.txt` | Completion marker of a chunk (written last), or the error message of a chunk that failed |
| `posteriors/tag_<Tag>_posterior_samples_<BATCH_TS>_paths.<csv\|rds\|feather>` | Posterior path file of each tag (format in `dp_global/README.md`; an empty `BATCH_TS` gives `tag_11_posterior_samples__paths.csv`) |
| `posteriors/.staging/` | Raw per-sample files of the engines; empty after a completed run |
| `bio_pars_report.pdf` | Plots of the estimated parameters (`main_cpp.R`) |
| `tag_<WHICH_TAG>_realism_summary_rcpp.csv`, `_realism_by_tag_rcpp.csv`, `_realism_tuning_suggestions_rcpp.csv` | Realism report (`main_cpp.R`, `RUN_REALISM_REPORT`) |
| `simulated_all_transition_cost_sweeps_rcpp.rds`, `simulated_all_transition_cost_sweeps.csv`, `simulated_all_transition_cost_sweep_jumps_rcpp.csv`, `simulated_all_transition_cost_jumps_rcpp.csv` / `.rds`, `simulated_all_transition_cost_sweeps_rcpp.pdf` | Sensitivity sweeps (`main_cpp.R`, `SENSITIVITY_MODE` with write; the PDF with `run+write+pdf`) |
| `k_sweep_join_vs_split_demo_rcpp.pdf` | k-sweep plot (`main_cpp.R`, `RUN_K_SWEEP_DEMO`) |

The chunked driver writes no combined RDS, no run-level PDF and no realism
report: it never holds the whole result in memory.

`basal_area_uncertainty.R` writes to `<RUN_DIR>/`:

| File | Content |
| --- | --- |
| `basal_area_tag_census.csv` | Per tag and census: `Tag`, `CensusID`, `Year` (mean ExactDate), `TotalBA_m2` (Σ π/4 (DBH/100)² over the measured stems), `NumStems` |
| `basal_area_tag_change.csv` | Per tag and census pair: `CensusID_from`, `CensusID_to`, `Interval_yr`; exported reconstruction `Growth_BA` (surviving stems), `Loss_BA` (deaths, negative), `Gain_BA` (recruits), `DeltaBA_total` (identity-free), `NumSurvivors`, `NumDeaths`, `NumRecruits`; posterior `*_mean`, `*_sd`, `*_q025`, `*_q975` of growth, loss and gain; `DeltaBA_check` (weighted mean of the total change, equal to `DeltaBA_total`); `NumSurvivors_mean`, `NumDeaths_mean`, `NumRecruits_mean`; `NumPaths` |
| `basal_area_figures.pdf` | Summary pages, per-tag panels (BA trajectory, stem count, decomposition with uncertainty whiskers, demographics), a 3 × 3 page of the exported trajectory against eight posterior paths drawn by weight, per-stem BA trajectories, uncertainty histograms and posterior densities (pooled and per interval) |

## What it does

### `main_cpp.R`

1. **Load and prepare.** Reads `INPUT_FILE`; sets the parameter set of each
   row in `species` (`ensure_species_column()`): one set for all rows with
   `FORCE_ONE_SPECIES_PARAMETERS = TRUE` (label `FORCED_SPECIES_LABEL`,
   `"all"`), else the column `SPECIES_COL` or the first of the usual species
   column names.
2. **Parameters.** Estimates one parameter set per species
   (`estimate_bio_pars()`, censuses from `ANCHOR_START_CENSUS` on, fixed hard
   bounds and soft penalties, growth increments outside the fixed bounds and
   recruits above `RECRUIT_MAX_FIXED` dropped before the fit) and attaches the
   `Bio_*` columns (`attach_bio_columns()`).
3. **Engine.** Runs `run_dp_one_group()` on each (Tag, species) group: one tag
   (`RUN_ALL_TAGS = FALSE`, `WHICH_TAG`, matched exactly, so `"084555"` keeps
   its leading zero) or all tags in parallel over `MC_CORES` workers.
   `run_dp_one_group()` calls the DP with the driver settings (pruning bounds
   `PRUNE_BOUND_FACTOR` × the fixed limits, `prune_use_bio_bounds = FALSE`,
   `prune_use_bio_recruit = FALSE`, `temperature = 1`), sends the group to the
   probabilistic matcher when the DP raises an error (`DP_FallbackReason =
   "error:<message>"`), and applies the script-level pin sweep (below).
4. **Helpers.** On the assembled result: `maybe_add_posterior_bins()`, then the
   post-engine helper chain and `finalize_posterior_paths()` (below).
5. **Outputs** (table above), plus the optional sensitivity sweeps
   (`SENSITIVITY_MODE`), realism report (`RUN_REALISM_REPORT`, for a
   representative species), k-sweep plot (`RUN_K_SWEEP_DEMO`) and the
   parameter report.

`PLOT_PDF_ONE_TAG_ONLY` is set from the default `RUN_ALL_TAGS` before the
command line is applied, so a run with `--RUN_ALL_TAGS=TRUE` also needs
`--PLOT_PDF_ONE_TAG_ONLY=FALSE` to put every tag in the PDF. `main_cpp.R`
processes all groups in one call and cannot resume: a stopped run starts
again; use `main_cpp_chunk.R` when checkpoints are needed.

### `main_cpp_chunk.R`

The same steps, in chunks, to keep peak memory low:

1. Groups are `unique(xrun[, .(Tag, species)])`, sorted, split into chunks of
   `DP_CHUNK_SIZE` groups (`DP_CHUNK_START` / `DP_CHUNK_END` limit the range).
2. For each chunk the groups run in parallel (`parallel::mclapply()`, each
   worker with `data.table::setDTthreads(1)`); a group without any DBH or
   CensusID is returned as is (`skipped_no_data`).
3. The chunk result (`out_chunk`, with `DP_Chunk` and `run_out_dir`) goes
   through `maybe_add_posterior_bins()`, the post-engine helper chain and
   `finalize_posterior_paths()`, is appended to the CSV and written as chunk
   files, and is removed from memory before the next chunk.
4. The `_done.txt` marker is written last. With `DP_CHUNK_RESUME = TRUE` a
   chunk with a marker is skipped; a chunk without one is run again, even if a
   partial chunk file exists. `DP_CHUNK_OVERWRITE = TRUE` reruns finished
   chunks. A chunk that raises an error leaves `_failed.txt` and the run goes
   on.

The chunked driver always runs every group: `RUN_ALL_TAGS` and `WHICH_TAG`
only enter the folder name. It runs the DP for every `DP_MODE` (`"+bins"`
adds `DP_PosteriorBin`). It has no sensitivity sweeps or realism report: both
need the full result in memory.

Memory advice for large data: keep the incremental CSV (`WRITE_DP_CSV`, which
grows on disk) and the per-chunk RDS files (`WRITE_DP_RDS`, to rebuild the CSV
or reload chunks later), and avoid assembling the whole result in memory.

**After a run.** `merge_chunks_to_csv(<run_dir>, prefer = "rds")` (defined in
`main_cpp_chunk.R`; `prefer = "feather"` for Feather chunks) streams the chunk
files one at a time into `stem_reconstruction_dp_global_rcpp_merged.csv`;
`merge_chunk_rds_to_csv()` and `merge_chunk_feathers_to_csv()` do it for one
kind. On a machine with enough memory:

```r
library(data.table)
files <- list.files("<run_dir>", pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
all <- rbindlist(lapply(files, readRDS), use.names = TRUE, fill = TRUE)
```

When a run is resumed with `--OUT_DIR_OVERRIDE` and the CSV exists, the new
chunks are appended to it. All settings other than the output ones must match
the original run, so the resumed chunks are processed in the same way.

### Anchor scoping and pins

- With observations after `ANCHOR_START_CENSUS`, the DP is solved up to the
  anchor and the later rows are appended: rows with a `TrueStemID` (measured or
  not) get it (`given`), the others `none_after_anchor` (rules in
  `dp_global/README.md`, "Anchor census").
- With `ALLOW_PROVISIONAL_DP_ANCHOR = TRUE`, an anchor census whose measured
  rows lack pins gets provisional IDs (`provisional_dp`) instead of a fallback.
- With `PIN_TRUESTEMID = TRUE` the rows with a `TrueStemID` before the anchor
  are constrained to it in the DP and the matcher; with `FALSE` the engines are
  free to give them another stem.

### Pin sweep, `SweepAuditOverride` and duplicate-aware pinning

With `PIN_TRUESTEMID = TRUE` an idempotent sweep sets every row with a
`TrueStemID` to `ReconstructedStemID = TrueStemID` and
`ReconstructionMethod = "given"`. It runs inside `finalize_out()` (the DP's
fallback paths and resprout segment splits), inside
`match_stems_probabilistic()`, and in the drivers' `run_dp_one_group()` on
every output (the only sweep for tags the DP solves directly, since the DP
enforces the pins through its state enumeration). It covers rows the engines
never visit (unmeasured terminal rows pinned by the BCI steps, re-inserted
missing-from-field rows, matcher edge cases).

- **`SweepAuditOverride`**: `TRUE` where the engine had assigned another stem
  and the sweep replaced it; an `[audit]` log line names the tag and the count.
  On such rows the `DP_PosteriorTop*` columns describe the engine's choice, not
  the final stem; treat the row as observed (probability 1).
- **`ReconstructedStemID_PreSweep`**: the engine's stem, taken once by the
  first sweep that runs; equal to `ReconstructedStemID` where nothing was
  overridden.
- **Duplicate-aware pinning** (`run_dp_one_group()`): the rows are processed in
  stable order; before setting a pin, the sweep checks whether that value is
  already on another row of the same (Tag, CensusID) in the working column
  (engine values plus the pins already set). If so it keeps the engine's
  stem, leaves `ReconstructionMethod` unchanged and flags
  `SweepRollbackToPreSweep = TRUE` (with an `[audit]` log line). This happens
  when the database gives two stems of one census the same StemID.
  `SweepAuditOverride` records the disagreement, `SweepRollbackToPreSweep` the
  resolution.

The three columns are created on every output, also for tags without pins.

### Post-engine helper chain

Every driver applies, in this order (all in `dp_global/R/dp_global_main.R`;
rules in `dp_global/README.md`, "Post-engine helpers"):

1. `maybe_add_posterior_bins()` (`main_cpp.R`: on the assembled result;
   `main_cpp_chunk.R`: per chunk);
2. `apply_pin_track_rejoin()`: unpinned rows the sweep left behind rejoin
   their track's pin (`pin_track`);
3. `apply_carried_terminal_backfill()`: rows with no stem, no DBH and status
   dead, stem dead, broken below or missing take, after sorting by Tag, source
   ID and CensusID, the latest earlier stem of the same (Tag, source ID) (LOCF;
   source ID `StemID`, else `OriginalStemID`) (`carried_terminal`);
4. `apply_orphan_stem_backfill()`: rows with no stem, no pin, no DBH and a
   source ID take the source ID (`given_orphan`; e.g. a new StemID first
   recorded broken below without a DBH);
5. `apply_terminal_to_host()`: unmeasured terminal records that precede the
   life of their stem move to the stem that ended (`terminal_to_host`);
6. `apply_broken_below_invariants()`: rules R1 and R2;
7. `renumber_engine_minted_ids(out, posterior_top_k = DP_POSTERIOR_TOP_K, …)`:
   stems numbered 1..N per tag (first census, then largest DBH, then old ID;
   never zero or negative); returns the Tag / old_id / new_id mapping (its
   `posterior_samples_path` and `mapping_format` arguments are not used);
8. `finalize_posterior_paths(out, posterior_samples_path = out_dir, mapping = …)`:
   for each staging file, translates the labels with the mapping, applies the
   pins (`apply_pins_to_samples()`) and R1/R2 (`apply_bb_invariants_to_samples()`)
   to each sample, computes `path_sig`, `path_count`, `path_prob` and `recon`,
   writes the path file and deletes the staging file. `path_prob` is the count
   share when the samples carry no `logp` (matcher) and the share of
   `exp(logp)` when they do (DP).

In `main_cpp.R` these are steps 5.5–5.5e of `run_main()`; in
`main_cpp_bci.R` section 9b.

### `main_cpp_bci.R`

A debug driver for single BCI tags. It sources `main_cpp.R` for its functions
and command line, overrides the input and a few defaults, estimates no
parameters (the input carries the `Bio_*` columns) and sets pins before the
engine from the BCI database conventions:

- **Step 1(a)**: a row with a `StemTag` (a stem tagged in the field) is pinned
  to its `OriginalStemID`.
- **Step 1(b)**: a row from census 7 on (2010, when every stem was tagged) is
  pinned to its `OriginalStemID`, except rows without a DBH whose `Status` is
  dead, stem dead or broken below (the end of a trajectory, placed by the
  engine and the backfills).
- **Step 2**: within each (Tag, `OriginalStemID`), the rows after the last DBH
  form the terminal phase (2a). Their unpinned rows take a pin of the same
  phase, forward then backward (2b).
- **Step 3a**: a row with a DBH that is broken below or carries a resprout code
  (`R`, `RP`, `RF`, `RT`, `QR`, `OR`) is pinned to its `OriginalStemID` (the
  start of a new trajectory); death and break rows without a DBH are not.
- **Step 3a.5**: an `OriginalStemID` with no pin yet, at least one alive
  measured row and at least one dead / stem-dead / broken-below row without a
  DBH gets its alive measured rows pinned to it. In the BCI data a terminal
  record without a DBH nearly always continues its own StemID (99.8% of dead
  and stem-dead records, 55.6% of broken-below records have an earlier record
  of the same StemID); without this step the DP could give the unpinned alive
  history to a competing new stem and duplicate stems at later censuses.
- **Step 3b**: within each (Tag, `OriginalStemID`), when the measured rows hold
  exactly one pin value, it fills the unpinned rows; groups with conflicting
  pins are left alone and counted.

Rows without an unambiguous identity stay free for the engine. The production
driver has further pin rules: re-keyed field tags (1a'), Step 3a only for
StemIDs with another record, the gap and fast-link stops (3b, 3c) and hard
pins (`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/README.md`).

Defaults that differ from `main_cpp.R`: `INPUT_FILE` (above),
`WHICH_TAG = "080297"`, `FORCE_ONE_SPECIES_PARAMETERS = FALSE`,
`MAX_GROWTH_FIXED = 5`, `MAX_SHRINK_FIXED = -0.5`, `DP_MAX_STATES = 1039`,
`PRUNE_BOUND_FACTOR = 5`, `USE_MEASUREMENT_ERROR = FALSE`,
`POSTERIOR_SAMPLES = 0` (`--POSTERIOR_SAMPLES=<n>` to sample),
`POSTERIOR_SAMPLE_SEED = 123` (the seed of the matcher's draws, so runs are
reproducible), `PROB_N_SIGMA_ME = 3`. Output goes to
`dp_global/output/<timestamp>_BCI_tag<WHICH_TAG>_…`.

### `basal_area_uncertainty.R`

The basal area (BA) of a tag in a census, the sum of π/4 (DBH/100)² over its
measured stems (m²), does not depend on which stem a measurement belongs to.
The split of each change into growth of surviving stems, loss by mortality and
gain by recruitment does: different posterior paths make different stems
survivors, deaths or recruits.

1. Reads the reconstruction and the posterior path files of the run; without
   path files it gives the exported reconstruction only. DBH is in cm. Only CSV
   path files are read (`POSTERIOR_SAMPLES_FORMAT = "csv"`, the default of
   `main_cpp.R` and `main_cpp_chunk.R`).
2. Sums the BA of the measured stems per tag and census.
3. Splits each change between censuses into growth, loss and gain for the
   exported reconstruction (`decompose_ba_change()`: a stem in both censuses
   is a survivor, in the first only a death, in the second only a recruit).
4. Does the same for every posterior path. The paths hold the observations up
   to the anchor; an observation outside the census window of the paths
   (after the anchor, or before the paths' first census) joins, in every path,
   the path stem holding its exported stem-mate at the nearest window edge
   (the last in-path observation of the same exported stem for later censuses,
   the first for earlier ones); an exported stem without in-path observations
   keeps its own stem. Every census pair is thus decomposed from one consistent
   set of stems. When a measured observation inside the window is missing from
   the paths, the tag's paths are replaced by the exported reconstruction. The
   script reports both counts.
5. Summarises the posterior per tag and census pair, weighting each path by its
   share of the tag's samples, `path_count / sum(path_count)` (the `path_prob`
   column is not used: for the DP it counts each sample's probability twice):
   weighted mean, SD and 2.5% / 97.5% quantiles.
6. Writes the tables and the figures.

## Settings

| Setting | `main_cpp.R` | `main_cpp_chunk.R` | What it does |
| --- | --- | --- | --- |
| `INPUT_FILE` | `data_simulation/data/simulated_data_1.csv` | same | Input CSV |
| `RUN_ALL_TAGS`, `WHICH_TAG` | `FALSE`, `"20"` | not used for selection (`"0"`) | Tags to run |
| `FORCE_ONE_SPECIES_PARAMETERS`, `FORCED_SPECIES_LABEL`, `SPECIES_COL` | `TRUE`, `"all"`, auto | same | Parameter sets |
| `USE_MEASUREMENT_ERROR` | `TRUE` | `TRUE` | Measurement-error mixture in the estimation and the DP growth likelihood (the matcher has none) |
| `MAX_GROWTH_HARD_SOURCE`, `MAX_GROWTH_FIXED` | `"fixed"`, `7.5` | `"fixed"`, `5` | Hard upper growth bound (cm/yr) |
| `MAX_SHRINK_HARD_SOURCE`, `MAX_SHRINK_FIXED` | `"fixed"`, `-0.5` | same | Hard lower growth bound (cm/yr) |
| `K_SHRINK_SOURCE`, `K_SHRINK_FIXED`, `K_GROWTH_SOURCE`, `K_GROWTH_FIXED` | `"fixed"`, `0`, `"fixed"`, `0` | same | Soft penalties (0 = off) |
| `RECRUIT_MAX_SOURCE`, `RECRUIT_MAX_FIXED` | `"fixed"`, `MAX_GROWTH_FIXED × 5 + 0.9999` | same | Largest DBH of a new stem (cm) |
| `PRUNE_BOUND_FACTOR` | `1.25` | `1.25` | Pruning bounds = factor × the hard bounds and the recruit limit |
| `ANCHOR_START_CENSUS` | `7` | `7` | Anchor census |
| `DP_MODE` | `"marginals+bins"` | same | `"map"` and `"marginals"` run the same DP (MAP plus marginals) and differ in the folder label; `"marginals+bins"` adds `DP_PosteriorBin`; `"none"` skips the DP (`main_cpp.R` only) |
| `DP_VERBOSE` | `TRUE` | `TRUE` | Progress messages (matcher warnings go to stderr; redirect both streams to log them) |
| `DP_POSTERIOR_TOP_K` | `2` | `2` | Number of `DP_PosteriorTop<k>ID` / `Prob` pairs |
| `DP_MAX_STATES` | `40000` | `40000` | State limit per census (and `max_states²` per census pair); with K = n + 1 the DP solves up to 6 stems per census exactly (tables in `dp_global/README.md`) |
| `DP_MAX_TRACKS` | auto | auto | Track cap (auto: largest stem count of a census up to the anchor, plus 1) |
| `DP_SLACK_TRACKS` | `1` | `1` | Slack tracks |
| `DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE`, `DP_SLACK_REQUIRE_ANCHOR_EPS` | `TRUE`, `1e-6` | same | Slack tracks only when an anchor stem is no larger than the recruit limit |
| `ALLOW_PROVISIONAL_DP_ANCHOR` | `TRUE` | `TRUE` | Provisional anchor IDs instead of a fallback |
| `PIN_TRUESTEMID` | `TRUE` | `TRUE` | Pins in the DP enumeration, the matcher and the sweep |
| `DP_FALLBACK_GROWTH_FORMS` | none | none | `growth_form` values sent to the matcher (comma- or semicolon-separated; matched exactly) |
| `PROB_SPECIES` | none | none | `Species` values sent to the matcher |
| `NON_TAPER_CORRECTED_GROWTH_FORMS` | palm, strangler_fig, tree_fern | same | Growth forms with their own pruning bounds and HOM widening |
| `NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH`, `_MAX_GROWTH` | `PRUNE_BOUND_FACTOR` × the hard bounds | same | Their pruning bounds (equal to the general ones by default, so only the HOM widening differs) |
| `HOM_TOLERANCE_SCALE` | `2.0` | `2.0` | HOM widening per census pair, `scale × max(abs(HOM − 1.3)) / interval` (cm per m; 0 = off) |
| `PROB_N_SAMPLES` | `200` | `200` | Matcher samples |
| `PROB_LOOKAHEAD_WEIGHT` | `1` | `1` | Matcher look-ahead weight (0 = off; active with K ≥ 4) |
| `PROB_N_SIGMA_ME` | `3` | engine default `3` | Matcher cumulative-shrinkage cut (`Inf` = off) |
| `USE_BIO_HARD_SHRINK_IN_PROB` | `TRUE` | `TRUE` | Matcher: `Bio_Max_Shrink` gate on links and the cumulative-shrinkage cut; `FALSE` allows large shrinkage (soft `k_shrink` only). The DP is not affected |
| `USE_BIO_HARD_GROWTH_IN_PROB` | `TRUE` | `TRUE` | Matcher: `Bio_Max_Growth` gate on links. The DP is not affected |
| `POSTERIOR_SAMPLES` | `200` | `200` | Posterior samples per tag (0 = none) |
| `POSTERIOR_SAMPLES_FORMAT` | `"csv"` | `"csv"` | Path-file format (`feather` needs `arrow`, else `rds`) |
| `POSTERIOR_SAMPLES_PATH` | run folder | run folder | Base folder of `posteriors/` (a path ending in `posteriors` is trimmed) |
| `POSTERIOR_SAMPLE_SEED` | `123` when sampling | same | Seed of the samples and of the matcher's draws |
| `WRITE_DP_CSV`, `WRITE_DP_RDS`, `WRITE_DP_FEATHER`, `WRITE_DP_PDF` | `TRUE`, `TRUE`, `FALSE`, `TRUE` | same | Output files (Feather needs `arrow`) |
| `WRITE_DP_PDF_PER_CHUNK` | — | `TRUE` | Per-chunk PDFs (with `WRITE_DP_PDF`) |
| `DP_PDF_INCLUDE_REFERENCE` | `TRUE` | `TRUE` | Second PDF panel with the trajectories grouped by `OriginalStemID` |
| `PLOT_PDF_ONE_TAG_ONLY` | `TRUE` with one tag | — | PDF of `WHICH_TAG` only |
| `DP_CHUNK_SIZE` | — | `7` | Groups per chunk |
| `DP_CHUNK_RESUME`, `DP_CHUNK_OVERWRITE` | — | `TRUE`, `FALSE` | Skip finished chunks; rerun them |
| `DP_CHUNK_START`, `DP_CHUNK_END` | — | none | Chunk range (tests) |
| `MANUAL_CORES`, `MANUAL_CORES_VALUE` | `TRUE`, `1` | same | Workers; `FALSE` uses all cores but one (forked with `mclapply()`, more than one only on Unix-alikes) |
| `OUT_DIR_OVERRIDE`, `CONFIG_NAME`, `BATCH_TS`, `PROJECT_ROOT` | none, none, `""`, project root | same | Output folder, its name and its root |
| `SENSITIVITY_MODE` | `"none"` | — | `"run"`, `"run+write"`, `"run+write+pdf"`: sensitivity sweeps, written, with PDF |
| `RUN_REALISM_REPORT`, `RUN_K_SWEEP_DEMO` | `FALSE`, `FALSE` | — | Realism report; k-sweep plot |

The engine defaults used by these drivers without a setting: birth-death
assignment in the matcher (`prob_birth_death = TRUE`) and no DBH classes
(`dbh_round_censuses`); the BCI driver exposes both (`PROB_BIRTH_DEATH`,
`DBH_ROUND_CENSUSES`).

## Checks

- `POSTERIOR_SAMPLES` and `POSTERIOR_SAMPLE_SEED` must be non-negative
  integers and `DP_MODE` and `POSTERIOR_SAMPLES_FORMAT` allowed values; the
  drivers stop otherwise. A missing `INPUT_FILE` stops the run.
- `[audit]` lines in `run_log.txt` report pin overrides and duplicate
  rollbacks per tag.
- `run_finished.txt`, the `_done.txt` / `_failed.txt` markers and an empty
  `posteriors/.staging/` show a completed run.
- `basal_area_uncertainty.R` reports the spliced observations and the tags that
  fell back to the exported reconstruction; `DeltaBA_check` equals
  `DeltaBA_total`.

## Requirements

R packages `data.table`, `here`, `Rcpp`, `igraph`, `MASS` (loaded by
`dp_global_main.R`), `lpSolve` (matcher), `parallel`; optional `arrow`
(Feather), `ggplot2` and `cowplot` (PDF plots), `reshape2` and `scales`
(parameter report, k-sweep plot). `basal_area_uncertainty.R` needs
`data.table`, `ggplot2`, `patchwork` and `scales` (`here` resolves a relative
`RUN_DIR` when installed). A C++ compiler for the engine.
