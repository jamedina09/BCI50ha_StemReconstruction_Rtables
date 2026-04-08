# dp_global/scripts

This document describes the current behavior of the `dp_global` driver scripts:
- `dp_global/scripts/main_cpp.R` — the interactive/CLI driver for single-tag or targeted runs.
- `dp_global/scripts/main_cpp_chunk.R` — the chunked driver optimized for large runs.
- `dp_global/scripts/main_cpp_bci.R` — the BCI debug driver for single-tag runs on BCI census data.
- `dp_global/scripts/basal_area_uncertainty.R` — posterior-based basal area uncertainty quantification.

Both `main_cpp.R` and `main_cpp_chunk.R` accept command-line overrides of defaults using `--KEY=VALUE` flags. Keys are case-insensitive and may use `-` or `_` as separators. `main_cpp_bci.R` inherits the same CLI interface from `main_cpp.R`.

---

## Overview ✨

`main_cpp.R` is the central driver for the `dp_global` pipeline. It:

- Loads input tree census data (CSV via `input_file`).
- Ensures species information is present and optionally forces a single species label.
- Estimates biological parameters per species (via `estimate_bio_pars`).
- Runs Dynamic Programming (DP) reconstruction logic (via `run_dp_one_group` which calls the Rcpp-based DP functions).
- Optionally runs sensitivity sweeps and realism checks.
- Writes outputs (CSV, RDS, PDF) to an automatically created output directory.

You can run either script directly via Rscript:

- Single-tag or small dataset: `Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=... --RUN_ALL_TAGS=TRUE`
- Large dataset (chunked, recommended): `Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=...`

(When running from the `dp_global` folder: `Rscript scripts/main_cpp.R`.)

---

## Important behavior notes ⚠️

- `main_cpp.R` defaults to single-tag mode (`RUN_ALL_TAGS = FALSE`). Use `--RUN_ALL_TAGS=TRUE` to run across all tags.
- For large datasets, prefer `main_cpp_chunk.R` which divides groups (Tag + species) into chunks (`DP_CHUNK_SIZE`) and writes chunk outputs incrementally to disk to keep peak RAM low.
- When chunking, the combined `out` is not assembled in memory. Instead:
  - Each chunk is appended to `stem_reconstruction_dp_global_rcpp.csv` (if `WRITE_DP_CSV=TRUE`),
  - Each chunk is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds` (if `WRITE_DP_RDS=TRUE`) and serves as a resume/completion marker.
- `maybe_add_posterior_bins()` is applied per-chunk when chunking; in non-chunked runs it is applied after assembling the full `out`.
- The chunked driver disables or comments out features that are not applicable to chunked runs (for example: full-run sensitivity sweeps and realism report generation).

---

## Key CLI flags & defaults 🔧

Common flags used by both drivers (case-insensitive, but use capital letters to keep it clean):

- `INPUT_FILE` — default: `data_simulation/data/simulated_data_1.csv`.
- `FORCE_ONE_SPECIES_PARAMETERS` — default: `TRUE`.
- `FORCED_SPECIES_LABEL` — default: `"all"` (used when forcing one species label).
- `SPECIES_COL` — default: `NULL` (auto-detected if not set).
- `USE_MEASUREMENT_ERROR` — default: `TRUE`.

DP / reconstruction option
- `DP_MODE` — default: `"marginals+bins"`. Allowed: 'none' 'marginals' 'marginals+bins' 'map'
- `WHICH_TAG` — character; used for single-tag runs (relevant to `main_cpp.R`). Must match the `Tag` column exactly (e.g., `"084555"` preserves leading zeros). The chunked runner processes groups (`Tag`, `species`) and does not rely on `WHICH_TAG`.
- `ANCHOR_START_CENSUS` — default: `7L`.
- `ALLOW_PROVISIONAL_DP_ANCHOR` — default: `TRUE` — when `TRUE` the DP can assign provisional anchor IDs at the last observed DBH census if the requested anchor census lacks `TrueStemID` but has DBH; set to `FALSE` to require an explicit anchored census or to fall back to the probabilistic matcher.
- `DP_VERBOSE` — default: `TRUE`.
- `DP_POSTERIOR_TOP_K` — default: `2L` - top=k posterior reconstructions to track. 
- `DP_MAX_TRACKS` — default: `NULL` (auto-computed per data when `NULL`) — optionally force max tracks; if `NULL` auto-computed.
- `DP_MAX_STATES` — default: `40000L`. Maximum number of injective assignment states allowed per census. Also controls the inter-census transition limit (`max_edges = max_states²`). The number of states at a census with $n$ observed stems and $K$ tracks is $P(K,n) = K!/(K-n)!$. When any census exceeds `max_states`, or when the cross-product of states between two adjacent censuses exceeds `max_states²`, the solver falls back to the probabilistic matcher. **Practical limits with default 40,000:** DP handles up to 6 observed stems per census exactly (P(7,6) = 5,040 < 40,000); 7+ stems trigger fallback (P(8,7) = 40,320 > 40,000). With `DP_MAX_STATES = 1,000`: max 5 stems. With `DP_MAX_STATES = 20,000`: max 6 stems. See `dp_global/README.md` for detailed tables and how to choose a value.
- `DP_SLACK_TRACKS` — default: `1L` - slack (additional) tracks allowed for DP.
- `DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE` — default: `TRUE` - require anchor to be recruitable (if DBH less than max recruitment size) before granting slack.
- `DP_SLACK_REQUIRE_ANCHOR_EPS` — default: `1e-6`.
- `PROB_N_SAMPLES` — default: `200L`. Number of Gumbel-noise stochastic samples drawn by the probabilistic greedy matching fallback. Higher values produce more accurate posterior estimates at the cost of computation time.
- `PROB_LOOKAHEAD_WEIGHT` — default: `0.5`. Weight for sequential backward conditioning in the probabilistic matcher. When > 0 and $K \geq 4$, the cost matrix for each census pair is conditioned on the already-resolved forward assignment, improving path continuity. Set to `0` to disable conditioning.
- `DP_FALLBACK_GROWTH_FORMS` — default: `character(0)`; comma- or
  semicolon-separated list of values in the `growth_form` column that should
  trigger an immediate probabilistic fallback and prevent the DP solver from running
  on that tag. The driver automatically splits the string into a vector.
- `NON_TAPER_CORRECTED_GROWTH_FORMS` — default: `c("palm", "strangler_fig", "tree_fern")`;
  growth forms whose DBH is NOT taper-corrected. These forms exhibit both real
  biological DBH growth (palms: 1–3 cm/yr; strangler figs: variable as they
  encircle hosts) and apparent DBH shifts when the measurement height (HOM)
  changes between censuses. Because their trunk geometry does not allow taper
  correction, any HOM change produces an uncompensated apparent DBH change.
  The DP solver replaces the general pruning bounds with wide base bounds
  (see below) and optionally applies HOM-proportional widening on top.
- `NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH` — default: `-0.625` cm/year
  (`1.25 × MAX_SHRINK_FIXED`); lower annual growth bound used during pruning
  for non-taper-corrected growth forms. Replaces the general effective pruning
  minimum when the tag's `growth_form` matches `NON_TAPER_CORRECTED_GROWTH_FORMS`.
  Note: in the default driver configuration the general `prune_min_growth` is
  also `1.25 × MAX_SHRINK_FIXED`, so this override is only active when
  `prune_use_bio_bounds = TRUE` or narrower general bounds are set.
- `NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH` — default: `6.25` cm/year
  (`1.25 × MAX_GROWTH_FIXED`); upper annual growth bound for non-taper-corrected
  growth forms during pruning. Same interplay with the general bounds as above.
- `HOM_TOLERANCE_SCALE` — default: `2.0` (cm/year per meter of HOM deviation);
  when a non-taper-corrected tag has a `hom` (or `HOM`) column in its data,
  the prune bounds are widened for each census pair by
  `hom_tolerance_scale × max(|HOM − 1.3|) / interval_years`.
  This is where the effective differentiation between regular trees and
  non-taper-corrected forms occurs in the default configuration. NA HOM values
  are treated as 1.3 (zero contribution). Set to `0` to disable
  HOM widening while keeping the wide base bounds.

**Anchor scoping and post-anchor preservation:** If observations exist after the requested `ANCHOR_START_CENSUS`, the DP is scoped to censuses <= `ANCHOR_START_CENSUS` and post-anchor rows are preserved and appended to the output. Post-anchor rows with non-NA `DBH` and a `TrueStemID` that was used by the DP are set to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. Remaining post-anchor rows without DP assignments receive `ReconstructionMethod = "none_after_anchor"`. If scoping removes all pre-anchor observations, the original rows are returned and observed `TrueStemID` values are treated as `given` while other rows are labeled `none_after_anchor`.

**Provisional anchor behavior:** When a requested anchor census lacks `TrueStemID` but contains DBH observations and `ALLOW_PROVISIONAL_DP_ANCHOR=TRUE`, the DP will assign provisional anchor IDs at the last-observed DBH census and mark those anchor rows with `ReconstructionMethod = "provisional_dp"`.

Posterior sampling:
- `POSTERIOR_SAMPLES` — default: `200L` (set to `0` to disable sampling). When `>0`, the DP engine will draw full-path posterior samples and write them into a per-run `posteriors/` subdirectory.
- `POSTERIOR_SAMPLES_FORMAT` — default: `"csv"` (options: `rds`, `feather`, `csv`).
- `POSTERIOR_SAMPLES_PATH` — default: `NULL`. If `NULL` the run's `out_dir` is used and the DP writes to `<out_dir>/posteriors/`. If you supply a path that itself ends in `posteriors`, the script strips that suffix to avoid creating nested `posteriors/posteriors` folders (the DP will create the `posteriors/` subdirectory itself).
- `POSTERIOR_SAMPLE_SEED` — default: `NULL`. If sampling is enabled and the seed is unset, the chunked runner defaults the seed to `123L` to improve reproducibility; you can override this with `--POSTERIOR_SAMPLE_SEED=<int>`.

Example posterior output files (per tag/run): `tag_11_posterior_samples__summary.csv`, `tag_11_posterior_samples__paths.csv`.

Output controls:
- `WRITE_DP_CSV` — default: `TRUE` - write incremental/combined CSV output - memory heavy.
- `WRITE_DP_RDS` — default: `TRUE` - write per-chunk RDS or combined for non-chunk runs.
- `WRITE_DP_FEATHER` — default: `FALSE` (requires the `arrow` package) - write per-chunk feather (.feather) files or combined for non-chunk runs.
- `WRITE_DP_PDF` — default: `TRUE` - write pdf visualizations per tag - memory heavy.
- `DP_PDF_INCLUDE_REFERENCE` — `main_cpp.R` default: `TRUE`; `main_cpp_chunk.R` default: `TRUE` — include biologically-informed reference lines in PDFs (useful with simulated data).
- `WRITE_DP_PDF_PER_CHUNK` — default in `main_cpp_chunk.R`: `TRUE` (controls per-chunk PDFs).

Parallel & chunking controls (chunked runner specific):
- `DP_CHUNK_SIZE` — default in `main_cpp_chunk.R`: `7L` (set `<= 0` to disable chunking behavior when applicable).
- `DP_CHUNK_RESUME` — default: `TRUE` (skip chunks whose `_done.txt` completion marker exists) — allows stopping and resuming runs. A chunk is considered complete only when its `_done.txt` file is present; partial RDS files from interrupted runs are re-processed.
- `OUT_DIR_OVERRIDE` — default: `NULL`. When set, bypasses automatic output directory creation and writes into the specified path directly. Use this to resume into an existing output directory (e.g., `--OUT_DIR_OVERRIDE=dp_global/output/<previous_run_dir>`).
- `DP_CHUNK_OVERWRITE` — default: `FALSE` (when `TRUE`, overwrite existing chunk RDS files).
- `DP_CHUNK_START`, `DP_CHUNK_END` — default: `NULL` (limit chunk range for tests).
- `RUN_ALL_TAGS` — default: `FALSE` (when `TRUE` run across all tags; chunking is recommended for large runs).
- `MANUAL_CORES` & `MANUAL_CORES_VALUE` — default: `TRUE` and `1L` respectively.

Notes on CLI differences:
- `main_cpp_chunk.R` exposes a reduced `CLI_REFERENCE` relative to `main_cpp.R` (it omits `WHICH_TAG` and leaves sensitivity/realism flags commented out) to reflect the chunked runner's intent; however it still accepts CLI overrides for many run-level parameters when invoked via Rscript.

Helpful post-run utilities:
- Merge chunk RDS/Feather files into a single CSV (run this in R or source the script and call the helper):

```r
# from R in project root
source("dp_global/scripts/main_cpp_chunk.R")
merge_chunks_to_csv("dp_global/output/<your_run_dir>")
```

This streams each chunk file to a single CSV to avoid loading the full dataset into memory.

`main_cpp.R` runs the non-chunked workflow (single-tag or parallelized tags) and does not perform per-chunk writing.

Misc:
- `USE_MEASUREMENT_ERROR` (default: `TRUE`) — enable measurement-error-aware parameter estimation.

Output directory & naming:
- `PROJECT_ROOT` (default: project root via `here::here()`) — override to set a different project root and thus change where `dp_global/output/` is created.
- `base_out_dir` (default: `dp_global/output`) — base directory where run-specific output directories are created.
- `CONFIG_NAME` (default: `NULL`) — optional string used when assembling the run-specific output directory name.
- The final `out_dir` is automatically constructed from timestamp, config name, DP mode, and other key parameters. The script writes a `run_parameters_full.txt` file into `out_dir` documenting the run configuration.

Files produced as run markers/logs:
- `run_started.txt` and `run_finished.txt` — small timestamp files written at start and finish to allow job watchers to detect progress.
- `run_parameters_full.txt` — text file capturing all important run-level variables for reproducibility.
- `run_log.txt` — appended by `log_msg()` throughout the run; writes performed via `maybe_write()` ensure directories exist and the script records success/failure messages here (for example: `Wrote RDS chunk 2: <path>`).

PDF & plotting controls:
- `WRITE_DP_PDF` (default: `TRUE`) — control whether PDFs are generated via `plot_tag_to_pdf()`.
- `DP_PDF_INCLUDE_REFERENCE` (default: `TRUE`) — include biologically-informed reference lines in PDFs.
- `PLOT_PDF_ONE_TAG_ONLY` (main: `TRUE` when `RUN_ALL_TAGS=FALSE`; not used by `main_cpp_chunk.R`) — when `TRUE` produce PDFs only for `WHICH_TAG` (useful for single-tag runs).

Sensitivity & realism flags (available in `main_cpp.R`):
- `SENSITIVITY_MODE` (default: `"none"`) — Options: `"none"`, `"run"`, `"run+write"`, `"run+write+pdf"`. Controls whether sensitivity sweeps are executed and if results are written.
- `WRITE_OUTPUTS` (derived from `SENSITIVITY_MODE`) — internal flag to control writing sensitivity outputs when requested.
- `MAKE_ALL_SWEEPS_PDF` (derived) — whether to render all sweeps to PDF when `SENSITIVITY_MODE="run+write+pdf"`.
- `RUN_REALISM_REPORT` (default: `FALSE`) — when `TRUE` the script will generate a realism report for a representative species.
- `RUN_K_SWEEP_DEMO` (default: `FALSE`) — optional demo mode for k-sweep visualizations.

Note: the chunked runner (`main_cpp_chunk.R`) disables or comments out these options because per-chunk processing does not assemble a full `out` object for full-run sensitivity/realism processing.

Biological realism settings (defaults in script):
- `MAX_GROWTH_HARD_SOURCE = "fixed"`, `MAX_GROWTH_FIXED = 7.5` (`main_cpp.R`) / `5` (`main_cpp_chunk.R`)
- `MAX_SHRINK_HARD_SOURCE = "fixed"`, `MAX_SHRINK_FIXED = -0.5`
- `K_SHRINK_SOURCE = "fixed"`, `K_SHRINK_FIXED = 0`
- `K_GROWTH_SOURCE = "fixed"`, `K_GROWTH_FIXED = 0`
- `RECRUIT_MAX_SOURCE = "fixed"`, `RECRUIT_MAX_FIXED = (MAX_GROWTH_FIXED * 5) + 0.9999`
- `USE_MEASUREMENT_ERROR = TRUE`

Notes about chunking & downstream outputs:
- When chunking is active (`DP_CHUNK_SIZE > 0`), the script processes and writes chunk outputs incrementally and intentionally sets the in-memory `out` object to `NULL` to avoid excessive memory use.
- Because `out` is not assembled in memory for chunked runs, downstream steps that expect a combined `out` (e.g., writing a single combined RDS `stem_reconstruction_dp_global_rcpp.rds`, generating per-run PDFs from a combined `out`, or creating the realism report from `out`) will be skipped. Instead, you can work with the incremental CSV or per-chunk RDS files produced by the run.

---

## Chunking details & recommendations ✅

- When `RUN_ALL_TAGS=TRUE` and `DP_CHUNK_SIZE > 0`, groups are constructed from `unique(xrun[, .(Tag, species)])` and split into chunks of size `DP_CHUNK_SIZE`.
- For each chunk `ci`:
  - The script runs the DP in parallel over the groups in the chunk (each child sets `data.table::setDTthreads(1L)` to limit thread contention).
  - The chunk's results are combined into `out_chunk` and annotated with `DP_Chunk = ci`.
  - `maybe_add_posterior_bins(out_chunk)` is applied to add `DP_PosteriorBin` if requested.
- The chunk rows include a `run_out_dir` column set to the basename of the run `out_dir`.
- If `WRITE_DP_CSV=TRUE`, `out_chunk` is appended to `stem_reconstruction_dp_global_rcpp.csv` (the script writes the header only on the first write).
- If `WRITE_DP_RDS=TRUE`, `out_chunk` is saved as `stem_reconstruction_dp_global_rcpp_chunk_###.rds` (used as a resume marker).
- If `WRITE_DP_PDF=TRUE` and `WRITE_DP_PDF_PER_CHUNK=TRUE`, the script will attempt to generate a per-chunk PDF `stem_reconstruction_dp_global_rcpp_chunk_###.pdf` using `plot_tag_to_pdf()`; PDF generation errors are logged but will not abort the run.

Memory-saving recommendations:
- Prefer the chunking + incremental CSV approach for very large datasets (keeps peak RAM low).
- Keep `WRITE_DP_CSV=TRUE` so you get a single on-disk CSV that grows incrementally (append is memory-friendly).
- Keep per-chunk RDS files (`WRITE_DP_RDS=TRUE`) as reliable completion markers for resume and for reproducibility. These RDS files contain chunk results and can be merged later using `data.table::rbindlist(lapply(chunk_files, readRDS), use.names=TRUE, fill=TRUE)` on a machine with enough RAM or processed in streaming fashion.
- Avoid assembling a full `out` in memory if your dataset is large — the script intentionally does `out <- NULL` when chunking to prevent a memory spike.

---

## Where `maybe_add_posterior_bins` is applied

- If the script runs in **non-chunked** mode (i.e., `RUN_ALL_TAGS=TRUE` but `DP_CHUNK_SIZE <= 0` or chunking disabled), the full `out` is assembled in memory and `maybe_add_posterior_bins(out)` is applied immediately after assembly.
- If **chunking** is active, `maybe_add_posterior_bins()` is applied to each `out_chunk` before the chunk is written to disk.

This ensures posterior-bin computation is done while memory per-chunk is small and avoids double-processing.

Runner integration

- `main_cpp.R` and `main_cpp_chunk.R` are the primary entrypoints and are designed to be invoked directly with `Rscript`. External orchestrators can build CLI flags using the canonical names in `CLI_REFERENCE` (defined in each script); keys are case-insensitive and can use `-` or `_`.

## Post-reconstruction assertions

`run_dp_one_group()` (in both `main_cpp.R` and `main_cpp_chunk.R`) includes a **TrueStemID preservation assertion** that runs after every group, regardless of whether the DP or probabilistic pathway was used. The check verifies that every row carrying a non-NA `TrueStemID` has `ReconstructedStemID == TrueStemID`. Violations are logged as `WARN`-level messages via `log_msg()`, which writes to both stderr and `run_log.txt`.

**Warning messages from the probabilistic matcher** (sample-level repair counts, ME cumulative-shrinkage breaks, safety-net repair warnings) are emitted via `message()` on stderr and are also printed to stdout via `cat()` when `DP_VERBOSE=TRUE`. To capture all warnings in a log file, redirect both streams: `Rscript ... > log.txt 2>&1`.


---

## BCI Debug Driver (`main_cpp_bci.R`)

`main_cpp_bci.R` is a lightweight wrapper around `main_cpp.R` designed for debugging single tags from the BCI (Barro Colorado Island) multi-stem census dataset. It:

1. Sources `main_cpp.R` to load all helper functions and infrastructure modules.
2. Loads `bci_data/bci_multistem_xrun_debug.rds` (an RDS file, not tracked by git) as input.
3. Maps BCI-specific column names to the standard DP schema (e.g., `stemID` → `TrueStemID`, `dbh` → `DBH`, `censusID` → `CensusID`).
4. Runs the full DP pipeline for a single tag (default or `--WHICH_TAG=<id>`).

### BCI-specific defaults

- Input: `bci_data/bci_multistem_xrun_debug.rds` (loaded via `readRDS()`)
- Output: written to `dp_global/output/` under a BCI-specific run directory
- Uses `withr::with_dir()` for bundle sourcing compatibility
- `ANCHOR_START_CENSUS` defaults per the script configuration

### Example invocations

```bash
# Default tag
Rscript dp_global/scripts/main_cpp_bci.R

# Specific tag
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=123375

# Quiet mode
Rscript dp_global/scripts/main_cpp_bci.R --WHICH_TAG=187064 --DP_VERBOSE=FALSE
```

### Prerequisites

- R packages: `data.table`, `here`, `withr` (in addition to standard DP prerequisites)
- `bci_data/bci_multistem_xrun_debug.rds` must exist (not tracked by git; prepared separately)

---

## Example invocations

- Run a single tag interactively:

```bash
Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WHICH_TAG=20
```

- Run all tags with chunking (use the chunked driver):

Note: `dp_global/scripts/main_cpp_chunk.R` is the chunked driver intended for large runs. You can edit configuration variables at the top of the file or pass overrides via CLI flags (e.g., `--DP_CHUNK_SIZE=7`) when invoking it with `Rscript`. Run it directly via:

```
Rscript dp_global/scripts/main_cpp_chunk.R
```

- Run all tags without chunking (not recommended for very large datasets):

```
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=0
```

---

## Post-run: combining chunk files or troubleshooting

- By default, the merge helpers write to `stem_reconstruction_dp_global_rcpp_merged.csv` in the run `out_dir`. They locate chunk files matching `stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds` or `stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.feather` (you can use `prefer = 'feather'` with `merge_chunks_to_csv()` to prefer Feather sources).

To merge per-chunk RDS files into a single RDS (only do this on a machine with enough RAM):

```r
library(data.table)
chunk_files <- list.files("<out_dir>", pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
all <- rbindlist(lapply(chunk_files, readRDS), use.names = TRUE, fill = TRUE)
```

- If `stem_reconstruction_dp_global_rcpp.csv` already exists when you resume with `DP_CHUNK_RESUME=TRUE`, the script will append new chunk rows and skip chunks whose `_done.txt` completion marker exists. Use `--OUT_DIR_OVERRIDE=<path>` to point at the existing output directory.

---

## Basal Area Uncertainty (`basal_area_uncertainty.R`)

Post-processing script that uses posterior path samples from a completed run to quantify how identity uncertainty propagates into basal area (BA) estimates.

### Usage

```bash
Rscript dp_global/scripts/basal_area_uncertainty.R \
  --RUN_DIR=dp_global/output/<run_dir>
```

The script reads:
- `<RUN_DIR>/stem_reconstruction_dp_global_rcpp.csv` (main reconstruction)
- `<RUN_DIR>/posteriors/tag_*_posterior_samples__paths.csv` (posterior paths)

### Outputs

Three CSV files written to `<RUN_DIR>/`:

| File | Rows | Description |
|------|------|-------------|
| `basal_area_uncertainty_tag.csv` | Tag × Census | Total BA with posterior mean, SD, 95% CI |
| `basal_area_uncertainty_stem.csv` | Stem × Census | Per-stem BA posterior: mean, SD, median, 95% CI, MAP |
| `basal_area_uncertainty_growth.csv` | Stem × Census pair | BA growth rate ($\Delta$BA/$\Delta t$) posterior |

### Key insight

Tag-level total BA per census is **invariant** to identity assignment — the same DBH values are summed regardless of which stem identity each observation receives. All uncertainty is at the **per-stem level**: different identity assignments allocate different DBH values to each stem, producing different BA trajectories and growth rates.

---