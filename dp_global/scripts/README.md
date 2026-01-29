# dp_global/scripts

This document describes the current behavior of the `dp_global` driver scripts:
- `dp_global/scripts/main_cpp.R` — the interactive/CLI driver for single-tag or targeted runs, and
- `dp_global/scripts/main_cpp_chunk.R` — the chunked driver optimized for large runs.

Both scripts accept command-line overrides of defaults using `--KEY=VALUE` flags. Keys are case-insensitive and may use `-` or `_` as separators.

---

## Overview ✨

`main_cpp.R` is the central driver for the `dp_global` pipeline. It:

- Loads input tree census data (CSV via `input_file`).
- Ensures species information is present and optionally forces a single species label.
- Estimates biological parameters per species (via `estimate_bio_pars`).
- Runs Dynamic Programming (DP) reconstruction logic (via `run_dp_one_group` which calls the Rcpp-based DP functions).
- Optionally runs sensitivity sweeps and realism checks.
- Writes outputs (CSV, RDS, PDF) to an automatically created output directory.

You can run the script directly via Rscript or by using the provided orchestrator wrapper:

- Direct: `Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=... --RUN_ALL_TAGS=TRUE` 
- Or via the orchestrator: `bin/run_dp_future_single.R` (this wrapper constructs CLI flags and runs the driver)

(When running in the workspace, a canonical example is `Rscript scripts/main_cpp.R` from the `dp_global` folder.)

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
- `WHICH_TAG` — used for single-tag runs (relevant to `main_cpp.R`); the chunked runner processes groups (`Tag`, `species`) and does not rely on `WHICH_TAG`.
- `ANCHOR_START_CENSUS` — default: `7L`.
- `ALLOW_PROVISIONAL_DP_ANCHOR` — default: `TRUE` — when `TRUE` the DP can assign provisional anchor IDs at the last observed DBH census if the requested anchor census lacks `TrueStemID` but has DBH; set to `FALSE` to require an explicit anchored census or to fall back to the igraph matcher.
- `DP_VERBOSE` — default: `TRUE`.
- `DP_POSTERIOR_TOP_K` — default: `2L` - top=k posterior reconstructions to track. 
- `DP_MAX_TRACKS` — default: `NULL` (auto-computed per data when `NULL`) — optionally force max tracks; if `NULL` auto-computed.
- `DP_MAX_STATES` — default: `40000L` - limit DP max states to control memory/CPU.
- `DP_SLACK_TRACKS` — default: `1L` - slack (additional) tracks allowed for DP.
- `DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE` — default: `TRUE` - require anchor to be recruitable (if DBH less than max recruitment size) before granting slack.
- `DP_SLACK_REQUIRE_ANCHOR_EPS` — default: `1e-6`.

**Anchor scoping and post-anchor preservation:** If observations exist after the requested `ANCHOR_START_CENSUS`, the DP is scoped to censuses <= `ANCHOR_START_CENSUS` and post-anchor rows are preserved and appended to the output. Post-anchor rows with non-NA `DBH` and a `TrueStemID` that was used by the DP are set to `ReconstructedStemID = TrueStemID` and `ReconstructionMethod = "given"`. Remaining post-anchor rows without DP assignments receive `ReconstructionMethod = "none_after_anchor"`. If scoping removes all pre-anchor observations, the original rows are returned and observed `TrueStemID` values are treated as `given` while other rows are labeled `none_after_anchor`.

**Provisional anchor behavior:** When a requested anchor census lacks `TrueStemID` but contains DBH observations and `ALLOW_PROVISIONAL_DP_ANCHOR=TRUE`, the DP will assign provisional anchor IDs at the last-observed DBH census and mark those anchor rows with `ReconstructionMethod = "provisional_dp"`. If the DP cannot anchor or a fallback is used, an igraph fallback can assign provisional IDs and mark them with `ReconstructionMethod = "provisional_igraph"`.

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
- `DP_PDF_INCLUDE_REFERENCE` — `main_cpp.R` default: `TRUE`; `main_cpp_chunk.R` default: `FALSE` — include biologically-informed reference lines in PDFs (useful with simulated data).
- `WRITE_DP_PDF_PER_CHUNK` — default in `main_cpp_chunk.R`: `TRUE` (controls per-chunk PDFs).

Parallel & chunking controls (chunked runner specific):
- `DP_CHUNK_SIZE` — default in `main_cpp_chunk.R`: `7L` (set `<= 0` to disable chunking behavior when applicable).
- `DP_CHUNK_RESUME` — default: `TRUE` (skip chunks with existing RDS files) - allows to stop runs and continue later.
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
- `PLOT_PDF_ONE_TAG_ONLY` (main: `TRUE` when `RUN_ALL_TAGS=FALSE`; not used by `main_cpp_chunk.R`) — when `TRUE` produce PDFs only for `which_tag` (useful for single-tag runs).

Sensitivity & realism flags (available in `main_cpp.R`):
- `SENSITIVITY_MODE` (default: `"none"`) — Options: `"none"`, `"run"`, `"run+write"`, `"run+write+pdf"`. Controls whether sensitivity sweeps are executed and if results are written.
- `WRITE_OUTPUTS` (derived from `SENSITIVITY_MODE`) — internal flag to control writing sensitivity outputs when requested.
- `MAKE_ALL_SWEEPS_PDF` (derived) — whether to render all sweeps to PDF when `SENSITIVITY_MODE="run+write+pdf"`.
- `RUN_REALISM_REPORT` (default: `FALSE`) — when `TRUE` the script will generate a realism report for a representative species.
- `RUN_K_SWEEP_DEMO` (default: `FALSE`) — optional demo mode for k-sweep visualizations.

Note: the chunked runner (`main_cpp_chunk.R`) disables or comments out these options because per-chunk processing does not assemble a full `out` object for full-run sensitivity/realism processing.

Biological realism settings (defaults in script):
- `MAX_GROWTH_HARD_SOURCE = "data"`, `MAX_GROWTH_FIXED = 7.5`
- `MAX_SHRINK_HARD_SOURCE = "data"`, `MAX_SHRINK_FIXED = -0.5`
- `K_SHRINK_SOURCE = "data"`, `K_SHRINK_FIXED = 0`
- `K_GROWTH_SOURCE = "data"`, `K_GROWTH_FIXED = 0`
- `RECRUIT_MAX_SOURCE = "fixed"`, `RECRUIT_MAX_FIXED = 5`

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

- When an orchestrator (`bin/run_dp_future_single.R` / `bin/run_dp_future.R`) runs DP, it attempts to detect the run `out_dir` by matching `BATCH_TS` and `CONFIG_NAME`. If the `out_dir` is located, the orchestrator will include the main run's `run_log.txt` in the per-config worker log and record the `main_out_dir` in the joblog CSV for easier discovery of artifacts and diagnostics.

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

- If `stem_reconstruction_dp_global_rcpp.csv` already exists when you resume with `DP_CHUNK_RESUME=TRUE`, the script will append new chunk rows and skip chunks with existing RDS files.

---