# dp_global/scripts

This document describes the `main_cpp.R` driver script (located at `dp_global/scripts/main_cpp.R`), what it does, how to run it, and the meaning and defaults of every major parameter and CLI flag. It also documents the chunking/resume behavior implemented to avoid out-of-memory issues when running reconstruction on large datasets.

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

- By default the script runs in a single-tag mode (`RUN_ALL_TAGS = FALSE`). To run all tags, set `--RUN_ALL_TAGS=TRUE`.
- Running `RUN_ALL_TAGS=TRUE` on very large datasets can cause out-of-memory (OOM) errors. To mitigate this, the script supports a **chunking** mode that processes groups (Tag + species) in chunks of `DP_CHUNK_SIZE` and writes chunk outputs to disk incrementally.
- When chunking, the script will NOT assemble a full in-memory `out` table (to avoid OOM). Instead it:
  - Appends each chunk's rows to the main CSV (`stem_reconstruction_dp_global_rcpp.csv`) so a single CSV on disk is produced incrementally, and
  - Writes a per-chunk RDS file `stem_reconstruction_dp_global_rcpp_chunk_###.rds` for reproducibility and resume markers.
- Posterior bin computation (`add_dp_posterior_bins`) is applied per-chunk in the chunking loop and for the non-chunking path it is applied after the full `out` is assembled. This avoids double-processing.

---

## Key CLI flags & variables (with defaults) 🔧

Below are the canonical keys exposed via the CLI and their defaults in the script. Keys are case-insensitive and accept `-`/`_` separators.

- `INPUT_FILE` (default: `data_simulation/data/simulated_data_1.csv`) — path to input CSV containing tree observations.
- `FORCE_ONE_SPECIES_PARAMETERS` (default: `TRUE`) — if `TRUE`, all trees are treated as one species using `FORCED_SPECIES_LABEL`.
- `FORCED_SPECIES_LABEL` (default: `"all"`) — species label used when forcing one species.
- `SPECIES_COL` (default: `NULL`) — explicit species column name to use in input; otherwise the script attempts to find a candidate.

DP / reconstruction options:
- `DP_MODE` (default: `"marginals+bins"`) — DP mode. Allowed: `none`, `marginals`, `marginals+bins`, `map`.
- `which_tag` (default: `19`) — Tag to process when not running all tags.
- `anchor_start_census` (default: `7`) — census used as anchor start.
- `DP_VERBOSE` (default: `TRUE`) — verbose DP logging.
- `DP_POSTERIOR_TOP_K` (default: `2`) — top-k posterior reconstructions to track.
- `dp_max_tracks` (default: `NULL`) — optionally force max tracks; if `NULL` auto-computed.
- `dp_max_states` (default: `40000`) — limit DP max states to control memory/CPU.
- `dp_slack_tracks` (default: `1`) — slack tracks allowed for DP.
- `dp_slack_require_anchor_recruitable` (default: `TRUE`) — require anchor to be recruitable before granting slack.
- `dp_slack_require_anchor_eps` (default: `1e-6`) — numerical tolerance for anchor recruitability checks.

Posterior sampling:
- `POSTERIOR_SAMPLES` (default: `200`) — number of full-path reconstructions to draw; `0` disables.
- `POSTERIOR_SAMPLES_FORMAT` (default: `"csv"`) — `rds`,`feather`, or `csv`.
- `POSTERIOR_SAMPLES_PATH` (default: `NULL`) — path for posterior sample files; default to `<out_dir>/posteriors` when sampling.
- `POSTERIOR_SAMPLE_SEED` (default: `NULL`) — seed for reproducible posterior sampling.

Output controls:
- `WRITE_DP_CSV` (default: `TRUE`) — write incremental/combined CSV output.
- `WRITE_DP_RDS` (default: `TRUE`) — write per-chunk RDS (or combined RDS for non-chunk runs).
- `WRITE_DP_FEATHER` (default: `FALSE`) — write per-chunk Feather files (`.feather`) for faster IO. Requires the **arrow** R package.
- `WRITE_DP_PDF` (default: `TRUE`) — write PDF visualizations (`plot_tag_to_pdf`).
- `DP_PDF_INCLUDE_REFERENCE` (default: `TRUE`) — include reference lines in PDF output.

Parallel & chunking controls:
- `RUN_ALL_TAGS` (default: `FALSE`) — run all tags instead of a single `which_tag`.
- `MANUAL_CORES` (default: `TRUE`) — set to `TRUE` to use `MANUAL_CORES_VALUE`.
- `MANUAL_CORES_VALUE` (default: `1`) — number of cores to use when `MANUAL_CORES=TRUE`.
Note: There is a dedicated chunked driver script `dp_global/scripts/main_cpp_chunk.R` that implements chunked processing. It is a standalone script intended to be run manually (no CLI argument parsing); edit parameters at the top of `main_cpp_chunk.R` to change behavior, then run it via:

- `Rscript dp_global/scripts/main_cpp_chunk.R`

Helpful post-run utilities:
- Merge chunk RDS/Feather files into a single CSV (run this in R or source the script and call the helper):

```r
# from R in project root
source("dp_global/scripts/main_cpp_chunk.R")
merge_chunks_to_csv("dp_global/output/<your_run_dir>")
```

This streams each chunk file to a single CSV to avoid loading the full dataset into memory.

For backwards compatibility, the original `main_cpp.R` runs the non-chunked workflow (single-tag or parallelized tags) and does not perform per-chunk writing.

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

PDF & plotting controls:
- `WRITE_DP_PDF` (default: `TRUE`) — control whether PDFs are generated via `plot_tag_to_pdf()`.
- `DP_PDF_INCLUDE_REFERENCE` (default: `TRUE`) — include biologically-informed reference lines in PDFs.
- `PLOT_PDF_ONE_TAG_ONLY` (default: `TRUE` when `RUN_ALL_TAGS=FALSE`) — when `TRUE` produce PDFs only for `which_tag` (useful for single-tag runs).

Sensitivity & realism flags:
- `SENSITIVITY_MODE` (default: `"none"`) — Options: `"none"`, `"run"`, `"run+write"`, `"run+write+pdf"`. Controls whether sensitivity sweeps are executed and if results are written.
- `WRITE_OUTPUTS` (derived from `SENSITIVITY_MODE`) — internal flag to control writing sensitivity outputs when requested.
- `MAKE_ALL_SWEEPS_PDF` (derived) — whether to render all sweeps to PDF when `SENSITIVITY_MODE="run+write+pdf"`.
- `RUN_REALISM_REPORT` (default: `FALSE`) — when `TRUE` the script will generate a realism report for a representative species.
- `RUN_K_SWEEP_DEMO` (default: `FALSE`) — optional demo mode for k-sweep visualizations.

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

- Run one tag (interactive):

```
Rscript dp_global/scripts/main_cpp.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WHICH_TAG=19
```

- Run all tags with chunking (use the chunked driver):

Note: `dp_global/scripts/main_cpp_chunk.R` is a **standalone, manual** chunk driver that intentionally does not parse the main CLI. Edit the configuration block at the top of `main_cpp_chunk.R` (e.g., `DP_CHUNK_SIZE`, `DP_CHUNK_START`, `DP_CHUNK_END`, `WRITE_DP_PDF_PER_CHUNK`, `POSTERIOR_SAMPLES`, `MANUAL_CORES_VALUE`) and then run it directly via:

```
Rscript dp_global/scripts/main_cpp_chunk.R
```



- Run all tags without chunking (not recommended for very large datasets):

```
Rscript dp_global/scripts/main_cpp.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=0
```

---

## Post-run: combining chunk files or troubleshooting

- To merge per-chunk RDS files into a single RDS (only do this on a machine with enough RAM):

```r
library(data.table)
chunk_files <- list.files("<out_dir>", pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
all <- rbindlist(lapply(chunk_files, readRDS), use.names = TRUE, fill = TRUE)
```

- If `stem_reconstruction_dp_global_rcpp.csv` already exists when you resume with `DP_CHUNK_RESUME=TRUE`, the script will append new chunk rows and skip chunks with existing RDS files.

---