# Runner quick reference — preferred: `run_dp_future.R`

IMPORTANT: DRY_RUN first

- Before executing real runs, always use the `--DRY_RUN` flag to verify the constructed commands and logging paths. Example:

  ```bash
  ./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
  ```

Purpose

- The recommended entrypoint for running experiments on a single machine is `bin/run_dp_future.R` (concurrent runner using `future` + `progressr`).
- Legacy serial helpers (`bin/run_dp_full_cpp.sh`, `bin/run_dp_full.sh`) still exist for running configurations one-at-a-time, but `run_dp_future.R` is preferred for reproducible concurrent runs and better logging.

Quick start

1. Make executable (if needed):

   chmod +x bin/run_dp_full_cpp.sh

2. Dry-run (fast, no DP work):

   bin/run_dp_full_cpp.sh --config=fixed --MANUAL_CORES_VALUE=16 --SENSITIVITY_MODE=none --DP_MODE=none --DRY_RUN

3. Run all configurations (serial):

   bin/run_dp_full_cpp.sh

Flags & options

- --config=<name>        : run a single config (e.g., `--config=fixed`)
- --MANUAL_CORES_VALUE=N : set per-job cores for resource heuristics (default 15)
- --DRY_RUN              : prints Rscript invocations instead of executing them
- --dp_max_states=N      : set the DP enumerator maximum states per census (default 40000). Use a larger value for more complex tags, but be mindful of memory; see Notes below for how to estimate a dataset-appropriate value.
- Any other `--NAME=VALUE` args are passed through to the R driver.

Examples

- Run a single config with 15 cores and DRY_RUN:

  bin/run_dp_full_cpp.sh --config=fixed --MANUAL_CORES_VALUE=16 --SENSITIVITY_MODE=none --DP_MODE=none --DRY_RUN

- Run all configs (serial):

  bin/run_dp_full_cpp.sh

Notes & troubleshooting

- The script still respects an exported `BATCH_TS` (if you want to group outputs by timestamp).
- If you need parallel execution, run concurrent invocations manually and manage `BATCH_TS` and resource coordination yourself; the project no longer provides the GNU-parallel wrapper.
- The `make smoke` target now runs the DRY_RUN and also executes `Rscript dp_global/scripts/run_smoke.R` to verify core functions load.
- Choosing `--dp_max_states`: the appropriate maximum states depends on dataset complexity. Use the estimator script `dp_global/R/test_complexity_estimator.R` to get per-tag complexity summaries and a recommended max-states guidance (it prints ranked tags and exports a CSV summary to `data_simulation/data/`). Run it with:

  ```bash
  Rscript dp_global/R/test_complexity_estimator.R
  ```

  and inspect the exported `report_run_<dataset>.csv` for per-tag `max_states` and other helpful diagnostics.

Concurrent runs with `future` + `progressr` (not working yet)

- Purpose: Run multiple `bin/run_dp_full_cpp.sh` configs concurrently in a robust, observable way.
- Script: `bin/run_dp_future.R` (R script using `future`, `future.apply` and `progressr`.)

Installation (if needed)

```r
install.packages(c("future", "future.apply", "progressr"))
```

Usage examples

- Dry-run (safe; prints commands and writes per-config dry logs):

  ./bin/run_dp_future.R --workers 3 --cores-per-job 5 --configs "fixed data_hard" -- --DRY_RUN

- Real run (executes the configs concurrently):

  ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_soft"

- Full run

  ./bin/run_dp_full_cpp.sh --MANUAL_CORES_VALUE=16

Options

- `--workers N` (or `-j N`): number of concurrent tasks (future workers). Default: 3
- `--cores-per-job N`: assumed CPU cores each task uses (used for safety check). Default: 5
- `--configs "c1 c2 ..."`: list of configs to run (space- or comma-separated)
- `--joblog FILE`: output CSV summary (default `parallel_future.log`)
- `--force`: override oversubscription safety check
- Any args after `--` are forwarded to `bin/run_dp_full_cpp.sh` (e.g., `--DRY_RUN`, `--SENSITIVITY_MODE=none`)

Behavior & outputs

- Oversubscription safety: script checks `workers * cores_per_job <= available logical CPUs` and will abort unless `--force` used.
- Per-config logs: `tests/parallel_future_logs/<config>.log` (created if not present)
- Joblog CSV: (default) `parallel_future.log` summarizing start/end, exit status, and log file path
- Progress: displayed with `progressr` (text progress bar by default)

Caveats & tips

- I/O caution: on OneDrive or other networked filesystems concurrent writes can be slow or cause conflicts. Use lower concurrency or write to a fast local temp directory and move outputs after completion.
- Oversubscription: ensure `workers * cores-per-job <= physical cores` (or use `--force`) to avoid thrashing.
- DRY_RUN: use first to verify the exact commands that will be run before executing.

Example DRY_RUN with forwarded args:

  ./bin/run_dp_future.R --workers 2 --cores-per-job 8 --configs "fixed data_soft" -- --SENSITIVITY_MODE=none --DP_MODE=none --DRY_RUN
