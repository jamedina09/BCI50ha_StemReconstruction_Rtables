# Runner quick reference — `run_dp_future.R` and `run_dp_future_single.R`

**IMPORTANT: DRY_RUN first ✅**

- Always verify constructed commands and logging paths before running real experiments. Use the `-- --DRY_RUN` flag to print the fully sanitized Rscript invocations and write per-config dry logs. Examples:

  ```bash
  ./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
  ./bin/run_dp_future_single.R --workers 1 --cores-per-job 1 -- --DRY_RUN
  ```

Purpose

- Recommended entrypoints for running experiments on a single machine:
  - `bin/run_dp_future.R` — concurrent runner using `future` + `progressr` (preferred for multi-config runs).
  - `bin/run_dp_future_single.R` — helper to run a fixed configuration set with the same concurrent/future orchestration semantics (useful for testing or single-config experiments).
- For single-config serial execution, use `bin/run_dp_future_single.R` or run `run_dp_future.R` with `--workers 1` for equivalent behavior and consistent logging.

Quick start

1. Dry-run (fast; no heavy computation):

   ```bash
   ./bin/run_dp_future.R --workers 2 --cores-per-job 4 --configs "fixed data_hard" -- --DRY_RUN
   ```

2. Real run (executes configs concurrently):

   ```bash
   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_soft"
   ```

3. Single-config runner (dry-run or real):

   ```bash
   ./bin/run_dp_future_single.R --workers 1 --cores-per-job 4 -- --DRY_RUN
   ```

Flags & options

- `--config=<name>`        : run a single config (e.g., `--config=fixed`) (used by serial helpers; for `run_dp_future.R` use `--configs`)
- `--DRY_RUN`               : prints Rscript invocations instead of executing them (pass after `--` to forward to `main_cpp.R`)
- `--dp_max_states=N`       : set the DP enumerator maximum states per census (default 40000). Use a larger value for more complex tags, but be mindful of memory.
- Any other `--NAME=VALUE` args are passed through to the R driver (`dp_global/scripts/main_cpp.R`).

Examples

- Dry-run a single config (via `run_dp_future_single.R`):

  ```bash
  ./bin/run_dp_future_single.R --workers 1 --cores-per-job 4 -- --DRY_RUN
  ```

Notes & troubleshooting

- The runners respect an exported `BATCH_TS` (useful to group outputs by timestamp).
- Choosing `--dp_max_states`: the appropriate maximum depends on dataset complexity. Use `dp_global/R/test_complexity_estimator.R` to evaluate per-tag complexity and guide the choice. Example:

  ```bash
  Rscript dp_global/R/test_complexity_estimator.R
  ```

- Prefer passing `interval_years` explicitly to the DP or including a per-row interval column such as `Bio_IntervalYears` in your input. To enable automatic per-pair interval detection in DP functions set `interval_years = NULL` (functions will search for candidate columns when `NULL`).

Concurrent runs with `future` + `progressr` (working)

- Purpose: Run multiple named experiment configs concurrently in a robust, observable way.
- Script: `bin/run_dp_future.R` (R script using `future`, `future.apply` and `progressr`).

Installation (if needed)

```r
install.packages(c("future", "future.apply", "progressr"))
```

Usage examples

- Dry-run (safe; prints commands and writes per-config dry logs):

  ```bash
  ./bin/run_dp_future.R --workers 3 --cores-per-job 5 --configs "fixed data_hard" -- --DRY_RUN
  ```

- Real run (executes the configs concurrently):

  ```bash
  ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_soft"
  ```

Options

- `--workers N` (or `-j N`): number of concurrent tasks (future workers). Default: 3
- `--cores-per-job N`: assumed CPU cores each task uses (used for safety check). Default: 5
- `--configs "c1 c2 ..."`: list of configs to run (space- or comma-separated)
- `--joblog FILE`: output CSV summary (default `parallel_future.log`)
- `--force`: override oversubscription safety check
- Any args after `--` are forwarded to `dp_global/scripts/main_cpp.R` (e.g., `--DRY_RUN`, `--SENSITIVITY_MODE=none`)

Behavior & outputs

- Oversubscription safety: script checks `workers * cores_per_job <= available logical CPUs` and will abort unless `--force` used.
- Per-config logs: `tests/parallel_future_logs/<config>.log` (created if not present)
- Joblog CSV: (default) `parallel_future.log` summarizing start/end, exit status, and log file path
- Progress: displayed with `progressr` (text progress bar by default)

Caveats & tips

- I/O caution: on OneDrive or other networked filesystems concurrent writes can be slow or cause conflicts. Use lower concurrency or write to a local fast directory and move outputs after completion.
- Oversubscription: ensure `workers * cores-per-job <= physical cores` (or use `--force`) to avoid thrashing.
- DRY_RUN: use first to verify the exact commands that will be run before executing.
