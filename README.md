# STEM_IDENTIFICATION_TEST

## Overview
A collection of scripts and R modules for running and testing the DP stem-identification workflow.

Quick pointers:

- See `dp_global/README.md` for details on the DP implementation and R driver. 🔧
- Run `bin/run_dp_future.R` to execute experiments on a single machine (this is the recommended entrypoint). For single-config runs, use `bin/run_dp_future_single.R` or run `run_dp_future.R` with `--workers 1`. Please consult `bin/README.md` for examples and options.
- See `bin/README.md` for a quick reference to the runners and usage notes.

Quickstart
---------

1. Install dependencies: R.
2. Quick smoke checks:

   - Dry-run the recommended concurrent runner (safe; prints commands and writes per-config dry logs):

     ```bash
     ./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
     ```

   - Or run the lightweight smoke test that verifies core functions load:

     ```bash
     make smoke
     # which runs the DRY_RUN above and also: Rscript dp_global/scripts/run_smoke.R
     ```

Notes
-----
- Output directories are written under `dp_global/output/` by default. If you're using OneDrive or similar sync services, consider excluding this directory from sync to avoid I/O delays and accidental commits. ⚠️

Conventions
-----------
- Top-level helper scripts intended for users live in `bin/` (e.g., `bin/run_dp_future.R` — the recommended entrypoint for running experiments). Module internals and helpers live under `dp_global/`.
- Directory names are lowercase snake_case (e.g., `data_simulation`, `dp_global`). Most scripts and documentation follow this convention.

Testing & CI
------------
- Run `make smoke` to perform a lightweight smoke test: it executes the serial runner in `--DRY_RUN` mode and runs `Rscript dp_global/scripts/run_smoke.R` to verify core functions load without performing heavy computations.

