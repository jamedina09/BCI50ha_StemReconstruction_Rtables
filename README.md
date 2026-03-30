# STEM_IDENTIFICATION_TEST

## Overview
A collection of scripts and R modules for running and testing the DP stem-identification workflow.

Quick pointers:

- See `dp_global/README.md` for details on the DP implementation and R modules.
- Run `dp_global/scripts/main_cpp.R` for single-tag or small-dataset runs. For large datasets, use `dp_global/scripts/main_cpp_chunk.R` which writes outputs incrementally and supports resume.
- See `dp_global/scripts/README.md` for a full CLI flag reference, chunking details, and example invocations.
- You can force the DP solver to skip particular growth forms by passing
  `--DP_FALLBACK_GROWTH_FORMS=<form1,form2>` (for example, `tree` or `fig`) to
  the driver; the flag accepts a comma‑ or semicolon‑separated list. Any tag
  containing a row with a matching `growth_form` value will be reconstructed
  with the igraph matcher.

## Prerequisites

R (packages: `data.table`, `igraph`).

## Quickstart

1. Install R and required packages (see `## Prerequisites`).

2. Run a single tag to verify the pipeline:

```bash
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=20
```

3. Run the smoke test to verify core functionality:

```bash
make smoke
```

Notes
-----
- Output directories are written under `dp_global/output/` by default. If you're using OneDrive or similar sync services, consider excluding this directory from sync to avoid I/O delays and accidental commits. ⚠️

Conventions
-----------
- Driver scripts live under `dp_global/scripts/` (`main_cpp.R` for single-tag/small runs; `main_cpp_chunk.R` for large chunked runs). Module internals and helpers live under `dp_global/R/`.
- Directory names are lowercase snake_case (e.g., `data_simulation`, `dp_global`). Most scripts and documentation follow this convention.

Testing & CI
------------
- Run `make smoke` to perform a lightweight smoke test: sources `dp_global/R/dp_global_main.R` to verify all R modules load without errors.

## TODO / Future improvements

Short actionable items to stabilize the DP and posterior-attachment workflow and make downstream usage easier:

- **Add a convenience alias `ObsRowID` → `obs_row_id` on final outputs** to simplify joins with posterior `paths.csv` (which use `ObsRowID`). Keep `obs_row_id` as the canonical internal name; create `ObsRowID` before export and add a unit test to assert equality.

- **Convert ad-hoc validation scripts into formal tests (testthat) and add to CI**:
  - Add tests for anchor-scoping semantics (anchor at first census, anchor > last observed census, post-anchor preservation, and provisional anchor behavior).
  - Add tests for posterior-path attachment (`attach_paths_to_output()`) covering `ObsRowID` mapping, alternative reconstruction formats, expand/aggregate correctness, and failure modes.
  - Add a regression test asserting the `DP_PruneInfo` attribute exists (defensive init) and that early-return branches do not error.

- **Document `obs_row_id`/`ObsRowID` expectations** in `dp_global/README.md` (how posterior reconstructions map back to output rows).

- **Other improvements / robustness**:
  - Ensure posterior sampling output (feather/csv/rds) preserves `ObsRowID` when requested.
  - Consider adding tests for pruning diagnostics and timing under realistic small examples.
