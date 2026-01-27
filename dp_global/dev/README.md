# dp_global dev

Development tools (tests, benchmarks, profiling) for the dp_global stem-matching workflow.

Notes:
- For per-pair interval handling in tests and examples, prefer `interval_years = NULL` (this enables detection of per-row interval columns like `Bio_IntervalYears`) or pass an explicit scalar `interval_years` when needed.
- Use `--WRITE_DP_RDS=TRUE` to produce `.rds` artifacts alongside CSV/PDF outputs for R-native validation in tests.

## Quickstart

Run the main tests (from project root):

```bash
Rscript dp_global/dev/test_transition_cost_tracks_bio_batch.R
Rscript dp_global/dev/test_match_stems_dp_global_backward_marginals_batch.R
```

## Prerequisites

R (packages: `data.table`, `Rcpp`); access to project root and synthetic data in `data_simulation/data/` for profiling.

## Scripts Overview

### Regression Tests

These scripts validate the correctness and functionality of the dp_global implementation.

- **`test_transition_cost_tracks_bio_batch.R`**: 
  - **Purpose**: Validates the batch transition cost computation against scalar (single-transition) results.
  - **Tests**: Ensures the biological cost model correctly handles multiple candidate transitions. Compares R batch implementation against C++ batch (if `Rcpp` available).
  - **Inputs**: Synthetic track data with various DBH transitions (growth, shrinkage, mortality, recruitment).
  - **Outputs**: Success message if all comparisons pass; errors if discrepancies found.
  - **Usage**: `Rscript dp_global/dev/test_transition_cost_tracks_bio_batch.R`

- **`test_match_stems_dp_global_backward_marginals_batch.R`**: 
  - **Purpose**: Validates that the full batch DP algorithm executes without errors on synthetic data.

- **`test_estimate_bio_pars_intervals.R`**:
  - **Purpose**: Tests `estimate_bio_pars()` interval handling (scalar vs per-row/per-pair interval columns), verifies inference and diagnostic outputs, and writes per-pair interval diagnostics for inspection.
  - **Usage**: `Rscript dp_global/dev/test_estimate_bio_pars_intervals.R`
  - **Tests**: Runs the DP on minimal synthetic tree data. Tests both R and C++ implementations (if available).
  - **Inputs**: Synthetic `tree_data` with multiple censuses, anchor census, and populated biological parameters.
  - **Outputs**: Success messages for R and C++ versions; errors if the DP fails to run.
  - **Usage**: `Rscript dp_global/dev/test_match_stems_dp_global_backward_marginals_batch.R`

### Benchmarking and Profiling

- **`profiling_code.R`**:
  - **Purpose**: Profile the C++-accelerated DP pipeline to locate hotspots (CPP-only profiler).
  - **Inputs**: Synthetic dataset from `data_simulation/data/simulated_data_1.csv`.
  - **Outputs**: Rprof profiling file `dp_global_cpp.prof` written to `dp_global/dev`.
  - **Usage**:
    - Run the profiler: `Rscript --vanilla dp_global/dev/profiling_code.R`
    - Disable tie-breaking for performance tests: `DP_EPS_TIEBREAK=0 Rscript --vanilla dp_global/dev/profiling_code.R`

## Output Files

- **`dp_global_cpp.prof`**: Profiling output for the C++-only profiler. This file contains R profiling data that can be analyzed with `summaryRprof()` or visualization tools.

## Usage

All scripts should be run from the project root directory to ensure correct relative paths.

### Running Tests

Execute the regression tests to verify functionality:

```bash
# Test transition cost batch computation
Rscript dp_global/dev/test_transition_cost_tracks_bio_batch.R

# Test full DP algorithm
Rscript dp_global/dev/test_match_stems_dp_global_backward_marginals_batch.R
```

### Running Benchmarks

Use the profiling script to compare performance:

```bash
# Benchmark all variants
Rscript --vanilla dp_global/dev/profiling_code.R

# Profile specific variant
PROFILE_VARIANT=cpp Rscript --vanilla dp_global/dev/profiling_code.R
```

## Notes

- C++ features require the `Rcpp` package and on-the-fly compilation of `src/transition_cost_rcpp.cpp`.
- Benchmarks measure execution time; profiling identifies bottlenecks.
- Synthetic data ensures reproducible testing without relying on real datasets.
- If tests fail, check for missing dependencies or incorrect column names in synthetic data.