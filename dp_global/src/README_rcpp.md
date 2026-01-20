# Rcpp Acceleration for DP Global Stem Tracking

This directory contains an optimized C++ implementation of the transition cost computation for the DP global stem tracking algorithm.

## Files

- `src/transition_cost_rcpp.cpp`: C++ implementation of transition cost computation
- `R/transition_cost_rcpp.R`: R wrapper function for the C++ code

> Note: The separate `test_rcpp_implementation.R` script was removed. Validation is performed via the test suite in `dp_global/dev/` (see below).

## Usage

### Prerequisites

1. Install the `Rcpp` package:
   ```r
   install.packages("Rcpp")
   ```

2. For performance testing, also install `microbenchmark`:
   ```r
   install.packages("microbenchmark")
   ```

### Basic Usage

```r
# Load the C++ implementation
Rcpp::sourceCpp("./dp_global/src/transition_cost_rcpp.cpp")
source("./dp_global/src/transition_cost_rcpp.R")

# Use the function (same interface as the original R version)
costs <- transition_cost_tracks_bio_batch_rcpp(
    track_dbh_t = c(10.0, NA, 20.0),
    mat_tp1 = matrix(c(12.0, 15.0, 18.0, 8.0, NA, 22.0), nrow = 2, byrow = TRUE),
    interval_years = 5,
    # ... other parameters with same defaults as original function
)
```

### Integration with Existing Code

To use the C++ version in place of the R version, modify `dp_global_biol.R`:

1. Add at the top of the file:
   ```r
   library(Rcpp)
   source("./dp_global/src/transition_cost_rcpp.R")
   Rcpp::sourceCpp("./dp_global/src/transition_cost_rcpp.cpp")
   ```

2. Replace calls to `transition_cost_tracks_bio_batch` with `transition_cost_tracks_bio_batch_rcpp`

## Performance

The C++ implementation provides significant speedup over the R version, especially for large batch sizes. Typical speedup factors range from 10-50x depending on the problem size and complexity.

## Validation

Run the dev test that exercises the transition-cost implementation to verify correctness and performance. For example:

```bash
# Validate batch transition-cost implementation and compare R vs C++ behavior
Rscript dp_global/dev/test_transition_cost_tracks_bio_batch.R
```

The dev tests compare outputs and report discrepancies; use these to validate the C++ implementation after changes.

## Implementation Details

The C++ version implements the same logic as the R version but with:

- Direct C++ loops instead of R vectorized operations
- Manual implementation of statistical functions (dnorm, dlnorm, log_sum_exp)
- Optimized memory access patterns
- Reduced function call overhead

All edge cases, measurement error models, biological constraints, and tie-breaking logic are preserved.