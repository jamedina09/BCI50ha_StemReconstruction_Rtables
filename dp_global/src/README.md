# Rcpp Acceleration for DP Global Stem Tracking

C++ implementation and R wrapper for Rcpp-accelerated transition-cost computations used by the DP workflow.

## Files

- `src/transition_cost_rcpp.cpp`: C++ implementation of transition cost computation
- `R/transition_cost_rcpp.R`: R wrapper function for the C++ code

> Note: Validation of the C++ implementation is performed via the test suite in `dp_global/dev/` (see below).

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

### Quickstart & API

Load and call from R (project root):

```r
library(Rcpp)
Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
source("dp_global/src/transition_cost_rcpp.R")

# Minimal example
defaults <- list(mu_const = 0, mu_gamma = 0, sigma0 = 1, sigma1 = 0,
                 max_shrink = -Inf, k_shrink = 0, max_growth = Inf,
                 max_growth_soft = Inf, k_growth = 0, use_measurement_error = FALSE,
                 meas_sd1_a = 0.0062, meas_sd1_b = 0.0904, meas_sd2 = 4.64,
                 meas_p_big = 0.05, h0 = 0, beta = 0,
                 recruit_meanlog = 0, recruit_sdlog = 1, recruit_max_dbh = 200,
                 recruit_lambda = 0, eps_tiebreak = 1e-6)

track_dbh_t <- c(10.0, NA, 20.0)
mat_tp1 <- matrix(c(12.0, 15.0, 18.0, 8.0, NA, 22.0), nrow = 2, byrow = TRUE)

costs <- do.call(transition_cost_tracks_bio_batch_rcpp, c(list(track_dbh_t = track_dbh_t, track_dbh_tp1 = mat_tp1, interval_years = 5), defaults))
print(costs)
```

See `R/transition_cost_rcpp.R` for detailed parameter descriptions and expected types.
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

**Notes:**
- Use the `dp_global/dev/` test suite to validate the C++ implementation and to exercise end-to-end behavior through `main_cpp.R` when needed.
- When running end-to-end checks, pass `--WRITE_DP_RDS=TRUE` to the run command to ensure `.rds` outputs are produced and inspected by R-based tests (preferred over parsing CSVs for validation).

If you edit the C++ implementation or the R wrapper, re-run the dev validation tests and the profiler to ensure correctness and to check performance.

## Troubleshooting

- If `Rcpp::sourceCpp()` fails on macOS, ensure Xcode command-line tools are installed (`xcode-select --install`) and that your `R` can find clang/clang++.
- If you encounter small numeric differences versus the R implementation, verify `eps_tiebreak` and floating-point tolerances.
- For reproducible CI, consider packaging the C++ code into an R package to avoid on-the-fly compilation variability.

## Implementation Details

The C++ version implements the same logic as the R version but with:

- Direct C++ loops instead of R vectorized operations
- Manual implementation of statistical functions (dnorm, dlnorm, log_sum_exp)
- Optimized memory access patterns
- Reduced function call overhead

All edge cases, measurement error models, biological constraints, and tie-breaking logic are preserved.