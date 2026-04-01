# Rcpp Acceleration for DP Global Stem Tracking

C++ implementation and R wrapper for Rcpp-accelerated transition-cost and phase-feasibility computations used by the DP workflow.

## Files

- `transition_cost_rcpp.cpp`: C++ implementation of transition cost computation and phase-feasibility batch checks
- `transition_cost_rcpp.R`: R wrapper function for the C++ code

## Usage

### Prerequisites

R (packages: `Rcpp`; optional for benchmarking: `microbenchmark`) 

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

## Performance

The C++ implementation is the sole backend for transition-cost computation in the DP workflow. It is substantially faster than an equivalent pure-R loop implementation would be, due to direct C++ iteration, manual statistical functions, and reduced function-call overhead.

## Validation

To validate the C++ implementation after changes, run an end-to-end check through `main_cpp_chunk.R`:

```bash
Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WRITE_DP_RDS=TRUE
```

Compare chunk RDS output against a saved reference to detect regressions. The `--WRITE_DP_RDS=TRUE` flag ensures `.rds` outputs are produced for programmatic comparison.

If you edit the C++ implementation or the R wrapper, rerun this check and inspect outputs before merging.

## Troubleshooting

- If `Rcpp::sourceCpp()` fails on macOS, ensure Xcode command-line tools are installed (`xcode-select --install`) and that your `R` can find clang/clang++.
- If you encounter unexpected numeric differences between run configurations, verify `eps_tiebreak` and floating-point tolerances.
- For reproducible CI, consider packaging the C++ code into an R package to avoid on-the-fly compilation variability.

## Implementation Details

The C++ implementation uses:

- Direct C++ loops for per-track and per-batch iteration
- Manual implementation of statistical functions (dnorm, dlnorm, log_sum_exp)
- Batch phase-feasibility checking (`derive_phase_prev_batch_rcpp`) for all (i, j) assignment pairs
- Optimized memory access patterns
- Reduced function call overhead relative to equivalent R code

All edge cases, measurement error models, biological constraints, and tie-breaking logic are implemented directly in C++.