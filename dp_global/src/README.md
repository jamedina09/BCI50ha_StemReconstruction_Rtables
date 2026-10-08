# Rcpp Acceleration for DP Global Stem Tracking

C++ implementation and R wrappers for Rcpp-accelerated transition-cost and phase-feasibility computations used by the DP workflow, and the exact assignment solver used by the probabilistic matcher.

## Files

- `transition_cost_rcpp.cpp`: C++ implementation of transition cost computation, phase-feasibility batch checks and minimum-cost assignment
- `transition_cost_rcpp.R`: R wrapper functions for the two transition-cost functions (`transition_cost_tracks_bio_batch_rcpp()`, `transition_cost_paired_rcpp()`)

## Usage

### Prerequisites

R (package: `Rcpp`) and a C++ compiler.

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

See `transition_cost_rcpp.R` (this folder) for detailed parameter descriptions and expected types. `hard_penalty` defaults to `1e6` in the R wrappers.

Optional arguments `round_t`, `round_tp1` (default `FALSE`), `round_max_dbh` (5.5 cm) and `round_width` (0.5 cm): when census t (or t+1) recorded DBH in classes rounded down, a DBH below `round_max_dbh` is read as a true size in `[d, d + round_width)` — the growth likelihood uses the class mid-points and adds `round_width^2 / 12` per rounded value to the variance of the DBH difference (divided by `interval_years^2` for the annual growth), while the hard growth limits stay on the measured values. With both flags `FALSE` rounding has no effect on the costs. The DP sets the flags from `dbh_round_censuses` (see `dp_global/README.md`, *DBH recorded in classes*).

## Stem Identity Renumbering

All drivers in the dp_global workflow use a universal post-engine helper chain: `maybe_add_posterior_bins()`, `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, and finally `finalize_posterior_paths()`. After these steps, all `ReconstructedStemID` values are renumbered **sequentially from 1 to N within each tag**, ordered by the earliest census in which each stem appears. If multiple stems first appear in the same census, the largest DBH at that census gets the lower ID, with ties broken by original ID. **Negative or zero IDs are never produced.**

## Performance

The C++ implementation is the sole backend for transition-cost computation in the DP workflow: the backward pass calls `derive_phase_prev_batch_rcpp()` and then `transition_cost_paired_rcpp()` once per census pair (`dp_global/R/dp_global_dp.R`). Per-track and per-pair iteration runs in C++ loops, with the statistical functions written out in C++, so no R function is called per transition.

## Validation

To validate the C++ implementation after changes, run an end-to-end check through `main_cpp_chunk.R`:

```bash
Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WRITE_DP_RDS=TRUE
```

Compare chunk RDS output against a saved reference to detect regressions. `--WRITE_DP_RDS=TRUE` (the default of `main_cpp_chunk.R`) writes the `.rds` outputs used for programmatic comparison.

If you edit the C++ implementation or the R wrapper, rerun this check and inspect outputs before merging.

## Troubleshooting

- If `Rcpp::sourceCpp()` fails on macOS, ensure Xcode command-line tools are installed (`xcode-select --install`) and that your `R` can find clang/clang++.
- If you encounter unexpected numeric differences between run configurations, verify `eps_tiebreak` and floating-point tolerances.
- For reproducible CI, consider packaging the C++ code into an R package to avoid on-the-fly compilation variability.

## Implementation Details

The C++ implementation provides four exported functions:

- `transition_cost_tracks_bio_batch_rcpp_cpp()`: Computes transition costs for a batch of candidate next-states given a single current state (R wrapper `transition_cost_tracks_bio_batch_rcpp()`; the DP itself calls the paired function)
- `transition_cost_paired_rcpp_cpp()`: Computes transition costs for pre-paired (current, next) state matrices — used for the batched backward pass, one call per census pair (R wrapper `transition_cost_paired_rcpp()`)
- `derive_phase_prev_batch_rcpp()`: Batch phase-feasibility checking with integrated hard pruning (growth bounds and recruit size) — determines valid phase transitions for all (i, j) assignment pairs and returns only feasible pairs
- `hungarian_min_rcpp()`: Exact minimum-cost assignment of a square cost matrix (Kuhn-Munkres, O(n³)) — used by the birth-death mode of the probabilistic matcher (`greedy_assignment_gumbel()` in `dp_global/R/dp_probabilistic_matching.R`)

The implementation uses:

- Direct C++ loops for per-track and per-batch iteration
- Manual implementation of statistical functions (dnorm, dlnorm, log_sum_exp)
- Batch phase-feasibility checking (`derive_phase_prev_batch_rcpp`) for all (i, j) assignment pairs
- Reduced function call overhead relative to equivalent R code

The measurement error model, DBH rounding, biological constraints, and tie-breaking logic are implemented directly in C++.
