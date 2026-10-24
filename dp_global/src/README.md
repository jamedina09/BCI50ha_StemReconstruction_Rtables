# dp_global/src

C++ functions of the `dp_global` engine: the transition costs and the
feasibility checks of the DP, and the exact assignment solver of the
probabilistic matcher. `dp_global/R/dp_global_main.R` compiles them when it
is sourced; the DP cannot run without them.

## Contents

| File | What it does |
| --- | --- |
| `transition_cost_rcpp.cpp` | The four C++ functions below |
| `transition_cost_rcpp.R` | R wrappers of the two cost functions, `transition_cost_tracks_bio_batch_rcpp()` and `transition_cost_paired_rcpp()` (default `hard_penalty = 1e6`) |

## How to run

From the project root:

```r
Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
source("dp_global/src/transition_cost_rcpp.R")

pars <- list(mu_const = 0, mu_gamma = 0, sigma0 = 1, sigma1 = 0,
             max_shrink = -Inf, k_shrink = 0, max_growth = Inf,
             max_growth_soft = Inf, k_growth = 0, use_measurement_error = FALSE,
             meas_sd1_a = 0.0062, meas_sd1_b = 0.0904, meas_sd2 = 4.64,
             meas_p_big = 0.05, h0 = 0, beta = 0,
             recruit_meanlog = 0, recruit_sdlog = 1, recruit_max_dbh = 200,
             recruit_lambda = 0, eps_tiebreak = 1e-6)
track_dbh_t <- c(10.0, NA, 20.0)                     # DBH per track at census t (cm)
mat_tp1 <- matrix(c(12.0, 15.0, 18.0, 8.0, NA, 22.0), nrow = 2, byrow = TRUE)  # two candidate states at t+1
do.call(transition_cost_tracks_bio_batch_rcpp,
        c(list(track_dbh_t = track_dbh_t, track_dbh_tp1 = mat_tp1, interval_years = 5), pars))
```

`transition_cost_rcpp.R` documents every argument.

## What it does

| Function | What it does |
| --- | --- |
| `transition_cost_paired_rcpp_cpp()` | Cost of every (state at t, state at t+1) pair of a census pair, in one call; used by the DP backward pass |
| `transition_cost_tracks_bio_batch_rcpp_cpp()` | Cost from one state at t to a batch of candidate states at t+1 |
| `derive_phase_prev_batch_rcpp()` | Feasibility of every state pair: phase rules, resprout rule and pruning bounds (growth, recruit size); returns the feasible pairs |
| `hungarian_min_rcpp()` | Exact minimum-cost assignment of a square matrix (Kuhn–Munkres, O(n³)); used by the birth-death mode of the matcher |

The cost is the negative log-likelihood of the transition summed over the
tracks: recruitment, mortality, the growth likelihood (Gaussian, or the
four-component measurement-error mixture), the hard bounds (`hard_penalty`),
the soft penalties and a rank tie-break (`dp_global/README.md`, *Transition
cost*). The statistical functions are written in C++, so no R function is
called per transition.

DBH recorded in classes: with `round_t` / `round_tp1` (default `FALSE`), a DBH
below `round_max_dbh` (5.5 cm) at census t / t+1 is read as a size in
[d, d + `round_width`) (0.5 cm): the growth likelihood uses the class
mid-points and adds `round_width`² / 12 per rounded value to the variance,
while the hard bounds stay on the measured values. The DP sets the flags from
`dbh_round_censuses`.

## Tests

`dp_global/tests/test_probabilistic_birth_death.R` tests `hungarian_min_rcpp()`.
After a change to the C++ code, rerun a driver on a fixed input and compare its
output with a saved reference:

```sh
Rscript dp_global/scripts/main_cpp_chunk.R --INPUT_FILE=data_simulation/data/simulated_data_1.csv --WRITE_DP_RDS=TRUE
```

## Requirements

R package `Rcpp` and a C++ compiler (macOS: `xcode-select --install`). Compile
on each machine; do not copy compiled files between machines.
