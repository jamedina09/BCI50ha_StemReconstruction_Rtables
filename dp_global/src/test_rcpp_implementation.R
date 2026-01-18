#' Standalone test script to validate C++ implementation of transition_cost_tracks_bio_batch
#'
#' This script compares the R and C++ implementations to ensure they produce
#' identical results.

# Load required packages
library(Rcpp)
library(microbenchmark)

if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this script.")
}
library(here)

# Source the original R function
source(here("dp_global","R","dp_global_biol.R"))

# Source the C++ wrapper
source(here("dp_global","src","transition_cost_rcpp.R"))

# Compile the C++ code
Rcpp::sourceCpp(here("dp_global","src","transition_cost_rcpp.cpp"))

# Compile the C++ code
Rcpp::sourceCpp("./dp_global/src/transition_cost_rcpp.cpp")

# Test parameters (using defaults from the original function)
test_params <- list(
    interval_years = 5,
    mu_const = 0.5,
    mu_gamma = 0.1,
    sigma0 = 0.5,
    sigma1 = 0.05,
    max_shrink = -2.0,
    k_shrink = 1.0,
    max_growth = 5.0,
    max_growth_soft = 3.0,
    k_growth = 0.5,
    use_measurement_error = TRUE,
    meas_sd1_a = 0.02,
    meas_sd1_b = 0.1,
    meas_sd2 = 1.0,
    meas_p_big = 0.1,
    h0 = 0.01,
    beta = 0.05,
    recruit_meanlog = 2.0,
    recruit_sdlog = 0.5,
    recruit_max_dbh = 200.0,
    recruit_lambda = 0.1,
    eps_tiebreak = 0.001,
    hard_penalty = 1000.0
)

# Test case 1: Simple case with 3 tracks, 2 candidates
cat("Test 1: Simple case\n")
track_dbh_t <- c(10.0, NA, 20.0)
mat_tp1 <- matrix(c(
    12.0, 15.0, 18.0,  # candidate 1
    8.0, NA, 22.0       # candidate 2
), nrow = 2, byrow = TRUE)

# Compute with R version
cost_r <- transition_cost_tracks_bio_batch(
    track_dbh_t = track_dbh_t,
    track_dbh_tp1 = mat_tp1,
    interval_years = test_params$interval_years,
    mu_const = test_params$mu_const,
    mu_gamma = test_params$mu_gamma,
    sigma0 = test_params$sigma0,
    sigma1 = test_params$sigma1,
    max_shrink = test_params$max_shrink,
    k_shrink = test_params$k_shrink,
    max_growth = test_params$max_growth,
    max_growth_soft = test_params$max_growth_soft,
    k_growth = test_params$k_growth,
    use_measurement_error = test_params$use_measurement_error,
    meas_sd1_a = test_params$meas_sd1_a,
    meas_sd1_b = test_params$meas_sd1_b,
    meas_sd2 = test_params$meas_sd2,
    meas_p_big = test_params$meas_p_big,
    h0 = test_params$h0,
    beta = test_params$beta,
    recruit_meanlog = test_params$recruit_meanlog,
    recruit_sdlog = test_params$recruit_sdlog,
    recruit_max_dbh = test_params$recruit_max_dbh,
    recruit_lambda = test_params$recruit_lambda,
    eps_tiebreak = test_params$eps_tiebreak,
    hard_penalty = test_params$hard_penalty
)

# Compute with C++ version
cost_cpp <- transition_cost_tracks_bio_batch_rcpp(
    track_dbh_t = track_dbh_t,
    mat_tp1 = mat_tp1,
    interval_years = test_params$interval_years,
    mu_const = test_params$mu_const,
    mu_gamma = test_params$mu_gamma,
    sigma0 = test_params$sigma0,
    sigma1 = test_params$sigma1,
    max_shrink = test_params$max_shrink,
    k_shrink = test_params$k_shrink,
    max_growth = test_params$max_growth,
    max_growth_soft = test_params$max_growth_soft,
    k_growth = test_params$k_growth,
    use_measurement_error = test_params$use_measurement_error,
    meas_sd1_a = test_params$meas_sd1_a,
    meas_sd1_b = test_params$meas_sd1_b,
    meas_sd2 = test_params$meas_sd2,
    meas_p_big = test_params$meas_p_big,
    h0 = test_params$h0,
    beta = test_params$beta,
    recruit_meanlog = test_params$recruit_meanlog,
    recruit_sdlog = test_params$recruit_sdlog,
    recruit_max_dbh = test_params$recruit_max_dbh,
    recruit_lambda = test_params$recruit_lambda,
    eps_tiebreak = test_params$eps_tiebreak,
    hard_penalty = test_params$hard_penalty
)

cat("R version results:", cost_r, "\n")
cat("C++ version results:", cost_cpp, "\n")
max_diff <- max(abs(cost_r - cost_cpp))
cat("Max difference:", max_diff, "\n")

if (max_diff < 1e-10) {
    cat("✓ Test 1 PASSED\n\n")
} else {
    cat("✗ Test 1 FAILED\n\n")
    stop("Test 1 failed - results don't match")
}

# Test case 2: Larger batch with more complex data
cat("Test 2: Larger batch\n")
set.seed(42)  # For reproducible results
K <- 5  # 5 tracks
n_batch <- 20  # 20 candidates

track_dbh_t <- c(15.0, NA, 25.0, 8.0, 30.0)

# Generate random candidate matrices
mat_tp1 <- matrix(NA, nrow = n_batch, ncol = K)
for (i in 1:n_batch) {
    for (j in 1:K) {
        if (runif(1) < 0.8) {  # 80% chance of having a DBH
            mat_tp1[i, j] <- rnorm(1, mean = 20, sd = 10)
            mat_tp1[i, j] <- max(1, mat_tp1[i, j])  # Ensure positive
        }
    }
}

# Compute with both versions
cost_r2 <- transition_cost_tracks_bio_batch(
    track_dbh_t = track_dbh_t,
    track_dbh_tp1 = mat_tp1,
    interval_years = test_params$interval_years,
    mu_const = test_params$mu_const,
    mu_gamma = test_params$mu_gamma,
    sigma0 = test_params$sigma0,
    sigma1 = test_params$sigma1,
    max_shrink = test_params$max_shrink,
    k_shrink = test_params$k_shrink,
    max_growth = test_params$max_growth,
    max_growth_soft = test_params$max_growth_soft,
    k_growth = test_params$k_growth,
    use_measurement_error = test_params$use_measurement_error,
    meas_sd1_a = test_params$meas_sd1_a,
    meas_sd1_b = test_params$meas_sd1_b,
    meas_sd2 = test_params$meas_sd2,
    meas_p_big = test_params$meas_p_big,
    h0 = test_params$h0,
    beta = test_params$beta,
    recruit_meanlog = test_params$recruit_meanlog,
    recruit_sdlog = test_params$recruit_sdlog,
    recruit_max_dbh = test_params$recruit_max_dbh,
    recruit_lambda = test_params$recruit_lambda,
    eps_tiebreak = test_params$eps_tiebreak,
    hard_penalty = test_params$hard_penalty
)

cost_cpp2 <- transition_cost_tracks_bio_batch_rcpp(
    track_dbh_t = track_dbh_t,
    mat_tp1 = mat_tp1,
    interval_years = test_params$interval_years,
    mu_const = test_params$mu_const,
    mu_gamma = test_params$mu_gamma,
    sigma0 = test_params$sigma0,
    sigma1 = test_params$sigma1,
    max_shrink = test_params$max_shrink,
    k_shrink = test_params$k_shrink,
    max_growth = test_params$max_growth,
    max_growth_soft = test_params$max_growth_soft,
    k_growth = test_params$k_growth,
    use_measurement_error = test_params$use_measurement_error,
    meas_sd1_a = test_params$meas_sd1_a,
    meas_sd1_b = test_params$meas_sd1_b,
    meas_sd2 = test_params$meas_sd2,
    meas_p_big = test_params$meas_p_big,
    h0 = test_params$h0,
    beta = test_params$beta,
    recruit_meanlog = test_params$recruit_meanlog,
    recruit_sdlog = test_params$recruit_sdlog,
    recruit_max_dbh = test_params$recruit_max_dbh,
    recruit_lambda = test_params$recruit_lambda,
    eps_tiebreak = test_params$eps_tiebreak,
    hard_penalty = test_params$hard_penalty
)

max_diff2 <- max(abs(cost_r2 - cost_cpp2))
cat("Max difference:", max_diff2, "\n")

if (max_diff2 < 1e-10) {
    cat("✓ Test 2 PASSED\n\n")
} else {
    cat("✗ Test 2 FAILED\n\n")
    stop("Test 2 failed - results don't match")
}

# Performance comparison
cat("Performance comparison:\n")
bench_result <- microbenchmark(
    R_version = transition_cost_tracks_bio_batch(
        track_dbh_t = track_dbh_t,
        track_dbh_tp1 = mat_tp1,
        interval_years = test_params$interval_years,
        mu_const = test_params$mu_const,
        mu_gamma = test_params$mu_gamma,
        sigma0 = test_params$sigma0,
        sigma1 = test_params$sigma1,
        max_shrink = test_params$max_shrink,
        k_shrink = test_params$k_shrink,
        max_growth = test_params$max_growth,
        max_growth_soft = test_params$max_growth_soft,
        k_growth = test_params$k_growth,
        use_measurement_error = test_params$use_measurement_error,
        meas_sd1_a = test_params$meas_sd1_a,
        meas_sd1_b = test_params$meas_sd1_b,
        meas_sd2 = test_params$meas_sd2,
        meas_p_big = test_params$meas_p_big,
        h0 = test_params$h0,
        beta = test_params$beta,
        recruit_meanlog = test_params$recruit_meanlog,
        recruit_sdlog = test_params$recruit_sdlog,
        recruit_max_dbh = test_params$recruit_max_dbh,
        recruit_lambda = test_params$recruit_lambda,
        eps_tiebreak = test_params$eps_tiebreak,
        hard_penalty = test_params$hard_penalty
    ),
    CPP_version = transition_cost_tracks_bio_batch_rcpp(
        track_dbh_t = track_dbh_t,
        mat_tp1 = mat_tp1,
        interval_years = test_params$interval_years,
        mu_const = test_params$mu_const,
        mu_gamma = test_params$mu_gamma,
        sigma0 = test_params$sigma0,
        sigma1 = test_params$sigma1,
        max_shrink = test_params$max_shrink,
        k_shrink = test_params$k_shrink,
        max_growth = test_params$max_growth,
        max_growth_soft = test_params$max_growth_soft,
        k_growth = test_params$k_growth,
        use_measurement_error = test_params$use_measurement_error,
        meas_sd1_a = test_params$meas_sd1_a,
        meas_sd1_b = test_params$meas_sd1_b,
        meas_sd2 = test_params$meas_sd2,
        meas_p_big = test_params$meas_p_big,
        h0 = test_params$h0,
        beta = test_params$beta,
        recruit_meanlog = test_params$recruit_meanlog,
        recruit_sdlog = test_params$recruit_sdlog,
        recruit_max_dbh = test_params$recruit_max_dbh,
        recruit_lambda = test_params$recruit_lambda,
        eps_tiebreak = test_params$eps_tiebreak,
        hard_penalty = test_params$hard_penalty
    ),
    times = 10
)

print(bench_result)

# Calculate speedup
r_times <- summary(bench_result)$mean[summary(bench_result)$expr == "R_version"]
cpp_times <- summary(bench_result)$mean[summary(bench_result)$expr == "CPP_version"]
speedup <- r_times / cpp_times
cat("\nSpeedup factor:", speedup, "\n")

cat("\n✓ All tests passed! C++ implementation is working correctly.\n")
