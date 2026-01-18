############################################################
### test_transition_cost_tracks_bio_batch.R
############################################################

############################################################
# Run with:
# Rscript STEM_IDENTIFICATION_TEST/dp_global/dev/test_transition_cost_tracks_bio_batch.R
#
# Tests:
# - Batch transition cost matches scalar transition cost
# - C++ batch transition cost matches R batch (if Rcpp available)
############################################################

find_project_root <- function(start_dir) {
    d <- normalizePath(start_dir)
    for (i in 0:6) {
        cand <- d
        if (dir.exists(file.path(cand, "STEM_IDENTIFICATION_TEST"))) {
            return(cand)
        }
        d2 <- dirname(d)
        if (identical(d2, d)) break
        d <- d2
    }
    stop("Could not find project root containing STEM_IDENTIFICATION_TEST/ starting from: ", start_dir)
}

# Make the script runnable from anywhere
root_dir <- find_project_root(getwd())
setwd(file.path(root_dir, "STEM_IDENTIFICATION_TEST"))

if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this script.")
}
library(here)

source(here("dp_global","R","dp_global_biol.R"))

# Load C++ implementation for additional testing
if (requireNamespace("Rcpp", quietly = TRUE)) {
    library(Rcpp)
    source(here("dp_global","src","transition_cost_rcpp.R"))
    Rcpp::sourceCpp(here("dp_global","src","transition_cost_rcpp.cpp"))
    has_rcpp <- TRUE
} else {
    has_rcpp <- FALSE
    cat("Rcpp not available; skipping C++ tests.\n")
}

set.seed(1)

# Small synthetic example with NAs and mixed cases
K <- 6
n_batch <- 25

track_dbh_t <- c(10, 20, NA, 5, 30, NA)

# Candidate next states (rows)
mat_tp1 <- matrix(NA_real_, nrow = n_batch, ncol = K)
for (i in seq_len(n_batch)) {
    row <- track_dbh_t

    # Randomly flip some tracks to NA / DBH
    for (k in seq_len(K)) {
        if (is.na(row[k])) {
            # NA->NA or NA->DBH
            if (runif(1) < 0.4) {
                row[k] <- NA_real_
            } else {
                row[k] <- exp(rnorm(1, log(2), 0.4))
            }
        } else {
            # DBH->DBH or DBH->NA
            if (runif(1) < 0.2) {
                row[k] <- NA_real_
            } else {
                row[k] <- row[k] + rnorm(1, mean = 0.5, sd = 2.0)
            }
        }
    }

    mat_tp1[i, ] <- row
}

pars <- list(
    interval_years = 5,
    mu_const = 0.8,
    mu_gamma = 0.1,
    sigma0 = 0.2,
    sigma1 = 0.01,
    max_shrink = -1,
    k_shrink = 25,
    max_growth = 5,
    max_growth_soft = 1.0,
    k_growth = 2,
    use_measurement_error = TRUE,
    meas_sd1_a = 0.0062,
    meas_sd1_b = 0.0904,
    meas_sd2 = 4.64,
    meas_p_big = 0.05,
    h0 = 0.01,
    beta = 0.0,
    recruit_meanlog = log(2),
    recruit_sdlog = 0.5,
    recruit_max_dbh = 10,
    recruit_lambda = 0.1,
    eps_tiebreak = 1e-6,
    hard_penalty = 1e6
)

# Scalar baseline (transition_cost_tracks_bio does not have hard_penalty)
pars_scalar <- pars
pars_scalar$hard_penalty <- NULL

scalar <- vapply(
    seq_len(n_batch),
    function(i) {
        do.call(
            transition_cost_tracks_bio,
            c(list(track_dbh_t = track_dbh_t, track_dbh_tp1 = mat_tp1[i, ]), pars_scalar)
        )
    },
    numeric(1L)
)

# Batched
batched <- do.call(
    transition_cost_tracks_bio_batch,
    c(list(track_dbh_t = track_dbh_t, track_dbh_tp1 = mat_tp1), pars)
)

max_abs_diff <- max(abs(scalar - batched))
cat("max_abs_diff=", format(max_abs_diff, digits = 10), "\n", sep = "")

stopifnot(is.finite(max_abs_diff), max_abs_diff < 1e-8)
cat("OK: transition_cost_tracks_bio_batch matches transition_cost_tracks_bio on this test.\n")

# Test C++ implementation if available
if (has_rcpp) {
    cpp <- do.call(
        transition_cost_tracks_bio_batch_rcpp,
        c(list(track_dbh_t = track_dbh_t, mat_tp1 = mat_tp1), pars)
    )
    
    max_abs_diff_cpp <- max(abs(batched - cpp))
    cat("max_abs_diff_cpp=", format(max_abs_diff_cpp, digits = 10), "\n", sep = "")
    
    stopifnot(is.finite(max_abs_diff_cpp), max_abs_diff_cpp < 1e-8)
    cat("OK: transition_cost_tracks_bio_batch_rcpp matches transition_cost_tracks_bio_batch on this test.\n")
} else {
    cat("Skipping C++ test (Rcpp not available).\n")
}
