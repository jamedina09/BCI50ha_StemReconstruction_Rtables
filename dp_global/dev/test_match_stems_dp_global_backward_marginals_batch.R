############################################################
### test_match_stems_dp_global_backward_marginals_batch.R
############################################################

############################################################
# Run with:
# Rscript STEM_IDENTIFICATION_TEST/dp_global/dev/test_match_stems_dp_global_backward_marginals_batch.R
#
# Tests:
# - Batch DP runs without error (R version)
# - Batch DP runs without error (C++ version, if available)
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

# Minimal synthetic dataset that satisfies the DP's expectations:
# - multiple censuses
# - anchor census has TrueStemID present for all observations
# - Bio_* columns populated (constants are OK for a regression test)

tree_data <- data.table::data.table(
    Tag = rep(1L, 3L),
    species = rep("sp1", 3L),
    CensusID = 1:3,
    DBH = c(10, 11, 12),
    TrueStemID = c(
        NA_integer_,
        101L,
        NA_integer_
    )
)

# Populate Bio_* columns (constant across rows for the test)
tree_data[, `:=`(
    Bio_IntervalYears = 5,
    Bio_Mu_Growth = 0.0,
    Bio_Gamma_Growth = 0.0,
    Bio_Sigma0_Growth = 1.0,
    Bio_Sigma1_Growth = 0.00,
    Bio_Max_Shrink = -2.0,
    Bio_K_Shrink = 0.0,
    Bio_Max_Growth = 10.0,
    Bio_Max_Growth_Soft = 10.0,
    Bio_K_Growth = 0.0,
    Bio_H0_Mortality = 0.0,
    Bio_Beta_Mortality = 0.0,
    Bio_Recruit_Meanlog = log(2.0),
    Bio_Recruit_Sdlog = 0.3,
    Bio_Recruit_MaxDBH_unit = 6.0,
    Bio_Recruitment_lambda = 0.0
)]

args <- list(
    min_growth = -Inf,
    max_growth = Inf,
    anchor_start = 2L,
    max_tracks = 3L,
    slack_tracks = 1L,
    max_states = 50000L,
    temperature = 1.0,
    posterior_top_k = 2L,
    eps_tiebreak = 0,
    use_measurement_error = FALSE,
    verbose = FALSE
)

# Test with R implementation
out_batch <- do.call(
    match_stems_dp_global_backward_marginals_batch,
    c(list(tree_data = data.table::copy(tree_data)), args)
)

cat("OK: batch DP (R) runs without error on this test.\n")

# Test with C++ implementation if available
if (has_rcpp) {
    # Redefine transition_cost_tracks_bio_batch to use C++
    original_transition_cost_tracks_bio_batch <- transition_cost_tracks_bio_batch
    transition_cost_tracks_bio_batch <- function(
      track_dbh_t,
      track_dbh_tp1,
      interval_years,
      # -----------------------
      # GROWTH MODEL PARAMETERS
      # -----------------------
      mu_const = Bio_Mu_Growth_unit,
      mu_gamma = 0,
      sigma0 = Bio_Sigma0_unit,
      sigma1 = Bio_Sigma1_unit,
      max_shrink = Bio_max_shrink_unit,
      k_shrink = Bio_k_shrink_unit,
      max_growth = Inf,
      max_growth_soft = Inf,
      k_growth = 0,
      # -----------------
      # MEASUREMENT ERROR (optional)
      # -----------------
      use_measurement_error = FALSE,
      meas_sd1_a = 0.0062,
      meas_sd1_b = 0.0904,
      meas_sd2 = 4.64,
      meas_p_big = 0.05,
      # -------------------------
      # MORTALITY MODEL PARAMETERS
      # -------------------------
        h0 = Bio_H0_Mortality,
        beta = Bio_Beta_Mortality,
      # ----------------------------
      # RECRUITMENT MODEL PARAMETERS
      # ----------------------------
      recruit_meanlog = Bio_Recruit_Meanlog_unit,
      recruit_sdlog = Bio_Recruit_Sdlog_unit,
      recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
      recruit_lambda = Bio_Recruitment_lambda_unit,
      # -----------------
      # DETERMINISTIC TIE-BREAK
      # -----------------
      eps_tiebreak = 1e-6,
      hard_penalty = 1e6
    ) {
        if (is.list(track_dbh_tp1)) {
            mat_tp1 <- do.call(rbind, track_dbh_tp1)
        } else {
            mat_tp1 <- as.matrix(track_dbh_tp1)
        }
        result <- transition_cost_tracks_bio_batch_rcpp(
            track_dbh_t = track_dbh_t,
            mat_tp1 = mat_tp1,
            interval_years = interval_years,
            mu_const = mu_const,
            mu_gamma = mu_gamma,
            sigma0 = sigma0,
            sigma1 = sigma1,
            max_shrink = max_shrink,
            k_shrink = k_shrink,
            max_growth = max_growth,
            max_growth_soft = max_growth_soft,
            k_growth = k_growth,
            use_measurement_error = use_measurement_error,
            meas_sd1_a = meas_sd1_a,
            meas_sd1_b = meas_sd1_b,
            meas_sd2 = meas_sd2,
            meas_p_big = meas_p_big,
            h0 = h0,
            beta = beta,
            recruit_meanlog = recruit_meanlog,
            recruit_sdlog = recruit_sdlog,
            recruit_max_dbh = recruit_max_dbh,
            recruit_lambda = recruit_lambda,
            eps_tiebreak = eps_tiebreak,
            hard_penalty = hard_penalty
        )
        return(result)
    }

    out_batch_cpp <- do.call(
        match_stems_dp_global_backward_marginals_batch,
        c(list(tree_data = data.table::copy(tree_data)), args)
    )

    cat("OK: batch DP (C++) runs without error on this test.\n")

    # Restore original
    transition_cost_tracks_bio_batch <- original_transition_cost_tracks_bio_batch
}
