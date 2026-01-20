#!/usr/bin/env Rscript
library(data.table)
# Source the DP code (only functional defs)
source(here::here("dp_global","R","dp_global_biol.R"))
# Minimal synthetic dataset
tree_data <- data.table::data.table(
    Tag = rep(1L, 3L),
    species = rep("sp1", 3L),
    CensusID = 1:3,
    DBH = c(10, 11, 12),
    TrueStemID = c(NA_integer_, 101L, NA_integer_)
)
# Populate Bio_* constants (we'll jitter per-row below)
tree_data[, `:=`(
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
# Add jitter per-census
set.seed(123)
means_by_c <- c(5, 3, 7)
for (ci in unique(tree_data$CensusID)) {
    idx <- which(tree_data$CensusID == ci)
    tree_data[idx, Bio_IntervalYears := means_by_c[as.integer(ci)] + rnorm(length(idx), 0, 1e-3)]
}
cat("Bio_IntervalYears per row:\n")
print(tree_data$Bio_IntervalYears)
# Run DP (R version) - use marginal batch entrypoint for coverage
res <- match_stems_dp_global_backward_marginals_batch(
    tree_data = data.table::copy(tree_data),
    min_growth = -Inf,
    max_growth = Inf,
    anchor_start = 2L,
    max_tracks = 3L,
    max_states = 50000L,
    slack_tracks = 1L,
    temperature = 1.0,
    posterior_top_k = 2L,
    use_measurement_error = FALSE,
    verbose = TRUE
)
cat("DP result (R):\n")
print(res)
# If Rcpp available, test C++ path by redefining the batch function
if (requireNamespace("Rcpp", quietly = TRUE)) {
    source(here::here("dp_global","src","transition_cost_rcpp.R"))
    Rcpp::sourceCpp(here::here("dp_global","src","transition_cost_rcpp.cpp"))
    original_transition_cost_tracks_bio_batch <- transition_cost_tracks_bio_batch
    transition_cost_tracks_bio_batch <- function(
        track_dbh_t,
        track_dbh_tp1,
        interval_years,
        mu_const = Bio_Mu_Growth_unit,
        mu_gamma = 0,
        sigma0 = Bio_Sigma0_unit,
        sigma1 = Bio_Sigma1_unit,
        max_shrink = Bio_max_shrink_unit,
        k_shrink = Bio_k_shrink_unit,
        max_growth = Inf,
        max_growth_soft = Inf,
        k_growth = 0,
        use_measurement_error = FALSE,
        meas_sd1_a = 0.0062,
        meas_sd1_b = 0.0904,
        meas_sd2 = 4.64,
        meas_p_big = 0.05,
        h0 = Bio_H0_Mortality,
        beta = Bio_Beta_Mortality,
        recruit_meanlog = Bio_Recruit_Meanlog_unit,
        recruit_sdlog = Bio_Recruit_Sdlog_unit,
        recruit_max_dbh = Bio_Recruit_MaxDBH_unit,
        recruit_lambda = Bio_Recruitment_lambda_unit,
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
    res_cpp <- match_stems_dp_global_backward_marginals_batch(
        tree_data = data.table::copy(tree_data),
        min_growth = -Inf,
        max_growth = Inf,
        anchor_start = 2L,
        max_tracks = 3L,
        max_states = 50000L,
        slack_tracks = 1L,
        temperature = 1.0,
        posterior_top_k = 2L,
        use_measurement_error = FALSE,
        verbose = TRUE
    )
    cat("DP result (C++):\n")
    print(res_cpp)
}
cat("Done test\n")
