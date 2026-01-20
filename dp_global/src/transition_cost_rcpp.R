#' @title Transition Cost Computation (Rcpp Accelerated)
#' @description Computes transition costs for biological stem tracking using batched processing
#' @param track_dbh_t Numeric vector of DBH values at time t (NA for unobserved)
#' @param mat_tp1 Matrix of candidate DBH values at time t+1 (rows = candidates, cols = tracks)
#' @param interval_years Time interval in years
#' @param mu_const Growth model intercept
#' @param mu_gamma Growth model slope for log(DBH)
#' @param sigma0 Growth model baseline SD
#' @param sigma1 Growth model SD coefficient for DBH
#' @param max_shrink Maximum allowed shrinkage (cm/year, negative)
#' @param k_shrink Soft penalty coefficient for shrinkage
#' @param max_growth Maximum allowed growth (cm/year)
#' @param max_growth_soft Soft cap for growth (cm/year)
#' @param k_growth Soft penalty coefficient for excessive growth
#' @param use_measurement_error Whether to use measurement error model
#' @param meas_sd1_a Measurement error SD coefficient for DBH
#' @param meas_sd1_b Measurement error SD intercept
#' @param meas_sd2 Measurement error SD for large errors
#' @param meas_p_big Probability of large measurement error
#' @param h0 Mortality hazard intercept
#' @param beta Mortality hazard coefficient for DBH
#' @param recruit_meanlog Recruitment log-mean
#' @param recruit_sdlog Recruitment log-SD
#' @param recruit_max_dbh Maximum DBH for recruitment
#' @param recruit_lambda Recruitment rate parameter
#' @param eps_tiebreak Tie-breaking coefficient
#' @param hard_penalty Hard constraint penalty
#' @return Numeric vector of transition costs (one per candidate)
transition_cost_tracks_bio_batch_rcpp <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  mu_const,
  mu_gamma,
  sigma0,
  sigma1,
  max_shrink,
  k_shrink,
  max_growth,
  max_growth_soft,
  k_growth,
  use_measurement_error,
  meas_sd1_a,
  meas_sd1_b,
  meas_sd2,
  meas_p_big,
  h0,
  beta,
  recruit_meanlog,
  recruit_sdlog,
  recruit_max_dbh,
  recruit_lambda,
  eps_tiebreak,
  hard_penalty = 1e6
) {
    # cat("C++ version called\n")
    if (is.list(track_dbh_tp1)) {
        mat_tp1 <- do.call(rbind, track_dbh_tp1)
    } else {
        mat_tp1 <- as.matrix(track_dbh_tp1)
    }
    result <- transition_cost_tracks_bio_batch_rcpp_cpp(
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