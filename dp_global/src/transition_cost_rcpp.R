#' @title Transition Cost Computation (Rcpp Accelerated)
#' @description Compute per-candidate negative-log-likelihood "costs" for possible
#' transitions between observed tracks at time t and candidate states at time t+1.
#' This function evaluates the following cases for each track-candidate pair:
#' (1) NA -> NA (no observation), (2) NA -> DBH (recruitment likelihood),
#' (3) DBH -> NA (mortality probability), and (4) DBH -> DBH (growth/mismatch likelihood).
#' It supports an optional measurement-error mixture model for growth likelihoods,
#' hard biological guardrails (hard_penalty applied when constraints violated),
#' and soft penalties for shrinkage and extreme growth. A small deterministic
#' rank-based tie-break cost may be added via `eps_tiebreak` to ensure reproducible
#' ordering of equal-cost candidates.
#'
#' All DBH units are in cm and growth rates are interpreted as cm/year (i.e.,
#' (d1 - d0) / interval_years). The returned vector has length equal to the number
#' of candidate rows in `mat_tp1` and contains the computed cost for each candidate.
#'
#' @param track_dbh_t Numeric vector of DBH values at time t (length K). Use NA for unobserved tracks.
#' @param mat_tp1 Numeric matrix or list-of-vectors of candidate DBH values at time t+1. Each row corresponds to one candidate and each column to a track (nrow = candidates, ncol = K).
#' @param interval_years Numeric; time interval in years between censuses (used to annualize growth).
#' @param mu_const Numeric; intercept for the size-dependent expected growth model E[g] = mu_const + mu_gamma * log(DBH). Default behavior when mu_gamma==0 is size-independent mean.
#' @param mu_gamma Numeric; slope for log(DBH) in the growth mean model (units: change in cm per unit log(DBH)).
#' @param sigma0 Numeric; baseline growth SD (cm/year).
#' @param sigma1 Numeric; additional growth SD proportional to DBH (cm/year per cm).
#' @param max_shrink Numeric; hard lower bound on annualized growth (cm/year, typically negative for shrink). Violations incur `hard_penalty`.
#' @param k_shrink Numeric; quadratic soft penalty coefficient applied when observed shrinkage occurs (units: 1/(cm^2)).
#' @param max_growth Numeric; hard upper bound on annualized growth (cm/year). Violations incur `hard_penalty`.
#' @param max_growth_soft Numeric; soft cap for positive growth (cm/year). Growth above this cap incurs a quadratic penalty controlled by `k_growth`.
#' @param k_growth Numeric; quadratic penalty coefficient for excess positive growth (units: 1/(cm^2)).
#' @param use_measurement_error Logical; if TRUE, evaluate the growth likelihood under a 4-component measurement-error mixture model (small/large errors) combined with process variance; otherwise use a Gaussian growth likelihood.
#' @param meas_sd1_a Numeric; slope for measurement error SD as function of DBH (small-error component), units: SD ≈ meas_sd1_a * DBH + meas_sd1_b.
#' @param meas_sd1_b Numeric; intercept for measurement error SD (small-error component).
#' @param meas_sd2 Numeric; SD for the large-error component (cm/year after conversion; used in mixture components).
#' @param meas_p_big Numeric; probability of a large measurement error (between 0 and 1) used to build the 4-component mixture weights.
#' @param h0 Numeric; mortality hazard intercept (hazard per year at DBH = 0). Used to compute death probability over `interval_years` via 1 - exp(-hazard * interval_years).
#' @param beta Numeric; mortality hazard coefficient multiplying DBH in the exponential hazard (hazard = h0 * exp(beta * DBH)).
#' @param recruit_meanlog Numeric; mean (on log scale) of recruit DBH distribution (log-normal) used for NA->DBH candidates.
#' @param recruit_sdlog Numeric; SD (on log scale) of recruit DBH distribution.
#' @param recruit_max_dbh Numeric; maximum DBH allowed for recruits (DBH > recruit_max_dbh is treated as impossible and penalized).
#' @param recruit_lambda Numeric; Poisson recruitment rate per year used to compute the probability of recruitment over `interval_years`.
#' @param eps_tiebreak Numeric; small coefficient multiplied by a deterministic rank-distance tie-break (default small positive value to break exact ties). Set to 0 to disable tie-breaking.
#' @param hard_penalty Numeric; large cost applied when hard biological constraints are violated (default 1e6).
#' @return Numeric vector of transition costs (one per candidate row of `mat_tp1`).
#' @export
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

#' @title Paired Transition Cost Computation (Rcpp Accelerated)
#' @description Compute transition costs for n_pairs paired source/destination
#' track-state vectors in a single C++ call. Both source and destination vary
#' per pair (unlike the batch version which fixes one source).
#' @param tdbh0_mat Numeric matrix [n_pairs x K] — source DBH per pair/track.
#' @param tdbh1_mat Numeric matrix [n_pairs x K] — destination DBH per pair/track.
#' @inheritParams transition_cost_tracks_bio_batch_rcpp
#' @return Numeric vector of length n_pairs — total NLL cost per pair.
#' @export
transition_cost_paired_rcpp <- function(
  tdbh0_mat,
  tdbh1_mat,
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
    transition_cost_paired_rcpp_cpp(
        tdbh0_mat = tdbh0_mat,
        tdbh1_mat = tdbh1_mat,
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
}