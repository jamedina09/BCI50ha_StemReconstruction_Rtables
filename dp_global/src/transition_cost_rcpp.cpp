#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

// ---------------------------------------------------------------------------
// transition_cost_tracks_bio_batch_rcpp_cpp
//
// Computes the total negative-log-likelihood transition cost for a batch of
// candidate track-state vectors relative to a current track-state vector.
//
// Arguments:
//   track_dbh_t    : NumericVector of length K — current DBH values per track
//                    (NA = track not currently observed).
//   mat_tp1        : NumericMatrix of dimensions n_batch × K — candidate next-
//                    state DBH values; each row is one candidate assignment.
//   interval_years : Length of the census interval in years.
//   mu_const       : Intercept of the linear growth mean (cm/year).
//   mu_gamma       : Slope on log(DBH) for growth mean (0 = constant mean).
//   sigma0, sigma1 : Intercept and slope of the process SD model:
//                      SD(D) = sigma0 + sigma1 * D.
//   max_shrink     : Hard lower bound on annual growth (cm/year); transitions
//                    below this rate receive hard_penalty.
//   k_shrink       : Soft quadratic penalty weight for shrinkage below 0.
//   max_growth     : Hard upper bound on annual growth (cm/year); transitions
//                    above this rate receive hard_penalty.
//   max_growth_soft: Soft threshold for excess growth penalty.
//   k_growth       : Soft quadratic penalty weight for growth exceeding
//                    max_growth_soft.
//   use_measurement_error : When TRUE, growth likelihood is computed as a
//                    4-component mixture (process + 3 measurement-error
//                    components) rather than a single Gaussian.
//   meas_sd1_a, meas_sd1_b : Parameters of the small-error SD model:
//                              SD1(D) = meas_sd1_a * D + meas_sd1_b.
//   meas_sd2       : SD of the large measurement-error component.
//   meas_p_big     : Mixing weight for the large measurement-error component.
//   h0, beta       : Hazard model parameters for mortality probability:
//                      p_die = 1 − exp(−h0 * exp(beta * D) * interval_years).
//   recruit_meanlog, recruit_sdlog : Log-normal parameters for recruit DBH.
//   recruit_max_dbh : Hard cap on recruit DBH; exceeding this gives hard_penalty.
//   recruit_lambda : Poisson recruitment rate (per year) for p_recruit.
//   eps_tiebreak   : Small deterministic rank-sum term added to each row cost
//                    to break ties in a reproducible way.
//   hard_penalty   : Cost assigned to biologically impossible transitions
//                    (default 1e6).
//
// Per-track cases evaluated for each candidate row:
//   NA → NA   : log(1 − p_recruit) — track absent in both censuses.
//   NA → DBH  : log(p_recruit) + dlnorm(DBH | recruit_meanlog, recruit_sdlog)
//               − cap at recruit_max_dbh (hard_penalty if exceeded).
//   DBH → NA  : log(p_die) — observed tree goes absent (mortality).
//   DBH → DBH : Growth likelihood (Gaussian or 4-component mixture) plus soft
//               and hard growth-rate penalties.
//
// Returns: NumericVector of length n_batch — total (summed across tracks)
//          negative-log-likelihood cost for each candidate, plus eps_tiebreak
//          scaled by the rank sum of current DBH values for tie-breaking.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector transition_cost_tracks_bio_batch_rcpp_cpp(
    const Rcpp::NumericVector& track_dbh_t,
    const Rcpp::NumericMatrix& mat_tp1,
    double interval_years,
    double mu_const,
    double mu_gamma,
    double sigma0,
    double sigma1,
    double max_shrink,
    double k_shrink,
    double max_growth,
    double max_growth_soft,
    double k_growth,
    bool use_measurement_error,
    double meas_sd1_a,
    double meas_sd1_b,
    double meas_sd2,
    double meas_p_big,
    double h0,
    double beta,
    double recruit_meanlog,
    double recruit_sdlog,
    double recruit_max_dbh,
    double recruit_lambda,
    double eps_tiebreak,
    double hard_penalty
) {
    int K = track_dbh_t.size();
    int n_batch = mat_tp1.nrow();

    if (mat_tp1.ncol() != K) {
        Rcpp::stop("track_dbh_tp1 must have K columns");
    }

    Rcpp::NumericVector cost(n_batch, 0.0);

    // Recruitment probability over interval
    double p_recruit = 1.0 - std::exp(-recruit_lambda * interval_years);
    p_recruit = std::max(1e-12, std::min(1.0 - 1e-12, p_recruit));

    // Lambda function for growth mean
    auto mu_growth = [&](double d) -> double {
        if (!std::isfinite(mu_gamma) || mu_gamma == 0.0 || !std::isfinite(d) || d <= 0.0) {
            return mu_const;
        }
        return mu_const + mu_gamma * std::log(d);
    };

    // Lambda function for measurement error SD
    auto meas_sd1 = [&](double d) -> double {
        return std::max(1e-6, meas_sd1_a * d + meas_sd1_b);
    };

    // Lambda function for log-sum-exp
    auto log_sum_exp = [](const std::vector<double>& x) -> double {
        if (x.empty()) return -INFINITY;
        double max_val = *std::max_element(x.begin(), x.end());
        double sum = 0.0;
        for (double val : x) {
            sum += std::exp(val - max_val);
        }
        return max_val + std::log(sum);
    };

    // Lambda function for normal log density
    auto dnorm_log = [](double x, double mean, double sd) -> double {
        double diff = x - mean;
        return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * diff * diff / (sd * sd);
    };

    // Lambda function for log-normal log density
    auto dlnorm_log = [](double x, double meanlog, double sdlog) -> double {
        if (x <= 0.0) return -INFINITY;
        double log_x = std::log(x);
        double diff = log_x - meanlog;
        return -std::log(x * sdlog * std::sqrt(2.0 * M_PI)) - 0.5 * diff * diff / (sdlog * sdlog);
    };

    // Loop over tracks (vectorized across candidates)
    for (int k = 0; k < K; k++) {
        double d0 = track_dbh_t[k];
        Rcpp::NumericVector d1_vec = mat_tp1(Rcpp::_, k);

        // CASE 1 + 2: NA -> *
        if (Rcpp::NumericVector::is_na(d0)) {
            for (int i = 0; i < n_batch; i++) {
                double d1 = d1_vec[i];

                if (Rcpp::NumericVector::is_na(d1)) {
                    // NA -> NA
                    cost[i] -= std::log(1.0 - p_recruit);
                } else {
                    // NA -> DBH
                    bool hard = (!std::isfinite(d1)) || (d1 <= 0.0) || (d1 > recruit_max_dbh);
                    if (hard) {
                        cost[i] += hard_penalty;
                    } else {
                        cost[i] -= std::log(p_recruit) + dlnorm_log(d1, recruit_meanlog, recruit_sdlog);
                    }
                }
            }
            continue;
        }

        // CASE 3: DBH -> NA
        for (int i = 0; i < n_batch; i++) {
            double d1 = d1_vec[i];
            if (!Rcpp::NumericVector::is_na(d1)) continue;

            double hazard = h0 * std::exp(beta * d0);
            double p_death = 1.0 - std::exp(-hazard * interval_years);
            p_death = std::max(1e-12, std::min(1.0 - 1e-12, p_death));
            cost[i] -= std::log(p_death);
        }

        // CASE 4: DBH -> DBH
        for (int i = 0; i < n_batch; i++) {
            double d1 = d1_vec[i];
            if (Rcpp::NumericVector::is_na(d1)) continue;

            double g = (d1 - d0) / interval_years;

            // Hard biological constraints
            bool hard = false;
            if (std::isfinite(max_shrink) && (g < max_shrink)) hard = true;
            if (std::isfinite(max_growth) && (g > max_growth)) hard = true;

            if (hard) {
                cost[i] += hard_penalty;
                continue;
            }

            // Size-dependent growth variance
            double sigma_d = sigma0 + sigma1 * d0;
            sigma_d = std::max(sigma_d, 1e-6);
            double mu = mu_growth(d0);

            if (use_measurement_error) {
                // Measurement-error-aware likelihood
                double s_small0 = meas_sd1(d0);
                double s_small1 = meas_sd1(d1);
                double s_big = meas_sd2;
                double w_small = 1.0 - meas_p_big;
                double w_big = meas_p_big;

                // Build 4-component mixture
                std::vector<double> sd_meas_mix = {
                    std::sqrt(s_small0*s_small0 + s_small1*s_small1) / interval_years,
                    std::sqrt(s_small0*s_small0 + s_big*s_big) / interval_years,
                    std::sqrt(s_big*s_big + s_small1*s_small1) / interval_years,
                    std::sqrt(s_big*s_big + s_big*s_big) / interval_years
                };

                std::vector<double> wt_meas_mix = {
                    w_small * w_small,
                    w_small * w_big,
                    w_big * w_small,
                    w_big * w_big
                };

                std::vector<double> sd_tot(4);
                for (int j = 0; j < 4; j++) {
                    sd_tot[j] = std::sqrt(sigma_d*sigma_d + sd_meas_mix[j]*sd_meas_mix[j]);
                }

                std::vector<double> ll(4);
                for (int j = 0; j < 4; j++) {
                    ll[j] = std::log(wt_meas_mix[j]) + dnorm_log(g, mu, sd_tot[j]);
                }

                cost[i] -= log_sum_exp(ll);
            } else {
                // Gaussian growth likelihood
                double diff = g - mu;
                cost[i] += diff*diff / (2.0 * sigma_d*sigma_d) + std::log(sigma_d) + 0.5 * std::log(2.0 * M_PI);
            }

            // Soft penalty for shrinkage
            if (std::isfinite(k_shrink) && k_shrink > 0.0 && d1 < d0) {
                double dd = d0 - d1;
                cost[i] += k_shrink * dd * dd;
            }

            // Soft penalty for extreme positive growth
            if (std::isfinite(max_growth_soft) && std::isfinite(k_growth) && k_growth > 0.0) {
                double d1_soft_cap = d0 + max_growth_soft * interval_years;
                if (std::isfinite(d1_soft_cap) && d1 > d1_soft_cap) {
                    double dd = d1 - d1_soft_cap;
                    cost[i] += k_growth * dd * dd;
                }
            }
        }
    }

    // Deterministic tie-break
    if (eps_tiebreak > 0.0) {
        // Compute ranks for track_dbh_t
        std::vector<std::pair<double, int>> ranked_t;
        for (int k = 0; k < K; k++) {
            if (!Rcpp::NumericVector::is_na(track_dbh_t[k])) {
                ranked_t.emplace_back(track_dbh_t[k], k);
            }
        }
        std::sort(ranked_t.begin(), ranked_t.end());
        std::vector<int> r0(K, 0);
        for (size_t i = 0; i < ranked_t.size(); i++) {
            r0[ranked_t[i].second] = i + 1;
        }

        for (int i = 0; i < n_batch; i++) {
            // Compute ranks for mat_tp1[i, ]
            std::vector<std::pair<double, int>> ranked_tp1;
            for (int k = 0; k < K; k++) {
                if (!Rcpp::NumericVector::is_na(mat_tp1(i, k))) {
                    ranked_tp1.emplace_back(mat_tp1(i, k), k);
                }
            }
            std::sort(ranked_tp1.begin(), ranked_tp1.end());
            std::vector<int> r1(K, 0);
            for (size_t j = 0; j < ranked_tp1.size(); j++) {
                r1[ranked_tp1[j].second] = j + 1;
            }

            // Sum absolute differences for tracks that are observed in both
            double tie_break = 0.0;
            for (int k = 0; k < K; k++) {
                if (!Rcpp::NumericVector::is_na(track_dbh_t[k]) && !Rcpp::NumericVector::is_na(mat_tp1(i, k))) {
                    tie_break += std::abs(r0[k] - r1[k]);
                }
            }
            cost[i] += eps_tiebreak * tie_break;
        }
    }

    return cost;
}

// ---------------------------------------------------------------------------
// transition_cost_paired_rcpp_cpp
//
// Computes the total negative-log-likelihood transition cost for n_pairs
// paired source/destination track-state vectors.  Unlike the batch version
// (which fixes one source and varies destinations), both source and
// destination change per pair.
//
// Arguments:
//   tdbh0_mat      : NumericMatrix [n_pairs × K] — source DBH per pair/track.
//   tdbh1_mat      : NumericMatrix [n_pairs × K] — destination DBH per pair/track.
//   (all other arguments identical to transition_cost_tracks_bio_batch_rcpp_cpp)
//
// Returns: NumericVector of length n_pairs — total NLL cost per pair.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector transition_cost_paired_rcpp_cpp(
    const Rcpp::NumericMatrix& tdbh0_mat,
    const Rcpp::NumericMatrix& tdbh1_mat,
    double interval_years,
    double mu_const,
    double mu_gamma,
    double sigma0,
    double sigma1,
    double max_shrink,
    double k_shrink,
    double max_growth,
    double max_growth_soft,
    double k_growth,
    bool use_measurement_error,
    double meas_sd1_a,
    double meas_sd1_b,
    double meas_sd2,
    double meas_p_big,
    double h0,
    double beta,
    double recruit_meanlog,
    double recruit_sdlog,
    double recruit_max_dbh,
    double recruit_lambda,
    double eps_tiebreak,
    double hard_penalty
) {
    int K = tdbh0_mat.ncol();
    int n_pairs = tdbh0_mat.nrow();

    if (tdbh1_mat.ncol() != K || tdbh1_mat.nrow() != n_pairs) {
        Rcpp::stop("tdbh0_mat and tdbh1_mat must have identical dimensions");
    }

    Rcpp::NumericVector cost(n_pairs, 0.0);

    // Recruitment probability over interval
    double p_recruit = 1.0 - std::exp(-recruit_lambda * interval_years);
    p_recruit = std::max(1e-12, std::min(1.0 - 1e-12, p_recruit));

    // Lambda: growth mean
    auto mu_growth = [&](double d) -> double {
        if (!std::isfinite(mu_gamma) || mu_gamma == 0.0 || !std::isfinite(d) || d <= 0.0)
            return mu_const;
        return mu_const + mu_gamma * std::log(d);
    };

    // Lambda: measurement error SD
    auto meas_sd1 = [&](double d) -> double {
        return std::max(1e-6, meas_sd1_a * d + meas_sd1_b);
    };

    // Lambda: log-sum-exp
    auto log_sum_exp = [](const std::vector<double>& x) -> double {
        if (x.empty()) return -INFINITY;
        double max_val = *std::max_element(x.begin(), x.end());
        double sum = 0.0;
        for (double val : x) sum += std::exp(val - max_val);
        return max_val + std::log(sum);
    };

    // Lambda: normal log density
    auto dnorm_log = [](double x, double mean, double sd) -> double {
        double diff = x - mean;
        return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * diff * diff / (sd * sd);
    };

    // Lambda: log-normal log density
    auto dlnorm_log = [](double x, double meanlog, double sdlog) -> double {
        if (x <= 0.0) return -INFINITY;
        double log_x = std::log(x);
        double diff = log_x - meanlog;
        return -std::log(x * sdlog * std::sqrt(2.0 * M_PI)) - 0.5 * diff * diff / (sdlog * sdlog);
    };

    // Main loop: iterate over pairs, then over tracks within each pair
    for (int i = 0; i < n_pairs; i++) {
        for (int k = 0; k < K; k++) {
            double d0 = tdbh0_mat(i, k);
            double d1 = tdbh1_mat(i, k);

            if (Rcpp::NumericVector::is_na(d0)) {
                if (Rcpp::NumericVector::is_na(d1)) {
                    // NA -> NA
                    cost[i] -= std::log(1.0 - p_recruit);
                } else {
                    // NA -> DBH (recruitment)
                    bool hard = (!std::isfinite(d1)) || (d1 <= 0.0) || (d1 > recruit_max_dbh);
                    if (hard) {
                        cost[i] += hard_penalty;
                    } else {
                        cost[i] -= std::log(p_recruit) + dlnorm_log(d1, recruit_meanlog, recruit_sdlog);
                    }
                }
            } else if (Rcpp::NumericVector::is_na(d1)) {
                // DBH -> NA (mortality)
                double hazard = h0 * std::exp(beta * d0);
                double p_death = 1.0 - std::exp(-hazard * interval_years);
                p_death = std::max(1e-12, std::min(1.0 - 1e-12, p_death));
                cost[i] -= std::log(p_death);
            } else {
                // DBH -> DBH (growth)
                double g = (d1 - d0) / interval_years;

                // Hard biological constraints
                bool hard = false;
                if (std::isfinite(max_shrink) && (g < max_shrink)) hard = true;
                if (std::isfinite(max_growth) && (g > max_growth)) hard = true;

                if (hard) {
                    cost[i] += hard_penalty;
                } else {
                    double sigma_d = std::max(sigma0 + sigma1 * d0, 1e-6);
                    double mu = mu_growth(d0);

                    if (use_measurement_error) {
                        double s_small0 = meas_sd1(d0);
                        double s_small1 = meas_sd1(d1);
                        double s_big = meas_sd2;
                        double w_small = 1.0 - meas_p_big;
                        double w_big = meas_p_big;

                        std::vector<double> sd_meas_mix = {
                            std::sqrt(s_small0*s_small0 + s_small1*s_small1) / interval_years,
                            std::sqrt(s_small0*s_small0 + s_big*s_big) / interval_years,
                            std::sqrt(s_big*s_big + s_small1*s_small1) / interval_years,
                            std::sqrt(s_big*s_big + s_big*s_big) / interval_years
                        };

                        std::vector<double> wt_meas_mix = {
                            w_small * w_small,
                            w_small * w_big,
                            w_big * w_small,
                            w_big * w_big
                        };

                        std::vector<double> sd_tot(4);
                        for (int j = 0; j < 4; j++)
                            sd_tot[j] = std::sqrt(sigma_d*sigma_d + sd_meas_mix[j]*sd_meas_mix[j]);

                        std::vector<double> ll(4);
                        for (int j = 0; j < 4; j++)
                            ll[j] = std::log(wt_meas_mix[j]) + dnorm_log(g, mu, sd_tot[j]);

                        cost[i] -= log_sum_exp(ll);
                    } else {
                        double diff = g - mu;
                        cost[i] += diff*diff / (2.0 * sigma_d*sigma_d) + std::log(sigma_d) + 0.5 * std::log(2.0 * M_PI);
                    }

                    // Soft penalty for shrinkage
                    if (std::isfinite(k_shrink) && k_shrink > 0.0 && d1 < d0) {
                        double dd = d0 - d1;
                        cost[i] += k_shrink * dd * dd;
                    }

                    // Soft penalty for extreme positive growth
                    if (std::isfinite(max_growth_soft) && std::isfinite(k_growth) && k_growth > 0.0) {
                        double d1_soft_cap = d0 + max_growth_soft * interval_years;
                        if (std::isfinite(d1_soft_cap) && d1 > d1_soft_cap) {
                            double dd = d1 - d1_soft_cap;
                            cost[i] += k_growth * dd * dd;
                        }
                    }
                }
            }
        }

        // Deterministic tie-break
        if (eps_tiebreak > 0.0) {
            // Compute ranks for source
            std::vector<std::pair<double, int>> ranked0;
            for (int k = 0; k < K; k++) {
                double d0 = tdbh0_mat(i, k);
                if (!Rcpp::NumericVector::is_na(d0))
                    ranked0.emplace_back(d0, k);
            }
            std::sort(ranked0.begin(), ranked0.end());
            std::vector<int> r0(K, 0);
            for (size_t j = 0; j < ranked0.size(); j++)
                r0[ranked0[j].second] = j + 1;

            // Compute ranks for destination
            std::vector<std::pair<double, int>> ranked1;
            for (int k = 0; k < K; k++) {
                double d1 = tdbh1_mat(i, k);
                if (!Rcpp::NumericVector::is_na(d1))
                    ranked1.emplace_back(d1, k);
            }
            std::sort(ranked1.begin(), ranked1.end());
            std::vector<int> r1(K, 0);
            for (size_t j = 0; j < ranked1.size(); j++)
                r1[ranked1[j].second] = j + 1;

            // Sum absolute rank differences for tracks observed in both
            double tie_break = 0.0;
            for (int k = 0; k < K; k++) {
                if (!Rcpp::NumericVector::is_na(tdbh0_mat(i, k)) &&
                    !Rcpp::NumericVector::is_na(tdbh1_mat(i, k))) {
                    tie_break += std::abs(r0[k] - r1[k]);
                }
            }
            cost[i] += eps_tiebreak * tie_break;
        }
    }

    return cost;
}

// ---------------------------------------------------------------------------
// derive_phase_prev_batch_rcpp
//
// Vectorized C++ implementation of the R `derive_phase_prev()` logic.
// Checks phase-transition feasibility for every (current_assignment i,
// next_full_state j) pair and returns the indices of feasible pairs plus
// the derived phase_t vector for each pair.
//
// Arguments:
//   tdbh0_mat     : numeric matrix [n_cc  × K] — track DBH at t   per current assignment
//   tdbh1_mat     : numeric matrix [n_next × K] — track DBH at t+1 per next assignment
//   phase_tp1_mat : integer matrix [n_next × K] — phase at t+1 per next full-state
//   resprout_mat  : logical matrix [n_next × K] — resprout flag per next full-state
//                   (pass a 0-row matrix when no resprouts present)
//   prune_hard    : logical — whether to apply hard growth-rate pruning
//   interval_val  : numeric — census interval in years (NA / non-finite → skip pruning)
//   eff_min_grow  : numeric — minimum allowed annual DBH growth (cm/yr)
//   eff_max_grow  : numeric — maximum allowed annual DBH growth (cm/yr)
//   eff_recruit_max: numeric — maximum DBH for a recruit (NA / non-finite → skip)
//
// Returns a list with:
//   from_i    : integer vector (1-based) — current-assignment index for each feasible pair
//   to_j      : integer vector (1-based) — next-full-state index for each feasible pair
//   phase_t   : integer matrix [n_feasible × K] — derived phase at t for each feasible pair
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List derive_phase_prev_batch_rcpp(
    const Rcpp::NumericMatrix&  tdbh0_mat,
    const Rcpp::NumericMatrix&  tdbh1_mat,
    const Rcpp::IntegerMatrix&  phase_tp1_mat,
    const Rcpp::LogicalMatrix&  resprout_mat,
    bool   prune_hard,
    double interval_val,
    double eff_min_grow,
    double eff_max_grow,
    double eff_recruit_max
) {
    const int n_cc   = tdbh0_mat.nrow();
    const int n_next = tdbh1_mat.nrow();
    const int K      = tdbh0_mat.ncol();
    const bool has_resprout = (resprout_mat.nrow() == n_next) && (resprout_mat.ncol() == K);
    const bool do_prune = prune_hard && std::isfinite(interval_val) && interval_val > 0.0;

    // Worst-case pre-allocation (all pairs feasible).
    std::vector<int>   from_i_vec;  from_i_vec.reserve(n_cc * n_next);
    std::vector<int>   to_j_vec;    to_j_vec.reserve(n_cc * n_next);
    // phase_t stored row-major: [n_feasible × K]
    std::vector<int>   phase_t_flat; phase_t_flat.reserve(n_cc * n_next * K);

    std::vector<int>   phase_t_buf(K);

    for (int i = 0; i < n_cc; ++i) {
        for (int j = 0; j < n_next; ++j) {

            // ---- 1. Phase-transition feasibility check ----
            bool feasible = true;

            // Quick global guard: alive_tp1 must have phase 1 and dead_tp1 must not.
            for (int k = 0; k < K && feasible; ++k) {
                double d1      = tdbh1_mat(j, k);
                int    ph_tp1  = phase_tp1_mat(j, k);
                bool   alive1  = !Rcpp::NumericVector::is_na(d1);
                if (alive1 && ph_tp1 != 1) { feasible = false; break; }
                if (!alive1 && ph_tp1 == 1) { feasible = false; break; }
            }
            if (!feasible) continue;

            // Derive phase_t per track
            for (int k = 0; k < K && feasible; ++k) {
                double d0      = tdbh0_mat(i, k);
                double d1      = tdbh1_mat(j, k);
                int    ph_tp1  = phase_tp1_mat(j, k);
                bool   alive0  = !Rcpp::NumericVector::is_na(d0);
                bool   alive1  = !Rcpp::NumericVector::is_na(d1);
                bool   resp    = has_resprout && resprout_mat(j, k);

                if (alive1) {
                    if (resp) {
                        // Resprout barrier: track was unborn at t.
                        if (alive0) { feasible = false; break; }
                        phase_t_buf[k] = 0;
                    } else {
                        phase_t_buf[k] = alive0 ? 1 : 0;
                    }
                } else {
                    // dead or never-born at t+1
                    if (ph_tp1 == 0) {
                        if (alive0) { feasible = false; break; }
                        phase_t_buf[k] = 0;
                    } else if (ph_tp1 == 2) {
                        phase_t_buf[k] = alive0 ? 1 : 2;
                    } else {
                        // ph_tp1 == 1 but !alive1 — impossible (caught above)
                        feasible = false; break;
                    }
                }
            }
            if (!feasible) continue;

            // Final consistency checks
            for (int k = 0; k < K && feasible; ++k) {
                double d0 = tdbh0_mat(i, k);
                bool alive0 = !Rcpp::NumericVector::is_na(d0);
                if (alive0  && phase_t_buf[k] != 1) { feasible = false; }
                if (!alive0 && phase_t_buf[k] == 1) { feasible = false; }
            }
            if (!feasible) continue;

            // ---- 2. Hard growth-rate pruning ----
            if (do_prune) {
                for (int k = 0; k < K && feasible; ++k) {
                    double d0 = tdbh0_mat(i, k);
                    double d1 = tdbh1_mat(j, k);
                    bool alive0 = !Rcpp::NumericVector::is_na(d0);
                    bool alive1 = !Rcpp::NumericVector::is_na(d1);
                    if (alive0 && alive1) {
                        double g = (d1 - d0) / interval_val;
                        if (g < eff_min_grow || g > eff_max_grow) { feasible = false; }
                    } else if (!alive0 && alive1) {
                        if (std::isfinite(eff_recruit_max) && d1 > eff_recruit_max) { feasible = false; }
                    }
                }
            }
            if (!feasible) continue;

            // ---- 3. Record feasible pair ----
            from_i_vec.push_back(i + 1);  // convert to 1-based for R
            to_j_vec.push_back(j + 1);
            for (int k = 0; k < K; ++k) phase_t_flat.push_back(phase_t_buf[k]);
        }
    }

    const int n_feasible = (int)from_i_vec.size();

    Rcpp::IntegerVector r_from(from_i_vec.begin(), from_i_vec.end());
    Rcpp::IntegerVector r_to(to_j_vec.begin(), to_j_vec.end());
    Rcpp::IntegerMatrix r_phase(n_feasible, K);
    for (int r = 0; r < n_feasible; ++r)
        for (int k = 0; k < K; ++k)
            r_phase(r, k) = phase_t_flat[r * K + k];

    return Rcpp::List::create(
        Rcpp::Named("from_i")  = r_from,
        Rcpp::Named("to_j")    = r_to,
        Rcpp::Named("phase_t") = r_phase
    );
}
