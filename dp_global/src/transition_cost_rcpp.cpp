#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

// ---------------------------------------------------------------------------
// transition_cost_tracks_bio_batch_rcpp_cpp
//
// Transition cost (negative log-likelihood) from one track-state vector at
// census t to each of n_batch candidate track-state vectors at census t+1.
// DBH is in cm, growth in cm/year and the interval in years.
//
// Arguments:
//   track_dbh_t    : NumericVector of length K — DBH per track at t
//                    (NA = no stem on the track).
//   mat_tp1        : NumericMatrix [n_batch × K] — candidate DBH per track at
//                    t+1; each row is one candidate assignment.
//   interval_years : Length of the census interval in years.
//   mu_const       : Intercept of the growth mean (cm/year).
//   mu_gamma       : Slope of the growth mean on log(DBH):
//                      mean(D) = mu_const + mu_gamma * log(D).
//                    The mean is mu_const when mu_gamma is 0 or non-finite, or
//                    D is not positive.
//   sigma0, sigma1 : Intercept and slope of the process SD model (cm/year):
//                      SD(D) = sigma0 + sigma1 * D, floored at 1e-6.
//   max_shrink     : Hard lower bound on annual growth (cm/year); transitions
//                    below this rate receive hard_penalty. Not applied when
//                    non-finite.
//   k_shrink       : Weight of the soft shrinkage penalty: a DBH decrease adds
//                    k_shrink * (D_t − D_t+1)^2 (difference in cm). Applied
//                    when k_shrink > 0.
//   max_growth     : Hard upper bound on annual growth (cm/year); transitions
//                    above this rate receive hard_penalty. Not applied when
//                    non-finite.
//   max_growth_soft: Soft growth cap (cm/year). A DBH at t+1 above
//                    D_t + max_growth_soft * interval_years adds
//                    k_growth * excess^2 (excess in cm).
//   k_growth       : Weight of the soft growth penalty. Applied when
//                    k_growth > 0 and max_growth_soft is finite.
//   use_measurement_error : When TRUE, the growth likelihood is a 4-component
//                    mixture (a small or a large measurement error at t, times
//                    a small or a large one at t+1, each added to the process
//                    variance) rather than a single Gaussian.
//   meas_sd1_a, meas_sd1_b : Parameters of the small-error SD model (cm):
//                              SD1(D) = meas_sd1_a * D + meas_sd1_b,
//                            floored at 1e-6.
//   meas_sd2       : SD of a large measurement error (cm).
//   meas_p_big     : Probability that one measurement carries a large error.
//   h0, beta       : Hazard model parameters for mortality probability:
//                      p_die = 1 − exp(−h0 * exp(beta * D) * interval_years).
//   recruit_meanlog, recruit_sdlog : Log-normal parameters for recruit DBH.
//   recruit_max_dbh : Hard cap on recruit DBH; exceeding this gives hard_penalty.
//   recruit_lambda : Recruitment rate per year of a track with no stem at t:
//                      p_recruit = 1 − exp(−recruit_lambda * interval_years).
//   eps_tiebreak   : Weight of the deterministic rank term added to each cost
//                    to break ties in a reproducible way (not added when
//                    eps_tiebreak <= 0).
//   hard_penalty   : Cost assigned to biologically impossible transitions.
//                    No default here; the R wrappers in transition_cost_rcpp.R
//                    default to 1e6.
//   round_t, round_tp1 : TRUE when census t (t+1) recorded DBH in classes,
//                    rounded down (default false).
//   round_max_dbh  : Only DBH below this value was rounded (cm, default 5.5).
//   round_width    : Width of the rounding classes (cm, default 0.5).
//
// Cost per track, summed over the K tracks of each candidate row:
//   NA → NA   : −log(1 − p_recruit) — no stem on the track at either census.
//   NA → DBH  : −log(p_recruit) − log dlnorm(DBH | recruit_meanlog,
//               recruit_sdlog); hard_penalty instead when DBH is not positive,
//               not finite or above recruit_max_dbh.
//   DBH → NA  : −log(p_die) — the stem is gone at t+1 (mortality).
//   DBH → DBH : hard_penalty when the annual growth
//               g = (D_t+1 − D_t) / interval_years is outside
//               [max_shrink, max_growth]; otherwise −log of the likelihood of
//               g (Gaussian or 4-component mixture) plus the soft penalties.
//               At a census flagged as rounded, see the comment in the code.
// p_recruit and p_die are clamped to [1e-12, 1 − 1e-12].
//
// Returns: NumericVector of length n_batch — the cost of each candidate summed
//          across tracks, plus eps_tiebreak × Σ |rank_t − rank_t+1| over the
//          tracks observed at both censuses, where a rank is the position of
//          the track's DBH among the observed DBHs of that census.
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
    double hard_penalty,
    bool round_t = false,
    bool round_tp1 = false,
    double round_max_dbh = 5.5,
    double round_width = 0.5
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

    auto mu_growth = [&](double d) -> double {
        if (!std::isfinite(mu_gamma) || mu_gamma == 0.0 || !std::isfinite(d) || d <= 0.0) {
            return mu_const;
        }
        return mu_const + mu_gamma * std::log(d);
    };

    // SD of a small measurement error at DBH d
    auto meas_sd1 = [&](double d) -> double {
        return std::max(1e-6, meas_sd1_a * d + meas_sd1_b);
    };

    auto log_sum_exp = [](const std::vector<double>& x) -> double {
        if (x.empty()) return -INFINITY;
        double max_val = *std::max_element(x.begin(), x.end());
        double sum = 0.0;
        for (double val : x) {
            sum += std::exp(val - max_val);
        }
        return max_val + std::log(sum);
    };

    auto dnorm_log = [](double x, double mean, double sd) -> double {
        double diff = x - mean;
        return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * diff * diff / (sd * sd);
    };

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

            // DBH recorded in classes of width round_width and rounded down
            // (round_t / round_tp1: the census is flagged; only stems below
            // round_max_dbh were rounded, e.g. BCI 1982 and 1985 below 55 mm).
            // The true DBH lies in [d, d + round_width): the expected growth moves
            // to the class mid-points and its variance gains round_width^2 / 12
            // per rounded measurement. Without flagged censuses nothing changes.
            bool r0 = round_t && d0 < round_max_dbh;
            bool r1 = round_tp1 && d1 < round_max_dbh;
            double g_mid = g + ((r1 ? 0.5 * round_width : 0.0) - (r0 ? 0.5 * round_width : 0.0)) / interval_years;
            double var_round = ((r0 ? 1.0 : 0.0) + (r1 ? 1.0 : 0.0)) * round_width * round_width / 12.0 /
                               (interval_years * interval_years);
            double d0_mid = d0 + (r0 ? 0.5 * round_width : 0.0);

            // Hard biological constraints (on the measured DBHs: rounding changes
            // how likely a link is, never which links are allowed)
            bool hard = false;
            if (std::isfinite(max_shrink) && (g < max_shrink)) hard = true;
            if (std::isfinite(max_growth) && (g > max_growth)) hard = true;

            if (hard) {
                cost[i] += hard_penalty;
                continue;
            }

            // Growth SD and mean at the DBH at t (class mid-point when rounded)
            double sigma_d = sigma0 + sigma1 * d0_mid;
            sigma_d = std::max(sigma_d, 1e-6);
            double mu = mu_growth(d0_mid);

            if (use_measurement_error) {
                double s_small0 = meas_sd1(d0);
                double s_small1 = meas_sd1(d1);
                double s_big = meas_sd2;
                double w_small = 1.0 - meas_p_big;
                double w_big = meas_p_big;

                // Measurement-error SD of the annual growth for the four error
                // combinations (t, t+1): small/small, small/large, large/small,
                // large/large; wt_meas_mix holds their probabilities
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
                    sd_tot[j] = (var_round > 0.0) ? std::sqrt(sigma_d*sigma_d + var_round + sd_meas_mix[j]*sd_meas_mix[j])
                                             : std::sqrt(sigma_d*sigma_d + sd_meas_mix[j]*sd_meas_mix[j]);
                }

                std::vector<double> ll(4);
                for (int j = 0; j < 4; j++) {
                    ll[j] = std::log(wt_meas_mix[j]) + dnorm_log(g_mid, mu, sd_tot[j]);
                }

                cost[i] -= log_sum_exp(ll);
            } else {
                // Gaussian growth likelihood (process SD plus rounding, if any)
                double sd_g = (var_round > 0.0) ? std::sqrt(sigma_d*sigma_d + var_round) : sigma_d;
                double diff = g_mid - mu;
                cost[i] += diff*diff / (2.0 * sd_g*sd_g) + std::log(sd_g) + 0.5 * std::log(2.0 * M_PI);
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

    // Deterministic tie-break: among candidates of equal cost, prefer the one
    // that keeps the size order of the stems. Ranks order the observed DBHs
    // within a census (1 = smallest; equal DBHs in track order).
    if (eps_tiebreak > 0.0) {
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
// destination change per pair. The per-track cost and the tie-break term are
// those of transition_cost_tracks_bio_batch_rcpp_cpp.
//
// Arguments:
//   tdbh0_mat      : NumericMatrix [n_pairs × K] — source DBH per pair/track
//                    (census t).
//   tdbh1_mat      : NumericMatrix [n_pairs × K] — destination DBH per
//                    pair/track (census t+1); same dimensions as tdbh0_mat.
//   (all other arguments identical to transition_cost_tracks_bio_batch_rcpp_cpp)
//
// Returns: NumericVector of length n_pairs — total NLL cost per pair,
//          including the tie-break term.
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
    double hard_penalty,
    bool round_t = false,
    bool round_tp1 = false,
    double round_max_dbh = 5.5,
    double round_width = 0.5
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

    auto mu_growth = [&](double d) -> double {
        if (!std::isfinite(mu_gamma) || mu_gamma == 0.0 || !std::isfinite(d) || d <= 0.0)
            return mu_const;
        return mu_const + mu_gamma * std::log(d);
    };

    // SD of a small measurement error at DBH d
    auto meas_sd1 = [&](double d) -> double {
        return std::max(1e-6, meas_sd1_a * d + meas_sd1_b);
    };

    auto log_sum_exp = [](const std::vector<double>& x) -> double {
        if (x.empty()) return -INFINITY;
        double max_val = *std::max_element(x.begin(), x.end());
        double sum = 0.0;
        for (double val : x) sum += std::exp(val - max_val);
        return max_val + std::log(sum);
    };

    auto dnorm_log = [](double x, double mean, double sd) -> double {
        double diff = x - mean;
        return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * diff * diff / (sd * sd);
    };

    auto dlnorm_log = [](double x, double meanlog, double sdlog) -> double {
        if (x <= 0.0) return -INFINITY;
        double log_x = std::log(x);
        double diff = log_x - meanlog;
        return -std::log(x * sdlog * std::sqrt(2.0 * M_PI)) - 0.5 * diff * diff / (sdlog * sdlog);
    };

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

                // DBH recorded in classes of width round_width and rounded down
                // (round_t / round_tp1: the census is flagged; only stems below
                // round_max_dbh were rounded, e.g. BCI 1982 and 1985 below 55 mm).
                // The true DBH lies in [d, d + round_width): the expected growth moves
                // to the class mid-points and its variance gains round_width^2 / 12
                // per rounded measurement. Without flagged censuses nothing changes.
                bool r0 = round_t && d0 < round_max_dbh;
                bool r1 = round_tp1 && d1 < round_max_dbh;
                double g_mid = g + ((r1 ? 0.5 * round_width : 0.0) - (r0 ? 0.5 * round_width : 0.0)) / interval_years;
                double var_round = ((r0 ? 1.0 : 0.0) + (r1 ? 1.0 : 0.0)) * round_width * round_width / 12.0 /
                                   (interval_years * interval_years);
                double d0_mid = d0 + (r0 ? 0.5 * round_width : 0.0);

                // Hard biological constraints (on the measured DBHs: rounding changes
                // how likely a link is, never which links are allowed)
                bool hard = false;
                if (std::isfinite(max_shrink) && (g < max_shrink)) hard = true;
                if (std::isfinite(max_growth) && (g > max_growth)) hard = true;

                if (hard) {
                    cost[i] += hard_penalty;
                } else {
                    double sigma_d = std::max(sigma0 + sigma1 * d0_mid, 1e-6);
                    double mu = mu_growth(d0_mid);

                    if (use_measurement_error) {
                        double s_small0 = meas_sd1(d0);
                        double s_small1 = meas_sd1(d1);
                        double s_big = meas_sd2;
                        double w_small = 1.0 - meas_p_big;
                        double w_big = meas_p_big;

                        // Error combinations (t, t+1): small/small, small/large,
                        // large/small, large/large
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
                            sd_tot[j] = (var_round > 0.0) ? std::sqrt(sigma_d*sigma_d + var_round + sd_meas_mix[j]*sd_meas_mix[j])
                                             : std::sqrt(sigma_d*sigma_d + sd_meas_mix[j]*sd_meas_mix[j]);

                        std::vector<double> ll(4);
                        for (int j = 0; j < 4; j++)
                            ll[j] = std::log(wt_meas_mix[j]) + dnorm_log(g_mid, mu, sd_tot[j]);

                        cost[i] -= log_sum_exp(ll);
                    } else {
                        double sd_g = (var_round > 0.0) ? std::sqrt(sigma_d*sigma_d + var_round) : sigma_d;
                        double diff = g_mid - mu;
                        cost[i] += diff*diff / (2.0 * sd_g*sd_g) + std::log(sd_g) + 0.5 * std::log(2.0 * M_PI);
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

        // Deterministic tie-break (rank term of the batch version)
        if (eps_tiebreak > 0.0) {
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
// Checks phase-transition feasibility for every (current_assignment i,
// next_full_state j) pair and returns the indices of feasible pairs plus
// the derived phase_t vector for each pair. Used by the DP backward pass,
// which knows the phases at census t+1 and derives those at census t.
//
// Phase of a track at a census: 0 = stem not yet recruited, 1 = alive (the
// track holds a DBH), 2 = stem dead, or track not used. Rules per track:
//   - at t+1, phase 1 exactly when the track holds a DBH;
//   - DBH at t+1: phase 1 at t when the track holds a DBH at t, else phase 0
//     (recruited in the interval); a resprout-flagged DBH at t+1 must
//     continue a DBH at t;
//   - phase 0 at t+1: no DBH at t, phase 0 at t;
//   - phase 2 at t+1: phase 1 at t when the track holds a DBH at t (death in
//     the interval), else phase 2.
//
// Arguments:
//   tdbh0_mat     : numeric matrix [n_cc  × K] — track DBH at t   per current assignment
//   tdbh1_mat     : numeric matrix [n_next × K] — track DBH at t+1 per next assignment
//   phase_tp1_mat : integer matrix [n_next × K] — phase at t+1 per next full-state
//   resprout_mat  : logical matrix [n_next × K] — resprout flag per next full-state
//                   (a matrix of any other shape, e.g. 0 rows, means no resprouts)
//   prune_hard    : logical — whether to apply hard growth-rate pruning
//   interval_val  : numeric — census interval in years (non-finite or ≤ 0 → skip pruning)
//   eff_min_grow  : numeric — minimum allowed annual DBH growth (cm/yr)
//   eff_max_grow  : numeric — maximum allowed annual DBH growth (cm/yr)
//   eff_recruit_max: numeric — maximum DBH for a recruit, applied with the
//                   growth pruning (non-finite → no limit)
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

            // At t+1 a track has phase 1 exactly when it holds a DBH.
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
                        // A resprout-flagged DBH at t+1 stays on the track of the stem
                        // it came from: the track must be alive at t (no recruitment).
                        if (!alive0) { feasible = false; break; }
                        phase_t_buf[k] = 1;
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

            // Phase 1 at t exactly on the tracks that hold a DBH at t
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

// ---------------------------------------------------------------------------
// hungarian_min_rcpp
//
// Exact minimum-cost assignment of a square cost matrix (Kuhn-Munkres with
// row/column potentials, O(n^3)). Used by the probabilistic matcher to solve
// each perturbed census pair exactly (greedy_assignment_gumbel(), birth-death
// mode). Forbidden cells must be passed as large finite costs, not Inf.
//
// Arguments:
//   cost : NumericMatrix [n x n] — cost of assigning row i to column j.
//
// Returns: IntegerVector of length n — the 1-based column assigned to each row.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::IntegerVector hungarian_min_rcpp(const Rcpp::NumericMatrix& cost) {
    const int n = cost.nrow();
    if (cost.ncol() != n) Rcpp::stop("hungarian_min_rcpp: cost must be a square matrix");
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (!std::isfinite(cost(i, j))) Rcpp::stop("hungarian_min_rcpp: costs must be finite");
    const double INF = std::numeric_limits<double>::infinity();
    // 1-based potentials u (rows), v (columns); p[j] = row matched to column j
    std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0), minv(n + 1);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);
    std::vector<char> used(n + 1);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::fill(minv.begin(), minv.end(), INF);
        std::fill(used.begin(), used.end(), 0);
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                const double cur = cost(i0 - 1, j - 1) - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do { const int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0 != 0);
    }
    Rcpp::IntegerVector ans(n);
    for (int j = 1; j <= n; ++j) if (p[j] > 0) ans[p[j] - 1] = j;
    return ans;
}
