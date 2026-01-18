#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

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
