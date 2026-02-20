# compare_costs.R
# Compare transition cost from 2.7 -> 22.2 for various soft penalty parameters
library(data.table)
source('dp_global/R/dp_global_bio.R')
source('dp_global/src/transition_cost_rcpp.R')

params <- list(
    mu_const = 0.350960664887803,
    mu_gamma = 0.314805038516982,
    sigma0 = 0.269453148265253,
    sigma1 = 0.0105915836625149,
    max_shrink = -0.5,
    k_shrink = 0,
    max_growth = 7.5,
    use_measurement_error = TRUE,
    meas_sd1_a = 0.0062,
    meas_sd1_b = 0.0904,
    meas_sd2 = 4.64,
    meas_p_big = 0.05,
    h0 = 0.00873479569505417,
    beta = 0.00834401015479029,
    recruit_meanlog = -0.250726454828908,
    recruit_sdlog = 0.838668739452498,
    recruit_max_dbh = 38.4999,
    recruit_lambda = 0.0552806029785688
)
interval_years <- as.numeric(as.Date("1995-08-21") - as.Date("1990-12-18")) / 365.25

pairs <- data.table(
    max_growth_soft = c(1.5, 3.97855492092496),
    k_growth = c(0, 0)
)

# evaluate cost for each row separately
results <- pairs[, .(
    cost = transition_cost_tracks_bio_components(
        track_dbh_t = 2.7, track_dbh_tp1 = 22.2, interval_years = interval_years,
        mu_const = params$mu_const, mu_gamma = params$mu_gamma,
        sigma0 = params$sigma0, sigma1 = params$sigma1,
        max_shrink = params$max_shrink, k_shrink = params$k_shrink,
        max_growth = params$max_growth, max_growth_soft = max_growth_soft,
        k_growth = k_growth, , use_measurement_error = params$use_measurement_error,
        meas_sd1_a = params$meas_sd1_a, meas_sd1_b = params$meas_sd1_b,
        meas_sd2 = params$meas_sd2, meas_p_big = params$meas_p_big,
        h0 = params$h0, beta = params$beta,
        recruit_meanlog = params$recruit_meanlog, recruit_sdlog = params$recruit_sdlog,
        recruit_max_dbh = params$recruit_max_dbh, recruit_lambda = params$recruit_lambda,
        eps_tiebreak = 1e-6
    )$total
), by = .(max_growth_soft, k_growth)]
print(results)

# Also compute via Rcpp batch for comparison
for (i in seq_len(nrow(pairs))) {
    cat(sprintf("Evaluating batch for max_growth_soft=%.3f, k_growth=%.3f...\n", pairs$max_growth_soft[i], pairs$k_growth[i]))
    batch <- transition_cost_tracks_bio_batch_rcpp(
        track_dbh_t = 2.7, track_dbh_tp1 = list(22.2), interval_years = interval_years,
        mu_const = params$mu_const, mu_gamma = params$mu_gamma,
        sigma0 = params$sigma0, sigma1 = params$sigma1,
        max_shrink = params$max_shrink, k_shrink = params$k_shrink,
        max_growth = params$max_growth, max_growth_soft = pairs$max_growth_soft[i],
        k_growth = pairs$k_growth[i], use_measurement_error = params$use_measurement_error,
        meas_sd1_a = params$meas_sd1_a, meas_sd1_b = params$meas_sd1_b,
        meas_sd2 = params$meas_sd2, meas_p_big = params$meas_p_big,
        h0 = params$h0, beta = params$beta,
        recruit_meanlog = params$recruit_meanlog, recruit_sdlog = params$recruit_sdlog,
        recruit_max_dbh = params$recruit_max_dbh, recruit_lambda = params$recruit_lambda,
        eps_tiebreak = 1e-6
    )
    print(batch)
}
