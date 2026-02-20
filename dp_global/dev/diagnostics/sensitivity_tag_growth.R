# Sensitivity sweep for the suspicious 2.7 -> 22.2 jump (Tag 2747)
# Sweeps: max_growth and k_growth

library(data.table)
source('dp_global/R/dp_global_bio.R')
source('dp_global/R/sensitivity_transition_cost_bio.R')
source('dp_global/src/transition_cost_rcpp.R')

# Baseline args taken from the earlier DP run for tag 2747 (approx.)
base_args <- list(
  mu_const = 0.350960664887803,
  mu_gamma = 0.314805038516982,
  sigma0 = 0.269453148265253,
  sigma1 = 0.0105915836625149,
  h0 = 0.00873479569505417,
  beta = 0.00834401015479029,
  recruit_meanlog = -0.250726454828908,
  recruit_sdlog = 0.838668739452498,
  recruit_max_dbh = 38.4999,
  recruit_lambda = 0.0552806029785688,
  max_shrink = -0.5,
  k_shrink = 0,
  max_growth = 7.5,
  max_growth_soft = 3.97855492092496,
  k_growth = 0,
  use_measurement_error = TRUE,
  meas_sd1_a = 0.0062,
  meas_sd1_b = 0.0904,
  meas_sd2 = 4.64,
  meas_p_big = 0.05
)

# Scenario: DBH at census t = 2.7  --> t+1 = 22.2 (interval years from data)
interval_years <- as.numeric(as.Date('1995-08-21') - as.Date('1990-12-18'))/365.25
scenarios <- list(tag2747_jump = list(name = 'Tag2747 2.7->22.2', t = c(2.7), tp1 = c(22.2)))

# Build grids: focus on ranges around the problematic growth rate (~4.18 cm/yr)
# max_growth: 0..7 ; k_growth: 0..2
mg_vals <- seq(0, 10, length.out = 201)
kg_vals <- seq(0, 10, length.out = 201)

dts <- list()
# Sweep max_growth while keeping k_growth fixed at 0
dt_mg <- sweep_transition_cost(
    track_dbh_t = scenarios$tag2747_jump$t, track_dbh_tp1 = scenarios$tag2747_jump$tp1,
    interval_years = interval_years, base_args = base_args, param = "max_growth", values = mg_vals,
    eps_tiebreak = 1e-6, hard_penalty = 1e6
)
# Sweep k_growth while keeping max_growth at baseline
dt_kg <- sweep_transition_cost(
    track_dbh_t = scenarios$tag2747_jump$t, track_dbh_tp1 = scenarios$tag2747_jump$tp1,
    interval_years = interval_years, base_args = base_args, param = "k_growth", values = kg_vals,
    eps_tiebreak = 1e-6, hard_penalty = 1e6
)

# out_pdf <- 'dp_global/dev/diagnostics/sensitivity_tag_2747_growth.pdf'
plot_sweep_components(dt_mg, title = 'Sensitivity: max_growth (2.7->22.2)', subtitle = sprintf('interval=%.3f yrs', interval_years))
plot_sweep_components(dt_kg, title = 'Sensitivity: k_growth (2.7->22.2)', subtitle = sprintf('interval=%.3f yrs', interval_years))

# Save small CSV summaries: where cost crosses the recruit cost (baseline recruit cost ~13.30453)
recruit_cost <- transition_cost_tracks_bio_components(
    track_dbh_t = NA, track_dbh_tp1 = 22.2,
    interval_years = interval_years, mu_const = base_args$mu_const, mu_gamma = base_args$mu_gamma, sigma0 = base_args$sigma0, sigma1 = base_args$sigma1, max_shrink = base_args$max_shrink, k_shrink = base_args$k_shrink, max_growth = base_args$max_growth, max_growth_soft = base_args$max_growth_soft, k_growth = base_args$k_growth, h0 = base_args$h0, beta = base_args$beta, recruit_meanlog = base_args$recruit_meanlog, recruit_sdlog = base_args$recruit_sdlog, recruit_max_dbh = base_args$recruit_max_dbh, recruit_lambda = base_args$recruit_lambda
)
recruit_total <- recruit_cost$total

# find minimal parameter values that make DBH->DBH cost >= recruit_total
mg_threshold <- min(dt_mg$value[dt_mg$total >= recruit_total], na.rm = TRUE)
kg_threshold <- min(dt_kg$value[dt_kg$total >= recruit_total], na.rm = TRUE)

fsummary <- data.table(param = c('recruit_total', 'mg_threshold', 'kg_threshold'), value = c(recruit_total, mg_threshold, kg_threshold))
# fwrite(dt_mg, 'dp_global/dev/diagnostics/sensitivity_tag_2747_max_growth_sweep.csv')
# fwrite(dt_kg, 'dp_global/dev/diagnostics/sensitivity_tag_2747_k_growth_sweep.csv')
# fwrite(fsummary, 'dp_global/dev/diagnostics/sensitivity_tag_2747_summary.csv')

cat('Wrote sweep CSVs and summary to dp_global/dev/diagnostics/ (sensitivity_tag_2747_*)\n')
cat('Suggested thresholds (approx):\n')
print(fsummary)
