# Parameter tuning (2D) to prevent DBH jump 2.7 -> 22.2
library(data.table)
library(ggplot2)
source('dp_global/R/dp_global_bio.R')
source('dp_global/src/transition_cost_rcpp.R')

# Baseline args (from earlier runs)
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

# Observed pair and interval
d0 <- 2.7; d1 <- 22.2
interval_years <- as.numeric(as.Date('1995-08-21') - as.Date('1990-12-18'))/365.25

# Baseline recruit cost (NA -> 22.2)
recruit_comp <- transition_cost_tracks_bio_components(track_dbh_t = NA, track_dbh_tp1 = d1, interval_years = interval_years, mu_const = base_args$mu_const, mu_gamma = base_args$mu_gamma, sigma0 = base_args$sigma0, sigma1 = base_args$sigma1, max_shrink = base_args$max_shrink, k_shrink = base_args$k_shrink, max_growth = base_args$max_growth, max_growth_soft = base_args$max_growth_soft, k_growth = base_args$k_growth, h0 = base_args$h0, beta = base_args$beta, recruit_meanlog = base_args$recruit_meanlog, recruit_sdlog = base_args$recruit_sdlog, recruit_max_dbh = base_args$recruit_max_dbh, recruit_lambda = base_args$recruit_lambda)
recruit_total <- recruit_comp$total

# 2D grid: max_growth_soft x k_growth
mg_soft_vals <- seq(1.5, 4.5, length.out = 121)
kg_vals <- seq(0, 12, length.out = 121)

grid <- CJ(max_growth_soft = mg_soft_vals, k_growth = kg_vals)
setorder(grid, max_growth_soft, k_growth)

# compute DBH->DBH total for each pair
calc_row <- function(mg_soft, kg) {
  comp <- transition_cost_tracks_bio_components(track_dbh_t = d0, track_dbh_tp1 = d1, interval_years = interval_years, mu_const = base_args$mu_const, mu_gamma = base_args$mu_gamma, sigma0 = base_args$sigma0, sigma1 = base_args$sigma1, max_shrink = base_args$max_shrink, k_shrink = base_args$k_shrink, max_growth = base_args$max_growth, max_growth_soft = mg_soft, k_growth = kg, h0 = base_args$h0, beta = base_args$beta, recruit_meanlog = base_args$recruit_meanlog, recruit_sdlog = base_args$recruit_sdlog, recruit_max_dbh = base_args$recruit_max_dbh, recruit_lambda = base_args$recruit_lambda)
  return(as.numeric(comp$per_track$total_track))
}

# Vectorized loop (fast enough for 14k evaluations)
grid[, dbh_dbh_cost := mapply(calc_row, max_growth_soft, k_growth)]

# Flag combos where DBH->DBH cost >= recruit_total
grid[, ok := dbh_dbh_cost >= recruit_total]

# Find minimal k_growth for each max_growth_soft such that ok==TRUE
min_kg <- grid[ok == TRUE, .(min_k = min(k_growth)), by = max_growth_soft]
if (nrow(min_kg) == 0L) {
  message('No (mg_soft, k_growth) within grid makes DBH->DBH >= recruit cost')
} else {
  # pick recommended point: smallest k overall (minimise penalty) where ok==TRUE
  rec <- min_kg[which.min(min_kg$min_k)]
  rec_pair <- list(max_growth_soft = rec$max_growth_soft, k_growth = rec$min_k)
  message('Recommended (soft cap, k_growth) ~', rec_pair$max_growth_soft, ',', rec_pair$k_growth)
}

# Save CSV and contour plot
fwrite(grid, 'dp_global/dev/diagnostics/tune_growth_grid_tag_2747.csv')

p <- ggplot(grid, aes(x = max_growth_soft, y = k_growth, z = dbh_dbh_cost)) +
  geom_raster(aes(fill = dbh_dbh_cost)) +
  geom_contour(color = 'white', breaks = c(recruit_total)) +
  scale_fill_viridis_c(name = 'DBH->DBH cost') +
  geom_point(data = data.table(max_growth_soft = base_args$max_growth_soft, k_growth = base_args$k_growth), aes(x = max_growth_soft, y = k_growth), color = 'red', size = 3) +
  labs(title = 'Cost(DBH->DBH) for 2.7 -> 22.2; contour = recruit cost', x = 'max_growth_soft (cm/yr)', y = 'k_growth (1/cm^2)') +
  theme_minimal()

ggsave('dp_global/dev/diagnostics/tune_growth_tag_2747_contour.pdf', p, width = 8, height = 6)

# Write recommendation CSV
if (exists('rec_pair')) {
  reco_dt <- data.table(max_growth_soft = rec_pair$max_growth_soft, k_growth = rec_pair$k_growth)
  fwrite(reco_dt, 'dp_global/dev/diagnostics/tune_growth_tag_2747_recommendation.csv')
} else {
  fwrite(data.table(), 'dp_global/dev/diagnostics/tune_growth_tag_2747_recommendation.csv')
}

cat('Grid + contour + recommendation written to dp_global/dev/diagnostics/\n')
