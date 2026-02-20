source('dp_global/dev/diagnostics/utils.R'); compile_transition_cost(); source('dp_global/src/transition_cost_rcpp.R')

# DBHs at Census 3 and Census 4 for Tag 2747
dbh3 <- c(2.7, 1.4, 1.0)
dbh4 <- c(22.2, 2.7)
# exact interval between dates in dataset
interval <- as.numeric(as.Date('1995-08-21') - as.Date('1990-12-18')) / 365.25

params <- list(mu_const=0.350960664887803,
  mu_gamma=0.314805038516982,
  sigma0=0.269453148265253,
  sigma1=0.0105915836625149,
  max_shrink=-0.5,
  k_shrink=0,
  max_growth=4.0, # lowered hard cap
  max_growth_soft=3.97855492092496,
  k_growth=0,
  use_measurement_error=TRUE,
  meas_sd1_a=0.0062,
  meas_sd1_b=0.0904,
  meas_sd2=4.64,
  meas_p_big=0.05,
  h0=0.00873479569505417,
  beta=0.00834401015479029,
  recruit_meanlog=-0.250726454828908,
  recruit_sdlog=0.838668739452498,
  recruit_max_dbh=38.4999,
  recruit_lambda=0.0552806029785688,
  eps_tiebreak=1e-6)

for (i in seq_along(dbh3)) {
  cat(sprintf("--- From Census3 DBH=%.2f (interval=%.2f yrs) ---\n", dbh3[i], interval))
  costs <- transition_cost_tracks_bio_batch_rcpp(track_dbh_t = dbh3[i], track_dbh_tp1 = as.list(as.data.frame(t(dbh4))), interval_years = interval, mu_const = params$mu_const, mu_gamma = params$mu_gamma, sigma0 = params$sigma0, sigma1 = params$sigma1, max_shrink = params$max_shrink, k_shrink = params$k_shrink, max_growth = params$max_growth, max_growth_soft = params$max_growth_soft, k_growth = params$k_growth, use_measurement_error = params$use_measurement_error, meas_sd1_a = params$meas_sd1_a, meas_sd1_b = params$meas_sd1_b, meas_sd2 = params$meas_sd2, meas_p_big = params$meas_p_big, h0 = params$h0, beta = params$beta, recruit_meanlog = params$recruit_meanlog, recruit_sdlog = params$recruit_sdlog, recruit_max_dbh = params$recruit_max_dbh, recruit_lambda = params$recruit_lambda, eps_tiebreak = params$eps_tiebreak)
  print(data.frame(to_dbh=dbh4, cost=costs))
}
