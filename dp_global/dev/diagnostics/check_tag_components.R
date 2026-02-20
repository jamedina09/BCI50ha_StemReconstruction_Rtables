source('dp_global/dev/diagnostics/utils.R'); compile_transition_cost(); source('dp_global/R/dp_global_bio.R')

interval <- as.numeric(as.Date('1995-08-21') - as.Date('1990-12-18'))/365.25
params_row <- list(
  Bio_Mu_Growth=0.350960664887803,
  Bio_Gamma_Growth=0.314805038516982,
  Bio_Sigma0_Growth=0.269453148265253,
  Bio_Sigma1_Growth=0.0105915836625149,
  Bio_Max_Shrink=-0.5,
  Bio_K_Shrink=0,
  Bio_Max_Growth=7.5,
  Bio_Max_Growth_Soft=3.97855492092496,
  Bio_K_Growth=0,
  Bio_H0_Mortality=0.00873479569505417,
  Bio_Beta_Mortality=0.00834401015479029,
  Bio_Recruit_Meanlog=-0.250726454828908,
  Bio_Recruit_Sdlog=0.838668739452498,
  Bio_Recruit_MaxDBH_unit=38.4999,
  Bio_Recruitment_lambda=0.0552806029785688
)

cat('\n-- components: NA -> 22.2 (recruit) --\n')
comp <- transition_cost_tracks_bio_components(track_dbh_t = NA, track_dbh_tp1 = 22.2, interval_years = interval, mu_const = params_row$Bio_Mu_Growth, mu_gamma = params_row$Bio_Gamma_Growth, sigma0 = params_row$Bio_Sigma0_Growth, sigma1 = params_row$Bio_Sigma1_Growth, max_shrink = params_row$Bio_Max_Shrink, k_shrink = params_row$Bio_K_Shrink, max_growth = params_row$Bio_Max_Growth, max_growth_soft = params_row$Bio_Max_Growth_Soft, k_growth = params_row$Bio_K_Growth, h0 = params_row$Bio_H0_Mortality, beta = params_row$Bio_Beta_Mortality, recruit_meanlog = params_row$Bio_Recruit_Meanlog, recruit_sdlog = params_row$Bio_Recruit_Sdlog, recruit_max_dbh = params_row$Bio_Recruit_MaxDBH_unit, recruit_lambda = params_row$Bio_Recruitment_lambda)
print(comp$per_track)
cat('per-track total =', comp$total, '\n')

cat('\n-- components: 21.8 -> NA (mortality / missing) --\n')
comp <- transition_cost_tracks_bio_components(track_dbh_t = 21.8, track_dbh_tp1 = NA, interval_years = interval, mu_const = params_row$Bio_Mu_Growth, mu_gamma = params_row$Bio_Gamma_Growth, sigma0 = params_row$Bio_Sigma0_Growth, sigma1 = params_row$Bio_Sigma1_Growth, max_shrink = params_row$Bio_Max_Shrink, k_shrink = params_row$Bio_K_Shrink, max_growth = params_row$Bio_Max_Growth, max_growth_soft = params_row$Bio_Max_Growth_Soft, k_growth = params_row$Bio_K_Growth, h0 = params_row$Bio_H0_Mortality, beta = params_row$Bio_Beta_Mortality, recruit_meanlog = params_row$Bio_Recruit_Meanlog, recruit_sdlog = params_row$Bio_Recruit_Sdlog, recruit_max_dbh = params_row$Bio_Recruit_MaxDBH_unit, recruit_lambda = params_row$Bio_Recruitment_lambda)
print(comp$per_track)
cat('per-track total =', comp$total, '\n')

cat('\n-- components: 21.8 -> 22.2 (growth across 2->4 if present) --\n')
comp <- transition_cost_tracks_bio_components(track_dbh_t = 21.8, track_dbh_tp1 = 22.2, interval_years = as.numeric(as.Date('1995-08-21') - as.Date('1985-06-21'))/365.25, mu_const = params_row$Bio_Mu_Growth, mu_gamma = params_row$Bio_Gamma_Growth, sigma0 = params_row$Bio_Sigma0_Growth, sigma1 = params_row$Bio_Sigma1_Growth, max_shrink = params_row$Bio_Max_Shrink, k_shrink = params_row$Bio_K_Shrink, max_growth = params_row$Bio_Max_Growth, max_growth_soft = params_row$Bio_Max_Growth_Soft, k_growth = params_row$Bio_K_Growth, h0 = params_row$Bio_H0_Mortality, beta = params_row$Bio_Beta_Mortality, recruit_meanlog = params_row$Bio_Recruit_Meanlog, recruit_sdlog = params_row$Bio_Recruit_Sdlog, recruit_max_dbh = params_row$Bio_Recruit_MaxDBH_unit, recruit_lambda = params_row$Bio_Recruitment_lambda)
print(comp$per_track)
cat('per-track total =', comp$total, '\n')
