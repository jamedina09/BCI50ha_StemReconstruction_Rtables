# Utility helpers for diagnostics scripts
suppressPackageStartupMessages({
  library(data.table)
})

# Load dataset CSV (path defaults to diagnostics diagnostics_dataset.csv if present)
load_dataset <- function(path = NULL) {
  if (is.null(path)) {
    candidate <- file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv')
    if (file.exists(candidate)) path <- candidate
    else stop('No dataset path provided and default diagnostics dataset not found at ', candidate)
  }
  if (!file.exists(path)) stop('Missing input CSV at ', path)
  Dt <- fread(path)
  if ('ExactDate' %in% names(Dt)) Dt[, ExactDate := as.Date(ExactDate, format='%m/%d/%y')]
  Dt[, Tag := as.integer(Tag)]; Dt[, CensusID := as.integer(CensusID)]; Dt[, DBH := as.numeric(DBH)]
  # coerce Bio_* columns to numeric to avoid integer truncation warnings when assigning doubles
  bio_cols <- grep('^Bio_', names(Dt), value = TRUE)
  for (bc in bio_cols) Dt[, (bc) := as.numeric(get(bc))]
  if ('TrueStemID' %in% names(Dt)) Dt[, TrueStemID := as.integer(TrueStemID)]
  return(Dt)
}

# Safe compile/sourcing of transition cost C++/R wrappers
compile_transition_cost <- function() {
  cpp_path <- file.path('dp_global','src','transition_cost_rcpp.cpp')
  r_wrap <- file.path('dp_global','src','transition_cost_rcpp.R')
  tryCatch({
    if (file.exists(cpp_path)) Rcpp::sourceCpp(cpp_path)
    if (file.exists(r_wrap)) sys.source(r_wrap, envir = globalenv())
    TRUE
  }, error = function(e) { message('[compile] failed: ', e$message); FALSE })
}

# Load DP function and helpers into a new env and return the match function
load_dp_fn <- function() {
  safe_source <- function(p) if (file.exists(p)) sys.source(p, envir = globalenv())
  safe_source(file.path('dp_global','R','dp_global_utils.R'))
  safe_source(file.path('dp_global','R','dp_global_states.R'))
  safe_source(file.path('dp_global','R','dp_global_matchers.R'))
  safe_source(file.path('dp_global','R','dp_global_bio.R'))
  safe_source(file.path('dp_global','R','dp_global_diag.R'))

  env_test <- new.env(parent = globalenv())
  sys.source(file.path('dp_global','R','dp_global_dp.R'), envir = env_test)
  fn <- get('match_stems_dp_global_backward_marginals_batch', envir = env_test)
  fn_env <- environment(fn)
  # copy wrappers if available
  if (exists('transition_cost_tracks_bio_batch_rcpp_cpp', mode='function', inherits = TRUE)) assign('transition_cost_tracks_bio_batch_rcpp_cpp', transition_cost_tracks_bio_batch_rcpp_cpp, envir = fn_env)
  if (exists('transition_cost_tracks_bio_batch_rcpp', mode='function', inherits = TRUE)) assign('transition_cost_tracks_bio_batch_rcpp', transition_cost_tracks_bio_batch_rcpp, envir = fn_env)
  helpers <- c('count_injective_states','enumerate_states_injective','state_key','state_to_track_dbh','add_constraint_violation')
  for (h in helpers) if (!exists(h, envir = fn_env, mode='function') && exists(h, mode='function', inherits = TRUE)) assign(h, get(h, inherits = TRUE), envir = fn_env)
  return(fn)
}

# Run DP on a data.table with common defaults
run_dp_on_dt <- function(Dt, anchor_start = 7, max_tracks = 100, slack_tracks = 0, posterior_top_k = 3, verbose = FALSE, min_growth = -2.5, max_growth = 7.5) {
  fn <- load_dp_fn()
  res <- fn(copy(Dt), anchor_start = anchor_start, min_growth = min_growth, max_growth = max_growth, max_tracks = max_tracks, slack_tracks = slack_tracks, posterior_top_k = posterior_top_k, verbose = verbose)
  setDT(res)
  return(res)
}

# Compute cost components between two DBHs or NA using available R function
compute_cost_components <- function(track_dbh_t, track_dbh_tp1, params_row, interval_years = 5) {
  if (!exists('transition_cost_tracks_bio_components', mode='function')) {
    source(file.path('dp_global','R','dp_global_bio.R'))
  }
  comps <- transition_cost_tracks_bio_components(track_dbh_t = track_dbh_t, track_dbh_tp1 = track_dbh_tp1, interval_years = interval_years, mu_const = params_row$Bio_Mu_Growth, mu_gamma = params_row$Bio_Gamma_Growth, sigma0 = params_row$Bio_Sigma0_Growth, sigma1 = params_row$Bio_Sigma1_Growth, max_shrink = params_row$Bio_Max_Shrink, k_shrink = params_row$Bio_K_Shrink, max_growth = params_row$Bio_Max_Growth, max_growth_soft = params_row$Bio_Max_Growth_Soft, k_growth = params_row$Bio_K_Growth, h0 = params_row$Bio_H0_Mortality, beta = params_row$Bio_Beta_Mortality, recruit_meanlog = params_row$Bio_Recruit_Meanlog, recruit_sdlog = params_row$Bio_Recruit_Sdlog, recruit_max_dbh = params_row$Bio_Recruit_MaxDBH_unit, recruit_lambda = params_row$Bio_Recruitment_lambda)
  return(comps)
}

# Convert two costs into a recruit probability via softmin
p_recruit_from_costs <- function(cost_recruit, cost_link_prev) {
  costs <- c(recruit = cost_recruit, link_prev = cost_link_prev)
  # handle non-finite
  if (!is.finite(costs['recruit'])) return(0)
  m <- max(-costs)
  probs <- exp(-costs - m)
  probs <- probs / sum(probs)
  return(as.numeric(probs['recruit']))
}

# Helper: write tidy CSV report
save_report <- function(dt, fname) {
  dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
  fwrite(dt, fname)
  message('[save_report] wrote ', fname)
  invisible(TRUE)
}
