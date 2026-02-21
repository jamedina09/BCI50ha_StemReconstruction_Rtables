library(data.table)
source(file.path('dp_global','R','dp_global_main.R'))

# Construct the user's example
x <- data.table(
  Tag = rep('014993', 5),
  StemID = rep('10989', 5),
  TrueStemID = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 10989L),
  CensusID = 1:5,
  DBH = c(46, 48, 50, 52, 54),
  ExactDate = as.IDate(c('2000-01-01','2005-01-01','2010-01-01', '2015-01-01', '2020-01-01'))
)

x[, species := 'sp']
# Minimal bio params
x[, Bio_Mu_Growth := 1]
x[, Bio_Gamma_Growth := 0]
x[, Bio_Sigma0_Growth := 0.1]
x[, Bio_Sigma1_Growth := 0]
x[, Bio_H0_Mortality := 0.01]
x[, Bio_Beta_Mortality := 0.01]
x[, Bio_Recruit_Meanlog := -0.25]
x[, Bio_Recruit_Sdlog := 0.8]
x[, Bio_Recruit_MaxDBH_unit := 10]
x[, Bio_Recruitment_lambda := 0.1]
x[, Bio_Max_Shrink := -0.5]
x[, Bio_K_Shrink := 0.1]
x[, Bio_Max_Growth := 7.5]

res <- match_stems_dp_global_backward_marginals_batch(
  x,
  min_growth = -0.5,
  max_growth = 7.5,
  anchor_start = 5L,
  max_tracks = 5L,
  max_states = 5000L,
  posterior_samples = 0L,
  fallback_growth_forms = "tree",
  verbose = TRUE
)

# --- growth form fallback --------------------------------------------------
x2 <- copy(x)
x2[, growth_form := 'fig']

res2 <- match_stems_dp_global_backward_marginals_batch(
  x2,
  min_growth = -0.5,
  max_growth = 7.5,
  anchor_start = 5L,
  max_tracks = 5L,
  max_states = 5000L,
  posterior_samples = 0L,
  fallback_growth_forms = "fig",
  verbose = TRUE
)

if (!any(grepl('igraph', res2$ReconstructionMethod, ignore.case = TRUE))) stop('Expected igraph-based reconstruction method for growth form')
if (!('DP_FallbackReason' %in% names(res2))) stop('Missing DP_FallbackReason column in growth form test')
if (all(is.na(res2$DP_FallbackReason))) stop('DP_FallbackReason should be set when growth_form triggers fallback')
if (!all(res2$DP_FallbackReason %in% 'growth_form_forced')) stop('Unexpected fallback reasons in growth form test')
cat('OK: growth_form fallback test passed; DP_FallbackReason =', unique(na.omit(res2$DP_FallbackReason)), '\n')
