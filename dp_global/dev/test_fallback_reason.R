#!/usr/bin/env Rscript
library(data.table)
source(file.path('dp_global','R','dp_global_main.R'))

# Construct the user's example
x <- data.table(
  Tag = rep('014993', 3),
  StemID = rep('10989', 3),
  TrueStemID = c(NA_integer_, 1L, NA_integer_),
  CensusID = 1:3,
  DBH = c(58.6, 45.2, NA_real_),
  ExactDate = as.IDate(c('2000-01-01','2005-01-01','2010-01-01'))
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
  anchor_start = 2L,
  max_tracks = 5L,
  max_states = 5000L,
  posterior_samples = 0L,
  verbose = TRUE
)

# Assert we fell back to igraph (some rows should have ReconstructionMethod either 'igraph' or 'provisional_igraph')
if (!any(grepl('igraph', res$ReconstructionMethod, ignore.case = TRUE))) stop('Expected igraph-based reconstruction method')
# Assert fallback reason present
if (!('DP_FallbackReason' %in% names(res))) stop('Missing DP_FallbackReason column')
if (all(is.na(res$DP_FallbackReason))) stop('DP_FallbackReason should be set when falling back')

cat('OK: fallback reason test passed; DP_FallbackReason values:', unique(na.omit(res$DP_FallbackReason)), '\n')
