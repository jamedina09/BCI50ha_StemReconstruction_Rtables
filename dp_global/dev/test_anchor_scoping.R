#!/usr/bin/env Rscript
library(data.table)
source(file.path("dp_global","R","dp_global_main.R"))

# Construct test dataset: Tag 999 with anchor at 7, post-anchor DBH present and NA
x <- data.table(
  Species = rep("sp1", 3),
  Tag = rep(999L, 3),
  OriginalStemID = rep(1L, 3),
  TrueStemID = c(1L, 1L, NA_integer_),
  CensusID = c(7L, 8L, 9L),
  DBH = c(1.1, 2.1, NA_real_),
  ExactDate = as.IDate(c("2010-01-18", "2015-02-10", "2022-02-15"))
)

x[, species := "all"]

# Define minimal fixed constants used by the main runner
MAX_GROWTH_FIXED <- 7.5
MAX_SHRINK_FIXED <- -0.5
RECRUIT_MAX_FIXED <- (MAX_GROWTH_FIXED * 5) - 0.9999

# For this small test, construct minimal biological parameters (avoid needing many observations)
bio_pars <- list()
bio_pars[["all"]] <- list(
    growth = list(mu = 1, gamma = 0, sigma0 = 0.1, sigma1 = 0, guardrails = list(hard = list(value = MAX_GROWTH_FIXED), soft = list(value = MAX_GROWTH_FIXED * 0.5)), penalties = list(soft = list(k = 0.1))),
    shrinkage = list(max_shrink = MAX_SHRINK_FIXED, k_shrink = 0.1),
    mortality = list(h0 = 0.01, beta = 0.01),
    recruitment = list(meanlog = -0.25, sdlog = 0.8, recruit_max_dbh = RECRUIT_MAX_FIXED, lambda = 0.1),
    measurement_error = list(sd1_a = 0.0062, sd1_b = 0.0904, sd2 = 4.64, p_big = 0.05)
)
# Manually attach minimal Bio_* columns used by DP (avoid sourcing runner helpers)
x[, Bio_Mu_Growth := bio_pars[["all"]]$growth$mu]
x[, Bio_Gamma_Growth := bio_pars[["all"]]$growth$gamma]
x[, Bio_Sigma0_Growth := bio_pars[["all"]]$growth$sigma0]
x[, Bio_Sigma1_Growth := bio_pars[["all"]]$growth$sigma1]
x[, Bio_H0_Mortality := bio_pars[["all"]]$mortality$h0]
x[, Bio_Beta_Mortality := bio_pars[["all"]]$mortality$beta]
x[, Bio_Recruit_Meanlog := bio_pars[["all"]]$recruitment$meanlog]
x[, Bio_Recruit_Sdlog := bio_pars[["all"]]$recruitment$sdlog]
x[, Bio_Recruit_MaxDBH_unit := bio_pars[["all"]]$recruitment$recruit_max_dbh]
x[, Bio_Recruitment_lambda := bio_pars[["all"]]$recruitment$lambda]
x[, Bio_Max_Shrink := bio_pars[["all"]]$shrinkage$max_shrink]
x[, Bio_K_Shrink := bio_pars[["all"]]$shrinkage$k_shrink]
x[, Bio_Max_Growth := MAX_GROWTH_FIXED]
x[, Bio_Max_Growth_Soft := MAX_GROWTH_FIXED * 0.5]
x[, Bio_K_Growth := bio_pars[["all"]]$growth$penalties$soft$k]


res <- match_stems_dp_global_backward_marginals_batch(
  x,
  min_growth = -0.5,
  max_growth = 7.5,
  anchor_start = 7L,
  max_tracks = 10L,
  max_states = 100L,
  slack_tracks = 1L,
  posterior_samples = 0L,
  verbose = TRUE
)

# Assertions
stopifnot(nrow(res) == nrow(x))
# rows with CensusID <= anchor should have ReconstructedStemID possibly non-NA
pre_anchor_idx <- which(res$CensusID <= 7)
post_anchor_idx <- which(res$CensusID > 7)
stopifnot(length(pre_anchor_idx) > 0)
# At least one pre-anchor row should have a reconstructed stem id (not NA)
stopifnot(any(!is.na(res$ReconstructedStemID[pre_anchor_idx])))
# Post-anchor rows: rows with DBH and TrueStemID that were used by DP should be marked 'given'
post_idx_dbh_true <- post_anchor_idx[!is.na(res$TrueStemID[post_anchor_idx]) & !is.na(res$DBH[post_anchor_idx])]
used_ids <- unique(res$ReconstructedStemID[!is.na(res$ReconstructedStemID)])
if (length(post_idx_dbh_true) > 0 && length(used_ids) > 0) {
  for (i in post_idx_dbh_true) {
    if (res$TrueStemID[i] %in% used_ids) {
      stopifnot(res$ReconstructedStemID[i] == res$TrueStemID[i])
      stopifnot(!is.na(res$ReconstructionMethod[i]) && res$ReconstructionMethod[i] == "given")
    } else {
      stopifnot(is.na(res$ReconstructionMethod[i]) || res$ReconstructionMethod[i] == "none_after_anchor")
    }
  }
}

# Other post-anchor rows (no DBH or no TrueStemID) should be labelled or left as none_after_anchor
other_post_idx <- setdiff(post_anchor_idx, post_idx_dbh_true)
if (length(other_post_idx) > 0) {
  for (i in other_post_idx) {
    stopifnot(is.na(res$ReconstructionMethod[i]) || res$ReconstructionMethod[i] == "none_after_anchor")
  }
}

cat("OK: anchor scoping test updated and passed\n")
