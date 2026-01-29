#!/usr/bin/env Rscript
library(data.table)
source(file.path("dp_global","R","dp_global_main.R"))

# Integration test: ensure post-anchor rows with DBH + TrueStemID that were
# used by DP are annotated ReconstructedStemID == TrueStemID and
# ReconstructionMethod == "given".

# Construct test dataset for Tag 1001
x <- data.table(
  Tag = rep(1001L, 5),
  Species = rep("sp", 5),
  CensusID = 5:9,
  DBH = c(NA, 10, 12, 14, 16),
  ExactDate = as.IDate(c("2000-01-01","2001-01-01","2002-01-01","2003-01-01","2004-01-01")),
  TrueStemID = c(NA, 1L, 1L, 2L, 2L)
)

# Minimal bio columns to satisfy DP
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
x[, Bio_Max_Growth_Soft := 3.75]
x[, Bio_K_Growth := 0.1]

# Run DP with anchor_start = 7 (so CensusID > 7 are post-anchor)
res <- match_stems_dp_global_backward_marginals_batch(
  x,
  min_growth = -0.5,
  max_growth = 7.5,
  anchor_start = 7L,
  max_tracks = 5L,
  max_states = 100L,
  slack_tracks = 1L,
  posterior_samples = 0L,
  verbose = FALSE
)

stopifnot(nrow(res) == nrow(x))

# Identify pre/post anchor indices
pre_anchor_idx <- which(res$CensusID <= 7)
post_anchor_idx <- which(res$CensusID > 7)
stopifnot(length(pre_anchor_idx) > 0)
stopifnot(length(post_anchor_idx) > 0)

# Ensure at least one pre-anchor row has a reconstructed stem id
stopifnot(any(!is.na(res$ReconstructedStemID[pre_anchor_idx])))

# Identify post-anchor rows with DBH + TrueStemID
post_dbh_true_idx <- post_anchor_idx[!is.na(res$DBH[post_anchor_idx]) & !is.na(res$TrueStemID[post_anchor_idx])]

# Gather the set of ReconstructedStemIDs used in pre-anchor
used_ids <- unique(res$ReconstructedStemID[!is.na(res$ReconstructedStemID) & res$CensusID <= 7])

# If there are such rows and the DP used some ids, assert the 'given' propagation
if (length(post_dbh_true_idx) > 0 && length(used_ids) > 0) {
  for (i in post_dbh_true_idx) {
    if (res$TrueStemID[i] %in% used_ids) {
      stopifnot(!is.na(res$ReconstructionMethod[i]) && res$ReconstructionMethod[i] == "given")
      stopifnot(res$ReconstructedStemID[i] == res$TrueStemID[i])
    } else {
      stopifnot(is.na(res$ReconstructionMethod[i]) || res$ReconstructionMethod[i] == "none_after_anchor")
    }
  }
}

# Also check that other post-anchor rows are either NA or marked none_after_anchor
other_post_idx <- setdiff(post_anchor_idx, post_dbh_true_idx)
if (length(other_post_idx) > 0) {
  for (i in other_post_idx) {
    stopifnot(is.na(res$ReconstructionMethod[i]) || res$ReconstructionMethod[i] == "none_after_anchor")
  }
}

cat("OK: integration test 'post_anchor_given' passed\n")
