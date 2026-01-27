# Simulate datasets for DP benchmark
# Produces: data_simulation/data/simulated_for_benchmark.csv

library(data.table)
library(here)

# Config
OUT_PATH <- here::here("data_simulation", "data", "simulated_for_benchmark.csv")
ANCHOR_START <- 7L
BIO_VALS <- list(
  Bio_Mu_Growth = 0.01452756,
  Bio_Gamma_Growth = 0.07163091,
  Bio_Sigma0_Growth = 0.003578623,
  Bio_Sigma1_Growth = 0.006697133,
  Bio_H0_Mortality = 0.04685938,
  Bio_Beta_Mortality = -0.02755028,
  Bio_Recruit_Meanlog = 0.4339083,
  Bio_Recruit_Sdlog = 0.4847472,
  Bio_Recruit_MaxDBH_unit = 36.5001,
  Bio_Recruitment_lambda = 0.08044933,
  Bio_Max_Shrink = -1,
  Bio_K_Shrink = 0,
  Bio_Max_Growth = 7.5,
  Bio_Max_Growth_Soft = 0.8861551,
  Bio_K_Growth = 0
)

# Threshold for stopping (TransitionComputations)
TC_STOP <- 2e6
MAX_STEMS <- 7L

# Helper: compute transition computations for tags (copied logic)
compute_transition_table <- function(xs, anchor_start = ANCHOR_START) {
  xs <- copy(xs)
  setDT(xs)
  key <- paste(xs$Tag, xs$CensusID, sep = ":::")
  key_non_na <- key[!is.na(xs$DBH)]
  cnts <- collapse::fcount(key_non_na)

  if (is.data.frame(cnts)) {
    keys <- as.character(cnts$x)
    vals <- as.integer(cnts$N)
  } else {
    keys <- names(cnts)
    vals <- if (is.list(cnts)) as.integer(unlist(cnts)) else as.integer(cnts)
  }

  spl <- strsplit(keys, ":::")
  split_k <- do.call(rbind, spl)
  dt_counts <- data.table(Tag = split_k[,1], CensusID = as.integer(split_k[,2]), n = vals)
  dt_counts_anchor <- dt_counts[CensusID <= anchor_start]
  wide <- data.table::dcast(dt_counts_anchor, Tag ~ CensusID, value.var = "n", fill = 0)

  census_cols <- as.character(seq_len(anchor_start))
  missing_cols <- setdiff(census_cols, setdiff(names(wide), "Tag"))
  for (col in missing_cols) wide[, (col) := 0L]
  setcolorder(wide, c("Tag", census_cols))

  mat <- as.matrix(wide[, ..census_cols])
  K_mat <- matrix(as.integer(do.call(pmax, as.data.frame(mat))), nrow = nrow(mat), ncol = ncol(mat))
  n_obs_mat <- mat
  n_states_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  n_states_mat[n_obs_mat == 0] <- 1
  ok_mask <- (n_obs_mat != 0) & (n_obs_mat <= K_mat)
  if (any(ok_mask)) {
    n_states_mat[ok_mask] <- exp(lgamma(K_mat[ok_mask] + 1) - lgamma(K_mat[ok_mask] - n_obs_mat[ok_mask] + 1))
  }
  total_states <- rowSums(n_states_mat)
  transition_computations <- rowSums(n_states_mat[, -ncol(n_states_mat), drop = FALSE] * n_states_mat[, -1, drop = FALSE])
  data.table(Tag = wide$Tag, TransitionComputations = transition_computations)
}

# Build synthetic tags
rows <- list()
tag_id <- 1L
results <- list()
for (n_stems in seq_len(MAX_STEMS)) {
  tag <- paste0("nstem_", n_stems)
  # For simplicity, we observe all stems every census (max complexity)
  for (c in seq_len(ANCHOR_START)) {
    for (i in seq_len(n_stems)) {
      dbh <- round(10 + i * 0.5 + (c - 1) * 0.2, 2)
      rows[[length(rows) + 1]] <- data.table(
        Tag = tag,
        CensusID = as.integer(c),
        DBH = dbh,
        TrueStemID = if (c == ANCHOR_START) as.integer(i) else NA_integer_
      )
    }
  }
  # temporary table and compute TC
  tmp <- rbindlist(rows)
  # compute TC for all tags currently present
  tc_tab <- compute_transition_table(tmp, anchor_start = ANCHOR_START)
  # find TC for this tag
  tc_val <- tc_tab[Tag == tag, TransitionComputations]
  results[[length(results) + 1]] <- list(Tag = tag, NStems = n_stems, TransitionComputations = tc_val)
  # Stop when we have generated up to the requested max stems
  if (n_stems >= MAX_STEMS) break
}

# Final dataset: attach bio vals and write
Dt <- rbindlist(rows)
for (nm in names(BIO_VALS)) Dt[, (nm) := as.numeric(BIO_VALS[[nm]])]
# Ensure required columns are present for DP
if (!("Species" %in% names(Dt))) Dt[, Species := "synth"]
if (!("OriginalStemID" %in% names(Dt))) Dt[, OriginalStemID := NA_integer_]
if (!("ExactDate" %in% names(Dt))) Dt[, ExactDate := as.IDate(Sys.Date())]

fwrite(Dt, OUT_PATH)
cat("Wrote synthetic dataset to:", OUT_PATH, "\n")
cat("Summary of generated tags (Tag, NStems, TransitionComputations):\n")
print(rbindlist(lapply(results, as.data.table)))

invisible(NULL)
