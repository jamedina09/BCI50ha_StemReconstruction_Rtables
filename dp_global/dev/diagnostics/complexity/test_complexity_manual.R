rm(list = ls())

# Minimal complexity diagnostic for DP decisions

library(data.table)
library(here)
library(collapse)

# Configuration
INPUT_FILE <- here::here("data_simulation/data/simulated_for_benchmark.csv")
ANCHOR_START <- 7L
SLACK_TRACKS <- 1L
# Default DP cap (per-census max states) — matches main defaults
DP_MAX_STATES <- as.integer(Sys.getenv("DP_MAX_STATES", "40000"))
# Optionally specify target runtime (hours) via DP_SUGGEST_HOURS; leave unset to skip suggestions
DP_SUGGEST_HOURS <- as.numeric(Sys.getenv("DP_SUGGEST_HOURS", NA))

# Output directory for diagnostics (inside the 'dev' folder)
outdir <- here::here("dp_global", "dev", "diagnostics", "complexity")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# Read data and compute per-Tag observed counts up to anchor
xs_all <- as.data.table(fread(INPUT_FILE))
xs_all[, species := Species]
base_dt <- xs_all[, .(CensusID, Tag, TrueStemID, DBH)][, .SD, by = Tag]
base_dt[CensusID < ANCHOR_START, TrueStemID := NA_integer_]

key <- paste(base_dt$Tag, base_dt$CensusID, sep = ":::")
key_non_na <- key[!is.na(base_dt$DBH)]
cnts <- collapse::fcount(key_non_na)

if (is.data.frame(cnts)) {
  keys <- as.character(cnts$x); vals <- as.integer(cnts$N)
} else {
  keys <- names(cnts); vals <- if (is.list(cnts)) as.integer(unlist(cnts)) else as.integer(cnts)
}

spl <- strsplit(keys, ":::")
split_k <- do.call(rbind, spl)
dt_counts <- data.table(Tag = split_k[, 1], CensusID = as.integer(split_k[, 2]), n = vals)

dt_counts_anchor <- dt_counts[CensusID <= ANCHOR_START]
wide <- data.table::dcast(dt_counts_anchor, Tag ~ CensusID, value.var = "n", fill = 0)

# Matrix of obs counts per census
census_cols <- as.character(seq_len(ANCHOR_START))
for (col in setdiff(census_cols, setdiff(names(wide), "Tag"))) wide[, (col) := 0L]
setcolorder(wide, c("Tag", census_cols))
mat <- as.matrix(wide[, ..census_cols])

# Basic per-tag metrics
max_obs_any_census <- do.call(pmax, as.data.frame(mat))
# births
if (ncol(mat) >= 2L) {
  if (ncol(mat) == 2L) {
    births_needed <- as.integer(pmax(0L, mat[, 2] - mat[, 1]))
  } else if (nrow(mat) > 1L) {
    diffs <- mat[, -1, drop = FALSE] - mat[, -ncol(mat), drop = FALSE]
    diffs[diffs < 0] <- 0L
    births_needed <- rowSums(diffs)
  } else {
    births_needed <- as.integer(sum(pmax(0L, mat[1, -1, drop = TRUE] - mat[1, -ncol(mat), drop = TRUE])))
  }
} else {
  births_needed <- integer(nrow(mat))
}
K_from_counts <- as.integer(mat[, 1] + births_needed)
K_base <- pmax(0L, max_obs_any_census, K_from_counts)
K <- K_base
if (SLACK_TRACKS > 0L) K[K_base == max_obs_any_census] <- K[K_base == max_obs_any_census] + SLACK_TRACKS

# Estimate number of states per census using permutation count P(K, n_obs)
K_mat <- matrix(K, nrow = nrow(mat), ncol = ncol(mat))
n_obs_mat <- mat
n_states_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))

n_states_mat[n_obs_mat == 0] <- 1
n_states_mat[n_obs_mat > K_mat] <- 0
ok_mask <- (n_obs_mat != 0) & (n_obs_mat <= K_mat)
if (any(ok_mask)) {
  n_states_mat[ok_mask] <- exp(lgamma(K_mat[ok_mask] + 1) - lgamma(K_mat[ok_mask] - n_obs_mat[ok_mask] + 1))
}

TotalStates <- rowSums(n_states_mat)
MaxStatesPerCensus <- apply(n_states_mat, 1, max)
TransitionComputations <- rowSums(n_states_mat[, -ncol(n_states_mat), drop = FALSE] * n_states_mat[, -1, drop = FALSE])

res <- data.table(Tag = wide$Tag, K = as.integer(K), MaxObs = as.integer(max_obs_any_census), TotalStates = TotalStates, MaxStatesPerCensus = as.integer(MaxStatesPerCensus), TransitionComputations = TransitionComputations)
setorder(res, -TransitionComputations)

# DP gating checks
res[, dp_max_states := as.integer(DP_MAX_STATES)]
res[, WillRunByDP := MaxStatesPerCensus <= dp_max_states]
res[, WillRunByTransition := TransitionComputations <= dp_max_states]

cat(sprintf("DP_MAX_STATES = %d (per-census max states cap)\n", DP_MAX_STATES))
cat(sprintf("WillRunByDP: %d run / %d skip\n", sum(res$WillRunByDP), sum(!res$WillRunByDP)))
cat(sprintf("WillRunByTransition: %d run / %d skip\n\n", sum(res$WillRunByTransition), sum(!res$WillRunByTransition)))

# Predict runtime helper
predict_time <- function(N, unit = c("sec", "min", "hour")) {
  unit <- match.arg(unit)
  logN <- log10(N)
  logT <- -1.45817 + (-0.07313 * logN) + (0.14164 * logN^2)
  secs <- 10^logT
  switch(unit, sec = secs, min = secs / 60, hour = secs / 3600)
}

res[, PredictedTime_hour := predict_time(TransitionComputations, "hour")]

# Suggest dp_max_states for a target runtime in hours
DP_SUGGEST_HOURS <- as.numeric(Sys.getenv("DP_SUGGEST_HOURS", NA))
if (!is.na(DP_SUGGEST_HOURS)) {
  invert_predict_time <- function(target_hours) {
    a <- -1.45817; b <- -0.07313; c <- 0.14164
    logN_min <- -b / (2 * c)
    f_log <- function(logN) predict_time(10^logN, "hour") - target_hours
    lo <- max(logN_min + 1e-6, 0); hi <- 16
    if (f_log(hi) <= 0) return(10^hi)
    if (f_log(lo) > 0) return(10^lo)
    10^uniroot(f_log, lower = lo, upper = hi)$root
  }
  N_thresh <- invert_predict_time(DP_SUGGEST_HOURS)
  n_obs_per_tag <- rowSums(mat > 0)
  n_intervals_per_tag <- pmax(1L, as.integer(n_obs_per_tag - 1L))
  res[, SuggestedDPMaxStates := as.integer(ceiling(sqrt(N_thresh / pmax(1, n_intervals_per_tag))))]
  res[, SuggestedPredictedTime_hour := predict_time(pmax(1, n_intervals_per_tag) * (SuggestedDPMaxStates^2), "hour")]

  medS <- as.integer(median(res$SuggestedDPMaxStates, na.rm = TRUE))
  p90S <- as.integer(quantile(res$SuggestedDPMaxStates, 0.9, na.rm = TRUE))
  cat(sprintf("\nTarget runtime: %g hours (%g minutes)\n", DP_SUGGEST_HOURS, DP_SUGGEST_HOURS * 60))
  cat(sprintf("TransitionComputations threshold ~ %g\n", N_thresh))
  cat(sprintf("Suggested dp_max_states (median=%d, 90th=%d)\n", medS, p90S))
  cat(sprintf("If DP_MAX_STATES = %d -> %d/%d tags run (%.1f%%)\n", medS, sum(res$MaxStatesPerCensus <= medS), nrow(res), 100 * sum(res$MaxStatesPerCensus <= medS) / nrow(res)))

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  data.table::fwrite(res[order(-SuggestedDPMaxStates), .(Tag, K, MaxObs, TotalStates, TransitionComputations, MaxStatesPerCensus, SuggestedDPMaxStates, SuggestedPredictedTime_hour)], file = file.path(outdir, sprintf("suggested_dp_max_states_%ghours.csv", DP_SUGGEST_HOURS)))
  cat(sprintf("Wrote per-tag suggestions CSV to: %s\n", file.path(outdir, sprintf("suggested_dp_max_states_%ghours.csv", DP_SUGGEST_HOURS))))
}

# Write res summary CSV
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
res_out <- res[, .(Tag, K, MaxObs, TotalStates, MaxStatesPerCensus, TransitionComputations, dp_max_states, WillRunByDP, WillRunByTransition, PredictedTime_hour)]
data.table::fwrite(res_out, file = file.path(outdir, "res_summary.csv"))
cat(sprintf("Wrote res summary to: %s\n", file.path(outdir, "res_summary.csv")))

# Plot Predicted time vs DP max states (median & 90th)
S_max <- max(1000L, as.integer(max(res$SuggestedDPMaxStates, na.rm = TRUE)), DP_MAX_STATES)
S_max <- min(S_max, 1e6)
S_grid <- unique(ceiling(10^seq(log10(1), log10(S_max), length.out = 200)))

n_intervals_per_tag <- pmax(1L, as.integer(rowSums(mat > 0) - 1L))
time_mat <- sapply(S_grid, function(S) {
  N_per_tag <- pmax(1, n_intervals_per_tag) * (S^2)
  predict_time(N_per_tag, "hour")
})
median_time <- apply(time_mat, 2, median, na.rm = TRUE)
q90_time <- apply(time_mat, 2, function(x) quantile(x, 0.9, na.rm = TRUE))

png(file.path(outdir, "PredictedTime_vs_DPMaxStates.png"), width = 900, height = 600)
plot(S_grid, median_time, log = "xy", type = "l", lwd = 2, col = "blue", xlab = "DP max states per census (S)", ylab = "Predicted Time (hours)", main = "Predicted Runtime vs DP max states")
lines(S_grid, q90_time, col = "red", lty = 2, lwd = 2)
abline(h = c(1/60, 1, 3), col = c("grey50", "black", "grey50"), lty = 3)
legend("bottomleft", legend = c("Median across tags", "90th percentile across tags", "refs: 1 min, 1 hr, 3 hr"), col = c("blue", "red", "black"), lty = c(1, 2, 3), bty = "n")
if (!is.na(DP_SUGGEST_HOURS)) points(res$SuggestedDPMaxStates, res$SuggestedPredictedTime_hour, pch = 20, col = "darkgreen")
dev.off()

cat(sprintf("Wrote plot to: %s\n", file.path(outdir, "PredictedTime_vs_DPMaxStates.png")))

cat("Completed test_complexity_manual minimal run.\n")