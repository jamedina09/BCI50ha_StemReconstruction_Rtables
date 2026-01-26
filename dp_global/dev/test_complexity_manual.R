rm(list = ls())

library(data.table)
library(here)
library(collapse)

# Config
INPUT_FILE <- here::here("data_simulation/data/simulated_for_benchmark.csv")
anchor_start <- 7L
slack_tracks <- 1L

xs_all <- as.data.table(fread(INPUT_FILE))
xs_all[, species := Species]

# Choose a source tag (use first tag or a known one)

base_tag <- xs_all$Tag # [1:10]
base_dt <- xs_all[Tag %in% base_tag, .(CensusID, Tag, TrueStemID, DBH, species)]

xs <- base_dt
setDT(xs)

# Prepare TrueStemID like main script
# xs[, TrueStemID := as.character(TrueStemID)]
xs[CensusID < anchor_start, TrueStemID := NA_integer_]

# Use collapse to compute Tag x Census counts
key <- paste(xs$Tag, xs$CensusID, sep = ":::")
key_non_na <- key[!is.na(xs$DBH)]
cnts <- collapse::fcount(key_non_na)

# Parse output
if (is.data.frame(cnts)) {
  keys <- as.character(cnts$x)
  vals <- as.integer(cnts$N)
} else {
  keys <- names(cnts)
  vals <- if (is.list(cnts)) as.integer(unlist(cnts)) else as.integer(cnts)
}

spl <- strsplit(keys, ":::")
len <- lengths(spl)
if (!all(len == 2L)) stop("unexpected keys")
split_k <- do.call(rbind, spl)
dt_counts <- data.table(Tag = split_k[, 1], CensusID = as.integer(split_k[, 2]), n = vals)

# Keep counts up to anchor and pivot
dt_counts_anchor <- dt_counts[CensusID <= anchor_start]
wide <- data.table::dcast(dt_counts_anchor, Tag ~ CensusID, value.var = "n", fill = 0)

message("Wide counts:")
print(wide)

# Ensure census columns
census_cols <- as.character(seq_len(anchor_start))
missing_cols <- setdiff(census_cols, setdiff(names(wide), "Tag"))
for (col in missing_cols) wide[, (col) := 0L]
setcolorder(wide, c("Tag", census_cols))

# Compute simple metrics (K roughly) and print
mat <- as.matrix(wide[, ..census_cols])
max_obs_any_census <- do.call(pmax, as.data.frame(mat))

# births (handle various shapes robustly)
if (ncol(mat) >= 2L) {
  if (ncol(mat) == 2L) {
    # two-census case: single diff per row
    diffs_vec <- mat[, 2] - mat[, 1]
    births_needed <- as.integer(pmax(0L, diffs_vec))
  } else if (nrow(mat) > 1L) {
    diffs <- mat[, -1, drop = FALSE] - mat[, -ncol(mat), drop = FALSE]
    # avoid pmax() turning matrix into vector; do elementwise floor at zero
    diffs[diffs < 0] <- 0L
    births_needed <- rowSums(diffs)
  } else {
    # single-row multi-census case
    diffr <- pmax(0L, mat[1, -1, drop = TRUE] - mat[1, -ncol(mat), drop = TRUE])
    births_needed <- as.integer(sum(diffr))
  }
} else {
  births_needed <- integer(nrow(mat))
}

K_from_counts <- as.integer(mat[, 1] + births_needed)

K_base <- pmax(0L, max_obs_any_census, K_from_counts)
K <- K_base
if (slack_tracks > 0L) K[K_base == max_obs_any_census] <- K[K_base == max_obs_any_census] + slack_tracks

# n_states via lgamma
K_mat <- matrix(K, nrow = nrow(mat), ncol = ncol(mat))
n_obs_mat <- mat
n_states_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))

n_states_mat[n_obs_mat == 0] <- 1
n_states_mat[n_obs_mat > K_mat] <- 0
ok_mask <- (n_obs_mat != 0) & (n_obs_mat <= K_mat)
if (any(ok_mask)) {
  n_states_mat[ok_mask] <- exp(lgamma(K_mat[ok_mask] + 1) - lgamma(K_mat[ok_mask] - n_obs_mat[ok_mask] + 1))
}

total_states <- rowSums(n_states_mat)
max_states_per_census <- apply(n_states_mat, 1, max)
transition_computations <- rowSums(n_states_mat[, -ncol(n_states_mat), drop = FALSE] * n_states_mat[, -1, drop = FALSE])

res <- data.table(Tag = wide$Tag, K = as.integer(K), MaxObs = as.integer(max_obs_any_census), TotalStates = total_states, TransitionComputations = transition_computations)

## sort res by transition computations descending
setorder(res, -TransitionComputations)

# DP_MAX_STATES parameter: set the threshold (on TransitionComputations) above which tags will be
# skipped. You can override it via environment variable e.g. Sys.setenv(DP_MAX_STATES = "500000").
DP_MAX_STATES <- as.numeric(Sys.getenv("DP_MAX_STATES", "1e6"))

# Add per-tag columns:
# - dp_max_states: the configured cap applied to each Tag (equal to DP_MAX_STATES)
# - WillRun: logical indicating whether this Tag would be run (TRUE when TransitionComputations <= DP_MAX_STATES)
res[, dp_max_states := as.numeric(DP_MAX_STATES)]
res[, WillRun := TransitionComputations <= DP_MAX_STATES]

message(sprintf("DP_MAX_STATES = %g; %d tags would run; %d tags would be skipped.",
                DP_MAX_STATES, sum(res$WillRun, na.rm = TRUE), sum(!res$WillRun, na.rm = TRUE)))

# show top skipped tags (if any)
if (any(!res$WillRun)) {
  message("Top skipped tags by TransitionComputations:")
  print(head(res[WillRun == FALSE][order(-TransitionComputations)], n = 10))
}


# predict_time(N, unit = c("sec","min","hour")) - Predict runtime from TransitionComputations.
# - N: numeric scalar or vector of TransitionComputations
# - unit: one of "sec" (seconds), "min" (minutes), or "hour" (hours). Default is "sec".
predict_time <- function(N, unit = c("sec", "min", "hour")) {
  unit <- match.arg(unit)
  logN <- log10(N)
  logT <- -1.45817 +
    (-0.07313 * logN) +
    (0.14164 * logN^2)
  secs <- 10^logT
  out <- switch(unit,
    sec = secs,
    min = secs / 60,
    hour = secs / 3600
  )
  out
}

res[, PredictedTime_hour := predict_time(TransitionComputations, "hour")]

round(unique(res[TotalStates <= 1000]$PredictedTime_hour), 3)

# --- Helper: suggest DP max_states matching a target runtime (hours) ---
# Usage examples:
#  - interactively: suggest_dp_max_states(1)  # returns list with per-tag suggestions
#  - via env var: set DP_SUGGEST_HOURS=1 and run the script to add columns to `res` and print summary

# Add MaxStatesPerCensus and observed-intervals to `res` for use in suggestions
res[, MaxStatesPerCensus := as.integer(max_states_per_census)]
# Number of observed census per tag (>=1)
n_obs_per_tag <- rowSums(mat > 0)
n_intervals_per_tag <- pmax(1L, as.integer(n_obs_per_tag - 1L))
res[, N_Intervals := n_intervals_per_tag]

# Invert the predict_time() curve numerically to map target hours -> TransitionComputations threshold
invert_predict_time <- function(target_hours) {
  if (!is.numeric(target_hours) || length(target_hours) != 1 || is.na(target_hours) || target_hours <= 0) stop("target_hours must be a single positive number")
  f_log <- function(logN) predict_time(10^logN, "hour") - target_hours
  lo <- -6; hi <- 16 # search bounds (10^-6 .. 10^16)
  if (f_log(hi) <= 0) return(10^hi)
  if (f_log(lo) > 0) return(10^lo)
  10^uniroot(f_log, lower = lo, upper = hi)$root
}

suggest_dp_max_states <- function(target_hours) {
  N_thresh <- invert_predict_time(target_hours)
  # rough per-tag dp_max_states estimate using upper bound TransitionComputations <= sum S_t * S_{t+1} <= (n_intervals) * S^2
  # => S ≈ sqrt(N_thresh / n_intervals)
  suggested_S <- ceiling(sqrt(N_thresh / pmax(1, n_intervals_per_tag)))
  list(
    target_hours = target_hours,
    TransitionComputationsThreshold = N_thresh,
    suggested_dp_max_states_per_tag = suggested_S,
    suggested_dp_max_states_median = as.integer(median(suggested_S, na.rm = TRUE)),
    suggested_dp_max_states_90 = as.integer(quantile(suggested_S, 0.9, na.rm = TRUE))
  )
}

# If user sets DP_SUGGEST_HOURS env var, compute suggestions and add to `res`
DP_SUGGEST_HOURS <- as.numeric(Sys.getenv("DP_SUGGEST_HOURS", NA))
if (!is.na(DP_SUGGEST_HOURS)) {
  s <- suggest_dp_max_states(DP_SUGGEST_HOURS)
  res[, SuggestedTransitionComputations := as.numeric(s$TransitionComputationsThreshold)]
  res[, SuggestedDPMaxStates := as.integer(s$suggested_dp_max_states_per_tag)]
  message(sprintf("Suggested TransitionComputations threshold for %g hours: %g", DP_SUGGEST_HOURS, s$TransitionComputationsThreshold))
  message(sprintf("Suggested dp_max_states (per-tag) - median: %d; 90th pct: %d", s$suggested_dp_max_states_median, s$suggested_dp_max_states_90))
}

# Ensure output directory exists before saving plot
with(
  res[is.finite(PredictedTime_hour), ],
  plot(TransitionComputations, PredictedTime_hour,
    log = "xy",
    xlab = "Transition Computations", ylab = "Predicted Time (hours)",
    main = "Predicted Runtime vs Transition Computations",
    type = "l"
  )
)

# --- New plot: Predicted time vs DP max_states ---
# Summary across tags: for a grid of dp_max_states (S), estimate TransitionComputations <= n_intervals * S^2
# and convert to predicted time using predict_time(). We plot median and 90th percentile across tags.
S_guess_median <- if ("SuggestedDPMaxStates" %in% names(res)) median(res$SuggestedDPMaxStates, na.rm = TRUE) else NA_real_
S_max <- max(1000L, as.integer(na.omit(c(max(res$MaxStatesPerCensus, na.rm = TRUE), DP_MAX_STATES, ifelse(is.finite(S_guess_median), S_guess_median * 10L, 0L)))))
S_max <- min(S_max, 1e6)
S_grid <- unique(ceiling(10^seq(log10(1), log10(S_max), length.out = 200)))

# Ensure n_intervals_per_tag is available (defined earlier)
if (!exists("n_intervals_per_tag")) n_intervals_per_tag <- pmax(1L, as.integer(rowSums(mat > 0) - 1L))

# Build a matrix of predicted times: columns = S_grid, rows = tags
time_mat <- sapply(S_grid, function(S) {
  N_per_tag <- pmax(1, n_intervals_per_tag) * (S^2)
  predict_time(N_per_tag, "hour")
})

median_time <- apply(time_mat, 2, median, na.rm = TRUE)
q90_time <- apply(time_mat, 2, function(x) quantile(x, 0.9, na.rm = TRUE))

plot(S_grid, median_time,
  log = "xy",
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "DP max states per census (S)",
  ylab = "Predicted Time (hours)",
  main = "Predicted Runtime vs DP max states"
)
lines(S_grid, q90_time, col = "red", lty = 2, lwd = 2)
# horizontal reference lines: 1 min, 1 hour, 3 hours
abline(h = c(1/60, 1, 3), col = c("grey50", "black", "grey50"), lty = 3)
legend("bottomleft", legend = c("Median across tags", "90th percentile across tags", "refs: 1 min, 1 hr, 3 hr"), col = c("blue", "red", "black"), lty = c(1, 2, 3), bty = "n")

# Add per-tag suggested points (if present)
if ("SuggestedDPMaxStates" %in% names(res)) {
  pts_x <- res$SuggestedDPMaxStates
  pts_y <- predict_time(pmax(1, n_intervals_per_tag) * (pts_x^2), "hour")
  points(pts_x, pts_y, pch = 20, col = "darkgreen")
}

