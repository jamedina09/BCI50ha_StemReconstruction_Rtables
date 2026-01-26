rm(list = ls())

library(data.table)
library(here)
library(collapse)

# Config
INPUT_FILE <- here::here("./data_simulation/data/simulated_data_1.csv")
anchor_start <- 7L
slack_tracks <- 1L

xs_all <- as.data.table(read.csv(INPUT_FILE))

# Choose a source tag (use first tag or a known one)

base_tag <- xs_all$Tag#[1:10]
base_dt <- xs_all[Tag %in% base_tag, .(CensusID, Tag, TrueStemID, DBH, Species)]

xs <- base_dt
setDT(xs)

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
dt_counts <- data.table(Tag = split_k[,1], CensusID = as.integer(split_k[,2]), n = vals)

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
} else births_needed <- integer(nrow(mat))

K_from_counts <- as.integer(mat[,1] + births_needed)

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

