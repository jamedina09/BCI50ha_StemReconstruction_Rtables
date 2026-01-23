library(data.table)

# Use the same stubbed transition cost as the test script (fast)
transition_cost_tracks_bio_batch_rcpp <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  ...,
  eps_tiebreak = 0,
  hard_penalty = 1e6
) {
    if (is.list(track_dbh_tp1)) {
        return(rep(0, length(track_dbh_tp1)))
    }
    mat <- as.matrix(track_dbh_tp1)
    return(rep(0, nrow(mat)))
}

# Source project internals
src_dir <- file.path("dp_global", "R")
source(file.path(src_dir, "dp_global_states.R"))
source(file.path(src_dir, "dp_global_utils.R"))
source(file.path(src_dir, "dp_global_diag.R"))
source(file.path(src_dir, "dp_global_matchers.R"))
source(file.path(src_dir, "dp_global_dp.R"))

# Load the example dataset from the other test (re-create in script to be self-contained)
rows <- list()
for (cc in 1:4) {
    rows[[length(rows) + 1]] <- list(
        Tag = 1L, species = "sp",
        CensusID = cc, DBH = 10 + (cc - 1) * 1, ExactDate = as.Date("2000-01-01") + (cc - 1) * 365, TrueStemID = NA_integer_
    )
    rows[[length(rows) + 1]] <- list(
        Tag = 1L, species = "sp",
        CensusID = cc, DBH = 20 + (cc - 1) * 1, ExactDate = as.Date("2000-01-01") + (cc - 1) * 365, TrueStemID = NA_integer_
    )
}
rows[[length(rows) + 1]] <- list(
    Tag = 1L, species = "sp",
    CensusID = 7L, DBH = NA_real_, ExactDate = as.Date("2000-01-01") + 6 * 365, TrueStemID = NA_integer_
)
rows[[length(rows) + 1]] <- list(
    Tag = 1L, species = "sp",
    CensusID = 7L, DBH = NA_real_, ExactDate = as.Date("2000-01-01") + 6 * 365, TrueStemID = NA_integer_
)

dt <- rbindlist(lapply(rows, as.data.table))
## add rowid
dt[, RowID := .I]
dt[RowID == 7, TrueStemID := 1L]
dt[RowID == 8, TrueStemID := 2L]

bio_cols <- list(
    Bio_Mu_Growth = 0.5,
    Bio_Gamma_Growth = 0,
    Bio_Sigma0_Growth = 0.1,
    Bio_Sigma1_Growth = 0,
    Bio_Max_Shrink = -5,
    Bio_K_Shrink = 0,
    Bio_Max_Growth = 50,
    Bio_Max_Growth_Soft = 50,
    Bio_K_Growth = 0,
    Bio_H0_Mortality = 0.01,
    Bio_Beta_Mortality = 0,
    Bio_Recruit_Meanlog = log(5),
    Bio_Recruit_Sdlog = 0.5,
    Bio_Recruit_MaxDBH_unit = 10,
    Bio_Recruitment_lambda = 0.1
)
for (nm in names(bio_cols)) dt[[nm]] <- bio_cols[[nm]]
setorder(dt, CensusID)

# Baseline parameters
baseline_min <- -5
baseline_max <- 15

baseline <- match_stems_dp_global_backward_marginals_batch(
    tree_data = dt,
    min_growth = baseline_min,
    max_growth = baseline_max,
    anchor_start = 7L,
    max_tracks = 6L,
    slack_tracks = 1L,
    max_states = 50000L,
    temperature = 1.0,
    prune_hard = TRUE,
    verbose = FALSE
)

# baseline vector to compare: ReconstructedStemID for CensusID 1..4, ordered by RowID
baseline_vec <- baseline[CensusID %in% 1:4, .(RowID, ReconstructedStemID)][order(RowID), ReconstructedStemID]

# Grid to search
min_vals <- seq(-10, 0, by = 1)
max_vals <- seq(5, 25, by = 1)

res <- list()
count <- 0
start_time <- proc.time()
for (mng in min_vals) {
    for (mxg in max_vals) {
        count <- count + 1
        out <- match_stems_dp_global_backward_marginals_batch(
            tree_data = dt,
            min_growth = mng,
            max_growth = mxg,
            anchor_start = 7L,
            max_tracks = 6L,
            slack_tracks = 1L,
            max_states = 50000L,
            temperature = 1.0,
            prune_hard = TRUE,
            verbose = FALSE
        )
        vec <- out[CensusID %in% 1:4, .(RowID, ReconstructedStemID)][order(RowID), ReconstructedStemID]
        same <- identical(vec, baseline_vec)
        any_igraph <- any(out$ReconstructionMethod == "igraph", na.rm = TRUE)
        prune_info <- attr(out, "DP_PruneInfo")
        total_examined <- if (!is.null(prune_info)) prune_info$total_examined else NA
        total_pruned <- if (!is.null(prune_info)) prune_info$total_pruned else NA
        res[[count]] <- list(min_growth = mng, max_growth = mxg, same = same, any_igraph = any_igraph, total_examined = total_examined, total_pruned = total_pruned)
    }
}
end_time <- proc.time()

res_dt <- rbindlist(res)
# summary: which parameter pairs keep same
same_dt <- res_dt[same == TRUE]
cat(sprintf("Baseline min=%g max=%g\n", baseline_min, baseline_max))
cat(sprintf("Grid search points: %d\n", nrow(res_dt)))
cat(sprintf("Number of pairs keeping same path as baseline: %d\n", nrow(same_dt)))

# For each min_growth, show contiguous range of max_growth that keep same
summary_by_min <- same_dt[, .(min_mx = min(max_growth), max_mx = max(max_growth), n = .N), by = min_growth][order(min_growth)]
print(summary_by_min)

# Write full results and the filtered pairs
fwrite(res_dt, "./growth_bounds_sweep_results.csv")
if (nrow(same_dt) > 0) fwrite(same_dt, "./growth_bounds_pairs_same.csv")

cat(sprintf("Elapsed time: %.2fs\n", (end_time - start_time)[3]))

print("Done. Results written to growth_bounds_sweep_results.csv and growth_bounds_pairs_same.csv (if any).")
