library(data.table)

# Simple stub for transition cost so we don't require Rcpp compilation during testing
transition_cost_tracks_bio_batch_rcpp <- function(
  track_dbh_t,
  track_dbh_tp1,
  interval_years,
  ...,
  eps_tiebreak = 0,
  hard_penalty = 1e6
) {
    # track_dbh_tp1 may be a list of vectors or a matrix
    if (is.list(track_dbh_tp1)) {
        return(rep(0, length(track_dbh_tp1)))
    }
    mat <- as.matrix(track_dbh_tp1)
    return(rep(0, nrow(mat)))
}

# Source necessary internals
src_dir <- file.path("dp_global", "R")
source(file.path(src_dir, "dp_global_states.R"))
source(file.path(src_dir, "dp_global_utils.R"))
source(file.path(src_dir, "dp_global_diag.R"))
source(file.path(src_dir, "dp_global_dp.R"))

# Construct a small simulated dataset
# Two stems present in Census 1..4 and missing in 5..6; anchor at 7 has two anchor trees.
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

# Census 5 and 6: no observations for those two stems (they died)
# anchor census 7: two anchor observations with TrueStemID defined
rows[[length(rows) + 1]] <- list(
    Tag = 1L, species = "sp",
    CensusID = 7L, DBH = NA_real_, ExactDate = as.Date("2000-01-01") + 6 * 365, TrueStemID = NA_integer_
)

rows[[length(rows) + 1]] <- list(
    Tag = 1L, species = "sp",
    CensusID = 7L, DBH = NA_real_, ExactDate = as.Date("2000-01-01") + 6 * 365, TrueStemID = NA_integer_
)

# Build data.table and add required bio columns
dt <- rbindlist(lapply(rows, as.data.table))
## add rowid
dt[, RowID := .I]
dt[RowID == 7, TrueStemID := 1L]
dt[RowID == 8, TrueStemID := 2L]

# Append required bio columns with simple constants
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
print("Input data:")
print(dt)

# Run the DP with anchor_start = 7
out <- match_stems_dp_global_backward_marginals_batch(
    tree_data = dt,
    min_growth = -Inf,
    max_growth = Inf,
    anchor_start = 7L,
    max_tracks = 6L,
    slack_tracks = 1L,
    max_states = 50000L,
    temperature = 1.0,
    verbose = TRUE, 
    posterior_samples = 10, 
    posterior_sample_seed = 42L,
    posterior_samples_format = "csv", 
    posterior_samples_path = "./"

)

print("Output summary (first rows):")
print(out)

# Summarize how many rows got ReconstructionMethod == 'dp' or 'given' or 'igraph'
print("ReconstructionMethod counts:")
print(table(out$ReconstructionMethod, useNA = "ifany"))

# Show ReconstructedStemID for census 1..4
print("Census 1-4 rows:")
print(out[CensusID %in% 1:4, .(CensusID, DBH, TrueStemID, ReconstructedStemID, ReconstructionMethod)])

# Show DP diagnostic columns
print("DP diagnostics:")
print(unique(out[, .(DP_KUsed, DP_MaxStatesPerCensus, DP_MaxStatesCensusID)]))

fwrite(out, "./out.csv")
