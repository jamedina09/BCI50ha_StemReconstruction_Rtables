rm(list = ls())

# Example: demonstrate DP behavior when census 5 has no observations
# Places an anchor at census 7 and shows reconstructed IDs and phases

library(data.table)

cur <- normalizePath(getwd())
dp_r_dir <- file.path(cur, "dp_global", "R")

if (dir.exists(dp_r_dir)) {
    # Only source core DP function files to avoid running examples/tests with side-effects
    whitelist <- c(
        "dp_global_bio.R",
        "dp_global_matchers.R",
        "dp_global_states.R",
        "dp_global_utils.R",
        "dp_global_dp.R",
        "dp_global_main.R",
        "dp_global_diag.R"
    )

    dp_files_all <- list.files(dp_r_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
    dp_files <- dp_files_all[basename(dp_files_all) %in% whitelist]

    for (f in sort(dp_files)) {
        message("Sourcing dp_global (core): ", basename(f))
        tryCatch(
            source(f, local = new.env()),
            error = function(e) {
                message("(non-fatal) error sourcing ", basename(f), ": ", e$message)
                NULL
            }
        )
    }
}

# Ensure required packages are available
if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it (install.packages('igraph')).")
}
library(igraph)

# Build a minimal tree_data example for a single Tag with two stems
# Census: 1..7 (no rows for census 5)
cdates <- as.Date(c("2000-01-01", "2001-01-01", "2002-01-01", "2003-01-01", "2004-01-01", "2005-01-01", "2006-01-01"))
rows <- list(
    # Census 1
    list(Tag = "T1", species = "sp1", CensusID = 1L, DBH = 10, ExactDate = cdates[1], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 1L, DBH = 8, ExactDate = cdates[1], TrueStemID = NA_integer_),
    # Census 2
    list(Tag = "T1", species = "sp1", CensusID = 2L, DBH = 11, ExactDate = cdates[2], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 2L, DBH = 9, ExactDate = cdates[2], TrueStemID = NA_integer_),
    # Census 3
    list(Tag = "T1", species = "sp1", CensusID = 3L, DBH = 12, ExactDate = cdates[3], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 3L, DBH = 10, ExactDate = cdates[3], TrueStemID = NA_integer_),
    # Census 4
    list(Tag = "T1", species = "sp1", CensusID = 4L, DBH = 13, ExactDate = cdates[4], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 4L, DBH = 11, ExactDate = cdates[4], TrueStemID = NA_integer_),
    # Census 5
    list(Tag = "T1", species = "sp1", CensusID = 5L, DBH = NA_real_, ExactDate = cdates[5], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 5L, DBH = NA_real_, ExactDate = cdates[5], TrueStemID = NA_integer_),
    # Census 6 (no observations recorded at census 5)
    list(Tag = "T1", species = "sp1", CensusID = 6L, DBH = 3.0, ExactDate = cdates[6], TrueStemID = NA_integer_),
    list(Tag = "T1", species = "sp1", CensusID = 6L, DBH = 2.5, ExactDate = cdates[6], TrueStemID = NA_integer_),
    # Anchor census 7 with known TrueStemID values (pins the reconstruction)
    list(Tag = "T1", species = "sp1", CensusID = 7L, DBH = 3.2, ExactDate = cdates[7], TrueStemID = 1L),
    list(Tag = "T1", species = "sp1", CensusID = 7L, DBH = 2.9, ExactDate = cdates[7], TrueStemID = 2L)
)

# Convert to data.table and add required bio parameters (one value per dataset expected by functions)
tree_data <- rbindlist(rows)

tree_data_incomplete <- copy(tree_data)
tree_data_incomplete[CensusID == 5, ExactDate := NA]
tree_data_incomplete[CensusID == 5, CensusID := NA]

# Bio params (simple illustrative values)
tree_data[, `:=`(
    Bio_Mu_Growth = 0.5,
    Bio_Gamma_Growth = 0,
    Bio_Sigma0_Growth = 0.1,
    Bio_Sigma1_Growth = 0,
    Bio_Max_Shrink = -10,
    Bio_K_Shrink = 0,
    Bio_Max_Growth = 10,
    Bio_Max_Growth_Soft = 5,
    Bio_K_Growth = 0,
    Bio_H0_Mortality = 0.01,
    Bio_Beta_Mortality = 0,
    Bio_Recruit_Meanlog = log(2),
    Bio_Recruit_Sdlog = 0.2,
    Bio_Recruit_MaxDBH_unit = 10,
    Bio_Recruitment_lambda = 0.2
)]

tree_data_incomplete[, `:=`(
    Bio_Mu_Growth = 0.5,
    Bio_Gamma_Growth = 0,
    Bio_Sigma0_Growth = 0.1,
    Bio_Sigma1_Growth = 0,
    Bio_Max_Shrink = -10,
    Bio_K_Shrink = 0,
    Bio_Max_Growth = 10,
    Bio_Max_Growth_Soft = 5,
    Bio_K_Growth = 0,
    Bio_H0_Mortality = 0.01,
    Bio_Beta_Mortality = 0,
    Bio_Recruit_Meanlog = log(2),
    Bio_Recruit_Sdlog = 0.2,
    Bio_Recruit_MaxDBH_unit = 10,
    Bio_Recruitment_lambda = 0.2
)]

# Run DP with anchor at census 7
message("Running DP with anchor_start=7 (verbose)")
res <- match_stems_dp_global_backward_marginals_batch(
    tree_data = copy(tree_data),
    anchor_start = 7L,
    max_tracks = 4L,
    slack_tracks = 0L,
    prune_hard = TRUE,
    verbose = TRUE,
    posterior_top_k = 2L
)

# NOTE: Even when no DBH, censuses should be continuous (no gaps)
message("Running DP with anchor_start=7 on data with missing census 5 (verbose)")
res_incomplete <- match_stems_dp_global_backward_marginals_batch(
    tree_data = copy(tree_data_incomplete),
    anchor_start = 7L,
    max_tracks = 4L,
    slack_tracks = 0L,
    prune_hard = TRUE,
    verbose = TRUE,
    posterior_top_k = 2L
)

# message("--- DP result (rows) ---")
# print(res[, .(CensusID, DBH, TrueStemID, ReconstructedStemID, ReconstructionMethod)])

# message("--- DP diagnostics ---")
# print(attr(res, "DP_PruneInfo"))

# message("Example finished. Note: census 5 had no observations; DP creates an empty state at that census and enforces life-cycle constraints (no return from Phase 2 -> Phase 1).")
