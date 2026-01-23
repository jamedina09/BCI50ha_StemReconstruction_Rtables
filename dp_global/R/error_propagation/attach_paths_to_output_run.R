
#!/usr/bin/env Rscript

# attach_paths_to_output_run.R
# ---------------------------
# Demonstration script that attaches DP posterior path reconstructions to the
# main reconstruction output file using the per-row identifier `obs_row_id`.
#
# This script is intentionally simple and documented to make the workflow clear:
# - Load posterior paths (via `load_posterior_paths()`)
# - Attach chosen paths to the main reconstruction output (`attach_paths_to_output()`)
# - Write out updated CSV (appended filename with `_with_paths.csv`)
# - Create a small `simp` data.frame and plot original vs. attached DP reconstructions
#
# Usage:
#   Rscript dp_global/R/error_propagation/attach_paths_to_output_run.R
#
# Files (expected to already exist in this workspace):
# - paths_file: the DP posterior paths summary file (CSV produced by DP)
# - out_file: the main per-row reconstruction output CSV

# Load core package (kept minimal, no startup messages here)
library(data.table)

# Source the tiny helper module that contains the functions used below:
# - load_posterior_paths()
# - attach_paths_to_output()
source("dp_global/R/error_propagation/reconstruct_and_propagate.R")

# -----------------------------------------------------------------------------
# Configuration: change these to point at the files you want to process
# -----------------------------------------------------------------------------
paths_file <- "dp_global/output/20260123_005933_fixed_T20_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/posteriors/tag_20_posterior_samples_20260123_005948_paths.csv"
out_file <- "dp_global/output/20260123_005933_fixed_T20_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv"

# Quick existence check to fail early with a helpful message
if (!file.exists(paths_file) || !file.exists(out_file)) {
    stop("Required files not found in expected locations: check 'paths_file' and 'out_file' variables at the top of this script")
}

# -----------------------------------------------------------------------------
# Load posterior paths and attach to the main output
# -----------------------------------------------------------------------------
# `load_posterior_paths()` will return a data.table with parsed reconstructions
paths_dt <- load_posterior_paths(paths_file)

# Attach the top-2 most-probable unique paths to the main output by obs_row_id.
# This will add columns to the output table:
#   DP_PathSig_1, DP_ReconstructedStemID_1, DP_PathSig_2, DP_ReconstructedStemID_2
# and will write a new CSV to disk with suffix `_with_paths.csv`.
out_dt <- attach_paths_to_output(paths_dt, out_file, which = "top_n", n = 2, write_out = TRUE)
cat("Attached top 2 paths; new columns:", paste0("DP_PathSig_", 1:2), "and", paste0("DP_ReconstructedStemID_", 1:2), "\n")

# Show a short head of the updated table so you can inspect the mapping quickly
print(head(out_dt[, .(DBH, CensusID, obs_row_id, DP_PathSig_1, DP_ReconstructedStemID_1, DP_PathSig_2, DP_ReconstructedStemID_2)]))

# -----------------------------------------------------------------------------
# Small tidy dataset for plotting and exploration
# -----------------------------------------------------------------------------
# Create a simplified view containing only the columns we need for plotting
simp <- out_dt[, .(ReconstructedStemID, DBH, CensusID, obs_row_id, DP_PathSig_1, DP_ReconstructedStemID_1, DP_PathSig_2, DP_ReconstructedStemID_2)]

# Example behavior: for censuses 7..9, if the DP-reconstructed mapping is NA we
# fall back to the (existing) `ReconstructedStemID` value for plotting. This is
# purely a visualization convenience and does not modify the main output on disk.
simp <- simp[CensusID %in% c(7, 8, 9), DP_ReconstructedStemID_1 := fifelse(is.na(DP_ReconstructedStemID_1), ReconstructedStemID, DP_ReconstructedStemID_1)]
simp <- simp[CensusID %in% c(7, 8, 9), DP_ReconstructedStemID_2 := fifelse(is.na(DP_ReconstructedStemID_2), ReconstructedStemID, DP_ReconstructedStemID_2)]

# -----------------------------------------------------------------------------
# Plotting: original vs DP-mapped trajectories
# -----------------------------------------------------------------------------
# Use ggplot2 for quick visual checks; if it's not available we print a short note
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' not available; skipping example plots. Install it to enable plots.")
} else {
    library(ggplot2)

    # Original (ReconstructedStemID) trajectories
    p1 <- ggplot(simp, aes(y = DBH, x = as.factor(CensusID), color = as.factor(ReconstructedStemID), group = as.factor(ReconstructedStemID))) +
        geom_point(size = 2, alpha = 0.7) +
        geom_line() +
        theme_minimal() +
        labs(title = "Original ReconstructedStemID trajectories", x = "CensusID", color = "ReconstructedStemID")
    print(p1)

    # DP mapping: top path 1
    p2 <- ggplot(simp, aes(y = DBH, x = as.factor(CensusID), color = as.factor(DP_ReconstructedStemID_1), group = as.factor(DP_ReconstructedStemID_1))) +
        geom_point(size = 2, alpha = 0.7) +
        geom_line() +
        theme_minimal() +
        labs(title = "DP Path 1 mapping", x = "CensusID", color = "DP_ReconstructedStemID_1")
    print(p2)

    # DP mapping: top path 2
    p3 <- ggplot(simp, aes(y = DBH, x = as.factor(CensusID), color = as.factor(DP_ReconstructedStemID_2), group = as.factor(DP_ReconstructedStemID_2))) +
        geom_point(size = 2, alpha = 0.7) +
        geom_line() +
        theme_minimal() +
        labs(title = "DP Path 2 mapping", x = "CensusID", color = "DP_ReconstructedStemID_2")
    print(p3)
}

# End of script: use invisibly return for convenience when sourced
invisible(NULL)

paths_file <- "dp_global/output/20260123_005933_fixed_T20_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/posteriors/tag_20_posterior_samples_20260123_005948_paths.csv"
out_file <- "dp_global/output/20260123_005933_fixed_T20_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv"

if (!file.exists(paths_file) || !file.exists(out_file)) stop("Required files not found in expected locations")

paths_dt <- load_posterior_paths(paths_file)
# Attach top 2 paths to the main output and write to disk
out_dt <- attach_paths_to_output(paths_dt, out_file, which = "top_n", n = 2, write_out = TRUE)
cat("Attached top 2 paths; new columns:", paste0("DP_PathSig_", 1:2), "and", paste0("DP_ReconstructedStemID_", 1:2), "\n")
print(head(out_dt[, .(DBH, CensusID, obs_row_id, DP_PathSig_1, DP_ReconstructedStemID_1, DP_PathSig_2, DP_ReconstructedStemID_2)]))

simp <- out_dt[, .(ReconstructedStemID, DBH, CensusID, obs_row_id, DP_PathSig_1, DP_ReconstructedStemID_1, DP_PathSig_2, DP_ReconstructedStemID_2)]
simp <- simp[CensusID %in% c(7, 8, 9), DP_ReconstructedStemID_1 := fifelse(is.na(DP_ReconstructedStemID_1), ReconstructedStemID, DP_ReconstructedStemID_1)]
simp <- simp[CensusID %in% c(7, 8, 9), DP_ReconstructedStemID_2 := fifelse(is.na(DP_ReconstructedStemID_2), ReconstructedStemID, DP_ReconstructedStemID_2)]

library(ggplot2)

ggplot(simp, aes(
    y = DBH, x = as.factor(CensusID),
    color = as.factor(ReconstructedStemID),
    group = as.factor(ReconstructedStemID)
)) +
    geom_point(, size = 2, alpha = 0.7) +
    geom_line() +
    theme_minimal()

ggplot(simp, aes(
    y = DBH, x = as.factor(CensusID),
    color = as.factor(DP_ReconstructedStemID_1),
    group = as.factor(DP_ReconstructedStemID_1)
)) +
    geom_point(, size = 2, alpha = 0.7) +
    geom_line() +
    theme_minimal()

ggplot(simp, aes(
    y = DBH, x = as.factor(CensusID),
    color = as.factor(DP_ReconstructedStemID_2),
    group = as.factor(DP_ReconstructedStemID_2)
)) +
    geom_point(, size = 2, alpha = 0.7) +
    geom_line() +
    theme_minimal()
