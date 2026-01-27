# attach_paths_to_output_run.R
# ---------------------------------------
# Refactored runner: clearly separated functions and a small runner that
# demonstrates attachment of DP posterior paths and how to expand draws.
#
# Notes:
# - Keep `process_posteriors.R` as the canonical implementation file
#   (parsing, attaching). This runner sources it and calls its exported
#   helpers (`load_posterior_paths`, `attach_paths_to_output`).
# - This file is safe to source interactively or run via `Rscript`.

rm(list = ls())

library(data.table)
source("dp_global/R/error_propagation/process_posteriors.R")

vcat <- function(...) cat(..., "\n")

# -----------------------------------------------------------------------------
# Small library of convenience functions
# -----------------------------------------------------------------------------
check_inputs <- function(paths_file, out_file) {
    if (!file.exists(paths_file)) stop("paths file not found: ", paths_file)
    if (!file.exists(out_file)) stop("output file not found: ", out_file)
    invisible(TRUE)
}

# Use `load_posterior_paths()` from `process_posteriors.R` (no local wrapper).

# Use `attach_paths_to_output()` directly from `process_posteriors.R` (no wrapper in the runner).

# Use `expand_draws()` from `process_posteriors.R` to expand summary+paths into per-draw long table.

# Use `aggregate_draws()` from `process_posteriors.R` to aggregate expanded draws into per-ObsRowID probabilities.

# Use `check_map_in_paths()` from `process_posteriors.R` for MAP diagnostics.

# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------
# This script is intended to be sourced in RStudio or run interactively via
# `source()`.

# Configuration defaults (edit interactively as needed)
default_paths <- "./dp_global/output/20260123_184403_unknown_T19_DP_MB_ME_gD_sD_kgD_ksD_rcpp/posteriors/posteriors/tag_19_posterior_samples_20260123_184616_paths.csv"
default_out <- "./dp_global/output/20260123_184403_unknown_T19_DP_MB_ME_gD_sD_kgD_ksD_rcpp/stem_reconstruction_dp_global_rcpp.csv"

# Main interactive runner (explicit arguments, easy to call from R)
run_attach_and_expand <- function(paths_file = default_paths,
                                  out_file = default_out,
                                  attach_n = 2L,
                                  write_out = TRUE,
                                  run_expand_draws = FALSE,
                                  expand_N = 1000L,
                                  expand_out_csv = "./dp_global/output/expanded_draws_sample_1000.csv") {
    check_inputs(paths_file, out_file)
    paths_dt <- load_posterior_paths(paths_file)

    out_dt <- attach_paths_to_output(paths_dt, out_file, which = "top_n", n = as.integer(attach_n), write_out = write_out)

    # diagnostic
    m <- check_map_in_paths(paths_dt, out_dt)
    if (isTRUE(m$found)) {
        vcat("MAP joint path was found among sampled unique paths (idx):", m$idx)
    } else if (!is.null(m$best)) {
        vcat(sprintf("MAP not found; best partial idx=%d frac=%.3f prob=%.5f", m$best$idx, m$best$frac, m$best$prob))
        vcat("running more samples may help capture MAP path")
    }

    if (run_expand_draws) {
        summary_file <- sub("_paths\\.csv$", "_summary.csv", paths_file)
        if (!file.exists(summary_file)) stop("summary file not found: ", summary_file)
        summary_dt <- fread(summary_file)
        res_dt <- expand_draws(summary_dt, paths_dt, N = as.integer(expand_N))
        agg <- aggregate_draws(res_dt)
        vcat("Writing aggregated draws to:", expand_out_csv)
        fwrite(agg, expand_out_csv)
        vcat("Wrote", nrow(agg), "rows of aggregated frequencies.")
    }

    invisible(list(paths_dt = paths_dt, out_dt = out_dt))
}

# Run automatically when sourced interactively with defaults
if (interactive()) {
    vcat("Interactive: running `run_attach_and_expand()` with defaults (no expand_draws)")
    run_attach_and_expand()
}
