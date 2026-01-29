rm(list = ls())
# run_attach_tag41.R
# -----------------------------------------------------------------------------
# Small script to:
#  1) Run a DP chunk (defaults to chunk 6) with posterior sampling,
#  2) find the output directory produced,
#  3) read the paths summary for a target Tag (defaults to 41),
#  4) attach top-N posterior paths to the tag's rows and merge back into the
#     full DP output CSV, and
#  5) optionally expand draws and aggregate per-ObsRowID frequencies.
#
# Usage:
#  - Source interactively in RStudio to step through line-by-line, or run non-
#    interactively as `Rscript dp_global/dev/run_attach_tag41.R` (defaults will run).
#  - Edit the configuration block below to change behavior.
# -----------------------------------------------------------------------------

# --- Configuration (edit interactively when sourcing) ------------------------
CHUNK_START <- 6L
CHUNK_END   <- 6L
POSTERIOR_SAMPLES <- 200L
POSTERIOR_SAMPLES_FORMAT <- "csv"
DP_POSTERIOR_TOP_K <- 3L
MANUAL_CORES <- TRUE
MANUAL_CORES_VALUE <- 1L
WRITE_DP_FEATHER <- FALSE
WRITE_DP_PDF <- FALSE

TARGET_TAGS <- c(41L)
ATTACH_N <- 7L        # top-N per tag
ATTACH_WHICH <- "top_n" # options: "top_n", "map", "indices", "sample"
RUN_EXPAND <- TRUE
EXPAND_N <- 200L
EXPAND_OUT_DIR <- file.path("dp_global", "output")

# -----------------------------------------------------------------------------

# Safety: working directory should be project root
stopifnot(file.exists("dp_global/scripts/main_cpp_chunk.R"))

# Load required code and helpers
message("Sourcing runner and posterior helpers...")
source("dp_global/scripts/main_cpp_chunk.R")
source("dp_global/R/error_propagation/process_posteriors.R")

# Export configuration into the global environment used by the runner
assign("DP_CHUNK_START", as.integer(CHUNK_START), envir = .GlobalEnv)
assign("DP_CHUNK_END", as.integer(CHUNK_END), envir = .GlobalEnv)
assign("POSTERIOR_SAMPLES", as.integer(POSTERIOR_SAMPLES), envir = .GlobalEnv)
assign("POSTERIOR_SAMPLES_FORMAT", POSTERIOR_SAMPLES_FORMAT, envir = .GlobalEnv)
assign("DP_POSTERIOR_TOP_K", as.integer(DP_POSTERIOR_TOP_K), envir = .GlobalEnv)
assign("MANUAL_CORES", as.logical(MANUAL_CORES), envir = .GlobalEnv)
assign("MANUAL_CORES_VALUE", as.integer(MANUAL_CORES_VALUE), envir = .GlobalEnv)
assign("WRITE_DP_FEATHER", as.logical(WRITE_DP_FEATHER), envir = .GlobalEnv)
assign("WRITE_DP_PDF", as.logical(WRITE_DP_PDF), envir = .GlobalEnv)

message(sprintf("Running DP chunk %d..%d with posterior samples=%d (format=%s)", CHUNK_START, CHUNK_END, POSTERIOR_SAMPLES, POSTERIOR_SAMPLES_FORMAT))

# Run chunked DP (will create output directory + posteriors)
run_main_chunked()

# Locate the most recent output directory under dp_global/output
out_dirs <- list.dirs("dp_global/output", full.names = TRUE, recursive = FALSE)
if (length(out_dirs) == 0L) stop("No output directories found in dp_global/output")
latest_out <- tail(sort(out_dirs), 1)
message("Using output dir: ", latest_out)

# Load the main chunk CSV produced by the runner
out_file <- file.path(latest_out, "stem_reconstruction_dp_global_rcpp.csv")
if (!file.exists(out_file)) stop("DP output CSV not found: ", out_file)
out_dt <- data.table::fread(out_file)

# Attach paths for each requested tag
# for (tag in TARGET_TAGS) {
tag <- TARGET_TAGS[1]
    message("Processing Tag=", tag)
    post_dir <- file.path(latest_out, "posteriors")
    if (!dir.exists(post_dir)) stop("Posterior directory not found: ", post_dir)

    # Heuristic: find the most recently written "paths" file for the tag
    pattern <- paste0("^tag_", tag, "_posterior_samples.*_paths\\.(csv|feather|rds)$")
    cand <- list.files(post_dir, pattern = pattern, full.names = TRUE)
    if (length(cand) == 0L) stop("No posterior paths file found for Tag=", tag, " in ", post_dir)
    paths_file <- tail(sort(cand), 1)
    message("Found paths file: ", paths_file)

    # Load paths
    paths_dt <- load_posterior_paths(paths_file)

    # Subset the DP output for the tag
    out_sub <- out_dt[Tag == tag]
    if (nrow(out_sub) == 0L) {
        warning("No rows found for Tag=", tag, " in DP output; skipping")
        next
    }

    # Attach selected paths to subset
    message("Attaching paths to Tag=", tag, " (which=", ATTACH_WHICH, ", n=", ATTACH_N, ")")
    out_sub_res <- attach_paths_to_output(paths_dt, out_sub, which = ATTACH_WHICH, n = as.integer(ATTACH_N), write_out = FALSE)

    # Merge newly created DP_ReconstructedStemID_* and DP_PathSig_* columns back into full output
    new_cols <- grep("^DP_ReconstructedStemID_|^DP_PathSig_", names(out_sub_res), value = TRUE)
    if (length(new_cols) > 0L) {
        # Ensure obs_row_id exists in both tables
        if (!("obs_row_id" %in% names(out_dt)) || !("obs_row_id" %in% names(out_sub_res))) stop("Missing obs_row_id in output or subset; cannot merge")
        # Map rows by obs_row_id and assign values by position to avoid NSE pitfalls
        ix_out <- match(out_sub_res$obs_row_id, out_dt$obs_row_id)
        for (c in new_cols) {
            out_dt[ix_out, (c) := out_sub_res[[c]]]
        }
        message("Merged ", length(new_cols), " new columns into full output for Tag=", tag)
    } else {
        message("No new DP_* columns were produced for Tag=", tag)
    }

    # Optionally expand draws and aggregate
    if (isTRUE(RUN_EXPAND)) {
        # Find summary file (matching the same tag)
        summary_pattern <- paste0("^tag_", tag, "_posterior_samples.*_summary\\.(csv|feather|rds)$")
        summary_cand <- list.files(post_dir, pattern = summary_pattern, full.names = TRUE)
        if (length(summary_cand) == 0L) {
            warning("No posterior summary file found for Tag=", tag, "; skipping expand_draws")
        } else {
            summary_file <- tail(sort(summary_cand), 1)
            message("Found summary file: ", summary_file, "; expanding N=", EXPAND_N)
            summary_dt <- data.table::fread(summary_file)
            res_dt <- expand_draws(summary_dt, paths_dt, N = as.integer(EXPAND_N))
            agg_dt <- aggregate_draws(res_dt)
            expand_out_file <- file.path(latest_out, paste0("expanded_draws_tag", tag, "_N", EXPAND_N, ".csv"))
            data.table::fwrite(agg_dt, expand_out_file)
            message("Wrote aggregated draws to: ", expand_out_file)
        }
    }

    # Write an updated output file for this tag merge
    out_file_new <- file.path(latest_out, paste0("stem_reconstruction_dp_global_rcpp_with_paths_tag", tag, ".csv"))
    data.table::fwrite(out_dt, out_file_new)
    message("Wrote merged DP output with attached paths: ", out_file_new)
# }

message("Done.")
