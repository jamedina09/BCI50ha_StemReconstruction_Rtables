############################################################
### test_complexity_estimator.R
### Rank tags by expected DP run time before batch submission.
###
### Parameters match main_cpp_chunk.R defaults.
### Run from the project root:
###   Rscript dp_global/R/complexity/test_complexity_estimator.R
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=bci_data/bci_multistem_xrun_debug.rds
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=path/to/file.csv --ANCHOR_START=7 --DP_MAX_STATES=1100
############################################################
rm(list = ls())

library(here)
library(data.table)

source(here("dp_global", "R", "complexity", "estimate_dp_complexity_function.R"))

# --- CLI argument parsing ---
.args <- commandArgs(trailingOnly = TRUE)
.overrides <- list()
for (.a in .args) {
    if (grepl("^--[A-Za-z].*=", .a)) {
        .kv  <- strsplit(sub("^--", "", .a), "=", fixed = TRUE)[[1L]]
        .key <- toupper(gsub("-", "_", .kv[1L]))
        .val <- paste(.kv[-1L], collapse = "=")
        if (grepl("^[+-]?[0-9]+$", .val))            .val <- as.integer(.val)
        else if (grepl("^[+-]?[0-9]*\\.[0-9]+$", .val)) .val <- as.numeric(.val)
        else if (tolower(.val) %in% c("true","false")) .val <- as.logical(tolower(.val))
        .overrides[[.key]] <- .val
    }
}

# --- Parameters: keep in sync with main_cpp_chunk.R ---
DATA_PATH     <- here("data_simulation", "data", "simulated_data_1.csv")
ANCHOR_START  <- 7L
DP_MAX_STATES <- 1100L   # --DP_MAX_STATES in main_cpp_chunk.R
SLACK_TRACKS  <- 1L
MIN_GROWTH    <- -0.5    # MAX_SHRINK_FIXED
MAX_GROWTH    <- 5.0     # MAX_GROWTH_FIXED
RECRUIT_MAX   <- (MAX_GROWTH * 5) + 0.9999  # RECRUIT_MAX_FIXED

# Apply CLI overrides
if (!is.null(.overrides$INPUT_FILE))   DATA_PATH     <- .overrides$INPUT_FILE
if (!is.null(.overrides$ANCHOR_START)) ANCHOR_START  <- as.integer(.overrides$ANCHOR_START)
if (!is.null(.overrides$DP_MAX_STATES)) DP_MAX_STATES <- as.integer(.overrides$DP_MAX_STATES)
if (!is.null(.overrides$SLACK_TRACKS)) SLACK_TRACKS  <- as.integer(.overrides$SLACK_TRACKS)
if (!is.null(.overrides$MIN_GROWTH))   MIN_GROWTH    <- as.numeric(.overrides$MIN_GROWTH)
if (!is.null(.overrides$MAX_GROWTH))   MAX_GROWTH    <- as.numeric(.overrides$MAX_GROWTH)
if (!is.null(.overrides$RECRUIT_MAX))  RECRUIT_MAX   <- as.numeric(.overrides$RECRUIT_MAX)
# Resolve relative paths from project root
if (!file.exists(DATA_PATH)) {
    candidate <- file.path(here::here(), DATA_PATH)
    if (file.exists(candidate)) DATA_PATH <- candidate
}
RECRUIT_MAX <- (MAX_GROWTH * 5) + 0.9999   # recompute if MAX_GROWTH changed and RECRUIT_MAX not set
if (!is.null(.overrides$RECRUIT_MAX)) RECRUIT_MAX <- as.numeric(.overrides$RECRUIT_MAX)

# -------------------------------------------------------
cat(sprintf("[complexity] Input: %s\n", DATA_PATH))
cat(sprintf("[complexity] anchor=%d  max_states=%d  min_growth=%.2f  max_growth=%.2f  recruit_max=%.2f\n\n",
    ANCHOR_START, DP_MAX_STATES, MIN_GROWTH, MAX_GROWTH, RECRUIT_MAX))

complexity <- estimate_dp_complexity(
    data                  = DATA_PATH,
    anchor_start          = ANCHOR_START,
    slack_tracks          = SLACK_TRACKS,
    max_states            = DP_MAX_STATES,
    min_growth            = MIN_GROWTH,
    max_growth            = MAX_GROWTH,
    prune_min_growth      = MIN_GROWTH * 2.5,
    prune_max_growth      = MAX_GROWTH * 1.5,
    prune_use_bio_bounds  = FALSE,
    recruit_max_dbh       = RECRUIT_MAX * 1.2,
    prune_use_bio_recruit = FALSE,
    fast = TRUE
)

# -------------------------------------------------------
# Summary counts
n_dp     <- sum(!complexity$estimated_fallback)
n_igraph <- sum( complexity$estimated_fallback)
total_h  <- round(sum(complexity$predicted_hours, na.rm = TRUE), 2)

cat(sprintf("Tags via DP     : %d\n", n_dp))
cat(sprintf("Tags via igraph : %d\n", n_igraph))
cat(sprintf("Total predicted : %.2f hours\n\n", total_h))

# -------------------------------------------------------
# Clean ranked display (already sorted slowest first by estimate_dp_complexity)
display_cols <- c("Tag", "K", "max_obs", "n_censuses",
                  "max_states_per_census", "estimated_edges_unpruned",
                  "estimated_fallback", "predicted_seconds", "predicted_hours")

disp <- complexity[, display_cols[display_cols %in% names(complexity)], with = FALSE]
disp[, predicted_seconds := round(predicted_seconds, 1)]
disp[, predicted_hours   := round(predicted_hours,   4)]

cat("=== All tags ranked slowest first ===\n")
print(disp, topn = nrow(disp))

cat(sprintf("\n=== Top 10 slowest tags ===\n"))
print(disp[seq_len(min(10L, nrow(disp)))], topn = 10L)

# -------------------------------------------------------
# Export
output_path <- here("data_simulation", "data", "report_run_simulated_data_1.csv")
fwrite(complexity, output_path)
cat(sprintf("\nResults exported to: %s\n", output_path))
