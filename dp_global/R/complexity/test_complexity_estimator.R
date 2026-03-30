############################################################
### test_complexity_estimator.R
### Predict DP running time and show how parameter changes
### affect runtime.
###
### Produces:
###   1. Ranked tag list with predicted run time
###   2. Overall runtime summary (DP vs igraph, total hours)
###   3. Parameter sensitivity table showing how tightening
###      pruning or lowering max_states changes runtime
###      (skipped when --NO_SWEEP is passed)
###
### Parameters match main_cpp_chunk.R defaults.
### Run from the project root:
###   Rscript dp_global/R/complexity/test_complexity_estimator.R
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=bci_data/bci_multistem_xrun_debug.rds --NO_SWEEP
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=path/to/file.csv --ANCHOR_START=7 --DP_MAX_STATES=10000 --TOP_N=30
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=bci_data/bci_multistem_xrun_debug.rds --DP_MAX_STATES=10000 --NO_SWEEP --TOP_N=50
###   Rscript dp_global/R/complexity/test_complexity_estimator.R --INPUT_FILE=bci_data/bci_multistem_xrun_debug.rds --DP_MAX_STATES=10000 --NO_SWEEP --TOP_N=50 --SAMPLE_N=10

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
    } else if (grepl("^--[A-Za-z]", .a)) {
        # valueless flags like --NO_SWEEP treated as TRUE
        .key <- toupper(gsub("^--", "", .a))
        .overrides[[.key]] <- TRUE
    }
}

# ==================================================================
# Parameters — synced with main_cpp_chunk.R defaults
# ==================================================================
DATA_PATH        <- here("data_simulation", "data", "simulated_data_1.csv")
ANCHOR_START     <- 7L
DP_MAX_STATES    <- 40000L
SLACK_TRACKS     <- 1L
MAX_GROWTH_FIXED <- 5.0
MAX_SHRINK_FIXED <- -0.5
RECRUIT_MAX      <- (MAX_GROWTH_FIXED * 5) + 0.9999
PRUNE_MARGIN     <- 1.25   # prune bounds = fixed bounds * this margin
TOP_N            <- 30L    # number of slowest tags to display
RUN_SWEEP        <- TRUE   # set FALSE via --NO_SWEEP to skip the slow parameter sweep
SAMPLE_N         <- NULL   # if set, randomly sample this many tags (for quick checks)

# Apply CLI overrides
if (!is.null(.overrides$INPUT_FILE))      DATA_PATH        <- .overrides$INPUT_FILE
if (!is.null(.overrides$ANCHOR_START))    ANCHOR_START     <- as.integer(.overrides$ANCHOR_START)
if (!is.null(.overrides$DP_MAX_STATES))   DP_MAX_STATES    <- as.integer(.overrides$DP_MAX_STATES)
if (!is.null(.overrides$SLACK_TRACKS))    SLACK_TRACKS     <- as.integer(.overrides$SLACK_TRACKS)
if (!is.null(.overrides$MAX_GROWTH))      MAX_GROWTH_FIXED <- as.numeric(.overrides$MAX_GROWTH)
if (!is.null(.overrides$MIN_GROWTH))      MAX_SHRINK_FIXED <- as.numeric(.overrides$MIN_GROWTH)
if (!is.null(.overrides$RECRUIT_MAX))     RECRUIT_MAX      <- as.numeric(.overrides$RECRUIT_MAX)
if (!is.null(.overrides$PRUNE_MARGIN))    PRUNE_MARGIN     <- as.numeric(.overrides$PRUNE_MARGIN)
if (!is.null(.overrides$TOP_N))           TOP_N            <- as.integer(.overrides$TOP_N)
if (!is.null(.overrides$SAMPLE_N))        SAMPLE_N         <- as.integer(.overrides$SAMPLE_N)
if (isTRUE(.overrides$NO_SWEEP))          RUN_SWEEP        <- FALSE

# Recompute RECRUIT_MAX if MAX_GROWTH was overridden but RECRUIT_MAX was not
if (!is.null(.overrides$MAX_GROWTH) && is.null(.overrides$RECRUIT_MAX)) {
    RECRUIT_MAX <- (MAX_GROWTH_FIXED * 5) + 0.9999
}

# Resolve relative paths from project root
if (!file.exists(DATA_PATH)) {
    candidate <- file.path(here::here(), DATA_PATH)
    if (file.exists(candidate)) DATA_PATH <- candidate
}

# Derived prune bounds (match main_cpp_chunk.R)
PRUNE_MIN <- MAX_SHRINK_FIXED * PRUNE_MARGIN
PRUNE_MAX <- MAX_GROWTH_FIXED * PRUNE_MARGIN
PRUNE_REC <- RECRUIT_MAX * PRUNE_MARGIN

# ==================================================================
# BCI data preprocessing
# Detect if input is a BCI-style RDS (has raw mm DBH and StemID columns)
# and normalize to the standard columns expected by estimate_dp_complexity.
# ==================================================================
.is_bci_raw <- grepl("\\.rds$", DATA_PATH, ignore.case = TRUE)
if (.is_bci_raw) {
    cat("  Detected RDS input — loading and preprocessing BCI data...\n")
    .raw <- as.data.table(readRDS(DATA_PATH))

    # BCI-specific column mapping
    if ("dbh_with_best_candidate_taper_corrected" %in% names(.raw) &&
        !("DBH" %in% names(.raw))) {
        .raw[, DBH := dbh_with_best_candidate_taper_corrected / 10]  # mm -> cm
    }
    if ("StemID" %in% names(.raw) && !("TrueStemID" %in% names(.raw))) {
        .raw[, TrueStemID := as.integer(as.character(StemID))]
    }
    if ("Mnemonic" %in% names(.raw) && !("species" %in% names(.raw))) {
        .raw[, species := as.character(Mnemonic)]
    }
    .raw[, Tag := as.character(Tag)]

    # Keep only multi-stemmed tags (those that need DP)
    if ("single_stem_tags" %in% names(.raw)) {
        .n_before <- uniqueN(.raw$Tag)
        .raw <- .raw[single_stem_tags == FALSE]
        cat(sprintf("  Kept %d multi-stemmed tags (dropped %d single-stemmed).\n",
            uniqueN(.raw$Tag), .n_before - uniqueN(.raw$Tag)))
    }
    # NOTE: do NOT filter by CensusID here — the estimator needs all historical
    # censuses (not just >= anchor_start) to compute transition complexity.

    cat(sprintf("  Tags after filtering: %d\n", uniqueN(.raw$Tag)))
    DATA_PATH <- .raw   # pass the prepared data.table directly
    rm(.raw)
}
# Subsample tags (for quick smoke-tests)
if (!is.null(SAMPLE_N) && is.data.table(DATA_PATH)) {
    .all_tags <- unique(DATA_PATH$Tag)
    set.seed(42L)
    .keep <- sample(.all_tags, min(SAMPLE_N, length(.all_tags)))
    DATA_PATH <- DATA_PATH[Tag %in% .keep]
    cat(sprintf("  [SAMPLE_N=%d] Subsampled to %d tags.\n", SAMPLE_N, uniqueN(DATA_PATH$Tag)))
}
.data_path_label <- if (is.character(DATA_PATH)) DATA_PATH else
    sprintf("[BCI data.table: %d rows]", nrow(DATA_PATH))
cat("============================================================\n\n")
cat(sprintf("  Input file    : %s\n", .data_path_label))
cat(sprintf("  Anchor census : %d\n", ANCHOR_START))
cat(sprintf("  Max states    : %s\n", format(DP_MAX_STATES, big.mark = ",")))
cat(sprintf("  Growth bounds : [%.2f, %.2f] cm/yr\n", MAX_SHRINK_FIXED, MAX_GROWTH_FIXED))
cat(sprintf("  Prune bounds  : [%.2f, %.2f] cm/yr  (margin=%.2fx)\n", PRUNE_MIN, PRUNE_MAX, PRUNE_MARGIN))
cat(sprintf("  Recruit max   : %.2f cm  (prune: %.2f cm)\n", RECRUIT_MAX, PRUNE_REC))
cat("\n")

# ==================================================================
# 1. Run complexity estimation at current settings
# ==================================================================
cat("--- Running complexity estimation at current settings ---\n\n")

complexity <- estimate_dp_complexity(
    data                  = DATA_PATH,
    anchor_start          = ANCHOR_START,
    slack_tracks          = SLACK_TRACKS,
    max_states            = DP_MAX_STATES,
    min_growth            = MAX_SHRINK_FIXED,
    max_growth            = MAX_GROWTH_FIXED,
    prune_min_growth      = PRUNE_MIN,
    prune_max_growth      = PRUNE_MAX,
    prune_use_bio_bounds  = FALSE,
    recruit_max_dbh       = PRUNE_REC,
    prune_use_bio_recruit = FALSE,
    fast                  = TRUE
)

# ==================================================================
# 2. Overall Summary
# ==================================================================
n_total  <- nrow(complexity)
n_dp     <- sum(!complexity$estimated_fallback)
n_igraph <- sum( complexity$estimated_fallback)
total_h  <- sum(complexity$predicted_hours, na.rm = TRUE)
dp_only  <- complexity[estimated_fallback == FALSE]

cat("\n")
cat("============================================================\n")
cat("  OVERALL SUMMARY\n")
cat("============================================================\n")
cat(sprintf("  Total tags       : %d\n", n_total))
cat(sprintf("  Tags via DP      : %d  (%.1f%%)\n", n_dp, 100 * n_dp / max(1, n_total)))
cat(sprintf("  Tags via igraph  : %d  (%.1f%%)\n", n_igraph, 100 * n_igraph / max(1, n_total)))
cat(sprintf("  Total predicted  : %.2f hours\n", total_h))
if (nrow(dp_only) > 0L) {
    cat(sprintf("  Slowest tag      : %s (%.1f min)\n", dp_only$Tag[1L], dp_only$predicted_seconds[1L] / 60))
    cat(sprintf("  Median tag time  : %.1f sec\n", median(dp_only$predicted_seconds, na.rm = TRUE)))
}
cat("\n")

# ==================================================================
# 3. Top N Slowest Tags
# ==================================================================
display_cols <- c("Tag", "Species", "K", "max_obs", "n_censuses",
                  "max_states_per_census", "estimated_edges_unpruned",
                  "estimated_edges_pruned",
                  "estimated_fallback", "predicted_seconds", "predicted_hours")
display_cols <- display_cols[display_cols %in% names(complexity)]

n_show <- min(TOP_N, nrow(complexity))
disp <- complexity[seq_len(n_show), ..display_cols]
disp[, predicted_seconds := round(predicted_seconds, 1)]
disp[, predicted_hours   := round(predicted_hours, 4)]

cat("============================================================\n")
cat(sprintf("  TOP %d SLOWEST TAGS\n", n_show))
cat("============================================================\n")
print(disp, topn = n_show)
cat("\n")

# ==================================================================
# 4. Fallback breakdown (why tags use igraph)
# ==================================================================
if (n_igraph > 0L) {
    fb <- complexity[estimated_fallback == TRUE, .N, by = fallback_reason]
    setorder(fb, -N)
    cat("============================================================\n")
    cat("  IGRAPH FALLBACK REASONS\n")
    cat("============================================================\n")
    print(fb)
    cat("\n")
}

# ==================================================================
# 5. Parameter Sensitivity — how changing settings affects runtime
# ==================================================================
if (!isTRUE(RUN_SWEEP)) {
    cat("  (Parameter sweep skipped — pass without --NO_SWEEP to enable)\n\n")
} else {
cat("  PARAMETER SENSITIVITY\n")
cat("  How changing pruning bounds and max_states affects runtime\n")
cat("============================================================\n\n")

# Build scenarios. Each row overrides one parameter at a time from the current
# settings so you can see the marginal effect.
scenarios <- rbindlist(list(
    # Baseline (current settings)
    data.table(label = "current settings",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = PRUNE_MIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),

    # --- Vary max_states ---
    data.table(label = "max_states 5,000",
        max_states = 5000L, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = PRUNE_MIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),
    data.table(label = "max_states 10,000",
        max_states = 10000L, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = PRUNE_MIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),
    data.table(label = "max_states 100,000",
        max_states = 100000L, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = PRUNE_MIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),

    # --- Vary max_growth (tighter/wider) ---
    data.table(label = "max_growth 3 cm/yr",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = 3.0,
        prune_min_growth = MAX_SHRINK_FIXED * PRUNE_MARGIN, prune_max_growth = 3.0 * PRUNE_MARGIN,
        recruit_max_dbh = ((3.0 * 5) + 0.9999) * PRUNE_MARGIN),
    data.table(label = "max_growth 7.5 cm/yr",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = 7.5,
        prune_min_growth = MAX_SHRINK_FIXED * PRUNE_MARGIN, prune_max_growth = 7.5 * PRUNE_MARGIN,
        recruit_max_dbh = ((7.5 * 5) + 0.9999) * PRUNE_MARGIN),
    data.table(label = "max_growth 10 cm/yr",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = 10.0,
        prune_min_growth = MAX_SHRINK_FIXED * PRUNE_MARGIN, prune_max_growth = 10.0 * PRUNE_MARGIN,
        recruit_max_dbh = ((10.0 * 5) + 0.9999) * PRUNE_MARGIN),

    # --- Vary shrinkage ---
    data.table(label = "min_growth -0.25 cm/yr",
        max_states = DP_MAX_STATES, min_growth = -0.25, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = -0.25 * PRUNE_MARGIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),
    data.table(label = "min_growth -1.0 cm/yr",
        max_states = DP_MAX_STATES, min_growth = -1.0, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = -1.0 * PRUNE_MARGIN, prune_max_growth = PRUNE_MAX,
        recruit_max_dbh = PRUNE_REC),

    # --- Tighter prune margin ---
    data.table(label = "prune margin 1.0x (no margin)",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = MAX_SHRINK_FIXED * 1.0, prune_max_growth = MAX_GROWTH_FIXED * 1.0,
        recruit_max_dbh = RECRUIT_MAX * 1.0),
    data.table(label = "prune margin 1.5x",
        max_states = DP_MAX_STATES, min_growth = MAX_SHRINK_FIXED, max_growth = MAX_GROWTH_FIXED,
        prune_min_growth = MAX_SHRINK_FIXED * 1.5, prune_max_growth = MAX_GROWTH_FIXED * 1.5,
        recruit_max_dbh = RECRUIT_MAX * 1.5)
), fill = TRUE)

sweep <- sweep_dp_complexity(
    data        = DATA_PATH,
    scenarios   = scenarios,
    base_params = list(
        anchor_start       = ANCHOR_START,
        slack_tracks       = SLACK_TRACKS,
        prune_use_bio_bounds  = FALSE,
        prune_use_bio_recruit = FALSE
    )
)

# Display sweep results
show_cols <- c("label", "n_tags_dp", "n_tags_igraph", "pct_dp",
               "total_hours", "slowest_min", "median_sec")
show_cols <- show_cols[show_cols %in% names(sweep)]
cat("\n")
print(sweep[, ..show_cols], topn = nrow(sweep))

cat("\n")
cat("  Notes:\n")
cat("  - 'total_hours' is the sum of predicted DP time for all tags\n")
cat("  - 'slowest_min' is the predicted time for the single slowest tag\n")
cat("  - 'n_tags_igraph' = tags that exceed max_states and fall back\n")
cat("  - Tighter pruning reduces edges but does not change the state count;\n")
cat("    lowering max_states pushes more tags to igraph fallback\n")
cat("\n")

} # end if (RUN_SWEEP)

# ==================================================================
# 6. Export results
# ==================================================================
output_path <- if (is.character(DATA_PATH)) {
    sub("\\.[^.]+$", "_complexity_report.csv", DATA_PATH)
} else {
    file.path(here("dp_global", "R", "complexity"), "complexity_report.csv")
}
fwrite(complexity, output_path)
cat(sprintf("  Full per-tag results exported to: %s\n", output_path))

if (isTRUE(RUN_SWEEP)) {
    sweep_path <- if (is.character(DATA_PATH)) {
        sub("\\.[^.]+$", "_complexity_sweep.csv", DATA_PATH)
    } else {
        file.path(here("dp_global", "R", "complexity"), "complexity_sweep.csv")
    }
    fwrite(sweep, sweep_path)
    cat(sprintf("  Parameter sweep results exported to: %s\n\n", sweep_path))
}
