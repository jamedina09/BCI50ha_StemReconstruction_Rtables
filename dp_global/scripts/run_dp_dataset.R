############################################################
### run_dp_dataset.R — Run the DP on a pre-prepared dataset
############################################################
# This script loads a CSV that already contains observations and all
# Bio_* parameter columns, runs the backward DP per Tag×species group,
# and writes the reconstruction output.
#
# Usage:
#   Rscript dp_global/scripts/run_dp_dataset.R
#   Rscript dp_global/scripts/run_dp_dataset.R --INPUT_FILE=path/to/data.csv
#   Rscript dp_global/scripts/run_dp_dataset.R --WHICH_TAG=159367

############################################################
### 0) Housekeeping
############################################################
if (sys.nframe() == 0L) rm(list = ls())

############################################################
### 1) CLI parsing
############################################################
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    overrides <- list()
    for (a in args) {
        if (grepl("^--", a)) {
            kv <- sub("^--", "", a)
            eq <- regexpr("=", kv, fixed = TRUE)
            if (eq > 0L) {
                key <- toupper(substr(kv, 1L, eq - 1L))
                val <- substr(kv, eq + 1L, nchar(kv))
                overrides[[key]] <- val
            }
        }
    }
    overrides
}
cli <- parse_args()

############################################################
### 2) Defaults
############################################################
INPUT_FILE <- "data_simulation/data/stem_reconstruction_dp_global_rcpp.csv"
OUTPUT_DIR <- NULL            # auto-generated if NULL
WHICH_TAG  <- NULL            # NULL = run all tags
ANCHOR_START_CENSUS <- 7L

# DP tuning
MAX_GROWTH_FIXED   <- 7.5
MAX_SHRINK_FIXED   <- -0.5
RECRUIT_MAX_FIXED  <- (MAX_GROWTH_FIXED * 5) + 0.9999
DP_MAX_STATES      <- 40000L
DP_MAX_TRACKS      <- NULL    # auto
DP_SLACK_TRACKS    <- 1L
DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE <- TRUE
DP_SLACK_REQUIRE_ANCHOR_EPS <- 1e-6
DP_POSTERIOR_TOP_K <- 2L
USE_MEASUREMENT_ERROR <- TRUE
ALLOW_PROVISIONAL_DP_ANCHOR <- TRUE
DP_VERBOSE         <- TRUE

# Posterior sampling
POSTERIOR_SAMPLES         <- 200L
POSTERIOR_SAMPLES_FORMAT  <- "csv"
POSTERIOR_SAMPLE_SEED     <- 123L

# Output
WRITE_DP_CSV <- TRUE
WRITE_DP_RDS <- TRUE

############################################################
### 3) Apply CLI overrides
############################################################
coerce <- function(val, current) {
    if (is.integer(current))  return(suppressWarnings(as.integer(val)))
    if (is.numeric(current))  return(suppressWarnings(as.numeric(val)))
    if (is.logical(current))  return(as.logical(toupper(val)))
    val
}
for (nm in names(cli)) {
    if (exists(nm, inherits = FALSE)) {
        assign(nm, coerce(cli[[nm]], get(nm)))
    }
}

############################################################
### 4) Dependencies
############################################################
if (!requireNamespace("here", quietly = TRUE)) {
    here <- function(...) file.path(getwd(), ...)
} else {
    library(here)
}
setwd(here())

source(here("dp_global", "R", "dp_global_main.R"))

############################################################
### 5) Load data
############################################################
library(data.table)
cat("[run_dp_dataset] Loading:", INPUT_FILE, "\n")
xrun <- fread(INPUT_FILE)

# Ensure ExactDate is Date type (fread may read non-ISO formats as character)
if (is.character(xrun$ExactDate)) {
    # Try ISO first, then common US format
    parsed <- as.Date(xrun$ExactDate, format = "%Y-%m-%d")
    if (all(is.na(parsed))) {
        parsed <- as.Date(xrun$ExactDate, format = "%m/%d/%y")
    }
    xrun[, ExactDate := parsed]
    cat("[run_dp_dataset] Converted ExactDate from character to Date\n")
}

# Ensure species column
if (!"species" %in% names(xrun) && "Species" %in% names(xrun)) {
    setnames(xrun, "Species", "species")
}

# Filter to single tag if requested
if (!is.null(WHICH_TAG)) {
    WHICH_TAG <- suppressWarnings(as.integer(WHICH_TAG))
    if (!(WHICH_TAG %in% unique(xrun$Tag))) {
        stop("Requested WHICH_TAG=", WHICH_TAG, " not found in data.")
    }
    xrun <- xrun[Tag == WHICH_TAG]
    cat("[run_dp_dataset] Filtered to Tag=", WHICH_TAG, "\n")
}

cat("[run_dp_dataset] Rows:", nrow(xrun), " Tags:", length(unique(xrun$Tag)), "\n")

############################################################
### 6) Output directory
############################################################
if (is.null(OUTPUT_DIR)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    tag_label <- if (!is.null(WHICH_TAG)) paste0("T", WHICH_TAG) else "Tall"
    OUTPUT_DIR <- file.path("dp_global", "output",
                            paste0(ts, "_dataset_", tag_label, "_rcpp"))
}
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
posterior_path <- file.path(OUTPUT_DIR, "posteriors")
if (!dir.exists(posterior_path)) dir.create(posterior_path, recursive = TRUE)
cat("[run_dp_dataset] Output directory:", OUTPUT_DIR, "\n")

############################################################
### 7) Auto max tracks
############################################################
if (is.null(DP_MAX_TRACKS)) {
    if (exists("auto_dp_max_tracks", mode = "function")) {
        dp_max_tracks <- as.integer(auto_dp_max_tracks(xrun))
    } else {
        dp_max_tracks <- 20L
    }
} else {
    dp_max_tracks <- as.integer(DP_MAX_TRACKS)
}

############################################################
### 8) DP runner (per Tag×species group)
############################################################
run_dp_one_group <- function(dtg) {
    tag_label <- as.character(dtg$Tag[[1]])
    if (all(is.na(dtg$DBH)) || all(is.na(dtg$CensusID))) {
        cat("[run_dp_dataset] Skipping Tag=", tag_label, ": all DBH or CensusID NA\n")
        return(NULL)
    }

    match_stems_dp_global_backward_marginals_batch(
        tree_data         = copy(dtg),
        min_growth        = MAX_SHRINK_FIXED,
        max_growth        = MAX_GROWTH_FIXED,
        anchor_start      = ANCHOR_START_CENSUS,
        max_tracks        = dp_max_tracks,
        max_states        = DP_MAX_STATES,
        slack_tracks      = DP_SLACK_TRACKS,
        slack_require_anchor_recruitable = DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE,
        slack_require_anchor_eps         = DP_SLACK_REQUIRE_ANCHOR_EPS,
        temperature       = 1,
        posterior_top_k   = DP_POSTERIOR_TOP_K,
        fallback_growth_forms = character(0),
        posterior_samples        = POSTERIOR_SAMPLES,
        posterior_samples_format = POSTERIOR_SAMPLES_FORMAT,
        posterior_samples_path   = posterior_path,
        posterior_sample_seed    = POSTERIOR_SAMPLE_SEED,
        use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
        prune_hard            = TRUE,
        prune_min_growth      = MAX_SHRINK_FIXED * 2.5,
        prune_max_growth      = MAX_GROWTH_FIXED * 1.5,
        prune_use_bio_bounds  = FALSE,
        prune_recruit_max_dbh = RECRUIT_MAX_FIXED * 1.2,
        prune_use_bio_recruit = FALSE,
        allow_provisional_anchor = isTRUE(ALLOW_PROVISIONAL_DP_ANCHOR),
        verbose = isTRUE(DP_VERBOSE)
    )
}

############################################################
### 9) Run DP per group
############################################################
groups <- unique(xrun[, .(Tag, species)])
setorder(groups, Tag, species)

cat("[run_dp_dataset] Running DP on", nrow(groups), "group(s)...\n")

res_list <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    cat("\n========== Tag=", g$Tag, " species=", g$species, " ==========\n")
    dtg <- xrun[Tag == g$Tag & species == g$species]
    res_list[[i]] <- run_dp_one_group(dtg)
}

res_list <- Filter(Negate(is.null), res_list)
out <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

# Add posterior bins if available
if (exists("add_dp_posterior_bins", mode = "function")) {
    out <- add_dp_posterior_bins(out)
}

############################################################
### 10) Write output
############################################################
csv_path <- file.path(OUTPUT_DIR, "stem_reconstruction_dp_global_rcpp.csv")
rds_path <- file.path(OUTPUT_DIR, "stem_reconstruction_dp_global_rcpp.rds")

if (isTRUE(WRITE_DP_CSV)) {
    fwrite(out, csv_path)
    cat("[run_dp_dataset] Wrote CSV:", csv_path, "\n")
}
if (isTRUE(WRITE_DP_RDS)) {
    saveRDS(out, rds_path)
    cat("[run_dp_dataset] Wrote RDS:", rds_path, "\n")
}

############################################################
### 11) Summary
############################################################
cat("\n========== RECONSTRUCTION SUMMARY ==========\n")
print(out[, .(Tag, CensusID, DBH, ListOfTSM, TrueStemID, ReconstructedStemID, ReconstructionMethod)])
cat("\nDone.\n")
