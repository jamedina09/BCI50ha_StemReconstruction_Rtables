############################################################
### main_cpp_chunk_bci.R — dp_global BCI chunked driver
############################################################
# Goal
#   BCI-specific variant of the chunked DP_GLOBAL pipeline. Sources the DP
#   solver from a self-contained bundle (dp_bundle_path) via withr::with_dir
#   so bundle-internal relative paths resolve correctly when running from
#   another directory. Groups (Tag + species) are processed in chunks of
#   DP_CHUNK_SIZE and outputs are written incrementally to disk.
#
# BCI-specific behaviour (vs. main_cpp_chunk.R)
# - Loads data from an RDS file (readRDS) instead of a CSV (fread).
# - Estimates biological parameters from post-anchor census data, separately
#   for trees/shrubs, palms/tree-ferns, figs, strangler figs, and unknowns.
# - Splits tags into single-stemmed (bypass) and multi-stemmed (DP) subsets.
# - FORCE_ONE_SPECIES_PARAMETERS=FALSE — uses real BCI species identities.
# - Uses withr::with_dir to source dp_global_main.R from the bundle.
#
# Note for orchestrators
# - This script accepts CLI overrides of internal variables via --KEY=VALUE.
# - See the `CLI_REFERENCE` variable below for the canonical keys used by
#   external orchestrators; they should construct flags matching these names
#   (case-insensitive, '-' or '_' allowed).
#
# Differences from main_cpp.R
# - run_main_chunked() writes each chunk to disk and sets out <- NULL after
#   each chunk, keeping peak memory proportional to chunk size.
# - Sensitivity sweeps and realism reports are disabled (require a full out).
# - PLOT_PDF_ONE_TAG_ONLY is not used; per-chunk PDFs are controlled by
#   WRITE_DP_PDF_PER_CHUNK.
#
# Table of Contents
#  0) Housekeeping       — workspace reset guard
#  1) CLI parsing        — parse_args(); overrides applied later in section 4
#  2) Dependencies       — package checks and imports
#  3) Defaults           — editable run defaults
#    3.1) Biological parameter estimation settings
#    3.2) DP solver settings
#    3.3) Chunking & posterior sampling settings
#    3.4) Parallelism settings
#    3.5) Output & path settings
#  4) CLI reference & override mapping — CLI_REFERENCE list; apply overrides
#  5) Input validation   — validate INPUT_FILE; derive boolean flags; compute out_dir
#  6) Source project code — load dp_global R modules
#  7) Helpers            — filesystem, logging, data-manipulation utilities
#  8) Core DP functions  — run_dp_one_group(); maybe_add_posterior_bins()
#  9) Main pipeline      — run_main_chunked(); chunk loop; incremental writes
# 10) Post-run utilities — merge_chunk_rds_to_csv(); merge_chunks_to_csv()
# 11) Entrypoint         — execute run_main_chunked() when called via Rscript

############################################################
### 0) Housekeeping
############################################################
# Avoid nuking the user's interactive environment when sourcing this file.
# Only clear the workspace when running as a top-level script.
if (sys.nframe() == 0L) {
    rm(list = ls())
}

############################################################
### 1) CLI parsing — Command-line parsing & overrides
############################################################
# Parse command-line arguments to override defaults.
# Usage: Rscript main_cpp_chunk.R --RUN_ALL_TAGS=TRUE --DP_CHUNK_SIZE=7 --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=4
# Supported args: any config variable name prefixed with --

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    overrides <- list()
    for (arg in args) {
        # Short help flag
        if (identical(arg, "-h") || identical(arg, "--help")) {
            overrides[["help"]] <- TRUE
            next
        }

        if (grepl("^--", arg)) {
            # --key=value or --FLAG
            if (grepl("=", arg)) {
                kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
                key <- kv[1]
                val <- kv[2]
                # Try to convert to appropriate type (handle booleans, integers, floats, including negatives)
                # Keys whose values must stay character (e.g. Tag IDs with leading zeros)
                .char_keys <- c("WHICH_TAG", "PROB_SPECIES", "DP_FALLBACK_GROWTH_FORMS",
                                "NON_TAPER_CORRECTED_GROWTH_FORMS", "CONFIG_NAME",
                                "INPUT_FILE", "POSTERIOR_SAMPLES_FORMAT", "SPECIES_COL")
                if (tolower(val) %in% c("true", "false")) {
                    val <- as.logical(tolower(val))
                } else if (!(toupper(key) %in% .char_keys) && grepl("^[+-]?[0-9]+$", val)) {
                    val <- as.integer(val)
                } else if (!(toupper(key) %in% .char_keys) && grepl("^[+-]?[0-9]*\\.[0-9]+$", val)) {
                    val <- as.numeric(val)
                }
            } else {
                key <- sub("^--", "", arg)
                # Valueless flags are treated as TRUE
                val <- TRUE
            }
            overrides[[key]] <- val
        } else if (grepl("^-[A-Za-z]$", arg)) {
            # Short single-letter flags (-h handled above); treat others as boolean TRUE
            key <- sub("^-", "", arg)
            overrides[[key]] <- TRUE
        }
    }
    overrides
}

overrides <- parse_args()

############################################################
### 2) Dependencies — package imports
############################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install the 'data.table' package to run this script.")
}
library(data.table)

if (!requireNamespace("here", quietly = TRUE)) {
    stop("Please install the 'here' package to run this script.")
}
library(here)

if (!requireNamespace("withr", quietly = TRUE)) {
    stop("Please install the 'withr' package to run this script.")
}
library(withr)

############################################################
### 3) Defaults — editable run defaults
############################################################

############################################################
### 3.1) Biological parameter estimation settings
############################################################
# Input data and species handling
INPUT_FILE <- here("DATA", "PROCESSED", "6_ViewFullTable_taper_corrected_growth_forms.rds")
FORCE_ONE_SPECIES_PARAMETERS <- FALSE
if (isTRUE(FORCE_ONE_SPECIES_PARAMETERS)) {
    FORCED_SPECIES_LABEL <- "all"
    message("[dp_global main_cpp_chunk_bci.R] FORCE_ONE_SPECIES_PARAMETERS=TRUE: using single species label '", FORCED_SPECIES_LABEL, "' for all trees.")
} else {
    message("[dp_global main_cpp_chunk_bci.R] FORCE_ONE_SPECIES_PARAMETERS=FALSE: using species column from data for parameter estimation.")
}
SPECIES_COL <- NULL

# Biological parameter sources and fixed fallback values.
# _SOURCE controls whether the bound is estimated from data ("data") or fixed ("fixed").
# _FIXED is the fallback value used when _SOURCE = "fixed" or data are too sparse.

# NOTE: DEFINE VERSION OF DP BUNDLE - The Code
dp_bundle_path <- here("2_STEM_IDENTIFICATION", "dpglobal_bundle_full_20260408_183305")

USE_MEASUREMENT_ERROR <- TRUE
MAX_GROWTH_HARD_SOURCE <- "fixed"
MAX_GROWTH_FIXED <- 5
MAX_SHRINK_HARD_SOURCE <- "fixed"
MAX_SHRINK_FIXED <- -0.5
K_SHRINK_SOURCE <- "fixed"
K_SHRINK_FIXED <- 0 # 0 to disable soft penalty
K_GROWTH_SOURCE <- "fixed"
K_GROWTH_FIXED <- 0 # 0 to disable soft penalty
RECRUIT_MAX_SOURCE <- "fixed"
RECRUIT_MAX_FIXED <- (MAX_GROWTH_FIXED * 5) + 0.9999

############################################################
### 3.2) DP solver settings
############################################################
# Controls the DP reconstruction: mode, anchor census, state-space limits,
# slack tracks, growth-form fallbacks, and non-taper-corrected pruning bounds.
DP_MODE <- "marginals+bins" # Options: "none", "marginals", "marginals+bins"
# WHICH_TAG is not used for group selection in the chunked runner (all groups
# are processed), but build_out_dir_name() reads it when RUN_ALL_TAGS=FALSE.
# Keep it defined as 0 (meaning "all") for directory naming purposes.
WHICH_TAG <- "0"
ANCHOR_START_CENSUS <- 7L
DP_VERBOSE <- TRUE
DP_POSTERIOR_TOP_K <- 2L
DP_MAX_TRACKS <- NULL # auto (computed from data)
# DP_MAX_STATES <- 1039L # maximum 30 minutes of runtime per tag check "./2_STEM_IDENTIFICATION/test_complexity_manual.R"
DP_MAX_STATES <- 10000L # maximum 30 minutes of runtime per tag check "./2_STEM_IDENTIFICATION/test_complexity_manual.R"
# DP_MAX_STATES <- 40000L # maximum 30 minutes of runtime per tag check "./2_STEM_IDENTIFICATION/test_complexity_manual.R"
DP_SLACK_TRACKS <- 1L
# NOTE: Optionally require that slack be granted only if an anchor DBH is recruitable
# (i.e., DBH <= Bio_Recruit_MaxDBH_unit + eps). Set FALSE to preserve current behavior.
DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE <- TRUE
# Tolerance (cm) used when comparing anchor DBH to recruit_max_dbh
DP_SLACK_REQUIRE_ANCHOR_EPS <- 1e-6
# Growth forms forcing probabilistic fallback. See main_cpp.R for details.
DP_FALLBACK_GROWTH_FORMS <- character(0)
# Non-taper-corrected growth forms (palms, strangler figs, tree ferns):
# These show real DBH growth plus large apparent variation when HOM changes.
# Wide base prune bounds prevent spurious pruning.
# HOM tolerance adds per-census-pair widening when a HOM column is present.
# Units: cm/year.  Set to NULL to disable the override.
PRUNE_BOUND_FACTOR <- 1.25
NON_TAPER_CORRECTED_GROWTH_FORMS <- c("palm", "strangler_fig", "tree_fern")
NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH <- PRUNE_BOUND_FACTOR * MAX_SHRINK_FIXED
NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH <- PRUNE_BOUND_FACTOR * MAX_GROWTH_FIXED
# HOM tolerance scale: cm of annual DBH tolerance per meter of HOM deviation
# from 1.3 m.  Set 0 to disable HOM widening.
HOM_TOLERANCE_SCALE <- 2.0

# Posterior sampling defaults (disabled by default)
# - POSTERIOR_SAMPLES: number of full-path reconstructions to draw from the DP posterior
# - POSTERIOR_SAMPLES_FORMAT: output format forwarded to DP ('rds','feather','csv')
# - POSTERIOR_SAMPLES_PATH: optional path to write posterior files; when NULL DP writes to out_dir/posteriors
# - POSTERIOR_SAMPLE_SEED: integer seed for reproducible posterior sampling. If NULL, this script defaults to 123L
#   when sampling is enabled; pass --POSTERIOR_SAMPLE_SEED=<int> to override.
POSTERIOR_SAMPLES <- 200L
POSTERIOR_SAMPLES_FORMAT <- "feather" # options: 'rds', 'feather', 'csv'
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- NULL
# Option: allow DP to use a provisional anchor at the last observed DBH census when no TrueStemID exists
ALLOW_PROVISIONAL_DP_ANCHOR <- TRUE

# Number of stochastic samples drawn by the probabilistic greedy matcher
PROB_N_SAMPLES <- 200L

# Species that should bypass DP and go directly to the probabilistic greedy
# matcher. Provide a character vector of Species column values.
PROB_SPECIES <- c("oenoma", "bactma")

# Lookahead weight for probabilistic matcher (0 = disabled, 0.5 = default)
PROB_LOOKAHEAD_WEIGHT <- 0.5

############################################################
### 3.3) Chunking & posterior sampling settings
############################################################
# Controls incremental chunk processing and posterior path sampling.
DP_CHUNK_SIZE <- 15L
DP_CHUNK_RESUME <- TRUE
DP_CHUNK_OVERWRITE <- FALSE
# Optional: limit chunks to a specific range for testing (NULL means all)
DP_CHUNK_START <- NULL
DP_CHUNK_END <- NULL

############################################################
### 3.4) Parallelism settings
############################################################
RUN_ALL_TAGS <- TRUE
MANUAL_CORES <- TRUE # Flag to manually define cores instead of auto-detecting
MANUAL_CORES_VALUE <- 15L # Number of cores to use if MANUAL_CORES=TRUE

############################################################
### 3.5) Output & path settings
############################################################
## create output directory within project
# Base output directory
# default to the project workspace unless overridden.
# Users can set BASE_OUT_DIR via the CLI (e.g. --BASE_OUT_DIR=/some/path) to
# redirect output anywhere (including the home directory) without editing this
# file.
# Pass --BASE_OUT_DIR=/some/path on the CLI to redirect output anywhere
# (e.g. a home-directory folder on a remote machine).
base_out_dir <- here("2_STEM_IDENTIFICATION", "output")
message("[dp_global main_cpp_chunk_bci.R] here root: ", here::here())
message("[dp_global main_cpp_chunk_bci.R] base_out_dir (raw): ", base_out_dir)
base_out_dir <- normalizePath(base_out_dir, winslash = "/", mustWork = FALSE)
message("[dp_global main_cpp_chunk_bci.R] base_out_dir (normalized): ", base_out_dir)

# Optional: explicitly set a subdirectory name for outputs.
# If NULL, an automatic name based on timestamp + key config flags is used.
# OUT_DIR_NAME <- NULL
# CONFIG_NAME is set by the orchestrator (e.g., run_dp_future) to identify the
# experimental configuration; default to NULL so override parsing treats it as
# a valid, known variable rather than an unknown override.
CONFIG_NAME <- NULL
# When non-NULL, OUT_DIR_OVERRIDE bypasses build_out_dir_name() and uses this
# path directly as out_dir. Pass --OUT_DIR_OVERRIDE=/path/to/previous/run on
# the command line to resume into an existing output directory.
OUT_DIR_OVERRIDE <- NULL

# Output path helpers (encode_num, build_out_dir_name) are provided by
# dp_global/R/naming_helpers.R; sourced at the end of this section.

WRITE_DP_CSV <- FALSE
WRITE_DP_RDS <- FALSE
WRITE_DP_FEATHER <- TRUE
WRITE_DP_PDF_PER_CHUNK <- WRITE_DP_PDF <- FALSE
# Set to FALSE when input data have no TrueStemID reference lines to plot.
DP_PDF_INCLUDE_REFERENCE <- FALSE

# Default project root so --PROJECT_ROOT=/path overrides are accepted by the CLI parser
PROJECT_ROOT <- here::here()
# Batch timestamp can be provided by orchestrators; default empty so overrides like --BATCH_TS=... are accepted without warnings
BATCH_TS <- ""
# Naming helpers (encode_num, build_out_dir_name) live in a separate helper
# file to keep the main script concise. Source it early so it's available
# when we compute `out_dir` below.
with_dir(dp_bundle_path, source(file.path("dp_global", "R", "naming_helpers.R")))

############################################################
### 4) CLI reference & override mapping
############################################################
# CLI_REFERENCE maps canonical uppercase CLI flag names to internal variable
# names. Keys are case-insensitive on the command line; use '-' or '_'.
# Orchestrators should keep their flag list in sync with this table.
CLI_REFERENCE <- list(
    INPUT_FILE = "INPUT_FILE",
    FORCE_ONE_SPECIES_PARAMETERS = "FORCE_ONE_SPECIES_PARAMETERS",
    DP_MODE = "DP_MODE",
    ANCHOR_START_CENSUS = "ANCHOR_START_CENSUS",
    DP_VERBOSE = "DP_VERBOSE",
    RUN_ALL_TAGS = "RUN_ALL_TAGS",
    MANUAL_CORES = "MANUAL_CORES",
    MANUAL_CORES_VALUE = "MANUAL_CORES_VALUE",
    WRITE_DP_CSV = "WRITE_DP_CSV",
    WRITE_DP_RDS = "WRITE_DP_RDS",
    WRITE_DP_FEATHER = "WRITE_DP_FEATHER",
    WRITE_DP_PDF = "WRITE_DP_PDF",
    WRITE_DP_PDF_PER_CHUNK = "WRITE_DP_PDF_PER_CHUNK",
    DP_PDF_INCLUDE_REFERENCE = "DP_PDF_INCLUDE_REFERENCE",
    DP_MAX_STATES = "DP_MAX_STATES",
    DP_FALLBACK_GROWTH_FORMS = "DP_FALLBACK_GROWTH_FORMS",
    NON_TAPER_CORRECTED_GROWTH_FORMS = "NON_TAPER_CORRECTED_GROWTH_FORMS",
    NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH = "NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH",
    NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH = "NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH",
    PRUNE_BOUND_FACTOR = "PRUNE_BOUND_FACTOR",
    HOM_TOLERANCE_SCALE = "HOM_TOLERANCE_SCALE",
    DP_SLACK_TRACKS = "DP_SLACK_TRACKS",
    DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE = "DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE",
    DP_SLACK_REQUIRE_ANCHOR_EPS = "DP_SLACK_REQUIRE_ANCHOR_EPS",
    POSTERIOR_SAMPLES = "POSTERIOR_SAMPLES",
    POSTERIOR_SAMPLES_FORMAT = "POSTERIOR_SAMPLES_FORMAT",
    POSTERIOR_SAMPLES_PATH = "POSTERIOR_SAMPLES_PATH",
    POSTERIOR_SAMPLE_SEED = "POSTERIOR_SAMPLE_SEED",
    DP_CHUNK_SIZE = "DP_CHUNK_SIZE",
    DP_CHUNK_RESUME = "DP_CHUNK_RESUME",
    DP_CHUNK_OVERWRITE = "DP_CHUNK_OVERWRITE",
    DP_CHUNK_START = "DP_CHUNK_START",
    DP_CHUNK_END = "DP_CHUNK_END",
    PROJECT_ROOT = "PROJECT_ROOT",
    BATCH_TS = "BATCH_TS",
    CONFIG_NAME = "CONFIG_NAME",
    USE_MEASUREMENT_ERROR = "USE_MEASUREMENT_ERROR",
    ALLOW_PROVISIONAL_DP_ANCHOR = "ALLOW_PROVISIONAL_DP_ANCHOR",
    DP_MAX_TRACKS = "DP_MAX_TRACKS",
    OUT_DIR_OVERRIDE = "OUT_DIR_OVERRIDE",
    PROB_N_SAMPLES = "PROB_N_SAMPLES",
    PROB_SPECIES = "PROB_SPECIES",
    PROB_LOOKAHEAD_WEIGHT = "PROB_LOOKAHEAD_WEIGHT"
)

# Sensitivity and realism flags are not applicable to the chunked runner
# (they require a fully assembled out object). Kept as comments for reference.
# SENSITIVITY_MODE <- "none"
# RUN_REALISM_REPORT <- FALSE

############################################################
### 4.1) Help & override application
############################################################
print_help <- function() {
    cat("Usage: Rscript scripts/main_cpp_chunk_bci.R [--KEY=VALUE] [--FLAG]\n")
    cat("Common keys and defaults:\n")

    # Prefer printing the canonical CLI keys from CLI_REFERENCE so the help is
    # always in sync with the internal mapping.
    for (key in names(CLI_REFERENCE)) {
        varname <- CLI_REFERENCE[[key]]
        if (exists(varname, envir = globalenv())) {
            val <- get(varname, envir = globalenv())
        } else {
            val <- "<not set>"
        }
        cat(sprintf("  --%s = %s\n", key, as.character(val)))
    }

    cat("\nFlags without =value are treated as boolean TRUE (e.g., --DRY_RUN).\n")
}

# If user asked for help, print and exit (do this before applying overrides)
if (isTRUE(overrides$help) || isTRUE(overrides$h)) {
    print_help()
    quit(save = "no", status = 0)
}

# Apply command-line overrides with validation and warnings for unknown keys
# This logic is intentionally flexible: users can pass keys such as
#   --POSTERIOR-SAMPLES=200  or  --posterior-samples=200  or  --POSTERIOR_SAMPLES=200
# All of these will match the canonical in-memory variable names such as
# `POSTERIOR_SAMPLES` or `POSTERIOR_SAMPLES_FORMAT` or `which_tag` regardless of
# case or punctuation.

# Helper: normalize an override key to uppercase letters and underscores
normalize_key <- function(k) {
    toupper(gsub("[- ]", "_", k))
}

# Helper: find an existing variable name that matches the normalized key (case-insensitive)
find_matching_var <- function(norm_key) {
    # Normalize current environment variable names for comparison (case-insensitive)
    norm_var <- function(v) toupper(gsub("[^A-Z0-9_]", "", toupper(v)))
    vars <- ls(envir = globalenv())
    matches <- sapply(vars, function(v) norm_var(v) == norm_key)
    if (any(matches)) {
        return(vars[which(matches)[1]])
    }
    NULL
}

for (name in names(overrides)) {
    if (name %in% c("help", "h")) next

    norm <- normalize_key(name)
    match_var <- find_matching_var(norm)

    if (is.null(match_var)) {
        warning(sprintf("[dp_global main_cpp_chunk_bci.R] Unknown override '%s' (ignored).\n", name))
        next
    }

    # Coerce type where appropriate (basic numeric/logical conversions were
    # already attempted in parse_args). Here we ensure that integer-like
    # numeric values are stored as integers when the default was integer.
    old_val <- get(match_var, inherits = FALSE)
    new_val <- overrides[[name]]

    if (is.integer(old_val) && is.numeric(new_val)) {
        new_val <- as.integer(new_val)
    }

    assign(match_var, new_val, envir = globalenv())
    message("[dp_global main_cpp_chunk_bci.R] Overriding ", match_var, " = ", as.character(new_val))
}

# Backwards-compatibility aliases removed. Use canonical ALL-CAPS variables (e.g., WHICH_TAG, INPUT_FILE) everywhere; update scripts that relied on lowercase globals.

# Post-override validation: Check a few key options for allowed values and types
if (!DP_MODE %in% c("none", "marginals", "marginals+bins", "map")) {
    stop("Invalid DP_MODE: ", DP_MODE, ". Allowed: 'none','marginals','marginals+bins','map'.")
}
if (!POSTERIOR_SAMPLES_FORMAT %in% c("rds", "feather", "csv")) {
    stop("Invalid POSTERIOR_SAMPLES_FORMAT: ", POSTERIOR_SAMPLES_FORMAT, ". Allowed: 'rds','feather','csv'.")
}
if (!is.null(POSTERIOR_SAMPLES) && (!is.numeric(POSTERIOR_SAMPLES) || as.integer(POSTERIOR_SAMPLES) < 0L)) {
    stop("POSTERIOR_SAMPLES must be a non-negative integer or 0 to disable.")
}
if (!is.null(POSTERIOR_SAMPLE_SEED) && (!is.numeric(POSTERIOR_SAMPLE_SEED) || as.integer(POSTERIOR_SAMPLE_SEED) < 0L)) {
    stop("POSTERIOR_SAMPLE_SEED must be a non-negative integer or NULL.")
}

# Recompute derived values that depend on overridable inputs so CLI overrides take effect
# - MC_CORES depends on MANUAL_CORES and MANUAL_CORES_VALUE
MC_CORES <- if (exists("MANUAL_CORES") && isTRUE(MANUAL_CORES)) {
    as.integer(MANUAL_CORES_VALUE)
} else if (requireNamespace("parallel", quietly = TRUE)) {
    max(1L, parallel::detectCores(logical = TRUE) - 1L)
} else {
    1L
}

# Optional override: explicitly set base output directory directly.
# If --BASE_OUT_DIR=/path was supplied, the override loop already set
# base_out_dir; re-normalize here to handle tilde expansion.
# Usage: --BASE_OUT_DIR=/absolute/path/to/output
if ("BASE_OUT_DIR" %in% toupper(gsub("[- ]", "_", names(overrides)))) {
    base_out_dir <- normalizePath(base_out_dir, winslash = "/", mustWork = FALSE)
    message("[dp_global main_cpp_chunk_bci.R] BASE_OUT_DIR override applied: ", base_out_dir)
}

# Optional override: explicitly set project root (useful when running under different working dirs or job launchers)
# Usage: --PROJECT_ROOT=/absolute/path/to/project
if ("PROJECT_ROOT" %in% toupper(gsub("[- ]", "_", names(overrides)))) {
    message("[dp_global main_cpp_chunk_bci.R] Using PROJECT_ROOT override: ", PROJECT_ROOT)
    base_out_dir <- normalizePath(file.path(PROJECT_ROOT, "2_STEM_IDENTIFICATION", "output"), winslash = "/", mustWork = FALSE)
    message("[dp_global main_cpp_chunk_bci.R] base_out_dir overridden to: ", base_out_dir)
}

############################################################
### 5) Input validation
############################################################
# Validate INPUT_FILE, derive boolean flags from mode strings, compute out_dir,
# and define filesystem/logging helpers used throughout the script.

if (!exists("INPUT_FILE") || is.null(INPUT_FILE) || !nzchar(INPUT_FILE)) {
    stop(
        "No INPUT_FILE specified. ",
        "Provide --INPUT_FILE=... on the command line or define INPUT_FILE explicitly."
    )
}

if (!file.exists(INPUT_FILE)) {
    stop("INPUT_FILE does not exist: ", INPUT_FILE)
}

# Derive booleans from modes
RUN_DP <- DP_MODE != "none"
ADD_DP_POSTERIOR_BINS <- DP_MODE == "marginals+bins"
# Not used by chunked runner; leave commented out to avoid confusion
# RUN_SENSITIVITY <- SENSITIVITY_MODE != "none"
# WRITE_OUTPUTS <- SENSITIVITY_MODE %in% c("run+write", "run+write+pdf")
# MAKE_ALL_SWEEPS_PDF <- SENSITIVITY_MODE == "run+write+pdf"

# Final output directory for this run (created at runtime in run_main_chunked()).
# OUT_DIR_OVERRIDE takes precedence: set it to resume into an existing output
# directory (--OUT_DIR_OVERRIDE=/path/to/previous/run).
if (!is.null(OUT_DIR_OVERRIDE) && nzchar(as.character(OUT_DIR_OVERRIDE))) {
    out_dir <- normalizePath(OUT_DIR_OVERRIDE, winslash = "/", mustWork = FALSE)
    message("[dp_global main_cpp_chunk_bci.R] OUT_DIR_OVERRIDE set — using existing dir: ", out_dir)
} else {
    out_dir <- file.path(base_out_dir, build_out_dir_name())
    message("[dp_global main_cpp_chunk_bci.R] out_dir (computed): ", out_dir)
}
message("[dp_global main_cpp_chunk_bci.R] getwd(): ", getwd())

# Define filesystem and logging helpers here (before source() calls) so they
# are available in the input validation block above and in section 6.
ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        message("[dp_global main_cpp_chunk_bci.R] Created directory: ", path)
    } else {
        message("[dp_global main_cpp_chunk_bci.R] Directory already exists: ", path)
    }
    invisible(path)
}

log_msg <- function(msg, level = "INFO") {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full <- sprintf("[%s] %s: %s", ts, level, msg)
    message(full)
    if (exists("out_dir") && nzchar(out_dir)) {
        tryCatch(write(full, file = file.path(out_dir, "run_log.txt"), append = TRUE), error = function(e) NULL)
    }
    invisible(full)
}

# Centralized DP naming and path helpers (same helpers as main_cpp.R)
DP_BASE <- "stem_reconstruction_dp_global_rcpp"
make_out_path <- function(base = DP_BASE, ext = "csv", dir = out_dir) {
    file.path(dir, paste0(base, ".", ext))
}
out_path <- function(target, ext = NULL) {
    if (!is.null(ext) && nzchar(ext)) {
        ext_use <- ext
    } else {
        ext_use <- "csv"
    }
    switch(target,
        dp = make_out_path(DP_BASE, ext_use),
        dp_csv = make_out_path(DP_BASE, "csv"),
        dp_rds = make_out_path(DP_BASE, "rds"),
        dp_feather = make_out_path(DP_BASE, "feather"),
        dp_pdf = make_out_path(DP_BASE, "pdf"),
        posteriors = file.path(out_dir, "posteriors"),
        stop("Unknown out_path target: ", target)
    )
}

maybe_write <- function(flag, path, write_expr, msg = NULL) {
    if (!isTRUE(flag)) {
        return(invisible(FALSE))
    }
    ensure_dir(dirname(path))
    tryCatch(
        {
            write_expr()
            log_msg(sprintf("Wrote %s: %s", if (!is.null(msg)) msg else basename(path), path))
            TRUE
        },
        error = function(e) {
            log_msg(sprintf("Failed to write %s: %s", path, conditionMessage(e)), "ERROR")
            FALSE
        }
    )
}

DP_CSV_FILE <- out_path("dp_csv")
DP_RDS_FILE <- out_path("dp_rds")
DP_FEATHER_FILE <- out_path("dp_feather")
DP_PDF_FILE <- out_path("dp_pdf")
# Note: chunk-level per-file outputs (per-chunk RDS/Feather/PDF) still use chunk-specific names; these top-level globals point to the combined run-level names.

# If posterior sampling is requested, provide the DP with the run's base out_dir.
# The DP itself will create the 'posteriors/' subdirectory (avoids double-nesting).
# These defaults can still be overridden via CLI (e.g., --POSTERIOR_SAMPLES_PATH=... or --POSTERIOR_SAMPLE_SEED=...)
if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L) {
    if (is.null(POSTERIOR_SAMPLES_PATH) || !nzchar(POSTERIOR_SAMPLES_PATH)) {
        # Provide the DP with the run's base out_dir; the DP will create the
        # 'posteriors/' subdirectory itself, avoiding double-nesting.
        POSTERIOR_SAMPLES_PATH <- out_dir
    }
    # If user supplied a path that ends in 'posteriors', remove that suffix so
    # DP's internal creation of '<base>/posteriors' doesn't create nested folders.
    if (basename(POSTERIOR_SAMPLES_PATH) == "posteriors") {
        POSTERIOR_SAMPLES_PATH <- dirname(POSTERIOR_SAMPLES_PATH)
    }
    # normalize path for consistency; don't require it to exist yet (created in run_main)
    POSTERIOR_SAMPLES_PATH <- normalizePath(POSTERIOR_SAMPLES_PATH, winslash = "/", mustWork = FALSE)
    if (is.null(POSTERIOR_SAMPLE_SEED)) {
        POSTERIOR_SAMPLE_SEED <- as.integer(123L)
    } else {
        POSTERIOR_SAMPLE_SEED <- as.integer(POSTERIOR_SAMPLE_SEED)
    }
} else {
    # Ensure disabled sampling leaves a NULL path to avoid accidental writes
    POSTERIOR_SAMPLES_PATH <- NULL
}

# Coerce POSTERIOR_SAMPLE_SEED whenever explicitly provided via CLI, even when
# posterior sampling is disabled (POSTERIOR_SAMPLES=0).  The seed is also used
# by the probabilistic matcher's Gumbel-noise draws, so honouring it here
# makes fallback runs reproducible.
if (!is.null(POSTERIOR_SAMPLE_SEED)) {
    POSTERIOR_SAMPLE_SEED <- as.integer(POSTERIOR_SAMPLE_SEED)
}

############################################################
### 6) Source project code
############################################################
# Load dp_global R modules: DP solver, biological parameter estimation,
# sensitivity and realism helpers, naming utilities.
with_dir(dp_bundle_path, source(file.path("dp_global", "R", "dp_global_main.R")))
# NOTE: sensitivity_transition_cost_bio.R, realism_calibration.R, and
# k_tuning_viz.R are NOT sourced here — they require a fully assembled
# output object and are not applicable to the chunked runner.
# naming_helpers.R is already sourced in section 3.5.

############################################################
### 7) Helpers — data-manipulation utilities
############################################################
# These functions support parameter estimation and input preparation.
# Filesystem/logging helpers (ensure_dir, log_msg, maybe_write) are defined
# in section 5 so they are available before project code is sourced.

# Unit reminder for soft-penalty tuning:
# - Soft penalties operate on DBH differences (cm) over the census interval.
# - Hard guardrails operate on annualized growth (cm/year).
# - Soft penalty is quadratic: soft_cost = k * delta_cm^2
#   To have delta_cm = D contribute cost C, set k = C / D^2.
ensure_species_column <- function(x) {
    if (isTRUE(FORCE_ONE_SPECIES_PARAMETERS)) {
        x[, species := FORCED_SPECIES_LABEL]
        return(x)
    }

    if ("species" %in% names(x)) {
        x[, species := as.character(species)]
        return(x)
    }

    species_col <- SPECIES_COL
    if (is.null(species_col) || !nzchar(species_col) || !(species_col %in% names(x))) {
        candidates <- c(
            "Species", "species", "SP", "sp", "spcode", "sp_code",
            "taxon", "Taxon", "spname", "Sp"
        )
        found <- candidates[candidates %in% names(x)]
        if (length(found) > 0L) {
            species_col <- found[[1L]]
            message("[dp_global main_cpp_chunk_bci.R] Using '", species_col, "' as species column. Set SPECIES_COL to override.")
        }
    }

    if (is.null(species_col) || !nzchar(species_col) || !(species_col %in% names(x))) {
        stop(
            "FORCE_ONE_SPECIES_PARAMETERS=FALSE but no species column found. ",
            "Add a 'species' column to the input, or set SPECIES_COL to an existing column name."
        )
    }

    x[, species := as.character(get(species_col))]
    x
}

get_growth_mu_const <- function(growth_list) {
    # Growth mean model: mu(DBH) = alpha + gamma*log(DBH)
    # Backward-compat: if alpha/gamma absent, use constant mean mu and gamma=0.
    if (!is.null(growth_list$alpha)) {
        return(growth_list$alpha)
    }
    growth_list$mu
}

get_nested_numeric <- function(x, expr, fallback = NULL) {
    v <- tryCatch(eval(expr, envir = x), error = function(e) NULL)
    if (!is.null(v) && length(v) == 1L && is.finite(v)) {
        return(as.numeric(v))
    }
    fallback
}

attach_bio_columns <- function(xrun, bio_pars) {
    xrun[, `:=`(
        Bio_Mu_Growth = get_growth_mu_const(bio_pars[[.BY$species]]$growth),
        Bio_Gamma_Growth = {
            g <- bio_pars[[.BY$species]]$growth
            if (!is.null(g$gamma)) g$gamma else 0
        },
        Bio_Sigma0_Growth = bio_pars[[.BY$species]]$growth$sigma0,
        Bio_Sigma1_Growth = bio_pars[[.BY$species]]$growth$sigma1,
        Bio_H0_Mortality = bio_pars[[.BY$species]]$mortality$h0,
        Bio_Beta_Mortality = bio_pars[[.BY$species]]$mortality$beta,
        Bio_Recruit_Meanlog = bio_pars[[.BY$species]]$recruitment$meanlog,
        Bio_Recruit_Sdlog = bio_pars[[.BY$species]]$recruitment$sdlog,
        Bio_Recruit_MaxDBH_unit = bio_pars[[.BY$species]]$recruitment$recruit_max_dbh,
        Bio_Recruitment_lambda = bio_pars[[.BY$species]]$recruitment$lambda,
        Bio_Max_Shrink = {
            sh <- bio_pars[[.BY$species]]$shrinkage
            get_nested_numeric(sh, quote(guardrails$hard$value), fallback = sh$max_shrink)
        },
        Bio_K_Shrink = {
            sh <- bio_pars[[.BY$species]]$shrinkage
            get_nested_numeric(sh, quote(penalties$soft$k), fallback = sh$k_shrink)
        },
        Bio_Max_Growth = {
            g <- bio_pars[[.BY$species]]$growth
            get_nested_numeric(g, quote(guardrails$hard$value), fallback = g$max_growth)
        },
        Bio_Max_Growth_Soft = {
            g <- bio_pars[[.BY$species]]$growth
            get_nested_numeric(g, quote(guardrails$soft$value), fallback = g$max_growth_soft)
        },
        Bio_K_Growth = {
            g <- bio_pars[[.BY$species]]$growth
            get_nested_numeric(g, quote(penalties$soft$k), fallback = g$k_growth)
        }
    ), by = species]
    xrun
}

auto_dp_max_tracks <- function(xrun) {
    max_obs_any_tag_census <- xrun[
        CensusID <= ANCHOR_START_CENSUS & !is.na(DBH),
        .N,
        by = .(Tag, CensusID)
    ][, max(N, na.rm = TRUE)]
    if (!is.finite(max_obs_any_tag_census)) max_obs_any_tag_census <- 0L
    as.integer(max_obs_any_tag_census + 1L)
}

############################################################
### 8) Core DP functions
############################################################
# run_dp_one_group()       — runs the DP solver for a single (Tag, species) group
# maybe_add_posterior_bins() — applies posterior binning when DP_MODE includes bins

run_dp_one_group <- function(dtg, dp_max_tracks, chunk_id = NULL) {
    tag_label <- tryCatch(as.character(unique(dtg$Tag)[1]), error = function(e) "<unknown>")

    out <- tryCatch(
    match_stems_dp_global_backward_marginals_batch(
        tree_data = data.table::copy(dtg),
        chunk_id = chunk_id,
        min_growth = MAX_SHRINK_FIXED,
        max_growth = MAX_GROWTH_FIXED,
        anchor_start = ANCHOR_START_CENSUS,
        max_tracks = dp_max_tracks,
        max_states = DP_MAX_STATES,
        slack_tracks = DP_SLACK_TRACKS,
        slack_require_anchor_recruitable = DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE,
        slack_require_anchor_eps = DP_SLACK_REQUIRE_ANCHOR_EPS,
        temperature = 1,
        posterior_top_k = DP_POSTERIOR_TOP_K,
        fallback_growth_forms = DP_FALLBACK_GROWTH_FORMS,
        non_taper_corrected_growth_forms = NON_TAPER_CORRECTED_GROWTH_FORMS,
        non_taper_corrected_prune_min_growth = NON_TAPER_CORRECTED_PRUNE_MIN_GROWTH,
        non_taper_corrected_prune_max_growth = NON_TAPER_CORRECTED_PRUNE_MAX_GROWTH,
        hom_tolerance_scale = HOM_TOLERANCE_SCALE,
        # posterior sampling controls (disabled by default)
        posterior_samples = POSTERIOR_SAMPLES,
        posterior_samples_format = POSTERIOR_SAMPLES_FORMAT,
        posterior_samples_path = POSTERIOR_SAMPLES_PATH,
        posterior_sample_seed = POSTERIOR_SAMPLE_SEED,
        use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
        # prune controls
        # NOTE: You can always define very wide based on the parameter data you have.
        prune_hard = TRUE,
        prune_min_growth = MAX_SHRINK_FIXED * PRUNE_BOUND_FACTOR, # very wide fixed bounds
        prune_max_growth = MAX_GROWTH_FIXED * PRUNE_BOUND_FACTOR, # very wide fixed bounds
        prune_use_bio_bounds = FALSE, # use fixed prune bounds instead of biological ones
        prune_recruit_max_dbh = RECRUIT_MAX_FIXED * PRUNE_BOUND_FACTOR, # very high recruit max dbh
        prune_use_bio_recruit = FALSE, # FALSE = use prune_recruit_max_dbh instead of biological (and margin) one, TRUE, set prune_recruit_max_dbh as min(prune_recruit_max_dbh, bio_recruit_max_dbh * 1.2)
        allow_provisional_anchor = isTRUE(ALLOW_PROVISIONAL_DP_ANCHOR),
        verbose = isTRUE(DP_VERBOSE),
        prob_n_samples = PROB_N_SAMPLES,
        prob_species = PROB_SPECIES,
        prob_lookahead_weight = PROB_LOOKAHEAD_WEIGHT
        ),
        error = function(e) {
            tag_val <- tryCatch(unique(dtg$Tag)[1], error = function(e2) NA)
            msg <- conditionMessage(e)
            log_msg(sprintf("DP error for Tag=%s: %s — falling back to probabilistic", tag_val, msg), "WARN")
            # Probabilistic fallback so the tag is not lost
            out <- match_stems_probabilistic(
                tree_data = data.table::copy(dtg),
                min_growth = MAX_SHRINK_FIXED,
                max_growth = MAX_GROWTH_FIXED,
                anchor_start = ANCHOR_START_CENSUS,
                n_samples = PROB_N_SAMPLES,
                temperature = 1,
                posterior_top_k = DP_POSTERIOR_TOP_K,
                posterior_samples_path = POSTERIOR_SAMPLES_PATH,
                posterior_samples_format = POSTERIOR_SAMPLES_FORMAT,
                prune_min_growth = MAX_SHRINK_FIXED * PRUNE_BOUND_FACTOR,
                prune_max_growth = MAX_GROWTH_FIXED * PRUNE_BOUND_FACTOR,
                prune_recruit_max_dbh = RECRUIT_MAX_FIXED * PRUNE_BOUND_FACTOR,
                prob_lookahead_weight = PROB_LOOKAHEAD_WEIGHT,
                verbose = isTRUE(DP_VERBOSE)
            )
            if (!("DP_FallbackReason" %in% names(out))) out[, DP_FallbackReason := NA_character_]
            out[, DP_FallbackReason := paste0("error:", substr(msg, 1, 200))]
            # R-boundary splitting: sever tracks that cross live R-coded censuses
            .r_regex_eh <- "\\b(R|RP|RF|RT|QR|OR)\\b"
            if ("ListOfTSM" %in% names(out)) {
                .pre_cc_eh <- sort(unique(out$CensusID[out$CensusID <= ANCHOR_START_CENSUS]))
                for (.cc_eh in .pre_cc_eh) {
                    .lr_eh <- which(out$CensusID == .cc_eh & !is.na(out$DBH))
                    if (length(.lr_eh) == 0L) next
                    .tsm_eh <- out$ListOfTSM[.lr_eh]
                    if (!any(!is.na(.tsm_eh) & grepl(.r_regex_eh, .tsm_eh, perl = TRUE))) next
                    .before_eh <- .pre_cc_eh[.pre_cc_eh < .cc_eh]
                    if (length(.before_eh) == 0L) next
                    .ids_bef <- unique(out$ReconstructedStemID[out$CensusID %in% .before_eh & !is.na(out$ReconstructedStemID)])
                    .ids_aft <- unique(out$ReconstructedStemID[out$CensusID >= .cc_eh & !is.na(out$ReconstructedStemID)])
                    .cross <- intersect(.ids_bef, .ids_aft)
                    .mx_eh <- suppressWarnings(max(out$ReconstructedStemID, na.rm = TRUE))
                    if (!is.finite(.mx_eh)) .mx_eh <- 0L
                    for (.old_eh in .cross) {
                        .mx_eh <- .mx_eh + 1L
                        out[CensusID %in% .before_eh & ReconstructedStemID == .old_eh, ReconstructedStemID := as.integer(.mx_eh)]
                    }
                }
            }
            out
        }
    )

    # --- TrueStemID preservation assertion --------------------------------
    # After both DP and probabilistic pathways, verify that every row with
    # a known TrueStemID has ReconstructedStemID == TrueStemID.
    if ("TrueStemID" %in% names(out) && "ReconstructedStemID" %in% names(out)) {
        .given <- out[!is.na(TrueStemID)]
        if (nrow(.given) > 0L) {
            .bad <- .given[as.integer(TrueStemID) != as.integer(ReconstructedStemID)]
            if (nrow(.bad) > 0L) {
                log_msg(sprintf(
                    "WARNING Tag=%s: %d row(s) with TrueStemID != ReconstructedStemID after reconstruction (censuses: %s)",
                    tag_label, nrow(.bad),
                    paste(unique(.bad$CensusID), collapse = ",")), "WARN")
            }
        }
    }

    out
}

maybe_add_posterior_bins <- function(out) {
    if (isTRUE(ADD_DP_POSTERIOR_BINS) && !is.null(out)) {
        out <- add_dp_posterior_bins(
            out,
            confident_prob = 0.95,
            unlinked_prob = 0.5,
            use_reconstructed_prob = TRUE,
            out_col = "DP_PosteriorBin"
        )
    }
    out
}

############################################################
### 9) Main pipeline — run_main_chunked()
############################################################
# Loads input data, estimates biological parameters, processes groups in
# chunks of DP_CHUNK_SIZE, writes incremental CSV/RDS/PDF outputs per chunk,
# and records run markers (run_started.txt, run_finished.txt, run_log.txt).
run_main_chunked <- function() {
    ensure_dir(out_dir)
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_started.txt"))
        },
        error = function(e) {
            message("[dp_global main_cpp_chunk_bci.R] Warning writing run_started marker: ", conditionMessage(e))
        }
    )
    log_msg("Started run")

    # Create posteriors subdirectory (DP writes its files into <base>/posteriors)
    if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && !is.null(POSTERIOR_SAMPLES_PATH) && nzchar(POSTERIOR_SAMPLES_PATH)) {
        ensure_dir(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"))
        log_msg(paste("Ensured posterior samples path:", file.path(POSTERIOR_SAMPLES_PATH, "posteriors")))
    }

    # 5.1 Load data
    cat("[dp_global main_cpp_chunk_bci.R] Loading input data from: ", INPUT_FILE, "\n")
    xraw <- as.data.table(readRDS(INPUT_FILE))

    xraw[, StemID := as.integer(as.character(StemID))] # ensure StemID is integer (handles cases where it might be read as numeric or factor)

    ## ---- Clean raw data ------------------------------------------------------
    xraw[
        is.na(Mnemonic) | is.infinite(Mnemonic),
        Mnemonic := "unknown"
    ][, Tag := as.character(Tag)]

    ## ---- Anchor data ---------------------------------------------------------
    anchor_data <- xraw[
        CensusID >= ANCHOR_START_CENSUS &
            !is.na(dbh_with_best_candidate_taper_corrected),
        .(
            Tag,
            CensusID,
            TrueStemID = StemID,
            Mnemonic,
            # NOTE: BCI data is by default in mm
            DBH = dbh_with_best_candidate_taper_corrected / 10, # mm → cm
            ExactDate,
            growth_form
        )
    ]

    ## ---- Mean dates ----------------------------------------------------------
    mean_date_tag_census <- anchor_data[
        , .(MeanExactDate = as.IDate(mean(as.numeric(ExactDate)))),
        keyby = .(Tag, CensusID)
    ]

    mean_date_census <- anchor_data[
        , .(MeanExactDate_Census = as.IDate(mean(as.numeric(ExactDate)))),
        keyby = CensusID
    ]

    ## ---- Complete grid -------------------------------------------------------
    # Get unique Tag-TrueStemID combinations (preserves which stems belong to which tags)
    tag_stem_combos <- unique(anchor_data[, .(Tag, TrueStemID)])

    full_anchor_data <- CJ(
        Tag_TrueStemID = tag_stem_combos[, paste(Tag, TrueStemID, sep = "_")],
        CensusID = unique(anchor_data[CensusID >= ANCHOR_START_CENSUS, CensusID]),
        sorted = FALSE
    )

    # Split IDs
    full_anchor_data[, c("Tag", "TrueStemID") := tstrsplit(Tag_TrueStemID, "_", fixed = TRUE)]
    full_anchor_data[, Tag_TrueStemID := NULL]
    full_anchor_data[, TrueStemID := as.integer(as.character(TrueStemID))]

    ## ---- Join dates ----------------------------------------------------------
    full_anchor_data[mean_date_tag_census, MeanExactDate := i.MeanExactDate,
        on = .(Tag, CensusID)
    ]
    full_anchor_data[mean_date_census, MeanExactDate_Census := i.MeanExactDate_Census,
        on = .(CensusID)
    ]

    # Fill missing with census means
    full_anchor_data[is.na(MeanExactDate), MeanExactDate := MeanExactDate_Census]

    ## ---- Join anchor values --------------------------------------------------
    full_anchor_data[
        anchor_data,
        `:=`(DBH = i.DBH, Mnemonic = i.Mnemonic),
        on = .(Tag, TrueStemID, CensusID)
    ]

    ## ---- Fill species from Tag lookup ----------------------------------------
    tag_mnemonic <- anchor_data[!is.na(Mnemonic), .(Mnemonic = first(Mnemonic)), keyby = .(Tag, growth_form)]
    full_anchor_data[tag_mnemonic, Mnemonic := i.Mnemonic, on = "Tag"]
    # add growth form to full_anchor_data
    full_anchor_data[tag_mnemonic, growth_form := i.growth_form, on = "Tag"]

    ## ---- Finalize ------------------------------------------------------------
    full_anchor_data[
        , `:=`(
            ExactDate = as.IDate(MeanExactDate),
            MeanExactDate = NULL,
            MeanExactDate_Census = NULL
        )
    ]

    setorder(full_anchor_data, Tag, TrueStemID, CensusID)

    # # ---- Check all tags x censusid are complete ----
    # # Compute expected vs actual counts per Tag and identify any gaps (tags with missing censuses)
    # # using the taper-corrected table. Add missing Tag × CensusID rows below if needed.
    # xraw_unique <- unique(full_anchor_data[, .(Tag, CensusID)])
    # # Get the range per tag
    # tag_ranges <- xraw_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]
    # # Add expected count (how many censuses should exist)
    # tag_ranges[, expected_count := max_c - min_c + 1L]
    # # Get actual count per tag
    # actual_counts <- xraw_unique[, .(actual_count = .N), by = Tag]
    # # Merge and compare
    # tag_check <- tag_ranges[actual_counts, on = "Tag"]
    # tag_check[, complete := actual_count == expected_count]
    # # Summary
    # tag_check[, .N, by = complete]

    # # See tags with missing censuses
    # missing_tags <- tag_check[complete == FALSE]
    # missing_tags[, gap := expected_count - actual_count]

    # # Summary stats
    # cat("Total tags:", nrow(tag_check), "\n")
    # cat("Complete tags:", tag_check[complete == TRUE, .N], "\n")
    # cat("Tags with gaps:", tag_check[complete == FALSE, .N], "\n")
    # if (nrow(missing_tags) > 0) {
    #   cat("Total missing observations:", sum(missing_tags$gap), "\n")
    # }

    ## ---- Select species with sufficient repeated DBH -------------------------
    to_do <- copy(full_anchor_data)

    setorder(to_do, Tag, TrueStemID, CensusID)

    # previous DBH within Tag x Stem
    to_do[
        , shift_DBH := shift(DBH),
        by = .(Tag, TrueStemID)
    ]

    # valid consecutive DBH pairs
    to_do[
        , valid_growth_pair := !is.na(shift_DBH) & !is.na(DBH)
    ]

    # per-species summaries (single grouped pass)
    to_do[
        , `:=`(
            n_valid_growth_pair = sum(valid_growth_pair),
            n_tags = uniqueN(Tag)
        ),
        by = Mnemonic
    ]

    ## ---- Filter: require minimum stems per log-DBH bin ----------------------
    ## Species must have at least `min_per_bin` valid-pair stems in each of
    ## `n_bins` equal-width bins on the log(DBH) scale.  This ensures that
    ## growth-mean, variance, and guardrail regressions are anchored across
    ## the full size range, preventing extrapolation artefacts for large stems.
    ##
    ## Species that fail are excluded from per-species estimation and fall
    ## back to the pooled "all_tree_shrub" parameters.
    has_dbh_coverage <- function(dt, n_bins = 4L, min_per_bin = 3L) {
        # Step 1: Extract valid DBH values
        d <- dt[!is.na(DBH) & is.finite(DBH) & DBH > 0, DBH]
        # Step 2: Early exit if insufficient data
        if (length(d) < n_bins * min_per_bin) {
            return(FALSE)
        }
        # Step 3: Log-transform (because tree diameters are log-normally distributed)
        log_d <- log(d)
        # Step 4: Create bin edges
        breaks <- seq(min(log_d), max(log_d), length.out = n_bins + 1L)
        ## Widen edges slightly so min/max fall inside
        breaks[1] <- breaks[1] - 1e-6 # Ensure min value is included
        breaks[n_bins + 1] <- breaks[n_bins + 1] + 1e-6 # Ensure max value is included
        # Step 5: Assign each observation to a bin
        bin_id <- findInterval(log_d, breaks, rightmost.closed = TRUE)
        # Step 6: Count observations per bin - Every bin must meet the minimum count
        tab <- tabulate(bin_id, nbins = n_bins)
        # Step 7: Check if all bins meet minimum threshold
        all(tab >= min_per_bin)
    }

    ## ---- Helper: prepare data for bio parameter estimation -------------------
    # Returns a filtered data.table ready for run_bio_par_estimation().
    # - growth_forms:  character vector of growth_form values to include
    # - species_label: NULL   -> per-species mode (species := Mnemonic, coverage by Mnemonic)
    #                  string -> pooled mode (species := species_label, coverage on whole subset)
    # - min_tags:  minimum n_tags required; NULL to skip this filter
    # - min_pairs: minimum n_valid_growth_pair required
    # - n_bins, min_per_bin: forwarded to has_dbh_coverage()
    prepare_bio_data <- function(dt, growth_forms, species_label = NULL,
                                 min_tags = 30L, min_pairs = 10L,
                                 n_bins = 4L, min_per_bin = 4L) {
        out <- dt[growth_form %in% growth_forms]
        if (is.null(species_label)) {
            # Per-species mode: compute coverage by Mnemonic, filter, then label
            out[, has_size_coverage := has_dbh_coverage(.SD, n_bins = n_bins, min_per_bin = min_per_bin),
                by = Mnemonic
            ]
            keep <- out$n_valid_growth_pair >= min_pairs & out$has_size_coverage == TRUE
            if (!is.null(min_tags)) keep <- keep & out$n_tags >= min_tags
            out <- out[keep][, species := Mnemonic][]
        } else {
            # Pooled mode: assign label, compute coverage on whole subset, filter
            out[, species := species_label]
            out[, has_size_coverage := has_dbh_coverage(.SD, n_bins = n_bins, min_per_bin = min_per_bin),
                by = species
            ]
            keep <- out$n_valid_growth_pair >= min_pairs & out$has_size_coverage == TRUE
            if (!is.null(min_tags)) keep <- keep & out$n_tags >= min_tags
            out <- out[keep][]
        }
        out
    }

    ## ---- Diagnostic: sweep parameter grid for species coverage ----------------
    # For each growth_form group, sweep over candidate values of min_tags,
    # min_pairs, n_bins, and min_per_bin.  Returns a data.table with one row
    # per parameter combination showing how many (and which fraction of)
    # species pass all filters.
    sweep_bio_thresholds <- function(dt,
                                     growth_forms,
                                     min_tags_vals = c(5L, 10L, 15L, 20L, 30L, 50L),
                                     min_pairs_vals = c(3L, 5L, 10L, 15L, 20L),
                                     n_bins_vals = c(3L, 4L, 5L),
                                     min_per_bin_vals = c(3L, 4L, 5L, 8L)) {
        sub <- dt[growth_form %in% growth_forms]
        if (nrow(sub) == 0L) {
            warning("No rows match growth_forms: ", paste(growth_forms, collapse = ", "))
            return(data.table())
        }

        mnemonics <- unique(sub$Mnemonic)
        n_total <- length(mnemonics)

        # Pre-compute per-species scalars (these don't depend on n_bins / min_per_bin)
        sp_stats <- sub[, .(
            n_tags              = first(n_tags),
            n_valid_growth_pair = first(n_valid_growth_pair)
        ), by = Mnemonic]

        # Pre-compute has_dbh_coverage for every (species, n_bins, min_per_bin) combo
        coverage_grid <- CJ(
            Mnemonic = mnemonics,
            n_bins = n_bins_vals,
            min_per_bin = min_per_bin_vals,
            sorted = FALSE
        )

        coverage_grid[, passes_coverage := {
            has_dbh_coverage(sub[Mnemonic == .BY$Mnemonic],
                n_bins = .BY$n_bins,
                min_per_bin = .BY$min_per_bin
            )
        }, by = .(Mnemonic, n_bins, min_per_bin)]

        # Full parameter grid
        param_grid <- CJ(
            min_tags = min_tags_vals,
            min_pairs = min_pairs_vals,
            n_bins = n_bins_vals,
            min_per_bin = min_per_bin_vals,
            sorted = FALSE
        )

        # For each parameter combo, count passing species
        results <- param_grid[,
            {
                cov <- coverage_grid[n_bins == .BY$n_bins & min_per_bin == .BY$min_per_bin]
                merged <- sp_stats[cov, on = "Mnemonic"]
                pass <- merged[n_tags >= .BY$min_tags &
                    n_valid_growth_pair >= .BY$min_pairs &
                    passes_coverage == TRUE, Mnemonic]
                .(
                    n_species = length(pass),
                    pct_species = round(100 * length(pass) / n_total, 1),
                    species_list = list(pass)
                )
            },
            by = .(min_tags, min_pairs, n_bins, min_per_bin)
        ]

        setattr(results, "n_total_species", n_total)
        setattr(results, "growth_forms", growth_forms)
        setorder(results, -n_species, min_tags, min_pairs, n_bins, min_per_bin)
        results
    }

    # # Sweep for trees/shrubs
    # sweep_ts <- sweep_bio_thresholds(to_do, growth_forms = "trees_shrubs")
    # # Sweep for palms
    # sweep_palm <- sweep_bio_thresholds(to_do, growth_forms = "palm")
    # # Sweep for figs
    # sweep_fig <- sweep_bio_thresholds(to_do, growth_forms = "fig")
    # # Sweep for strangler figs
    # sweep_strangler_fig <- sweep_bio_thresholds(to_do, growth_forms = "strangler_fig")
    # # Sweep for all growth forms combined
    # sweep_all <- sweep_bio_thresholds(to_do, growth_forms = unique(to_do$growth_form))

    min_tags <- 20
    min_pairs <- 20
    n_bins <- 4
    min_per_bin <- 3

    ## ---- Get Bio Pars per species within growth_form -------------------------
    # Per-species estimates for trees/shrubs (species with sufficient data)
    to_do_tree_shrub_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = "trees_shrubs",
        species_label = NULL,
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # Pooled fallback for trees/shrubs
    to_do_tree_shrub_all_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = c("trees_shrubs"), # include unknown growth form in pooled tree/shrub group
        species_label = "all_trees_shrubs",
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # per species estimates for palms (species with sufficient data)
    to_do_palm_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = "palm",
        species_label = NULL,
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # Pooled palms + tree ferns (NOTE: no tree fern data in census 7, 8, 9)
    to_do_palm_tree_fern_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = c("palm", "tree_fern"),
        species_label = "all_palm_tree_fern",
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # per species estimates for figs (species with sufficient data)
    to_do_fig_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = "fig",
        species_label = NULL,
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # Pooled figs (relax tag requirement; rely on size-coverage filter)
    to_do_fig_all_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = "fig",
        species_label = "all_fig",
        min_tags = min_tags, min_pairs = min_pairs,
        n_bins = n_bins, min_per_bin = min_per_bin
    )

    # # Inspect which species pass for a specific combo
    # # sweep_strangler_fig[min_tags == 20 & min_pairs == 20 & n_bins == 4 & min_per_bin == 3, species_list]
    # sweep_strangler_fig[min_tags == 5 & min_pairs == 3 & n_bins == 3 & min_per_bin == 3, species_list]

    # Pooled strangler figs (relax tag requirement; rely on size-coverage filter)
    to_do_strangler_fig_all_sp <- prepare_bio_data(
        dt = to_do,
        growth_forms = "strangler_fig",
        species_label = "all_strangler_fig",
        min_tags = 10, min_pairs = 5,
        n_bins = 3, min_per_bin = 3
    )

    # Function to run bio parameter estimation for a given data.table and return a
    # list of results keyed by species. This allows us to run the estimation
    # separately for trees, shrubs, palms, etc., and then combine the results as
    # needed.
    run_bio_par_estimation <- function(dt, verbose = FALSE) {
        res <- vector("list", length = uniqueN(dt$species))
        names(res) <- unique(dt$species)
        sp_names <- names(res)
        for (i in seq_along(sp_names)) {
            sp <- sp_names[i]
            res[[sp]] <- estimate_bio_pars(
                dt[species == sp],
                use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
                # Hard shrink guardrail (max_shrink)
                # - "data": estimated from observed shrink tail (with measurement-error support)
                # - "fixed": use a fixed constant bound (cm/year)
                max_shrink_source = MAX_SHRINK_HARD_SOURCE,
                max_shrink_fixed = MAX_SHRINK_FIXED,
                # Soft shrinkage penalty strength (k_shrink)
                # - "data": estimate from measurement-error scale (preferred) or from data variance
                # - "fixed": use a fixed constant (units: 1/cm^2)
                k_shrink_source = K_SHRINK_SOURCE,
                k_shrink_fixed = K_SHRINK_FIXED,
                # Soft extreme-growth penalty strength (k_growth), analogous to k_shrink
                # - "data": estimate from measurement-error scale (preferred) or from data variance
                # - "fixed": use a fixed constant (units: 1/cm^2); set 0 to disable soft penalty
                k_growth_source = K_GROWTH_SOURCE,
                k_growth_fixed = K_GROWTH_FIXED,
                # Hard growth guardrail (max_growth)
                # - "data": estimated from observed extreme-growth tail
                # - "fixed": use a fixed constant bound (cm/year)
                max_growth_source = MAX_GROWTH_HARD_SOURCE,
                max_growth_fixed = MAX_GROWTH_FIXED,
                # Quantiles used to set conservative guardrails
                # the lowest value of the probability function to get lowest shrink from measurement error
                shrink_hard_prob = 1e-4,
                # the lowest value of the empirical quantile to get lowest shrink from data
                shrink_data_quantile = 0.001,
                # if masurement error, then, lowest shrink is the min between the two for hard shrink guardrail
                #################
                # Extreme-growth guardrails (upper tail)
                # - growth_hard_prob is the *upper-tail* probability (e.g., 1e-4 means 99.99th percentile)
                # - growth_data_quantile is the empirical upper quantile used as a guardrail
                # - growth_soft_quantile sets a softer threshold used for a quadratic penalty
                # the highest (1 - growth_hard_prob) value of the probability function to get highest growth from measurement error
                growth_hard_prob = 1e-4,
                # Upper quantile for hard growth guardrail from empirical data
                growth_data_quantile = 0.999,
                # if masurement error, then, highest growth is the max between the two for hard growth guardrail
                #################
                # to etimate the growth soft penalty k_growth - its used if it becomes the minimum between max grwoth from measurement error or fixed or data
                growth_soft_quantile = 0.99,
                # Recruitment max DBH (upper bound for recruits dbh at first census)
                recruit_max_quantile = 0.999,
                recruit_max_source = RECRUIT_MAX_SOURCE,
                recruit_max_fixed = as.numeric(RECRUIT_MAX_FIXED),
                # -----------------------------------------------------------------
                # Optional enforcement of user-specified growth/recruit bounds
                # Units: growth bounds in cm/year; recruit max in cm.
                # If 'enforce_growth_bounds' is TRUE, observations outside the provided fixed
                # bounds ('growth_min_fixed' and/or 'growth_max_fixed') are dropped before estimation.
                enforce_growth_bounds = TRUE,
                growth_min_fixed = MAX_SHRINK_FIXED,
                growth_max_fixed = MAX_GROWTH_FIXED,
                # If 'enforce_recruit_max' is TRUE, recruits with DBH > 'recruit_max_fixed' are dropped
                # before fitting the recruitment-size lognormal.
                enforce_recruit_max = TRUE
            )
            if (verbose) {
                cat("Estimated bio parameters for species:", sp, "\n")
            }
        }
        res
    }

    bio_pars_tree_shrub_sp <- run_bio_par_estimation(to_do_tree_shrub_sp, verbose = TRUE)
    bio_pars_tree_shrub_all_sp <- run_bio_par_estimation(to_do_tree_shrub_all_sp, verbose = TRUE)
    bio_pars_palm_sp <- run_bio_par_estimation(to_do_palm_sp, verbose = TRUE)
    bio_pars_palm_tree_fern_sp <- run_bio_par_estimation(to_do_palm_tree_fern_sp, verbose = TRUE)
    bio_pars_fig_sp <- run_bio_par_estimation(to_do_fig_sp, verbose = TRUE)
    bio_pars_fig_all_sp <- run_bio_par_estimation(to_do_fig_all_sp, verbose = TRUE)
    bio_pars_strangler_fig_all_sp <- run_bio_par_estimation(to_do_strangler_fig_all_sp, verbose = TRUE)
    # For unknown species, use the tree all species parameters as a fallback, but
    # with a distinct name to avoid confusion in downstream analyses.
    bio_pars_unknown_sp <- bio_pars_tree_shrub_all_sp
    names(bio_pars_unknown_sp) <- "all_unknown"

    # Combine bio pars
    bio_pars <- c(
        # tree+shrub species parameters
        bio_pars_tree_shrub_sp,
        bio_pars_tree_shrub_all_sp,
        # palm parameters
        bio_pars_palm_sp,
        # palm + tree fern parameters (these are pooled together since tree ferns have no data in recent censuses)
        bio_pars_palm_tree_fern_sp,
        # fig parameters
        bio_pars_fig_sp,
        bio_pars_fig_all_sp,
        # strangler fig parameters
        bio_pars_strangler_fig_all_sp,
        # unkonwn species parameters (use tree all sp as fallback, but with a distinct name)
        bio_pars_unknown_sp
    )

    # Normalize Mnemonic
    # stop if any Mnemonic is NA or infinite in xraw, since that would cause issues downstream
    if (any(is.na(xraw$Mnemonic) | is.infinite(xraw$Mnemonic))) {
        stop("Error: 'Mnemonic' column contains NA or infinite values. Please clean the data before proceeding.")
    }

    # if species tree name in the species specific tree/shrub estimates, keep mnemonic
    # if not tree or shrub, set to non_tree_shrub (these will be run with igraph and have no species-specific bio pars)
    xraw[
        growth_form == "trees_shrubs",
        species := fifelse(
            Mnemonic %in% names(bio_pars_tree_shrub_sp),
            Mnemonic, "all_trees_shrubs"
        )
    ]

    # unique(xraw[is.na(species)]$growth_form)

    xraw[
        growth_form == "palm",
        species := fifelse(
            Mnemonic %in% names(bio_pars_palm_sp),
            Mnemonic, "all_palm_tree_fern"
        )
    ]

    # unique(xraw[is.na(species)]$growth_form)

    xraw[
        growth_form %in% c("tree_fern"),
        species := "all_palm_tree_fern"
    ]

    # unique(xraw[is.na(species)]$growth_form)

    xraw[
        growth_form == "fig",
        species := fifelse(
            Mnemonic %in% names(bio_pars_fig_sp),
            Mnemonic, "all_fig"
        )
    ]

    # unique(xraw[is.na(species)]$growth_form)

    xraw[
        growth_form == "strangler_fig",
        species := "all_strangler_fig"
    ]

    # unique(xraw[is.na(species)]$growth_form)

    xraw[
        growth_form == "unknown",
        species := "all_unknown"
    ]

    # Convert and normalize columns
    xraw_single_stems <- xraw[single_stem_tags == TRUE]
    xraw_single_stems[, `:=`(
        ExactDate = as.IDate(ExactDate),
        DBH_mm_original_backup = DBH,
        DBH = dbh_with_best_candidate_taper_corrected / 10,
        CensusID = as.integer(CensusID),
        StemID = as.integer(as.character(StemID)),
        TrueStemID = StemID
    )]

    xraw_multi_stems <- xraw[single_stem_tags == FALSE]
    xraw_multi_stems[, `:=`(
        ExactDate = as.IDate(ExactDate),
        DBH_mm_original_backup = DBH,
        DBH = dbh_with_best_candidate_taper_corrected / 10,
        CensusID = as.integer(CensusID),
        StemID = as.integer(as.character(StemID)) # ,
        # TrueStemID = as.integer(as.character(fifelse(CensusID >= ANCHOR_START_CENSUS & !is.na(DBH), StemID, NA)))
        # NOTE: this gives truestemid to NA dbh, which is fine, i have safeguards. but I dont want that now
        # TrueStemID = as.character(fifelse(CensusID >= 7, StemID, NA))
    )]

    # The truestemid for multistem trees are those stemid after census 7 (2010) and
    # the stemid for those stems with a stemtag (since those have always been
    # identified). This means that if a multistem tree has a stemid before census 7, but
    # no stemtag, this stemid was given in the computer and need to be checked.
    # However, if a stem has a stemid before census 7 and also a stemtag, we don't
    # need to track those as those have always had id.
    xraw_multi_stems[
        ,
        TrueStemID := fcase(
            is.na(DBH), NA_integer_, # No measurement → no reliable ID
            !is.na(StemTag), StemID, # Has physical tag → trust it
            CensusID >= 7, StemID, # Census 7+ → systematic tracking
            default = NA_integer_ # Early census, no tag → unreliable
        )
    ]

    # SAFETY CODE:
    # fill missing CensusDate with mean per CensusID
    xraw_multi_stems[, mean_date_census :=
        as.IDate(mean(as.numeric(ExactDate), na.rm = TRUE),
            origin = "1970-01-01"
        ),
    by = CensusID
    ]

    xraw_multi_stems[
        is.na(ExactDate),
        ExactDate := as.IDate(mean_date_census)
    ][, mean_date_census := NULL]

    xraw_multi_stems <- ensure_species_column(xraw_multi_stems)
    # print(inspectdf::inspect_na(xraw), n = 60)

    xrun <- data.table::copy(xraw_multi_stems)

    # ------------------------------------------------------------------
    # free memory: the remaining objects are only needed to build xrun/
    # bio_pars, which we now have; drop the big tables before the DP work
    # ------------------------------------------------------------------
    rm(
        xraw,
        anchor_data, mean_date_tag_census, mean_date_census, full_anchor_data,
        tag_stem_combos, tag_mnemonic,
        to_do,
        to_do_tree_shrub_sp, to_do_tree_shrub_all_sp,
        to_do_palm_tree_fern_sp, to_do_fig_sp, to_do_fig_all_sp,
        to_do_strangler_fig_all_sp,
        bio_pars_tree_shrub_sp, bio_pars_tree_shrub_all_sp,
        bio_pars_palm_sp, bio_pars_palm_tree_fern_sp,
        bio_pars_fig_sp, bio_pars_fig_all_sp,
        bio_pars_strangler_fig_all_sp,
        bio_pars_unknown_sp
    )
    invisible(gc())


    # now continue with parameter logging / attachment, etc.
    # Write a small text file recording the parameters used to build the
    # run-specific output directory name so runs are reproducible.
    # Gather all parameters in the environment for full run context

    params <- list()
    params$TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    params$OUT_DIR <- out_dir
    params$GENERATED_DIR_NAME <- basename(out_dir)
    # Record BATCH_TS explicitly: use provided BATCH_TS if present; otherwise
    # extract the timestamp prefix used in the generated directory name (YYYYmmdd_HHMMSS)
    if (exists("BATCH_TS") && nzchar(BATCH_TS)) {
        params$BATCH_TS <- BATCH_TS
    } else {
        m <- regexpr("^[0-9]{8}_[0-9]{6}", basename(out_dir))
        params$BATCH_TS <- if (m[1] == -1) "" else regmatches(basename(out_dir), m)
    }

    # Add called parameters (command-line overrides)
    params$CALLED_PARAMETERS <- overrides

    # Add tag values and unique tags/species in data
    # params$UNIQUE_TAGS <- if (exists("xrun")) unique(xrun$Tag) else NA
    params$UNIQUE_SPECIES <- if (exists("xrun")) unique(xrun$species) else NA
    # Record the exact command line used to invoke this run
    # This is critical for full reproducibility
    params$COMMAND_LINE <- paste(commandArgs(), collapse = " ")

    # Add a readable version of bio_pars to the log
    if (exists("bio_pars")) {
        params$BIO_PARAMETERS_STR <- paste(capture.output(str(bio_pars)), collapse = "\n")
    }

    format_param <- function(x) {
        if (is.null(x)) {
            "NULL"
        } else if (is.function(x)) {
            "<function>"
        } else if (is.environment(x)) {
            "<environment>"
        } else if (is.list(x)) {
            "<list>"
        } else if (length(x) == 0) {
            "<empty>"
        } else if (is.logical(x)) {
            if (isTRUE(x)) "TRUE" else "FALSE"
        } else if (is.atomic(x) && length(x) == 1 && is.na(x)) {
            "NA"
        } else if (is.atomic(x)) {
            paste(as.character(x), collapse = ", ")
        } else {
            paste0("<", typeof(x), ">")
        }
    }
    params_lines <- unlist(
        Map(function(name, value) {
            paste0(name, ": ", format_param(value))
        }, names(params), params)
    )
    params_file <- file.path(out_dir, "run_parameters_full.txt")
    writeLines(params_lines, con = params_file)

    # 5.3 Attach Bio_* columns (DP reads parameters from columns)
    xrun <- attach_bio_columns(xrun, bio_pars)

    # NOTE: Return the same mnemonic names to species
    xrun[, species := Mnemonic]

    # bci_multistem_xrun_debug <- xrun

    # # rename needed columns
    # bci_multistem_xrun_debug <- bci_multistem_xrun_debug[, .(
    #     Species = Mnemonic, Tag, OriginalStemID = StemID, TrueStemID, CensusID, ExactDate,
    #     DBH, ListOfTSM, HOM,
    #     Bio_Mu_Growth, Bio_Gamma_Growth, Bio_Sigma0_Growth, Bio_Sigma1_Growth,
    #     Bio_H0_Mortality, Bio_Beta_Mortality,
    #     Bio_Recruit_Meanlog, Bio_Recruit_Sdlog, Bio_Recruit_MaxDBH_unit, Bio_Recruitment_lambda,
    #     Bio_Max_Shrink, Bio_K_Shrink,
    #     Bio_Max_Growth, Bio_Max_Growth_Soft, Bio_K_Growth
    # )]

    # vars_to_export <- c(
    #     "Species", "Tag", "OriginalStemID", "TrueStemID", "CensusID", "ExactDate",
    #     "DBH", "ListOfTSM", "HOM",
    #     grep("^Bio_", names(xrun), value = TRUE)
    # )

    # bci_multistem_xrun_debug <- bci_multistem_xrun_debug[, ..vars_to_export]
    # saveRDS(bci_multistem_xrun_debug, "./2_STEM_IDENTIFICATION/bci_multistem_xrun_debug.rds")

    #     # # select Mnemonic and all bio parameter columns
    #     # cat("Exporting bio parameters by species to CSV for DP use...\n")
    #     # cols_to_select <- c("growth_form", "Mnemonic", grep("^Bio_", names(xrun), value = TRUE))
    #     # write.csv(unique(xrun[, ..cols_to_select]), file.path(out_dir, "bio_parameters_by_species.csv"), row.names = FALSE)

    # # !
    # check_tags <- as.data.table(read.csv("./2_STEM_IDENTIFICATION/all_tags_to_re_process.csv", stringsAsFactors = FALSE))
    # xrun <- xrun[Tag %in% unique(check_tags$Tag)]
    # # !

    # 5.4 DP meta settings
    dp_max_tracks_local <- if (is.null(DP_MAX_TRACKS)) auto_dp_max_tracks(xrun) else as.integer(DP_MAX_TRACKS)
    dp_max_tracks_local <- as.integer(dp_max_tracks_local)

    ## Create chunks based on unique (Tag, species) groups
    if (!requireNamespace("parallel", quietly = TRUE)) stop("Package not available: parallel.")
    groups <- unique(xrun[, .(Tag, species)])
    data.table::setorder(groups, Tag, species)
    group_idx <- seq_len(nrow(groups))
    chunks <- split(group_idx, ceiling(seq_along(group_idx) / as.integer(DP_CHUNK_SIZE)))
    first_chunk <- !file.exists(DP_CSV_FILE)

    # Optionally limit to a subset of chunks for testing
    start_ci <- if (exists("DP_CHUNK_START") && !is.null(DP_CHUNK_START)) as.integer(DP_CHUNK_START) else 1L
    end_ci <- if (exists("DP_CHUNK_END") && !is.null(DP_CHUNK_END)) as.integer(DP_CHUNK_END) else length(chunks)

    # Handle edge cases: empty chunk list or out-of-range values
    if (length(chunks) == 0L) {
        log_msg("No groups/chunks to process — exiting.")
        tryCatch(
            {
                writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt"))
            },
            error = function(e) NULL
        )
        return(invisible(list(xrun = xrun, bio_pars = bio_pars)))
    }

    # Clamp values to valid range
    start_ci <- max(1L, min(length(chunks), start_ci))
    end_ci <- max(1L, min(length(chunks), end_ci))

    if (start_ci > end_ci) {
        stop(sprintf(
            "Invalid chunk range: DP_CHUNK_START=%s, DP_CHUNK_END=%s after clamping => start=%d > end=%d",
            if (exists("DP_CHUNK_START") && !is.null(DP_CHUNK_START)) as.character(DP_CHUNK_START) else "NULL",
            if (exists("DP_CHUNK_END") && !is.null(DP_CHUNK_END)) as.character(DP_CHUNK_END) else "NULL",
            start_ci, end_ci
        ))
    }

    for (ci in seq_len(length(chunks))) {
        if (ci < start_ci || ci > end_ci) {
            next
        }
        chunk_rds <- file.path(out_dir, sprintf(paste0(DP_BASE, "_chunk_%03d.rds"), ci))
        chunk_done <- file.path(out_dir, sprintf(paste0(DP_BASE, "_chunk_%03d_done.txt"), ci))

        if (isTRUE(DP_CHUNK_RESUME) && file.exists(chunk_done) && !isTRUE(DP_CHUNK_OVERWRITE)) {
            log_msg(sprintf("Skipping chunk %d/%d — already completed (resume enabled)", ci, length(chunks)))
            first_chunk <- !file.exists(DP_CSV_FILE)
            next
        }

        groups_ci <- groups[chunks[[ci]]]
        log_msg(sprintf("Chunk %d/%d — %d groups", ci, length(chunks), nrow(groups_ci)))

        # Run chunk with error isolation so one bad chunk doesn't kill the run
        status <- tryCatch(
            {
                res <- parallel::mclapply(seq_len(nrow(groups_ci)), function(j) {
                    data.table::setDTthreads(1L)
                    g <- groups_ci[j]
                    dtg <- xrun[Tag == g$Tag & species == g$species]

                    # Return rows as-is for groups with missing DBH or CensusID (all NA)
                    if (!("DBH" %in% names(dtg)) || !("CensusID" %in% names(dtg)) || all(is.na(dtg$DBH)) || all(is.na(dtg$CensusID))) {
                        log_msg(sprintf("Skipping Tag=%s species=%s in chunk %d: missing or all NA DBH/CensusID; returning rows as-is", g$Tag, g$species, ci), "WARN")
                        dtg[, ReconstructionMethod := "skipped_no_data"]
                        return(dtg)
                    }

                    run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local, chunk_id = ci)
                }, mc.cores = MC_CORES)

                # Remove NULLs (skipped groups) before binding
                res <- Filter(Negate(is.null), res)
                if (length(res) == 0L) {
                    out_chunk <- data.table::data.table()
                } else {
                    out_chunk <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
                }

                if (nrow(out_chunk) > 0L) {
                    out_chunk[, DP_Chunk := ci]
                    out_chunk <- maybe_add_posterior_bins(out_chunk)
                    # Record run output directory (basename) in each row to avoid variable/column name collision
                    out_chunk[, run_out_dir := basename(out_dir)]

                    if (isTRUE(WRITE_DP_CSV)) {
                        maybe_write(isTRUE(WRITE_DP_CSV), DP_CSV_FILE, function() {
                            data.table::fwrite(out_chunk, file = DP_CSV_FILE, append = !first_chunk)
                        }, sprintf("DP CSV chunk %d", ci))
                    }

                    if (isTRUE(WRITE_DP_FEATHER)) {
                        feather_path <- file.path(out_dir, sprintf(paste0(DP_BASE, "_chunk_%03d.feather"), ci))
                        maybe_write(isTRUE(WRITE_DP_FEATHER), feather_path, function() {
                            if (!requireNamespace("arrow", quietly = TRUE)) stop("'arrow' package not available; skipping feather output")
                            arrow::write_feather(out_chunk, feather_path)
                        }, sprintf("Feather for chunk %d", ci))
                    }

                    if (isTRUE(WRITE_DP_RDS)) {
                        maybe_write(isTRUE(WRITE_DP_RDS), chunk_rds, function() {
                            saveRDS(out_chunk, file = chunk_rds)
                        }, sprintf("RDS chunk %d", ci))
                    }

                    if (isTRUE(WRITE_DP_PDF_PER_CHUNK) && isTRUE(WRITE_DP_PDF)) {
                        pdf_path <- file.path(out_dir, sprintf(paste0(DP_BASE, "_chunk_%03d.pdf"), ci))
                        maybe_write(isTRUE(WRITE_DP_PDF_PER_CHUNK) && isTRUE(WRITE_DP_PDF), pdf_path, function() {
                            plot_tag_to_pdf(out_chunk, pdf_file = pdf_path, include_reference = DP_PDF_INCLUDE_REFERENCE)
                        }, sprintf("PDF chunk %d", ci))
                    }
                } else {
                    log_msg(sprintf("Chunk %d returned no rows.", ci))
                    if (isTRUE(WRITE_DP_RDS)) {
                        maybe_write(isTRUE(WRITE_DP_RDS), chunk_rds, function() {
                            saveRDS(out_chunk, file = chunk_rds)
                        }, sprintf("RDS chunk %d (empty)", ci))
                    }
                }

                # Write a completion marker so resume can skip this chunk
                # even when WRITE_DP_RDS=FALSE.
                tryCatch(
                    writeLines(as.character(Sys.time()), con = chunk_done),
                    error = function(e) log_msg(sprintf("Warning writing chunk done marker: %s", conditionMessage(e)), "WARN")
                )

                TRUE
            },
            error = function(e) {
                # Write a small error marker for the chunk and continue
                err_file <- file.path(out_dir, sprintf(paste0(DP_BASE, "_chunk_%03d_failed.txt"), ci))
                writeLines(conditionMessage(e), con = err_file)
                log_msg(sprintf("Chunk %d failed: %s", ci, conditionMessage(e)), "ERROR")
                FALSE
            }
        )

        first_chunk <- FALSE
        # Safe cleanup: only rm() variables that were successfully assigned
        for (.v in c("res", "out_chunk", "groups_ci")) {
            if (exists(.v, inherits = FALSE)) rm(list = .v)
        }
        invisible(gc())
    }

    # Write a small finished marker so users and wrappers can detect job completion
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt"))
            log_msg("Finished run")
        },
        error = function(e) {
            # message("[dp_global main_cpp_chunk.R] Warning writing run_finished marker: ", conditionMessage(e))
            log_msg(sprintf("Warning writing run_finished marker: %s", conditionMessage(e)), "WARN")
        }
    )

    invisible(list(xrun = xrun, bio_pars = bio_pars))
}

############################################################
### 10) Post-run utilities
############################################################
# Helpers for merging per-chunk RDS or Feather files into a single CSV.
# Call merge_chunks_to_csv(out_dir) after a run to produce a combined flat
# file without loading all chunks into memory simultaneously.

# Merge helpers: combine per-chunk RDS or Feather files into a single CSV
merge_chunk_rds_to_csv <- function(out_dir, out_csv = file.path(out_dir, paste0(DP_BASE, "_merged.csv"))) {
    files <- list.files(out_dir, pattern = paste0(DP_BASE, "_chunk_\\d{3}\\.rds$"), full.names = TRUE)
    if (length(files) == 0L) stop("No chunk RDS files found in ", out_dir)
    first <- TRUE
    for (f in sort(files)) {
        log_msg(paste("Merging RDS", basename(f)))
        dt <- readRDS(f)
        if (nrow(dt) == 0L) next
        if (first) {
            maybe_write(TRUE, out_csv, function() data.table::fwrite(dt, file = out_csv), sprintf("Merged RDS first write: %s", basename(out_csv)))
            first <- FALSE
        } else {
            maybe_write(TRUE, out_csv, function() data.table::fwrite(dt, file = out_csv, append = TRUE), sprintf("Merged RDS append: %s", basename(out_csv)))
        }
        rm(dt)
        invisible(gc())
    }
    log_msg(paste("Merged", length(files), "RDS chunks to", out_csv))
    out_csv
}

merge_chunk_feathers_to_csv <- function(out_dir, out_csv = file.path(out_dir, paste0(DP_BASE, "_merged.csv"))) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow package required to read feather files")
    files <- list.files(out_dir, pattern = paste0(DP_BASE, "_chunk_\\d{3}\\.feather$"), full.names = TRUE)
    if (length(files) == 0L) stop("No chunk Feather files found in ", out_dir)
    first <- TRUE
    for (f in sort(files)) {
        log_msg(paste("Merging Feather", basename(f)))
        dt <- as.data.table(arrow::read_feather(f))
        if (nrow(dt) == 0L) next
        if (first) {
            maybe_write(TRUE, out_csv, function() data.table::fwrite(dt, file = out_csv), sprintf("Merged Feather first write: %s", basename(out_csv)))
            first <- FALSE
        } else {
            maybe_write(TRUE, out_csv, function() data.table::fwrite(dt, file = out_csv, append = TRUE), sprintf("Merged Feather append: %s", basename(out_csv)))
        }
        rm(dt)
        invisible(gc())
    }
    log_msg(paste("Merged", length(files), "Feather chunks to", out_csv))
    out_csv
}

merge_chunks_to_csv <- function(out_dir, prefer = c("rds", "feather")) {
    prefer <- match.arg(prefer)
    if (prefer == "feather") {
        chars <- list.files(out_dir, pattern = paste0(DP_BASE, "_chunk_\\d{3}\\.feather$"), full.names = TRUE)
        if (length(chars) > 0L) {
            return(merge_chunk_feathers_to_csv(out_dir))
        }
        return(merge_chunk_rds_to_csv(out_dir))
    }
    # prefer rds
    rds <- list.files(out_dir, pattern = paste0(DP_BASE, "_chunk_\\d{3}\\.rds$"), full.names = TRUE)
    if (length(rds) > 0L) {
        return(merge_chunk_rds_to_csv(out_dir))
    }
    return(merge_chunk_feathers_to_csv(out_dir))
}

############################################################
### 11) Entrypoint
############################################################
# Executed only when called directly via Rscript; safe to source() without
# triggering a run (sys.nframe() > 0 when sourced interactively).
if (sys.nframe() == 0L) {
    message("[dp_global main_cpp_chunk_bci.R] Starting chunked run_main_chunked()")
    run_main_chunked()
}

############################################################
### Usage guide
############################################################
#
# PURPOSE
#   Memory-efficient chunked DP pipeline. Splits all (Tag, species) groups into
#   batches of DP_CHUNK_SIZE, processes each batch independently, and writes a
#   per-chunk RDS file plus an incrementally-appended CSV. Use this script
#   instead of main_cpp.R when the full data set does not fit comfortably in RAM
#   or when you want fault-tolerant checkpointing.
#
# BASIC RUN — all chunks, default settings
#   Rscript dp_global/scripts/main_cpp_chunk.R \
#     --INPUT_FILE=data/my_stems.csv
#
# COMMON OVERRIDES
#   --DP_CHUNK_SIZE=7             Groups per chunk (increase for faster runs with enough RAM)
#   --DP_MAX_STATES=1100          Max DP states per track
#   --MANUAL_CORES=TRUE           Use a fixed core count instead of auto-detect
#   --MANUAL_CORES_VALUE=8        Number of cores (requires MANUAL_CORES=TRUE)
#   --WRITE_DP_PDF=TRUE/FALSE     Whether to produce per-chunk PDF plots
#   --WRITE_DP_PDF_PER_CHUNK=TRUE Per-chunk PDFs in addition to any run-level PDF
#   --POSTERIOR_SAMPLES=250       Draw N posterior samples per group (0 = disabled)
#   --POSTERIOR_SAMPLES_FORMAT=csv|rds|feather
#   --USE_MEASUREMENT_ERROR=TRUE  Enable measurement-error model for bio params
#   --DP_FALLBACK_GROWTH_FORMS="fig,tree"
#                                 Comma-separated list of growth forms that bypass
#                                 species-specific bio params-falls back to probabilistic
#   --DP_CHUNK_START=3            Start from chunk N (skip earlier chunks)
#   --DP_CHUNK_END=9              Stop after chunk N
#
# STOPPING A RUN
#   Send SIGINT (Ctrl-C in the terminal) or SIGTERM to the Rscript process.
#   The completion marker (_done.txt) is written as the very last step of each
#   chunk, after all outputs (RDS, CSV, feather, PDF) have been flushed.
#   A chunk without a _done.txt is treated as incomplete and will be re-run
#   on resume, even if a partial _chunk_NNN.rds exists from the interrupted run.
#   Files in the output directory:
#     stem_reconstruction_dp_global_rcpp_chunk_NNN.rds   — chunk data
#     stem_reconstruction_dp_global_rcpp_chunk_NNN_done.txt — completion flag
#                                                             (only present if chunk finished)
#
# RESUMING A STOPPED RUN
#   Pass the path of the existing output directory via --OUT_DIR_OVERRIDE and
#   enable the resume flag. The script will skip every chunk that already has a
#   _done.txt marker and continue from where it stopped.
#
#   Rscript dp_global/scripts/main_cpp_chunk.R \
#     --INPUT_FILE=data/my_stems.csv \
#     --OUT_DIR_OVERRIDE=dp_global/output/20260329_220720_unknown_T0_DP_MB_ME_g5_sm0p5_kg0_ks0_rcpp \
#     --DP_CHUNK_RESUME=TRUE \
#     --MANUAL_CORES=TRUE \
#     --MANUAL_CORES_VALUE=8
#
#   Important: all non-output parameters (DP_MAX_STATES, POSTERIOR_SAMPLES,
#   DP_FALLBACK_GROWTH_FORMS, etc.) must match the original run so that the
#   resumed chunks are processed identically to the ones already completed.
#
# MERGING CHUNK FILES AFTER A RUN
#   If the incremental CSV is missing or incomplete, rebuild it from the RDS
#   checkpoints without re-running the DP:
#
#   source("dp_global/scripts/main_cpp_chunk.R")   # defines helpers only
#   merge_chunk_rds_to_csv("dp_global/output/<run_dir>")
#
# HELP
#   Rscript dp_global/scripts/main_cpp_chunk.R --help