############################################################
### main_cpp.R — dp_global driver
############################################################
# Goal
#   One place to run the DP_GLOBAL workflow end-to-end
#
# Note for orchestrators
# - This script accepts CLI overrides of internal variables via --KEY=VALUE.
# - See the `CLI_REFERENCE` variable below for the canonical keys used by
#   external orchestrators; they should construct flags matching these canonical
#   names (case-insensitive, '-' or '_' allowed).
#   Keep the orchestrator in sync with `CLI_REFERENCE`.
#
# Table of Contents
#  0) Housekeeping           — safe top-level behavior
#  1) CLI parsing            — parse and coerce command-line overrides
#  2) Dependencies           — package checks and imports
#  3) Defaults
#    3.1) Biological parameter estimation settings
#    3.2) DP solver settings
#    3.3) Parallelism settings
#    3.4) Output & path settings
#  4) CLI reference & override mapping
#    4.1) Sensitivity & realism settings
#    4.2) Help & override application
#  5) Input validation       — abort early on missing/invalid inputs
#  6) Source project code    — load dp_global R modules
#  7) Helpers                — data-manipulation utilities
#  8) Core DP functions      — run_dp_one_group() and helpers
#  9) Main pipeline          — run_main()
# 10) Entrypoint             — execute when invoked via Rscript
#
# Use the numbered sections to quickly scan the file. Each section contains a
# short header and concise responsibilities to make navigation quick and clear.

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
# Usage: Rscript main_cpp.R --WHICH_TAG=2 --RUN_ALL_TAGS=TRUE --DP_MODE=marginals+bins --SENSITIVITY_MODE=run --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=8
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
                .char_keys <- c(
                    "WHICH_TAG", "PROB_SPECIES", "DP_FALLBACK_GROWTH_FORMS",
                    "NON_TAPER_CORRECTED_GROWTH_FORMS", "CONFIG_NAME",
                    "INPUT_FILE", "POSTERIOR_SAMPLES_FORMAT", "SPECIES_COL"
                )
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

############################################################
### 3) Defaults
############################################################
# Editable run defaults and output-naming variables shared across all sections.
INPUT_FILE <- here("data_simulation", "data", "simulated_data_1.csv")
FORCE_ONE_SPECIES_PARAMETERS <- TRUE
if (isTRUE(FORCE_ONE_SPECIES_PARAMETERS)) {
    FORCED_SPECIES_LABEL <- "all"
    message("[dp_global main_cpp.R] FORCE_ONE_SPECIES_PARAMETERS=TRUE: using single species label '", FORCED_SPECIES_LABEL, "' for all trees.")
} else {
    message("[dp_global main_cpp.R] FORCE_ONE_SPECIES_PARAMETERS=FALSE: using species column from data for parameter estimation.")
}
SPECIES_COL <- NULL

############################################################
### 3.1) Biological parameter estimation settings
############################################################
# Growth rates, shrinkage limits, and soft-penalty k values.
# NOTE: You can define them with parameter data from your species of interest.
USE_MEASUREMENT_ERROR <- TRUE
MAX_GROWTH_HARD_SOURCE <- "fixed"
MAX_GROWTH_FIXED <- 7.5
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
# DP algorithm parameters: mode, anchoring, state budget, slack, posterior sampling.
DP_MODE <- "marginals+bins" # Options: "none", "marginals", "marginals+bins"
WHICH_TAG <- "20"
ANCHOR_START_CENSUS <- 7L
DP_VERBOSE <- TRUE
DP_POSTERIOR_TOP_K <- 2L
DP_MAX_TRACKS <- NULL # auto (computed from data)
DP_MAX_STATES <- 40000L
DP_SLACK_TRACKS <- 1L
# NOTE: Optionally require that slack be granted only if an anchor DBH is recruitable
# (i.e., DBH <= Bio_Recruit_MaxDBH_unit + eps). Set FALSE to preserve current behavior.
DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE <- TRUE
# Tolerance (cm) used when comparing anchor DBH to recruit_max_dbh
DP_SLACK_REQUIRE_ANCHOR_EPS <- 1e-6
# Growth forms which should bypass DP and force probabilistic matcher
# - character vector; values correspond to entries in the `growth_form`
#   column of the input dataset.  Pass to the DP function via
#   `fallback_growth_forms` argument.
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
# - POSTERIOR_SAMPLE_SEED: integer seed used to make posterior sampling reproducible. If NULL sampling is not
#   deterministically seeded. Pass --POSTERIOR_SAMPLE_SEED=<int> for reproducible CLI runs; the chunked runner
#   (main_cpp_chunk.R) defaults to 123L when sampling is enabled and no seed is set.
POSTERIOR_SAMPLES <- 200L
POSTERIOR_SAMPLES_FORMAT <- "csv" # options: 'rds', 'feather', 'csv'
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- NULL
# Option: allow DP to use a provisional anchor at the last observed DBH census when no TrueStemID exists
ALLOW_PROVISIONAL_DP_ANCHOR <- TRUE

# Number of stochastic samples drawn by the probabilistic greedy matcher
# when DP falls back due to intractable state spaces.
PROB_N_SAMPLES <- 200L

# Species that should bypass DP and go directly to the probabilistic greedy
# matcher. Provide a character vector of Species column values (e.g.,
# c("Oenocarpus mapora", "Socratea exorrhiza")). Empty vector disables.
PROB_SPECIES <- character(0)

# Lookahead weight for probabilistic matcher: controls how much the
# cost matrix at pair (i) is influenced by the assignment at pair (i+1).
# 0 = disabled (original independent sampling), 0.5 = default.
PROB_LOOKAHEAD_WEIGHT <- 1

# Bio hard bounds control for probabilistic matcher:
# When TRUE (default), use bio-estimated hard shrink/growth guardrails (strict).
# When FALSE, rely only on prune bounds (relaxed, for continuity rescue).
USE_BIO_HARD_SHRINK_IN_PROB <- TRUE
USE_BIO_HARD_GROWTH_IN_PROB <- TRUE

# ME cumulative-shrinkage threshold for probabilistic matcher (Layer 2 repair).
# Trajectories where cumulative shrinkage exceeds n_sigma_me * sqrt(SD(d_start)^2 + SD(d_curr)^2)
# are severed. Lower values = sever sooner (less negative growth). 0 = sever at first decline.
# Only active when USE_BIO_HARD_SHRINK_IN_PROB = TRUE.
PROB_N_SIGMA_ME <- 3

# Pin observations with known TrueStemID at non-anchor censuses to their
# correct track.  Reduces state space and prevents misidentification.
# Set FALSE to revert to pre-pinning behavior.
PIN_TRUESTEMID <- TRUE

############################################################
### 3.3) Parallelism settings
############################################################
RUN_ALL_TAGS <- FALSE
MANUAL_CORES <- TRUE # Flag to manually define cores instead of auto-detecting
MANUAL_CORES_VALUE <- 1L # Number of cores to use if MANUAL_CORES=TRUE

############################################################
### 3.4) Output & path settings
############################################################
# Base output directory
base_out_dir <- here("dp_global", "output")
message("[dp_global main_cpp.R] here root: ", here::here())
message("[dp_global main_cpp.R] base_out_dir (raw): ", base_out_dir)
base_out_dir <- normalizePath(base_out_dir, winslash = "/", mustWork = FALSE)
message("[dp_global main_cpp.R] base_out_dir (normalized): ", base_out_dir)

# Optional: explicitly set a subdirectory name for outputs.
# If NULL, an automatic name based on timestamp + key config flags is used.
# OUT_DIR_NAME <- NULL
# CONFIG_NAME is set by the orchestrator (e.g., run_dp_future) to identify the
# experimental configuration; default to NULL so override parsing treats it as
# a valid, known variable rather than an unknown override.
CONFIG_NAME <- NULL
# When non-NULL, OUT_DIR_OVERRIDE bypasses build_out_dir_name() and uses this
# path directly as out_dir. Pass --OUT_DIR_OVERRIDE=/path/to/existing/run on
# the command line to reuse an existing output directory.
OUT_DIR_OVERRIDE <- NULL

# Output path helpers (`encode_num`, `build_out_dir_name`) live in naming_helpers.R,
# sourced near the end of this section.

WRITE_DP_CSV <- TRUE
WRITE_DP_RDS <- TRUE
WRITE_DP_FEATHER <- FALSE
WRITE_DP_PDF <- TRUE
DP_PDF_INCLUDE_REFERENCE <- TRUE

if (!isTRUE(RUN_ALL_TAGS)) {
    PLOT_PDF_ONE_TAG_ONLY <- TRUE
} else {
    PLOT_PDF_ONE_TAG_ONLY <- FALSE
}

# Default project root so --PROJECT_ROOT=/path overrides are accepted by the CLI parser
PROJECT_ROOT <- here::here()
# Batch timestamp can be provided by orchestrators; default empty so overrides like --BATCH_TS=... are accepted without warnings
BATCH_TS <- ""
# Naming helpers (encode_num, build_out_dir_name) live in a separate helper
# file to keep the main script concise. Source it early so it's available
# when we compute `out_dir` below.
source(here("dp_global", "R", "naming_helpers.R"))

############################################################
### 4) CLI reference & override mapping
############################################################
# Canonical CLI flags for external orchestrators.
# Keys are case-insensitive and may use '-' or '_' as separators.
CLI_REFERENCE <- list(
    INPUT_FILE = "INPUT_FILE",
    FORCE_ONE_SPECIES_PARAMETERS = "FORCE_ONE_SPECIES_PARAMETERS",
    DP_MODE = "DP_MODE",
    WHICH_TAG = "WHICH_TAG",
    ANCHOR_START_CENSUS = "ANCHOR_START_CENSUS",
    DP_VERBOSE = "DP_VERBOSE",
    RUN_ALL_TAGS = "RUN_ALL_TAGS",
    MANUAL_CORES = "MANUAL_CORES",
    MANUAL_CORES_VALUE = "MANUAL_CORES_VALUE",
    WRITE_DP_CSV = "WRITE_DP_CSV",
    WRITE_DP_RDS = "WRITE_DP_RDS",
    WRITE_DP_FEATHER = "WRITE_DP_FEATHER",
    WRITE_DP_PDF = "WRITE_DP_PDF",
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
    SENSITIVITY_MODE = "SENSITIVITY_MODE",
    RUN_REALISM_REPORT = "RUN_REALISM_REPORT",
    RUN_K_SWEEP_DEMO = "RUN_K_SWEEP_DEMO",
    PROJECT_ROOT = "PROJECT_ROOT",
    BATCH_TS = "BATCH_TS",
    CONFIG_NAME = "CONFIG_NAME",
    USE_MEASUREMENT_ERROR = "USE_MEASUREMENT_ERROR",
    ALLOW_PROVISIONAL_DP_ANCHOR = "ALLOW_PROVISIONAL_DP_ANCHOR",
    DP_MAX_TRACKS = "DP_MAX_TRACKS",
    OUT_DIR_OVERRIDE = "OUT_DIR_OVERRIDE",
    PROB_N_SAMPLES = "PROB_N_SAMPLES",
    PROB_SPECIES = "PROB_SPECIES",
    PROB_LOOKAHEAD_WEIGHT = "PROB_LOOKAHEAD_WEIGHT",
    USE_BIO_HARD_SHRINK_IN_PROB = "USE_BIO_HARD_SHRINK_IN_PROB",
    USE_BIO_HARD_GROWTH_IN_PROB = "USE_BIO_HARD_GROWTH_IN_PROB",
    PROB_N_SIGMA_ME = "PROB_N_SIGMA_ME",
    PIN_TRUESTEMID = "PIN_TRUESTEMID"
)

############################################################
### 4.1) Sensitivity & realism settings
############################################################
SENSITIVITY_MODE <- "none" # Options: "none", "run", "run+write", "run+write+pdf"
RUN_K_SWEEP_DEMO <- FALSE
RUN_REALISM_REPORT <- FALSE

############################################################
### 4.2) Help & override application
############################################################
print_help <- function() {
    cat("Usage: Rscript scripts/main_cpp.R [--KEY=VALUE] [--FLAG]\n")
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
        warning(sprintf("[dp_global main_cpp.R] Unknown override '%s' (ignored).\n", name))
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
    message("[dp_global main_cpp.R] Overriding ", match_var, " = ", as.character(new_val))
}

# Use canonical ALL-CAPS variables (e.g., WHICH_TAG, INPUT_FILE) everywhere.

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

# Optional override: explicitly set project root (useful when running under different working dirs or job launchers)
# Usage: --PROJECT_ROOT=/absolute/path/to/project
if (exists("PROJECT_ROOT") && !is.null(PROJECT_ROOT) && nzchar(PROJECT_ROOT)) {
    message("[dp_global main_cpp.R] Using PROJECT_ROOT override: ", PROJECT_ROOT)
    base_out_dir <- normalizePath(file.path(PROJECT_ROOT, "dp_global", "output"), winslash = "/", mustWork = FALSE)
    message("[dp_global main_cpp.R] base_out_dir overridden to: ", base_out_dir)
}

############################################################
### 5) Input validation
############################################################
# Abort early on missing or unreadable inputs; derive boolean mode flags.

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
RUN_SENSITIVITY <- SENSITIVITY_MODE != "none"
WRITE_OUTPUTS <- SENSITIVITY_MODE %in% c("run+write", "run+write+pdf")
MAKE_ALL_SWEEPS_PDF <- SENSITIVITY_MODE == "run+write+pdf"

# Final output directory for this run (created at runtime in run_main()).
# OUT_DIR_OVERRIDE takes precedence: set it to reuse an existing output directory.
if (!is.null(OUT_DIR_OVERRIDE) && nzchar(as.character(OUT_DIR_OVERRIDE))) {
    out_dir <- normalizePath(OUT_DIR_OVERRIDE, winslash = "/", mustWork = FALSE)
    message("[dp_global main_cpp.R] OUT_DIR_OVERRIDE set — using existing dir: ", out_dir)
} else {
    out_dir <- file.path(base_out_dir, build_out_dir_name())
    message("[dp_global main_cpp.R] out_dir (computed): ", out_dir)
}
message("[dp_global main_cpp.R] getwd(): ", getwd())

# Centralized DP naming and path helpers 🔧
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

# If posterior sampling is requested, default to the run base `out_dir` so the
# DP can create a single `posteriors/` subdirectory. If users supply a path
# that already ends in 'posteriors', strip that suffix to avoid double-nesting.
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

    # normalize path for consistency; don't create it yet (created when needed)
    POSTERIOR_SAMPLES_PATH <- normalizePath(POSTERIOR_SAMPLES_PATH, winslash = "/", mustWork = FALSE)
    if (is.null(POSTERIOR_SAMPLE_SEED)) {
        POSTERIOR_SAMPLE_SEED <- as.integer(123L)
    } else {
        POSTERIOR_SAMPLE_SEED <- as.integer(POSTERIOR_SAMPLE_SEED)
    }
} else {
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
source(here("dp_global", "R", "dp_global_main.R"))
source(here("dp_global", "R", "sensitivity_transition_cost_bio.R"))
source(here("dp_global", "R", "realism_calibration.R"))
source(here("dp_global", "R", "k_tuning_viz.R"))
# naming_helpers.R is already sourced in section 3.4; not repeated here.

############################################################
### 7) Helpers — data-manipulation utilities
############################################################

# Soft penalties vs hard guardrails (unit reminder)
# - Soft penalties operate on DBH differences (cm) over the interval.
# - Hard guardrails operate on annualized growth (cm/year).
#
# Choosing k from a reference excess:
# - Soft penalty is quadratic: soft_cost = k * (delta_cm^2)
# - If you want delta_cm = D to contribute cost C, set k = C / (D^2).

ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        message("[dp_global main_cpp.R] Created directory: ", path)
    } else {
        message("[dp_global main_cpp.R] Directory already exists: ", path)
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
            message("[dp_global main_cpp.R] Using '", species_col, "' as species column. Set SPECIES_COL to override.")
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
# run_dp_one_group(dtg, dp_max_tracks)  — wraps match_stems_dp_global_backward_marginals_batch()
# maybe_add_posterior_bins(out)          — optionally annotates output with posterior bin labels

run_dp_one_group <- function(dtg, dp_max_tracks) {
    # Safety: skip groups with missing DBH or CensusID (all NA or missing columns)
    tag_label <- if ("Tag" %in% names(dtg) && length(dtg$Tag) > 0) as.character(dtg$Tag[[1]]) else "<unknown>"
    species_label <- if ("species" %in% names(dtg) && length(dtg$species) > 0) as.character(dtg$species[[1]]) else "<unknown>"

    if (nrow(dtg) == 0L) {
        log_msg(sprintf("Skipping Tag=%s species=%s: 0 rows; nothing to process", tag_label, species_label), "WARN")
        return(data.table::copy(dtg))
    }

    if (!("DBH" %in% names(dtg)) || !("CensusID" %in% names(dtg))) {
        log_msg(sprintf("Skipping Tag=%s species=%s: missing DBH or CensusID column; returning rows as-is", tag_label, species_label), "WARN")
        out <- data.table::copy(dtg)
        out[, ReconstructionMethod := "skipped_no_data"]
        return(out)
    }
    if (all(is.na(dtg$DBH)) || all(is.na(dtg$CensusID))) {
        log_msg(sprintf("Skipping Tag=%s species=%s: all DBH or all CensusID are NA; returning rows as-is", tag_label, species_label), "WARN")
        out <- data.table::copy(dtg)
        out[, ReconstructionMethod := "skipped_no_data"]
        return(out)
    }

    out <- tryCatch(
        match_stems_dp_global_backward_marginals_batch(
            tree_data = data.table::copy(dtg),
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
            # growth-form bypass list
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
            prob_lookahead_weight = PROB_LOOKAHEAD_WEIGHT,
            use_bio_hard_shrink_in_prob = isTRUE(USE_BIO_HARD_SHRINK_IN_PROB),
            use_bio_hard_growth_in_prob = isTRUE(USE_BIO_HARD_GROWTH_IN_PROB),
            prob_n_sigma_me = as.numeric(PROB_N_SIGMA_ME),
            pin_truestemid = isTRUE(PIN_TRUESTEMID)
        ),
        error = function(e) {
            msg <- conditionMessage(e)
            log_msg(sprintf("DP error for Tag=%s species=%s: %s — falling back to probabilistic", tag_label, species_label, msg), "WARN")
            # Scope to pre-anchor rows only — mirrors do_fallback() inside the DP.
            # Post-anchor rows are appended afterward with proper given/none_after_anchor labels.
            .pre_anchor_eh <- dtg[CensusID <= ANCHOR_START_CENSUS]
            .post_anchor_eh <- dtg[CensusID > ANCHOR_START_CENSUS]
            out <- match_stems_probabilistic(
                tree_data = data.table::copy(.pre_anchor_eh),
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
                use_bio_hard_shrink_in_prob = isTRUE(USE_BIO_HARD_SHRINK_IN_PROB),
                use_bio_hard_growth_in_prob = isTRUE(USE_BIO_HARD_GROWTH_IN_PROB),
                n_sigma_me = as.numeric(PROB_N_SIGMA_ME),
                pin_truestemid = isTRUE(PIN_TRUESTEMID),
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
            # Append post-anchor rows with proper labeling (mirrors finalize_out / propagate_post_anchor_given)
            if (nrow(.post_anchor_eh) > 0L) {
                .post <- data.table::copy(.post_anchor_eh)
                if (!("ReconstructedStemID" %in% names(.post))) .post[, ReconstructedStemID := NA_integer_]
                if (!("ReconstructionMethod" %in% names(.post))) .post[, ReconstructionMethod := NA_character_]
                if (!("ConstraintViolation" %in% names(.post))) .post[, ConstraintViolation := NA]
                # NOTE: no !is.na(DBH) guard here.  Terminal NA-DBH post-
                # anchor rows with a known TrueStemID (rows anchored by
                # the pre-DP terminal-event propagation in main_cpp_bci.R
                # Step 2 / Step 3) must also be honoured to satisfy the
                # hard invariant.
                .post[!is.na(TrueStemID), `:=`(
                    ReconstructedStemID = as.integer(TrueStemID),
                    ReconstructionMethod = "given"
                )]
                .post[is.na(ReconstructionMethod), ReconstructionMethod := "none_after_anchor"]
                .post[, DP_FallbackReason := paste0("error:", substr(msg, 1, 200))]
                out <- data.table::rbindlist(list(out, .post), use.names = TRUE, fill = TRUE)
            }
            out
        }
    )

    # ---- TrueStemID HARD-INVARIANT backstop sweep (script-level) ----------
    # Final, idempotent enforcement of TrueStemID == ReconstructedStemID for
    # every row with a non-NA TrueStemID, regardless of which engine path
    # produced this output (DP success, probabilistic fallback, MF
    # re-insertion).  Mirrors the equivalent backstop in main_cpp_chunk.R
    # so both single-tag and chunked drivers offer the same guarantee.
    if (isTRUE(PIN_TRUESTEMID) && !is.null(out) &&
        all(c("TrueStemID", "ReconstructedStemID", "ReconstructionMethod") %in% names(out))) {
        .ts_rows <- which(!is.na(out$TrueStemID))
        if (length(.ts_rows) > 0L) {
            # ---- Script-level audit: detect engine-vs-pin disagreements ----
            # Last-line defense.  At this point both finalize_out and the
            # probabilistic matcher have already run their own audits, so any
            # override caught here means a row leaked through both inner
            # sweeps -- worth surfacing loudly.
            if (!("SweepAuditOverride" %in% names(out))) {
                out[, SweepAuditOverride := FALSE]
            } else {
                # Per-row backfill — mirrors dp_global_dp.R finalize_out.
                out[is.na(SweepAuditOverride), SweepAuditOverride := FALSE]
            }
            # Snapshot the engine's pre-sweep ReconstructedStemID.
            # Per-row backfill (see dp_global_dp.R::finalize_out): create the
            # column if absent, otherwise populate only NA cells so rows that
            # were appended after an inner sweep (post-anchor block, sub-
            # segment merges) also get their pre-sweep value captured.
            if (!("ReconstructedStemID_PreSweep" %in% names(out))) {
                out[, ReconstructedStemID_PreSweep := ReconstructedStemID]
            } else {
                out[
                    is.na(ReconstructedStemID_PreSweep),
                    ReconstructedStemID_PreSweep := ReconstructedStemID
                ]
            }
            .pre_recon <- out$ReconstructedStemID_PreSweep[.ts_rows]
            .true_int <- as.integer(out$TrueStemID[.ts_rows])
            .override_local <- !is.na(.pre_recon) & .pre_recon != .true_int
            if (any(.override_local)) {
                out[.ts_rows[.override_local], SweepAuditOverride := TRUE]
                log_msg(sprintf("[audit] script-level sweep overrode %d engine-assigned ReconstructedStemID value(s) for tag=%s (rows flagged via SweepAuditOverride=TRUE)", sum(.override_local), tag_label), "WARN")
            }

            # ---- Fix 2: duplicate-aware pinning -------------------------
            # Per-row decision: tentatively pin Recon := TrueStemID, but
            # only commit if doing so does NOT create a duplicate
            # ReconstructedStemID at the same (Tag, CensusID) given the
            # current state of the column (engine-assigned baseline +
            # any pins already committed earlier in this pass). If a pin
            # would collide, RESPECT THE ID GIVEN BY THE FULL ENGINE
            # RECONSTRUCTION — keep the PreSweep value, flag the row via
            # SweepRollbackToPreSweep=TRUE, and leave its
            # ReconstructionMethod untouched. SweepAuditOverride above
            # already records the disagreement; the rollback flag
            # records the chosen resolution.
            #
            # Rows are processed in stable index order so the outcome is
            # deterministic. Tag 258411 C6 is the canonical case: a
            # retag-campaign reuse of OriginalStemID 995110 would force
            # two rows at (Tag=258411, CensusID=6) onto Recon=995110;
            # the engine had already given the second row a fresh ID
            # (995113), which we now retain.
            if (!("SweepRollbackToPreSweep" %in% names(out))) {
                out[, SweepRollbackToPreSweep := FALSE]
            } else {
                out[is.na(SweepRollbackToPreSweep), SweepRollbackToPreSweep := FALSE]
            }

            .working_recon <- as.integer(out$ReconstructedStemID)
            .tag_vec <- out$Tag
            .cen_vec <- out$CensusID
            .true_vec_all <- as.integer(out$TrueStemID)

            # Pre-build (Tag, CensusID) -> row index map for fast peer lookup.
            .grp_key <- paste(.tag_vec, .cen_vec, sep = "\u0001")
            .grp_map <- split(seq_len(nrow(out)), .grp_key)

            .pin_idx <- integer(0)
            .rollback_idx <- integer(0)
            for (.r in .ts_rows) {
                .v <- .true_vec_all[.r]
                .peers <- .grp_map[[.grp_key[.r]]]
                .peers <- .peers[.peers != .r]
                if (length(.peers) > 0L) {
                    .peer_vals <- .working_recon[.peers]
                    if (any(!is.na(.peer_vals) & .peer_vals == .v)) {
                        # Collision: skip the pin, keep engine's value.
                        .rollback_idx <- c(.rollback_idx, .r)
                        next
                    }
                }
                .working_recon[.r] <- .v
                .pin_idx <- c(.pin_idx, .r)
            }

            if (length(.pin_idx) > 0L) {
                out[.pin_idx, `:=`(
                    ReconstructedStemID  = .true_vec_all[.pin_idx],
                    ReconstructionMethod = "given"
                )]
            }
            if (length(.rollback_idx) > 0L) {
                # Revert ReconstructedStemID back to the engine's pre-sweep
                # value (in case an inner finalize_out sweep already pinned
                # the row before script-level sweep runs). Leave
                # ReconstructionMethod untouched — it reflects the engine's
                # original assignment path.
                .pre_for_rollback <- out$ReconstructedStemID_PreSweep[.rollback_idx]
                out[.rollback_idx, `:=`(
                    ReconstructedStemID = .pre_for_rollback,
                    SweepRollbackToPreSweep = TRUE
                )]
                log_msg(sprintf(
                    "[audit] script-level sweep skipped %d TrueStemID pin(s) for tag=%s to avoid duplicate ReconstructedStemID at same (Tag,CensusID); engine's PreSweep id retained (rows flagged via SweepRollbackToPreSweep=TRUE)",
                    length(.rollback_idx), tag_label
                ), "WARN")
            }
        }
        # Unconditional materialisation of the audit columns: tags with NO
        # non-NA TrueStemID never enter the .ts_rows branch above, so the
        # columns would otherwise be missing from their output and downstream
        # consumers would have to special-case schema differences.  Create
        # them here with the no-override defaults if absent.
        if (!("SweepAuditOverride" %in% names(out))) {
            out[, SweepAuditOverride := FALSE]
        }
        if (!("ReconstructedStemID_PreSweep" %in% names(out))) {
            out[, ReconstructedStemID_PreSweep := ReconstructedStemID]
        }
        if (!("SweepRollbackToPreSweep" %in% names(out))) {
            out[, SweepRollbackToPreSweep := FALSE]
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
### 9) Main pipeline — run_main()
############################################################
# run_main() — loads data, estimates bio parameters, runs DP over all tags,
# writes CSV/RDS/PDF outputs, and returns list(xrun, bio_pars).
run_main <- function() {
    ensure_dir(out_dir)
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_started.txt"))
        },
        error = function(e) {
            message("[dp_global main_cpp.R] Warning writing run_started marker: ", conditionMessage(e))
        }
    )
    log_msg("Started run")

    # Create posteriors subdirectory (DP writes its files into <base>/posteriors)
    if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && !is.null(POSTERIOR_SAMPLES_PATH) && nzchar(POSTERIOR_SAMPLES_PATH)) {
        ensured <- file.path(POSTERIOR_SAMPLES_PATH, "posteriors")
        ensure_dir(ensured)
        log_msg(paste("Ensured posterior samples path:", ensured))
    }

    # 5.1 Load data
    xraw <- data.table::fread(INPUT_FILE)
    xraw <- ensure_species_column(xraw)
    xrun <- data.table::copy(xraw)

    # 5.2 Estimate biological parameters
    bio_pars <- list()

    for (sp in unique(xrun$species)) {
        bio_pars[[sp]] <- estimate_bio_pars(
            xrun[species == sp],
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
            # if measurement error, then lowest shrink is the min between the two for hard shrink guardrail
            #################
            # Extreme-growth guardrails (upper tail)
            # - growth_hard_prob is the *upper-tail* probability (e.g., 1e-4 means 99.99th percentile)
            # - growth_data_quantile is the empirical upper quantile used as a guardrail
            # - growth_soft_quantile sets a softer threshold used for a quadratic penalty
            # the highest (1 - growth_hard_prob) value of the probability function to get highest growth from measurement error
            growth_hard_prob = 1e-4,
            # Upper quantile for hard growth guardrail from empirical data
            growth_data_quantile = 0.999,
            # if measurement error, then highest growth is the max between the two for hard growth guardrail
            #################
            # to estimate the growth soft penalty k_growth - used if it becomes the minimum between max growth from measurement error or fixed or data
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
    }

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
    params$UNIQUE_TAGS <- if (exists("xrun")) unique(xrun$Tag) else NA
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
    # 5.4 DP meta settings
    dp_max_tracks_local <- if (is.null(DP_MAX_TRACKS)) auto_dp_max_tracks(xrun) else as.integer(DP_MAX_TRACKS)
    dp_max_tracks_local <- as.integer(dp_max_tracks_local)

    # 5.5 DP reconstruction
    out <- NULL
    if (isTRUE(RUN_DP)) {
        if (!isTRUE(RUN_ALL_TAGS)) {
            if (!(WHICH_TAG %in% unique(xrun$Tag))) {
                stop("Requested WHICH_TAG=", WHICH_TAG, " not found in data. Set WHICH_TAG to an existing Tag or enable RUN_ALL_TAGS=TRUE.")
            }
            dt_tag <- xrun[Tag == WHICH_TAG]
            if (!("DBH" %in% names(dt_tag)) || !("CensusID" %in% names(dt_tag)) || all(is.na(dt_tag$DBH)) || all(is.na(dt_tag$CensusID))) {
                log_msg(sprintf("Skipping WHICH_TAG=%s: all DBH or all CensusID missing for this Tag; returning rows as-is", WHICH_TAG), "WARN")
                out <- data.table::copy(dt_tag)
                out[, ReconstructionMethod := "skipped_no_data"]
            } else {
                out <- xrun[Tag == WHICH_TAG, run_dp_one_group(.SD, dp_max_tracks = dp_max_tracks_local), by = .(Tag, species)]
            }
        } else {
            if (!requireNamespace("parallel", quietly = TRUE)) {
                stop("Package not available: parallel. (It normally ships with R.)")
            }
            groups <- unique(xrun[, .(Tag, species)])
            data.table::setorder(groups, Tag, species)

            res_list <- parallel::mclapply(seq_len(nrow(groups)), function(i) {
                data.table::setDTthreads(1L)
                g <- groups[i]
                dtg <- xrun[Tag == g$Tag & species == g$species]
                # Return rows as-is for groups with missing DBH or CensusID (all NA)
                if (!("DBH" %in% names(dtg)) || !("CensusID" %in% names(dtg)) || all(is.na(dtg$DBH)) || all(is.na(dtg$CensusID))) {
                    log_msg(sprintf("Skipping Tag=%s species=%s: all DBH or all CensusID missing; returning rows as-is", g$Tag, g$species), "WARN")
                    dtg[, ReconstructionMethod := "skipped_no_data"]
                    return(dtg)
                }
                run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local)
            }, mc.cores = MC_CORES)

            # Remove NULLs (skipped groups) before binding
            res_list <- Filter(Negate(is.null), res_list)
            if (length(res_list) == 0L) {
                out <- NULL
            } else {
                out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
            }
        }
    }
    out <- maybe_add_posterior_bins(out)

    # 5.5b Post-engine backfill of orphan terminal-event rows.
    #      For rows with Status %in% {"dead","stem dead","broken below"} AND
    #      NA DBH that the engine could not match (ReconstructedStemID NA),
    #      copy the ReconstructedStemID from the most recent prior row in the
    #      same (Tag, OriginalStemID) group (LOCF). Sets
    #      ReconstructionMethod = "carried_terminal" on filled rows.
    #      Shared helper defined in dp_global/R/dp_global_main.R; mirrors
    #      Step 9b in main_cpp_bci.R and the per-chunk call in main_cpp_chunk.R.
    out <- apply_carried_terminal_backfill(out)

    # 5.5c Born-orphan stem backfill. Rows with NA Recon, NA TrueStemID,
    #      NA DBH, and a non-NA source-id (StemID / OriginalStemID) are stems
    #      whose first appearance had no measurement (e.g. born broken-below
    #      at C7+); the engine cannot reach them. Fill Recon = source-id and
    #      tag ReconstructionMethod = "given_orphan". Mirrors Step 9c in
    #      main_cpp_bci.R and the per-chunk call in main_cpp_chunk.R.
    out <- apply_orphan_stem_backfill(out)

    # 5.5d Broken-below invariant pass. Enforce R1 (split-on-break) and R2
    #      (terminate-on-stump) per `apply_broken_below_invariants` in
    #      dp_global/R/dp_global_main.R. Mirrors Step 9d in main_cpp_bci.R
    #      and the per-chunk call in main_cpp_chunk.R.
    out <- apply_broken_below_invariants(out)

    # 5.5e Reverse the direction of ReconstructedStemID numbering so
    #      smaller integers correspond to earlier stem appearances.
    #      Must run AFTER apply_broken_below_invariants(). Writes per-tag
    #      companion mapping files alongside posterior path files.
    #      See dp_global/improvements.md for the full algorithm.
    if (!is.null(out)) {
        .renum <- renumber_engine_minted_ids(
            out,
            posterior_top_k = DP_POSTERIOR_TOP_K,
            posterior_samples_path = out_dir
        )
        out <- .renum$out
        # Finalize posterior path files in the renumbered ID space
        # (recommended architecture, see dp_global/improvements.md).
        finalize_posterior_paths(
            out,
            posterior_samples_path = out_dir,
            mapping = .renum$mapping
        )
    }

    # Record run output directory (basename) in each row to avoid variable/column name collision
    if (!is.null(out)) {
        out[, run_out_dir := basename(out_dir)]
    }

    # 5.6 Optional: write DP outputs (centralized helpers)
    maybe_write(
        isTRUE(WRITE_DP_CSV) && !is.null(out), out_path("dp_csv"),
        function() data.table::fwrite(out, out_path("dp_csv")),
        "DP CSV"
    )

    maybe_write(
        isTRUE(WRITE_DP_RDS) && !is.null(out), out_path("dp_rds"),
        function() saveRDS(out, file = out_path("dp_rds")),
        "DP RDS"
    )

    maybe_write(
        isTRUE(WRITE_DP_FEATHER) && !is.null(out), out_path("dp_feather"),
        function() {
            if (!requireNamespace("arrow", quietly = TRUE)) stop("'arrow' package not available; skipping feather output")
            arrow::write_feather(out, out_path("dp_feather"))
        },
        "DP Feather"
    )

    maybe_write(
        isTRUE(WRITE_DP_PDF) && !is.null(out), out_path("dp_pdf"),
        function() {
            plot_tag_to_pdf(
                out,
                pdf_file = out_path("dp_pdf"),
                include_reference = DP_PDF_INCLUDE_REFERENCE,
                tag = if (isTRUE(PLOT_PDF_ONE_TAG_ONLY)) WHICH_TAG else NULL
            )
        },
        "DP PDF"
    )

    # 5.7 Optional: realism report
    if (isTRUE(RUN_REALISM_REPORT) && !is.null(out)) {
        sp0 <- unique(out$species)
        sp0 <- sp0[!is.na(sp0) & nzchar(sp0)]
        sp0 <- if (length(sp0) > 0L) sp0[[1L]] else FORCED_SPECIES_LABEL

        base_args0 <- bio_pars_to_transition_args(bio_pars[[sp0]])
        #* TODO: Define interval_years per row if needed
        rep0 <- realism_report_from_reconstruction(out, interval_years = 5, base_args = base_args0)

        # Add out_dir to realism outputs
        rep0$summary[, out_dir := basename(out_dir)]
        rep0$by_group[, out_dir := basename(out_dir)]
        rep0$suggestions[, out_dir := basename(out_dir)]

        # Write realism outputs via maybe_write for consistent logs and directory creation
        sum_path <- file.path(out_dir, paste0("tag_", WHICH_TAG, "_realism_summary_rcpp.csv"))
        by_path <- file.path(out_dir, paste0("tag_", WHICH_TAG, "_realism_by_tag_rcpp.csv"))
        sug_path <- file.path(out_dir, paste0("tag_", WHICH_TAG, "_realism_tuning_suggestions_rcpp.csv"))

        maybe_write(isTRUE(RUN_REALISM_REPORT), sum_path, function() {
            data.table::fwrite(rep0$summary, file = sum_path)
        }, "Realism summary")

        maybe_write(isTRUE(RUN_REALISM_REPORT), by_path, function() {
            data.table::fwrite(rep0$by_group, file = by_path)
        }, "Realism by-tag")

        maybe_write(isTRUE(RUN_REALISM_REPORT), sug_path, function() {
            data.table::fwrite(rep0$suggestions, file = sug_path)
        }, "Realism tuning suggestions")
    }

    # 5.8 Optional: sensitivity sweeps
    if (isTRUE(RUN_SENSITIVITY)) {
        sp_sens <- unique(xrun$species)
        sp_sens <- sp_sens[!is.na(sp_sens) & nzchar(sp_sens)]
        sp_sens <- if (length(sp_sens) > 0L) sp_sens[[1L]] else FORCED_SPECIES_LABEL

        base <- bio_pars_to_transition_args(bio_pars[[sp_sens]])
        #* TODO: Define interval_years per row if needed
        sc <- make_demo_scenarios(base, interval_years = 5)
        param_grids <- default_param_grids(base, n = 200)

        all_sweeps <- build_all_sweeps(
            scenarios = sc,
            #* TODO: Define interval_years per row if needed
            interval_years = 5,
            base_args = base,
            grids = param_grids,
            abs_jump = 1000
        )

        dts <- all_sweeps$dts
        dt_all <- all_sweeps$all
        dt_jumps <- all_sweeps$jumps

        # Add out_dir to sensitivity outputs
        dt_jumps[, out_dir := basename(out_dir)]
        dt_all[, out_dir := basename(out_dir)]
        for (key in names(dts)) {
            dts[[key]][, out_dir := basename(out_dir)]
        }
        # Update all_sweeps with modified components
        all_sweeps$dts <- dts
        all_sweeps$all <- dt_all

        # example_key <- "growth_ok__sigma1"
        # if (example_key %in% names(dts)) {
        #     print(plot_sweep_components(dts[[example_key]]), yscale = "delta")
        # }

        if (isTRUE(WRITE_OUTPUTS)) {
            sweeps_rds <- file.path(out_dir, "simulated_all_transition_cost_sweeps_rcpp.rds")
            sweeps_csv <- file.path(out_dir, "simulated_all_transition_cost_sweeps.csv")
            jumps_csv <- file.path(out_dir, "simulated_all_transition_cost_sweep_jumps_rcpp.csv")
            jumps_rds <- file.path(out_dir, "simulated_all_transition_cost_jumps_rcpp.rds")
            jumps_csv2 <- file.path(out_dir, "simulated_all_transition_cost_jumps_rcpp.csv")

            maybe_write(isTRUE(WRITE_OUTPUTS), sweeps_rds, function() {
                saveRDS(all_sweeps, file = sweeps_rds)
            }, "Sensitivity sweeps (RDS)")

            # Optionally write the combined sweep table
            maybe_write(isTRUE(WRITE_OUTPUTS), sweeps_csv, function() {
                data.table::fwrite(dt_all, file = sweeps_csv)
            }, "Sensitivity sweeps (CSV)")

            maybe_write(isTRUE(WRITE_OUTPUTS), jumps_csv, function() {
                data.table::fwrite(dt_jumps, file = jumps_csv)
            }, "Sensitivity jump summaries (CSV)")

            maybe_write(isTRUE(WRITE_OUTPUTS), jumps_rds, function() {
                saveRDS(dt_jumps, file = jumps_rds)
            }, "Sensitivity jump summaries (RDS)")

            maybe_write(isTRUE(WRITE_OUTPUTS), jumps_csv2, function() {
                data.table::fwrite(dt_jumps, file = jumps_csv2)
            }, "Sensitivity jump summaries (CSV alt)")
        }

        maybe_write(isTRUE(MAKE_ALL_SWEEPS_PDF), file.path(out_dir, "simulated_all_transition_cost_sweeps_rcpp.pdf"), function() {
            plot_all_sweeps_to_pdf(
                all_sweeps,
                pdf_file = file.path(out_dir, "simulated_all_transition_cost_sweeps_rcpp.pdf"),
                y_scale = "delta",
                subtitle = basename(out_dir)
            )
        }, "Sensitivity sweeps PDF")
    }

    # Write a small finished marker so users and wrappers can detect job completion
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt"))
            log_msg("Finished run")
        },
        error = function(e) {
            log_msg(sprintf("Warning writing run_finished marker: %s", conditionMessage(e)), "WARN")
        }
    )

    invisible(list(out = out, xrun = xrun, bio_pars = bio_pars))
}

############################################################
### 10) Entrypoint
############################################################
# Execute run_main() when called via Rscript; skip when sourced interactively.

# When you run this file with Rscript, sys.nframe()==0 and we execute.
# When you source() this file from another script/session, we only define helpers.
if (sys.nframe() == 0L) {
    res <- run_main()
}

############################################################
### Optional demo: k sweep join-vs-split analysis
############################################################
# Usage (run after run_main() so that `bio_pars` exists):

if (sys.nframe() == 0L && isTRUE(RUN_K_SWEEP_DEMO)) {
    xrun <- res$xrun
    bio_pars <- res$bio_pars
    sp0 <- if (isTRUE(FORCE_ONE_SPECIES_PARAMETERS)) {
        FORCED_SPECIES_LABEL
    } else {
        # Prefer keys from `bio_pars` (guaranteed to be valid list indices)
        nms <- names(bio_pars)
        if (length(nms) > 0L && nzchar(nms[[1L]])) nms[[1L]] else unique(as.character(xrun$species))[[1L]]
    }
    dt_sweep <- k_sweep_join_vs_split(
        scenarios = data.frame(
            d0 = c(20, 20, 40, 40, 10, 10, 50, 30, 30),
            d1 = c(18, 10, 45, 70, 9, 30, 55, 32, 28),
            label = c(
                "small shrink", "big shrink", "normal growth",
                "extreme growth", "small tree shrink", "small tree extreme growth",
                "large tree growth", "moderate growth", "moderate shrink"
            )
        ),
        #* TODO: Define interval_years per row if needed
        interval_years = 5,
        bio = bio_pars[[sp0]],
        temperature = 1L,
        which_k = "auto",
        # optional: include candidate-pruning bounds from the DP enumerator
        prune_min_annual_growth = MAX_SHRINK_FIXED,
        prune_max_annual_growth = MAX_GROWTH_FIXED,
        subtitle = basename(out_dir)
    )

    print(k_sweep_crosspoints(dt_sweep))

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        pp <- plot_k_sweep_join_vs_split(
            dt_sweep,
            k_max = 1000,
            out_path = here(out_dir, "k_sweep_join_vs_split_demo_rcpp.pdf"),
            subtitle = basename(out_dir)
        )
    } else {
        message("ggplot2 not available; install ggplot2 to visualize k sweeps.")
    }
}

# Export bio-parameter report if `res` and `res$bio_pars` are available.
if (exists("res") && !is.null(res) && !is.null(res$bio_pars)) {
    source(here("dp_global", "R", "check_functions.r"))
    tryCatch(
        {
            export_bio_pars_report(res$bio_pars,
                species = NULL,
                interval_years = 5,
                out_file = file.path(out_dir, "bio_pars_report.pdf")
            )
        },
        error = function(e) {
            message("[dp_global main_cpp.R] Warning exporting bio_pars_report: ", conditionMessage(e))
        }
    )
} else {
    message("[dp_global main_cpp.R] Skipping bio_pars report: 'res$bio_pars' not available (chunked run or earlier error).")
}

############################################################
### Usage guide
############################################################
#
# PURPOSE
#   Single-pass DP pipeline. Loads the full data set into memory, runs the
#   dynamic-programming stem reconstruction over all (or one) Tag x species
#   groups, and writes a combined CSV / RDS / PDF output. Also supports
#   optional sensitivity sweeps, realism reports, and k-sweep demos.
#   Use main_cpp_chunk.R instead when memory is limited.
#
# BASIC RUN — all tags, default settings
#   Rscript dp_global/scripts/main_cpp.R \
#     --INPUT_FILE=data/my_stems.csv
#
# SINGLE-TAG MODE
#   Rscript dp_global/scripts/main_cpp.R \
#     --INPUT_FILE=data/my_stems.csv \
#     --RUN_ALL_TAGS=FALSE \
#     --WHICH_TAG=2747
#
# COMMON OVERRIDES
#   --DP_MAX_STATES=1100          Max DP states per track (higher = slower, more accurate)
#   --MANUAL_CORES=TRUE           Use a fixed core count instead of auto-detect
#   --MANUAL_CORES_VALUE=8        Number of cores (requires MANUAL_CORES=TRUE)
#   --WRITE_DP_PDF=TRUE/FALSE     Whether to produce the per-tag PDF plot
#   --POSTERIOR_SAMPLES=250       Draw N posterior samples per group (0 = disabled)
#   --POSTERIOR_SAMPLES_FORMAT=csv|rds|feather
#   --USE_MEASUREMENT_ERROR=TRUE  Enable measurement-error model for bio params
#   --DP_FALLBACK_GROWTH_FORMS="fig,tree"
#                                 Comma-separated list of growth forms that bypass
#                                 species-specific bio params
#   --SENSITIVITY_MODE=run+write+pdf
#                                 Run transition-cost sensitivity sweep and write PDF
#   --RUN_REALISM_REPORT=TRUE     Write per-tag realism / tuning-suggestion tables
#   --RUN_K_SWEEP_DEMO=TRUE       Produce k-sweep join-vs-split plot after the run
#
# STOPPING & RESTARTING
#   This script processes all groups in a single call and does not support
#   incremental checkpointing. If stopped mid-run, restart from scratch.
#   For large data sets where checkpointing is needed, use main_cpp_chunk.R.
#
# HELP
#   Rscript dp_global/scripts/main_cpp.R --help

# Rscript dp_global/scripts/main_cpp.R \
#     --INPUT_FILE=/Users/medinaja/GDrive_Science/STRI/STEM_IDENTIFICATION_TEST/data_simulation/data/simulated_data_1.csv \
#     --RUN_ALL_TAGS=FALSE \
#     --DP_MAX_STATES=2 \
#     --WHICH_TAG=20
