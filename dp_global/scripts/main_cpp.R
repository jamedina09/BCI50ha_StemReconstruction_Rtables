############################################################
### main.R — dp_global driver
############################################################
# Goal
#   One place to run the DP_GLOBAL workflow end-to-end

############################################################
### 0) Housekeeping
############################################################
# Avoid nuking the user's interactive environment when sourcing this file.
# Only clear the workspace when running as a top-level script.
if (sys.nframe() == 0L) {
    rm(list = ls())
}

############################################################
### Command-line argument parsing
############################################################
# Parse command-line arguments to override defaults.
# Usage: Rscript main_cpp.R --which_tag=2 --RUN_ALL_TAGS=TRUE --DP_MODE=marginals+bins --SENSITIVITY_MODE=run --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=8
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
                # Try to convert to appropriate type
                if (tolower(val) %in% c("true", "false")) {
                    val <- as.logical(tolower(val))
                } else if (grepl("^[0-9]+$", val)) {
                    val <- as.integer(val)
                } else if (grepl("^[0-9]+\\.[0-9]+$", val)) {
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
### 1) Dependencies
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
### Manual overrides (uncomment and edit for interactive/manual runs)
############################################################
## 2.1 Input data and species handling
input_file <- here("data_simulation", "data", "simulated_data_1.csv")
FORCE_ONE_SPECIES_PARAMETERS <- TRUE
if (isTRUE(FORCE_ONE_SPECIES_PARAMETERS)) {
    FORCED_SPECIES_LABEL <- "all"
    message("[dp_global main.R] FORCE_ONE_SPECIES_PARAMETERS=TRUE: using single species label '", FORCED_SPECIES_LABEL, "' for all trees.")
} else {
    message("[dp_global main.R] FORCE_ONE_SPECIES_PARAMETERS=FALSE: using species column from data for parameter estimation.")
}
SPECIES_COL <- NULL

############################################################
### 2.2 Parameter estimation settings
############################################################
# All settings related to parameter estimation and biological realism
USE_MEASUREMENT_ERROR <- TRUE
MAX_GROWTH_HARD_SOURCE <- "data"
MAX_GROWTH_FIXED <- 7.5
MAX_SHRINK_HARD_SOURCE <- "data"
MAX_SHRINK_FIXED <- -0.5
K_SHRINK_SOURCE <- "data"
K_SHRINK_FIXED <- 0 # 0 to disable soft penalty
K_GROWTH_SOURCE <- "data"
K_GROWTH_FIXED <- 0 # 0 to disable soft penalty
RECRUIT_MAX_SOURCE <- "fixed"
RECRUIT_MAX_FIXED <- 5

############################################################
### 2.3 DP running settings
############################################################
DP_MODE <- "marginals+bins" # Options: "none", "marginals", "marginals+bins"
which_tag <- 20L # 20 to compare outputs
anchor_start_census <- 7L
DP_VERBOSE <- TRUE
DP_POSTERIOR_TOP_K <- 2L
dp_max_tracks <- NULL # auto (computed from data)
dp_max_states <- 40000L
dp_slack_tracks <- 1L
# NOTE: Optionally require that slack be granted only if an anchor DBH is recruitable
# (i.e., DBH <= Bio_Recruit_MaxDBH_unit + eps). Set FALSE to preserve current behavior.
dp_slack_require_anchor_recruitable <- TRUE
# Tolerance (cm) used when comparing anchor DBH to recruit_max_dbh
dp_slack_require_anchor_eps <- 1e-6

# Posterior sampling defaults (disabled by default)
POSTERIOR_SAMPLES <- 0L
POSTERIOR_SAMPLES_FORMAT <- "rds" # options: 'rds', 'feather', 'csv'
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- NULL

############################################################
### 2.4 Parallel and output settings
############################################################
RUN_ALL_TAGS <- FALSE
MANUAL_CORES <- TRUE # Flag to manually define cores instead of auto-detecting
MANUAL_CORES_VALUE <- 1L # Number of cores to use if MANUAL_CORES=TRUE

############################################################
### 2.5 CPP
############################################################

## create output directory within project
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

# Encode numeric values for directory-safe names
# -0.5 -> m0p5, 7.5 -> 7p5
encode_num <- function(x) {
    if (is.null(x) || is.na(x)) {
        return("NA")
    }
    s <- as.character(x)
    s <- gsub("-", "m", s)
    s <- gsub("\\.", "p", s)
    s
}

build_out_dir_name <- function() {
    # # Use explicit OUT_DIR_NAME if set
    # if (!is.null(OUT_DIR_NAME) && nzchar(OUT_DIR_NAME)) {
    #     return(OUT_DIR_NAME)
    # }

    # Timestamp: use BATCH_TS if provided; else fallback to current date+time
    ts <- if (exists("BATCH_TS") && nzchar(BATCH_TS)) BATCH_TS else format(Sys.time(), "%Y%m%d_%H%M%S")

    # Config name (for output directory label)
    config_part <- if (exists("CONFIG_NAME") && !is.null(CONFIG_NAME)) {
        CONFIG_NAME
    } else {
        "unknown"
    }

    # Tag info
    tag_part <- if (isTRUE(RUN_ALL_TAGS)) {
        "allT"
    } else {
        paste0("T", as.integer(which_tag))
    }

    # DP mode label
    dp_part <- switch(DP_MODE,
        "none" = "NO_DP",
        "map" = "DP_S",
        "marginals" = "DP_M",
        "marginals+bins" = "DP_MB",
        "DP_U"
    )

    # Measurement error label
    me_part <- if (isTRUE(USE_MEASUREMENT_ERROR)) "ME" else "NME"

    # Encode numeric values for directory-safe names
    encode_num <- function(x) {
        if (is.null(x) || is.na(x)) {
            return("NA")
        }
        s <- as.character(x)
        s <- gsub("-", "m", s)
        s <- gsub("\\.", "p", s)
        s
    }

    max_growth_hard_ <- switch(MAX_GROWTH_HARD_SOURCE,
        "fixed" = paste0("g", encode_num(MAX_GROWTH_FIXED)),
        "data"  = "gD",
        "gU"
    )

    max_shrink_hard_ <- switch(MAX_SHRINK_HARD_SOURCE,
        "fixed" = paste0("s", encode_num(MAX_SHRINK_FIXED)),
        "data"  = "sD",
        "sU"
    )

    soft_growth_ <- switch(K_GROWTH_SOURCE,
        "fixed" = paste0("kg", encode_num(K_GROWTH_FIXED)),
        "data"  = "kgD",
        "kgU"
    )

    soft_shrink_ <- switch(K_SHRINK_SOURCE,
        "fixed" = paste0("ks", encode_num(K_SHRINK_FIXED)),
        "data"  = "ksD",
        "ksU"
    )

    # Assemble final directory name
    dir_name <- paste(
        ts,
        config_part,
        tag_part,
        paste0(dp_part, "_", me_part),
        max_growth_hard_,
        max_shrink_hard_,
        soft_growth_,
        soft_shrink_,
        "rcpp",
        sep = "_"
    )

    return(dir_name)
}

WRITE_DP_CSV <- TRUE
WRITE_DP_RDS <- TRUE
WRITE_DP_PDF <- TRUE
DP_PDF_INCLUDE_REFERENCE <- TRUE

if (!isTRUE(RUN_ALL_TAGS)) {
    PLOT_PDF_ONE_TAG_ONLY <- TRUE
} else {
    PLOT_PDF_ONE_TAG_ONLY <- FALSE
}

# Default project root so --PROJECT_ROOT=/path overrides are accepted by the CLI parser
PROJECT_ROOT <- here::here()
############################################################
### 2.5 Sensitivity analysis settings
############################################################
SENSITIVITY_MODE <- "none" # Options: "none", "run", "run+write", "run+write+pdf"
RUN_K_SWEEP_DEMO <- FALSE

############################################################
### 2.6 Realism report settings
############################################################
RUN_REALISM_REPORT <- FALSE

print_help <- function() {
    cat("Usage: Rscript scripts/main_cpp.R [--KEY=VALUE] [--FLAG]\n")
    cat("Common keys and defaults:\n")
    keys <- c(
        "input_file", "FORCE_ONE_SPECIES_PARAMETERS", "DP_MODE", "which_tag",
        "anchor_start_census", "DP_VERBOSE", "RUN_ALL_TAGS",
        "MANUAL_CORES", "MANUAL_CORES_VALUE", "WRITE_DP_CSV", "WRITE_DP_RDS", "WRITE_DP_PDF",
        "dp_max_states", "dp_slack_tracks", "dp_slack_require_anchor_recruitable", "dp_slack_require_anchor_eps", "POSTERIOR_SAMPLES", "POSTERIOR_SAMPLES_FORMAT", "POSTERIOR_SAMPLES_PATH", "POSTERIOR_SAMPLE_SEED", "SENSITIVITY_MODE", "RUN_K_SWEEP_DEMO", "RUN_REALISM_REPORT", "PROJECT_ROOT",
        "BATCH_TS", "CONFIG_NAME", "USE_MEASUREMENT_ERROR"
    )
    for (k in keys) {
        val <- if (exists(k, inherits = FALSE)) get(k) else "<not set>"
        cat(sprintf("  --%s = %s\n", k, as.character(val)))
    }
    cat("\nFlags without =value are treated as boolean TRUE (e.g., --DRY_RUN).\n")
}

# If user asked for help, print and exit (do this before applying overrides)
if (isTRUE(overrides$help) || isTRUE(overrides$h)) {
    print_help()
    quit(save = "no", status = 0)
}

# Apply command-line overrides with validation and warnings for unknown keys
for (name in names(overrides)) {
    if (name %in% c("help", "h")) next
    if (!exists(name, inherits = FALSE)) {
        warning(sprintf("[dp_global main_cpp.R] Unknown override '%s' (ignored).\n", name))
        next
    }
    assign(name, overrides[[name]])
    message("[dp_global main_cpp.R] Overriding ", name, " = ", overrides[[name]])
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
### Input file validation (must be explicit)
############################################################

if (!exists("input_file") || is.null(input_file) || !nzchar(input_file)) {
    stop(
        "No input_file specified. ",
        "Provide --input_file=... on the command line or define input_file explicitly."
    )
}

if (!file.exists(input_file)) {
    stop("input_file does not exist: ", input_file)
}

# Derive booleans from modes
RUN_DP <- DP_MODE != "none"
ADD_DP_POSTERIOR_BINS <- DP_MODE == "marginals+bins"
RUN_SENSITIVITY <- SENSITIVITY_MODE != "none"
WRITE_OUTPUTS <- SENSITIVITY_MODE %in% c("run+write", "run+write+pdf")
MAKE_ALL_SWEEPS_PDF <- SENSITIVITY_MODE == "run+write+pdf"

# Final output directory for this run (created at runtime in run_main())
out_dir <- file.path(base_out_dir, build_out_dir_name())
message("[dp_global main_cpp.R] out_dir (computed): ", out_dir)
message("[dp_global main_cpp.R] getwd(): ", getwd())

DP_CSV_FILE <- file.path(out_dir, "stem_reconstruction_dp_global_rcpp.csv")
DP_RDS_FILE <- file.path(out_dir, "stem_reconstruction_dp_global_rcpp.rds")
DP_PDF_FILE <- file.path(out_dir, "stem_reconstruction_dp_global_rcpp.pdf")

############################################################
### 3) Source project code
############################################################
source(here("dp_global", "R", "dp_global_main.R"))
source(here("dp_global", "R", "sensitivity_transition_cost_bio.R"))
source(here("dp_global", "R", "realism_calibration.R"))
source(here("dp_global", "R", "k_tuning_viz.R"))

############################################################
### 4) Helpers
############################################################

############################################################
### 4.1) Optional tuning / inspection helpers
############################################################

# Soft penalties vs hard guardrails (unit reminder)
# - Soft penalties operate on DBH differences (cm) over the interval.
# - Hard guardrails operate on annualized growth (cm/year).
#
# Choosing k from a reference excess:
# - Soft penalty is quadratic: soft_cost = k * (delta_cm^2)
# - If you want delta_cm = D to contribute cost C, set k = C / (D^2).

# Helper: compute the *actual* soft-penalty cost for a given k and delta.
# - delta_cm is in cm over the interval (NOT cm/year).
# - temperature is the marginal-DP temperature; weight multiplier is exp(-soft_cost / temperature).
# soft_cost_from_k <- function(delta_cm, k, temperature = 1) {
#     delta_cm <- as.numeric(delta_cm)
#     k <- as.numeric(k)
#     temperature <- as.numeric(temperature)
#     cost <- k * (delta_cm^2)
#     data.frame(
#         delta_cm = delta_cm,
#         k = k,
#         soft_cost = cost,
#         temperature = temperature,
#         weight_multiplier = exp(-cost / temperature)
#     )
# }

# soft_cost_from_k(delta_cm = seq(-10, 10, by = 1), k = 20, temperature = 1)

ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        message("[dp_global main_cpp.R] Created directory: ", path)
    } else {
        message("[dp_global main_cpp.R] Directory already exists: ", path)
    }
    invisible(path)
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
            message("[DP_GLOBAL main.R] Using '", species_col, "' as species column. Set SPECIES_COL to override.")
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
        CensusID <= anchor_start_census & !is.na(DBH),
        .N,
        by = .(Tag, CensusID)
    ][, max(N, na.rm = TRUE)]
    if (!is.finite(max_obs_any_tag_census)) max_obs_any_tag_census <- 0L
    as.integer(max_obs_any_tag_census + 1L)
}

run_dp_one_group <- function(dtg, dp_max_tracks) {
    match_stems_dp_global_backward_marginals_batch(
        tree_data = data.table::copy(dtg),
        min_growth = MAX_SHRINK_FIXED,
        max_growth = MAX_GROWTH_FIXED,
        anchor_start = anchor_start_census,
        max_tracks = dp_max_tracks,
        max_states = dp_max_states,
        slack_tracks = dp_slack_tracks,
        slack_require_anchor_recruitable = dp_slack_require_anchor_recruitable,
        slack_require_anchor_eps = dp_slack_require_anchor_eps,
        temperature = 1,
        posterior_top_k = DP_POSTERIOR_TOP_K,
        # posterior sampling controls (disabled by default)
        posterior_samples = get0("POSTERIOR_SAMPLES", ifnotfound = 0L),
        posterior_samples_format = get0("POSTERIOR_SAMPLES_FORMAT", ifnotfound = "rds"),
        posterior_samples_path = get0("POSTERIOR_SAMPLES_PATH", ifnotfound = NULL),
        posterior_sample_seed = get0("POSTERIOR_SAMPLE_SEED", ifnotfound = NULL),
        use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
        verbose = isTRUE(DP_VERBOSE)
    )
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
### 5) Main pipeline
############################################################
run_main <- function() {
    ensure_dir(out_dir)

    # Write a small startup marker so parallel runs can be observed immediately
    # (helps verify jobs start concurrently before heavy computation)
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_started.txt"))
        },
        error = function(e) {
            message("[dp_global main_cpp.R] Warning writing run_started marker: ", conditionMessage(e))
        }
    )

    # 5.1 Load data
    xraw <- data.table::fread(input_file)
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
            recruit_max_source = get0("RECRUIT_MAX_SOURCE", ifnotfound = "data"),
            recruit_max_fixed = as.numeric(get0("RECRUIT_MAX_FIXED", ifnotfound = 5))
        )
    }

    # Write a small text file recording the parameters used to build the
    # run-specific output directory name so runs are reproducible.
    # Gather all parameters in the environment for full run context

    params <- list()
    params$TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    params$OUT_DIR <- out_dir
    params$GENERATED_DIR_NAME <- basename(out_dir)

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
    dp_max_tracks_local <- if (is.null(dp_max_tracks)) auto_dp_max_tracks(xrun) else as.integer(dp_max_tracks)
    dp_max_tracks_local <- as.integer(dp_max_tracks_local)

    # 5.5 DP reconstruction
    out <- NULL
    if (isTRUE(RUN_DP)) {
        if (!isTRUE(RUN_ALL_TAGS)) {
            if (!(which_tag %in% unique(xrun$Tag))) {
                stop("Requested which_tag=", which_tag, " not found in data. Set which_tag to an existing Tag or enable RUN_ALL_TAGS=TRUE.")
            }
            out <- xrun[Tag == which_tag, run_dp_one_group(.SD, dp_max_tracks = dp_max_tracks_local), by = .(Tag, species)]
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
                run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local)
            }, mc.cores = MC_CORES)

            out <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
        }
    }
    out <- maybe_add_posterior_bins(out)

    # Add output directory name as a column for reference
    if (!is.null(out)) {
        out[, out_dir := basename(out_dir)]
    }

    # 5.6 Optional: write DP outputs
    if (isTRUE(WRITE_DP_CSV) && !is.null(out)) {
        data.table::fwrite(out, file = DP_CSV_FILE)
    }
    if (isTRUE(WRITE_DP_RDS) && !is.null(out)) {
        saveRDS(out, file = DP_RDS_FILE)
    }
    if (isTRUE(WRITE_DP_PDF) && !is.null(out)) {
        plot_tag_to_pdf(
            out,
            pdf_file = DP_PDF_FILE,
            include_reference = DP_PDF_INCLUDE_REFERENCE,
            tag = if (isTRUE(PLOT_PDF_ONE_TAG_ONLY)) which_tag else NULL
        )
    }

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

        data.table::fwrite(rep0$summary, file = file.path(out_dir, paste0("tag_", which_tag, "_realism_summary_rcpp.csv")))
        data.table::fwrite(rep0$by_group, file = file.path(out_dir, paste0("tag_", which_tag, "_realism_by_tag_rcpp.csv")))
        data.table::fwrite(rep0$suggestions, file = file.path(out_dir, paste0("tag_", which_tag, "_realism_tuning_suggestions_rcpp.csv")))
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
            saveRDS(all_sweeps, file = file.path(out_dir, "simulated_all_transition_cost_sweeps_rcpp.rds"))
            # data.table::fwrite(dt_all, file = file.path(out_dir, "simulated_all_transition_cost_sweeps.csv"))
            data.table::fwrite(dt_jumps, file = file.path(out_dir, "simulated_all_transition_cost_sweep_jumps_rcpp.csv"))
            saveRDS(dt_jumps, file = file.path(out_dir, "simulated_all_transition_cost_jumps_rcpp.rds"))
            data.table::fwrite(dt_jumps, file = file.path(out_dir, "simulated_all_transition_cost_jumps_rcpp.csv"))
        }

        if (isTRUE(MAKE_ALL_SWEEPS_PDF)) {
            plot_all_sweeps_to_pdf(
                all_sweeps,
                pdf_file = file.path(out_dir, "simulated_all_transition_cost_sweeps_rcpp.pdf"),
                y_scale = "delta",
                subtitle = basename(out_dir)
            )
        }
    }

    # Write a small finished marker so users and wrappers can detect job completion
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt"))
        },
        error = function(e) {
            message("[dp_global main_cpp.R] Warning writing run_finished marker: ", conditionMessage(e))
        }
    )

    invisible(list(out = out, xrun = xrun, bio_pars = bio_pars))
}

############################################################
### 6) Script entrypoint
############################################################

# When you run this file with Rscript, sys.nframe()==0 and we execute.
# When you source() this file from another script/session, we only define helpers.
if (sys.nframe() == 0L) {
    res <- run_main()
}

############################################################
### EOptional demo: k sweep join-vs-split analysis
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

source(here("dp_global", "R", "check_functions.r"))
export_bio_pars_report(res$bio_pars,
    species = NULL,
    interval_years = 5,
    out_file = file.path(out_dir, "bio_pars_report.pdf")
)
