# FIXME - RUN PROFILE TO dp_global/scripts/main_cpp_chunk.R
########################################################################
### main_cpp_chunk.R — dp_global driver (standalone chunked version)
###
### Standalone chunked driver that does NOT source `main_cpp.R`.
### It uses the same helpers in `dp_global/R/` but contains its own
### CLI parsing, defaults, and chunk-only run loop.
########################################################################
############################################################
### 0) Housekeeping
############################################################
# Avoid nuking the user's interactive environment when sourcing this file.
# Only clear the workspace when running as a top-level script.
if (sys.nframe() == 0L) {
    rm(list = ls())
}

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

# Small helper to create directories safely
ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        message("[main_cpp_chunk.R] Created directory: ", path)
    }
    invisible(path)
}

# Simple logging helper (prints and appends to run_log.txt in out_dir when available)
log_msg <- function(msg, level = "INFO") {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full <- sprintf("[%s] %s: %s", ts, level, msg)
    message(full)
    # Append to run log if out_dir exists
    if (exists("out_dir") && nzchar(out_dir)) {
        tryCatch(
            write(full, file = file.path(out_dir, "run_log.txt"), append = TRUE),
            error = function(e) NULL
        )
    }
} 

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
MAX_GROWTH_HARD_SOURCE <- "fixed"
MAX_GROWTH_FIXED <- 7.5
MAX_SHRINK_HARD_SOURCE <- "fixed"
MAX_SHRINK_FIXED <- -0.5
K_SHRINK_SOURCE <- "fixed"
K_SHRINK_FIXED <- 0 # 0 to disable soft penalty
K_GROWTH_SOURCE <- "fixed"
K_GROWTH_FIXED <- 0 # 0 to disable soft penalty
RECRUIT_MAX_SOURCE <- "fixed"
RECRUIT_MAX_FIXED <- (MAX_GROWTH_FIXED  * 5) + 0.9999

############################################################
### 2.3 DP running settings
############################################################
DP_MODE <- "marginals+bins" # Options: "none", "marginals", "marginals+bins"
anchor_start_census <- 7L
DP_VERBOSE <- TRUE
DP_POSTERIOR_TOP_K <- 1L
dp_max_tracks <- NULL # auto (computed from data)
dp_max_states <- 40000L
dp_slack_tracks <- 1L
# NOTE: Optionally require that slack be granted only if an anchor DBH is recruitable
# (i.e., DBH <= Bio_Recruit_MaxDBH_unit + eps). Set FALSE to preserve current behavior.
dp_slack_require_anchor_recruitable <- TRUE
# Tolerance (cm) used when comparing anchor DBH to recruit_max_dbh
dp_slack_require_anchor_eps <- 1e-6

# Chunk-specific defaults
DP_CHUNK_SIZE <- 7L
DP_CHUNK_RESUME <- TRUE
DP_CHUNK_OVERWRITE <- FALSE
# Optional: limit chunks to a specific range for testing (NULL means all)
DP_CHUNK_START <- NULL
DP_CHUNK_END <- NULL

# Posterior sampling defaults (disabled by default)
# - POSTERIOR_SAMPLES: number of full-path reconstructions to draw from the DP posterior
# - POSTERIOR_SAMPLES_FORMAT: output format forwarded to DP ('rds','feather','csv')
# - POSTERIOR_SAMPLES_PATH: optional path to write posterior files; when NULL DP writes to out_dir/posteriors
# - POSTERIOR_SAMPLE_SEED: integer seed used to make posterior sampling reproducible. If NULL sampling is not deterministically seeded; runners
#   (e.g., bin/run_dp_future_single.R) may auto-generate a seed when running in batch/parallel to avoid RNG misuse warnings and ensure
#   reproducible sampling across tasks. If you want reproducible CLI runs, pass --POSTERIOR_SAMPLE_SEED explicitly.
POSTERIOR_SAMPLES <- 200L
POSTERIOR_SAMPLES_FORMAT <- "feather" # options: 'rds', 'feather', 'csv'
POSTERIOR_SAMPLES_PATH <- NULL
POSTERIOR_SAMPLE_SEED <- NULL

############################################################
### 2.4 Parallel and output settings
############################################################
RUN_ALL_TAGS <- TRUE
MANUAL_CORES <- TRUE # Flag to manually define cores instead of auto-detecting
MANUAL_CORES_VALUE <- 7L # Number of cores to use if MANUAL_CORES=TRUE

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

WRITE_DP_CSV <- FALSE
WRITE_DP_RDS <- FALSE
WRITE_DP_FEATHER <- TRUE
WRITE_DP_PDF_PER_CHUNK <- WRITE_DP_PDF <- TRUE
DP_PDF_INCLUDE_REFERENCE <- TRUE
PROJECT_ROOT <- here::here()

# Validate a few things
if (!POSTERIOR_SAMPLES_FORMAT %in% c("rds", "feather", "csv")) stop("Invalid POSTERIOR_SAMPLES_FORMAT: ", POSTERIOR_SAMPLES_FORMAT)

MC_CORES <- if (exists("MANUAL_CORES") && isTRUE(MANUAL_CORES)) as.integer(MANUAL_CORES_VALUE) else if (requireNamespace("parallel", quietly = TRUE)) max(1L, parallel::detectCores(logical = TRUE) - 1L) else 1L

############################################################
### 2.6 Realism report settings
############################################################
RUN_REALISM_REPORT <- FALSE

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

# Final output directory for this run (created at runtime in run_main())
out_dir <- file.path(base_out_dir, build_out_dir_name())
message("[dp_global main_cpp.R] out_dir (computed): ", out_dir)
message("[dp_global main_cpp.R] getwd(): ", getwd())

DP_CSV_FILE <- file.path(out_dir, "stem_reconstruction_dp_global_rcpp.csv")

# If posterior sampling is requested, provide the DP with the run's base out_dir.
# The DP itself will create the 'posteriors/' subdirectory (avoids double-nesting).
# These defaults can still be overridden via CLI (e.g., --POSTERIOR_SAMPLES_PATH=... or --POSTERIOR_SAMPLE_SEED=...)
if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L) {
    if (is.null(POSTERIOR_SAMPLES_PATH) || !nzchar(POSTERIOR_SAMPLES_PATH)) {
        # Provide the DP with the run's base out_dir; the DP will create the
        # 'posteriors/' subdirectory itself, avoiding double-nesting.
        POSTERIOR_SAMPLES_PATH <- out_dir
    }
    # normalize path for consistency; don't require it to exist yet (created in run_main)
    POSTERIOR_SAMPLES_PATH <- normalizePath(POSTERIOR_SAMPLES_PATH, winslash = "/", mustWork = FALSE)
    if (is.null(POSTERIOR_SAMPLE_SEED)) {
        POSTERIOR_SAMPLE_SEED <- as.integer(123L)
    } else {
        POSTERIOR_SAMPLE_SEED <- as.integer(POSTERIOR_SAMPLE_SEED)
    }
}

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

# (using the ensure_dir helper defined near the top of this file)

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
        posterior_samples = POSTERIOR_SAMPLES,
        posterior_samples_format = POSTERIOR_SAMPLES_FORMAT,
        posterior_samples_path = POSTERIOR_SAMPLES_PATH,
        posterior_sample_seed = POSTERIOR_SAMPLE_SEED,
        use_measurement_error = isTRUE(USE_MEASUREMENT_ERROR),
        # prune controls
        prune_hard = TRUE,
        prune_min_growth = MAX_SHRINK_FIXED * 2.5, # very wide fixed bounds
        prune_max_growth = MAX_GROWTH_FIXED * 1.5, # very wide fixed bounds
        prune_use_bio_bounds = FALSE, # use fixed prune bounds instead of biological ones
        prune_recruit_max_dbh = RECRUIT_MAX_FIXED * 1.2, # very high recruit max dbh
        prune_use_bio_recruit = FALSE, # FALSE = use prune_recruit_max_dbh instead of biological (and margin) one, TRUE, set prune_recruit_max_dbh as min(prune_recruit_max_dbh, bio_recruit_max_dbh * 1.2)
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
run_main_chunked <- function() {
    ensure_dir(out_dir)

    # Create posteriors subdirectory (DP writes its files into <base>/posteriors)
    if (!is.null(POSTERIOR_SAMPLES) && as.integer(POSTERIOR_SAMPLES) > 0L && !is.null(POSTERIOR_SAMPLES_PATH) && nzchar(POSTERIOR_SAMPLES_PATH)) {
        ensure_dir(file.path(POSTERIOR_SAMPLES_PATH, "posteriors"))
    }

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
            recruit_max_fixed = as.numeric(get0("RECRUIT_MAX_FIXED", ifnotfound = (7.5 * 5) + 0.9999))
        )
    }

    # Write a small text file recording the parameters used to build the
    # run-specific output directory name so runs are reproducible.
    # Gather all parameters in the environment for full run context

    params <- list()
    params$TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    params$OUT_DIR <- out_dir
    params$GENERATED_DIR_NAME <- basename(out_dir)

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
        tryCatch({ writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt")) }, error = function(e) NULL)
        return(invisible(list(xrun = xrun, bio_pars = bio_pars)))
    }

    # Clamp values to valid range
    start_ci <- max(1L, min(length(chunks), start_ci))
    end_ci <- max(1L, min(length(chunks), end_ci))

    if (start_ci > end_ci) {
        stop(sprintf("Invalid chunk range: DP_CHUNK_START=%s, DP_CHUNK_END=%s after clamping => start=%d > end=%d", 
                     if (exists("DP_CHUNK_START") && !is.null(DP_CHUNK_START)) as.character(DP_CHUNK_START) else "NULL",
                     if (exists("DP_CHUNK_END") && !is.null(DP_CHUNK_END)) as.character(DP_CHUNK_END) else "NULL",
                     start_ci, end_ci))
    }

    for (ci in seq_len(length(chunks))) {
        if (ci < start_ci || ci > end_ci) {
            next
        }
        chunk_rds <- file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.rds", ci))

        if (isTRUE(DP_CHUNK_RESUME) && file.exists(chunk_rds) && !isTRUE(DP_CHUNK_OVERWRITE)) {
            log_msg(sprintf("Skipping chunk %d/%d — chunk RDS exists (resume enabled)", ci, length(chunks)))
            first_chunk <- !file.exists(DP_CSV_FILE)
            next
        }

        groups_ci <- groups[chunks[[ci]]]
        log_msg(sprintf("Chunk %d/%d — %d groups", ci, length(chunks), nrow(groups_ci)))

        # Run chunk with error isolation so one bad chunk doesn't kill the run
        status <- tryCatch({
            res <- parallel::mclapply(seq_len(nrow(groups_ci)), function(j) {
                data.table::setDTthreads(1L)
                g <- groups_ci[j]
                dtg <- xrun[Tag == g$Tag & species == g$species]
                run_dp_one_group(dtg, dp_max_tracks = dp_max_tracks_local)
            }, mc.cores = MC_CORES)

            out_chunk <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)

            if (nrow(out_chunk) > 0L) {
                out_chunk[, DP_Chunk := ci]
                out_chunk <- maybe_add_posterior_bins(out_chunk)
                # Record run output directory (basename) in each row to avoid variable/column name collision
                out_chunk[, run_out_dir := basename(out_dir)]

                if (isTRUE(WRITE_DP_CSV)) {
                    data.table::fwrite(out_chunk, file = DP_CSV_FILE, append = !first_chunk)
                    log_msg(sprintf("Wrote CSV for chunk %d (nrow=%d)", ci, nrow(out_chunk)))
                }

                if (isTRUE(WRITE_DP_FEATHER)) {
                    if (!requireNamespace("arrow", quietly = TRUE)) {
                        log_msg("'arrow' package not available; skipping feather output", "WARN")
                    } else {
                        arrow::write_feather(out_chunk, file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.feather", ci)))
                        log_msg(sprintf("Wrote Feather for chunk %d", ci))
                    }
                }

                if (isTRUE(WRITE_DP_RDS)) {
                    saveRDS(out_chunk, file = chunk_rds)
                    log_msg(sprintf("Wrote RDS for chunk %d", ci))
                }

                if (isTRUE(WRITE_DP_PDF_PER_CHUNK) && isTRUE(WRITE_DP_PDF)) {
                    tryCatch({
                        plot_tag_to_pdf(out_chunk, pdf_file = file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d.pdf", ci)), include_reference = DP_PDF_INCLUDE_REFERENCE)
                        log_msg(sprintf("Wrote PDF for chunk %d", ci))
                    }, error = function(e) log_msg(sprintf("Failed to write PDF for chunk %d: %s", ci, conditionMessage(e)), "ERROR"))
                }
            } else {
                log_msg(sprintf("Chunk %d returned no rows.", ci))
                if (isTRUE(WRITE_DP_RDS)) saveRDS(out_chunk, file = chunk_rds)
            }

            TRUE
        }, error = function(e) {
            # Write a small error marker for the chunk and continue
            err_file <- file.path(out_dir, sprintf("stem_reconstruction_dp_global_rcpp_chunk_%03d_failed.txt", ci))
            writeLines(conditionMessage(e), con = err_file)
            log_msg(sprintf("Chunk %d failed: %s", ci, conditionMessage(e)), "ERROR")
            FALSE
        })

        first_chunk <- FALSE
        rm(res, out_chunk, groups_ci)
        invisible(gc())
    }

    # Write a small finished marker so users and wrappers can detect job completion
    tryCatch(
        {
            writeLines(as.character(Sys.time()), con = file.path(out_dir, "run_finished.txt"))
        },
        error = function(e) {
            message("[dp_global main_cpp_chunk.R] Warning writing run_finished marker: ", conditionMessage(e))
        }
    )

    invisible(list(xrun = xrun, bio_pars = bio_pars))
}

# Merge helpers: combine per-chunk RDS or Feather files into a single CSV
merge_chunk_rds_to_csv <- function(out_dir, out_csv = file.path(out_dir, "stem_reconstruction_dp_global_rcpp_merged.csv")) {
    files <- list.files(out_dir, pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
    if (length(files) == 0L) stop("No chunk RDS files found in ", out_dir)
    first <- TRUE
    for (f in sort(files)) {
        log_msg(paste("Merging RDS", basename(f)))
        dt <- readRDS(f)
        if (nrow(dt) == 0L) next
        if (first) {
            data.table::fwrite(dt, file = out_csv)
            first <- FALSE
        } else {
            data.table::fwrite(dt, file = out_csv, append = TRUE)
        }
        rm(dt)
        invisible(gc())
    }
    log_msg(paste("Merged", length(files), "RDS chunks to", out_csv))
    out_csv
}

merge_chunk_feathers_to_csv <- function(out_dir, out_csv = file.path(out_dir, "stem_reconstruction_dp_global_rcpp_merged.csv")) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow package required to read feather files")
    files <- list.files(out_dir, pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.feather$", full.names = TRUE)
    if (length(files) == 0L) stop("No chunk Feather files found in ", out_dir)
    first <- TRUE
    for (f in sort(files)) {
        log_msg(paste("Merging Feather", basename(f)))
        dt <- as.data.table(arrow::read_feather(f))
        if (nrow(dt) == 0L) next
        if (first) {
            data.table::fwrite(dt, file = out_csv)
            first <- FALSE
        } else {
            data.table::fwrite(dt, file = out_csv, append = TRUE)
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
        chars <- list.files(out_dir, pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.feather$", full.names = TRUE)
        if (length(chars) > 0L) return(merge_chunk_feathers_to_csv(out_dir))
        return(merge_chunk_rds_to_csv(out_dir))
    }
    # prefer rds
    rds <- list.files(out_dir, pattern = "stem_reconstruction_dp_global_rcpp_chunk_\\d{3}\\.rds$", full.names = TRUE)
    if (length(rds) > 0L) return(merge_chunk_rds_to_csv(out_dir))
    return(merge_chunk_feathers_to_csv(out_dir))
}

############################################################
### 6) Script entrypoint
############################################################

# When you run this file with Rscript, sys.nframe()==0 and we execute.
# When you source() this file from another script/session, we only define helpers.
if (sys.nframe() == 0L) {
    message("[dp_global main_cpp_chunk.R] Starting chunked run_main_chunked()")
    run_main_chunked()
}
