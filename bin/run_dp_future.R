#!/usr/bin/env Rscript
# run_dp_future.R — concurrent experiment runner (future + progressr)
#
# Overview
# -------
# Runs multiple experimental configurations concurrently on a single machine
# using the 'future' package (multisession plan). Each experimental config is
# executed by invoking `Rscript dp_global/scripts/main_cpp.R` with a set of
# command-line arguments. This script provides experiment definitions (same
# configs as the serial runner), logging, per-config overrides, and safety
# checks to avoid CPU oversubscription.
#
# Note for maintainers/orchestrators
# - This runner constructs canonical CLI flags (e.g., --POSTERIOR_SAMPLES=200)
#   expected by `dp_global/scripts/main_cpp.R`. The mapping is case-insensitive
#   and accepts '-' or '_' separators. Keep this script's MAIN_CLI_KEYS in sync
#   with the main script's `CLI_REFERENCE` to avoid mismatches.
#
# Design & resource model
# ----------------------
# - Outer concurrency (W workers): `--workers N` controls how many R
#   worker processes are launched by `future::plan(multisession, workers = N)`.
# - Inner parallelism (C cores per job): `--cores-per-job M` sets
#   `--MANUAL_CORES_VALUE=M` which is forwarded into `main_cpp.R` so the
#   per-experiment driver uses M R cores (via its own parallel/forking logic).
# - Total cores requested = W * M. The script aborts (unless `--force`) when
#   this exceeds available logical CPUs to avoid thrashing.
# - To be conservative and avoid nested threading issues with BLAS/OpenMP,
#   each launched Rscript process runs with environment variables set to
#   single-threaded by default (OMP_NUM_THREADS=1, MKL_NUM_THREADS=1,
#   OPENBLAS_NUM_THREADS=1). You can change this behavior in the script if
#   you prefer BLAS-threaded inner jobs instead of R-level forks.
#
# Logging & outputs
# -----------------
# - Per-config logs are written to `tests/parallel_future_logs/<config>.log`.
# - A job summary CSV (`parallel_future.log` by default) lists start/end times,
#   exit statuses, log file paths, and the sanitized command used.
# - For privacy, logged commands are sanitized to replace the absolute project
#   root path with `.` so your full filesystem path is not leaked in logs.
#
# Overrides & experiment tuning
# -----------------------------
# - Global overrides: `--override KEY=VAL` adds KEY=VAL to every config.
# - Per-config overrides: --cfg-override fixed:KEY=VAL adds KEY=VAL for the fixed config.
# - Extras after -- are appended to the command and are passed verbatim to
#   main_cpp.R (same place the serial runner would receive them).
# - Precedence: BASE_ARGS < config default args < extras < --override < --cfg-override
#
# Posterior sampling note:
# - If you request posterior sampling via --posterior-samples but don't provide
#   --posterior-seed, the runner will auto-generate a reproducible integer seed,
#   forward it to main_cpp.R as --POSTERIOR_SAMPLE_SEED, and record it in the joblog.
#   This reduces RNG misuse warnings from `future` and helps reproducibility.
#
# Safety & recommended workflow
# ----------------------------
# 1) DRY_RUN first to verify commands and logging: pass `-- --DRY_RUN` after
#    run_dp_future.R's options.
# 2) Inspect `tests/parallel_future_logs` and `parallel_future.log` for the
#    constructed sanitized commands before doing a real run.
# 3) Start with conservative defaults (BLAS/OMP set to 1) and monitor `htop`.
#
# Examples
# --------
# Dry-run a 4×4 experiment set (4 sessions × 4 cores each = 16 cores):
#   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_hard_soft data_soft fixed_k50 fixed_k25 data_hard_k50 data_hard_k25" -- --DRY_RUN
# Real run (same configs):
#   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_hard_soft data_soft fixed_k50 fixed_k25 data_hard_k50 data_hard_k25"
# Global override example (set DP_MODE=none for all configs):
#   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --override DP_MODE=none --configs "fixed data_hard" -- --DRY_RUN
# Per-config override example (set K_GROWTH_FIXED=50 only for 'fixed'):
#   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --cfg-override fixed:K_GROWTH_FIXED=50 --configs "fixed" -- --DRY_RUN
#
# Help/usage
# ----------
# Run with -h or --help to print a short usage summary and exit.

# Usage examples:
#   ./bin/run_dp_future.R --workers 3 --cores-per-job 5 --configs "fixed data_hard" -- --DRY_RUN

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(progressr)
})

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package to run this script.")
}
library(here)

# Prepare logging dir
log_dir <- here("tests", "parallel_future_logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

############################################################
### CLI defaults and mapping (runner-side)
############################################################
# The runner accepts a small set of flags which it maps into canonical DP
# CLI options (the DP expects flags like --POSTERIOR_SAMPLES=... etc.). Keep
# these names in sync with `dp_global/scripts/main_cpp.R`.
MAIN_CLI_KEYS <- c(
  "POSTERIOR_SAMPLES", "POSTERIOR_SAMPLES_FORMAT", "POSTERIOR_SAMPLE_SEED", "POSTERIOR_SAMPLES_PATH", "PROJECT_ROOT", "BATCH_TS", "CONFIG_NAME",
  # support the growth-form fallback parameter
  "DP_FALLBACK_GROWTH_FORMS"
)

# Runner options (lowercase names for user-facing flags)
opt <- list(
  workers = 3L,
  cores_per_job = 5L,
  configs = NULL,
  joblog = "parallel_future.log",
  force = FALSE,
  overrides = character(0),
  cfg_overrides = list(),
  dry_run = FALSE,
  no_blas_limit = FALSE,
  posterior_samples = 0L, # maps -> POSTERIOR_SAMPLES
  posterior_format = NULL, # maps -> POSTERIOR_SAMPLES_FORMAT
  posterior_seed = NULL, # maps -> POSTERIOR_SAMPLE_SEED
  posterior_path = NULL # maps -> POSTERIOR_SAMPLES_PATH
)

extras <- c()

i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("-j", "--workers")) {
    opt$workers <- as.integer(args[i + 1])
    i <- i + 2
  } else if (a == "--cores-per-job") {
    opt$cores_per_job <- as.integer(args[i + 1])
    i <- i + 2
  } else if (a == "--configs") {
    opt$configs <- strsplit(args[i + 1], "[ ,]+")[[1]]
    i <- i + 2
  } else if (a == "--joblog") {
    opt$joblog <- args[i + 1]
    i <- i + 2
  } else if (a == "--force") {
    opt$force <- TRUE
    i <- i + 1
  } else if (grepl("^--override=", a)) {
    opt$overrides <- c(opt$overrides, sub("^--override=", "", a))
    i <- i + 1
  } else if (a == "--override") {
    opt$overrides <- c(opt$overrides, args[i + 1])
    i <- i + 2
  } else if (grepl("^--cfg-override=", a)) {
    co <- sub("^--cfg-override=", "", a)
    parts <- strsplit(co, ":", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      cfg <- parts[1]
      kv <- parts[2]
      opt$cfg_overrides[[cfg]] <- c(opt$cfg_overrides[[cfg]], kv)
    }
    i <- i + 1
  } else if (a == "--cfg-override") {
    co <- args[i + 1]
    parts <- strsplit(co, ":", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      cfg <- parts[1]
      kv <- parts[2]
      opt$cfg_overrides[[cfg]] <- c(opt$cfg_overrides[[cfg]], kv)
    }
    i <- i + 2
  } else if (a == "--help" || a == "-h") {
    cat("Usage: run_dp_future_single.R [--workers N] [--cores-per-job N] [--joblog file] [--force] [--dry-run] [--no-blas-limit] [--override KEY=VAL] [--cfg-override fixed:KEY=VAL] [--posterior-samples N] [--posterior-format rds|feather|csv] [--posterior-seed N] [--posterior-path /path/to/dir] -- [extra args passed to main_cpp.R]\n")
    cat("\nNote: If --posterior-samples is supplied but no --posterior-seed is provided, the runner will auto-generate a reproducible seed and forward it to the DP (recorded in the joblog).\n")
    cat("\nCanonical DP CLI keys this runner will produce (case-insensitive, '-' or '_' allowed):\n")
    cat(sprintf("  %s\n", paste(MAIN_CLI_KEYS, collapse = ", ")))
    q(status = 0)
  } else if (a == "--dry-run" || a == "--dry_run") {
    opt$dry_run <- TRUE
    i <- i + 1
  } else if (a == "--no-blas-limit") {
    # Do not set OMP/MKL/OPENBLAS env vars; allow BLAS threading
    opt$no_blas_limit <- TRUE
    i <- i + 1
  } else if (a == "--posterior-samples") {
    opt$posterior_samples <- as.integer(args[i + 1])
    i <- i + 2
  } else if (a == "--posterior-format") {
    opt$posterior_format <- args[i + 1]
    i <- i + 2
  } else if (a == "--posterior-seed") {
    opt$posterior_seed <- as.integer(args[i + 1])
    i <- i + 2
  } else if (a == "--posterior-seed-auto") {
    # Internal/testing flag: can be used to force auto-seed behavior explicitly
    opt$posterior_seed <- NULL
    i <- i + 1
  } else if (a == "--posterior-path") {
    opt$posterior_path <- args[i + 1]
    i <- i + 2
  } else if (a == "--") {
    extras <- args[(i + 1):length(args)]
    break
  } else {
    extras <- c(extras, a)
    i <- i + 1
  }
}

if (is.null(opt$configs) || length(opt$configs) == 0L) {
  opt$configs <- c("fixed", "data_hard", "data_hard_soft", "data_soft", "fixed_k50", "fixed_k25", "data_hard_k50", "data_hard_k25")
} 

# Auto-generate a reproducible posterior seed when posterior sampling is requested
# and no explicit seed was provided. This helps avoid the future RNG misuse warning
# and ensures reproducible sampling across parallel runs. The generated seed is
# forwarded as `--POSTERIOR_SAMPLE_SEED=` to `main_cpp.R` and recorded in the joblog.
if (!is.null(opt$posterior_samples) && as.integer(opt$posterior_samples) > 0L && is.null(opt$posterior_seed)) {
  opt$posterior_seed <- as.integer(sample.int(.Machine$integer.max - 1L, 1))
  cat(sprintf("[run_dp_future] Auto-generated posterior seed: %d\n", opt$posterior_seed))
  flush.console()
}

# Safety check: avoid oversubscription
avail_cores <- tryCatch(
  parallel::detectCores(logical = TRUE),
  error = function(e) 1L
)
required <- opt$workers * opt$cores_per_job
if (required > avail_cores && !isTRUE(opt$force)) {
  stop(sprintf("ERROR: Requested cores (%d = %d workers × %d cores/worker) exceed available logical CPUs (%d). Use --force to override.", required, opt$workers, opt$cores_per_job, avail_cores))
}

# Place joblog CSV inside the same log directory as per-config logs for consistent UX
joblog_path <- file.path(log_dir, opt$joblog)

# Experiment BASE_ARGS (mirrors bin/run_dp_full_cpp.sh defaults)
BATCH_TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
BASE_ARGS <- c(
  paste0("--INPUT_FILE=", here("data_simulation", "data", "simulated_data_1.csv")),
  "--FORCE_ONE_SPECIES_PARAMETERS=FALSE",
  "--DP_MODE=marginals+bins",
  "--WHICH_TAG=20",
  "--ANCHOR_START_CENSUS=7",
  "--DP_VERBOSE=TRUE",
  "--RUN_ALL_TAGS=TRUE",
  # Allow overriding DP enumerator state budget from the orchestrator
  # (main_cpp.R defines `dp_max_states` in the DP settings section)
  # Use 0 to run igraph for all experiments in this script
  "--DP_MAX_STATES=40000", # default in main_cpp.R
  "--MANUAL_CORES=TRUE",
  sprintf("--MANUAL_CORES_VALUE=%d", opt$cores_per_job),
  "--DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE=TRUE",
  "--DP_SLACK_TRACKS=1",
  "--WRITE_DP_CSV=TRUE",
  "--WRITE_DP_RDS=TRUE",
  "--WRITE_DP_PDF=TRUE", # true to generate PDFs - set false if many Tags to save time/disk
  "--DP_PDF_INCLUDE_REFERENCE=TRUE",
  "--PLOT_PDF_ONE_TAG_ONLY=FALSE",
  "--SENSITIVITY_MODE=none", # run+write+pdf
  "--RUN_K_SWEEP_DEMO=FALSE",
  "--RUN_REALISM_REPORT=FALSE",
  paste0("--PROJECT_ROOT=", here::here()),
  paste0("--BATCH_TS=", BATCH_TS)
) 

get_config_args <- function(cfg) {
  switch(cfg,
    fixed = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=fixed",
      "--MAX_GROWTH_FIXED=7.5",
      "--MAX_SHRINK_HARD_SOURCE=fixed",
      "--MAX_SHRINK_FIXED=-0.5",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=0",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=0",
      "--RECRUIT_MAX_SOURCE=fixed",
      "--RECRUIT_MAX_FIXED=38"
    ),
    data_hard = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=data",
      "--MAX_SHRINK_HARD_SOURCE=data",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=0",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=0",
      "--RECRUIT_MAX_SOURCE=data"
    ),
    data_hard_soft = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=data",
      "--MAX_SHRINK_HARD_SOURCE=data",
      "--K_SHRINK_SOURCE=data",
      "--K_GROWTH_SOURCE=data"
    ),
    data_soft = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=fixed",
      "--MAX_GROWTH_FIXED=7.5",
      "--MAX_SHRINK_HARD_SOURCE=fixed",
      "--MAX_SHRINK_FIXED=-0.5",
      "--K_SHRINK_SOURCE=data",
      "--K_GROWTH_SOURCE=data"
    ),
    fixed_k50 = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=fixed",
      "--MAX_GROWTH_FIXED=7.5",
      "--MAX_SHRINK_HARD_SOURCE=fixed",
      "--MAX_SHRINK_FIXED=-0.5",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=50",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=50"
    ),
    fixed_k25 = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=fixed",
      "--MAX_GROWTH_FIXED=7.5",
      "--MAX_SHRINK_HARD_SOURCE=fixed",
      "--MAX_SHRINK_FIXED=-0.5",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=25",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=25"
    ),
    data_hard_k50 = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=data",
      "--MAX_SHRINK_HARD_SOURCE=data",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=50",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=50"
    ),
    data_hard_k25 = c(
      "--USE_MEASUREMENT_ERROR=TRUE",
      "--MAX_GROWTH_HARD_SOURCE=data",
      "--MAX_SHRINK_HARD_SOURCE=data",
      "--K_SHRINK_SOURCE=fixed",
      "--K_SHRINK_FIXED=25",
      "--K_GROWTH_SOURCE=fixed",
      "--K_GROWTH_FIXED=25"
    ),
    stop("Unknown config: ", cfg)
  )
}

# Preflight check: ensure the dp_main driver exists and is reachable from the project root
if (!file.exists(here("dp_global", "scripts", "main_cpp.R"))) {
  stop("dp_global/scripts/main_cpp.R not found. Run this script from the project root or set PROJECT_ROOT accordingly.")
}

# Setup future plan
plan(multisession, workers = opt$workers)
# Use explicit handler registration that works both interactively and in
# non-interactive Rscript runs (e.g., terminal, cron, SLURM). In interactive
# sessions prefer global handlers; in non-interactive runs prefer a combination
# of txtprogressbar and progress handlers which print to stdout reliably.
tryCatch(
  {
    if (interactive()) {
      handlers(global = TRUE)
    } else {
      # Explicitly register handlers that print to stdout and avoid clearing so
      # output is visible in non-interactive terminal runs (Rscript).
      handlers(list(
        progressr::handler_txtprogressbar(clear = FALSE, file = stdout()),
        progressr::handler_progress()
      ))
      # Ensure progressr is enabled for non-interactive sessions
      options(progressr.enable = TRUE)
    }
  },
  error = function(e) {
    # As a last resort, try to register a global handler silently so we don't
    # abort the whole run due to progress handler issues.
    try(handlers(global = TRUE), silent = TRUE)
  }
)

# Canonicalize overrides into --KEY=VAL uppercase form so main_cpp.R receives consistent keys
canonicalize_override <- function(ov) {
  if (!nzchar(ov) || is.na(ov)) return(ov)
  if (grepl("=", ov, fixed = TRUE)) {
    parts <- strsplit(ov, "=", fixed = TRUE)[[1]]
    key <- toupper(gsub("[- ]", "_", parts[1]))
    val <- parts[2]
    return(paste0("--", key, "=", val))
  }
  key <- toupper(gsub("[- ]", "_", ov))
  paste0("--", key)
}
norm_overrides <- if (length(opt$overrides) > 0) vapply(opt$overrides, canonicalize_override, character(1)) else character(0)

# Build commands: BASE_ARGS + config-specific args + extras + global overrides + per-config overrides
commands <- lapply(opt$configs, function(cfg) {
  cfg_args <- get_config_args(cfg)
  cfg_specific_raw <- if (!is.null(opt$cfg_overrides[[cfg]])) opt$cfg_overrides[[cfg]] else character(0)
  cfg_specific <- if (length(cfg_specific_raw) > 0) vapply(cfg_specific_raw, canonicalize_override, character(1)) else character(0)
  # precedence: BASE_ARGS < cfg_args < extras < global overrides < cfg_specific
  cmd <- c(BASE_ARGS, cfg_args, extras, norm_overrides, cfg_specific, paste0("--CONFIG_NAME=", cfg))

  # Helper to form canonical flags expected by `main_cpp.R` (uppercase, underscores)
  make_main_flag <- function(k, v) paste0("--", toupper(gsub("[- ]", "_", k)), "=", v)

  # Append posterior sampling flags when requested. If posterior_path is NULL
  # the DP will write samples to its own out_dir/posteriors by default.
  if (!is.null(opt$posterior_samples) && as.integer(opt$posterior_samples) > 0L) {
    cmd <- c(cmd, make_main_flag("posterior-samples", as.integer(opt$posterior_samples)))
    if (!is.null(opt$posterior_format) && nzchar(opt$posterior_format)) {
      cmd <- c(cmd, make_main_flag("posterior-samples-format", opt$posterior_format))
    }
    if (!is.null(opt$posterior_seed)) {
      cmd <- c(cmd, make_main_flag("posterior-sample-seed", as.integer(opt$posterior_seed)))
    }
    if (!is.null(opt$posterior_path) && nzchar(opt$posterior_path)) {
      cmd <- c(cmd, make_main_flag("posterior-samples-path", opt$posterior_path))
    }
  }

  list(config = cfg, cmd = cmd)
})

# Validate overrides to catch obvious typos early (map KEY=VAL pairs and warn
# if the KEY doesn't resemble a known main CLI key; this helps keep runner and
# main in sync and surfaces user mistakes instead of silent mis-overrides).
validate_override_key <- function(k) {
  normalize <- function(x) toupper(gsub("[- ]", "_", x))
  nk <- normalize(k)
  known_keys <- toupper(c(
    MAIN_CLI_KEYS,
    "INPUT_FILE","FORCE_ONE_SPECIES_PARAMETERS","DP_MODE","WHICH_TAG","ANCHOR_START_CENSUS",
    "DP_VERBOSE","RUN_ALL_TAGS","MANUAL_CORES","MANUAL_CORES_VALUE",
    "WRITE_DP_CSV","WRITE_DP_RDS","WRITE_DP_PDF","WRITE_DP_FEATHER",
    "DP_MAX_STATES","DP_SLACK_TRACKS","DP_SLACK_REQUIRE_ANCHOR_RECRUITABLE","DP_SLACK_REQUIRE_ANCHOR_EPS",
    "USE_MEASUREMENT_ERROR"
  ))
  if (!(nk %in% known_keys)) {
    warning(sprintf("[run_dp_future] Override key '%s' does not match known main options; check spelling.", k))
  }
}

# Check global overrides
if (length(opt$overrides) > 0) {
  for (ov in opt$overrides) {
    if (grepl("=", ov)) {
      k <- strsplit(ov, "=", fixed = TRUE)[[1]][1]
      validate_override_key(k)
    }
  }
}

# Check cfg-specific overrides
if (length(opt$cfg_overrides) > 0) {
  for (cfg in names(opt$cfg_overrides)) {
    for (ov in opt$cfg_overrides[[cfg]]) {
      if (grepl("=", ov)) {
        k <- strsplit(ov, "=", fixed = TRUE)[[1]][1]
        validate_override_key(k)
      }
    }
  }
}
# Total tasks to launch (one future per config)
total_tasks <- length(commands)
cat(sprintf("Launching %d tasks...\n", total_tasks))
flush.console()

# Create futures (non-blocking)
futures_list <- lapply(seq_along(commands), function(i) {
  entry <- commands[[i]]
  future(
    {
    cfg <- entry$config
    cmd <- entry$cmd
      # Unique log filename to avoid overwriting previous runs
      log_suffix <- paste0(BATCH_TS, "_", format(Sys.time(), "%H%M%S"), "_", sprintf("%06d", sample.int(1e6, 1)))
      log_file <- file.path(log_dir, paste0(cfg, "_", log_suffix, ".log"))
    # Ensure the log directory exists even if removed or not visible to workers
    log_dir_local <- dirname(log_file)
    if (!dir.exists(log_dir_local)) dir.create(log_dir_local, recursive = TRUE, showWarnings = FALSE)

    start <- Sys.time()
    cat(sprintf("[%s] STARTED at %s\n", cfg, format(start, tz = "UTC")))
    flush.console()
      # Record basic info in the final log immediately so we have a trace even if
      # later steps fail inside the worker
      cat(sprintf("FUTURE: pid=%d, start=%s\n", Sys.getpid(), format(start, tz = "UTC")), file = log_file, append = TRUE)
      cat(sprintf("FUTURE: R.version=%s, R.home=%s\n", paste(R.version, collapse = " "), R.home()), file = log_file, append = TRUE)

    # Build the full command string for logging (project-relative)
    cmd_str_raw <- paste(c("Rscript", shQuote("dp_global/scripts/main_cpp.R"), vapply(cmd, shQuote, character(1))), collapse = " ")
    cmd_str <- gsub(here::here(), ".", cmd_str_raw, fixed = TRUE)

      # Top-level dry-run override or explicit --DRY_RUN in args means we don't execute
      is_dry <- isTRUE(opt$dry_run) || any(grepl("--DRY_RUN", cmd))

      if (is_dry) {
        out_header <- c(sprintf("DRY RUN: %s", cmd_str), sprintf("Generated log file: %s", log_file))
        writeLines(out_header, con = log_file)
      exit_status <- 0L
      cat(sprintf("[DRY_RUN] %s: %s\n", cfg, cmd_str))
      flush.console()
      # Early return for dry-run to avoid runtime errors in future worker and
      # to ensure the future resolves to a clean success (status = 0)
      return(list(config = cfg, start = start, end = start, status = 0L, log = log_file, cmd = cmd_str, main_out_dir = NA_character_))
    } else {
        # Prepare env vars to limit BLAS/OMP threading unless disabled by flag
        if (isTRUE(opt$no_blas_limit)) {
          env_vars <- character(0)
        } else {
        env_vars <- c("OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1")
        }

        # Run command and stream stdout/stderr to a temp file then prepend header into final log
        temp_log <- tempfile(pattern = paste0(cfg, "_"), tmpdir = tempdir(), fileext = ".log")
        header_lines <- c(sprintf("COMMAND: %s", cmd_str), sprintf("ENV: %s", paste(env_vars, collapse = " ")), sprintf("START: %s", format(start, tz = "UTC")))
        writeLines(header_lines, con = log_file)

        # Add diagnostics to temp_log to help debug worker environment issues
        diag_lines <- c(
          sprintf("worker getwd: %s", getwd()),
          sprintf("Rscript (Sys.which): %s", Sys.which("Rscript")),
          sprintf("R.home(bin): %s", file.path(R.home("bin"), "Rscript")),
          sprintf("main_cpp exists: %s", file.exists("dp_global/scripts/main_cpp.R"))
        )
        # Use cat() with append to reliably add diagnostics to the temp log
        cat(paste0(diag_lines, collapse = "\n"), file = temp_log, sep = "\n", append = TRUE)

        # Prefer an explicit Rscript binary (fall back to R.home if PATH not available in worker)
        rscript_bin <- Sys.which("Rscript")
        if (!nzchar(rscript_bin)) rscript_bin <- file.path(R.home("bin"), "Rscript")

        # Log the exact command we will call into the final log for easier debugging
        cat(sprintf("FUTURE: About to call system2: rscript_bin=%s, temp_log=%s\n", rscript_bin, temp_log), file = log_file, append = TRUE)

        status <- tryCatch(
          {
            system2(rscript_bin, args = c("dp_global/scripts/main_cpp.R", cmd), stdout = temp_log, stderr = temp_log, env = env_vars)
          },
          error = function(e) {
            # system2 may throw in rare cases; capture message in both temp_log and final log
            cat(paste0("ERROR: ", conditionMessage(e), "\n"), file = temp_log, append = TRUE)
            cat(paste0("ERROR: system2 threw in worker: ", conditionMessage(e), "\n"), file = log_file, append = TRUE)
            1L
          }
        )
        # Log the return status in the final log as well
        cat(sprintf("FUTURE: system2 returned status=%s\n", as.character(status)), file = log_file, append = TRUE)

        # Append temp_log contents to final log (header already written)
        if (file.exists(temp_log)) {
          tryCatch(
            {
              logs <- readLines(temp_log, warn = FALSE)
              if (length(logs) > 0) cat(paste0(logs, collapse = "\n"), file = log_file, sep = "\n", append = TRUE, useBytes = TRUE)

              # Parse main script timestamps if the log contains standardized messages
              main_start <- NA
              main_end <- NA
              started_lines <- grep("Started run", logs, value = TRUE)
              if (length(started_lines) > 0) {
                ts_str <- sub("^\\[(.*?)\\].*$", "\\1", started_lines[[1]])
                main_start <- tryCatch(as.POSIXct(ts_str, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), error = function(e) NA)
                if (!is.na(main_start)) cat(sprintf("FUTURE: detected main 'Started run' at %s\n", format(main_start, tz = "UTC")), file = log_file, append = TRUE)
              }
              finished_lines <- grep("Finished run", logs, value = TRUE)
              if (length(finished_lines) > 0) {
                ts_str2 <- sub("^\\[(.*?)\\].*$", "\\1", finished_lines[[1]])
                main_end <- tryCatch(as.POSIXct(ts_str2, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), error = function(e) NA)
                if (!is.na(main_end)) cat(sprintf("FUTURE: detected main 'Finished run' at %s\n", format(main_end, tz = "UTC")), file = log_file, append = TRUE)
              }

              # If main wrote an output dir with our BATCH_TS and config name, try to include its run_log.txt
              out_base <- file.path(here::here(), "dp_global", "output")
              if (dir.exists(out_base)) {
                pattern <- paste0(BATCH_TS, ".*", entry$config)
                dirs <- list.dirs(out_base, full.names = TRUE, recursive = FALSE)
                matched <- dirs[grepl(pattern, basename(dirs))]
                if (length(matched) > 0) {
                  m_info <- file.info(matched)$mtime
                  chosen <- matched[which.max(m_info)]
                  main_out_dir <- chosen
                  run_log_path <- file.path(chosen, "run_log.txt")
                  if (file.exists(run_log_path)) {
                    cat(sprintf("FUTURE: including %s into %s\n", run_log_path, log_file), file = log_file, append = TRUE)
                    run_log_lines <- readLines(run_log_path, warn = FALSE)
                    if (length(run_log_lines) > 0) cat(paste0("\n--- run_log.txt (from main run) ---\n", paste(run_log_lines, collapse = "\n"), "\n--- end run_log.txt ---\n"), file = log_file, append = TRUE)
                  }
                } else {
                  main_out_dir <- NA_character_
                }
              }

              unlink(temp_log)
              cat("FUTURE: appended temp_log to log_file and unlinked temp_log\n", file = log_file, append = TRUE)
            },
            error = function(e) {
              cat(sprintf("FUTURE: failed to append or unlink temp_log: %s\n", conditionMessage(e)), file = log_file, append = TRUE)
            }
          )
        }

        exit_status <- as.integer(status)
        # Prefer timestamps reported by the main run (if detected in logs) for accuracy
        if (exists("main_start") && !is.na(main_start)) {
          start <- main_start
    }
        if (exists("main_end") && !is.na(main_end)) {
          end_override <- main_end
        } else {
          end_override <- NULL
        }
      }

      end <- if (!is.null(end_override)) end_override else Sys.time()
    cat(sprintf("[%s] DONE at %s (status=%d) log=%s\n", cfg, format(end, tz = "UTC"), exit_status, log_file))
    flush.console()

    list(config = cfg, start = start, end = end, status = exit_status, log = log_file, cmd = cmd_str, main_out_dir = if (exists('main_out_dir')) main_out_dir else NA_character_)
    },
    seed = TRUE
  )
})

# Monitor futures and print ETA information as tasks complete
results <- vector("list", length(futures_list))
completed <- rep(FALSE, length(futures_list))
check_interval <- 5 # seconds between status prints
last_print <- Sys.time()
while (!all(completed)) {
  for (i in seq_along(futures_list)) {
    if (!completed[i] && resolved(futures_list[[i]])) {
      # Try to get value; on error, record failure info
      val <- tryCatch(value(futures_list[[i]]), error = function(e) {
        entry <- commands[[i]]
        # Fallback: find the most recent log file for this config if available
        logs_found <- list.files(log_dir, pattern = paste0("^", entry$config, "_.*\\.log$"), full.names = TRUE)
        log_file_fallback <- if (length(logs_found) > 0) tail(sort(logs_found), 1) else file.path(log_dir, paste0(entry$config, ".log"))
        # Record the future error message in the fallback log for diagnosability
        tryCatch(
          writeLines(c("ERROR: Future error in worker:", conditionMessage(e)), con = log_file_fallback, sep = "\n", useBytes = TRUE, append = TRUE),
          error = function(we) NULL
        )
        list(config = entry$config, start = NA, end = NA, status = 1L, log = log_file_fallback, cmd = "", error = conditionMessage(e))
      })
      results[[i]] <- val
      completed[i] <- TRUE

      # Compute ETA using mean duration of completed tasks
      finished_vals <- Filter(Negate(is.null), results[completed])
      durations <- sapply(finished_vals, function(r) {
        if (is.na(r$start) || is.na(r$end)) {
          return(NA_real_)
        }
        as.numeric(difftime(r$end, r$start, units = "secs"))
      })
      durations <- durations[!is.na(durations)]
      finished_count <- sum(completed)
      remaining <- total_tasks - finished_count
      if (length(durations) > 0 && remaining > 0) {
        avg <- mean(durations)
        est_remaining <- avg * remaining
        eta_time <- Sys.time() + est_remaining
        cat(sprintf(
          "[ETA] %d/%d done, avg=%.1fs, remaining=%d, est_remaining=%.1fs (finish at %s)\n",
          finished_count, total_tasks, avg, remaining, est_remaining, format(eta_time, tz = "UTC")
        ))
      } else {
        cat(sprintf("[PROGRESS] %d/%d done\n", finished_count, total_tasks))
      }
      flush.console()
    }
  }
  if (!all(completed) && difftime(Sys.time(), last_print, units = "secs") >= check_interval) {
    pending <- sum(!completed)
    cat(sprintf("[STATUS] %d tasks pending...\n", pending))
    flush.console()
    last_print <- Sys.time()
  }
  if (!all(completed)) Sys.sleep(1)
}

# Write joblog CSV
joblog_rows <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    BATCH_TS = BATCH_TS,
    workers = opt$workers,
    cores_per_job = opt$cores_per_job,
    dry_run = as.logical(opt$dry_run),
    config = r$config,
    start = format(r$start, tz = "UTC"),
    end = format(r$end, tz = "UTC"),
    status = r$status,
    log = r$log,
    cmd = r$cmd,
    main_out_dir = if (is.null(r$main_out_dir)) NA_character_ else as.character(r$main_out_dir),
    posterior_samples = opt$posterior_samples,
    posterior_format = if (is.null(opt$posterior_format)) NA_character_ else opt$posterior_format,
    posterior_seed = if (is.null(opt$posterior_seed)) NA_integer_ else as.integer(opt$posterior_seed),
    posterior_path = if (is.null(opt$posterior_path)) NA_character_ else opt$posterior_path,
    stringsAsFactors = FALSE
  )
}))
write.csv(joblog_rows, file = joblog_path, row.names = FALSE)

# Summary
total <- nrow(joblog_rows)
failed <- sum(joblog_rows$status != 0)
succeeded <- total - failed
cat(sprintf("\nParallel future run summary: %d tasks, %d succeeded, %d failed. Joblog: %s\n", total, succeeded, failed, joblog_path))

if (failed > 0) {
  cat("Failed tasks:\n")
  print(joblog_rows[joblog_rows$status != 0, c("config", "status", "log")])
  q(status = 1)
}

q(status = 0)

# Rscript bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard data_hard_soft data_soft fixed_k50 fixed_k25 data_hard_k50 data_hard_k25" -- --DRY_RUN

# Rscript bin/run_dp_future.R --workers 3 --cores-per-job 5 --configs "fixed data_hard data_hard_soft data_soft fixed_k50 fixed_k25 data_hard_k50 data_hard_k25"

# Rscript bin/run_dp_future.R --workers 1 --cores-per-job 14 --configs "fixed"

# ./bin/run_dp_future.R --workers 1 --cores-per-job 14 --configs "fixed" --DP_FALLBACK_GROWTH_FORMS=tree
# ./bin/run_dp_future.R --workers 1 --cores-per-job 14 --configs "fixed" --DP_FALLBACK_GROWTH_FORMS=tree,fig

## NOTE: If you want to run all experiments for one tag only, set --WHICH_TAG=N --RUN_ALL_TAGS=FALSE
# Rscript bin/run_dp_future.R --workers 8 --cores-per-job 1 --WHICH_TAG=20 --RUN_ALL_TAGS=FALSE  --configs "fixed data_hard data_hard_soft data_soft fixed_k50 fixed_k25 data_hard_k50 data_hard_k25"