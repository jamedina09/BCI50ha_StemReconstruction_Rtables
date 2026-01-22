#!/usr/bin/env Rscript
# run_dp_future_single.R — concurrent experiment runner for fixed configuration (future + progressr)
#
# Overview
# -------
# Runs the 'fixed' experimental configuration concurrently on a single machine
# using the 'future' package (multisession plan). The experimental config is
# executed by invoking `Rscript dp_global/scripts/main_cpp.R` with a set of
# command-line arguments. This script provides the fixed config definition,
# logging, per-config overrides, and safety checks to avoid CPU oversubscription.
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
# - Global overrides: --override KEY=VAL adds KEY=VAL to the config.
# - Per-config overrides: --cfg-override fixed:KEY=VAL adds KEY=VAL for the fixed config.
# - Extras after -- are appended to the command and are passed verbatim to
#   main_cpp.R (same place the serial runner would receive them).
# - Precedence: BASE_ARGS < config default args < extras < --override < --cfg-override
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
# Dry-run the fixed experiment set:
#   ./bin/run_dp_future_single.R --workers 4 --cores-per-job 4 -- --DRY_RUN
# Real run:
#   ./bin/run_dp_future_single.R --workers 4 --cores-per-job 4
# Global override example (set DP_MODE=none):
#   ./bin/run_dp_future_single.R --workers 4 --cores-per-job 4 --override DP_MODE=none -- --DRY_RUN
# Per-config override example (set K_GROWTH_FIXED=50):
#   ./bin/run_dp_future_single.R --workers 4 --cores-per-job 4 --cfg-override fixed:K_GROWTH_FIXED=50 -- --DRY_RUN
#
# Help/usage
# ----------
# Run with -h or --help to print a short usage summary and exit.

# Usage examples:
#   ./bin/run_dp_future_single.R --workers 3 --cores-per-job 5 -- --DRY_RUN

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
# Simple arg parsing
opt <- list(workers = 3L, cores_per_job = 5L, configs = NULL, joblog = "parallel_future.log", force = FALSE, overrides = character(0), cfg_overrides = list())
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
    cat("Usage: run_dp_future_single.R [--workers N] [--cores-per-job N] [--joblog file] [--force] [--override KEY=VAL] [--cfg-override fixed:KEY=VAL] -- [extra args passed to main_cpp.R]\n")
    q(status = 0)
  } else if (a == "--") {
    extras <- args[(i + 1):length(args)]
    break
  } else {
    extras <- c(extras, a)
    i <- i + 1
  }
}

if (is.null(opt$configs) || length(opt$configs) == 0L) {
  opt$configs <- c("fixed")
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
  paste0("--input_file=", here("data_simulation", "data", "simulated_data_1.csv")),
  "--FORCE_ONE_SPECIES_PARAMETERS=FALSE",
  "--DP_MODE=marginals+bins",
  "--which_tag=1",
  "--anchor_start_census=7",
  "--DP_VERBOSE=TRUE",
  "--RUN_ALL_TAGS=FALSE",
  # Allow overriding DP enumerator state budget from the orchestrator
  # (main_cpp.R defines `dp_max_states` in the DP settings section)
  # Use 0 to run igraph for all experiments in this script
  "--dp_max_states=40000", # default in main_cpp.R
  "--MANUAL_CORES=TRUE",
  sprintf("--MANUAL_CORES_VALUE=%d", opt$cores_per_job),
  "--WRITE_DP_CSV=TRUE",
  "--WRITE_DP_RDS=TRUE",
  "--WRITE_DP_PDF=TRUE", # true to generate PDFs - set false if many Tags to save time/disk
  "--DP_PDF_INCLUDE_REFERENCE=TRUE",
  "--PLOT_PDF_ONE_TAG_ONLY=FALSE",
  "--SENSITIVITY_MODE=none", #run+write+pdf
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
      "--RECRUIT_MAX_FIXED=6"
    ),
    stop("Unknown config: ", cfg)
  )
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

# Build commands: BASE_ARGS + config-specific args + extras + global overrides + per-config overrides
commands <- lapply(opt$configs, function(cfg) {
  cfg_args <- get_config_args(cfg)
  cfg_specific <- if (!is.null(opt$cfg_overrides[[cfg]])) opt$cfg_overrides[[cfg]] else character(0)
  # precedence: BASE_ARGS < cfg_args < extras < global overrides < cfg_specific
  cmd <- c(BASE_ARGS, cfg_args, extras, opt$overrides, cfg_specific, paste0("--CONFIG_NAME=", cfg))
  list(config = cfg, cmd = cmd)
})

# Launch each config as an individual future and monitor progress with simple
# ETA/summary messages (progress bar removed to avoid terminal incompatibilities).
# This approach provides explicit STARTED/DONE lines and an ETA after each
# completed task based on the mean duration of finished tasks.

total_tasks <- length(commands)
cat(sprintf("Launching %d tasks...\n", total_tasks))
flush.console()

# Create futures (non-blocking)
futures_list <- lapply(seq_along(commands), function(i) {
  entry <- commands[[i]]
  future({
    cfg <- entry$config
    cmd <- entry$cmd
    log_file <- file.path(log_dir, paste0(cfg, ".log"))
    # Ensure the log directory exists even if removed or not visible to workers
    log_dir_local <- dirname(log_file)
    if (!dir.exists(log_dir_local)) dir.create(log_dir_local, recursive = TRUE, showWarnings = FALSE)

    start <- Sys.time()
    cat(sprintf("[%s] STARTED at %s\n", cfg, format(start, tz = "UTC")))
    flush.console()

    # Build the full command string for logging (project-relative)
    cmd_str_raw <- paste(c("Rscript", shQuote("dp_global/scripts/main_cpp.R"), vapply(cmd, shQuote, character(1))), collapse = " ")
    cmd_str <- gsub(here::here(), ".", cmd_str_raw, fixed = TRUE)

    if (any(grepl("--DRY_RUN", cmd))) {
      out <- paste("DRY RUN:", cmd_str)
      exit_status <- 0L
      writeLines(out, con = log_file)
      cat(sprintf("[DRY_RUN] %s: %s\n", cfg, cmd_str))
      flush.console()
    } else {
      res <- tryCatch({
        env_vars <- c("OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1")
        out <- system2("Rscript", args = c("dp_global/scripts/main_cpp.R", cmd), stdout = TRUE, stderr = TRUE, env = env_vars)
        writeLines(out, con = log_file)
        list(status = 0L, out = out)
      }, error = function(e) {
        writeLines(c(paste0("ERROR: ", conditionMessage(e)), ""), con = log_file)
        list(status = 1L, out = conditionMessage(e))
      })
      exit_status <- as.integer(res$status)
    }

    end <- Sys.time()
    cat(sprintf("[%s] DONE at %s (status=%d) log=%s\n", cfg, format(end, tz = "UTC"), exit_status, log_file))
    flush.console()

    list(config = cfg, start = start, end = end, status = exit_status, log = log_file, cmd = cmd_str)
  }, future.seed = TRUE)
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
        list(config = entry$config, start = NA, end = NA, status = 1L, log = file.path(log_dir, paste0(entry$config, ".log")), cmd = "")
      })
      results[[i]] <- val
      completed[i] <- TRUE

      # Compute ETA using mean duration of completed tasks
      finished_vals <- Filter(Negate(is.null), results[completed])
      durations <- sapply(finished_vals, function(r) {
        if (is.na(r$start) || is.na(r$end)) return(NA_real_)
        as.numeric(difftime(r$end, r$start, units = "secs"))
      })
      durations <- durations[!is.na(durations)]
      finished_count <- sum(completed)
      remaining <- total_tasks - finished_count
      if (length(durations) > 0 && remaining > 0) {
        avg <- mean(durations)
        est_remaining <- avg * remaining
        eta_time <- Sys.time() + est_remaining
        cat(sprintf("[ETA] %d/%d done, avg=%.1fs, remaining=%d, est_remaining=%.1fs (finish at %s)\n",
                    finished_count, total_tasks, avg, remaining, est_remaining, format(eta_time, tz = "UTC")))
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
  data.frame(config = r$config, start = format(r$start, tz = "UTC"), end = format(r$end, tz = "UTC"), status = r$status, log = r$log, cmd = r$cmd, stringsAsFactors = FALSE)
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

# Rscript bin/run_dp_future_single.R --workers 4 --cores-per-job 4 -- --DRY_RUN

# Rscript bin/run_dp_future_single.R --workers 4 --cores-per-job 4

# Rscript bin/run_dp_future_single.R --workers 1 --cores-per-job 1
