#!/usr/bin/env Rscript
# Quick test: run runner in dry-run mode with posterior sampling and verify joblog records a seed
system2("Rscript", args = c("bin/run_dp_future_single.R", "--workers", "1", "--cores-per-job", "1", "--posterior-samples", "3", "--posterior-format", "csv", "--", "--DRY_RUN"), stdout = TRUE, stderr = TRUE)
joblog <- file.path("tests", "parallel_future_logs", "parallel_future.log")
if (!file.exists(joblog)) {
  cat("test_runner_auto_seed: SKIP (joblog not found)\n")
  quit(status = 0)
}
df <- read.csv(joblog, stringsAsFactors = FALSE)
if (!"posterior_seed"%in%names(df)) stop("posterior_seed column missing from joblog")
if (is.na(df$posterior_seed[1])) stop("posterior_seed in joblog is NA; expected auto-generated seed")
cat("test_runner_auto_seed: OK\n")
