#!/usr/bin/env Rscript
# compare_truestemid_pre_post.R
#
# Compare PRE-change vs POST-change output buckets and produce a diff report.
# Writes a summary CSV + a markdown report to
# dp_global/output/_truestemid_compare/PRE_POST_DIFF.{csv,md}.
#
# Usage:  Rscript dp_global/scripts/compare_truestemid_pre_post.R

suppressPackageStartupMessages({
    library(data.table)
    library(here)
})

ROOT <- here("dp_global", "output", "_truestemid_compare")
PRE_DIR <- file.path(ROOT, "PRE")
POST_DIR <- file.path(ROOT, "POST")

stopifnot(dir.exists(PRE_DIR), dir.exists(POST_DIR))

CSV_NAME <- "stem_reconstruction_dp_global_rcpp.csv"
CHUNK_GLOB <- "stem_reconstruction_dp_global_rcpp_chunk_*.rds"

# ----- helpers --------------------------------------------------------------

# Strip the leading 14-char timestamp + "_" so PRE and POST runs match by config.
strip_ts <- function(dirname) sub("^[0-9]{8}_[0-9]{6}_", "", dirname)

# Heuristic to extract DP_MAX_STATES suffix used during the run.  We embed it
# into the run-key by checking the runner log filename pattern that produced
# this dir.  Easier: pair PRE and POST by stripped-name (which encodes
# CONFIG_NAME, tag, etc., all the bits affected by build_out_dir_name()) AND
# by which logfile produced them (we stored that mapping nowhere).  Instead we
# detect chunked vs single-tag by file presence and treat each pair as one.
# DP mode is encoded outside the dirname; rely on log-derived info: we'll
# read DP_FallbackReason and DP_MaxStatesPerCensus from the data itself.

load_run_data <- function(run_dir) {
    csv_p <- file.path(run_dir, CSV_NAME)
    if (file.exists(csv_p)) {
        d <- fread(csv_p)
        return(d)
    }
    # chunked run: rbind all chunk RDS
    chunks <- Sys.glob(file.path(run_dir, CHUNK_GLOB))
    if (length(chunks) == 0L) return(NULL)
    parts <- lapply(chunks, function(f) {
        x <- readRDS(f)
        if (is.list(x) && !is.data.frame(x)) {
            # list of group results
            x <- rbindlist(x, use.names = TRUE, fill = TRUE)
        }
        as.data.table(x)
    })
    rbindlist(parts, use.names = TRUE, fill = TRUE)
}

# ----- pair runs -----------------------------------------------------------

pre_runs  <- list.dirs(PRE_DIR,  recursive = FALSE)
post_runs <- list.dirs(POST_DIR, recursive = FALSE)
pre_runs  <- pre_runs[!grepl("^_",  basename(pre_runs))]
post_runs <- post_runs[!grepl("^_", basename(post_runs))]

cat("PRE  runs:", length(pre_runs), "\n")
cat("POST runs:", length(post_runs), "\n")

pre_keys  <- strip_ts(basename(pre_runs))
post_keys <- strip_ts(basename(post_runs))

# A run-key may appear twice in PRE (DP_MAX_STATES=10000 vs =2 produce the
# same dirname suffix).  We must distinguish via the per-log mapping stored
# in `_logs/`.  Use lexicographic timestamp order to assign DP10000 first
# and DP2 second within each key, mirroring the runner script's order.
order_within_key <- function(paths) {
    paths[order(basename(paths))]
}

pre_by_key  <- split(pre_runs,  pre_keys)
post_by_key <- split(post_runs, post_keys)

per_run_summary <- list()
combined_long   <- list()

for (key in intersect(names(pre_by_key), names(post_by_key))) {
    p_paths <- order_within_key(pre_by_key[[key]])
    q_paths <- order_within_key(post_by_key[[key]])
    n_pair  <- min(length(p_paths), length(q_paths))
    for (i in seq_len(n_pair)) {
        pre  <- p_paths[i]
        post <- q_paths[i]
        cat(sprintf("Pairing [%d] %s\n  PRE : %s\n  POST: %s\n",
                    i, key, basename(pre), basename(post)))

        d_pre  <- load_run_data(pre)
        d_post <- load_run_data(post)
        if (is.null(d_pre) || is.null(d_post)) {
            cat("  (skip — could not load data)\n"); next
        }

        # Restrict to columns we always have
        keep <- intersect(c("Tag","CensusID","OriginalStemID","TrueStemID",
                            "DBH","Status","ListOfTSM","StemTag",
                            "ReconstructedStemID","ReconstructionMethod",
                            "DP_FallbackReason","DP_KUsed","DP_MaxStatesPerCensus"),
                          names(d_pre))
        d_pre  <- d_pre[, ..keep]
        d_post <- d_post[, intersect(keep, names(d_post)), with = FALSE]

        # Per-run summary
        rk <- paste0(key, "_pair", i)
        per_run_summary[[rk]] <- data.table(
            run_key                            = rk,
            n_rows_pre                         = nrow(d_pre),
            n_rows_post                        = nrow(d_post),
            n_truestemid_pre                   = sum(!is.na(d_pre$TrueStemID)),
            n_truestemid_post                  = sum(!is.na(d_post$TrueStemID)),
            n_violations_pre                   = sum(!is.na(d_pre$TrueStemID) &
                                                     (is.na(d_pre$ReconstructedStemID) |
                                                      d_pre$TrueStemID != d_pre$ReconstructedStemID)),
            n_violations_post                  = sum(!is.na(d_post$TrueStemID) &
                                                     (is.na(d_post$ReconstructedStemID) |
                                                      d_post$TrueStemID != d_post$ReconstructedStemID)),
            n_recon_na_pre                     = sum(is.na(d_pre$ReconstructedStemID)),
            n_recon_na_post                    = sum(is.na(d_post$ReconstructedStemID)),
            method_given_pre                   = sum(d_pre$ReconstructionMethod == "given", na.rm = TRUE),
            method_given_post                  = sum(d_post$ReconstructionMethod == "given", na.rm = TRUE),
            method_dp_pre                      = sum(d_pre$ReconstructionMethod == "dp", na.rm = TRUE),
            method_dp_post                     = sum(d_post$ReconstructionMethod == "dp", na.rm = TRUE),
            method_prob_pre                    = sum(d_pre$ReconstructionMethod == "probabilistic", na.rm = TRUE),
            method_prob_post                   = sum(d_post$ReconstructionMethod == "probabilistic", na.rm = TRUE),
            fallback_count_pre                 = if ("DP_FallbackReason" %in% names(d_pre))  sum(!is.na(d_pre$DP_FallbackReason))  else NA_integer_,
            fallback_count_post                = if ("DP_FallbackReason" %in% names(d_post)) sum(!is.na(d_post$DP_FallbackReason)) else NA_integer_
        )

        # Per-row diff (align by Tag+CensusID+OriginalStemID)
        join_keys <- intersect(c("Tag","CensusID","OriginalStemID"), names(d_pre))
        m <- merge(d_pre,  d_post,  by = join_keys, suffixes = c("_pre","_post"), all = TRUE)
        m[, recon_changed := !((is.na(ReconstructedStemID_pre) & is.na(ReconstructedStemID_post)) |
                               (!is.na(ReconstructedStemID_pre) & !is.na(ReconstructedStemID_post) &
                                ReconstructedStemID_pre == ReconstructedStemID_post))]
        m[, method_changed := !((is.na(ReconstructionMethod_pre) & is.na(ReconstructionMethod_post)) |
                                (!is.na(ReconstructionMethod_pre) & !is.na(ReconstructionMethod_post) &
                                 ReconstructionMethod_pre == ReconstructionMethod_post))]
        per_run_summary[[rk]][, n_recon_changed := sum(m$recon_changed)]
        per_run_summary[[rk]][, n_method_changed := sum(m$method_changed)]
        m[, run_key := rk]
        combined_long[[rk]] <- m[recon_changed | method_changed,
                                 c(join_keys, "TrueStemID_pre","TrueStemID_post",
                                   "ReconstructedStemID_pre","ReconstructedStemID_post",
                                   "ReconstructionMethod_pre","ReconstructionMethod_post",
                                   "run_key"), with = FALSE]
    }
}

summary_dt <- rbindlist(per_run_summary, use.names = TRUE, fill = TRUE)
diff_dt    <- rbindlist(combined_long,   use.names = TRUE, fill = TRUE)

out_summary <- file.path(ROOT, "PRE_POST_SUMMARY.csv")
out_diff    <- file.path(ROOT, "PRE_POST_DIFF_ROWS.csv")
fwrite(summary_dt, out_summary)
fwrite(diff_dt,    out_diff)

cat("\n--- SUMMARY ---\n")
print(summary_dt)
cat("\nWrote:\n  ", out_summary, "\n  ", out_diff, "\n", sep = "")

# Aggregate invariant check
cat("\n--- INVARIANT TOTALS (TrueStemID != ReconstructedStemID counts) ---\n")
cat("PRE  total violations:", sum(summary_dt$n_violations_pre), "\n")
cat("POST total violations:", sum(summary_dt$n_violations_post), "\n")
