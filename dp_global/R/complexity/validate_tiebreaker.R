############################################################
### validate_tiebreaker.R
### Runs the DP twice on simulated_data_1.csv (single core,
### deterministic) and checks that results are 100% identical.
### Also reports per-tag accuracy vs TrueStemID.
### Run from project root:
###   Rscript dp_global/R/complexity/validate_tiebreaker.R
############################################################
suppressPackageStartupMessages({
    library(here)
    library(data.table)
})

OUT1 <- here("dp_global", "R", "complexity", "_val_run1.csv")
OUT2 <- here("dp_global", "R", "complexity", "_val_run2.csv")

run_dp <- function(out_file) {
    cmd <- paste(
        "Rscript dp_global/scripts/main_cpp.R",
        "--RUN_ALL_TAGS=TRUE",
        "--WRITE_DP_CSV=TRUE --WRITE_DP_RDS=FALSE --WRITE_DP_PDF=FALSE",
        "--POSTERIOR_SAMPLES=0",
        "--MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=1",
        "2>&1"
    )
    message("\n=== Running DP ===\n", cmd)
    system(cmd)

    # Find the most recent output CSV
    dirs <- list.dirs(here("dp_global", "output"), recursive = FALSE)
    dirs <- sort(dirs, decreasing = TRUE)
    for (d in dirs) {
        cand <- file.path(d, "stem_reconstruction_dp_global_rcpp.csv")
        if (file.exists(cand)) {
            file.copy(cand, out_file, overwrite = TRUE)
            message("  Saved -> ", out_file)
            return(invisible(NULL))
        }
    }
    stop("Could not find DP output CSV")
}

run_dp(OUT1)
run_dp(OUT2)

r1 <- fread(OUT1)
r2 <- fread(OUT2)

# ----------------------------------------------------------------
# 1. Determinism: run1 == run2 row-by-row on ReconstructedStemID
# ----------------------------------------------------------------
key_cols <- c("Tag", "CensusID", "OriginalStemID")
setorderv(r1, key_cols)
setorderv(r2, key_cols)

same_ids <- identical(r1$ReconstructedStemID, r2$ReconstructedStemID)
n_diff <- sum(r1$ReconstructedStemID != r2$ReconstructedStemID, na.rm = TRUE) +
          sum(xor(is.na(r1$ReconstructedStemID), is.na(r2$ReconstructedStemID)))

cat("\n================================================\n")
cat("  DETERMINISM CHECK (run1 vs run2)\n")
cat("================================================\n")
cat(sprintf("  Total rows          : %d\n", nrow(r1)))
cat(sprintf("  Rows with same ID   : %d\n", nrow(r1) - n_diff))
cat(sprintf("  Rows with diff ID   : %d\n", n_diff))
cat(sprintf("  Identical?          : %s\n", if (same_ids) "YES ✓" else "NO ✗"))

if (n_diff > 0L) {
    diffs <- r1[r1$ReconstructedStemID != r2$ReconstructedStemID, .(Tag, CensusID, OriginalStemID,
        run1 = ReconstructedStemID)]
    diffs[, run2 := r2[r1$ReconstructedStemID != r2$ReconstructedStemID, ReconstructedStemID]]
    cat("\n  First 20 differences:\n")
    print(head(diffs, 20))
}

# ----------------------------------------------------------------
# 2. Accuracy vs TrueStemID (using run1)
# ----------------------------------------------------------------
if ("TrueStemID" %in% names(r1) && "ReconstructedStemID" %in% names(r1)) {
    acc_rows <- r1[!is.na(TrueStemID) & !is.na(ReconstructedStemID)]
    overall  <- mean(acc_rows$TrueStemID == acc_rows$ReconstructedStemID)

    # Separate M-coded test tags (901, 902, 903) from baseline
    m_tags   <- c("901", "902", "903")
    acc_m    <- acc_rows[Tag %in% m_tags]
    acc_base <- acc_rows[!(Tag %in% m_tags)]

    cat("\n================================================\n")
    cat("  ACCURACY vs TrueStemID\n")
    cat("================================================\n")
    cat(sprintf("  Overall accuracy    : %.1f%%  (%d / %d rows)\n",
        100 * overall, sum(acc_rows$TrueStemID == acc_rows$ReconstructedStemID), nrow(acc_rows)))

    if (nrow(acc_base) > 0L) {
        acc_b <- mean(acc_base$TrueStemID == acc_base$ReconstructedStemID)
        cat(sprintf("  Non-M tags          : %.1f%%  (%d rows)\n",
            100 * acc_b, nrow(acc_base)))
    }
    if (nrow(acc_m) > 0L) {
        acc_mv <- mean(acc_m$TrueStemID == acc_m$ReconstructedStemID)
        cat(sprintf("  M-coded tags (901-903): %.1f%%  (%d rows)\n",
            100 * acc_mv, nrow(acc_m)))
        # Per-tag breakdown
        per_tag <- acc_m[, .(correct = sum(TrueStemID == ReconstructedStemID), N = .N), by = Tag]
        per_tag[, pct := round(100 * correct / N, 1)]
        print(per_tag)
    }
}

cat("\n================================================\n")
cat("  Done.\n")
cat("================================================\n")
