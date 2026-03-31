#!/usr/bin/env Rscript
# Validation script: M-code implementation
# Phase 2 — consistency on non-M tags (baseline vs M-aware run)
# Phase 3 — M-tag accuracy
# Run from project root:
#   Rscript dp_global/R/complexity/validate_m_implementation.R

library(data.table)

# Two deterministic single-core runs (MANUAL_CORES_VALUE=1) on the same data with M tags present.
# These must be identical — any diff is a real non-determinism bug, not parallel tie-breaking.
BASELINE_DIR <- "dp_global/output/20260331_122650_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp"
NEW_DIR      <- "dp_global/output/20260331_122808_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp"

b <- fread(file.path(BASELINE_DIR, "stem_reconstruction_dp_global_rcpp.csv"))
n <- fread(file.path(NEW_DIR,      "stem_reconstruction_dp_global_rcpp.csv"))

# ==========================================================================
# Phase 2: consistency check on non-M tags
# Expected: 100% match between baseline (pre-M) and M-aware run
# ==========================================================================
shared_tags <- intersect(unique(b$Tag), unique(n$Tag))
b2 <- b[Tag %in% shared_tags]
n2 <- n[Tag %in% shared_tags]
setkey(b2, Tag, CensusID, DBH)
setkey(n2, Tag, CensusID, DBH)
merged <- merge(
    b2[, .(Tag, CensusID, DBH, base = ReconstructedStemID)],
    n2[, .(Tag, CensusID, DBH, new  = ReconstructedStemID)],
    by = c("Tag", "CensusID", "DBH")
)
match_rate <- mean(
    merged$base == merged$new | (is.na(merged$base) & is.na(merged$new)),
    na.rm = TRUE
)
cat("\n=== Phase 2: Consistency on non-M tags ===\n")
cat(sprintf("  Baseline run  : %s\n", BASELINE_DIR))
cat(sprintf("  New run       : %s\n", NEW_DIR))
cat(sprintf("  Shared tags   : %d\n", length(shared_tags)))
cat(sprintf("  Total rows    : %d\n", nrow(merged)))
cat(sprintf("  Match rate    : %.6f  [%s]\n\n",
    match_rate, if (match_rate == 1) "PASS" else "FAIL"))
if (match_rate < 1) {
    cat("Mismatching rows:\n")
    print(merged[base != new | (is.na(base) != is.na(new))])
}

# ==========================================================================
# Phase 3: M-tag accuracy (pre-anchor rows only)
# Ground truth = OriginalStemID (never used by the algorithm).
# ==========================================================================
cat("=== Phase 3: M-tag accuracy ===\n")
m_tags <- c(901L, 902L, 903L)
res_m  <- n[Tag %in% m_tags]
pre    <- res_m[is.na(TrueStemID)]   # pre-anchor rows; algorithm had to reconstruct these
pre[, correct := (ReconstructedStemID == OriginalStemID)]
acc_overall <- mean(pre$correct, na.rm = TRUE)
cat(sprintf("  M-tag pre-anchor rows: %d\n", nrow(pre)))
cat(sprintf("  Overall accuracy     : %.6f  [%s]\n\n",
    acc_overall, if (acc_overall == 1) "PASS" else "FAIL"))

for (tg in m_tags) {
    sub <- pre[Tag == tg]
    acc <- mean(sub$correct, na.rm = TRUE)
    cat(sprintf("  Tag %d (%d rows): %.6f  [%s]\n",
        tg, nrow(sub), acc, if (acc == 1) "PASS" else "FAIL"))
}
cat("\n")

fails <- pre[correct == FALSE]
if (nrow(fails) > 0) {
    cat("FAILED rows:\n")
    print(fails[, .(Tag, CensusID, DBH, OriginalStemID, ReconstructedStemID)])
} else {
    cat("All M-tagged pre-anchor rows reconstructed correctly.\n")
}

# ==========================================================================
# Summary: Overall accuracy on all pre-anchor rows in new run
# ==========================================================================
cat("\n=== Overall accuracy (new run, all pre-anchor rows) ===\n")
all_pre <- n[is.na(TrueStemID) & !is.na(OriginalStemID)]
all_pre[, correct := (ReconstructedStemID == OriginalStemID)]
cat(sprintf("  Tags : %d\n", uniqueN(all_pre$Tag)))
cat(sprintf("  Rows : %d\n", nrow(all_pre)))
cat(sprintf("  Accuracy : %.6f\n\n", mean(all_pre$correct, na.rm = TRUE)))

# Worst tags
worst <- all_pre[, .(acc = mean(correct), n = .N), by = Tag][order(acc)]
cat("Bottom 10 tags by accuracy:\n")
print(head(worst, 10))
