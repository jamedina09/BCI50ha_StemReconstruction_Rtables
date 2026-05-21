#!/usr/bin/env Rscript
# bb_audit.R — score a stem-reconstruction CSV for broken-below R1/R2 violations.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript bb_audit.R <csv> [<label>]")
csv <- args[[1L]]
label <- if (length(args) >= 2L) args[[2L]] else basename(csv)
d <- fread(csv)
# Use OriginalStemID as series key when StemTag is unreliable; mirror the
# diagnostic in `broken_below_tags.csv`. Default to OriginalStemID.
# StemTag mirrors the diagnostic in broken_below_tags.csv: the per-StemTag
# series can have a BB+DBH row that should split off from the prior row even
# though OriginalStemID has not changed.
key_col <- if ("StemTag" %in% names(d)) "StemTag" else if ("OriginalStemID" %in% names(d)) "OriginalStemID" else "StemID"
setorderv(d, c("Tag", key_col, "CensusID"))
d[, prev_recon := shift(ReconstructedStemID), by = c("Tag", key_col)]
d[, next_recon := shift(ReconstructedStemID, type = "lead"), by = c("Tag", key_col)]
d[, next_dbh   := shift(DBH, type = "lead"), by = c("Tag", key_col)]

bb <- d[Status == "broken below"]
# R1: BB+DBH must split (recon != prev_recon when prev exists)
sub_r1 <- bb[!is.na(DBH) & !is.na(prev_recon)]
sub_r1[, split_ok := ReconstructedStemID != prev_recon]
n_r1_total <- nrow(sub_r1); n_r1_bad <- sum(!sub_r1$split_ok)

# R2: BB+NA-DBH followed by alive (non-NA DBH) row with SAME recon
sub_r2 <- bb[is.na(DBH) & !is.na(next_recon) & !is.na(next_dbh)]
sub_r2[, term_ok := next_recon != ReconstructedStemID]
n_r2_total <- nrow(sub_r2); n_r2_bad <- sum(!sub_r2$term_ok)

cat(sprintf("[%s] tags=%d rows=%d  BB rows=%d (DBH=%d, NA=%d)\n",
            label, uniqueN(d$Tag), nrow(d), nrow(bb),
            sum(!is.na(bb$DBH)), sum(is.na(bb$DBH))))
cat(sprintf("  R1 BB+DBH should split: %d / %d violated (%.1f%%)\n",
            n_r1_bad, n_r1_total, 100*n_r1_bad/max(1,n_r1_total)))
cat(sprintf("  R2 BB+NA→alive same recon: %d / %d violated (%.1f%%)\n",
            n_r2_bad, n_r2_total, 100*n_r2_bad/max(1,n_r2_total)))

cat("\nR1 violations by ReconstructionMethod:\n")
print(sub_r1[(!split_ok), .N, by=ReconstructionMethod])
cat("\nR2 violations by ReconstructionMethod:\n")
print(sub_r2[(!term_ok), .N, by=ReconstructionMethod])

cat("\nMulti-BB at one census (Tag, OriginalStemID, CensusID with >=2 BB rows):\n")
mbb <- bb[, .N, by=c("Tag", key_col, "CensusID")][N>=2L]
cat(sprintf("  occurrences: %d\n", nrow(mbb)))
