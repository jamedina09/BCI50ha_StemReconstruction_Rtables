#!/usr/bin/env Rscript
# _sigma_stats.R  <csv_path> <sigma_value>
args <- commandArgs(trailingOnly = TRUE)
csv_path <- args[1]
sigma    <- args[2]

library(data.table)
x <- fread(csv_path)

prob_rows <- x[ReconstructionMethod == "probabilistic" & !is.na(ReconstructedStemID)]

if (nrow(prob_rows) == 0L) {
    cat(sprintf("  n_sigma_me=%s | no probabilistic rows\n", sigma))
    quit(status = 0)
}

setorder(prob_rows, ReconstructedStemID, CensusID)

# Compute annualised growth rate per consecutive pair within each track.
# Interval is in census units (1 unit ≈ 5 years); DBH in cm.
prob_rows[, prev_DBH := shift(DBH, 1L, type = "lag"), by = ReconstructedStemID]
prob_rows[, prev_CensusID := shift(CensusID, 1L, type = "lag"), by = ReconstructedStemID]
prob_rows <- prob_rows[!is.na(prev_DBH) & !is.na(DBH) & (CensusID - prev_CensusID) > 0]

# Use census gap * 5 as interval in years (BCI standard)
prob_rows[, interval_yr := (CensusID - prev_CensusID) * 5]
prob_rows[, growth_rate := (DBH - prev_DBH) / interval_yr]

neg <- prob_rows[growth_rate < 0]
n_links <- nrow(prob_rows)

cat(sprintf(
    "  n_sigma_me=%s | total_links=%d | neg_links=%d (%.1f%%) | median_neg=%.4f | min_neg=%.4f\n",
    sigma,
    n_links,
    nrow(neg),
    if (n_links > 0) 100 * nrow(neg) / n_links else 0,
    if (nrow(neg) > 0) median(neg$growth_rate) else 0,
    if (nrow(neg) > 0) min(neg$growth_rate) else 0
))
