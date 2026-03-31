library(data.table)

pre  <- fread("dp_global/output/20260331_181952_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv")
post <- fread("dp_global/output/20260331_183123_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv")
cat("pre rows:", nrow(pre), "  post rows:", nrow(post), "\n")

# Compare all non-M tags via Tag + CensusID + TrueStemID
pre_nm  <- pre[!Tag  %in% c(901L, 902L, 903L)]
post_nm <- post[!Tag %in% c(901L, 902L, 903L)]
# Use DBH as additional key to disambiguate rows with same CensusID but no TrueStemID
m <- merge(pre_nm[, .(Tag, CensusID, TrueStemID, DBH, pre = ReconstructedStemID)],
           post_nm[, .(Tag, CensusID, TrueStemID, DBH, post = ReconstructedStemID)],
           by = c("Tag", "CensusID", "TrueStemID", "DBH"), allow.cartesian = FALSE)
cat("Merged rows:", nrow(m), "\n")
diffs <- m[pre != post | (is.na(pre) != is.na(post))]
cat("Differing rows:", nrow(diffs), "\n")
if (nrow(diffs) > 0L) print(diffs[seq_len(min(10L, nrow(diffs)))])

# M-tag accuracy
m_tags <- post[Tag %in% c(901L, 902L, 903L) & !is.na(OriginalStemID)]
acc <- mean(m_tags$OriginalStemID == m_tags$ReconstructedStemID, na.rm = TRUE)
cat(sprintf("M-tag accuracy: %.4f  (%d rows)\n", acc, nrow(m_tags)))

# Placeholder so the rest of the old script doesn't run
q("no")

# Non-M regression — join on obs_row_id (stable across runs)
new_nm  <- new[!Tag  %in% c(901L, 902L, 903L)]
base_nm <- base[!Tag %in% c(901L, 902L, 903L)]
m <- merge(base_nm[, .(obs_row_id, base = ReconstructedStemID)],
           new_nm[,  .(obs_row_id, new  = ReconstructedStemID)],
           by = "obs_row_id")
match_rate <- mean(m$base == m$new | (is.na(m$base) & is.na(m$new)), na.rm = TRUE)
cat(sprintf("Non-M regression match: %.4f  (%d rows differ / %d total)\n",
            match_rate, sum(m$base != m$new, na.rm = TRUE), nrow(m)))

# M-tag accuracy (OriginalStemID = ground truth in simulated data)
m_tags <- new[Tag %in% c(901L, 902L, 903L) & !is.na(OriginalStemID)]
acc <- mean(m_tags$OriginalStemID == m_tags$ReconstructedStemID, na.rm = TRUE)
cat(sprintf("M-tag accuracy:         %.4f  (%d rows)\n", acc, nrow(m_tags)))
