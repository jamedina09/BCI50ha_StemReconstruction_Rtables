#!/usr/bin/env Rscript
library(data.table)

b <- fread("dp_global/output/20260331_122650_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv")
n <- fread("dp_global/output/20260331_122808_unknown_allT_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp.csv")

shared_tags <- intersect(unique(b$Tag), unique(n$Tag))
b2 <- b[Tag %in% shared_tags]
n2 <- n[Tag %in% shared_tags]
setkey(b2, Tag, CensusID)
setkey(n2, Tag, CensusID)

mismatched_tags <- c(2,3,4,8,11,16,20,21,26,30,37,38,61350,123315,229214)

cat("Checking whether Phase 2 mismatches are pure symmetric swaps...\n")
set_match_results <- sapply(mismatched_tags, function(tg) {
    bc <- b2[Tag == tg, .(CensusID, v = ReconstructedStemID)]
    nc <- n2[Tag == tg, .(CensusID, v = ReconstructedStemID)]
    all(sapply(unique(c(bc$CensusID, nc$CensusID)), function(cc) {
        setequal(sort(bc[CensusID == cc, v]), sort(nc[CensusID == cc, v]))
    }))
})
names(set_match_results) <- mismatched_tags
print(set_match_results)
cat(sprintf("\nAll mismatches are pure symmetric swaps: %s\n",
    if (all(set_match_results)) "YES -> Phase 2 PASS (non-determinism only, no regression)"
    else "NO -> real regression detected"))
