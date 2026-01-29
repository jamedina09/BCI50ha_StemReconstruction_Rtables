library(data.table)
chunk8 <- readRDS('dp_global/output/20260129_122104_unknown_T0_DP_MB_ME_g7p5_sm0p5_kg0_ks0_rcpp/stem_reconstruction_dp_global_rcpp_chunk_008.rds')
sel <- chunk8[Tag %in% c(51,52,53), .(Tag, CensusID, DBH, TrueStemID, ReconstructedStemID, ReconstructionMethod, ConstraintViolation, DP_PosteriorReconstructedProb)]
print(sel)
