#!/bin/sh
# Rerun of the engine run 20261024_213308_unknown_allT_DP_MB_NME_R12_BD_LT_CU_NS_MS_PO_PG_PF_HP_LF_g5_sm0p5_kg0_ks0_rcpp (started 2026-10-24 21:33:13).
# Run from the project root; the rerun writes a new run folder (new time stamp).
# On macOS, 'caffeinate -i' before Rscript keeps the computer awake.
#
# 1. Code of the run: git commit 1d918b3e62dd20904b170f85ba7dacafb04d6ca4
#      git checkout 1d918b3e62dd20904b170f85ba7dacafb04d6ca4
#    (run_summary.txt lists files that were not tracked by git.)
# 2. Input: BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
#    MD5 16bbf281cfd8fe7048c8b3db0c966930; stage 1 must give the same file.
# 3. Command:
Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES=oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe,ficutr,oreoca \
  --DP_FALLBACK_GROWTH_FORMS=strangler \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE \
  --BASE_OUT_DIR=/Users/medinaja/outputs_bci_stem_identification
#
# The copy of the driver in this folder runs the same way, with the engine
# code (dp_global/) of the same commit:
#   Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/run_records/20261024_213308_unknown_allT_DP_MB_NME_R12_BD_LT_CU_NS_MS_PO_PG_PF_HP_LF_g5_sm0p5_kg0_ks0_rcpp/1_main_cpp_chunk_bci.R \
#     --DP_MAX_STATES=10000 \
#     --PROB_SPECIES=oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe,ficutr,oreoca \
#     --DP_FALLBACK_GROWTH_FORMS=strangler \
#     --POSTERIOR_SAMPLE_SEED=42 \
#     --MANUAL_CORES=TRUE \
#     --MANUAL_CORES_VALUE=18 \
#     --DP_CHUNK_SIZE=18 \
#     --USE_MEASUREMENT_ERROR=FALSE \
#     --BASE_OUT_DIR=/Users/medinaja/outputs_bci_stem_identification
