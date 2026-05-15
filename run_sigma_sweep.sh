#!/bin/bash
set -euo pipefail

COMMON="--POSTERIOR_SAMPLES=200 --WHICH_TAG=255814 --DP_MAX_STATES=10000 \
  --PROB_SPECIES=oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe \
  --DP_FALLBACK_GROWTH_FORMS=strangler_fig \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=16 --DP_CHUNK_SIZE=16 \
  --USE_MEASUREMENT_ERROR=FALSE"

for sigma in 0 0.5 1.0 1.5 2.0 2.5 3.0; do
  echo "=== PROB_N_SIGMA_ME=${sigma} ==="
  Rscript dp_global/scripts/main_cpp_bci.R \
    ${COMMON} \
    --PROB_N_SIGMA_ME=${sigma} \
    2>&1 | grep -E "Sample-level repair|DP done"

  latest_csv=$(ls -t dp_global/output/*/stem_reconstruction_dp_global_rcpp.csv 2>/dev/null | head -1)
  if [ -n "$latest_csv" ]; then
    Rscript dp_global/scripts/_sigma_stats.R "$latest_csv" "$sigma" 2>/dev/null
  fi
done
