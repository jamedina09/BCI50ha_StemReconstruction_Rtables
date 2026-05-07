#!/usr/bin/env bash
# run_bb_sample.sh — run main_cpp_bci.R on a list of tags.
#
# Usage:
#   bash dp_global/scripts/run_bb_sample.sh <out_label> [tags_file]
# Default tags_file is bb_sample_tags.txt at the project root.
# Writes one combined CSV at:
#   dp_global/output/_bb_sample_<out_label>.csv
# and a per-tag run log at:
#   dp_global/output/_bb_sample_<out_label>.log

set -euo pipefail
LABEL="${1:-baseline}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TAGS_FILE="${2:-${ROOT}/bb_sample_tags.txt}"
OUT_DIR="${ROOT}/dp_global/output"
COMBINED="${OUT_DIR}/_bb_sample_${LABEL}.csv"
LOG="${OUT_DIR}/_bb_sample_${LABEL}.log"
TMPDIR_RUN="${OUT_DIR}/_bb_sample_${LABEL}_runs"
mkdir -p "${TMPDIR_RUN}"

: > "${LOG}"
: > "${COMBINED}"
HEADER_WRITTEN=0

i=0
N=$(wc -l < "${TAGS_FILE}" | tr -d ' ')
while IFS= read -r TAG; do
    i=$((i+1))
    [[ -z "${TAG}" ]] && continue
    echo "[$(date +%H:%M:%S)] (${i}/${N}) Tag=${TAG}" | tee -a "${LOG}"
    cd "${ROOT}"
    Rscript dp_global/scripts/main_cpp_bci.R \
        --WHICH_TAG="${TAG}" \
        --DP_VERBOSE=FALSE \
        --WRITE_DP_PDF=FALSE \
        >> "${LOG}" 2>&1 || { echo "  FAILED tag ${TAG}" | tee -a "${LOG}"; continue; }

    # Locate the latest run dir for this tag
    LAST_DIR=$(ls -dt ${OUT_DIR}/*_BCI_tag${TAG}_*_DP_MB_NME_*_rcpp 2>/dev/null | head -1)
    if [[ -z "${LAST_DIR}" || ! -f "${LAST_DIR}/stem_reconstruction_dp_global_rcpp.csv" ]]; then
        echo "  NO OUTPUT for tag ${TAG}" | tee -a "${LOG}"
        continue
    fi
    CSV="${LAST_DIR}/stem_reconstruction_dp_global_rcpp.csv"
    if [[ "${HEADER_WRITTEN}" -eq 0 ]]; then
        cat "${CSV}" >> "${COMBINED}"
        HEADER_WRITTEN=1
    else
        tail -n +2 "${CSV}" >> "${COMBINED}"
    fi
    # Move the run dir under the sample folder so it does not pollute output/
    mv "${LAST_DIR}" "${TMPDIR_RUN}/" || true
done < "${TAGS_FILE}"

echo "[done] combined CSV: ${COMBINED}" | tee -a "${LOG}"
wc -l "${COMBINED}" | tee -a "${LOG}"
