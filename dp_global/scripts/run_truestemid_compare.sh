#!/usr/bin/env bash
# run_truestemid_compare.sh
#
# Runs the TrueStemID PRE/POST empirical comparison batch.
#
# Usage:
#   bash dp_global/scripts/run_truestemid_compare.sh PRE
#   bash dp_global/scripts/run_truestemid_compare.sh POST
#
# Each invocation:
#   1) Snapshots the current set of output directories under dp_global/output/.
#   2) Runs simulated-data full set (RUN_ALL_TAGS=TRUE) at DP_MAX_STATES=10000
#      and DP_MAX_STATES=2.
#   3) Runs each of the 14 BCI debug tags at DP_MAX_STATES=10000 and 2.
#   4) Identifies newly-created dirs and moves them into
#      dp_global/output/_truestemid_compare/<LABEL>/ where LABEL = PRE or POST.

set -euo pipefail

LABEL="${1:?Usage: $0 PRE|POST}"
if [[ "$LABEL" != "PRE" && "$LABEL" != "POST" ]]; then
    echo "ERROR: LABEL must be PRE or POST (got: $LABEL)" >&2
    exit 1
fi

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"
echo "[runner] PROJECT_ROOT=$PROJECT_ROOT"
echo "[runner] LABEL=$LABEL"

OUT_DIR="dp_global/output"
BUCKET="$OUT_DIR/_truestemid_compare/$LABEL"
mkdir -p "$BUCKET"

LOG_DIR="$BUCKET/_logs"
mkdir -p "$LOG_DIR"

# 14 BCI tags, kept in same order as main_cpp_bci.R header comments
BCI_TAGS=(
    119453 115427 123375 115203 242799 246746 277120
    190932 171486 220311 204785 242114 001080 005558
)

DP_MAX_STATES_LIST=(10000 2)

# Snapshot existing output dirs (excluding the bucket itself) so we can
# detect new directories created by this batch.
SNAPSHOT_FILE="$BUCKET/_pre_run_snapshot.txt"
( cd "$OUT_DIR" && find . -mindepth 1 -maxdepth 1 -type d ! -path './_truestemid_compare' | sort ) > "$SNAPSHOT_FILE"
echo "[runner] Snapshot of existing output dirs: $SNAPSHOT_FILE ($(wc -l < "$SNAPSHOT_FILE") entries)"

run_one() {
    local desc="$1"; shift
    local logfile="$1"; shift
    echo "[runner] >>> $desc"
    echo "[runner]     log: $logfile"
    local t0=$(date +%s)
    if "$@" > "$logfile" 2>&1; then
        local dt=$(( $(date +%s) - t0 ))
        echo "[runner]     OK (${dt}s)"
    else
        local rc=$?
        local dt=$(( $(date +%s) - t0 ))
        echo "[runner]     FAILED rc=$rc after ${dt}s — see $logfile"
        # Continue to next; don't abort the whole batch.
    fi
}

# --- Simulated data: RUN_ALL_TAGS=TRUE, both DP_MAX_STATES values ---
for ms in "${DP_MAX_STATES_LIST[@]}"; do
    run_one "SIM all tags DP_MAX_STATES=$ms" \
        "$LOG_DIR/sim_DP${ms}.log" \
        Rscript dp_global/scripts/main_cpp_chunk.R \
            --RUN_ALL_TAGS=TRUE \
            --DP_MAX_STATES="$ms" \
            --DP_CHUNK_OVERWRITE=TRUE
done

# --- BCI tags: each tag, both DP_MAX_STATES values ---
for tag in "${BCI_TAGS[@]}"; do
    for ms in "${DP_MAX_STATES_LIST[@]}"; do
        run_one "BCI tag=$tag DP_MAX_STATES=$ms" \
            "$LOG_DIR/bci_${tag}_DP${ms}.log" \
            Rscript dp_global/scripts/main_cpp_bci.R \
                --WHICH_TAG="$tag" \
                --DP_MAX_STATES="$ms" \
                --DP_VERBOSE=FALSE
    done
done

# Identify and move new dirs into the bucket.
POST_FILE="$BUCKET/_post_run_listing.txt"
( cd "$OUT_DIR" && find . -mindepth 1 -maxdepth 1 -type d ! -path './_truestemid_compare' | sort ) > "$POST_FILE"

NEW_FILE="$BUCKET/_new_dirs.txt"
comm -13 "$SNAPSHOT_FILE" "$POST_FILE" > "$NEW_FILE" || true

echo "[runner] Moving $(wc -l < "$NEW_FILE") new dirs into $BUCKET"
while IFS= read -r rel; do
    [ -z "$rel" ] && continue
    src="$OUT_DIR/${rel#./}"
    if [ -d "$src" ]; then
        mv "$src" "$BUCKET/"
    fi
done < "$NEW_FILE"

echo "[runner] DONE. Bucket: $BUCKET"
echo "[runner] Contents:"
ls -1 "$BUCKET"
