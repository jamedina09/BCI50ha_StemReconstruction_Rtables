#!/usr/bin/env bash
set -euo pipefail

echo "Starting run_dp_full.sh with args: $@"

###############################################################################
# Batch-level timestamp (ONE per invocation)
# If BATCH_TS is already exported in the environment, preserve it so multiple
# runs can share the same timestamp.
###############################################################################
if [[ -z "${BATCH_TS:-}" ]]; then
  BATCH_TS=$(date +"%Y%m%d_%H%M%S")
  export BATCH_TS
  echo "Batch timestamp (generated): $BATCH_TS"
else
  echo "Batch timestamp (from environment): $BATCH_TS"
fi

###############################################################################
# Parse --config if provided (single config run)
###############################################################################
CONFIG=""
ARGS=()

for arg in "$@"; do
  if [[ $arg == --config=* ]]; then
    CONFIG="${arg#--config=}"
  else
    ARGS+=("$arg")
  fi
done

# Detect DRY_RUN flag in extra args or env var. When DRY_RUN=1, we print Rscript invocations
# instead of executing them. This is useful for CI/smoke tests and quick local checks.
DRY_RUN=0
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  DRY_RUN=1
fi
for a in "${ARGS[@]}"; do
  if [[ $a == --DRY_RUN ]]; then
    DRY_RUN=1
    echo "DRY_RUN enabled; will not execute Rscript calls for these configs."
    break
  fi
done

###############################################################################
# Script directory
###############################################################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "SCRIPT_DIR: $SCRIPT_DIR"
# Project root (one level up from this script): repository root
PROJ_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
echo "PROJ_ROOT: $PROJ_ROOT"

###############################################################################
# Base arguments (shared across all configs)
###############################################################################
BASE_ARGS=(
  "--input_file=XXX"
  "--FORCE_ONE_SPECIES_PARAMETERS=FALSE"

  "--DP_MODE=marginals+bins"
  "--which_tag=1"
  "--anchor_start_census=7"
  "--census_interval_years=5"
  "--DP_VERBOSE=TRUE"

  "--RUN_ALL_TAGS=TRUE"
  # Allow overriding DP enumerator state budget (see dp_max_states in main_cpp.R)
  "--dp_max_states=40000"
  "--MANUAL_CORES=TRUE"
  "--MANUAL_CORES_VALUE=15"

  "--WRITE_DP_CSV=TRUE"
  "--WRITE_DP_PDF=TRUE"
  "--DP_PDF_INCLUDE_REFERENCE=TRUE"
  "--PLOT_PDF_ONE_TAG_ONLY=FALSE"

  "--SENSITIVITY_MODE=run+write+pdf"
  "--RUN_K_SWEEP_DEMO=TRUE"

  "--RUN_REALISM_REPORT=TRUE"

  "--RECRUIT_MAX_SOURCE=data",
  "--RECRUIT_MAX_FIXED=5",

  "--PROJECT_ROOT=${PROJ_ROOT}"

  "--BATCH_TS=${BATCH_TS}" # ensure all R runs share the same timestamp
) 

# Launch the R driver (use project-root-relative path)

###############################################################################
# Config list
###############################################################################
configs=(
  fixed
  data_hard
  data_hard_soft
  data_soft
  fixed_k50
  fixed_k25
  data_hard_k50
  data_hard_k25
)

if [[ -n "$CONFIG" ]]; then
  configs=("$CONFIG")
  echo "Running single config: $CONFIG"
else
  echo "Running all configs"
fi

echo "Configs to run: ${configs[*]}"

# Determine requested cores per job (from ARGS override, env, or default)
# To change per-job cores, pass --MANUAL_CORES_VALUE=<N> (default 15)

CORES_PER_JOB="${MANUAL_CORES_VALUE:-}"
for a in "${ARGS[@]}"; do
  if [[ $a == --MANUAL_CORES_VALUE=* ]]; then
    CORES_PER_JOB="${a#--MANUAL_CORES_VALUE=}"
  fi
done
if [[ -z "${CORES_PER_JOB:-}" ]]; then
  CORES_PER_JOB=15
fi

###############################################################################
# Main loop
###############################################################################
for CONFIG in "${configs[@]}"; do
  echo "------------------------------------------------------------"
  echo "Running config: $CONFIG"

  case "$CONFIG" in
  fixed)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=fixed"
      "--MAX_GROWTH_FIXED=7.5"
      "--MAX_SHRINK_HARD_SOURCE=fixed"
      "--MAX_SHRINK_FIXED=-0.5"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=0"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=0"
      "--RECRUIT_MAX_SOURCE=fixed"
      "--RECRUIT_MAX_FIXED=6"
    )
    ;;
  data_hard)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=data"
      "--MAX_SHRINK_HARD_SOURCE=data"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=0"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=0"
    )
    ;;
  data_hard_soft)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=data"
      "--MAX_SHRINK_HARD_SOURCE=data"
      "--K_SHRINK_SOURCE=data"
      "--K_GROWTH_SOURCE=data"
    )
    ;;
  data_soft)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=fixed"
      "--MAX_GROWTH_FIXED=7.5"
      "--MAX_SHRINK_HARD_SOURCE=fixed"
      "--MAX_SHRINK_FIXED=-0.5"
      "--K_SHRINK_SOURCE=data"
      "--K_GROWTH_SOURCE=data"
    )
    ;;
  fixed_k50)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=fixed"
      "--MAX_GROWTH_FIXED=7.5"
      "--MAX_SHRINK_HARD_SOURCE=fixed"
      "--MAX_SHRINK_FIXED=-0.5"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=50"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=50"
    )
    ;;
  fixed_k25)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=fixed"
      "--MAX_GROWTH_FIXED=7.5"
      "--MAX_SHRINK_HARD_SOURCE=fixed"
      "--MAX_SHRINK_FIXED=-0.5"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=25"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=25"
    )
    ;;
  data_hard_k50)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=data"
      "--MAX_SHRINK_HARD_SOURCE=data"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=50"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=50"
    )
    ;;
  data_hard_k25)
    CONFIG_ARGS=(
      "--USE_MEASUREMENT_ERROR=TRUE"
      "--MAX_GROWTH_HARD_SOURCE=data"
      "--MAX_SHRINK_HARD_SOURCE=data"
      "--K_SHRINK_SOURCE=fixed"
      "--K_SHRINK_FIXED=25"
      "--K_GROWTH_SOURCE=fixed"
      "--K_GROWTH_FIXED=25"
    )
    ;;
  *)
    echo "Unknown config: $CONFIG"
    continue
    ;;
  esac

  # Always pass current config name to R for directory naming
  CONFIG_ARGS+=("--CONFIG_NAME=$CONFIG")

  # Combine base + config + extra CLI args
  ALL_ARGS=(
    "${BASE_ARGS[@]}"
    "${CONFIG_ARGS[@]}"
    "${ARGS[@]}"
  )

  echo "Launching Rscript for $CONFIG"
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    # Print a project-relative invocation with sanitized args (replace PROJ_ROOT with '.')
    SAN_ARGS=()
    for a in "${ALL_ARGS[@]}"; do
      s="${a//$PROJ_ROOT/.}"
      SAN_ARGS+=("$s")
    done
    echo "DRY RUN: Rscript \"dp_global/scripts/main.R\" \"${SAN_ARGS[@]}\""
  else
    Rscript "$PROJ_ROOT/dp_global/scripts/main.R" "${ALL_ARGS[@]}"
  fi

  echo "Completed config: $CONFIG"
done

echo "============================================================"
echo "All configs completed successfully."

## quick run
# Run a single config in DRY_RUN mode:
#   ./bin/run_dp_full.sh --config=fixed --SENSITIVITY_MODE=none --DP_MODE=none --DRY_RUN

# Run all configs serially (default behavior):
#   ./bin/run_dp_full.sh
# You may change per-job cores with --MANUAL_CORES_VALUE=<N> (default 15).

# Notes:
# - This script runs configs serially in a loop and is intended for single-machine use.
# - If you need concurrency, run concurrent invocations manually and ensure you set a shared BATCH_TS and handle resource coordination yourself.
