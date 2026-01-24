#!/usr/bin/env bash
set -euo pipefail

# package_bundle.sh
# Create a minimal timestamped tarball containing everything needed to run DP
# - auto-generates INSTALL.txt into the archive
# - includes dpglobal_bundle minimal files + files used by scripts/main_cpp.R
# - writes a checksum (.sha256) next to the tarball
# - cleans up temporary staging area

# Usage:
# From project root or anywhere: dp_global/R/dpglobal_bundle/package_bundle.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
OUT_DIR="${SCRIPT_DIR}/dist"
mkdir -p "${OUT_DIR}"
TS="$(date +%Y%m%d_%H%M%S)"
TARBALL_NAME="dpglobal_bundle_full_${TS}.tar.gz"
TARBALL_PATH="${OUT_DIR}/${TARBALL_NAME}"
SHA_PATH="${TARBALL_PATH}.sha256"
STAGE_DIR="$(mktemp -d)"  # will be removed on exit

# CLI flags
BUILD_BUNDLE=false
for arg in "$@"; do
  case "$arg" in
    --build-bundle) BUILD_BUNDLE=true ;;
    --help|-h) echo "Usage: package_bundle.sh [--build-bundle]"; exit 0 ;;
    *) echo "[package_bundle] Unknown arg: $arg" ;;
  esac
done

# Optionally build dpglobal_bundle.RData before packing
if [ "$BUILD_BUNDLE" = true ]; then
  echo "[package_bundle] Building dpglobal_bundle.RData using dpglobal_bundle_loader.R"
  Rscript -e "source('${PROJECT_ROOT}/dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"
fi

cleanup() {
  rm -rf "${STAGE_DIR}"
}
trap cleanup EXIT

echo "[package_bundle] Staging to ${STAGE_DIR}"

# Minimal bundle files (bundle-local)
if [ -f "${SCRIPT_DIR}/dpglobal_bundle.RData" ]; then
  cp -v "${SCRIPT_DIR}/dpglobal_bundle.RData" "${STAGE_DIR}/"
else
  echo "[package_bundle] NOTE: dpglobal_bundle.RData not present in ${SCRIPT_DIR}; it will not be included in the package. Run dpglobal_bundle_loader.R to create it if desired."
fi
if [ -f "${SCRIPT_DIR}/dpglobal_bundle_manifest.rds" ]; then
  cp -v "${SCRIPT_DIR}/dpglobal_bundle_manifest.rds" "${STAGE_DIR}/"
else
  echo "[package_bundle] NOTE: dpglobal_bundle_manifest.rds not present in ${SCRIPT_DIR}; it will not be included in the package."
fi
cp -v "${SCRIPT_DIR}/README.md" "${STAGE_DIR}/"
# copy verify script if present (optional helper)
if [ -f "${SCRIPT_DIR}/verify_bundle.R" ]; then cp -v "${SCRIPT_DIR}/verify_bundle.R" "${STAGE_DIR}/"; fi

# Copy R wrapper and C++ source from the primary dp_global/src location (keep bundle lean)
if [ -f "${PROJECT_ROOT}/dp_global/src/transition_cost_rcpp.R" ]; then
  # copy wrapper both to archive root (convenience) and to the bundle subfolder used by README
  cp -v "${PROJECT_ROOT}/dp_global/src/transition_cost_rcpp.R" "${STAGE_DIR}/"
  cp -v "${PROJECT_ROOT}/dp_global/src/transition_cost_rcpp.R" "${STAGE_DIR}/dp_global/R/dpglobal_bundle/transition_cost_rcpp.R" || true
else
  echo "[package_bundle] NOTE: transition_cost_rcpp.R not found in dp_global/src; skipping"
fi
mkdir -p "${STAGE_DIR}/dp_global/R/dpglobal_bundle/src"
# ensure the bundle subdir exists before copying wrapper
mkdir -p "${STAGE_DIR}/dp_global/R/dpglobal_bundle"
if [ -f "${PROJECT_ROOT}/dp_global/src/transition_cost_rcpp.cpp" ]; then
  cp -v "${PROJECT_ROOT}/dp_global/src/transition_cost_rcpp.cpp" "${STAGE_DIR}/dp_global/R/dpglobal_bundle/src/"
else
  echo "[package_bundle] NOTE: transition_cost_rcpp.cpp not found in dp_global/src; skipping"
fi

# Copy the full dp_global project tree required to run scripts/main_cpp.R
# Exclude large runtime outputs and generated bundle dist files to keep tarball reasonable.
if command -v rsync >/dev/null 2>&1; then
  echo "[package_bundle] Using rsync to copy dp_global (excluding output and bundle dist)"
  rsync -av --exclude 'output' --exclude 'R/dpglobal_bundle/dist' --exclude '.git' "${PROJECT_ROOT}/dp_global" "${STAGE_DIR}/"
else
  echo "[package_bundle] rsync not found; falling back to cp -r (then removing known large dirs)"
  cp -r "${PROJECT_ROOT}/dp_global" "${STAGE_DIR}/"
  rm -rf "${STAGE_DIR}/dp_global/output" || true
  rm -rf "${STAGE_DIR}/dp_global/R/dpglobal_bundle/dist" || true
fi

# Ensure data_simulation data is present (example inputs, may be used by main_cpp.R)
if [ -d "${PROJECT_ROOT}/data_simulation/data" ]; then
  mkdir -p "${STAGE_DIR}/data_simulation"
  cp -v -r "${PROJECT_ROOT}/data_simulation/data" "${STAGE_DIR}/data_simulation/" || true
fi

# Copy helper scripts and top-level runner utilities
mkdir -p "${STAGE_DIR}/bin"
cp -v "${PROJECT_ROOT}/bin/run_dp_future_single.R" "${STAGE_DIR}/bin/" || true
cp -v "${PROJECT_ROOT}/bin/run_dp_future.R" "${STAGE_DIR}/bin/" || true
# Skipping copying Makefile by default (not needed now)
# cp -v "${PROJECT_ROOT}/Makefile" "${STAGE_DIR}/" || true

# Sanity checks: ensure main driver exists in staged copy
if [ ! -f "${STAGE_DIR}/dp_global/scripts/main_cpp.R" ]; then
  echo "[package_bundle] WARNING: main_cpp.R not found in staged dp_global/scripts/ — package may be incomplete"
fi

# Build INSTALL.txt (auto-generated to explain how to run main_cpp.R)
INSTALL_FILE="${STAGE_DIR}/INSTALL.txt"
cat > "${INSTALL_FILE}" <<'EOF'
DP_GLOBAL Bundle - INSTALL

This archive contains a full copy of the project's runtime files required to run `scripts/main_cpp.R`.

Files included (relative to the archive root):
- dp_global/                (full project R code and scripts, excluding heavy outputs)
- dp_global/R/dpglobal_bundle/dpglobal_bundle.RData
- dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds
- dp_global/R/dpglobal_bundle/README.md
- dp_global/R/dpglobal_bundle/transition_cost_rcpp.R
- dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp
- bin/                      (run helpers such as run_dp_future_single.R and run_dp_future.R)
- data_simulation/data/      (example input files)
- Makefile                  [NOT INCLUDED BY DEFAULT]

Quick install & run (from the extracted directory):

1) Install required R packages (non-interactive one-liner):
Rscript -e "manifest <- readRDS('dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds'); pkgs <- unique(unlist(manifest$required_pkgs_by_file)); pkgs <- pkgs[!pkgs %in% installed.packages()[, 'Package']]; if (length(pkgs)) install.packages(pkgs)"

2) Build/compile the C++ acceleration (recommended for speed):
# from the extracted archive root (uses source in dp_global/src packaged in dp_global/R/dpglobal_bundle/src/)
Rscript -e "Rcpp::sourceCpp('dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')"

3) Load the bundle functions (example):
Rscript -e "load('dp_global/R/dpglobal_bundle/dpglobal_bundle.RData')"

4) Verify quickly using a dry-run with the provided runner (recommended):
# Print the Rscript commands that would run without executing
./bin/run_dp_future_single.R --workers 1 --cores-per-job 1 -- --DRY_RUN

5) Run main driver (example; adapt flags as needed):
# Example: run tag 2 in a single-process mode
Rscript dp_global/scripts/main_cpp.R --WHICH_TAG=2 --RUN_ALL_TAGS=TRUE --MANUAL_CORES=TRUE --MANUAL_CORES_VALUE=1

Notes:
- The compiled shared object is platform- and R-version-specific. Compile on the target machine using Rcpp::sourceCpp().
- If a function is missing after loading, rebuild the bundle on the source machine using:
  Rscript -e "source('./dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"

- To let the packaging script build the `dpglobal_bundle.RData` for you before packaging, run:
  sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

EOF

# Show staged files summary and run sanity checks
echo "[package_bundle] Staged contents preview (top-level):"
if command -v tree >/dev/null 2>&1; then
  tree -L 2 "${STAGE_DIR}"
else
  ls -1 "${STAGE_DIR}" | sed -e 's/^/  /'
fi

# Sanity checks for critical files
declare -a expected=(
  "${STAGE_DIR}/dp_global/scripts/main_cpp.R"
  "${STAGE_DIR}/dp_global/R/dpglobal_bundle/dpglobal_bundle.RData"
  "${STAGE_DIR}/dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds"
  "${STAGE_DIR}/dp_global/R/dpglobal_bundle/transition_cost_rcpp.R"
  "${STAGE_DIR}/bin/run_dp_future_single.R"
)
missing=()
for f in "${expected[@]}"; do
  if [ ! -e "$f" ]; then
    missing+=("$f")
  fi
done
if [ ${#missing[@]} -gt 0 ]; then
  echo "[package_bundle] WARNING: The following expected files are missing:"
  for m in "${missing[@]}"; do echo "  - ${m}"; done
  echo "[package_bundle] Please review staging; tarball will still be created but may be incomplete."
fi

# create the tarball
pushd "${STAGE_DIR}" >/dev/null
  echo "[package_bundle] Creating tarball: ${TARBALL_PATH}"
  tar -czf "${TARBALL_PATH}" .
popd >/dev/null

# checksum (sha256)
if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${TARBALL_PATH}" > "${SHA_PATH}"
elif command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${TARBALL_PATH}" > "${SHA_PATH}"
else
  echo "[package_bundle] Warning: no sha256 checksum tool found; skipping checksum generation." >&2
fi

echo "[package_bundle] Created: ${TARBALL_PATH}"
[[ -f "${SHA_PATH}" ]] && echo "[package_bundle] Checksum: ${SHA_PATH}"

# Stage cleaned up via trap

exit 0

# ./dp_global/R/dpglobal_bundle/package_bundle.sh