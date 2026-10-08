# dpglobal_bundle

Portable packaging of the `dp_global` stem-identification algorithm for deployment on other machines.

## Contents

| File | Purpose |
|------|---------|
| `package_bundle.sh` | Creates a timestamped `.tar.gz` in `dist/` containing everything needed to run the algorithm (bash script) |
| `dpglobal_bundle_loader.R` | Generates `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` (run on source machine before packaging) |
| `verify_bundle.R` | Smoke test for bundle integrity — verifies all R modules load and core functions exist |
| `dist/` | Generated tarballs and checksums. The tarball in the repository (`dpglobal_bundle_full_20260529_211653.tar.gz`) holds the code as it was on that date; rebuild before deploying |

## 1. Build the bundle (source machine)

From the project root:

```bash
bash dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle
```

`--build-bundle` regenerates `dpglobal_bundle.RData` and the manifest before packaging (it runs `dpglobal_bundle_loader.R`, which sources `dp_global_main.R` and therefore needs the R packages listed under *Full Workflow* and a C++ toolchain). Output: `dist/dpglobal_bundle_full_<timestamp>.tar.gz` with a SHA256 checksum.

The script uses `rsync` (falls back to `cp -r`) to copy the whole `dp_global/` tree and excludes `dp_global/output/`, `dist/`, and `.git/`. It also stages `data_simulation/data/` as example input and writes an `INSTALL.txt` inside the archive. `README.md`, `verify_bundle.R`, `dpglobal_bundle.RData` and the manifest are also copied to the archive root.

## 2. Deploy on the target machine

```bash
mkdir -p ~/dp_global_bundle
tar -xzf dpglobal_bundle_full_<timestamp>.tar.gz -C ~/dp_global_bundle
cd ~/dp_global_bundle
```

### Install R packages

```r
manifest <- readRDS("dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds")
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs)) install.packages(pkgs)
```

`lpSolve` is used by the probabilistic matcher but is not in the manifest; install it as well (`install.packages("lpSolve")`).

### Compile the C++ functions

```r
Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
```

The archive contains source code, not binaries — always compile on the target machine. Sourcing `dp_global_main.R` (next step) runs this compilation itself; the DP engine cannot run without the compiled functions, so a C++ toolchain is required.

### Run

```r
withr::with_dir("~/dp_global_bundle",
  source(file.path("dp_global", "R", "dp_global_main.R")))
```

`withr::with_dir()` temporarily sets the working directory and restores it afterward (`dp_global_main.R` resolves its files from the working directory). The bundled `INSTALL.txt` contains the same `source()` command; for the manual compilation it points to the copy of the C++ file at `dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp`.

`dp_global_main.R` defines the engine and the post-engine helpers. To run a full reconstruction use a driver script, e.g. `Rscript dp_global/scripts/main_cpp.R` (see `dp_global/scripts/README.md`).

## 3. Verify (optional)

```r
myenv <- new.env()
load("dp_global/R/dpglobal_bundle/dpglobal_bundle.RData", envir = myenv)
exists("estimate_bio_pars", envir = myenv, mode = "function")
exists("match_stems_probabilistic", envir = myenv, mode = "function")
```

## Notes

- `dpglobal_bundle.RData` is a convenience snapshot of pre-sourced R objects. It is **not required** — the `source()` path above loads everything directly from the R files. The snapshot does not hold the compiled C++ functions.
- The bundle includes the probabilistic matching module (`dp_probabilistic_matching.R`) which handles tags with very large state spaces.
- The bundle also includes `basal_area_uncertainty.R` for posterior-based uncertainty quantification.
- Record your R and package versions for reproducibility (`sessioninfo::session_info()`).

## Full Workflow (step-by-step)

### A) Source machine — create the bundle

1. Prerequisites:
   - R (≥ 4.0). Recommended: same or similar R minor version on target systems.
   - Development toolchain (macOS): `xcode-select --install`.
   - R packages: `Rcpp`, `here`, `data.table`, `MASS`, `igraph`, `lpSolve`. The driver scripts also use `arrow`, `ggplot2`, `cowplot`, `reshape2`, `scales` and `patchwork`; the commands in this file use `withr`.

2. Build and package (from project root):

```bash
# Build RData + manifest, then create tarball
bash dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

# Or just package (if RData already exists)
bash dp_global/R/dpglobal_bundle/package_bundle.sh

# Output:
# dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz
# dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz.sha256
```

### B) Target machine — deploy and verify

1. Transfer and unpack:

```bash
mkdir -p ~/projects/dp_global_bundle
tar -xzf dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz -C ~/projects/dp_global_bundle
cd ~/projects/dp_global_bundle
```

1. Install missing R packages:

```r
manifest <- readRDS("dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds")
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!is.na(pkgs) & pkgs != ""]
missing <- pkgs[!pkgs %in% installed.packages()[, 1]]
if (length(missing)) install.packages(missing)
```

1. Compile the C++ functions (also done when `dp_global_main.R` is sourced):

```r
Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
```

The archive contains source code, not binaries — always compile on the target machine.

1. Run the verification script:

```bash
Rscript dp_global/R/dpglobal_bundle/verify_bundle.R
```

1. Run your analysis:

```r
withr::with_dir("~/projects/dp_global_bundle",
  source(file.path("dp_global", "R", "dp_global_main.R")))
```

### Troubleshooting

- **Missing functions after loading RData**: Rebuild the bundle on the source machine. Check `dpglobal_bundle_loader.R` output for sourcing errors.
- **Missing compiled symbol**: Compile with `Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")`. Requires `Rcpp` package and a C++ toolchain.
- **Cross-platform shared objects**: Do not copy compiled `.so`/`.dll` files between machines. Always recompile on the target.
- **Reproducibility**: Record `sessioninfo::session_info()` for your R and package versions.

## Stem Identity Renumbering

All drivers in the dp_global workflow use a universal post-engine helper chain: `maybe_add_posterior_bins()`, `apply_pin_track_rejoin()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_terminal_to_host()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, and finally `finalize_posterior_paths()`. After these steps, all `ReconstructedStemID` values are renumbered **sequentially from 1 to N within each tag**, ordered by the earliest census in which each stem appears. If multiple stems first appear in the same census, the largest DBH at that census gets the lower ID, with ties broken by original ID. **Negative or zero IDs are never produced.**
