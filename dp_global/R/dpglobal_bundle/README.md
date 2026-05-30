# dpglobal_bundle

Portable packaging of the `dp_global` stem-identification algorithm for deployment on other machines.

## Contents

| File | Purpose |
|------|---------|
| `package_bundle.sh` | Creates a timestamped `.tar.gz` in `dist/` containing everything needed to run the algorithm |
| `dpglobal_bundle_loader.R` | Generates `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` (run on source machine before packaging) |
| `verify_bundle.R` | Smoke test for bundle integrity — verifies all R modules load and core functions exist |
| `dist/` | Generated tarballs and checksums (not tracked by git) |

## 1. Build the bundle (source machine)

From the project root:

```bash
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle
```

`--build-bundle` regenerates `dpglobal_bundle.RData` and the manifest before packaging. Output: `dist/dpglobal_bundle_full_<timestamp>.tar.gz` with a SHA256 checksum.

The script uses `rsync` (falls back to `cp -r`) and excludes `dp_global/output/`, `dist/`, and `.git/`. It also stages `data_simulation/data/` as example input and writes an `INSTALL.txt` inside the archive.

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

### Compile C++ acceleration (recommended)

```r
Rcpp::sourceCpp("dp_global/src/transition_cost_rcpp.cpp")
```

The archive contains source code, not binaries — always compile on the target machine.

### Run

## Stem Identity Renumbering (2024+)

All drivers in the dp_global workflow use a universal post-engine helper chain: `maybe_add_posterior_bins()`, `apply_carried_terminal_backfill()`, `apply_orphan_stem_backfill()`, `apply_broken_below_invariants()`, `renumber_engine_minted_ids()`, and finally `finalize_posterior_paths()`. After these steps, all `ReconstructedStemID` values are renumbered **sequentially from 1 to N within each tag**, ordered by the earliest census in which each stem appears. If multiple stems first appear in the same census, the largest DBH at that census gets the lower ID, with ties broken by original ID. **Negative or zero IDs are never produced.**

```r
withr::with_dir("~/dp_global_bundle",
  source(file.path("dp_global", "R", "dp_global_main.R")))
```

`withr::with_dir()` temporarily sets the working directory and restores it afterward. The bundled `INSTALL.txt` contains the same commands.

## 3. Verify (optional)

```r
myenv <- new.env()
load("dp_global/R/dpglobal_bundle/dpglobal_bundle.RData", envir = myenv)
exists("estimate_bio_pars", envir = myenv, mode = "function")
exists("match_stems_probabilistic", envir = myenv, mode = "function")
```

## Notes

- `dpglobal_bundle.RData` is a convenience snapshot of pre-sourced R objects. It is **not required** — the `source()` path above loads everything directly from the R files.
- The bundle includes the probabilistic matching module (`dp_probabilistic_matching.R`) which handles tags with very large state spaces.
- The bundle also includes `basal_area_uncertainty.R` for posterior-based uncertainty quantification.
- Record your R and package versions for reproducibility (`sessioninfo::session_info()`).

## Full Workflow (step-by-step)

### A) Source machine — create the bundle

1. Prerequisites:
   - R (≥ 4.0). Recommended: same or similar R minor version on target systems.
   - Development toolchain (macOS): `xcode-select --install`.
   - R packages: `Rcpp`, `here`, `data.table`, `MASS`, `igraph`.

2. Build and package (from project root):

```bash
# Build RData + manifest, then create tarball
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

# Or just package (if RData already exists)
sh dp_global/R/dpglobal_bundle/package_bundle.sh

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

1. Compile C++ acceleration (recommended):

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
