# dpglobal_bundle

Portable packaging of the `dp_global` stem-identification algorithm for deployment on other machines.

## Contents

| File | Purpose |
|------|---------|
| `package_bundle.sh` | Creates a timestamped `.tar.gz` in `dist/` containing everything needed to run the algorithm |
| `dpglobal_bundle_loader.R` | Generates `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` (run on source machine before packaging) |
| `verify_bundle.R` | Optional smoke test for bundle integrity on the target machine |
| `src/` | C++ source for transition-cost acceleration (compiled on target with `Rcpp::sourceCpp`) |

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
- Record your R and package versions for reproducibility (`sessioninfo::session_info()`).
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle
```

Package the project (full package) and create checksum

```bash
# Create a full package (includes a copy of dp_global/ excluding outputs)
sh dp_global/R/dpglobal_bundle/package_bundle.sh

# Create a full package and build `dpglobal_bundle.RData` first to include the RData
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

# The tarball and checksum are written to:
# dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz
# dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz.sha256
```
Notes & troubleshooting
- `dpglobal_bundle.RData` does NOT contain compiled binaries. Compile the C++ on the target with `Rcpp::sourceCpp("./dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp")` (or the shorter path if you're in the bundle dir).
- If a function is missing after loading, rebuild the bundle using `dpglobal_bundle_loader.R` on the source machine and re-copy both the `RData` and `manifest`.
- If compiled symbols are missing after compilation, run the `sourceCpp()` line interactively and re-check the `exists()` call above.


## Full, detailed workflow (step-by-step) 📋

Below is a comprehensive walkthrough — follow these steps when preparing a bundle and installing it on a remote or separate machine.

### A) Source machine — create the bundle

1. Prerequisites (source machine):
   - R (>= your project's R version). Recommended: same or similar R minor version on target systems for reproducibility.
   - Development toolchain (macOS): `xcode-select --install` (required to compile Rcpp locally on this machine).
   - R packages: `Rcpp`, `here`, and any packages used by the project (the manifest will list per-file requirements, e.g., `data.table`, `igraph`, `MASS`).

2. Generate the bundle artifacts (from project root):

```bash
Rscript -e "source('./dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"
# or
Rscript dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R
```

This produces:
- `dp_global/R/dpglobal_bundle/dpglobal_bundle.RData`
- `dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds`


### B) Target machine — prepare and verify the bundle

1. Transfer the tarball (or the bundle folder) to the target machine and unpack it:

```bash
# example: unpack into your working directory
mkdir -p ~/projects/dp_global_bundle
tar -xzf dpglobal_bundle_YYYYMMDD_HHMMSS.tar.gz -C ~/projects/dp_global_bundle
```

2. Inspect the manifest to see package requirements and whether compiled acceleration was present when the bundle was created:

```r
manifest <- readRDS('dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds')
print(manifest)
```

3. Install missing packages (recommended):

- Interactive/manual: open R and run the following snippet to install missing packages listed in the manifest:

```r
manifest <- readRDS('dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds')
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!is.na(pkgs) & pkgs != ""]
missing <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(missing)) install.packages(missing)
```

- Automated packaging: create a tarball with the provided packaging script. From the project root:

```bash
# Create package (does not auto-build dpglobal_bundle.RData unless it exists)
sh dp_global/R/dpglobal_bundle/package_bundle.sh

# Create package and build `dpglobal_bundle.RData` first to include the RData
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

# result: dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_YYYYMMDD_HHMMSS.tar.gz
# and a checksum file next to it
```

4. Compile the C++ acceleration to obtain runtime speedups:

```bash
# From project root (recommended):
Rscript -e "Rcpp::sourceCpp('dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')"
# Or within an interactive R session:
# Rcpp::sourceCpp('dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')
```

Notes:
- Compilation requires a C++ toolchain (e.g., Xcode on macOS). If compilation fails, the wrapper fallback still exists and code will run in pure R (slower).
- The compiled shared object is platform- and R-version-specific. Always compile on the target machine.

5. Load the functions (recommended into a new env):

```r
myenv <- new.env()
load('dp_global/R/dpglobal_bundle/dpglobal_bundle.RData', envir = myenv)
ls(envir = myenv)
# quick check
exists('estimate_bio_pars', envir = myenv)
```

6. Verify compiled C++ symbol (if you compiled):

```r
exists('transition_cost_tracks_bio_batch_rcpp_cpp', mode = 'function', inherits = TRUE)
```

7. Run a short smoke test or the provided verify script:

```bash
Rscript dp_global/R/dpglobal_bundle/verify_bundle.R
```


### Troubleshooting & FAQ ❓

Q: Some functions are missing after I load the RData. What now?

- A: Check `dpglobal_bundle_loader.R` output when the bundle was created; any warnings about files that failed to source will be printed. Verify that you have the required packages installed and the bundle was created with the same set of modules. If needed, rebuild the bundle on the source machine.

Q: I get an error about a missing compiled symbol.

- A: This means the C++ code hasn't been compiled in the current session. Compile it with `Rcpp::sourceCpp('./dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')` (requires `Rcpp` and a C++ toolchain).

Q: Is it safe to copy the compiled shared object from one machine to another?

- A: Not recommended. Shared objects are platform and R-version specific. Recompile on each target to avoid ABI/compatibility issues.


### Security & reproducibility notes 🔐

- For reproducibility, keep a record of the R version and package versions used to create the bundle (consider `sessioninfo::session_info()`).
- Avoid running untrusted compiled code; verify provenance of the bundle before compiling or loading it.


---