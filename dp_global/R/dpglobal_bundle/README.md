# dpglobal_bundle — Focused bundle folder (what's here)

This folder contains a compact set of runtime artifacts and helpers that make it easy
to run the `dp_global` codebase on another machine without checking out the full
repository. Typical contents (when present):

- `dpglobal_bundle_loader.R` — helper to build `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` on the source machine.
- `dpglobal_bundle.RData` — optional pre-sourced R objects (convenience; not required).
- `dpglobal_bundle_manifest.rds` — optional manifest listing suggested R packages.
- `src/` and `transition_cost_rcpp.R` — C++ source and R wrapper for the transition-cost acceleration (compile on target with Rcpp::sourceCpp).
- `package_bundle.sh` — packaging helper that creates a timestamped tarball in `dist/`.

Quick RUN (recommended and minimal) ✅

1) Extract the bundle somewhere convenient (or place the extracted folder into your project directory):

```bash
mkdir -p /path/to/extracted_bundle
tar -xzf dpglobal_bundle_full_YYYYmmdd_HHMMSS.tar.gz -C /path/to/extracted_bundle
```

2) Run the main entrypoint from R with a temporary working directory — this avoids editing paths and keeps your session tidy.

# Interactive R example (recommended)

```r
withr::with_dir('/path/to/extracted_bundle', source(file.path('dp_global', 'R', 'dp_global_main.R')))
```

Notes:
- The packaged `INSTALL.txt` contains this same minimal example (the packaging script writes this quick-run snippet).
- `withr::with_dir()` temporarily sets the working directory for the sourced code and automatically restores your original working directory after the call.
- `dpglobal_bundle.RData` may be loaded to provide quick access to pre-sourced helper objects, but it is not required to run the code.

Optional: compile C++ acceleration for speed (recommended if you need performance)

```bash
Rscript -e "Rcpp::sourceCpp('dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')"
```

Manifest & packages

- If present, the manifest at `dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds` lists recommended packages. Install missing packages with:

```r
manifest <- readRDS('dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds')
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!pkgs %in% installed.packages()[, 'Package']]
if (length(pkgs)) install.packages(pkgs)
```

Quick verification (optional)

```r
# load into a fresh environment to avoid polluting your session
myenv <- new.env()
load(file.path('/path/to/extracted_bundle', 'dp_global', 'R', 'dpglobal_bundle.RData'), envir = myenv)
exists('estimate_bio_pars', envir = myenv, mode = 'function')  # should be TRUE if RData was created
```

Rebuilding the bundle (source machine)

- Run `dpglobal_bundle_loader.R` to regenerate `dpglobal_bundle.RData` and the manifest, then re-run `package_bundle.sh`.

```bash
# from project root on the machine that has the full source
Rscript -e "source('./dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle
```

Why this is simpler

- No path-editing or tricky sourcing: extract the bundle, point `withr::with_dir()` to the extracted bundle root, and `source('dp_global/R/dp_global_main.R')`.
- Loading `dpglobal_bundle.RData` is optional and only a convenience to avoid sourcing many files manually.
- The C++ compile step is independent and only required for performance.

Security & reproducibility

- Record R and package versions used to build the bundle (e.g., `sessioninfo::session_info()`).
- Do not run compiled code from untrusted bundles without verifying provenance.

That's it — short and practical. If you'd like, I can also add a tiny `verify_bundle.R` script that demonstrates the one-liner `withr::with_dir(...)` pattern and performs a basic smoke test.

---

package_bundle.sh — detailed behavior and usage

Purpose
- Create a timestamped tarball that contains everything needed to run `dp_global/scripts/main_cpp.R` on another machine.

What it stages (full package, by default)
- A copy of the `dp_global/` tree (it uses `rsync` if available) with these exclusions: `dp_global/output/`, `dp_global/R/dpglobal_bundle/dist/`, and `.git/` (keeps the archive free of large outputs and VCS metadata).
- `data_simulation/data/` (example inputs) and `bin/` runner scripts (e.g., `run_dp_future_single.R`, `run_dp_future.R`) and `Makefile` from project root.
- The bundle-local files: `dp_global/R/dpglobal_bundle/*` (if present) including `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` when available.
- The R wrapper and C++ source for the transition-cost functions are copied from `dp_global/src/` into the staged `dp_global/R/dpglobal_bundle/src/` so the README paths are consistent and `Rcpp::sourceCpp()` can be run as described.
- An `INSTALL.txt` auto-generated inside the archive with exact copy-paste commands to install packages, compile the C++, and run `main_cpp.R`.

Key behavior and safeguards
- Uses `rsync -av` when available; otherwise falls back to `cp -r` and removes known large dirs.
- Prints a staged contents preview and runs simple sanity checks (warns if expected critical files are missing). The tarball will still be created but you are warned if something appears missing.
- Generates a SHA256 checksum (`.sha256`) for the created tarball (if `shasum` or `sha256sum` is available).
- Creates a temporary staging directory; the staging area is cleaned up automatically on exit (no changes to your repo files).

Options
- `--build-bundle` — before packaging, runs:
    Rscript -e "source('dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"
  This generates/updates `dpglobal_bundle.RData` and `dpglobal_bundle_manifest.rds` in `dp_global/R/dpglobal_bundle` and ensures the tarball contains the RData.
- `--help` or `-h` — prints usage and exits.

Where the tarball is written
- `dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_<timestamp>.tar.gz`
- Checksum: `dp_global/R/dpglobal_bundle/dist/dpglobal_bundle_full_<timestamp>.tar.gz.sha256`

How to produce a full package (copy-paste)
- From the project root (recommended):

```bash
# Full package (uses existing dpglobal_bundle.RData if present)
sh dp_global/R/dpglobal_bundle/package_bundle.sh

# Full package and build dpglobal_bundle.RData first (recommended to include RData)
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle
```

How to use the package on the target machine (copy-paste)
- Extract the tarball:

```bash
mkdir -p ~/projects/dp_global_bundle
tar -xzf dpglobal_bundle_full_<timestamp>.tar.gz -C ~/projects/dp_global_bundle
cd ~/projects/dp_global_bundle
```

- Follow `INSTALL.txt` inside the extracted archive or run these copy-paste commands directly (from the extracted root):

```bash
# Install R packages listed in the manifest (non-interactive):
Rscript -e "manifest <- readRDS('dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds'); pkgs <- unique(unlist(manifest$required_pkgs_by_file)); pkgs <- pkgs[!pkgs %in% installed.packages()[, 'Package']]; if (length(pkgs)) install.packages(pkgs)"

# Compile the C++ acceleration (recommended):
Rscript -e "Rcpp::sourceCpp('dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp')"

# Load the prebuilt functions (optional sample):
Rscript -e "load('dp_global/R/dpglobal_bundle/dpglobal_bundle.RData')"

# Dry-run example (print Rscript invocations without running DP):
./bin/run_dp_future_single.R --workers 1 --cores-per-job 1 -- --DRY_RUN
```

Notes and recommendations
- The archive contains compiled-source (**not** compiled binaries). Always run `Rcpp::sourceCpp()` on the target machine to build the shared object locally (platform- and R-version-dependent).

---


Minimal quickstart (copy-paste-ready commands)

From the project root (recommended):

```r
# 1) Load the R functions into the global session (or into an environment)
load("./dp_global/R/dpglobal_bundle/dpglobal_bundle.RData")

# 2) Source the R wrapper (R-side wrapper comes from dp_global/src and is staged into the bundle):
sys.source("./dp_global/R/dpglobal_bundle/transition_cost_rcpp.R", envir = globalenv())

# 3) Compile the C++ implementation (optional but recommended for speed):
Rcpp::sourceCpp("./dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp")
```

Or, if you've cd'ed into the bundle directory (shorter paths):

```r
load("dpglobal_bundle.RData")
sys.source("transition_cost_rcpp.R", envir = globalenv())
Rcpp::sourceCpp("src/transition_cost_rcpp.cpp")
```

Install packages listed in the manifest (one-liner)

```r
manifest <- readRDS("dp_global/R/dpglobal_bundle/dpglobal_bundle_manifest.rds")
pkgs <- unique(unlist(manifest$required_pkgs_by_file))
pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(pkgs)) install.packages(pkgs)
```

Quick verification (recommended)

```r
# load into a new environment to avoid polluting global env
myenv <- new.env()
load("./dp_global/R/dpglobal_bundle/dpglobal_bundle.RData", envir = myenv)
exists("estimate_bio_pars", envir = myenv, mode = "function")  # should be TRUE
# check compiled symbol (once you've run Rcpp::sourceCpp):
exists("transition_cost_tracks_bio_batch_rcpp_cpp", mode = "function", inherits = TRUE)
```

Rebuild the bundle (on the source machine)

```bash
# run from project root on the machine that has the full source code
# Option A: build the RData manually, then package
Rscript -e "source('./dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R')"
sh dp_global/R/dpglobal_bundle/package_bundle.sh --build-bundle

# Option B: have the packaging script build it for you
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

- Automated packaging: you can create a single tarball containing the bundle and runtime files with the provided packaging script. From the project root run:

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

7. Run a short smoke test or the provided tests:

```bash
# example test that exists in the project
Rscript dp_global/dev/test_transition_cost_tracks_bio_batch.R
```


### Troubleshooting & FAQ ❓

Q: Some functions are missing after I load the RData. What now?

- A: Check `dpglobal_bundle_loader.R` output when the bundle was created; any warnings about files that failed to source will be printed. Verify that you have the required packages installed and the bundle was created with the same set of modules. If needed, rebuild the bundle on the source machine.

Q: I get an error about a missing compiled symbol.

- A: This means the C++ code hasn't been compiled in the current session. Compile it with Rcpp::sourceCpp("./dp_global/R/dpglobal_bundle/src/transition_cost_rcpp.cpp") (requires `Rcpp` and a C++ toolchain).

Q: Is it safe to copy the compiled shared object from one machine to another?

- A: Not recommended. Shared objects are platform and R-version specific. Recompile on each target to avoid ABI/compatibility issues.


### Security & reproducibility notes 🔐

- For reproducibility, keep a record of the R version and package versions used to create the bundle (consider `sessioninfo::session_info()`).
- Avoid running untrusted compiled code; verify provenance of the bundle before compiling or loading it.


---