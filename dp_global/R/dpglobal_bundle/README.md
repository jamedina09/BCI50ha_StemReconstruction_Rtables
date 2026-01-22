# dpglobal_bundle — Portable function bundle for DP_GLOBAL 🔧

This directory bundles the R code and Rcpp helper wrappers needed to run the DP_GLOBAL workflow (the `main_cpp.R` driver). The goal is a portable, reproducible snapshot you can ship to another machine and load quickly.

---

## What is included

- **dpglobal_bundle_loader.R**
  - Script that collects and sources the project R modules into a dedicated environment and saves them as `dpglobal_bundle.RData`.
  - It tries to use the manifest from `dp_global/R/dp_global_main.R` (the `r_files` variable) to preserve the intended load order. If that manifest is not available it falls back to sourcing all `.R` files found in `dp_global/R/`.

- **dpglobal_bundle.RData**
  - The output file containing the R objects (functions, lists, data) that were sourced into the bundle environment.

- **src/transition_cost_rcpp.cpp**
  - The C++ implementation for Rcpp-accelerated functions. **This is source code only** and must be compiled on the target system.

- **transition_cost_rcpp.R** (recommended to add here)
  - Minimal R wrapper that calls `Rcpp::sourceCpp("src/transition_cost_rcpp.cpp")`. Place this in the bundle directory so users can compile the C++ code locally.

---

## How it works (key details) 💡

- The loader (`dpglobal_bundle_loader.R`) creates a new environment (`dpglobal_env`) and attempts to source modules in the order specified by `dp_global/R/dp_global_main.R`'s `r_files` manifest. This preserves internal dependencies.
- If `r_files` is not found or an error occurs, it falls back to sourcing every `.R`/`.r` file in `dp_global/R/` (excluding the bundle directory).
- The loader also sources the R wrapper `dp_global/src/transition_cost_rcpp.R` into `dpglobal_env` (if present) so wrapper functions that rely on compiled code are present in the bundle.
- Finally, everything in `dpglobal_env` is saved to `dpglobal_bundle.RData`.

---

## Quick start (generate & use) 🚀

From the project root (recommended):

1. Create the bundle RData:

```r
source("dp_global/R/dpglobal_bundle/dpglobal_bundle_loader.R")
# -> creates dp_global/R/dpglobal_bundle/dpglobal_bundle.RData
```

2. On the target system/session:

```r
# load functions (recommended into a new env to avoid global conflicts)
myenv <- new.env()
load("/path/to/dp_global/R/dpglobal_bundle/dpglobal_bundle.RData", envir = myenv)

# compile C++ code (must be done each session)
# If a wrapper exists in the bundle, you can run it; otherwise:
Rcpp::sourceCpp("/path/to/dp_global/src/transition_cost_rcpp.cpp")

# call functions (from global env if you loaded there) or via myenv:
myenv$estimate_bio_pars(...)
# or if loaded to global env: estimate_bio_pars(...)
```

---

## Verifying the bundle ✅

After creating or loading the bundle, verify core functions are present:

```r
core_fns <- c("estimate_bio_pars", "match_stems_dp_global_backward_marginals_batch", "add_dp_posterior_bins", "transition_cost_tracks_bio_batch_rcpp")
missing <- core_fns[!vapply(core_fns, function(s) exists(s, mode = "function", inherits = TRUE), logical(1L))]
if (length(missing)) {
  message("Missing functions: ", paste(missing, collapse = ", "))
} else {
  message("All core functions present.")
}
```

If you saved into a custom environment (`myenv`) use `exists("estimate_bio_pars", envir = myenv)` or check `ls(envir = myenv)`.

---

## Limitations & tips ⚠️

- **RData cannot contain compiled binaries.** The RData stores R objects only (functions, closures, lists, data). The compiled shared library must be built on the target system.
- Always run the Rcpp compilation step (via `transition_cost_rcpp.R` or `Rcpp::sourceCpp`) before calling functions that rely on C++.
- Prefer loading into a separate environment (`myenv`) to avoid masking or polluting the global environment.
- If functions are missing after loading, check the loader messages (it prints warnings/errors for any source failures) and rebuild the bundle.

---
