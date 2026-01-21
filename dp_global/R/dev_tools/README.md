# dev_tools

This folder contains developer helper scripts to create, verify, and use a packaged set of `dp_global` functions and the associated Rcpp source.

Files
- `save_dp_functions.R` — Creates (or refreshes) `dp_functions.rdata` by sourcing `dp_global` modules and (optionally) compiling the C++ acceleration. Run: `Rscript dp_global/R/dev_tools/save_dp_functions.R`.
- `test_load_dp_functions.R` — Verifies the saved `dp_functions.rdata` can be loaded, the functions attached into a fresh environment, and the embedded C++ source can be recompiled. Run: `Rscript dp_global/R/dev_tools/test_load_dp_functions.R`.
- `dp_functions.rdata` — The exported object created by `save_dp_functions.R`. Contains at minimum:
  - `dp_functions` — named list of function objects exported from `dp_global` R modules.
  - `cpp_source` — character vector of the `transition_cost_rcpp.cpp` source (if present).
  - `cpp_filename` — basename of the C++ file.
  - `compiled_ok` — logical indicating whether compilation succeeded during save.

Quick usage
1. Load into R and (optionally) compile the embedded C++:

```r
# Load the saved objects and compile the embedded C++ into the session
source("dp_global/R/dev_tools/load_dp_functions.R")  # provides helpers
res <- load_dp_functions_rdata("dp_global/R/dev_tools/dp_functions.rdata", compile_cpp = TRUE)
# Attach functions into your current session
list2env(res$functions, envir = .GlobalEnv)
```

2. Or use the convenience helper to attach and recompile in one step:

```r
source("dp_global/R/load_dp_functions.R")
attach_dp_functions_from_rdata("dp_global/R/dev_tools/dp_functions.rdata", compile_cpp = TRUE)
```

Notes
- The `.rdata` stores the raw C++ source (text) rather than a compiled binary; this avoids cross-platform binary issues. The loader writes a temp `.cpp` file and calls `Rcpp::sourceCpp()` to compile on demand.
- You need `Rcpp` installed in the session where you recompile. If `Rcpp` is not present compilation is skipped with a warning.
- Use `save_dp_functions.R` to refresh `dp_functions.rdata` after changing R or C++ source.
