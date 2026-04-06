#!/usr/bin/env Rscript
# verify_bundle.R — smoke tests to validate an extracted dp_global bundle
# Run from the extracted bundle root: Rscript dp_global/R/dpglobal_bundle/verify_bundle.R

cat("Running bundle verification...\n")

all_ok <- TRUE
fail <- function(msg) { cat("FAIL:", msg, "\n"); all_ok <<- FALSE }
pass <- function(msg) { cat("PASS:", msg, "\n") }

# 1. Check that core R modules exist
core_files <- c(
  "dp_global/R/dp_global_main.R",
  "dp_global/R/dp_global_bio.R",
  "dp_global/R/dp_global_dp.R",
  "dp_global/R/dp_global_states.R",
  "dp_global/R/dp_global_matchers.R",
  "dp_global/R/dp_probabilistic_matching.R",
  "dp_global/R/dp_global_utils.R",
  "dp_global/R/dp_global_diag.R"
)
for (f in core_files) {
  if (file.exists(f)) {
    pass(sprintf("File exists: %s", f))
  } else {
    fail(sprintf("Missing file: %s", f))
  }
}

# 2. Check driver scripts
drivers <- c(
  "dp_global/scripts/main_cpp.R",
  "dp_global/scripts/main_cpp_chunk.R",
  "dp_global/scripts/basal_area_uncertainty.R"
)
for (f in drivers) {
  if (file.exists(f)) {
    pass(sprintf("Driver exists: %s", f))
  } else {
    fail(sprintf("Missing driver: %s", f))
  }
}

# 3. Check C++ source
cpp_file <- "dp_global/src/transition_cost_rcpp.cpp"
if (file.exists(cpp_file)) {
  pass("C++ source exists")
} else {
  fail("Missing C++ source: transition_cost_rcpp.cpp")
}

# 4. Source all modules (verifies they parse without error)
cat("\nSourcing dp_global_main.R...\n")
tryCatch({
  source("dp_global/R/dp_global_main.R")
  pass("All R modules sourced successfully")
}, error = function(e) {
  fail(sprintf("Sourcing failed: %s", e$message))
})

# 5. Check key functions are available
required_fns <- c(
  "estimate_bio_pars",
  "match_stems_dp_global_backward_marginals_batch",
  "match_stems_optimal_backward",
  "match_stems_probabilistic",
  "enumerate_states_injective",
  "add_dp_posterior_bins"
)
for (fn in required_fns) {
  if (exists(fn, mode = "function", inherits = TRUE)) {
    pass(sprintf("Function available: %s", fn))
  } else {
    fail(sprintf("Function missing: %s", fn))
  }
}

# 6. Check bundle artifacts (optional)
if (file.exists("dp_global/R/dpglobal_bundle/dpglobal_bundle.RData")) {
  pass("dpglobal_bundle.RData exists")
} else {
  cat("NOTE: dpglobal_bundle.RData not present (optional — modules can be sourced directly)\n")
}

cat("\n")
if (all_ok) {
  cat("Bundle verification: OK\n")
  quit(status = 0)
} else {
  cat("Bundle verification: FAILED\n")
  quit(status = 1)
}
