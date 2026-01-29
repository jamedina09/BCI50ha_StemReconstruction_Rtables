#!/usr/bin/env Rscript
# verify_bundle.R — basic smoke tests to validate an extracted bundle
# Runs a small set of regression tests expected to exist in the bundle.

cat("Running bundle verification tests...\n")

tests <- c(
  "dp_global/dev/test_anchor_scoping.R",
  "dp_global/dev/test_integration_post_anchor_given.R"
)

all_ok <- TRUE
for (t in tests) {
  if (!file.exists(t)) {
    cat(sprintf("MISSING: %s\n", t))
    all_ok <- FALSE
    next
  }
  cat(sprintf("Running: %s ... ", t))
  rc <- system2("Rscript", args = t, stdout = "", stderr = "")
  if (rc == 0) {
    cat("PASS\n")
  } else {
    cat("FAIL (exit code:", rc, ")\n")
    all_ok <- FALSE
  }
}

if (all_ok) {
  cat("Bundle verification: OK\n")
  quit(status = 0)
} else {
  cat("Bundle verification: FAILED\n")
  quit(status = 1)
}
