# Run provisional anchor tests by sourcing package R files
r_files <- list.files('dp_global/R', full.names = TRUE)
# Only source regular .R files (skip subdirectories)
r_files <- r_files[endsWith(r_files, '.R')]
for (f in r_files) {
  tryCatch(source(f), error = function(e) stop(sprintf('Error sourcing %s: %s', f, e$message)))
}

library(testthat)
# Run only our test file
test_file('dp_global/tests/test-provisional-anchor.R', reporter = 'summary')
