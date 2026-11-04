# run_tests.R — runs the testthat files given as arguments, one after the
# other, from the project root, and stops with a non-zero exit status if a
# test fails or errors (or if a file is not found), so that `make test` fails
# when the code does.
#   Rscript run_tests.R dp_global/tests/test_enumerate_states.R ...
files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0L) stop("run_tests.R: give the test files to run")
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) {
    cat("❌ Test file not found:", paste(missing_files, collapse = ", "), "\n")
    quit(status = 1L)
}
suppressMessages(library(testthat))
totals <- c(tests = 0L, failed = 0L, skipped = 0L, errors = 0L)
bad <- character()
for (.i in seq_along(files)) {
    # the path is kept in .path: a test file may assign other variables
    .path <- files[[.i]]
    cat("\n== ", .path, "\n", sep = "")
    .res <- as.data.frame(testthat::test_file(.path, reporter = "summary", stop_on_failure = FALSE))
    .n <- c(tests = nrow(.res), failed = sum(.res$failed > 0), skipped = sum(.res$skipped), errors = sum(.res$error))
    totals <- totals + .n
    if (.n[["failed"]] > 0L || .n[["errors"]] > 0L) bad <- c(bad, files[[.i]])
    cat(sprintf("%s: tests %d | failed %d | skipped %d | errors %d\n", files[[.i]], .n[["tests"]], .n[["failed"]], .n[["skipped"]], .n[["errors"]]))
}
cat(sprintf(
    "\nTotal: tests %d | failed %d | skipped %d | errors %d\n",
    totals[["tests"]], totals[["failed"]], totals[["skipped"]], totals[["errors"]]
))
if (length(bad) > 0L) {
    cat("❌ Failing test files:", paste(bad, collapse = ", "), "\n")
    quit(status = 1L)
}
cat("✓ All tests passed\n")
