# Test that run_dp_future includes RECRUIT flags in DRY_RUN mode
# This test runs a dry-run of two configs and checks the printed commands.
# Test fixed config
cmd1 <- c("--workers", "1", "--cores-per-job", "1", "--configs", "fixed", "--", "--DRY_RUN")
out1 <- system2("Rscript", args = c("bin/run_dp_future.R", cmd1), stdout = TRUE, stderr = TRUE)
if (!any(grepl("\\[DRY_RUN\\]", out1))) stop("No DRY_RUN lines found for fixed config; run may have failed or DRY_RUN was not passed.")
out1_str <- paste(out1, collapse = " ")
if (!grepl("\\[DRY_RUN\\].*fixed:.*--RECRUIT_MAX_SOURCE=fixed", out1_str)) stop("Fixed config did not include --RECRUIT_MAX_SOURCE=fixed in DRY_RUN output")
if (!grepl("\\[DRY_RUN\\].*fixed:.*--RECRUIT_MAX_FIXED=6", out1_str)) stop("Fixed config did not include --RECRUIT_MAX_FIXED=6 in DRY_RUN output")

# Test data_hard config
cmd2 <- c("--workers", "1", "--cores-per-job", "1", "--configs", "data_hard", "--", "--DRY_RUN")
out2 <- system2("Rscript", args = c("bin/run_dp_future.R", cmd2), stdout = TRUE, stderr = TRUE)
if (!any(grepl("\\[DRY_RUN\\]", out2))) stop("No DRY_RUN lines found for data_hard config; run may have failed or DRY_RUN was not passed.")
out2_str <- paste(out2, collapse = " ")
if (!grepl("\\[DRY_RUN\\].*data_hard:.*--RECRUIT_MAX_SOURCE=data", out2_str)) stop("data_hard config did not include --RECRUIT_MAX_SOURCE=data in DRY_RUN output")
cat("test_run_dp_future_recruit_max: OK\n")
