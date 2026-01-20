# Test run_dp_full.sh produces RECRUIT flags in DRY_RUN output
out <- system2("bash", args = c("bin/run_dp_full.sh", "--config=fixed", "--DRY_RUN"), stdout = TRUE, stderr = TRUE)
out_str <- paste(out, collapse = " ")
if (!grepl("--RECRUIT_MAX_SOURCE=fixed", out_str)) stop("run_dp_full.sh: missing --RECRUIT_MAX_SOURCE=fixed in DRY_RUN output")
if (!grepl("--RECRUIT_MAX_FIXED=6", out_str)) stop("run_dp_full.sh: missing --RECRUIT_MAX_FIXED=6 in DRY_RUN output")
cat("test_run_dp_full_recruit_max: OK\n")
