#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("dp_global/R/dp_global_main.R")

input_path <- "/Users/medinaja/Library/CloudStorage/OneDrive-SmithsonianInstitution/STRI/STEM_TABLES_FORESTGEO/3_BCI/2_STEM_IDENTIFICATION/example.csv"
cat("Reading:", input_path, "\n")
d <- fread(input_path)
cat("Rows read:", nrow(d), "\n")

# Inject user-provided biological parameters (force to dataset so DP sees them)
bio_vals <- list(
  Bio_Mu_Growth = 0.01452756,
  Bio_Gamma_Growth = 0.07163091,
  Bio_Sigma0_Growth = 0.003578623,
  Bio_Sigma1_Growth = 0.006697133,
  Bio_H0_Mortality = 0.04685938,
  Bio_Beta_Mortality = -0.02755028,
  Bio_Recruit_Meanlog = 0.4339083,
  Bio_Recruit_Sdlog = 0.4847472,
  Bio_Recruit_MaxDBH_unit = 36.5001,
  Bio_Recruitment_lambda = 0.08044933,
  Bio_Max_Shrink = -1,
  Bio_K_Shrink = 0,
  Bio_Max_Growth = 7.5,
  Bio_Max_Growth_Soft = 0.8861551,
  Bio_K_Growth = 0
)
for (nm in names(bio_vals)) {
  d[, (nm) := bio_vals[[nm]]]
}
cat("Injected Bio parameters into dataset\n")

run_dp_safe <- function(dt, anchor, slack=0, label="run"){
  cat(sprintf("\n--- %s: anchor=%d slack=%d nrows=%d unique_census=%s ---\n",
              label, anchor, slack, nrow(dt), paste(sort(unique(dt$CensusID)), collapse=",")))
  res <- tryCatch({
    match_stems_dp_global_backward_marginals_batch(dt, anchor_start=anchor, slack_tracks=slack, verbose=TRUE)
  }, error=function(e){
    cat(sprintf("ERROR (%s): %s\n", label, e$message))
    return(NULL)
  })
  if (is.null(res)) return(NULL)
  cat(sprintf("Result class: %s\n", paste(class(res), collapse=",")))
  cat(sprintf("Unique DP_KUsed: %s\n", paste(unique(res$DP_KUsed), collapse=",")))
  cat(sprintf("Unique DP_MaxStatesPerCensus: %s\n", paste(unique(res$DP_MaxStatesPerCensus), collapse=",")))
  pi <- attr(res, "DP_PruneInfo")
  if (!is.null(pi)) {
    cat("DP_PruneInfo:\n")
    print(pi)
  }
  cat("Head (first 6 rows):\n")
  print(head(res, 6))
  cat("Tail (last 6 rows):\n")
  print(tail(res, 6))
  invisible(res)
}

# 1) Full dataset, anchor 7
res_full <- run_dp_safe(d, 7, slack=0, label="Full dataset")

# 2) Small sampled subset: keep census 6 fully, sample 5 rows from census 7
cat("Sampling small subset: census 6 (all), census 7 (first 5)\n")
d6 <- d[CensusID==6 & !is.na(DBH)]
d7 <- d[CensusID==7 & !is.na(DBH)][1:5]
d_small <- rbindlist(list(d6, d7))[order(CensusID)]
res_small <- run_dp_safe(d_small, 7, slack=0, label="Small sampled subset")

cat("\nRerunning small subset with debug handler to capture stack on error...\n")
res_small_debug <- tryCatch({
  match_stems_dp_global_backward_marginals_batch(d_small, anchor_start=7, slack_tracks=0, verbose=TRUE)
}, error=function(e){
  cat("DEBUG ERROR: ", e$message, "\n", sep="")
  cat("sys.calls():\n")
  print(sys.calls())
  cat("sessionInfo():\n")
  print(sessionInfo())
  NULL
})

cat('\nDone.\n')
