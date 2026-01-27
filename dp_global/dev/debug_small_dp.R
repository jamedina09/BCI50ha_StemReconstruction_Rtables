#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("dp_global/R/dp_global_main.R")

input_path <- "/Users/medinaja/Library/CloudStorage/OneDrive-SmithsonianInstitution/STRI/STEM_TABLES_FORESTGEO/3_BCI/2_STEM_IDENTIFICATION/example.csv"
d <- fread(input_path)
# construct small subset
 d6 <- d[CensusID==6 & !is.na(DBH)]
 d7 <- d[CensusID==7 & !is.na(DBH)][1:5]
 d_small <- rbindlist(list(d6, d7))[order(CensusID)]
cat("Running small subset DP with options(error=traceback)\n")
options(error = function(){
  cat("Error encountered. Traceback:\n")
  traceback()
  q(status = 1)
})
# This call should either complete or print a traceback
match_stems_dp_global_backward_marginals_batch(d_small, anchor_start=7, slack_tracks=0, verbose=TRUE)
cat("If you see this message, the small DP completed without error.\n")
