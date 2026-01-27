#!/usr/bin/env Rscript
# Parameterized synthetic dataset CSV generator for diagnostics
suppressPackageStartupMessages({ library(data.table); library(lubridate) })

args <- commandArgs(trailingOnly = TRUE)
# Usage: Rscript generate_synthetic_data.R [out_path] [n_tags] [censuses_comma] [seed] [recruit_tags_comma] [bio_recruit_maxdbhs_comma] [noise_sd]
out_path <- if (length(args) >= 1) args[1] else file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv')
n_tags <- if (length(args) >= 2) as.integer(args[2]) else 12
censuses <- if (length(args) >= 3) as.integer(strsplit(args[3],',')[[1]]) else c(5L,6L,7L)
seed <- if (length(args) >= 4) as.integer(args[4]) else 42
recruit_tags <- if (length(args) >= 5) as.integer(strsplit(args[5],',')[[1]]) else c(6L,2L)
# single value or comma list -> take first for default recruit maxDBH but we'll set per-row after
bio_recruit_maxdbhs <- if (length(args) >= 6) as.numeric(strsplit(args[6],',')[[1]]) else 5
noise_sd <- if (length(args) >= 7) as.numeric(args[7]) else 0.5

set.seed(seed)
Tags <- seq_len(n_tags)
rows <- list()
for (t in censuses) {
  for (tag in Tags) {
    base <- 5 + (tag %% 4) * 2
    dbh <- base + (t - min(censuses)) * runif(1, 0, 3) + rnorm(1, 0, noise_sd)
    rows[[length(rows)+1]] <- data.table(Tag = tag, CensusID = t, DBH = round(max(dbh, 0.5), 2), ExactDate = as.character(ymd('2000-01-01') + years(t*5)))
  }
}
Dt <- rbindlist(rows)
# Add anchors by default: set TrueStemID for selected tags at the last census
Dt[, TrueStemID := NA_integer_]
for (rt in recruit_tags) {
  Dt[Tag == rt & CensusID == max(censuses), TrueStemID := rt]
}

# Add Bio_ columns with sensible defaults; ensure numeric types to avoid integer truncation warnings
Dt[, Bio_Recruit_MaxDBH_unit := as.numeric(bio_recruit_maxdbhs[1])]
Dt[, Bio_Max_Shrink := as.numeric(-2.5)]
Dt[, Bio_K_Shrink := as.numeric(0)]
Dt[, Bio_K_Growth := as.numeric(0)]
Dt[, Bio_Mu_Growth := as.numeric(0.5)]
Dt[, Bio_Gamma_Growth := as.numeric(0.1)]
Dt[, Bio_Sigma0_Growth := as.numeric(0.2)]
Dt[, Bio_Sigma1_Growth := as.numeric(0.05)]
Dt[, Bio_Max_Growth := as.numeric(7.5)]
Dt[, Bio_Max_Growth_Soft := as.numeric(6)]
Dt[, Bio_H0_Mortality := as.numeric(0.01)]
Dt[, Bio_Beta_Mortality := as.numeric(0.001)]
Dt[, Bio_Recruit_Meanlog := as.numeric(log(2))]
Dt[, Bio_Recruit_Sdlog := as.numeric(0.5)]
Dt[, Bio_Recruitment_lambda := as.numeric(0.1)]

# Write file
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(Dt, out_path)
message('[fake] wrote synthetic dataset to ', out_path)

invisible(NULL)
