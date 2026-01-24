#!/usr/bin/env Rscript
# Diagnose anchor recruit plausibility across a recruit_max grid
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
dataset_path <- if (length(args) >= 1) args[1] else NULL
grid_vals <- if (length(args) >= 2) as.numeric(strsplit(args[2], ',')[[1]]) else c(0.5,1,2,5,10)

source(file.path('dp_global','dev','diagnostics','utils.R'))
Dt <- load_dataset(dataset_path)
# require a DP run to get anchor bio params and prior context; try reading an rds
rds_path <- file.path('dp_global','dev','diagnostics','diagnostics_dp_diag.rds')
if (file.exists(rds_path)) {
  res <- readRDS(rds_path); setDT(res)
} else {
  message('[diag_anchor] No DP RDS found; running DP on synthetic data')
  compile_transition_cost()
  res <- run_dp_on_dt(Dt)
  saveRDS(res, rds_path)
}

# find anchors
anchors <- unique(res[!is.na(TrueStemID), .(TrueStemID, CensusID, DBH)])
anchors <- anchors[, .SD[which.min(CensusID)], by=TrueStemID]

out_list <- list()
for (i in seq_len(nrow(anchors))) {
  a <- anchors[i]
  anchor_row <- res[TrueStemID == a$TrueStemID & CensusID == a$CensusID][1]
  best_prev_total <- Inf
  prev_rows <- res[CensusID == (a$CensusID - 1)]
  if (nrow(prev_rows) > 0) {
    for (j in seq_len(nrow(prev_rows))) {
      pr <- prev_rows[j]
      comps <- compute_cost_components(pr$DBH, a$DBH, pr)
      if (comps$total < best_prev_total) best_prev_total <- comps$total
    }
  }
  records <- list()
  for (rm in grid_vals) {
    comps <- compute_cost_components(NA_real_, a$DBH, anchor_row)
    # FIXME: override recruit cap in per-track comps by temporarily changing anchor_row field
    anchor_row$Bio_Recruit_MaxDBH_unit <- rm
    comps <- compute_cost_components(NA_real_, a$DBH, anchor_row)
    p_rec <- p_recruit_from_costs(comps$total, best_prev_total)
    records[[length(records)+1]] <- data.table(TrueStemID = a$TrueStemID, anchor_census = a$CensusID, recruit_max_dbh = rm, recruit_allowed = (comps$per_track$cost_hard[1]==0), recruit_total = comps$total, best_prev_total = best_prev_total, p_recruit = p_rec)
  }
  out_list[[length(out_list)+1]] <- rbindlist(records)
}

out_dt <- rbindlist(out_list)
save_report(out_dt, file.path('dp_global','dev','diagnostics','anchor_plausibility_grid.csv'))

invisible(NULL)
