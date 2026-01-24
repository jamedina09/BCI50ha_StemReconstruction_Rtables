#!/usr/bin/env Rscript
# Compute per-parameter per-component cost time series for a specified reconstructed track
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
track_probe <- if (length(args)>=1) as.integer(args[1]) else 6
param_name <- if (length(args)>=2) args[2] else 'Bio_Recruit_MaxDBH_unit'
grid_vals <- if (length(args)>=3) as.numeric(strsplit(args[3],',')[[1]]) else c(0.5,1,2,5,10)
book1_path <- if (length(args)>=4) args[4] else NULL

source(file.path('dp_global','examples','diagnostics','utils.R'))
Dt <- load_dataset(book1_path)
compile_transition_cost()

out_rows <- list()
for (v in grid_vals) {
  Dt2 <- copy(Dt)
  Dt2[, (param_name) := v]
  res <- run_dp_on_dt(Dt2)
  setDT(res)
  # find the probe's rows across censuses
  probe_rows <- res[ReconstructedStemID == track_probe]
  if (nrow(probe_rows) == 0) next
  # Sort by census and compute consecutive pair costs
  probe_rows <- probe_rows[order(CensusID)]
  for (k in seq_len(nrow(probe_rows)-1)) {
    a <- probe_rows[k]
    b <- probe_rows[k+1]
    comps <- compute_cost_components(a$DBH, b$DBH, a)
    out_rows[[length(out_rows)+1]] <- data.table(param = v, from = a$CensusID, to = b$CensusID, total = comps$total, hard = comps$per_track$cost_hard[1], soft = comps$per_track$soft_penalty[1], growth = comps$per_track$growth_penalty[1])
  }
}
out_dt <- rbindlist(out_rows)
save_report(out_dt, file.path('dp_global','examples','diagnostics',paste0('track_cost_timeseries_probe_',track_probe,'.csv')))

invisible(NULL)
