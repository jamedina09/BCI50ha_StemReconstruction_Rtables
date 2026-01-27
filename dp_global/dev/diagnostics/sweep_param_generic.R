#!/usr/bin/env Rscript
# Generic parameter sweep: either cost-only or full DP per grid value
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
param_name <- if (length(args)>=1) args[1] else 'Bio_Recruit_MaxDBH_unit'
grid_str <- if (length(args)>=2) args[2] else '0.5,1,2,5,10'
mode <- if (length(args)>=3) args[3] else 'cost' # 'cost' or 'full'
dataset_path <- if (length(args)>=4) args[4] else NULL

source(file.path('dp_global','dev','diagnostics','utils.R'))
Dt <- load_dataset(dataset_path)
vals <- as.numeric(strsplit(grid_str, ',')[[1]])
out <- data.table(param = vals, chosen_chain_total = NA_real_, alternative_chain_total = NA_real_, note = NA_character_)
for (i in seq_along(vals)) {
  v <- vals[i]
  Dt2 <- copy(Dt)
  Dt2[, (param_name) := v]
  if (mode == 'full') {
    compile_transition_cost()
    res <- run_dp_on_dt(Dt2)
    # compute a representative chosen chain total for probe (first TrueStemID anchor if present)
    probe_anchor <- res[!is.na(TrueStemID)][order(CensusID)][1]
    if (!is.null(probe_anchor) && nrow(probe_anchor)>0) {
      # use previous census row to compute chain total vs alt
      r_prev <- res[CensusID == (probe_anchor$CensusID - 1)][1]
      if (!is.null(r_prev) && !is.na(r_prev$DBH)) {
        chosen_total <- compute_cost_components(r_prev$DBH, probe_anchor$DBH, r_prev)$total
        alt_total <- compute_cost_components(NA_real_, probe_anchor$DBH, probe_anchor)$total + compute_cost_components(r_prev$DBH, NA_real_, r_prev)$total
        out[i, `:=`(chosen_chain_total = chosen_total, alternative_chain_total = alt_total)]
      } else out[i, note := 'no prev candidate']
    } else out[i, note := 'no anchors']
  } else {
    # cost-only: compute NA->anchor cost for first anchor row and best_prev_total
    compile_transition_cost()
    res <- run_dp_on_dt(Dt2)
    anchor_row <- res[!is.na(TrueStemID)][order(CensusID)][1]
    if (!is.null(anchor_row) && !is.na(anchor_row$DBH)) {
      comps <- compute_cost_components(NA_real_, anchor_row$DBH, anchor_row)
      best_prev <- Inf
      prev_rows <- res[CensusID == (anchor_row$CensusID - 1)]
      if (nrow(prev_rows)>0) for (j in seq_len(nrow(prev_rows))) best_prev <- min(best_prev, compute_cost_components(prev_rows$DBH[j], anchor_row$DBH, prev_rows[j])$total)
      out[i, `:=`(chosen_chain_total = NA_real_, alternative_chain_total = comps$total, note = NA_character_)]
    } else out[i, note := 'no anchors']
  }
}
save_report(out, file.path('dp_global','dev','diagnostics',paste0('sweep_generic_',gsub('[^A-Za-z0-9]','_',param_name),'.csv')))

invisible(NULL)
