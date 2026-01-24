#!/usr/bin/env Rscript
# Example runner that creates fake data and runs the main diagnostics
suppressPackageStartupMessages({ library(data.table) })

source(file.path('dp_global','examples','diagnostics','generate_synthetic_data.R'))
source(file.path('dp_global','examples','diagnostics','utils.R'))

# generate fake data
cat('[run_synth] generating synthetic dataset (diagnostics_dataset.csv)\n')
system2('Rscript', args = c('dp_global/examples/diagnostics/generate_synthetic_data.R', file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv')))

# run a quick DP and save
Dt <- load_dataset(file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv'))
compile_transition_cost()
res <- run_dp_on_dt(Dt, max_tracks = 50, slack_tracks = 0)
saveRDS(res, file.path('dp_global','examples','diagnostics','diagnostics_dp_fake.rds'))
message('[run_synth] saved synthetic DP RDS')

# run diagnostics
cat('[run_synth] running anchor plausibility grid\n')
system2('Rscript', args = c('dp_global/examples/diagnostics/diagnose_anchor_plausibility.R', file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv'), '0.5,1,2,5,10'))
cat('[run_synth] running generic sweep (cost-only)\n')
system2('Rscript', args = c('dp_global/examples/diagnostics/sweep_param_generic.R', 'Bio_Recruit_MaxDBH_unit', '0.5,1,2,5,10', 'cost', file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv')))
cat('[run_synth] running track cost timeseries for probe=6\n')
system2('Rscript', args = c('dp_global/examples/diagnostics/track_cost_timeseries.R', '6', 'Bio_Recruit_MaxDBH_unit', '0.5,1,2,5,10', file.path('dp_global','examples','diagnostics','diagnostics_dataset.csv')))

cat('[run_synth] done. Reports saved to dp_global/examples/diagnostics/\n')

invisible(NULL)
