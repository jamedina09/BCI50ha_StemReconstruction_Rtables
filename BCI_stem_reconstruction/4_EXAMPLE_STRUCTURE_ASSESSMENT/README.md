# BCI Stem Reconstruction — Example Structure Assessment

This folder contains three analysis scripts for the BCI 50-ha stem reconstruction project:
`biomass_stocks_fluxes.R` and `basal_area_uncertainty.R` write their results to `outputs/`
under this folder; `general_plot_information.R` prints its summaries to the console.

Inputs:

- the stage-3 R tables `DATA/RTABLES/bci.stem1.Rdata` … `bci.stem9.Rdata` (all scripts);
- the species table `DATA/RTABLES/bci.spptable.rdata` (all three scripts). Stage 3 does not write it: copy the ForestGEO-format
  table there (for example `data_paper_and_repo_publication/RTABLES/bci.spptable.rdata`);
- the posterior samples `DATA/POSTERIORS/posterior_sampled_paths.rds` and the stage-2 table
  `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds`
  (`basal_area_uncertainty.R`);
- the wood density file in `wd/` (`biomass_stocks_fluxes.R`).

In the R tables:

- `stemID` is the stem number within its tree (1, 2, ...), so a stem is identified by
  `treeID` + `stemID` (every stem key in these scripts uses both);
- `ExactDate` is the raw field date (`NA` where there is no record); each script completes
  it from `date` (days since 1960-01-01, filled in stage 3 for every row), so every
  stem-census row has a date;
- `dbh` is exported as recorded, so a dead (`G`/`D`) record can carry a DBH; each script
  ignores it, because stocks, fluxes and summaries use alive (`A`) rows only.

## Scripts

### `biomass_stocks_fluxes.R`

Purpose:

- Estimate aboveground biomass (AGB) stocks, productivity, mortality, and net AGB change
  across nine BCI stem censuses (1982–2022/3).

Key processing steps:

- Load the stage-3 R tables for censuses 1–9.
- Merge taxonomy and wood specific gravity (WSG) data from BCI species tables and the
  2026 Wright & Muller-Landau Dryad WSG dataset.
- Fill missing WSG hierarchically by genus, then family, then global mean.
- Apply optional strangler-fig removal and palm DBH correction.
- Convert DBH from mm to cm and apply Cushman et al. 2014 taper correction to estimate
  DBH at 1.3 m (`dbh_t`).
- Interpolate missing DBH for alive stems using `linear` (alternatives: `locf`, or `mean`).
- Compute AGB with Chave et al. 2014 allometry plus Martinez-Cano et al. 2019 height model,
  with Goodman et al. 2013 palm-specific allometry.
- Correct 1985 small-stem DBH rounding bias using a Census 3 reference.
- Compute growth, recruitment, and mortality fluxes per stem and aggregate them by quadrat
  and size class.
- Optionally apply Kohyama et al. 2019 bias correction to productivity and mortality at the
  quadrat level.

Outputs:

- `outputs/plot_agb_dynamics.png` — standing AGB, productivity/mortality, and net AGB change
  for the whole plot.
- `outputs/plot_agb_by_size.png` — the same summaries stratified by DBH size class.

### `basal_area_uncertainty.R`

Purpose:

- Propagate stem-identity uncertainty from posterior reconstruction paths produced by the
  `dp_global` engine into basal area (BA) stocks and fluxes for the BCI 50-ha plot.

Key processing steps:

- Load the exported reconstruction (R tables) and the posterior paths
  (`DATA/POSTERIORS/posterior_sampled_paths.rds`).
- Treat palms and strangler figs as `biomass_stocks_fluxes.R` does (flags
  `remove_strangler_figs` and `use_median_palm_dbh`, both `TRUE`): giant strangler figs
  (*Ficus costaricana*, *obtusifolia*, *popenoei*, *trigonata* with a DBH > 500 mm in any
  census) are removed with every stem of their tree and leave the posterior paths too; every
  measured alive stem of a palm species other than *Socratea* gets the species' median DBH over
  all censuses. A palm's basal area is then constant, so identity uncertainty changes only its
  Loss and Gain.
- Separate trees with a single path (no identity uncertainty) from trees with several paths.
- Complete every path with the stage-2 links of observations outside its window (splice).
- Run Monte Carlo realizations that draw one path per multi-path tree, independently, with
  weight `path_count / sum(path_count)`.
- Aggregate basal area stocks and fluxes at tree and quadrat scales.
- Report the exported reconstruction and the empirical MC uncertainty (95 % interval).

Posterior weights and the two engines:

- Every tree has 200 posterior draws, collapsed into unique paths. `path_count` is the number
  of draws that produced a path, so `path_count / 200` is its posterior probability and each
  draw is equally likely. `path_prob` is **not** used: for DP trees it re-weights each draw by
  its own probability, which counts that probability twice.
- DP trees are sampled by backward sampling from the exact DP posterior, so likely
  trajectories repeat (median 4 unique paths per multi-path tree, some drawn many times).
- Trees routed to the probabilistic engine (`dp_probabilistic_matching.R`: palms and other
  forced species, stranglers, trees the DP cannot solve because the state space is too large
  or no assignment is feasible) get approximate draws: a noisy assignment per census pair
  that respects the database pins (two stems with different pins are never joined),
  stitched and repaired for growth violations. Almost every draw differs somewhere, so nearly
  every path is unique and weighs 1/200 (0.005). These trees have no single most probable
  path; their exported reconstruction is the most representative draw.
- Both engines are sampled together and in the same way (one draw per tree per realization).
  The probabilistic draws are an approximation rather than a calibrated posterior; the
  diagnostics report both engines separately (trees, unique paths, share of paths drawn once,
  and how often the exported partition is among the sampled paths).

Key assumptions and scope:

- Palm diameters and strangler figs are corrected as in `biomass_stocks_fluxes.R`; the other
  methodological choices (taper correction aside, no buttress or Kohyama correction) are not
  applied here. On the current tables this removes 29 trees and lowers the plot BA stock by
  0.5–0.9 m² ha⁻¹ (giant figs and the DBH fluctuations of palms).
- Uncertainty applies only to pre-anchor intervals because post-anchor censuses have confirmed
  stem identity. The anchor census used in this script is `ANCHOR_START_CENSUS = 7`.
- Only identity uncertainty is quantified. Trees are drawn independently, so their variations
  largely cancel at plot level and the ribbon is narrow; it does not include sampling
  (quadrat) uncertainty or model/systematic error.
- The script writes several outputs to `outputs/`, including MAP feather tables, MC realization
  feather files, summary files, and diagnostic figures.

Outputs (written to `outputs/`):

MAP (deterministic) tables:

- `ba_map_change_treeID.feather` — MAP tree-level BA flux per census pair.
- `ba_map_change_quadrat.feather` — MAP quadrat-level BA flux per census pair.
- `ba_map_stock_quadrat.feather` — MAP quadrat-level BA stock per census.

Monte Carlo realizations and summaries:

- `ba_mc_realizations_quadrat.feather` — quadrat-level MC realization fluxes.
- `ba_mc_realizations_stock_quadrat.feather` — quadrat-level MC realization stocks.
- `ba_mc_summary_quadrat.feather` — empirical 95 % CI of quadrat fluxes.
- `ba_mc_summary_stock_quadrat.feather` — empirical 95 % CI of quadrat stocks.
- `ba_mc_realizations_treeID/` — one feather per MC realization with tree-level fluxes, written only
  when `write_tree_realizations <- TRUE` (default `FALSE`; not tracked by git).

Figures:

- `fig1_BA_stock.pdf` — forest-level BA stock: exported reconstruction vs MC.
- `fig2_BA_fluxes.pdf` — forest-level BA fluxes: exported reconstruction vs MC.

Diagnostics:

- `ba_mc_diagnostics.txt` — tree counts (single path, sampled, fallback, by engine), splice
  counts, the exported partition vs the posterior per engine, invariance checks and the MC
  error of the 95 % interval bounds.
- `fig3_BA_trajectories.pdf` — individual-tree BA trajectories for a selected subset.

### `general_plot_information.R`

Purpose:

- Descriptive summaries of forest structure and composition for the 50-ha plot across all
  censuses (census 1, 1982, is excluded from the temporal comparisons).

Key processing steps:

- Load the stage-3 R tables and the species table; count an individual as one `treeID` and
  an alive stem as `Rstatus == "A"`.
- Use the Cushman et al. 2014 taper-corrected DBH (`use_taper_for_ba`), give alive stems whose
  DBH was missed a DBH from their own measurements (the rule of the other two scripts), and
  place trees without coordinates at the centre of their quadrat.
- Summarise the most recent census (individuals, stems, basal area, species, genera, families,
  lifeforms, diversity), per-hectare estimates by bootstrap over quadrats, temporal trends and
  changes between censuses, diversity across censuses, and the DBH size-class distribution
  (census 2 vs census 9).

Outputs: printed to the console; no files are written.
