# BCI Stem Reconstruction — Biomass Stocks and Fluxes

This folder contains two analysis scripts for the BCI 50-ha stem reconstruction project.
All results are written to `outputs/` under this folder.

## Scripts

### `biomass_stocks_fluxes.R`

Purpose:

- Estimate aboveground biomass (AGB) stocks, productivity, mortality, and net AGB change
  across nine BCI stem censuses (1982–2022/3).

Key processing steps:

- Load raw BCI stem RTABLE files for censuses 1–9.
- Merge taxonomy and wood specific gravity (WSG) data from BCI species tables and the
  2026 Wright & Muller-Landau Dryad WSG dataset.
- Fill missing WSG hierarchically by genus, then family, then global mean.
- Apply optional strangler-fig removal and palm DBH correction.
- Convert DBH from mm to cm and apply Cushman et al. 2014 taper correction to estimate
  DBH at 1.3 m (`dbh_t`).
- Interpolate missing DBH for alive stems using `linear`, `locf`, or `mean`.
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

- Load posterior stem reconstruction data.
- Separate deterministic MAP-equivalent paths from ambiguous trees with multiple reconstruction
  paths.
- Run Monte Carlo realizations that sample one path per ambiguous tree proportional to path
  probability.
- Aggregate basal area stocks and fluxes at tree and quadrat scales.
- Report MAP estimates and empirical MC uncertainty (95% CI).

Key assumptions and scope:

- Uncertainty applies only to pre-anchor intervals because post-anchor censuses have confirmed
  stem identity. The anchor census used in this script is `ANCHOR_START_CENSUS = 7`.
- The script writes several outputs to `outputs/`, including MAP feather tables, MC realization
  feather files, summary files, and diagnostic figures.

Outputs:

- `outputs/fig1_BA_stock.pdf`
- `outputs/fig2_BA_fluxes.pdf`
- `outputs/fig3_BA_trajectories.pdf`
