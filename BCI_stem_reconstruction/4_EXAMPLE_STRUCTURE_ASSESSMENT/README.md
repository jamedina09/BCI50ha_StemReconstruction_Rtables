# 4_EXAMPLE_STRUCTURE_ASSESSMENT

Stage 4 of the pipeline. Three example analyses of the R tables: biomass
stocks and fluxes, basal-area stocks and fluxes with their stem-identity
uncertainty, and plot summaries.

## Contents

| File | What it does |
| --- | --- |
| `biomass_stocks_fluxes.R` | Aboveground biomass (AGB) stocks, productivity, mortality and net change, for the plot and by size class |
| `basal_area_uncertainty.R` | Basal-area (BA) stocks and fluxes of the exported reconstruction, with the uncertainty from the posterior samples |
| `general_plot_information.R` | Structure and composition summaries, printed to the console |
| `wd/` | Wood density dataset (input; not under version control) |
| `outputs/` | Outputs of the first two scripts |

## How to run

From the project root, after stage 3, in any order:

```sh
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/biomass_stocks_fluxes.R
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/general_plot_information.R
```

## Inputs

| File | Used by |
| --- | --- |
| `DATA/RTABLES/bci.stem1.Rdata` … `bci.stem9.Rdata` | all three scripts |
| `DATA/RTABLES/bci.spptable.rdata` | all three scripts |
| `DATA/POSTERIORS/posterior_sampled_paths.rds` | `basal_area_uncertainty.R` |
| `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` | `basal_area_uncertainty.R` |
| `wd/doi_10_5061_dryad_5qfttdzn3__v20260403/WD_species.txt` | `biomass_stocks_fluxes.R`; download from <https://datadryad.org/dataset/doi:10.5061/dryad.5qfttdzn3> |

## Outputs

| File | Content |
| --- | --- |
| `outputs/plot_agb_dynamics.png` | Standing AGB, productivity, mortality and net AGB change of the plot |
| `outputs/plot_agb_by_size.png` | The same by DBH size class |
| `outputs/ba_map_change_treeID.feather`, `ba_map_change_quadrat.feather` | BA fluxes per tree and per quadrat and census pair, exported reconstruction |
| `outputs/ba_map_stock_quadrat.feather` | BA stock per quadrat and census, exported reconstruction |
| `outputs/ba_mc_realizations_quadrat.feather`, `ba_mc_realizations_stock_quadrat.feather` | BA fluxes and stocks per quadrat in each Monte Carlo realization |
| `outputs/ba_mc_summary_quadrat.feather`, `ba_mc_summary_stock_quadrat.feather` | 95% interval of the quadrat fluxes and stocks |
| `outputs/ba_mc_realizations_treeID/` | One file per realization with tree-level fluxes, written only with `write_tree_realizations <- TRUE` (default `FALSE`; not under version control) |
| `outputs/fig1_BA_stock.pdf`, `fig2_BA_fluxes.pdf` | Plot BA stock and fluxes: exported reconstruction and Monte Carlo interval |
| `outputs/fig3_BA_trajectories.pdf` | BA trajectories of a sample of trees |
| `outputs/ba_mc_diagnostics.txt` | Tree counts by engine, splice counts, the exported reconstruction against the posterior, invariance checks and the Monte Carlo error of the interval bounds |

## What it does

The three scripts read the R tables in the same way:

- A stem is `treeID` + `stemID`.
- A missing `ExactDate` is completed from `date`, so every row has a date.
- Only alive (`A`) rows enter stocks, fluxes and summaries; the DBH that a
  `G` or `D` record can carry is ignored.
- A stump (an `A` record without a DBH whose `DFstatus` is "broken below") is
  out of the stock in that census: no DBH, no stem count, no basal area or
  biomass. Before the stem's first measurement the stem enters the stock at
  that measurement (a recruit). Between two measurements its basal area or
  biomass leaves the stock and comes back (a loss, then a gain); it is not a
  death.
- Any other `A` record without a DBH is a missed measurement and gets a DBH
  from the stem's own measurements.

### `biomass_stocks_fluxes.R`

1. Loads the R tables of censuses 1–9, the species table and the wood density
   (Wright and Muller-Landau); missing wood density is filled by genus, then
   family, then the plot mean.
2. Removes the giant strangler figs and gives every palm stem its species'
   median DBH (flags at the top of the script).
3. Converts DBH to cm and corrects it to 1.3 m with the taper model of Cushman
   et al. (2014).
4. Interpolates the DBH of missed measurements (`linear`; `locf` and `mean`
   are alternatives); stump censuses get none.
5. Computes AGB with Chave et al. (2014) and the height model of
   Martínez Cano et al. (2019), and with Goodman et al. (2013) for palms.
6. Corrects the 1985 small-stem DBH rounding with census 3 as reference.
7. Computes growth, recruitment and mortality per stem, aggregates them by
   quadrat and size class, and applies the Kohyama et al. (2019) correction to
   productivity and mortality per quadrat (optional).

### `basal_area_uncertainty.R`

1. Loads the R tables, the posterior samples and the stage-2 table.
2. Treats palms and strangler figs as `biomass_stocks_fluxes.R` does
   (`remove_strangler_figs`, `use_median_palm_dbh`, both `TRUE`): the giant
   strangler figs (*Ficus costaricana*, *obtusifolia*, *popenoei*,
   *trigonata* above 500 mm in any census) are removed with every stem of
   their tree, and every measured alive palm stem other than *Socratea* gets
   its species' median DBH.
3. Separates trees with one sampled path from trees with several, and
   completes every path with the stage-2 links of observations outside its
   window (splice).
4. Runs Monte Carlo realizations: each draws one path per tree with several
   paths, independently, with weight `path_count / sum(path_count)`.
5. Aggregates BA stocks and fluxes per tree and quadrat, and reports the
   exported reconstruction with the 95% Monte Carlo interval.

Posterior weights:

- Every tree has 200 draws, collapsed into unique paths; `path_count` is the
  number of draws of a path, so each draw weighs 1/200. `path_prob` is not
  used: for DP trees it re-weights each draw by its own probability, which
  counts that probability twice.
- DP trees are sampled from the exact DP posterior, so likely paths repeat.
  Trees of the probabilistic matcher get approximate draws, nearly all
  different (weight 1/200 each); their exported reconstruction is the most
  representative draw. Both engines are sampled in the same way.
- Only identity uncertainty is quantified, before the anchor census
  (`ANCHOR_START_CENSUS = 7`). Trees are drawn independently, so the plot
  interval is narrow; it holds no sampling (quadrat) or model uncertainty.
  Apart from the palm and fig treatment and the taper correction, the
  corrections of `biomass_stocks_fluxes.R` are not applied.

### `general_plot_information.R`

1. Loads the R tables and the species table. An individual is a `treeID`; an
   alive stem is `Rstatus == "A"`, except a stump census.
2. Uses the taper-corrected DBH (`use_taper_for_ba`), gives missed
   measurements a DBH from the stem's own measurements and places trees
   without coordinates at the centre of their quadrat.
3. Summarises the last census (individuals, stems, basal area, species,
   genera, families, life forms, diversity), per-hectare estimates by
   bootstrap over quadrats, trends and changes between censuses (census 1 is
   left out), diversity across censuses, and the DBH size classes of censuses
   2 and 9.

## Requirements

R packages: `data.table`, `truncnorm`, `lubridate`, `HDInterval`, `ggplot2`,
`scales`, `cowplot` (`biomass_stocks_fluxes.R`); `data.table`, `ggplot2`,
`patchwork`, `scales`, `collapse`, `arrow` (`basal_area_uncertainty.R`);
`data.table`, `boot`, `HDInterval` (`general_plot_information.R`).
