# 4_EXAMPLE_STRUCTURE_ASSESSMENT

Stage 4 of the pipeline. Three example analyses of the R tables: biomass
stocks and fluxes, basal-area stocks and fluxes with their stem-identity
uncertainty, and plot summaries.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [What it does](#what-it-does)
  - [`biomass_stocks_fluxes.R`](#biomass_stocks_fluxesr)
  - [`basal_area_uncertainty.R`](#basal_area_uncertaintyr)
  - [`general_plot_information.R`](#general_plot_informationr)
- [Settings](#settings)
- [Checks](#checks)
- [Requirements](#requirements)

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
| `DATA/RTABLES/bci.spptable.rdata` (from `3_PREPARE_R_TABLES/3_create_species_table_BCI.R`) | all three scripts |
| `DATA/POSTERIORS/posterior_sampled_paths.rds` | `basal_area_uncertainty.R` |
| `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` | `basal_area_uncertainty.R` (stage-2 links for the splice) |
| `wd/doi_10_5061_dryad_5qfttdzn3__v20260403/WD_species.txt` | `biomass_stocks_fluxes.R`; download from <https://datadryad.org/dataset/doi:10.5061/dryad.5qfttdzn3> (Wright and Muller-Landau) |

## Outputs

| File | Content |
| --- | --- |
| `outputs/plot_agb_dynamics.png` | Standing AGB, productivity, mortality and net AGB change of the plot |
| `outputs/plot_agb_by_size.png` | The same by DBH size class |
| `outputs/ba_map_change_treeID.feather`, `ba_map_change_quadrat.feather` | BA fluxes per tree and per quadrat and census pair, exported reconstruction |
| `outputs/ba_map_stock_quadrat.feather` | BA stock per quadrat and census, exported reconstruction |
| `outputs/ba_mc_realizations_quadrat.feather`, `ba_mc_realizations_stock_quadrat.feather` | BA fluxes and stocks per quadrat in each Monte Carlo realization (`write_quadrat_realizations`) |
| `outputs/ba_mc_summary_quadrat.feather`, `ba_mc_summary_stock_quadrat.feather` | 95% interval of the quadrat fluxes and stocks |
| `outputs/ba_mc_realizations_treeID/` | One file per realization with tree-level fluxes, written only with `write_tree_realizations <- TRUE` (not under version control) |
| `outputs/fig1_BA_stock.pdf`, `fig2_BA_fluxes.pdf` | Plot BA stock and fluxes: exported reconstruction and Monte Carlo interval |
| `outputs/fig3_BA_trajectories.pdf` | BA trajectories of a sample of trees: exported reconstruction and posterior paths |
| `outputs/ba_mc_diagnostics.txt` | Tree counts (one path, sampled, fallback; by engine), splice counts, the exported reconstruction against the posterior per engine, invariance checks and the Monte Carlo error of the interval bounds |

`general_plot_information.R` writes no files.

## What it does

The three scripts read the R tables in the same way:

- A stem is `treeID` + `stemID` (`stemID` is the stem number within its tree).
- A missing `ExactDate` is completed from `date` (days since 1960-01-01,
  filled by stage 3 on every row), so every row has a date.
- Only alive (`A`) rows enter stocks, fluxes and summaries; the DBH that a
  `G` or `D` record can carry (stage 3 exports DBH as recorded) is ignored.
- A **stump** (an `A` record without a DBH whose `DFstatus` is "broken
  below") is a stem alive with no measurable DBH at 1.3 m (a broken-below
  record with a field tag is that stem). It is out of the stock in that census:
  no DBH, no stem count, no basal area or biomass. Before the stem's first
  measurement the stem enters the stock at that measurement (a recruit).
  Between two measurements its basal area or biomass leaves the stock and comes
  back (a loss, then a gain); its `Rstatus` stays `A`, so it is not a death.
- Any other `A` record without a DBH is a missed measurement and gets a DBH
  from the stem's own measurements.

### `biomass_stocks_fluxes.R`

For each census interval c → c+1 it computes growth (AGB change of stems
alive at both ends), recruitment (AGB of a stem at its first alive census,
divided by the interval) and mortality (AGB at the last alive census, divided
by the interval), per stem, then aggregates them.

1. Loads the R tables of censuses 1–9, the species table and the wood density;
   a species without wood density takes the mean of its genus, else its
   family, else all species (each species counted once).
2. Removes the giant strangler figs (trees of *Ficus costaricana*,
   *obtusifolia*, *popenoei* or *trigonata* with a stem above 500 mm) and gives
   every palm stem, *Socratea* aside, its species' median DBH
   (`remove_strangler_figs`, `use_median_palm_dbh`).
3. Converts DBH to cm and corrects it to 1.3 m with the taper model of Cushman
   et al. (2014) (the recorded value is kept in `dbh_raw`).
4. Interpolates the DBH of missed measurements (`dbh_interp_method`: `linear`
   between flanking measurements, `locf` or `mean`); stump censuses get none.
5. Computes AGB with Chave et al. (2014) (`biomass_allometry`; `chave05`
   available) and the height model of Martínez Cano et al. (2019)
   (`use_local_height_allometry`), and with Goodman et al. (2013) for palms
   (which ignores wood density).
6. Corrects the 1985 small-stem DBH rounding (stems below 5.5 cm recorded in
   5 mm classes) with census 3 as reference.
7. Flags a growth interval as an outlier when the DBH growth is outside
   [`lower_ddbh_threshold`, `upper_ddbh_threshold`] (−2.5, 7.5 cm/yr) or the
   height of measurement moved by more than `hom_change_threshold` (1 m), and
   replaces its AGB change by the size-class relative rate (non-outlier change
   / AGB of the class, recruits excluded) × its own AGB.
8. Aggregates growth, recruitment and mortality by quadrat, interval and size
   class, and applies the Kohyama et al. (2019) bias correction to
   productivity and mortality per quadrat (`use_kohyama19_correction`).
9. Plots the plot-level and size-class dynamics.

### `basal_area_uncertainty.R`

1. Loads the R tables (the exported reconstruction), the posterior samples and
   the stage-2 table.
2. Treats palms and strangler figs as `biomass_stocks_fluxes.R` does: the
   giant strangler figs are removed with every stem of their tree (and leave
   the posterior paths too), and every measured alive palm stem other than
   *Socratea* gets its species' median DBH over all censuses. A palm's basal
   area is then constant, so identity uncertainty changes only its loss and
   gain.
3. Computes the stocks and fluxes of the exported reconstruction per tree and
   quadrat.
4. Separates trees with one sampled path (no identity uncertainty) from trees
   with several.
5. **Splice.** The paths cover the engine's window (at most up to the anchor
   census); every other measured observation of a sampled tree is added to
   every path with the stage-2 links: an observation outside the window joins
   the path stem that holds its exported stem-mate at the nearest window edge
   (the last in-path observation of the same exported stem for later censuses,
   the first one for earlier censuses); an exported stem without in-path
   observations keeps its own label. A tree with a measured observation missing
   inside its window, or without overlap, uses the exported reconstruction
   (fallback). Every interval and stock of a sampled tree is thus computed from
   one consistent set of stems.
6. **Monte Carlo.** `K_realizations` realizations (seed `mc_seed`); each draws
   one path per tree with several paths, independently, with weight
   `path_count / sum(path_count)`, and recomputes the stocks and fluxes.
7. Reports the exported reconstruction with the 95% Monte Carlo interval
   (centre line `mc_center`), the figures and the diagnostics.

Posterior weights and the two engines:

- Every tree has 200 draws, collapsed into unique paths; `path_count` is the
  number of draws of a path, so `path_count / 200` is its posterior
  probability and each draw is equally likely. `path_prob` is not used: for DP
  trees it re-weights each draw by its own probability, which counts that
  probability twice.
- DP trees are sampled by backward sampling from the exact DP posterior, so
  likely paths repeat.
- Trees of the probabilistic matcher (palm clumps and other species routed to
  it, strangler figs, trees the DP cannot solve because the state space is too
  large or no assignment is feasible) get approximate draws: a perturbed
  assignment per census pair that respects the database pins (two stems with
  different pins are never joined), stitched and repaired for growth
  violations. Nearly every draw differs, so nearly every path weighs 1/200;
  there is no single most probable path, and the exported reconstruction is
  the most representative draw.
- Both engines are sampled together, one draw per tree per realization. The
  matcher's draws approximate a posterior rather than calibrate it; the
  diagnostics report both engines separately (trees, unique paths, share of
  paths drawn once, how often the exported partition is among the sampled
  paths).

Scope:

- Only identity uncertainty is quantified, and only before the anchor census
  (`ANCHOR_START_CENSUS = 7`): after it identities are confirmed by the tags.
  Trees are drawn independently, so their variations largely cancel at plot
  level and the interval is narrow; it holds no sampling (quadrat) or model
  uncertainty.
- Apart from the palm and fig treatment and the taper correction, the choices
  of `biomass_stocks_fluxes.R` (outlier substitution, Kohyama correction) are
  not applied.

### `general_plot_information.R`

1. Loads the R tables and the species table. An individual is a `treeID`; an
   alive stem is `Rstatus == "A"`, except a stump census.
2. Uses the taper-corrected DBH (`use_taper_for_ba = TRUE`, so stand basal
   area equals the exported stock of `basal_area_uncertainty.R`; `FALSE` uses
   the DBH at the point of measurement), gives missed measurements a DBH from
   the stem's own measurements and places trees without coordinates at the
   centre of their quadrat.
3. Summarises the last census (individuals, stems, basal area, species,
   genera, families, life forms, diversity), per-hectare estimates by bootstrap
   over 20 m quadrats (`BOOT_R` replicates), trends and paired changes between
   censuses (census 1 left out), diversity across censuses, and the DBH size
   classes of censuses 2 and 9.

## Settings

| Setting | Script | Default | What it does |
| --- | --- | --- | --- |
| `remove_strangler_figs`, `use_median_palm_dbh` | biomass, basal area | `TRUE`, `TRUE` | Giant strangler figs removed; palm DBH set to the species median |
| `biomass_allometry`, `use_local_height_allometry` | biomass | `"chave14"`, `TRUE` | AGB allometry; Martínez Cano et al. (2019) height model |
| `use_kohyama19_correction` | biomass | `TRUE` | Kohyama et al. (2019) correction of productivity and mortality |
| `lower_ddbh_threshold`, `upper_ddbh_threshold`, `hom_change_threshold` | biomass | −2.5, 7.5 cm/yr; 1 m | Outlier growth intervals |
| `dbh_interp_method` | biomass | `"linear"` | DBH of missed measurements (`linear`, `locf`, `mean`) |
| `ANCHOR_START_CENSUS` | basal area | `7` | First census with confirmed identities |
| `K_realizations`, `mc_seed`, `mc_center` | basal area | `100`, `42`, `"median"` | Monte Carlo realizations, seed, centre line |
| `write_tree_realizations`, `write_quadrat_realizations` | basal area | `FALSE`, `TRUE` | Optional heavy outputs |
| `report_splice_methods` | basal area | `TRUE` | Report the stage-2 method of spliced observations |
| `first_plot_census` | basal area | `2` | First census in the figures |
| `use_taper_for_ba` | plot information | `TRUE` | Taper-corrected DBH for basal area and size classes |
| `BOOT_R` | plot information | `4999` | Bootstrap replicates |

## Checks

Each script prints ✓ or ❌ checks; a failed check prints ❌, raises a warning
and stops. `basal_area_uncertainty.R` checks, among others, that in every
realization growth + loss + gain equals the change of the plot stock, reports
the largest difference between a realization's plot stock and the exported
stock (only stems alive without a measurement can differ), and writes its
diagnostics to `ba_mc_diagnostics.txt`.

## Requirements

R packages: `data.table`, `truncnorm`, `lubridate`, `HDInterval`, `ggplot2`,
`scales`, `cowplot` (`biomass_stocks_fluxes.R`); `data.table`, `ggplot2`,
`patchwork`, `scales`, `collapse`, `arrow` (`basal_area_uncertainty.R`);
`data.table`, `boot`, `HDInterval` (`general_plot_information.R`).
