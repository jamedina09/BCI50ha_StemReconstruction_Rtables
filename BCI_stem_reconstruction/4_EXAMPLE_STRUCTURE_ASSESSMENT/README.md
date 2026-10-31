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
| `stem_dbh.R` | The diameter of a stem in a census, sourced by the three scripts: taper correction, stumps, the censuses to fill and their DBH, giant strangler figs |
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
| `outputs/plot_agb_fluxes_net.png` | Productivity and mortality (left axis) and the net AGB change (right axis) in one panel, with their 95% intervals across quadrats |
| `outputs/plot_agb_by_size.png` | The same by DBH size class |
| `outputs/ba_map_change_treeID.feather`, `ba_map_change_quadrat.feather` | BA fluxes per tree and per quadrat and census pair, exported reconstruction |
| `outputs/ba_map_stock_quadrat.feather` | BA stock per quadrat and census, exported reconstruction |
| `outputs/ba_mc_realizations_quadrat.feather`, `ba_mc_realizations_stock_quadrat.feather` | BA fluxes and stocks per quadrat in each Monte Carlo realization (`write_quadrat_realizations`) |
| `outputs/ba_mc_summary_quadrat.feather`, `ba_mc_summary_stock_quadrat.feather` | 95% interval of the quadrat fluxes and stocks |
| `outputs/ba_mc_realizations_treeID/` | One file per realization with tree-level fluxes, written only with `write_tree_realizations <- TRUE` (not under version control) |
| `outputs/fig1_BA_stock.pdf`, `fig2_BA_fluxes.pdf` | Plot BA stock and fluxes: exported reconstruction and Monte Carlo interval; the caption of `fig2` says why the exported line can lie outside the interval, with the shifts of the run |
| `outputs/fig3_BA_trajectories.pdf` | BA trajectories of a sample of trees: exported reconstruction and posterior paths |
| `outputs/ba_mc_diagnostics.txt` | Tree counts (one path, sampled, fallback; by engine), splice counts, the exported reconstruction against the posterior per engine, the exported reconstruction against the Monte Carlo (value, centre, interval and shift per interval and component; shift of the posterior mean and stems by engine), invariance checks and the Monte Carlo error of the interval bounds |

`biomass_stocks_fluxes.R` and `general_plot_information.R` write their checks
to `outputs/check_biomass_stocks_fluxes.txt` and
`outputs/check_general_plot_information.txt`; `general_plot_information.R`
writes no other file.

## What it does

The three scripts read the R tables in the same way, with the rules of
`stem_dbh.R`:

- A stem is `treeID` + `stemID` (`stemID` is the stem number within its tree).
- A missing `ExactDate` is completed from `date` (days since 1960-01-01,
  filled by stage 3 on every row), so every row has a date.
- **A stem is in the stand (stocks, fluxes, stem counts) in a census only
  when it has a DBH at 1.3 m**, measured or filled. The DBH that a `G` or `D`
  record can carry (stage 3 exports DBH as recorded) is ignored.
- **Taper.** Every DBH measured at another height than 1.3 m is corrected to
  1.3 m with Model 1 of Cushman et al. (2014), `D = d exp(b (h − 1.3))` (d
  measured at height h, in cm and m), with `log(b) = −2.0205 + species term −
  0.5053 log(d) + 0.3748 log(h)` from their Table 3; the species terms are
  those of their ten focal species, deviations from the base class of every
  other species (`use_taper_species_terms`). Palms keep the recorded DBH
  (`taper_palms = FALSE`): they have no buttress taper, and *Socratea* is
  measured above its stilt roots.
- **Taper for stocks and for fluxes.** The correction is set separately for
  the stocks (`use_taper_stocks`) and for the fluxes (`use_taper_fluxes`),
  both `TRUE` by default; without it the DBH is used as recorded, at its height
  of measurement. `general_plot_information.R` describes stocks only and has
  `use_taper_stocks` alone. With different settings the stocks and the fluxes
  are in different diameters: the change of the reported stock between two
  censuses is then not the net flux, and the scripts check their flux
  accounting, and apply the Kohyama correction, with the stock in the DBH of
  the fluxes.
- **Stump.** An `A` record without a DBH whose `DFstatus` is "broken below":
  the stem is alive but has no stem at 1.3 m (a broken-below record with a
  field tag is that stem). It gets no DBH and is not in the stand. It splits
  the stem's life into segments: the stem leaves the stand at the stump (its
  basal area or biomass is a loss; its `Rstatus` stays `A`, so it is not a
  death), and its next measurement is a recruit (a gain). An `A` record
  just before a stump with no sign of a stem at 1.3 m (no "alive" record, and
  no DBH on a record that is not dead: a dead, stem dead, missing or absent
  record, alive only because the stump comes later) is a stump too: the stem
  broke at a time not known between its last measurement and the stump; a DBH
  recorded on such a dead record is dropped, as for a dead stem.
- **Filled censuses.** An `A` record without a DBH, not a stump, that comes
  after a measurement of the stem with no stump in between, is a missed
  measurement: the stem stays in the stand, with a DBH from its own
  measurements of that segment, linear in time between two measurements, and
  after the last one the growth trend of the last two continued (a negative
  trend is not extrapolated; one measurement is kept).
- **Not in the stand** are, besides stumps, an `A` record before the stem's
  first measurement (the stem enters the stand when it is first measured),
  after a stump with no measurement since, and a stem never measured. A
  check stops a script if a census to fill lies between two measurements
  across a stump.

### `biomass_stocks_fluxes.R`

For each census interval c → c+1 it computes growth (AGB change of stems
alive at both ends), recruitment (AGB of a stem at its first alive census,
divided by the interval) and mortality (AGB at the last alive census, divided
by the interval), per stem, then aggregates them.

1. Loads the R tables of censuses 1–9, the species table and the wood density;
   a species without wood density takes the mean of its genus, else its
   family, else all species (each species counted once).
2. Removes the giant strangler figs (trees of *Ficus costaricana*,
   *obtusifolia*, *popenoei* or *trigonata* with a stem above 500 mm)
   (`remove_strangler_figs`).
3. Sets the DBH of the fluxes and the DBH of the stocks, each taper-corrected
   or as recorded (`use_taper_fluxes`, `use_taper_stocks`; above). The
   recorded value is kept in `dbh_raw`. Steps 4–6 apply to both. Growth,
   recruitment, mortality, the outlier rule, the Kohyama correction and the
   size classes of the fluxes use the DBH of the fluxes; the reported stocks
   (plot and size class) use the DBH of the stocks.
4. Converts DBH to cm and gives every alive palm stem, *Socratea* aside, its
   species' median DBH over all censuses (`use_median_palm_dbh`).
5. Fills the DBH of the censuses to fill (`dbh_interp_method`: `linear`, the
   rule of the other two scripts, `locf` or `mean`, within the stem's
   segment); nothing else gets a DBH.
6. Computes AGB with Chave et al. (2014) Eq. 4, `AGB = 0.0673 (ρ H D²)^0.976`
   (`biomass_allometry`; `chave05` available), with the height `H = 58.0
   D^0.73 / (21.8 + D^0.73)` of Martínez Cano et al. (2019) Eq. 7 (bias
   corrected for the back-transformation from log H;
   `use_local_height_allometry`), and with the diameter-only Arecaceae model
   of Goodman et al. (2013) for palms (which ignores wood density).
7. Corrects the 1985 small-stem DBH rounding (stems below 5.5 cm recorded in
   5 mm classes) with census 3 as reference.
8. Flags a growth interval as an outlier when the DBH growth is outside
   [`lower_ddbh_threshold`, `upper_ddbh_threshold`] (−2.5, 7.5 cm/yr) or the
   height of measurement moved by more than `hom_change_threshold` (1 m), and
   replaces its AGB change by the size-class relative rate (non-outlier change
   / AGB of the class, recruits excluded) × its own AGB.
9. Aggregates growth, recruitment and mortality by quadrat, interval and size
   class, and applies the Kohyama et al. (2019) bias correction to
   productivity and mortality per quadrat (`use_kohyama19_correction`).
10. Plots the plot-level and size-class dynamics.

### `basal_area_uncertainty.R`

1. Loads the R tables (the exported reconstruction), the posterior samples and
   the stage-2 table.
2. Treats palms, strangler figs and the taper correction as
   `biomass_stocks_fluxes.R` does: the giant strangler figs are removed with
   every stem of their tree (and leave the posterior paths too), the DBH of
   the fluxes and the DBH of the stocks are taper corrected or not
   (`use_taper_fluxes`, `use_taper_stocks`; palms aside), and every measured
   alive palm stem other than *Socratea* gets its species' median DBH over all
   censuses. A palm's basal area is then constant, so identity uncertainty
   changes only its loss and gain.
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

Why the exported line can lie outside the Monte Carlo interval (`fig2`,
before the anchor census, mostly in Gain and Loss):

- The bold line and the interval measure different things. The bold line is
  one reconstruction per tree: the most probable path of a DP tree, the most
  representative draw of a matcher tree. The interval is the spread of the
  plot total over realizations in which every tree takes a random posterior
  draw.
- The draws can split a stem history into two stems where the exported
  reconstruction keeps one stem, or keep one where it splits. A split turns
  one continuing stem into one that dies (its basal area counts as Loss) and
  a new one (its basal area counts as Gain); Growth loses only that stem's
  growth, and the net change stays the same. When the draws split more often
  than the exported reconstruction, the Monte Carlo Gain is higher and the
  Monte Carlo Loss more negative than the exported values, by about the same
  amount, while Growth and the net change barely move.
- The shift has the same sign in most trees, so it adds up over trees, while
  the random differences between draws cancel. The interval holds only that
  random spread, so it is narrower than the shift. A matcher tree's exported
  draw holds a small share of its posterior, so these trees can shift more
  than DP trees, whose exported path is usually the most probable one.
- After the anchor census every stem is tagged and the two agree.
- The caption of `fig2` gives the shifts of the run in percent of the
  exported value. `ba_mc_diagnostics.txt` gives, per interval and component,
  the exported value, the Monte Carlo centre and interval and the shift; the
  shift of the posterior mean by engine, computed exactly from the path
  weights; and the stems of the exported reconstruction against the posterior
  mean.
- The Monte Carlo centre is the expected flux under the model; the R tables
  carry one reconstruction per tree. Fluxes that should carry the identity
  uncertainty are computed over the posterior draws.

### `general_plot_information.R`

1. Loads the R tables and the species table. An individual is a `treeID`; the
   stems of a census are those with a DBH at 1.3 m (measured or filled).
2. Uses the taper-corrected DBH (`use_taper_stocks`; `FALSE`: the DBH as
   recorded) and the filled censuses of the other two scripts, and places
   trees without coordinates at the centre of their quadrat. It describes every
   tree of the plot: the giant strangler figs and the measured palm diameters
   stay, so its stand basal area is not the stock of
   `basal_area_uncertainty.R`.
3. Summarises the last census (individuals, stems, basal area, species,
   genera, families, life forms, diversity), per-hectare estimates by bootstrap
   over 20 m quadrats (`BOOT_R` replicates), trends and paired changes between
   censuses (census 1 left out), diversity across censuses, and the DBH size
   classes of censuses 2 and 9.

## Settings

| Setting | Script | Default | What it does |
| --- | --- | --- | --- |
| `remove_strangler_figs`, `use_median_palm_dbh` | biomass, basal area | `TRUE`, `TRUE` | Giant strangler figs removed; palm DBH set to the species median |
| `use_taper_stocks` | biomass, basal area, plot information | `TRUE` | Stocks from the taper-corrected DBH (`FALSE`: the DBH as recorded, at its height of measurement) |
| `use_taper_fluxes` | biomass, basal area | `TRUE` | Fluxes (growth, recruitment, mortality) from the taper-corrected DBH (`FALSE`: the DBH as recorded) |
| `use_taper_species_terms`, `taper_palms` | biomass, basal area, plot information | `TRUE`, `FALSE` | Taper model with the species terms of Cushman et al. (2014); palms keep the recorded DBH (`FALSE`, `TRUE`: the base class for every species and palms corrected too, as Piponiot et al. 2024) |
| `biomass_allometry`, `use_local_height_allometry` | biomass | `"chave14"`, `TRUE` | AGB allometry; Martínez Cano et al. (2019) height model |
| `use_kohyama19_correction` | biomass | `TRUE` | Kohyama et al. (2019) correction of productivity and mortality |
| `lower_ddbh_threshold`, `upper_ddbh_threshold`, `hom_change_threshold` | biomass | −2.5, 7.5 cm/yr; 1 m | Outlier growth intervals |
| `dbh_interp_method` | biomass | `"linear"` | DBH of the censuses to fill (`linear`, the rule of the other two scripts; `locf`, `mean`) |
| `ANCHOR_START_CENSUS` | basal area | `7` | First census with confirmed identities |
| `K_realizations`, `mc_seed`, `mc_center` | basal area | `100`, `42`, `"median"` | Monte Carlo realizations, seed, centre line |
| `write_tree_realizations`, `write_quadrat_realizations` | basal area | `FALSE`, `TRUE` | Optional heavy outputs |
| `report_splice_methods` | basal area | `TRUE` | Report the stage-2 method of spliced observations |
| `first_plot_census` | basal area | `2` | First census in the figures |
| `BOOT_R` | plot information | `4999` | Bootstrap replicates |

## Checks

Each script prints ✓ or ❌ checks; a failed check prints ❌, raises a warning,
writes to the script's check file and stops.

Each script checks that no census to fill lies between two measurements
across a stump, that every census to fill gets a DBH, and that no stump or
other record outside the stand gets one.

Each script first checks its inputs: the nine R tables and the species table
exist and are newer than the stage-2 table (otherwise rerun stage 3 after the
merge). `basal_area_uncertainty.R` also checks that the posterior samples are
newer than the stage-2 table and come from the same engine run (the
`engine_run` attribute of the posterior file against
`DATA/PROCESSED/engine_run.txt`).

`biomass_stocks_fluxes.R` and `basal_area_uncertainty.R` check that the
stocks and the fluxes hold the same stem-censuses, and that both use the same
DBH when their taper settings are equal. `basal_area_uncertainty.R` checks,
among others, that in every realization growth + loss + gain equals the
change of the plot stock in the DBH of the fluxes, reports
the largest difference between a realization's plot stock and the exported
stock (only stems alive without a measurement can differ), and writes its
diagnostics to `ba_mc_diagnostics.txt`.

## Requirements

R packages: `data.table`, `truncnorm`, `lubridate`, `HDInterval`, `ggplot2`,
`scales`, `cowplot` (`biomass_stocks_fluxes.R`); `data.table`, `ggplot2`,
`patchwork`, `scales`, `collapse`, `arrow` (`basal_area_uncertainty.R`);
`data.table`, `boot`, `HDInterval` (`general_plot_information.R`).
