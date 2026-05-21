## ============================================================
## biomass_uncertainty.R
## Biomass stocks and fluxes with posterior path uncertainty
## from multi-stem identification (DP reconstruction)
## ============================================================
##
## Stem identity in multi-stem trees is uncertain: different path
## realizations assign census observations to different physical
## stems, changing growth/mortality decomposition and thus AGB fluxes.
##
## Approach:
##   1. Pre-process all 9 censuses (same pipeline as BIOMASS_STOCKS_FLUXES.R):
##      taper correction → dbh_t, AGB allometry → agb_t.
##   2. Build an "observation lookup" from rows that carry an obs_row_id
##      (multi-stem trees flagged by the DP reconstruction).
##   3. Load posterior path samples (post_file).
##   4. Parse multi-path tag recons (vectorised) → per-(tag, path) stem
##      trajectories. Apply interpolation for any alive stem with NA dbh_t.
##   5. Pre-compute AGB stock and flux contributions per (tag, path_idx).
##   6. Single-stem trees and single-path multistem tags form the fixed
##      baseline (not subject to path uncertainty).
##   7. Draw K Monte Carlo realizations: for each multi-path tag, sample
##      one path proportional to path_prob; aggregate to plot level.
##   8. Summarise: mean ± 95% CI across K realizations.
##   9. Write CSV summaries and a diagnostic figure.
##
## Key simplifications relative to BIOMASS_STOCKS_FLUXES.R:
##   - No outlier substitution (path uncertainty is the focus here).
##   - No Kohyama et al. 2019 correction.
##   - No 1985 DBH-rounding bias correction.
##
## Outputs (written to ./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/):
##   agb_uncertainty_stock_realizations.csv — per-census × K realizations
##   agb_uncertainty_flux_realizations.csv  — per-interval × K realizations
##   agb_uncertainty_stock_summary.csv      — mean ± 95% CI for AGB stock
##   agb_uncertainty_flux_summary.csv       — mean ± 95% CI for AGB fluxes
##   agb_uncertainty_figure.png             — stock, productivity, mortality, net
## ============================================================

rm(list = ls())

library(data.table)
library(ggplot2)
library(cowplot)

## --- 0. Parameters -----------------------------------------------------------

K_realizations        <- 200      # Monte Carlo realizations
dbh_interp_method     <- "linear" # "linear", "locf", or "mean"
seed_val              <- 42
n_ha                  <- 50       # BCI plot area (ha)

# Same biological / QC flags as BIOMASS_STOCKS_FLUXES.R
remove_strangler_figs <- TRUE
use_median_palm_dbh   <- TRUE

set.seed(seed_val)

## --- 1. Load census data -----------------------------------------------------

bci_stem_nums <- as.character(1:9)

census_list <- lapply(bci_stem_nums, function(num) {
  fp <- paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", num, ".Rdata")
  if (!file.exists(fp)) stop("Missing census file: ", fp)
  load(fp)
  data.table(get(paste0("bci.stem", num)))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)

df_stem <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rm(census_list, bci_stem_nums)
df_stem[, CensusID := as.integer(CensusID)]

## --- 2. Taxonomy and wood density --------------------------------------------

bci_wd <- fread("./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/wd/doi_10_5061_dryad_5qfttdzn3__v20260403/WD_species.txt")
bci_wd <- bci_wd[, .(sp = tolower(sp6), wsg = wd100.mean)][!is.na(wsg)]

bci_spp <- fread("./BCI_stem_reconstruction/DATA/RTABLES/bci.spptable.csv")
bci_spp <- unique(bci_spp[, .(Family, sp, Genus, Latin)])

df_stem <- merge(df_stem, bci_spp, by = "sp", all.x = TRUE)
df_stem <- merge(df_stem, bci_wd,  by = "sp", all.x = TRUE)
rm(bci_wd, bci_spp)

# WSG: fill NAs hierarchically (genus → family → global mean)
df_stem[, wsg := ifelse(is.na(wsg), ave(wsg, Genus,  FUN = \(x) mean(x, na.rm = TRUE)), wsg)]
df_stem[, wsg := ifelse(is.na(wsg), ave(wsg, Family, FUN = \(x) mean(x, na.rm = TRUE)), wsg)]
df_stem[, wsg := ifelse(is.na(wsg), mean(wsg, na.rm = TRUE), wsg)]
df_stem[is.na(Family), Family := "Unknown"]
df_stem[is.na(Genus),  Genus  := "Unknown"]

## --- 3. Strangler fig removal ------------------------------------------------

if (remove_strangler_figs) {
  large_figs <- df_stem[
    dbh > 500 & Genus == "Ficus" &
      !is.na(Latin) & grepl("costaricana|obtusifolia|popenoei|trigonata", Latin)
  ]
  df_stem <- df_stem[
    !treeID %in% large_figs$treeID & !stemID %in% large_figs$stemID
  ]
  rm(large_figs)
}

## --- 4. DBH mm → cm, palm correction, taper correction → dbh_t, agb_t ------

# DBH: mm → cm (raw BCI units)
df_stem[, dbh := dbh / 10]

# Palm DBH: replace by species median (excluding Socratea; see BIOMASS_STOCKS_FLUXES.R)
if (use_median_palm_dbh) {
  df_stem[Family == "Arecaceae" & Genus != "Socratea",
          dbh := median(dbh, na.rm = TRUE), .(Latin)]
}

# Taper correction (Cushman et al. 2014): corrects DBH to 1.3 m equivalent.
# hom = NA stems get hom = 1.3 m → b ≈ 0 → dbh_t = dbh (no correction).
df_stem[, hom := ifelse(is.na(hom), 1.3, hom)]
df_stem[, b   := exp(-2.0205 - 0.5053 * log(dbh) + 0.3748 * log(hom))]
df_stem[!is.na(hom), dbh_t := dbh * exp(b * (hom - 1.3))]
df_stem[, b := NULL]

# AGB allometry:
#   Chave et al. 2014 (eq. 4) + Martinez-Cano et al. 2019 height model.
#   Palms: Goodman et al. 2013 (palm-specific).
agb_bci <- function(dbh, wsg, is_palm = NULL) {
  h   <- 58.0 * dbh^0.73 / (21.8 + dbh^0.73)
  agb <- (0.0673 * (wsg * h * dbh^2)^0.976) / 1000  # Mg
  if (!is.null(is_palm) && any(is_palm, na.rm = TRUE)) {
    agb[is_palm] <- 0.0417565 * dbh[is_palm]^2.7483 / 1000
  }
  agb
}

df_stem[, agb_t := agb_bci(dbh_t, wsg, is_palm = (Family == "Arecaceae"))]

## --- 5. Census date lookup ---------------------------------------------------

# Plot-median ExactDate per CensusID → CensusYear and inter-census interval dT.
# Convert ExactDate to numeric before median, then restore to Date to keep types consistent.
census_yr_lut <- df_stem[!is.na(ExactDate),
  .(date_med   = as.Date(median(as.numeric(ExactDate)), origin = "1970-01-01")),
  by = CensusID
]
census_yr_lut[, CensusYear := as.integer(format(date_med, "%Y"))]
setorder(census_yr_lut, CensusID)
census_yr_lut[, dT := c(as.numeric(diff(date_med), units = "days") / 365.25, NA_real_)]
census_yr_lut[, date_med := NULL]

max_census <- census_yr_lut[, max(CensusID)]

## --- 6. Separate multistem observations from single-stem trees ---------------

# obs_data: rows with obs_row_id — multistem observations assigned by the DP
# reconstruction. obs_row_id is an integer index within each tag (starts at 1).
# The posterior path maps obs_row_id → virtual StemID for each census.
obs_data <- df_stem[!is.na(obs_row_id), .(
  tag, treeID, CensusID, obs_row_id,
  dbh_t, agb_t, Rstatus, wsg, quadrat, Family
)]

# Single-stem trees (obs_row_id = NA): unambiguous identity; no path uncertainty.
# Keep only alive stems with valid AGB.
df_single <- df_stem[
  is.na(obs_row_id) & !is.na(dbh_t) & Rstatus == "A",
  .(treeID, stemID, CensusID, dbh_t, agb_t, quadrat)
]

rm(df_stem)
gc()

cat(sprintf("[UNCERT] Multistem obs: %d rows, %d tags\n",
  nrow(obs_data), uniqueN(obs_data$tag)))
cat(sprintf("[UNCERT] Single-stem alive obs: %d rows\n", nrow(df_single)))

## --- 7. Load posterior file --------------------------------------------------

post_file <- "./BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds"
post <- data.table(readRDS(post_file))
post[, path_idx := seq_len(.N), by = Tag]  # 1, 2, ... per tag

cat(sprintf("[UNCERT] Posterior: %d paths across %d tags\n",
  nrow(post), uniqueN(post$Tag)))

# Split into single-path (fixed) and multi-path (uncertain) tags.
n_per_tag       <- post[, .N, by = Tag]
single_path_tags <- n_per_tag[N == 1L, Tag]
multi_path_tags  <- n_per_tag[N >  1L, Tag]
cat(sprintf("[UNCERT] Single-path tags: %d | Multi-path tags: %d\n",
  length(single_path_tags), length(multi_path_tags)))

## --- 8. Vectorised recon parsing ---------------------------------------------

# Recon format: "obs_row_id:StemID;obs_row_id:StemID;..."
# obs_row_id is the within-tag index into obs_data (not globally unique).
# StemID is the virtual stem identifier assigned by the DP reconstruction.

parse_post_subset <- function(post_sub) {
  if (nrow(post_sub) == 0L) return(data.table())
  pairs_list <- strsplit(post_sub$recon, ";", fixed = TRUE)
  n_pairs    <- lengths(pairs_list)
  dt <- data.table(
    Tag       = rep(post_sub$Tag,       n_pairs),
    path_idx  = rep(post_sub$path_idx,  n_pairs),
    path_prob = rep(post_sub$path_prob, n_pairs),
    pair      = unlist(pairs_list, use.names = FALSE)
  )
  dt[, c("obs_row_id", "StemID") :=
    tstrsplit(pair, ":", fixed = TRUE, type.convert = TRUE)]
  dt[, pair := NULL]
  dt
}

cat("[UNCERT] Parsing multi-path recons...\n")
recon_multi <- parse_post_subset(post[Tag %in% multi_path_tags])
cat(sprintf("[UNCERT]   Multi-path: %d assignments\n", nrow(recon_multi)))

cat("[UNCERT] Parsing single-path recons...\n")
recon_single_path <- parse_post_subset(post[Tag %in% single_path_tags])
cat(sprintf("[UNCERT]   Single-path: %d assignments\n", nrow(recon_single_path)))

## --- 9. Join recon with observation data -------------------------------------

join_obs <- function(recon_dt) {
  if (nrow(recon_dt) == 0L) return(recon_dt)
  merge(
    recon_dt,
    obs_data[, .(Tag = tag, obs_row_id, CensusID, dbh_t, agb_t,
                 Rstatus, wsg, quadrat, Family)],
    by    = c("Tag", "obs_row_id"),
    all.x = TRUE
  )
}

recon_multi       <- join_obs(recon_multi)
recon_single_path <- join_obs(recon_single_path)

## --- 10. Interpolate missing dbh_t within reconstructed trajectories ---------

# After path assignment, a stem's trajectory may have gaps:
# a census observation was assigned to another StemID in the selected path,
# leaving this StemID with Rstatus = "A" but dbh_t = NA (alive but unmeasured
# in that census under this path assumption). Interpolate linearly between
# flanking measured values, consistent with BIOMASS_STOCKS_FLUXES.R.

interp_fill <- function(dbh_vec, method) {
  # Fills NA values that are between two measured values (internal gaps).
  # Called per (Tag, path_idx, StemID) ordered by CensusID.
  measured <- !is.na(dbh_vec)
  if (sum(measured) < 2L) return(dbh_vec)
  target <- is.na(dbh_vec)
  if (!any(target)) return(dbh_vec)
  idx_m <- which(measured)
  out   <- dbh_vec
  if (method == "linear") {
    yhat <- approx(idx_m, dbh_vec[idx_m], xout = seq_along(dbh_vec), rule = 1)$y
    out[target] <- yhat[target]
  } else if (method == "mean") {
    for (i in which(target)) {
      prev_i <- if (any(idx_m < i)) max(idx_m[idx_m < i]) else NA_integer_
      next_i <- if (any(idx_m > i)) min(idx_m[idx_m > i]) else NA_integer_
      vals <- c(if (!is.na(prev_i)) dbh_vec[prev_i],
                if (!is.na(next_i)) dbh_vec[next_i])
      if (length(vals)) out[i] <- mean(vals)
    }
  } else { # locf
    lv <- NA_real_
    for (i in seq_along(dbh_vec)) {
      if (!is.na(dbh_vec[i])) lv <- dbh_vec[i]
      else if (!is.na(lv))    out[i] <- lv
    }
  }
  out
}

interp_stems <- function(dt, method) {
  n_na <- dt[!is.na(Rstatus) & Rstatus == "A" & is.na(dbh_t), .N]
  if (n_na == 0L) {
    message("[INTERP] No alive stems with NA dbh_t in this subset.")
    return(invisible(dt))
  }
  setorder(dt, Tag, path_idx, StemID, CensusID)
  dt[, dbh_t := interp_fill(dbh_t, method = method), by = .(Tag, path_idx, StemID)]
  # Recompute agb_t where interpolation filled in a dbh_t.
  dt[!is.na(dbh_t) & is.na(agb_t),
    agb_t := agb_bci(dbh_t, wsg, is_palm = (Family == "Arecaceae"))
  ]
  n_filled <- n_na - dt[!is.na(Rstatus) & Rstatus == "A" & is.na(dbh_t), .N]
  message(sprintf("[INTERP] Filled %d / %d alive stems with NA dbh_t.", n_filled, n_na))
  invisible(dt)
}

interp_stems(recon_multi,       method = dbh_interp_method)
interp_stems(recon_single_path, method = dbh_interp_method)

## --- 11. Pre-compute AGB stock and flux per (Tag, path_idx) ------------------

# For each (tag, path): compute AGB stock per census and flux per interval.
# This pre-computation is done once; the realization loop then just sums up
# the sampled paths.
#
# Productivity convention: growth of survivors + recruitment, re-indexed to the
# STARTING census of each interval (= CensusID - 1). Consistent with the
# `build_demo()` convention in BIOMASS_STOCKS_FLUXES.R.
#
# Mortality convention: for each StemID, its last alive census in the path
# determines the mortality event. The mortality flux is agb_t[last_c] / dT[last_c],
# attributed to the interval starting at last_c.
#
# [NOTE] Negative growth (shrinkage) is included in the productivity sum,
# consistent with how Dagb_t_s is aggregated in build_demo(). No outlier
# substitution is applied here.

compute_contributions <- function(recon_dt, dT_lut = census_yr_lut) {
  # Restrict to alive observations with valid AGB.
  alive <- recon_dt[!is.na(agb_t) & (is.na(Rstatus) | Rstatus == "A")]
  setorder(alive, Tag, path_idx, StemID, CensusID)

  ## --- Stock ----------------------------------------------------------------
  stock <- alive[,
    .(agb_stock = sum(agb_t, na.rm = TRUE)),
    by = .(Tag, path_idx, path_prob, CensusID, quadrat)
  ]

  ## --- Productivity (growth + recruitment) ----------------------------------
  alive[, prev_agb   := shift(agb_t,   type = "lag"), by = .(Tag, path_idx, StemID)]
  alive[, prev_census := shift(CensusID, type = "lag"), by = .(Tag, path_idx, StemID)]

  # Merge dT at the END census of the interval.
  alive <- merge(alive, dT_lut[, .(CensusID, dT)], by = "CensusID", all.x = TRUE)

  # Annualised growth for survivors (stem present in consecutive censuses).
  alive[, Dagb := fifelse(
    !is.na(prev_agb) & prev_census == CensusID - 1L & !is.na(dT) & dT > 0,
    (agb_t - prev_agb) / dT,
    NA_real_
  )]
  # Recruitment gain: first appearance of StemID in this path.
  alive[is.na(prev_agb) & !is.na(dT) & dT > 0, Dagb := agb_t / dT]

  # Productivity: re-indexed to the STARTING census (CensusID - 1L).
  prod <- alive[!is.na(Dagb),
    .(prod = sum(Dagb, na.rm = TRUE)),
    by = .(Tag, path_idx, path_prob, quadrat, int_from = CensusID - 1L)
  ]

  ## --- Mortality ------------------------------------------------------------
  # Last alive census for each StemID in this path.
  # Mortality flux = agb_t[last_c] / dT[last_c], where dT is the time to
  # the NEXT plot-wide census. Stems surviving to the final census are NOT
  # counted as dead (last_c == max_census → excluded).
  stem_last <- alive[,
    {
      last_row <- .N  # already sorted by CensusID (ascending)
      .(last_c = CensusID[last_row], agb_last = agb_t[last_row],
        quadrat = quadrat[last_row])
    },
    by = .(Tag, path_idx, path_prob, StemID)
  ]
  stem_last <- stem_last[last_c < max_census]
  stem_last <- merge(stem_last, dT_lut[, .(CensusID, dT)],
    by.x = "last_c", by.y = "CensusID", all.x = TRUE)
  stem_last[, DagbM := fifelse(!is.na(dT) & dT > 0, agb_last / dT, NA_real_)]

  mort <- stem_last[!is.na(DagbM),
    .(mort = sum(DagbM, na.rm = TRUE)),
    by = .(Tag, path_idx, path_prob, quadrat, int_from = last_c)
  ]

  list(stock = stock, prod = prod, mort = mort)
}

cat("[UNCERT] Pre-computing multi-path contributions...\n")
multi_contrib  <- compute_contributions(recon_multi)

cat("[UNCERT] Pre-computing single-path contributions...\n")
single_contrib <- compute_contributions(recon_single_path)

# Free large recon tables — no longer needed.
rm(recon_multi, recon_single_path)
gc()

## --- 12. Fixed baseline: single-stem trees -----------------------------------

# Stock, productivity, and mortality from trees with no obs_row_id.
# These contributions are constant across all realizations.

setorder(df_single, treeID, stemID, CensusID)
df_single <- merge(df_single, census_yr_lut[, .(CensusID, dT)],
  by = "CensusID", all.x = TRUE)

df_single[, prev_agb    := shift(agb_t,   type = "lag"), by = .(treeID, stemID)]
df_single[, prev_census := shift(CensusID, type = "lag"), by = .(treeID, stemID)]
df_single[, Dagb := fifelse(
  !is.na(prev_agb) & prev_census == CensusID - 1L & !is.na(dT) & dT > 0,
  (agb_t - prev_agb) / dT, NA_real_
)]
df_single[is.na(prev_agb) & !is.na(dT) & dT > 0, Dagb := agb_t / dT]

baseline_stock <- df_single[!is.na(agb_t),
  .(agb_stock = sum(agb_t, na.rm = TRUE)), by = .(CensusID, quadrat)
]
baseline_prod <- df_single[!is.na(Dagb),
  .(prod = sum(Dagb, na.rm = TRUE)), by = .(quadrat, int_from = CensusID - 1L)
]

# Mortality from single-stem trees: last alive census approach.
single_last <- df_single[!is.na(agb_t),
  .(last_c = max(CensusID), agb_last = agb_t[which.max(CensusID)]),
  by = .(treeID, stemID, quadrat)
]
single_last <- single_last[last_c < max_census]
single_last <- merge(single_last, census_yr_lut[, .(CensusID, dT)],
  by.x = "last_c", by.y = "CensusID", all.x = TRUE)
single_last[, DagbM := fifelse(!is.na(dT) & dT > 0, agb_last / dT, NA_real_)]
baseline_mort <- single_last[!is.na(DagbM),
  .(mort = sum(DagbM, na.rm = TRUE)), by = .(quadrat, int_from = last_c)
]
rm(df_single, single_last)
gc()

## --- 13. Combine fixed baseline (single-stem + single-path multistem) --------

# Single-path multistem tags have only one path (path_prob = 1).
# Their contribution is fixed — aggregate once.
fixed_sp_stock <- single_contrib$stock[,
  .(agb_stock = sum(agb_stock, na.rm = TRUE)), by = .(CensusID, quadrat)
]
fixed_sp_prod <- single_contrib$prod[,
  .(prod = sum(prod, na.rm = TRUE)), by = .(quadrat, int_from)
]
fixed_sp_mort <- single_contrib$mort[,
  .(mort = sum(mort, na.rm = TRUE)), by = .(quadrat, int_from)
]
rm(single_contrib)

# Plot-level fixed totals (summed across all quadrats and converted to Mg ha⁻¹).
combine_and_sum <- function(a, b, by_col, val_col) {
  both <- rbindlist(list(a, b), use.names = TRUE, fill = TRUE)
  both[, lapply(.SD, sum, na.rm = TRUE), by = by_col, .SDcols = val_col]
}

fixed_stock_plot <- combine_and_sum(
  baseline_stock[, .(CensusID, agb_stock)],
  fixed_sp_stock[, .(CensusID, agb_stock)],
  "CensusID", "agb_stock"
)
fixed_prod_plot <- combine_and_sum(
  baseline_prod[, .(int_from, prod)],
  fixed_sp_prod[, .(int_from, prod)],
  "int_from", "prod"
)
fixed_mort_plot <- combine_and_sum(
  baseline_mort[, .(int_from, mort)],
  fixed_sp_mort[, .(int_from, mort)],
  "int_from", "mort"
)
rm(baseline_stock, baseline_prod, baseline_mort, fixed_sp_stock, fixed_sp_prod, fixed_sp_mort)

## --- 14. Variable contributions for multi-path tags -------------------------

# Indexed by (Tag, path_idx) for fast lookup in the realization loop.
var_stock <- multi_contrib$stock[,
  .(agb_stock = sum(agb_stock, na.rm = TRUE)), by = .(Tag, path_idx, CensusID)
]
var_prod <- multi_contrib$prod[,
  .(prod = sum(prod, na.rm = TRUE)), by = .(Tag, path_idx, int_from)
]
var_mort <- multi_contrib$mort[,
  .(mort = sum(mort, na.rm = TRUE)), by = .(Tag, path_idx, int_from)
]
rm(multi_contrib)
gc()

setkey(var_stock, Tag, path_idx)
setkey(var_prod,  Tag, path_idx)
setkey(var_mort,  Tag, path_idx)

## --- 15. K Monte Carlo realizations ------------------------------------------

# For each realization, sample one path per multi-path tag (proportional to
# path_prob), look up the pre-computed contributions, and add to the fixed
# baseline. Result is plot-level AGB stock and flux (Mg ha⁻¹ or Mg ha⁻¹ yr⁻¹).

cat(sprintf("[UNCERT] Running %d Monte Carlo realizations...\n", K_realizations))

path_opts <- post[Tag %in% multi_path_tags, .(Tag, path_idx, path_prob)]

# Pre-sample all paths for every tag once, then reuse these draws inside the loop.
cat("[UNCERT] Pre-sampling path choices for all realizations...\n")
sampled_paths <- path_opts[, .(
  path_idx    = sample(path_idx, K_realizations, replace = TRUE, prob = path_prob)
), by = Tag]
sampled_paths[, realization := seq_len(.N), by = Tag]
setkey(sampled_paths, realization)
rm(path_opts)
gc()

realization_stock <- vector("list", K_realizations)
realization_flux  <- vector("list", K_realizations)

for (k in seq_len(K_realizations)) {
  if (k %% 50L == 0L || k == 1L)
    cat(sprintf("  Realization %d / %d\n", k, K_realizations))

  # Use the pre-drawn path assignment for this realization.
  sampled <- sampled_paths[.(k), .(Tag, path_idx), on = "realization"]

  # Look up variable contributions for sampled paths
  sel <- sampled[, .(Tag, path_idx)]

  k_stock_var <- var_stock[sel, on = c("Tag", "path_idx"), nomatch = 0L][,
    .(agb_stock_var = sum(agb_stock, na.rm = TRUE)), by = CensusID
  ]
  k_prod_var <- var_prod[sel, on = c("Tag", "path_idx"), nomatch = 0L][,
    .(prod_var = sum(prod, na.rm = TRUE)), by = int_from
  ]
  k_mort_var <- var_mort[sel, on = c("Tag", "path_idx"), nomatch = 0L][,
    .(mort_var = sum(mort, na.rm = TRUE)), by = int_from
  ]

  # Total = fixed + variable (Mg ha⁻¹: divide by n_ha)
  total_stock <- merge(fixed_stock_plot, k_stock_var, by = "CensusID", all = TRUE)
  total_stock[is.na(agb_stock),   agb_stock   := 0]
  total_stock[is.na(agb_stock_var), agb_stock_var := 0]
  total_stock[, agb_ha := (agb_stock + agb_stock_var) / n_ha]
  total_stock[, realization := k]

  total_prod <- merge(fixed_prod_plot, k_prod_var, by = "int_from", all = TRUE)
  total_prod[is.na(prod), prod := 0]
  total_prod[is.na(prod_var), prod_var := 0]
  total_prod[, prod_ha := (prod + prod_var) / n_ha]

  total_mort <- merge(fixed_mort_plot, k_mort_var, by = "int_from", all = TRUE)
  total_mort[is.na(mort), mort := 0]
  total_mort[is.na(mort_var), mort_var := 0]
  total_mort[, mort_ha := (mort + mort_var) / n_ha]

  flux_k <- merge(total_prod[, .(int_from, prod_ha)],
                  total_mort[, .(int_from, mort_ha)],
                  by = "int_from", all = TRUE)
  flux_k <- merge(flux_k, census_yr_lut[, .(CensusID, CensusYear, dT)],
                  by.x = "int_from", by.y = "CensusID", all.x = TRUE)
  flux_k[, net_ha      := prod_ha - mort_ha]
  flux_k[, realization := k]

  realization_stock[[k]] <- total_stock[, .(CensusID, agb_ha, realization)]
  realization_flux[[k]]  <- flux_k[, .(int_from, CensusYear, dT,
                                        prod_ha, mort_ha, net_ha, realization)]
}

all_stock <- rbindlist(realization_stock)
all_flux  <- rbindlist(realization_flux)

## --- 16. Summarise across realizations ---------------------------------------

stock_summary <- all_stock[,
  .(agb_mean = mean(agb_ha),
    agb_sd   = sd(agb_ha),
    agb_lwr  = quantile(agb_ha, 0.025),
    agb_upr  = quantile(agb_ha, 0.975)),
  by = CensusID
]
stock_summary <- merge(stock_summary, census_yr_lut[, .(CensusID, CensusYear)],
  by = "CensusID", all.x = TRUE)
setorder(stock_summary, CensusID)

flux_summary <- all_flux[,
  .(prod_mean = mean(prod_ha),  prod_sd = sd(prod_ha),
    prod_lwr  = quantile(prod_ha, 0.025), prod_upr = quantile(prod_ha, 0.975),
    mort_mean = mean(mort_ha),  mort_sd = sd(mort_ha),
    mort_lwr  = quantile(mort_ha, 0.025), mort_upr = quantile(mort_ha, 0.975),
    net_mean  = mean(net_ha),   net_sd  = sd(net_ha),
    net_lwr   = quantile(net_ha, 0.025),  net_upr  = quantile(net_ha, 0.975)),
  by = .(int_from, CensusYear, dT)
]
setorder(flux_summary, int_from)
flux_summary[, x_mid := CensusYear + dT / 2]  # midpoint of the census interval

cat("\n===== AGB STOCK (Mg ha\u207b\u00b9) =====\n")
print(stock_summary[, .(CensusID, CensusYear,
  agb_mean = round(agb_mean, 1),
  agb_lwr  = round(agb_lwr,  1),
  agb_upr  = round(agb_upr,  1)
)])

cat("\n===== AGB FLUX (Mg ha\u207b\u00b9 yr\u207b\u00b9) =====\n")
print(flux_summary[, .(int_from, CensusYear,
  prod = round(prod_mean, 3),
  mort = round(mort_mean, 3),
  net  = round(net_mean,  3)
)])

## --- 17. Write CSV outputs ---------------------------------------------------

out_dir <- "./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES"

fwrite(all_stock,     file.path(out_dir, "agb_uncertainty_stock_realizations.csv"))
fwrite(all_flux,      file.path(out_dir, "agb_uncertainty_flux_realizations.csv"))
fwrite(stock_summary, file.path(out_dir, "agb_uncertainty_stock_summary.csv"))
fwrite(flux_summary,  file.path(out_dir, "agb_uncertainty_flux_summary.csv"))

cat("[UNCERT] CSV outputs written to", out_dir, "\n")

## --- 18. Figures -------------------------------------------------------------

# Figure A: AGB stock (Mg ha⁻¹)
p_stock <- ggplot(stock_summary, aes(x = CensusYear)) +
  geom_ribbon(aes(ymin = agb_lwr, ymax = agb_upr), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = agb_mean),  linewidth = 0.9, colour = "steelblue") +
  geom_point(aes(y = agb_mean), size = 2.5, colour = "steelblue") +
  labs(
    x = "Census year",
    y = expression("AGB (Mg ha"^{-1} * ")"),
    title = "AGB stock — path identity uncertainty (95% CI)"
  ) +
  theme_bw()

# Figure B: Productivity and mortality (Mg ha⁻¹ yr⁻¹)
p_flux <- ggplot(flux_summary, aes(x = x_mid)) +
  geom_ribbon(aes(ymin = prod_lwr, ymax = prod_upr),
              fill = "steelblue", alpha = 0.25) +
  geom_ribbon(aes(ymin = mort_lwr, ymax = mort_upr),
              fill = "firebrick", alpha = 0.25) +
  geom_line(aes(y = prod_mean, colour = "Productivity"), linewidth = 0.8) +
  geom_line(aes(y = mort_mean, colour = "Mortality"),    linewidth = 0.8) +
  geom_point(aes(y = prod_mean, colour = "Productivity"), size = 2) +
  geom_point(aes(y = mort_mean, colour = "Mortality"),    size = 2) +
  scale_colour_manual(values = c(Productivity = "steelblue", Mortality = "firebrick")) +
  labs(
    x = "Census year",
    y = expression("Flux (Mg ha"^{-1} ~ "yr"^{-1} * ")"),
    colour = NULL,
    title = "AGB fluxes — path identity uncertainty (95% CI)"
  ) +
  theme_bw()

# Figure C: Net AGB change (Mg ha⁻¹ yr⁻¹)
p_net <- ggplot(flux_summary, aes(x = x_mid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = net_lwr, ymax = net_upr), fill = "forestgreen", alpha = 0.25) +
  geom_line(aes(y = net_mean),  linewidth = 0.8, colour = "forestgreen") +
  geom_point(aes(y = net_mean), size = 2, colour = "forestgreen") +
  labs(
    x = "Census year",
    y = expression("Net AGB change (Mg ha"^{-1} ~ "yr"^{-1} * ")"),
    title = "Net AGB change — path identity uncertainty (95% CI)"
  ) +
  theme_bw()

fig <- plot_grid(p_stock, p_flux, p_net, ncol = 1, align = "v", labels = c("A", "B", "C"))

ggsave(
  plot = fig,
  filename = file.path(out_dir, "agb_uncertainty_figure.png"),
  width = 9, height = 12, units = "in", dpi = 300
)
cat("[UNCERT] Figure written.\n")
