## ============================================================
## BCI Aboveground Biomass Dynamics
## ============================================================
##
## Purpose: Estimate aboveground biomass (AGB) stock, productivity,
## mortality, and net change for the BCI 50-ha forest plot across
## 9 censuses (1982–2020).
##
## Key intermediate datasets:
##   df_stem_status  – standing AGB per alive stem per census
##   df_stem_prod    – AGB productivity per stem per census interval
##   df_stem_mort    – AGB mortality flux per stem at last alive census
##   df_stem_demo    – aggregated demographics per quadrat × size × interval
##   df_stem_demo_quadrat – same, all size classes pooled
##
## Indexing convention:
##   CensusID = c  is the INITIAL census of the interval  c → c+1.
##   Dagb  is lag-differenced at c+1 but re-indexed to c in df_stem_demo
##         so that Dagb, DagbM, and stock all share the same starting census.
##   DagbM is assigned at c (the last alive census before the stem dies).
##
## All AGB calculations use taper-corrected DBH (dbh_t, _t suffix).
## [EDGE CASE] notes flag non-trivial boundary conditions in the code.
## ============================================================

rm(list = ls())

## --- 0. Parameters -----------------------------------------------------------

remove_strangler_figs <- TRUE # exclude all Ficus strangler spp. (Rutishauser 2020)
use_median_palm_dbh <- TRUE # replace palm DBH with species median (except Socratea)
biomass_allometry <- "chave14" # "chave14" or "chave05"
use_local_height_allometry <- TRUE # use Martinez-Cano 2019 height model
use_kohyama19_correction <- TRUE # Kohyama et al. 2019 productivity correction

lower_ddbh_threshold <- -2.5 # cm yr⁻¹; minimum plausible DBH growth
upper_ddbh_threshold <- 7.5 # cm yr⁻¹; maximum plausible DBH growth
hom_change_threshold <- 1.0 # m;       flag if POM shifted by more than this

# DBH interpolation method for stems with Rstatus == "A" but dbh = NA.
# Change this and re-run for sensitivity analysis.
# One of: "linear" (default, interpolates between flanking measurements),
#         "locf"   (last observation carried forward, conservative),
#         "mean"   (mean of nearest flanking measured values)
dbh_interp_method <- "linear"

library(data.table)
library(truncnorm)
library(lubridate)

## --- 1. Load raw census data -------------------------------------------------

# Load all 9 census Rdata files and stack into one long table.
# idcol = "censusID" records which file each row came from.
bci_stem_nums <- as.character(1:9)

census_list <- lapply(bci_stem_nums, function(num) {
  filepath <- paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", num, ".Rdata")
  if (!file.exists(filepath)) stop("Missing census file: ", filepath)
  load(filepath)
  get(paste0("bci.stem", num))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)

df_stem <- data.table::rbindlist(census_list, fill = TRUE, idcol = "censusID")
rm(census_list, bci_stem_nums)

# CensusID as integer — essential for ordering and arithmetic (e.g. CensusID - 1L)
df_stem[, CensusID := as.integer(CensusID)]

## --- 2. Merge species taxonomy and wood density data ------------------------

# ## Download wood density data from Wright & Mulle-Landau 2026 Dryad
# url: "https://datadryad.org/dataset/doi:10.5061/dryad.5qfttdzn3"
bci_wd <- fread("./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/wd/doi_10_5061_dryad_5qfttdzn3__v20260403/WD_species.txt")
bci_wd <- bci_wd[, .(sp = tolower(sp6), wsg = wd100.mean)][!is.na(wsg)]

bci.spptable <- fread("./BCI_stem_reconstruction/DATA/RTABLES/bci.spptable.csv")
bci.spptable <- unique(bci.spptable[, .(Family, sp, Genus, Species = SpeciesName, Latin)])

df_stem <- merge(df_stem, bci.spptable, by = "sp", all.x = TRUE)
df_stem <- merge(df_stem, bci_wd, by = "sp", all.x = TRUE)

rm(bci_wd, bci.spptable)

# Fill missing wood specific gravity (WSG) hierarchically:
#   1st: species-level mean within the same genus
#   2nd: genus-level mean within the same family
#   3rd: global mean across all species
# [EDGE CASE] Species whose entire family has no WSG data will receive the
#             global mean (~0.58 g/cm³). Flag these if precision matters.
df_stem[, wsg := ifelse(is.na(wsg), ave(wsg, Genus, FUN = function(x) mean(x, na.rm = TRUE)), wsg)]
df_stem[, wsg := ifelse(is.na(wsg), ave(wsg, Family, FUN = function(x) mean(x, na.rm = TRUE)), wsg)]
df_stem[, wsg := ifelse(is.na(wsg), mean(wsg, na.rm = TRUE), wsg)]

# Replace NA taxonomy labels so downstream group-by operations produce clean groups
df_stem[, Family := ifelse(is.na(Family), "Unknown", Family)]
df_stem[, Genus := ifelse(is.na(Genus), "Unknown", Genus)]
df_stem[, Species := ifelse(is.na(Species), "Unknown", Species)]
df_stem[, Latin := ifelse(is.na(Latin), "Unknown", Latin)]

## --- 3. Strangler fig removal ------------------------------------------------
# Giant strangler Ficus (> 500 mm DBH = > 50 cm; raw BCI units before /10 conversion)
# are excluded following Rutishauser et al. 2020. They violate the standard
# allometric assumptions and inflate plot-level AGB disproportionately.
# Exclusion applies to ALL censuses of affected trees/stems.
large_strangler_figs <- df_stem[
  dbh > 500 &
    Genus == "Ficus" &
    Species %in% c("costaricana", "obtusifolia", "popenoei", "trigonata")
]

if (remove_strangler_figs) {
  df_stem <- df_stem[
    !treeID %in% large_strangler_figs$treeID &
      !stemID %in% large_strangler_figs$stemID
  ]
  message(sprintf("[STRANGLER] Removed %d trees / %d stems.",
    uniqueN(large_strangler_figs$treeID), uniqueN(large_strangler_figs$stemID)))
}
rm(large_strangler_figs)

## --- 4. DBH unit conversion: mm → cm ----------------------------------------
# Raw BCI RTABLE DBH is in mm. All allometric equations expect cm.
df_stem[, dbh := dbh / 10]

## --- 5. Palm DBH correction --------------------------------------------------
# Non-Socratea palms do not grow in diameter; observed DBH changes are measurement
# error. Replace each species' DBH with the species median across all censuses
# (Rutishauser et al. 2020).
if (use_median_palm_dbh) {
  df_stem[Family == "Arecaceae" & Genus != "Socratea", dbh := median(dbh, na.rm = TRUE), .(Latin)]
}

## --- 6. Taper correction -----------------------------------------------------
# DBH is corrected for taper using Cushman et al. 2014.
# Taper adjusts DBH to what it would be at 1.3 m when measured higher (e.g., above
# buttresses). The corrected value is stored as `dbh_t`; all downstream AGB uses _t.
#
# Applied universally: although ideal only for buttressed species, the equation
# returns dbh_t ≈ dbh when hom ≈ 1.3 m (b ≈ 0), so non-buttressed stems are
# unaffected.
#
# [EDGE CASE] Stems with hom = NA are imputed to 1.3 m (b = 0, dbh_t = dbh).
# [EDGE CASE] Interpolation of dbh_t (section 7) is done AFTER taper correction:
#             taper uses the raw measured dbh; filling missing dbh_t comes after.
df_stem[, hom := ifelse(is.na(hom), 1.3, hom)]
df_stem[, b := exp(-2.0205 - 0.5053 * log(dbh) + 0.3748 * log(hom))]
df_stem[!is.na(hom), dbh_t := dbh * exp(b * (hom - 1.3))]
df_stem[, b := NULL] # intermediate taper coefficient — no longer needed

## --- 7. Interpolate DBH for unmeasured-but-alive stems ---------------------

# Some stems have Rstatus == "A" (alive) but dbh == NA — they were assumed alive
# in the field and not re-measured (common in BCI between full re-censuses).
# These rows would otherwise produce NA growth / NA AGB and be silently dropped.
#
# We provide three interpolation methods. Set `dbh_interp_method` at the top of
# the script and re-run to assess sensitivity:
#   "linear" — linear interpolation between flanking MEASURED DBHs (preferred)
#   "locf"   — last observation carried forward (assumes no growth; conservative)
#   "mean"   — mean of nearest flanking measured DBHs
#
# Only rows where  Rstatus == "A"  AND  dbh is NA  are filled. Original raw
# DBH is preserved in `dbh_raw`; an integer flag `was_interpolated` marks
# the affected rows so downstream code (and tests) can audit them.
#
# [EDGE CASE] Stems with no measured DBH on at least one side cannot be
#             interpolated by "linear" or "mean" — they remain NA.
# [EDGE CASE] "locf" can fill a trailing gap (after the last measurement) only
#             when Rstatus == "A" persists; it cannot extrapolate growth.
# [EDGE CASE] If a stem has Rstatus = A throughout but never measured, no method
#             can fill it; it is excluded from analyses.

interpolate_dbh <- function(dt, method = c("linear", "locf", "mean"), var_to_interpolate = "dbh") {
  method <- match.arg(method)
  # Save raw values before interpolation
  if (!"dbh_raw" %in% names(dt)) {
    dt[, dbh_raw := get(var_to_interpolate)]
  }
  fill_one <- function(x, status) {
    target <- !is.na(status) & status == "A" & is.na(x)
    if (!any(target)) {
      return(x)
    }
    measured <- !is.na(x)
    if (sum(measured) < 1L) {
      return(x)
    }
    idx_meas <- which(measured)
    out <- x
    if (method == "linear" && sum(measured) >= 2L) {
      yhat <- approx(
        x = idx_meas, y = x[idx_meas], xout = seq_along(x),
        method = "linear", rule = 1
      )$y
      out[target] <- yhat[target]
    } else if (method == "mean") {
      for (i in which(target)) {
        prev_i <- if (any(idx_meas < i)) max(idx_meas[idx_meas < i]) else NA_integer_
        next_i <- if (any(idx_meas > i)) min(idx_meas[idx_meas > i]) else NA_integer_
        vals <- c(if (!is.na(prev_i)) x[prev_i], if (!is.na(next_i)) x[next_i])
        if (length(vals)) out[i] <- mean(vals)
      }
    } else { # "locf"
      last_val <- NA_real_
      for (i in seq_along(x)) {
        if (!is.na(x[i])) last_val <- x[i]
        if (target[i] && !is.na(last_val)) out[i] <- last_val
      }
    }
    out
  }
  # Apply interpolation using get() to reference column by name
  dt[, (var_to_interpolate) := fill_one(get(var_to_interpolate), Rstatus), by = .(treeID, stemID)]
  invisible(dt)
}

# Count alive stems that are missing taper-corrected DBH before interpolation.
# dbh_t is NA whenever dbh (raw) is NA, so the counts are equivalent, but we
# count dbh_t explicitly for clarity since that is the column being filled.
n_to_fill <- df_stem[!is.na(Rstatus) & Rstatus == "A" & is.na(dbh_t), .N]
message(sprintf(
  "[INTERP] %d stem-census rows are Rstatus=A & dbh_t=NA. Filling with method='%s'.",
  n_to_fill, dbh_interp_method
))
interpolate_dbh(df_stem, method = dbh_interp_method, var_to_interpolate = "dbh_t")
df_stem[, was_interpolated := is.na(dbh_raw) & !is.na(dbh_t)]
n_filled <- df_stem[was_interpolated == TRUE, .N]
message(sprintf("[INTERP] Filled %d / %d candidate rows.", n_filled, n_to_fill))
rm(n_to_fill, n_filled)

# Diagnostic: stems that are alive but still have NA dbh_t after interpolation.
# These cannot be filled (e.g. stem never measured on either side of the gap).
inc <- unique(df_stem[!is.na(Rstatus) & Rstatus == "A" & is.na(dbh_t)]$stemID)
if (length(inc) > 0L) {
  message(sprintf("[INTERP] %d stems remain with NA dbh_t after interpolation.", length(inc)))
}

## --- 8. Allometric functions -------------------------------------------------

# Two supported allometries: Chave et al. 2014 and Chave et al. 2005.
# Height is estimated with the multi-species model of Martinez-Cano et al. 2019
# (h = 58.0 * dbh^0.73 / (21.8 + dbh^0.73); asymptote ~58 m at BCI).
# Palms use the Goodman et al. 2013 palm-specific allometry.
# Results in Mg of dry biomass per stem.
agb_bci <- function(dbh, # dbh, in cm
                    wsg, # wood specific gravity values
                    palms = NULL, # logical: is the individual a palm?
                    method = "chave14",
                    use_height_allom = TRUE) {
  if (use_height_allom) {
    # with multi-species height allometry from Martinez Cano et al 2019
    h <- 58.0 * dbh^0.73 / (21.8 + dbh^0.73)
    if (method == "chave05") {
      # Chave et al 2005 - moist forests, with height (in Mg)
      agb <- 0.0509 * wsg * dbh^2 * h / 1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 4 with the BIOMASS package
      agb <- (0.0673 * (wsg * h * dbh^2)^0.976) / 1000
    }
  } else {
    # without any height information
    if (method == "chave05") {
      # Chave et al 2005 - moist forests, without height (in Mg)
      agb <- wsg * exp(-1.499 + 2.148 * log(dbh) +
        0.207 * log(dbh)^2 - 0.0281 * log(dbh)^3) / 1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 7 with the BIOMASS package (transform into kg)
      E <- 0.05176398
      agb <- exp(-2.023977 - 0.89563505 * E + 0.92023559 *
        log(wsg) + 2.79495823 * log(dbh) - 0.04606298 * (log(dbh)^2)) / 1000
    }
  }
  if (!is.null(palms)) {
    # palm specific allometry from Goodman et al. 2013
    agb[palms] <- 0.0417565 * dbh[palms]^2.7483 / 1000
  }
  return(agb)
}

## --- 9. Estimate AGB (taper-corrected) --------------------------------------

# AGB (Mg dry mass) from taper-corrected DBH (`dbh_t`) only.
# Allometry: Chave et al. 2014 (eq. 4) + Martinez-Cano et al. 2019 height model.
# Palms use the Goodman et al. 2013 palm-specific allometry.
df_stem[, agb_t := agb_bci(
  dbh = dbh_t,
  wsg = wsg,
  method = "chave14",
  use_height_allom = TRUE,
  palms = (Family == "Arecaceae")
)]

## --- Correction for 1985 DBH rounding bias ---------------------------------
#
# In the 1985 census (CensusID 2), stems with DBH < 5.5 cm were recorded in
# 5-mm intervals (rounded down to the nearest 5 mm). This compresses AGB at
# the START of the 1985→1990 interval and biases growth estimates upward.
#
# Correction:
#   1. Identify stems < 5.5 cm dbh_t in CensusID 2 (1985).
#   2. Assign each a 5-mm rounding class from their 1985 dbh_t:
#      class = floor(dbh_t / 0.5) * 0.5  (reproduces the field rounding).
#   3. Compute mean agb_t per rounding class from CensusID 3 (1990).
#      Using 1990 values as reference avoids carrying the rounding error forward.
#   4. Replace agb_t at BOTH CensusID 2 (1985) AND CensusID 3 (1990) with
#      the class mean. Setting both endpoints to the same value makes the
#      1985→1990 growth contribution from these stems ≈ 0, which is the
#      conservative but unbiased choice when the initial measurement is unreliable.
#
# [EDGE CASE] Stems present in 1985 but absent/unmeasured in 1990 do not
#             contribute to mean_agb; their CensusID 3 agb_t is left unchanged.
# [EDGE CASE] A rounding class with no CensusID 3 observations gets agb_t_m = NA
#             and the substitution is skipped for those stems (warning issued).

# Step 1: stems < 5.5 cm dbh_t in the 1985 census (CensusID 2)
small_stems_1985 <- df_stem[
  CensusID == 2L & !is.na(dbh_t) & dbh_t < 5.5,
  unique(stemID)
]
message(sprintf(
  "[ROUNDING] %d stems identified with dbh_t < 5.5 cm in CensusID 2 (1985).",
  length(small_stems_1985)
))

# Step 2: 5-mm rounding class from the 1985 (CensusID 2) taper-corrected DBH.
# One row per stemID — used as the join key for steps 3 and 4.
dbh_r_lut <- df_stem[
  stemID %in% small_stems_1985 & CensusID == 2L,
  .(stemID, dbh_r = floor(dbh_t / 0.5) * 0.5)
]

# Step 3: mean agb_t per rounding class from CensusID 3 (1990).
# Join the 1985 rounding class to the 1990 rows of the same stems, then
# aggregate. Only stems with a valid 1990 agb_t contribute to the mean.
mean_agb <- merge(
  df_stem[
    stemID %in% small_stems_1985 & CensusID == 3L & !is.na(agb_t),
    .(stemID, agb_t)
  ],
  dbh_r_lut,
  by = "stemID"
)[, .(agb_t_m = mean(agb_t, na.rm = TRUE)), by = dbh_r]

n_missing_class <- dbh_r_lut[!dbh_r %in% mean_agb$dbh_r, .N]
if (n_missing_class > 0L) {
  warning(sprintf("[ROUNDING] %d rounding classes have no CensusID 3 data; substitution skipped for those stems.", n_missing_class))
}

# Step 4: join the per-stem mean AGB and apply substitution at CensusID 2 AND 3.
# The join is on stemID only, so agb_t_m propagates to all census rows of the
# affected stems — but the assignment is restricted to the two target censuses,
# so all other censuses are untouched.
subst_lut <- merge(dbh_r_lut, mean_agb, by = "dbh_r")[, .(stemID, agb_t_m)]
df_stem <- merge(df_stem, subst_lut, by = "stemID", all.x = TRUE)
df_stem[!is.na(agb_t_m) & CensusID %in% c(2L, 3L), agb_t := agb_t_m]
df_stem[, agb_t_m := NULL]
rm(small_stems_1985, dbh_r_lut, mean_agb, subst_lut, n_missing_class)

## --- 10. Sort rows and assign DBH size classes -----------------------------

# Sort by stem and census — required for all lag/lead operations below.
data.table::setorder(df_stem, treeID, stemID, CensusID)

# One consistent size class definition used for outlier substitution AND aggregation.
# Breaks (cm): [0,10), [10,20), [20,50), [50,500)
# [EDGE CASE] Stems with NA DBH get size = NA and are excluded from
#             size-stratified analyses.
df_stem[, size := cut(dbh_t,
  breaks = c(0, 10, 20, 50, 500),
  include.lowest = TRUE, right = FALSE
)]

## --- 11. Annualised growth rates (recruits explicitly handled) -------------

# Convention: Ddbh and Dagb represent the annualised change from census c to
# c+1 of the SAME stem. The value is stored on the row of census c+1 (the END
# of the interval). When df_stem_demo is built, we re-index by CensusID - 1L
# so growth, mortality, and ntrees all refer to the INITIAL census c.
#
# Two sources contribute biomass GAIN:
#   (a) Growth of stems present in BOTH c and c+1
#       → Dagb = (agb[c+1] - agb[c]) / dT  (the standard lag-difference)
#   (b) Recruits: stems first observed alive at c+1 with no prior alive obs
#       → Dagb = agb[c+1] / dT_plot  where dT_plot is the quadrat census interval
#
# Recruit gain is allocated to interval (c → c+1) so it accumulates with growth
# at the same row. The `is_recruit` flag preserves provenance for tests.
#
# [EDGE CASE] A stem's first observation may legitimately be at census 1.
#             This is NOT a recruit — it is the start of monitoring. We
#             distinguish recruits as "first appearance at census >= 2".
# [EDGE CASE] dT ≤ 0 indicates a date recording error (duplicated census or
#             reversed dates). A warning is issued.
# [EDGE CASE] dT >= 10 yr is also flagged as suspicious (BCI intervals are ~5y).

# 11a. Time interval per stem between consecutive observations
data.table::setorder(df_stem, treeID, stemID, CensusID)

# Some rows lack ExactDate (stems from unidentified quadrats). We fill them with
# a two-step imputation so that every row gets a date and dT is never NA due to
# a missing date (only the first census of each stem legitimately has dT = NA
# because there is no prior row to difference against).
#
# Step 1: fill with the median date of all stems in the same quadrat × census.
# Step 2: fill any remaining NAs (e.g. quadrat itself is NA/unknown) with the
#         median date across the entire census (plot-wide).
#
# [EDGE CASE] If an entire census has no dated stems, date_plot_census is NA and
#             ExactDate remains NA for those rows. This is extremely unlikely with
#             BCI data; a warning is issued below if it occurs.
n_na_before <- df_stem[is.na(ExactDate), .N]
message(sprintf("[DATES] %d rows have NA ExactDate before imputation.", n_na_before))

df_stem[, date_quad_census := median(ExactDate, na.rm = TRUE), .(quadrat, CensusID)]
df_stem[, date_plot_census := median(ExactDate, na.rm = TRUE), .(CensusID)]

df_stem[, ExactDate := fifelse(is.na(ExactDate), date_quad_census, ExactDate)]
df_stem[is.na(ExactDate), ExactDate := date_plot_census]

df_stem[, `:=`(date_quad_census = NULL, date_plot_census = NULL)]

# Verify: how many ExactDate NAs remain?
n_na_after <- df_stem[is.na(ExactDate), .N]
if (n_na_after > 0L) {
  warning(sprintf(
    "[DATES] %d rows STILL have NA ExactDate after imputation (entire census undated?).",
    n_na_after
  ))
} else {
  message("[DATES] All ExactDate values filled — dT will be NA only for each stem's first census row (correct).")
}
rm(n_na_before, n_na_after)

df_stem[, prev_ExactDate := shift(ExactDate, type = "lag"), .(treeID, stemID)]
df_stem[, dT := as.numeric(difftime(ExactDate, prev_ExactDate, units = "days")) / 365.25]

n_bad_dT <- df_stem[!is.na(dT) & (dT <= 0 | dT >= 10), .N]
if (n_bad_dT > 0) {
  warning(sprintf("[EDGE CASE] %d stem-census rows have dT <= 0 or >= 10 yr. Inspect ExactDate.", n_bad_dT))
}
rm(n_bad_dT)

# 11b. Ensure rows are in order before lag/lead operations.
data.table::setorder(df_stem, treeID, stemID, CensusID)

# 11c. Lag base values per stem
df_stem[, prev_dbh_t := shift(dbh_t), .(treeID, stemID)]
df_stem[, prev_agb_t := shift(agb_t), .(treeID, stemID)]

# 11d. Recruit detection: first row per stem where the stem is alive (Rstatus=A,
#      valid dbh) AND there is no PRIOR alive observation. Recruits at CensusID 1
#      are excluded (they are start-of-monitoring, not new recruits).
df_stem[, n_alive_prior := cumsum(!is.na(dbh_t) & Rstatus == "A") -
  (!is.na(dbh_t) & Rstatus == "A"), .(treeID, stemID)]
df_stem[, is_recruit := !is.na(dbh_t) & Rstatus == "A" &
  n_alive_prior == 0L & CensusID >= 2L]

# 11e. Standard lag-difference growth (ongoing stems with both c and c+1 alive+measured)
df_stem[, Ddbh_t := fifelse(!is.na(dbh_t) & !is.na(prev_dbh_t), (dbh_t - prev_dbh_t) / dT, NA_real_)]
df_stem[, Dagb_t := fifelse(!is.na(agb_t) & !is.na(prev_agb_t), (agb_t - prev_agb_t) / dT, NA_real_)]

# 11f. Recruit gain: assign Dagb = agb / dT at the row where the recruit
#      first appears.
df_stem[is_recruit == TRUE, Dagb_t := fifelse(!is.na(dT) & dT > 0, agb_t / dT, NA_real_)]

# Remove temporary lag and helper columns
df_stem[, c(
  "prev_ExactDate", "n_alive_prior",
  "prev_dbh_t", "prev_agb_t"
) := NULL]

## --- 12. HOM change detection -----------------------------------------------
# Flag observations where the height of measurement (HOM/POM) shifted by more
# than `hom_change_threshold` (default: 1 m) between consecutive censuses.
# Such shifts produce unreliable DBH growth estimates and are treated as outliers
# regardless of the computed Ddbh value.
# [EDGE CASE] First observation per stem has prev_hom = NA → hom_change = FALSE.

data.table::setorder(df_stem, treeID, stemID, CensusID)
df_stem[, prev_hom := shift(hom), .(treeID, stemID)]
df_stem[, hom_change := !is.na(prev_hom) & abs(hom - prev_hom) > hom_change_threshold]
df_stem[, prev_hom := NULL]

## --- 13. Outlier detection and growth substitution -------------------------

# A growth interval is flagged as an outlier when:
#   (a) DBH growth exceeds the plausible range [lower_ddbh_threshold, upper_ddbh_threshold], OR
#   (b) the height of measurement changed between censuses (HOM shift).
#
# Outlier AGB/BA changes are substituted with the size-class relative rate:
#   Dagb_s = (sum non-outlier Dagb / sum non-outlier AGB) * own AGB
# This preserves the direction of the plot-level signal while removing extremes.
#
# [EDGE CASE] A size class with ALL stems as outliers will have tot_rawp = NaN
#             (0/0). Outlier substitution for those stems produces NaN Dagb_s.
# [EDGE CASE] Rows with NA Ddbh_t (first obs, both censuses missing) are NOT
#             flagged as outliers; they simply have no growth estimate.

df_stem[, outlier := !is.na(Ddbh_t) &
  (hom_change | Ddbh_t > upper_ddbh_threshold | Ddbh_t < lower_ddbh_threshold)]

# Size-class relative growth rate: used to substitute growth for outlier stems.
# Recruits are excluded from the rate denominator — their Dagb is biomass "birth",
# not growth, and would inflate the rate.
# lag_agb_t is the AGB at census c for each stem (the denominator of relative growth).
df_stem[, lag_agb_t := shift(agb_t), .(treeID, stemID)]
df_stem[, lag_agb_t := ifelse(is.na(lag_agb_t) & !is.na(Dagb_t) & is_recruit == FALSE, agb_t, lag_agb_t)]
df_stem[, tot_rawp_t := sum(Dagb_t[!outlier & is_recruit == FALSE & !is.na(Dagb_t)], na.rm = TRUE) /
  sum(lag_agb_t[!outlier & is_recruit == FALSE & !is.na(Dagb_t)], na.rm = TRUE), by = size]
df_stem[, lag_agb_t := NULL]

# Apply substitution: outlier rows get the size-class mean relative rate × their AGB.
# Non-outlier rows keep their measured value.
df_stem[, Dagb_t_s := ifelse(outlier, tot_rawp_t * agb_t, Dagb_t)]

# Remove substitution rate temporaries
df_stem[, c("tot_rawp_t") := NULL]

## --- 14. Mortality flux assignment ------------------------------------------

# A stem contributes a MORTALITY FLUX in the interval c → c+1 when:
#   (1) at census c the stem is ALIVE with a valid DBH (Rstatus == "A"), AND
#   (2) at census c+1 the stem is NOT alive (Rstatus is "D", "G", or absent).
#
# `last_census_alive` = last census where Rstatus == "A" and dbh is valid.
#  At that row, if a next census exists, mortality flux = agb / dT_mort
#  where dT_mort is the time from census c to census c+1.
#
# Rstatus codes used by the BCI RTABLES (per 3_PREPARE_R_TABLES/2_create_R_tables_BCI.R):
#   A = alive,  D = dead (whole stem),  G = stem dead but tree alive,
#   P = prior to first observation,  N = unresolved (rare; resolved to P/D/G)
#
# [EDGE CASE] Stems alive through the final census (9) have no forward
#             interval → dT_mort = NA → DagbM = NA. Correctly NOT counted as dead.
# [EDGE CASE] "Zombie" stems (A → D → A) have last_census_alive set to the
#             LAST alive census. Intermediate dead intervals are not tracked
#             as separate events here.
# [EDGE CASE] A stem alive in only one census is counted as a mortality
#             if its next census is D or G — ecologically correct.
# [EDGE CASE] Stem-level next_date can be NA if the stem has no row at c+1.
#             We fall back to plot-level dT_plot of CensusID c+1 in that case.

# Identify each stem's last census with a valid taper-corrected alive DBH.
# Using dbh_t for consistency: dbh_t is the corrected and interpolated
# value that feeds all downstream AGB calculations.
df_stem[
  !is.na(dbh_t) & Rstatus == "A",
  last_census_alive := max(CensusID),
  .(treeID, stemID)
]

# Status of the NEXT census (per stem)
df_stem[, next_Rstatus := shift(Rstatus, type = "lead"), .(treeID, stemID)]
df_stem[, next_date := shift(ExactDate, type = "lead"), .(treeID, stemID)]
df_stem[, dT_mort := as.numeric(difftime(next_date, ExactDate, units = "days")) / 365.25]

data.table::setorder(df_stem, treeID, stemID, CensusID)

# A stem is dead in the next census if next_Rstatus is D or G, OR if the stem
# row simply disappears after its last alive census (next_Rstatus == NA but
# last_census_alive < max plot census). The latter is rare in well-curated data.
max_census <- df_stem[, max(CensusID, na.rm = TRUE)]
df_stem[, dies_next := CensusID == last_census_alive &
  last_census_alive < max_census &
  (next_Rstatus %in% c("D", "G") | is.na(next_Rstatus))]
rm(max_census)

# Apply mortality flux at the last alive census whose next census is dead.
df_stem[
  dies_next == TRUE & !is.na(dT_mort) & dT_mort > 0,
  `:=`(
    DagbM_t = agb_t / dT_mort
  )
]

## ============================================================
## Dataset 1: Current Status
## ============================================================
##
## One row per ALIVE stem per census.
## Represents standing AGB at each measured census (the "stock").
## Use this to compute AGB at any given census snapshot.
##
# Filter on dbh_t (taper-corrected) rather than raw dbh so that only rows with
# a valid corrected measurement are included as "standing stock".
df_stem_status <- df_stem[
  !is.na(dbh_t) & Rstatus == "A",
  .(
    stemID, treeID, sp, Latin, Family, Genus, Species,
    quadrat, CensusID, dT, ExactDate,
    dbh_t, hom, agb_t, size, wsg
  )
]

## ============================================================
## Dataset 2: Productivity (growth + recruits)
## ============================================================
##
## One row per stem per census interval where AGB GAIN can be computed.
## CensusID here is c+1 (the END of the interval c → c+1).
## `is_recruit == TRUE` rows are first appearances; their Dagb = agb / dT_plot.
## `is_recruit == FALSE` rows are ongoing-stem growth (lag-difference).
##
## All growth values are taper-corrected (`_t` suffix):
##   Dagb_t   — raw lag-difference growth (before outlier substitution)
##   Dagb_t_s — outlier-substituted growth  ← USE THIS for analyses
##   Ddbh_t   — annualised taper-corrected DBH growth (cm yr⁻¹; diagnostics)
##
df_stem_prod <- df_stem[
  !is.na(Dagb_t_s) & !is.na(dT) & dT > 0,
  .(
    stemID, treeID, sp, quadrat, CensusID,
    dT, size, agb_t,
    Ddbh_t, Dagb_t, Dagb_t_s,
    outlier, hom_change, is_recruit, was_interpolated
  )
]

## ============================================================
## Dataset 3: Mortality
## ============================================================
##
## One row per stem at its LAST alive census (census c).
## DagbM = AGB / dT_mort = annualised biomass loss for the c → c+1 interval.
## CensusID here is c (the census at the START of the mortality interval).
##
## [EDGE CASE] Stems alive through census 9 are NOT in this dataset
##             (DagbM = NA because dT_mort = NA — correctly excluded).
##
df_stem_mort <- df_stem[
  !is.na(DagbM_t),
  .(
    stemID, treeID, sp, quadrat, CensusID,
    last_census_alive, dT_mort, size,
    agb_t, DagbM_t
  )
]

## ============================================================
## Dataset 4a: df_stem_demo  (per quadrat × size class × interval)
## Dataset 4b: df_stem_demo_quadrat (per quadrat × interval, all sizes pooled)
## ============================================================
##
## Row definition: quadrat q [× size class s], INITIAL census c of interval c→c+1.
## All components are indexed at c:
##   ntrees  — alive stems in census c  (stock)
##   agb     — total standing AGB at census c  (Mg)
##   dT      — mean interval length in years
##   Dagb    — total AGB productivity for the interval (Mg yr⁻¹; growth + recruits)
##   Dagb_growth  — productivity from growth of standing stems only
##   Dagb_recruit — productivity from new (recruit) stems only
##   DagbM   — total AGB mortality loss for the interval (Mg yr⁻¹)
##
## Productivity is stored at c+1 in df_stem_prod → subtract 1 to re-index to c.
## Mortality is already stored at c in df_stem_mort.
##
## [EDGE CASE] Census 9 has no forward interval → no Dagb / dT; rows dropped.
## [EDGE CASE] Stems with NA or empty quadrat are excluded from spatial aggregation.
## [EDGE CASE] size = NA (NA DBH) are excluded from df_stem_demo (size-stratified)
##             but ARE included in df_stem_demo_quadrat (pooled).

build_demo <- function(group_cols) {
  # Aggregate productivity at the INITIAL census (re-index by CensusID - 1).
  # Split into growth vs recruit components for diagnostics.
  prod_agg <- df_stem_prod[
    !is.na(quadrat) & quadrat != "",
    .(
      Dagb_growth = sum(Dagb_t_s[is_recruit == FALSE], na.rm = TRUE),
      Dagb_recruit = sum(Dagb_t_s[is_recruit == TRUE], na.rm = TRUE),
      Dagb = sum(Dagb_t_s, na.rm = TRUE),
      dT = mean(dT, na.rm = TRUE)
    ),
    by = group_cols
  ]
  prod_agg[, CensusID := CensusID - 1L] # re-index from c+1 to c

  mort_agg <- df_stem_mort[
    !is.na(quadrat) & quadrat != "",
    .(
      DagbM = sum(DagbM_t, na.rm = TRUE)
    ),
    by = group_cols
  ]

  stock_agg <- df_stem_status[
    !is.na(quadrat) & quadrat != "",
    .(
      ntrees = uniqueN(treeID[!is.na(dbh_t)]),
      agb_t = sum(agb_t, na.rm = TRUE)
    ),
    by = group_cols
  ]

  out <- merge(stock_agg, prod_agg, by = group_cols, all = TRUE)
  out <- merge(out, mort_agg, by = group_cols, all = TRUE)
  out[is.na(ntrees), ntrees := 0L]
  out[is.na(agb_t), agb_t := 0]
  out[is.na(Dagb_growth), Dagb_growth := 0]
  out[is.na(Dagb_recruit), Dagb_recruit := 0]
  out[is.na(Dagb), Dagb := 0]
  out[is.na(DagbM), DagbM := 0]
  # Drop final-census rows (no forward interval defined)
  max_c <- df_stem[, max(CensusID, na.rm = TRUE)]
  out <- out[CensusID < max_c]
  out[]
}

# 4a. Per quadrat × size class
df_stem_demo <- build_demo(c("quadrat", "size", "CensusID"))
data.table::setorder(df_stem_demo, quadrat, size, CensusID)

# 4b. Per quadrat (all size classes pooled)
df_stem_demo_quadrat <- build_demo(c("quadrat", "CensusID"))
data.table::setorder(df_stem_demo_quadrat, quadrat, CensusID)

## ============================================================
## Kohyama et al. 2019 correction for unrecorded growth
## ============================================================
##
## Stems that grow and die *between* two censuses contribute biomass gain
## that is never recorded. Kohyama et al. 2019
## derive bias-corrected productivity (G*) and mortality (M*) estimators.
##
## Notation (per quadrat × census-interval row):
##   B0  = agb_t  — stand AGB at the START of the interval (Mg per quadrat)
##   BT  = B0 + (Dagb − DagbM) × dT  — estimated AGB at END of interval
##   BS0 = B0 − DagbM × dT            — AGB of SURVIVING stems at the start
##
## Corrected fluxes (eqs. 7–8 in Kohyama et al. 2019):
##   G* = log(BT / BS0) × (BT − B0) / (dT × log(BT / B0))
##   M* = log(B0 / BS0) × (BT − B0) / (dT × log(BT / B0))
##
## Applied at the AGGREGATED quadrat level (NOT at stem level) to both
## df_stem_demo_quadrat and df_stem_demo.  Corrected values are stored in
## Dagb_k and DagbM_k.  Where the formula is undefined (BT == B0, or
## BS0 / BT ≤ 0) we fall back to the uncorrected value and log a count.
## When use_kohyama19_correction = FALSE the _k columns are plain aliases
## of Dagb / DagbM so all downstream code can always reference _k columns.

kohyama_correction <- function(stock, gain, loss, dT, output = "prod") {
  B0  <- stock
  BS0 <- B0 - loss * dT            # surviving AGB at interval start
  BT  <- B0 + (gain - loss) * dT   # estimated AGB at interval end
  denom <- dT * log(BT / B0)
  if (output == "prod") {
    log(BT / BS0) * (BT - B0) / denom
  } else if (output == "mort") {
    log(B0 / BS0) * (BT - B0) / denom
  } else {
    stop("'output' must be \"prod\" or \"mort\".")
  }
}

# Internal helper: applies the correction in-place; counts and replaces any
# undefined results (NaN / Inf) with the uncorrected value.
apply_kohyama <- function(dt, tag = "") {
  dt[, Dagb_k  := kohyama_correction(agb_t, Dagb, DagbM, dT, output = "prod")]
  dt[, DagbM_k := kohyama_correction(agb_t, Dagb, DagbM, dT, output = "mort")]
  n_bad_prod <- dt[!is.finite(Dagb_k)  | is.na(Dagb_k),  .N]
  n_bad_mort <- dt[!is.finite(DagbM_k) | is.na(DagbM_k), .N]
  dt[!is.finite(Dagb_k)  | is.na(Dagb_k),  Dagb_k  := Dagb]
  dt[!is.finite(DagbM_k) | is.na(DagbM_k), DagbM_k := DagbM]
  if (n_bad_prod + n_bad_mort > 0L) {
    message(sprintf(
      "[KOHYAMA%s] %d prod / %d mort rows fell back to uncorrected values (BT == B0 or BS0 <= 0).",
      tag, n_bad_prod, n_bad_mort
    ))
  }
  invisible(dt)
}

if (use_kohyama19_correction) {
  apply_kohyama(df_stem_demo_quadrat, tag = " (quadrat)")
  apply_kohyama(df_stem_demo,         tag = " (quadrat×size)")
  message("[KOHYAMA] Corrected fluxes stored in Dagb_k / DagbM_k.")
} else {
  # Alias _k columns to uncorrected values so downstream code is uniform.
  df_stem_demo_quadrat[, `:=`(Dagb_k = Dagb, DagbM_k = DagbM)]
  df_stem_demo[,         `:=`(Dagb_k = Dagb, DagbM_k = DagbM)]
  message("[KOHYAMA] Skipped (use_kohyama19_correction = FALSE). _k columns alias uncorrected values.")
}

## ============================================================
## TESTS — sanity checks for recruits, mortality, interpolation, conservation
## ============================================================
## Each test prints PASS / FAIL with a short message. Failures do NOT halt
## execution — they are diagnostic. Run them after every change to the pipeline.

run_tests <- function() {
  cat("\n========== PIPELINE SANITY TESTS ==========\n")
  pass <- function(msg) cat(sprintf("  [PASS] %s\n", msg))
  fail <- function(msg) cat(sprintf("  [FAIL] %s\n", msg))

  # T1. Recruits ARE included in productivity (not silently dropped)
  n_rec <- df_stem_prod[is_recruit == TRUE, .N]
  if (n_rec > 0) {
    pass(sprintf("T1 ingrowth/recruit present in df_stem_prod: %d rows", n_rec))
  } else {
    fail("T1 NO ingrowth/recruit in df_stem_prod — recruit logic likely broken")
  }

  # T2. Recruit Dagb is positive and finite
  bad_rec <- df_stem_prod[is_recruit == TRUE & (!is.finite(Dagb_t_s) | Dagb_t_s <= 0), .N]
  if (bad_rec == 0) {
    pass("T2 all ingrowth/recruit Dagb_t_s are positive and finite")
  } else {
    fail(sprintf("T2 %d ingrowth/recruit rows have non-positive or non-finite Dagb_t_s", bad_rec))
  }

  # T3. Recruit gain in df_stem_demo equals sum from df_stem_prod
  demo_rec <- df_stem_demo[, sum(Dagb_recruit, na.rm = TRUE)]
  prod_rec <- df_stem_prod[
    is_recruit == TRUE & !is.na(quadrat) & quadrat != "" & !is.na(size),
    sum(Dagb_t_s, na.rm = TRUE)
  ]
  if (isTRUE(all.equal(demo_rec, prod_rec))) {
    pass("T3 ingrowth/recruit Dagb conserved across aggregation")
  } else {
    fail(sprintf("T3 ingrowth/recruit Dagb mismatch: demo=%.4f prod=%.4f", demo_rec, prod_rec))
  }

  # T4. Mortality only at last_census_alive and never at the final census
  bad_mort <- df_stem_mort[CensusID != last_census_alive, .N]
  if (bad_mort == 0) {
    pass("T4 all mortality rows are at last_census_alive")
  } else {
    fail(sprintf("T4 %d mortality rows NOT at last_census_alive", bad_mort))
  }

  max_c <- df_stem[, max(CensusID, na.rm = TRUE)]
  bad_mort_last <- df_stem_mort[CensusID == max_c, .N]
  if (bad_mort_last == 0) {
    pass(sprintf("T4b no mortality at final census %d", max_c))
  } else {
    fail(sprintf("T4b %d stems flagged dead at final census", bad_mort_last))
  }

  # T5. Interpolation effects: number of filled rows and any negative DBHs
  n_int <- df_stem[was_interpolated == TRUE, .N]
  cat(sprintf(
    "  [INFO] %d stem rows had dbh interpolated (method='%s')\n",
    n_int, dbh_interp_method
  ))
  # Check taper-corrected DBH (dbh_t) — that is the column that was interpolated.
  bad_int <- df_stem[was_interpolated == TRUE & dbh_t <= 0, .N]
  if (bad_int == 0) {
    pass("T5 no non-positive interpolated dbh_t values")
  } else {
    fail(sprintf("T5 %d interpolated dbh_t values are <= 0", bad_int))
  }

  # T6. AGB stock should be non-decreasing then decreasing only when (mort > prod);
  #     at minimum, every interval should satisfy:  agb(c+1) ~ agb(c) + (Dagb - DagbM)*dT
  #     within rounding. Here we check the plot-wide identity.
  plot_check <- df_stem_demo_quadrat[
    , .(
      agb_t = sum(agb_t), Dagb = sum(Dagb), DagbM = sum(DagbM), dT = mean(dT)
    ),
    .(CensusID)
  ]
  data.table::setorder(plot_check, CensusID)
  plot_check[, agb_next := shift(agb_t, type = "lead")]
  plot_check[, predicted_next := agb_t + (Dagb - DagbM) * dT]
  plot_check[, abs_err := abs(agb_next - predicted_next)]
  worst <- plot_check[!is.na(abs_err), max(abs_err)]
  cat(sprintf("  [INFO] Plot-wide AGB conservation max abs err across intervals: %.3f Mg\n", worst))

  # T7. Quadrat coverage: per quadrat there should be 1 row per census (1–8 typically)
  q_counts <- df_stem_demo_quadrat[, .N, .(quadrat)][, range(N)]
  cat(sprintf(
    "  [INFO] df_stem_demo_quadrat rows per quadrat range: %d-%d\n",
    q_counts[1], q_counts[2]
  ))

  # T8. Date imputation completeness: dT should be NA ONLY for each stem's first
  #     census (prev_ExactDate = NA by design). Any additional NA dT on non-first
  #     rows indicates a date that could not be imputed.
  # "First census" rows: those where shift(CensusID) gives NA within each stem.
  non_first_rows <- df_stem[, .SD[-.1], .(treeID, stemID)] # drop first row per stem
  n_na_dt_non_first <- non_first_rows[is.na(dT), .N]
  if (n_na_dt_non_first == 0L) {
    pass("T8 dT is non-NA for all non-first-census rows (dates fully imputed)")
  } else {
    fail(sprintf(
      "T8 %d non-first-census rows still have NA dT — ExactDate imputation incomplete",
      n_na_dt_non_first
    ))
  }
  rm(non_first_rows)

  cat("========== END TESTS ==========\n\n")
  invisible(NULL)
}

run_tests()

## ============================================================
## Plot-level summary: Mg ha-1 (stock) and Mg ha-1 yr-1 (fluxes)
## ============================================================
##
## Design:
##   Stock  — aggregated directly from df_stem_status for ALL 9 censuses
##            (every census with alive stems is included naturally).
##   Fluxes — aggregated from df_stem_demo_quadrat, which covers intervals
##            c → c+1 indexed at c. Census 1 (interval 1→2) is excluded
##            at user request; the final interval (8→9) is included.
##
## Spatial replication: 1250 quadrats of 20×20 m (400 m2).
## ha⁻¹ conversion: multiply by 10000/400 = 25.
##
## Summary per census: plot mean ± 95% CI across quadrats.
## ============================================================

ha_factor <- 10000 / 400 # 400 m2 to ha-1

# ExactDate has been filled in section 11 (quadrat median → plot median fallback),
# so dT values are based on real measurement dates, not a mean census year.

# get median year per census for plotting (x-axis labels)
census_yr_lut <- df_stem_status[
  !is.na(ExactDate),
  .(CensusYear = round(median(year(ExactDate), na.rm = TRUE))),
  by = CensusID
]

## Bootstrap CI helper (percentile method, resampling quadrats with replacement)
## Uses .SD[[col]] (not get(col)) to ensure group-restricted access in data.table.
bootstrap_ci <- function(data,
                         group_cols,
                         value_cols,
                         n_boot = 10000,
                         ci_level = 0.95,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- 1 - ci_level
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  result <- data[
    ,
    {
      out <- list()
      for (col in value_cols) {
        vals <- na.omit(.SD[[col]]) # group-restricted via .SD
        col_mean <- mean(.SD[[col]], na.rm = TRUE)
        if (length(vals) >= 2L) {
          boots <- replicate(
            n_boot,
            mean(sample(vals, size = length(vals), replace = TRUE))
          )
          col_lwr <- quantile(boots, probs = lower_q, names = FALSE)
          col_upr <- quantile(boots, probs = upper_q, names = FALSE)
        } else {
          col_lwr <- NA_real_
          col_upr <- NA_real_
        }
        out[[paste0(col, "_mean")]] <- col_mean
        out[[paste0(col, "_lwr")]] <- col_lwr
        out[[paste0(col, "_upr")]] <- col_upr
      }
      out[["n_quadrats"]] <- .N
      out
    },
    by = group_cols,
    .SDcols = value_cols
  ]
  return(result)
}

## --- Stock: all 9 censuses (CensusID 2–9) ----------

stock_q <- df_stem_status[
  !is.na(quadrat) & quadrat != "",
  .(agb_ha = sum(agb_t, na.rm = TRUE) * ha_factor),
  .(quadrat, CensusID)
]
stock_q <- merge(stock_q, census_yr_lut, by = "CensusID")

stock_summary <- bootstrap_ci(
  stock_q[CensusID >= 2L],
  group_cols = c("CensusID", "CensusYear"),
  value_cols = "agb_ha",
  n_boot = 10000,
  ci_level = 0.95,
  seed = 42
)
setnames(
  stock_summary,
  c("agb_ha_mean", "agb_ha_lwr", "agb_ha_upr"),
  c("agb_mean", "agb_lwr", "agb_upr")
)
data.table::setorder(stock_summary, CensusID)

## --- Fluxes: intervals 2→3 through 8→9 (CensusID 2–8) ----------------------

flux_q <- merge(
  df_stem_demo_quadrat[CensusID >= 2L & !is.na(quadrat) & quadrat != ""],
  census_yr_lut,
  by = "CensusID"
)
# Kohyama-corrected fluxes (Dagb_k == Dagb when use_kohyama19_correction = FALSE).
# Net AGB change = productivity - mortality: positive = biomass sink, negative = source.
flux_q[, prod_ha := Dagb_k * ha_factor]
flux_q[, mort_ha := DagbM_k * ha_factor]
flux_q[, net_ha  := prod_ha - mort_ha]

flux_summary <- bootstrap_ci(
  flux_q,
  group_cols = c("CensusID", "CensusYear"),
  value_cols = c("prod_ha", "mort_ha", "net_ha", "dT"),
  n_boot = 10000,
  ci_level = 0.95,
  seed = 42
)
setnames(
  flux_summary,
  c(
    "prod_ha_mean", "prod_ha_lwr", "prod_ha_upr",
    "mort_ha_mean", "mort_ha_lwr", "mort_ha_upr",
    "net_ha_mean",  "net_ha_lwr",  "net_ha_upr"
  ),
  c(
    "prod_mean", "prod_lwr", "prod_upr",
    "mort_mean", "mort_lwr", "mort_upr",
    "net_mean",  "net_lwr",  "net_upr"
  )
)
data.table::setorder(flux_summary, CensusID)

## --- Print summary table -----------------------------------------------------
# Fluxes use the Kohyama et al. 2019 correction when use_kohyama19_correction = TRUE.
# Net = productivity - mortality (Mg ha-1 yr-1)

summary_tbl <- merge(
  stock_summary[, .(CensusID, CensusYear,
    agb_mean = round(agb_mean, 1),
    agb_lwr  = round(agb_lwr, 1),
    agb_upr  = round(agb_upr, 1)
  )],
  flux_summary[, .(CensusID, CensusYear,
    prod_mean = round(prod_mean, 3),
    mort_mean = round(mort_mean, 3),
    net_mean  = round(net_mean,  3)
  )],
  by = c("CensusID", "CensusYear"), all = TRUE
)
cat("\n========== PLOT SUMMARY (Mg ha\u207b\u00b9 | Mg ha\u207b\u00b9 yr\u207b\u00b9) ==========\n")
print(summary_tbl)
cat("Expected: AGB ~220 Mg ha\u207b\u00b9 | productivity 2\u20138 | mortality comparable | net \u00b11\u20132\n")

## --- Figure A: Standing AGB stock (Mg ha⁻¹) ---------------------------------

library(ggplot2)

p_stock <- ggplot(stock_summary, aes(x = CensusYear, y = agb_mean)) +
  geom_ribbon(aes(ymin = agb_lwr, ymax = agb_upr), fill = "grey70", alpha = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    x = NULL,
    y = expression("AGB (Mg ha"^
      {
        -1
      } * ")"),
    title = "Standing AGB  (mean \u00b1 95% CI across quadrats)"
  ) +
  theme_bw()

## --- Figure B: Productivity & mortality (Mg ha⁻¹ yr⁻¹) ---------------------
# Fluxes are plotted at the midpoint of each interval c → c+1.

df_flux_long <- rbind(
  flux_summary[, .(CensusID, CensusYear,
    x = CensusYear + dT_mean / 2,
    flux = "Productivity",
    mean = prod_mean, lwr = prod_lwr, upr = prod_upr
  )],
  flux_summary[, .(CensusID, CensusYear,
    x = CensusYear + dT_mean / 2,
    flux = "Mortality",
    mean = mort_mean, lwr = mort_lwr, upr = mort_upr
  )]
)

p_flux <- ggplot(df_flux_long, aes(x = x, y = mean, colour = flux, fill = flux)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(Productivity = "steelblue", Mortality = "firebrick")) +
  scale_fill_manual(values = c(Productivity = "steelblue", Mortality = "firebrick")) +
  labs(
    x = "Census year",
    y = expression("Flux (Mg ha"^{
      -1
    } ~ "yr"^
      {
        -1
      } * ")"),
    colour = NULL, fill = NULL,
    title = "AGB fluxes  (mean \u00b1 95% CI across quadrats)"
  ) +
  theme_bw()

## --- Figure C: Net AGB change (Mg ha⁻¹ yr⁻¹) --------------------------------
# Net change = productivity - mortality. Positive = biomass sink; negative = source.
# Plotted at the midpoint of each census interval (same x convention as Figure B).

p_net <- ggplot(flux_summary, aes(x = CensusYear + dT_mean / 2, y = net_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = net_lwr, ymax = net_upr), fill = "forestgreen", alpha = 0.25) +
  geom_line(linewidth = 0.8, colour = "forestgreen") +
  geom_point(size = 2, colour = "forestgreen") +
  labs(
    x = "Census year",
    y = expression("Net AGB change (Mg ha"^{-1} ~ "yr"^{-1} * ")"),
    title = "Net AGB change  (productivity \u2212 mortality, mean \u00b1 95% CI)"
  ) +
  theme_bw()

library(cowplot)
plots <- plot_grid(p_stock, p_flux, p_net, ncol = 1, align = "v", labels = c("A", "B", "C"))

ggsave(
  plot = plots,
  filename = "./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/plot_agb_dynamics.png",
  width = 9, height = 12, units = "in", dpi = 300
)

## ============================================================
## Figure C: AGB stock and fluxes by size class
## ============================================================
## Uses df_stem_demo (per quadrat × size × CensusID) — already built.
## Size classes: [0,10), [10,20), [20,50), [50,500) cm DBH.
## ha⁻¹ conversion and bootstrap CIs identical to the plot-level summaries.

size_stock_q <- df_stem_status[
  !is.na(quadrat) & quadrat != "" & !is.na(size),
  .(agb_ha = sum(agb_t, na.rm = TRUE) * ha_factor),
  .(quadrat, CensusID, size)
]
size_stock_q <- merge(size_stock_q, census_yr_lut, by = "CensusID")

size_stock_summary <- bootstrap_ci(
  size_stock_q[CensusID >= 2L],
  group_cols = c("CensusID", "CensusYear", "size"),
  value_cols = "agb_ha",
  n_boot = 5000,
  ci_level = 0.95,
  seed = 42
)
setnames(
  size_stock_summary,
  c("agb_ha_mean", "agb_ha_lwr", "agb_ha_upr"),
  c("agb_mean", "agb_lwr", "agb_upr")
)
data.table::setorder(size_stock_summary, size, CensusID)

p_stock_size <- ggplot(
  size_stock_summary,
  aes(x = CensusYear, y = agb_mean, colour = size, fill = size)
) +
  geom_ribbon(aes(ymin = agb_lwr, ymax = agb_upr), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    x = NULL,
    y = expression("AGB (Mg ha"^
      {
        -1
      } * ")"),
    colour = "DBH class (cm)", fill = "DBH class (cm)",
    title = "Standing AGB by size class  (mean \u00b1 95% CI)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

## --- Flux summary by size class ------------------------------------------
size_flux_q <- df_stem_demo[
  !is.na(quadrat) & quadrat != "" & !is.na(size) & CensusID >= 2L
]
size_flux_q <- merge(size_flux_q, census_yr_lut, by = "CensusID")
size_flux_q[, prod_ha := Dagb_k * ha_factor]
size_flux_q[, mort_ha := DagbM_k * ha_factor]
size_flux_q[, net_ha  := prod_ha - mort_ha]

size_flux_summary <- bootstrap_ci(
  size_flux_q,
  group_cols  = c("CensusID", "CensusYear", "size"),
  value_cols  = c("prod_ha", "mort_ha", "net_ha", "dT"),
  n_boot      = 5000,
  ci_level    = 0.95,
  seed        = 42
)
setnames(
  size_flux_summary,
  c(
    "prod_ha_mean", "prod_ha_lwr", "prod_ha_upr",
    "mort_ha_mean", "mort_ha_lwr", "mort_ha_upr",
    "net_ha_mean",  "net_ha_lwr",  "net_ha_upr"
  ),
  c(
    "prod_mean", "prod_lwr", "prod_upr",
    "mort_mean", "mort_lwr", "mort_upr",
    "net_mean",  "net_lwr",  "net_upr"
  )
)
data.table::setorder(size_flux_summary, size, CensusID)

# Flux midpoint on the real year axis (CensusYear, not CensusID).
df_size_flux_long <- rbind(
  size_flux_summary[, .(CensusID, CensusYear, size,
    x = CensusYear + dT_mean / 2,
    flux = "Productivity",
    mean = prod_mean, lwr = prod_lwr, upr = prod_upr
  )],
  size_flux_summary[, .(CensusID, CensusYear, size,
    x = CensusYear + dT_mean / 2,
    flux = "Mortality",
    mean = mort_mean, lwr = mort_lwr, upr = mort_upr
  )]
)

p_flux_size <- ggplot(
  df_size_flux_long,
  aes(x = x, y = mean, colour = flux, fill = flux)
) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~size, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(Productivity = "steelblue", Mortality = "firebrick")) +
  scale_fill_manual(values = c(Productivity = "steelblue", Mortality = "firebrick")) +
  labs(
    x = "Census year",
    y = expression("Flux (Mg ha"^{
      -1
    } ~ "yr"^
      {
        -1
      } * ")"),
    colour = NULL, fill = NULL,
    title = "AGB fluxes by size class  (mean \u00b1 95% CI)"
  ) +
  theme_bw()

## --- Figure D4: Net AGB change by size class ---------------------------------
p_net_size <- ggplot(
  size_flux_summary,
  aes(x = CensusYear + dT_mean / 2, y = net_mean, colour = size, fill = size)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = net_lwr, ymax = net_upr), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~size, scales = "free_y", ncol = 2) +
  labs(
    x = "Census year",
    y = expression("Net AGB change (Mg ha"^{-1} ~ "yr"^{-1} * ")"),
    colour = "DBH class (cm)", fill = "DBH class (cm)",
    title = "Net AGB change by size class  (mean \u00b1 95% CI)"
  ) +
  theme_bw()

size_plots <- plot_grid(
  p_stock_size, p_flux_size, p_net_size,
  ncol = 1, align = "v", labels = c("D1", "D2", "D3")
)
ggsave(
  plot = size_plots,
  filename = "./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/plot_agb_by_size.png",
  width = 9, height = 14, units = "in", dpi = 300
)