# =============================================================================
# stem_dbh.R — the diameter of a stem in a census, shared by the scripts of
# 4_EXAMPLE_STRUCTURE_ASSESSMENT (biomass_stocks_fluxes.R,
# basal_area_uncertainty.R, general_plot_information.R), so that all three
# use the same rules. Sourced by each of them; needs data.table.
#
# A stem is in the stand (stocks, fluxes, stem counts) in a census only when it
# has a DBH at 1.3 m: measured, or filled from its own measurements.
#   - Taper: a DBH measured at another height is corrected to 1.3 m with Model
#     1 of Cushman et al. (2014) (taper_2014()), for the stocks and for the
#     fluxes as each script's use_taper_stocks and use_taper_fluxes set; without
#     it the DBH is kept as recorded.
#   - Stump: an alive (A) record without a DBH whose raw status (DFstatus) is
#     "broken below". The stem is alive but has no stem at 1.3 m: it gets no
#     DBH, and it splits the stem's life into segments (`seg`). The stem
#     leaves the stand at a stump; a later measurement starts a new segment
#     (a recruit).
#   - Fill: an alive (A) census without a DBH, not a stump, that comes after a
#     measurement of the same segment (`fill`). Between two measurements its
#     DBH is linear in time; after the segment's last measurement it continues
#     the growth trend of the last two (dbh_from_own_measurements()). An alive
#     census before the stem's first measurement, after a stump with no
#     measurement since, or of a stem never measured gets no DBH.
# =============================================================================

# ---- Taper correction --------------------------------------------------------
# taper_2014(): DBH at `common_hom` from a DBH measured at height `hom`, with
# Model 1 of Cushman et al. (2014) and its parameter b from their Table 3:
#   b      = exp(-2.0205 + species term - 0.5053 * log(dbh_cm) + 0.3748 * log(hom))
#   dbh_at = dbh_cm * exp(b * (hom - common_hom))
# The species terms (TAPER_SPECIES_2014) are the deviations of their ten focal
# species from the base class (every other species), by BCI code; quaras and
# tet2pa are now named Quararibea stenophylla and Protium stevensonii. With
# sp = NULL every stem takes the base class. dbh_mm in mm, hom in m (NA is
# read as common_hom). Returns mm; NA where the DBH or the HOM is not
# positive. Stops when the vectors differ in length.
TAPER_SPECIES_2014 <- c(
    quaras = 0.4865, ceibpe = 0.0495, huracr = 0.0086, diptpa = 0.0150, alsebl = 0.1440,
    brosal = 0.3874, cavapl = 0.7982, anacex = 0.4371, tet2pa = -0.5141, tab2ar = 0.3576
)
taper_2014 <- function(dbh_mm, hom, sp = NULL, common_hom = 1.3) {
    if (length(dbh_mm) != length(hom) || (!is.null(sp) && length(sp) != length(hom))) {
        stop("'dbh_mm', 'hom' and 'sp' must have the same length")
    }
    dbh_mm <- as.numeric(dbh_mm)
    hom <- as.numeric(hom)
    hom[is.na(hom)] <- common_hom
    dbh_cm <- dbh_mm / 10
    dbh_cm[dbh_cm <= 0] <- NA_real_
    hom_for_log <- hom
    hom_for_log[hom_for_log <= 0] <- NA_real_
    sp_term <- if (is.null(sp)) 0 else fcoalesce(unname(TAPER_SPECIES_2014[sp]), 0)
    b <- exp(-2.0205 + sp_term - 0.5053 * log(dbh_cm) + 0.3748 * log(hom_for_log))
    out_mm <- 10 * dbh_cm * exp(b * (hom - common_hom))
    out_mm[is.na(out_mm) | is.infinite(out_mm)] <- NA_real_
    out_mm
}

# Taper-corrected DBH (mm) of the rows of `dt` (columns dbh in mm, hom in m,
# sp, Family): the species terms when species_terms is TRUE; palms keep the
# recorded DBH unless taper_palms is TRUE (they have no buttress taper, and
# Socratea is measured above its stilt roots).
stem_dbh_taper <- function(dt, species_terms = TRUE, taper_palms = FALSE) {
    out <- taper_2014(dt$dbh, dt$hom, if (species_terms) dt$sp else NULL)
    if (!taper_palms) out[dt$Family %in% "Arecaceae"] <- dt$dbh[dt$Family %in% "Arecaceae"]
    out
}

# ---- Giant strangler figs ----------------------------------------------------
# Trees of four strangler Ficus species with a stem above 500 mm DBH (recorded
# DBH, any census; Rutishauser et al. 2020): standard allometries do not
# describe them. Needs dbh (mm), Genus and Species.
GIANT_STRANGLER_SPECIES <- c("costaricana", "obtusifolia", "popenoei", "trigonata")
giant_strangler_trees <- function(dt) {
    unique(dt[dbh > 500 & Genus %in% "Ficus" & Species %in% GIANT_STRANGLER_SPECIES, treeID])
}

# ---- Stumps, segments and the censuses to fill -------------------------------
# Adds by reference, per stem (treeID + stemID) in census order:
#   stump : alive (A), no DBH, raw status broken below
#   seg   : stump segment (number of stumps up to this census)
#   fill  : alive (A), no DBH, not a stump, after a measurement of its segment
# `dbh_col` is the recorded DBH (NA on records that are not alive).
stump_segments <- function(dt, dbh_col = "dbh") {
    setorder(dt, treeID, stemID, CensusID)
    dt[, stump := Rstatus %in% "A" & is.na(get(dbh_col)) & DFstatus %in% "broken below"]
    dt[, seg := cumsum(stump), by = .(treeID, stemID)]
    dt[, fill := Rstatus %in% "A" & is.na(get(dbh_col)) & !stump &
        cumsum(!is.na(get(dbh_col))) > 0L, by = .(treeID, stemID, seg)]
    invisible(dt)
}

# Censuses to fill that lie between two measurements of their stem but not of
# their segment (a stump in between): none is expected; the scripts check it,
# because basal_area_uncertainty.R interpolates the gaps of the whole stem.
fill_across_stump <- function(dt, dbh_col = "dbh") {
    m <- dt[!is.na(get(dbh_col)), .(last_meas_stem = max(CensusID)), by = .(treeID, stemID)]
    s <- dt[!is.na(get(dbh_col)), .(last_meas_seg = max(CensusID)), by = .(treeID, stemID, seg)]
    x <- dt[fill == TRUE, .(treeID, stemID, seg, CensusID)]
    x <- m[x, on = .(treeID, stemID)]
    x <- s[x, on = .(treeID, stemID, seg)]
    x[CensusID > last_meas_seg & CensusID < last_meas_stem]
}

# ---- DBH of a census to fill -------------------------------------------------
# DBH at dates t_new from the measurements of the same segment (dates t_meas,
# DBHs d_meas): between two measurements linear in time; after the last one
# the growth trend of the last two continued in time (a negative trend,
# usually a measurement or POM problem, is not extrapolated; one measurement:
# kept). A date before the first measurement is not filled (stops).
dbh_from_own_measurements <- function(t_meas, d_meas, t_new) {
    o <- order(t_meas)
    t_meas <- t_meas[o]
    d_meas <- d_meas[o]
    n <- length(t_meas)
    if (any(t_new < t_meas[1L])) stop("dbh_from_own_measurements(): a date before the first measurement is not filled")
    if (n == 1L) {
        return(rep(d_meas, length(t_new)))
    }
    out <- approx(t_meas, d_meas, xout = t_new)$y
    after <- t_new > t_meas[n]
    if (any(after)) {
        growth <- max((d_meas[n] - d_meas[n - 1L]) / (t_meas[n] - t_meas[n - 1L]), 0)
        out[after] <- d_meas[n] + growth * (t_new[after] - t_meas[n])
    }
    out
}

# Fills `col` (a DBH) on the rows with fill == TRUE from the stem's measurements
# of the same segment, at the dates in `date_col`; returns the number of rows
# filled. Needs stump_segments() first.
fill_stem_dbh <- function(dt, col, date_col = "ExactDate") {
    dt[, (col) := {
        x <- get(col)
        target <- fill & is.na(x)
        if (any(target)) {
            t_num <- as.numeric(get(date_col))
            meas <- !is.na(x)
            x[target] <- dbh_from_own_measurements(t_num[meas], x[meas], t_num[target])
        }
        x
    }, by = .(treeID, stemID, seg)]
    invisible(dt[fill == TRUE & !is.na(get(col)), .N])
}
