# ==============================================================================
# BCI Stem Reconstruction — Basal Area: Stem-Identity Uncertainty
# ==============================================================================
#
# PURPOSE
# -------
# Propagates STEM-IDENTITY uncertainty (the posterior reconstruction paths of
# the dp_global engine) into basal area (BA) stocks and fluxes of the BCI 50-ha
# plot across nine stem censuses (1982–2022/3). Only identity uncertainty is
# quantified here; spatial sampling uncertainty is deliberately not included.
#
# WHAT IDENTITY UNCERTAINTY CAN AND CANNOT CHANGE
# -----------------------------------------------
# A reconstruction path says which measurement belongs to which stem across
# censuses. Every path contains the SAME measurements; paths differ only in how
# those measurements are linked into stems. Therefore:
#   • BA stock of measured stems is identical in every realization;
#   • the split of BA change into Growth (survivors), Loss (deaths) and Gain
#     (recruits) is what identity uncertainty changes;
#   • stock and net change can differ only through stems that are alive but
#     unmeasured in a census (gaps, interpolated in time; see fill_stem_gaps()).
# These properties are enforced with hard checks (bio_check()).
#
# DESIGN
# ------
# Exported reconstruction : the stem IDs of the R tables (DP trees: the DP's
#                           Viterbi decoding; probabilistic trees: the most
#                           representative of their draws; both plus stage-2
#                           post-processing). Used for every tree with no
#                           identity uncertainty and as the reference line in
#                           the figures.
# Posterior paths         : 200 posterior draws per tag, collapsed into unique
#                           paths with their sample counts
#                           (DATA/POSTERIORS/posterior_sampled_paths.rds).
#                           Path weights = path_count / sum(path_count), i.e.
#                           every one of a tree's 200 draws is equally likely.
#                           (path_prob is NOT used: for DP trees it re-weights
#                           each draw by its own probability and double counts.)
# Two engines, one sample : both engines' draws are used together, the same
#                           way. DP trees (dp_global_dp.R) are drawn by backward
#                           sampling from the exact DP posterior, so likely
#                           trajectories repeat (path_count > 1). Trees routed
#                           to the probabilistic engine (dp_probabilistic_matching.R:
#                           palms and other forced species, stranglers, trees
#                           whose state space is too large) are drawn
#                           approximately: a noisy assignment per census pair,
#                           stitched, growth-repaired and filtered by pins.
#                           Almost every such draw differs somewhere, so nearly
#                           every path has path_count = 1 (weight 1/200). Either
#                           way, one draw per tree per realization is taken with
#                           probability path_count / 200. The probabilistic
#                           draws are an approximation, not a calibrated
#                           posterior; their trees are counted separately in
#                           the diagnostics.
# Splice                  : a path covers the DP window of its tree. Measured
#                           observations outside that window were placed
#                           deterministically by stage 2 (DB StemID kept,
#                           broken-below splits, probabilistic / enumeration
#                           fallbacks) and keep those links: each joins the
#                           path stem that holds its exported stem-mate at the
#                           nearest window edge.
# Monte Carlo             : each realization draws one path per multi-path tree
#                           (independently, by weight) and aggregates BA at
#                           quadrat and plot level.
#
# ANCHOR CENSUSES AND SCOPE
# -------------------------
# The DP reconstructs identities backward from the anchor census
# (ANCHOR_START_CENSUS = Census 7). Post-anchor censuses (C8, C9) have confirmed
# stem IDs: their stocks and the C7→C8 and C8→C9 fluxes are the same in every
# realization, except for the rare stem whose unmeasured gap spans the anchor
# (its interpolated DBH depends on the path).
#
# OUTPUTS (written to BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs/)
# --------
#   ba_map_change_treeID.feather        exported reconstruction: tree-level flux
#   ba_map_change_quadrat.feather       exported reconstruction: quadrat-level flux
#   ba_map_stock_quadrat.feather        exported reconstruction: quadrat-level stock
#   ba_mc_realizations_quadrat.feather  MC realizations: quadrat-level fluxes
#   ba_mc_realizations_stock_quadrat.feather  MC realizations: quadrat-level stocks
#   ba_mc_summary_quadrat.feather       MC mean / sd / 95 % interval of quadrat fluxes
#   ba_mc_summary_stock_quadrat.feather MC mean / sd / 95 % interval of quadrat stocks
#   ba_mc_diagnostics.txt               tree counts, splice counts, checks, MC error
#   ba_mc_realizations_treeID/          optional per-realization tree-level fluxes
#   fig1_BA_stock.pdf                   BA stock: exported vs identity MC
#   fig2_BA_fluxes.pdf                  BA fluxes: exported vs identity MC
#   fig3_BA_trajectories.pdf            tree BA trajectories: exported vs posterior paths
#
# NOTES
# --------
# Additional methodological decisions are documented in biomass_stocks_fluxes.R.
# In this script, giant strangler ficus (> 500 mm DBH) are retained, no
# correction is applied for the buttress bias of the first census, the Kohyama
# et al. (2019) correction is not applied, and palm diameters are not modified.
# ==============================================================================

rm(list = ls())

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(collapse)
library(arrow)

workspace_root <- getwd()

# ============================================================
# SECTION 1: Configuration and data loading
# ============================================================

# ── Anchor census ──────────────────────────────────────────────────────────────
# First census with confirmed stem identities (2010 = Census 7). The DP samples
# identity paths backward from this anchor; censuses after it are deterministic.
ANCHOR_START_CENSUS <- 7L

# The first census is omitted from figures.
first_plot_census <- 2L

# Number of Monte Carlo realizations and seed (recorded in the diagnostics).
K_realizations <- 100L
mc_seed <- 42L

# Center of the MC distribution drawn as a dashed line ("mean" or "median").
mc_center <- "median"
mc_center <- match.arg(mc_center, c("mean", "median"))

# Optional heavy outputs.
write_tree_realizations <- FALSE # one tree-level feather per realization
write_quadrat_realizations <- TRUE # all quadrat-level realizations in one feather

# Report the stage-2 method of spliced observations.
report_splice_methods <- TRUE

out_dir <- file.path(workspace_root, "BCI_stem_reconstruction", "4_EXAMPLE_STRUCTURE_ASSESSMENT", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
diag_file <- file.path(out_dir, "ba_mc_diagnostics.txt")
if (file.exists(diag_file)) file.remove(diag_file)
diag_line <- function(...) cat(..., "\n", file = diag_file, append = TRUE, sep = "")

# ── Hard check that stays visible in interactive (line-by-line) runs ────────
# On failure: prints a ❌ line (count + examples), raises an immediate warning,
# writes the same text to the diagnostics file, and only then stops.
bio_check <- function(ok, msg, examples = NULL, n_bad = NULL) {
    if (isTRUE(all(ok))) {
        cat("✓", msg, "\n")
        diag_line("CHECK OK: ", msg)
        return(invisible(TRUE))
    }
    n_txt <- if (!is.null(n_bad)) sprintf(" [%d case(s)]", n_bad) else ""
    ex_txt <- if (length(examples) > 0L) {
        paste0(" | examples: ", paste(head(unique(examples), 10), collapse = ", "))
    } else {
        ""
    }
    full_msg <- paste0("CHECK FAILED: ", msg, n_txt, ex_txt)
    cat("❌", full_msg, "\n")
    warning(full_msg, call. = FALSE, immediate. = TRUE)
    diag_line(full_msg)
    stop(full_msg, call. = FALSE)
}

bci_stem_nums <- as.character(1:9)
census_list <- lapply(bci_stem_nums, function(num) {
    fp <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RTABLES", paste0("bci.stem", num, ".Rdata"))
    if (!file.exists(fp)) stop("Missing census file: ", fp)
    load(fp)
    get(paste0("bci.stem", num))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)
rec <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rec <- rec[!is.na(quadrat)]
rec[, CensusID := as.integer(CensusID)]
rm(census_list, bci_stem_nums)

# Stage 3 exports every DBH as recorded, so a dead (G/D) record can carry a
# DBH. It is not a living stem: stocks and fluxes use alive (A) rows only, so
# that DBH is ignored here (a measured row below means an alive measurement).
n_dead_dbh <- rec[Rstatus != "A" & !is.na(dbh), .N]
rec[Rstatus != "A" & !is.na(dbh), dbh := NA]
cat(sprintf("[BA] %d DBH values recorded on dead (G/D) records are not used\n", n_dead_dbh))
rm(n_dead_dbh)

# stemID is the stem number within its tree (1, 2, ...): treeID + stemID
# identify a stem, so every stem key below uses both.

# Stage 3 dates every alive (A) row, measured row and first G/D row; P rows
# and later G/D rows keep NA unless a date was recorded. Those get the modal
# field date of their tree in that census, else of their quadrat, else of the
# census (the rule stage 3 applied to every row before), so census dates and
# path gaps are computed as before. Same function in biomass_stocks_fluxes.R
# and general_plot_information.R.
date_mode_by <- function(dt, by_cols) {
    cnt <- dt[!is.na(ExactDate), .N, by = c(by_cols, "ExactDate")]
    setorderv(cnt, c(by_cols, "N", "ExactDate"), c(rep(1L, length(by_cols)), -1L, 1L))
    cnt[cnt[, .I[1L], by = by_cols]$V1, c(by_cols, "ExactDate"), with = FALSE]
}
fill_missing_dates <- function(dt) {
    m_tree <- date_mode_by(dt, c("treeID", "CensusID"))
    m_quad <- date_mode_by(dt[!is.na(quadrat)], c("quadrat", "CensusID"))
    m_cens <- date_mode_by(dt, "CensusID")
    dt[m_tree, on = .(treeID, CensusID), ExactDate := fcoalesce(ExactDate, i.ExactDate)]
    dt[m_quad, on = .(quadrat, CensusID), ExactDate := fcoalesce(ExactDate, i.ExactDate)]
    dt[m_cens, on = .(CensusID), ExactDate := fcoalesce(ExactDate, i.ExactDate)]
    invisible(dt)
}
bio_check(
    rec[Rstatus == "A" | !is.na(dbh), !anyNA(ExactDate)],
    "Every alive or measured stem-census row has an ExactDate (dated in stage 3)",
    n_bad = rec[(Rstatus == "A" | !is.na(dbh)) & is.na(ExactDate), .N]
)
n_undated <- rec[is.na(ExactDate), .N]
fill_missing_dates(rec)
bio_check(
    rec[, !anyNA(ExactDate)],
    sprintf("Every stem-census row has an ExactDate (%d P / later G/D rows dated from tree, quadrat or census)", n_undated),
    n_bad = rec[is.na(ExactDate), .N]
)
rm(n_undated)

# Cushman et al. 2014
taper_2014 <- function(dbh_mm, hom, common_hom = 1.3) {
    # Defensive checks
    if (length(dbh_mm) != length(hom)) {
        stop("'dbh_mm' and 'hom' must have the same length")
    }
    # copy inputs to avoid modifying caller's vectors
    dbh_mm <- as.numeric(dbh_mm)
    hom <- as.numeric(hom)
    # Replace NA heights with 1.3 m (do not modify valid measured heights)
    hom_na <- is.na(hom)
    hom[hom_na] <- common_hom
    # convert dbh from mm to cm for the model
    dbh_cm <- dbh_mm / 10
    # Protect against log(0) or negative inputs by coercing non-positive values to NA
    dbh_cm[dbh_cm <= 0] <- NA_real_
    hom_for_log <- hom
    hom_for_log[hom_for_log <= 0] <- NA_real_
    b <- exp(-2.0205 - 0.5053 * log(dbh_cm) + 0.3748 * log(hom_for_log))
    out <- dbh_cm / (exp(-b * (hom - common_hom)))
    # convert back to mm and set invalid values to NA
    out_mm <- out * 10
    out_mm[is.na(out_mm) | is.infinite(out_mm)] <- NA_real_
    return(out_mm)
}

# NOTE: dbh should be in cm for the equation.
rec[, hom := ifelse(is.na(hom), 1.3, hom)]
rec[, dbh_t := taper_2014(dbh_mm = dbh, hom = hom)]
rec[, dbh_raw := dbh]
rec[, dbh := fifelse(!is.na(dbh_t), dbh_t, dbh_raw)]

post_file <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "POSTERIORS", "posterior_sampled_paths.rds")
stage2_file <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "PROCESSED", "complete_dataset_final_with_reconstructed_stemids.rds")

# ============================================================
# SECTION 2: Helper functions
# ============================================================
# ba_m2(dbh_mm)     – stem basal area in m² (π/4 × (d/1000)²).
# fill_stem_gaps()  – adds a row for every census in which a stem is alive but
#                     unmeasured (between two of its measurements) with DBH
#                     interpolated linearly in time, flagged `interpolated`.
# decompose_ba()    – classifies each stem of an interval as survivor, death
#                     or recruit and returns Growth / Loss / Gain.
# summarise_flux()  – mean / sd / 95 % interval of fluxes across realizations.
# ============================================================

ba_m2 <- function(dbh_mm) pi / 4 * (dbh_mm / 1000)^2

# A stem measured at censuses a and b is alive at every census between them
# (stage-3 lifespan rule: alive later => never dead). Without a row for those
# censuses the flux decomposition would count a death followed by a new
# recruit. Missing rows are added with DBH interpolated linearly in time
# between the neighbouring measurements (measurement dates).
# Rstatus "A" without a DBH means the stem was alive and its measurement was
# missed in the field, so it stays in the stock until it is dead (D / G). A
# census in `alive` after the stem's last measurement gets the growth trend of
# its last two measurements continued in time (a negative trend is not
# extrapolated; a stem measured once keeps that value); one before its first
# measurement gets the first measurement. This is dbh_from_own_measurements()
# of general_plot_information.R and biomass_stocks_fluxes.R.
#   dt        : one row per stem × census with a measured `dbh` and its date `t`
#   stem_cols : columns identifying a stem (e.g. treeID, stemID)
#   gap_dates : date `t` of a missing census, joined on its other columns
#               (the stem's own row date, or the tree's date in that census)
#   alive     : stem_cols + CensusID of censuses in which the stem is alive
#               without a DBH (NULL: none)
# Added rows are flagged `interpolated` (not a measurement); those beyond the
# measured span are also flagged `carried`.
fill_stem_gaps <- function(dt, stem_cols, gap_dates, alive = NULL) {
    dt <- copy(dt)
    dt[, `:=`(interpolated = FALSE, carried = FALSE)]
    setkeyv(dt, c(stem_cols, "CensusID"))
    date_keys <- setdiff(names(gap_dates), "t")
    rng <- dt[, .(c_min = min(CensusID), c_max = max(CensusID), n = .N), by = stem_cols]
    added <- list()

    # 1. Censuses between two measurements: linear in time.
    gaps <- rng[c_max - c_min + 1L > n]
    if (nrow(gaps) > 0L) {
        full <- gaps[, .(CensusID = seq.int(c_min, c_max)), by = stem_cols]
        miss <- full[!dt, on = c(stem_cols, "CensusID")]
        prev <- dt[miss, on = c(stem_cols, "CensusID"), roll = Inf, .(d_prev = x.dbh, t_prev = x.t)]
        nxt <- dt[miss, on = c(stem_cols, "CensusID"), roll = -Inf, .(d_next = x.dbh, t_next = x.t)]
        miss[, c("d_prev", "t_prev") := prev]
        miss[, c("d_next", "t_next") := nxt]
        miss[gap_dates, on = date_keys, t := i.t]
        miss[, dbh := d_prev + (d_next - d_prev) * (t - t_prev) / (t_next - t_prev)]
        bad <- miss[is.na(dbh) | !(t_prev < t & t < t_next)]
        bio_check(
            nrow(bad) == 0L,
            "Every alive-but-unmeasured census lies in time between two measurements of its stem and gets an interpolated DBH",
            examples = bad$treeID, n_bad = nrow(bad)
        )
        added$gap <- miss[, c(stem_cols, "CensusID", "dbh", "t"), with = FALSE][, `:=`(interpolated = TRUE, carried = FALSE)]
    }

    # 2. Alive censuses outside the measured span.
    if (!is.null(alive) && nrow(alive) > 0L) {
        carry <- unique(alive[, c(stem_cols, "CensusID"), with = FALSE])
        carry <- rng[carry, on = stem_cols, nomatch = 0L][CensusID > c_max | CensusID < c_min]
        if (nrow(carry) > 0L) {
            carry[gap_dates, on = date_keys, t := i.t]
            # last measurement (after the span) and the one before it
            lastm <- dt[carry, on = c(stem_cols, "CensusID"), roll = Inf, .(c_l = x.CensusID, d_l = x.dbh, t_l = x.t)]
            carry[, c("c_l", "d_l", "t_l") := lastm]
            prevq <- carry[, c(stem_cols, "CensusID"), with = FALSE][, CensusID := carry$c_l - 1L]
            prevm <- dt[prevq, on = c(stem_cols, "CensusID"), roll = Inf, .(d_p = x.dbh, t_p = x.t)]
            carry[, c("d_p", "t_p") := prevm]
            # first measurement (before the span)
            carry[, d_f := dt[carry, on = c(stem_cols, "CensusID"), roll = -Inf, x.dbh]]
            carry[, growth := fifelse(is.na(d_p), 0, pmax((d_l - d_p) / (t_l - t_p), 0))]
            carry[, dbh := fifelse(CensusID > c_max, d_l + growth * (t - t_l), d_f)]
            bad <- carry[is.na(dbh) | (CensusID > c_max & !(t > t_l))]
            bio_check(
                nrow(bad) == 0L,
                "Every alive census outside a stem's measured span gets a DBH from the stem's own measurements",
                examples = bad$treeID, n_bad = nrow(bad)
            )
            added$carry <- carry[, c(stem_cols, "CensusID", "dbh", "t"), with = FALSE][, `:=`(interpolated = TRUE, carried = TRUE)]
        }
    }

    if (length(added) == 0L) {
        return(dt)
    }
    # added rows take the stem's tree-level attributes (quadrat) from its measured rows
    add <- rbindlist(added, use.names = TRUE)
    attr_dt <- unique(dt[, c(stem_cols, intersect(names(dt), "quadrat")), with = FALSE])
    add <- attr_dt[add, on = stem_cols]
    rbindlist(list(dt, add), use.names = TRUE, fill = TRUE)
}

decompose_ba <- function(m, by_cols) {
    m[, status := fifelse(
        !is.na(BA_from) & !is.na(BA_to), "survivor",
        fifelse(!is.na(BA_from), "death", "recruit")
    )]
    # fifelse(..., 0) keeps every component defined (0, not NA) when a tree
    # has no stem of that class in the interval.
    m[, .(
        Growth_BA     = fsum(fifelse(status == "survivor", BA_to - BA_from, 0)),
        Loss_BA       = -fsum(fifelse(status == "death", BA_from, 0)),
        Gain_BA       = fsum(fifelse(status == "recruit", BA_to, 0))
    ), by = by_cols]
}

summarise_flux <- function(dt, by_cols) {
    dt[,
        {
            qg <- fquantile(Growth_BA, c(0.025, 0.975))
            ql <- fquantile(Loss_BA, c(0.025, 0.975))
            qa <- fquantile(Gain_BA, c(0.025, 0.975))
            qd <- fquantile(DeltaBA_total, c(0.025, 0.975))
            .(
                Growth_mean = fmean(Growth_BA),     Growth_sd  = fsd(Growth_BA),
                Growth_lwr  = qg[1L],               Growth_upr = qg[2L],
                Loss_mean   = fmean(Loss_BA),       Loss_sd    = fsd(Loss_BA),
                Loss_lwr    = ql[1L],               Loss_upr   = ql[2L],
                Gain_mean   = fmean(Gain_BA),       Gain_sd    = fsd(Gain_BA),
                Gain_lwr    = qa[1L],               Gain_upr   = qa[2L],
                Delta_mean  = fmean(DeltaBA_total), Delta_sd   = fsd(DeltaBA_total),
                Delta_lwr   = qd[1L],               Delta_upr  = qd[2L]
            )
        },
        by = by_cols
    ]
}

# Per-interval decomposition of a stem table (stems identified by stem_cols;
# `tree_cols` are the grouping columns kept in the result).
decompose_intervals <- function(stems, pairs, stem_cols, tree_cols) {
    rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
        cf <- pairs$CensusID_from[i]
        ct <- pairs$CensusID_to[i]
        sf <- stems[CensusID == cf, c(stem_cols, "BA"), with = FALSE]
        st <- stems[CensusID == ct, c(stem_cols, "BA"), with = FALSE]
        setnames(sf, "BA", "BA_from")
        setnames(st, "BA", "BA_to")
        d <- decompose_ba(merge(sf, st, by = stem_cols, all = TRUE), by_cols = tree_cols)
        d[, `:=`(CensusID_from = cf, CensusID_to = ct)]
    }))
}

# ============================================================
# SECTION 3: Exported reconstruction
# ============================================================
# The stem IDs of the R tables are the exported reconstruction. Gaps (alive
# but unmeasured censuses) are filled before the decomposition so that they
# are not counted as a death plus a recruitment.
# Outputs:
#   map_tree_change    – tree-level flux per census pair
#   map_quadrat_change – quadrat-level flux
#   map_quadrat_stock  – quadrat-level BA stock
# (object names keep the historical "map_" prefix.)
# ============================================================

dates <- rec[, .(Date = median(ExactDate)), by = CensusID][order(CensusID)]
dates[, Year := as.integer(format(Date, "%Y"))]
bio_check(nrow(dates) >= 2L, "At least two censuses are available")

census_pairs <- data.table(
    CensusID_from = dates$CensusID[-nrow(dates)],
    CensusID_to   = dates$CensusID[-1L],
    Date_from     = dates$Date[-nrow(dates)],
    Date_to       = dates$Date[-1L]
)
census_pairs[, Date_mid := Date_from + (Date_to - Date_from) / 2]
census_pairs[, Year_mid := as.integer(format(Date_mid, "%Y"))]
census_pairs[, Interval_yr := as.numeric(difftime(Date_to, Date_from, units = "days")) / 365.25]

stem_obs <- rec[
    !is.na(dbh) & !is.na(treeID) & !is.na(stemID),
    .(quadrat, treeID, stemID, CensusID, StemPaths, dbh, t = as.numeric(ExactDate))
]
bio_check(
    stem_obs[, !anyDuplicated(stem_obs[, .(treeID, stemID, CensusID)])],
    "One measurement per stem and census in the exported reconstruction"
)

# Dates of a census in which a stem is alive but unmeasured: the stem's own
# row date (exported stems), or the tree's modal date in that census (stems
# of a posterior path, which have no row of their own).
stem_dates <- rec[, .(treeID, stemID, CensusID, t = as.numeric(ExactDate))]
tree_dates <- rec[, .N, by = .(treeID, CensusID, ExactDate)][order(-N)][
    , .(t = as.numeric(ExactDate[1L])),
    by = .(treeID, CensusID)
]

# Censuses in which a stem is alive (A) but its measurement was missed.
alive_nodbh <- rec[Rstatus == "A" & is.na(dbh), .(treeID, stemID, CensusID)]
exp_stems <- fill_stem_gaps(stem_obs, c("treeID", "stemID"), stem_dates, alive = alive_nodbh)
exp_stems[, BA := ba_m2(dbh)]
cat(sprintf(
    "[BA] exported: %d measured stem-censuses | alive without DBH: %d (interpolated %d, trend continued after the last measurement %d, never measured %d)\n",
    nrow(stem_obs), nrow(alive_nodbh), exp_stems[interpolated & !carried, .N], exp_stems[carried == TRUE, .N],
    alive_nodbh[!stem_obs, on = .(treeID, stemID), .N]
))

# An alive stem is in the stock in every census it is alive; a stem is never
# in the stock in a census where it is not alive.
alive_rows <- rec[Rstatus == "A", .(treeID, stemID, CensusID)][
    unique(stem_obs[, .(treeID, stemID)]),
    on = .(treeID, stemID), nomatch = 0L
]
not_in_stock <- alive_rows[!exp_stems, on = .(treeID, stemID, CensusID)]
bio_check(
    nrow(not_in_stock) == 0L,
    "Every alive (A) census of a measured stem is in the stock (measured, interpolated or trend-continued)",
    examples = not_in_stock[, paste(treeID, stemID, sep = "_")], n_bad = nrow(not_in_stock)
)
not_alive <- exp_stems[, .(treeID, stemID, CensusID)][!alive_rows, on = .(treeID, stemID, CensusID)]
bio_check(
    nrow(not_alive) == 0L,
    "No stem is in the stock in a census where it is not alive (A)",
    examples = not_alive[, paste(treeID, stemID, sep = "_")], n_bad = nrow(not_alive)
)
rm(alive_rows, not_in_stock, not_alive)

tree_census <- exp_stems[,
    .(TotalBA_m2 = fsum(BA), MeasuredBA_m2 = fsum(BA * !interpolated), NumStems = .N),
    by = .(quadrat, treeID, CensusID)
]
setorder(tree_census, treeID, CensusID)

map_change <- decompose_intervals(
    exp_stems, census_pairs,
    stem_cols = c("quadrat", "treeID", "stemID"), tree_cols = c("quadrat", "treeID")
)
map_change <- census_pairs[, .(CensusID_from, CensusID_to, Date_from, Date_to)][map_change, on = .(CensusID_from, CensusID_to)]
map_change[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
cat("[BA] exported decomposition:", nrow(map_change), "tree-intervals across", uniqueN(map_change$treeID), "treeIDs\n")

flux_cols <- c("Growth_BA", "Loss_BA", "Gain_BA")
map_tree_change <- map_change
map_quadrat_change <- map_tree_change[,
    lapply(.SD, fsum, na.rm = TRUE),
    .SDcols = flux_cols,
    by = .(quadrat, CensusID_from, CensusID_to)
]
map_quadrat_change[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
map_quadrat_stock <- tree_census[, .(TotalBA_m2 = fsum(TotalBA_m2)), by = .(quadrat, CensusID)]
cat("[BA] exported quadrat stock:", nrow(map_quadrat_stock), "quadrat×census rows\n")

# ============================================================
# SECTION 4: Identity uncertainty (Monte Carlo over posterior paths)
# ============================================================

# ---- 4.1 Posterior paths and weights ---------------------------------------
post_full <- as.data.table(readRDS(post_file))
post_full[, treeID := as.character(treeID)]
post_full[, n_paths := .N, by = treeID]
# Sample frequencies are the posterior probabilities of the unique paths, for
# both engines (DP draws often repeat; probabilistic draws are nearly all
# unique, so their paths weigh 1/200 each).
post_full[, w := path_count / sum(path_count), by = treeID]
bio_check(
    post_full[, abs(sum(w) - 1) < 1e-9, by = treeID][, all(V1)],
    "Path weights (path_count / sum(path_count)) sum to 1 within every tree"
)
cat(
    "[BA] posterior:", nrow(post_full), "paths for", uniqueN(post_full$treeID), "trees |",
    post_full[n_paths > 1L, uniqueN(treeID)], "trees with >1 path\n"
)

multi <- post_full[n_paths > 1L]
multi[, path_idx := seq_len(.N), by = treeID]
multi_trees <- unique(multi$treeID)

# Stage-2 reconstruction method of every observation (engine of each tree).
s2 <- as.data.table(readRDS(stage2_file))[, .(
    treeID = as.character(TreeID), StemPaths = as.integer(obs_row_id),
    method = as.character(ReconstructionMethod)
)]
# Both engines are sampled together; the engine is only used to report them.
prob_trees <- intersect(multi_trees, s2[method == "probabilistic", unique(treeID)])
# Single-stem tags (no reconstruction, hence no posterior), for the diagnostics.
single_stem_trees <- s2[method == "single_stem_tag_no_reconstructed", unique(treeID)]
cat(
    "[BA] multi-path trees:", length(multi_trees), "| DP engine:", length(multi_trees) - length(prob_trees),
    "| probabilistic engine:", length(prob_trees), "(both sampled the same way)\n"
)

# Parse every path into (observation, stem label). Path labels get a "p"
# prefix so they cannot collide with the labels of spliced stems ("x").
recon_split <- strsplit(multi$recon, ";", fixed = TRUE)
parts <- data.table(
    treeID = rep(multi$treeID, lengths(recon_split)),
    path_idx = rep(multi$path_idx, lengths(recon_split)),
    pair = unlist(recon_split, use.names = FALSE)
)
rm(recon_split)
parts[, c("StemPaths", "lab") := tstrsplit(pair, ":", fixed = TRUE)]
parts[, `:=`(StemPaths = as.integer(StemPaths), lab = paste0("p", lab), pair = NULL)]

# ---- 4.2 Observation universe of the multi-path trees ----------------------
# Every measured observation of these trees in ALL censuses. Paths cover the
# DP window (at most up to the anchor); post-anchor observations (confirmed
# identities) are spliced in below with the exported links, so every interval
# and stock of a sampled tree is computed from one consistent set of stems.
# `oid` identifies an observation independently of StemPaths.
U <- stem_obs[
    treeID %in% multi_trees,
    .(treeID, xstem = stemID, CensusID, StemPaths, dbh, t, quadrat)
]
# Trees without a quadrat (no coordinates) are not in `rec`, so they have no
# observations here and nothing to sample; they are reported, not used.
no_obs_trees <- setdiff(multi_trees, U$treeID)
multi_trees <- setdiff(multi_trees, no_obs_trees)
multi <- multi[treeID %in% multi_trees]
cat("[BA] multi-path trees without analysed observations (no quadrat):", length(no_obs_trees), "\n")
U[, oid := .I]
setkey(U, treeID, StemPaths)
in_story <- unique(parts[, .(treeID, StemPaths)])[U, on = .(treeID, StemPaths), nomatch = 0L]$oid
U[, story := oid %in% in_story]

# Story window per tree: census range of its measured observations in paths.
win <- U[story == TRUE, .(wmin = min(CensusID), wmax = max(CensusID)), by = treeID]
U <- win[U, on = "treeID"]
U[, position := fcase(
    story == TRUE, "in story",
    is.na(wmin), "tree without story overlap",
    CensusID < wmin, "before window",
    CensusID > wmax, "after window",
    default = "inside window, missing"
)]

# Trees that cannot be spliced safely fall back to the exported reconstruction.
fallback_trees <- unique(U[position %in% c("inside window, missing", "tree without story overlap"), treeID])

# ---- 4.3 Splice: complete every path with the stage-2 links ----------------
# An observation outside the story window joins the path stem holding its
# exported stem-mate at the nearest window edge: the last in-story observation
# of the same exported stem (for later censuses) or the first one (for earlier
# censuses). Exported stems without any in-story observation keep their own
# label ("x" + exported stem).
edge_obs <- U[story == TRUE, .(
    first_in = oid[which.min(CensusID)],
    last_in = oid[which.max(CensusID)]
), by = .(treeID, xstem)]
outside <- U[position %in% c("before window", "after window") & !treeID %in% fallback_trees]
outside <- edge_obs[outside, on = .(treeID, xstem)]
outside[, anchor_oid := fifelse(position == "after window", last_in, first_in)]

# Story labels per path, keyed by observation id (only measured observations).
story_lab <- parts[U[story == TRUE & !treeID %in% fallback_trees, .(treeID, StemPaths, oid, CensusID, dbh, t, quadrat)],
    on = .(treeID, StemPaths), nomatch = 0L
][, .(treeID, path_idx, oid, lab, CensusID, dbh, t, quadrat)]

anchored <- story_lab[, .(treeID, path_idx, anchor_oid = oid, lab)][
    outside[!is.na(anchor_oid), .(treeID, oid, anchor_oid, CensusID, dbh, t, quadrat)],
    on = .(treeID, anchor_oid), allow.cartesian = TRUE, nomatch = 0L
][, .(treeID, path_idx, oid, lab, CensusID, dbh, t, quadrat)]

own_stem <- unique(multi[!treeID %in% fallback_trees, .(treeID, path_idx)])[
    outside[is.na(anchor_oid), .(treeID, oid, lab = paste0("x", xstem), CensusID, dbh, t, quadrat)],
    on = "treeID", allow.cartesian = TRUE, nomatch = 0L
]

comp <- rbindlist(list(story_lab, anchored, own_stem), use.names = TRUE)
rm(anchored, own_stem, edge_obs)

# A spliced path must never hold two observations of one stem in one census.
dup_trees <- comp[, .N, by = .(treeID, path_idx, lab, CensusID)][N > 1L, unique(treeID)]
if (length(dup_trees) > 0L) {
    fallback_trees <- union(fallback_trees, dup_trees)
    comp <- comp[!treeID %in% dup_trees]
}
sampled_trees <- setdiff(multi_trees, fallback_trees)

# Every completed path holds exactly the tree's measured observations.
n_expect <- U[treeID %in% sampled_trees, .N, by = treeID]
n_have <- comp[, .N, by = .(treeID, path_idx)]
cov_chk <- n_expect[n_have, on = "treeID"][N != i.N]
bio_check(
    nrow(cov_chk) == 0L &&
        nrow(n_have) == multi[treeID %in% sampled_trees, .N],
    "Every completed path contains exactly its tree's measured observations (all censuses)",
    examples = cov_chk$treeID, n_bad = uniqueN(cov_chk$treeID)
)
rm(n_expect, n_have, cov_chk)

cat(sprintf(
    "[BA] identity MC: %d trees sampled | %d fallback (exported reconstruction) | spliced observations: %d before / %d after the story window\n",
    length(sampled_trees), length(fallback_trees),
    U[position == "before window" & treeID %in% sampled_trees, .N],
    U[position == "after window" & treeID %in% sampled_trees, .N]
))

# ---- 4.4 Gap filling, per-path stocks and fluxes ----------------------------
# An exported stem alive (A) without a DBH outside its measured span stays
# alive in every path: the census joins the path stem that holds the exported
# stem's last measurement (after the span) or first one (before it), as in the
# splice. Censuses between two measurements are gap-filled anyway.
xspan <- U[treeID %in% sampled_trees, .(
    first_oid = oid[which.min(CensusID)], last_oid = oid[which.max(CensusID)],
    c_min = min(CensusID), c_max = max(CensusID)
), by = .(treeID, xstem)]
path_alive <- xspan[alive_nodbh[, .(treeID, xstem = stemID, CensusID)], on = .(treeID, xstem), nomatch = 0L][
    CensusID > c_max | CensusID < c_min
][, anchor_oid := fifelse(CensusID > c_max, last_oid, first_oid)]
path_alive <- unique(comp[, .(treeID, path_idx, anchor_oid = oid, lab)][
    path_alive[, .(treeID, anchor_oid, CensusID)],
    on = .(treeID, anchor_oid), allow.cartesian = TRUE, nomatch = 0L
][, .(treeID, path_idx, lab, CensusID)])
rm(xspan)

comp_f <- fill_stem_gaps(comp, c("treeID", "path_idx", "lab"), tree_dates, alive = path_alive)
comp_f[, BA := ba_m2(dbh)]

# Measured stock per tree and census must be identical in every path.
meas_path <- comp_f[interpolated == FALSE, .(ba = fsum(BA)), by = .(treeID, path_idx, CensusID)]
meas_exp <- U[treeID %in% sampled_trees, .(ba_exp = fsum(ba_m2(dbh))), by = .(treeID, CensusID)]
inv <- meas_exp[meas_path, on = .(treeID, CensusID)][abs(ba - ba_exp) > 1e-9]
bio_check(
    nrow(inv) == 0L,
    "Measured BA stock per tree and census is identical in every path (identity cannot change measurements)",
    examples = inv$treeID, n_bad = uniqueN(inv$treeID)
)
rm(meas_path, meas_exp, inv)

# All intervals and censuses of sampled trees come from the completed paths.
# After the anchor, every observation is linked by the (confirmed) exported
# links, so post-anchor values equal the exported ones except where a gap
# spanning the anchor is interpolated from a path-specific neighbour.
post_decomp <- decompose_intervals(
    comp_f, census_pairs,
    stem_cols = c("quadrat", "treeID", "path_idx", "lab"),
    tree_cols = c("quadrat", "treeID", "path_idx")
)
setkey(post_decomp, treeID, path_idx)
tree_stock <- comp_f[, .(TotalBA_m2 = fsum(BA)), by = .(quadrat, treeID, path_idx, CensusID)]
setkey(tree_stock, treeID, path_idx)
cat("[BA] per-path decompositions:", nrow(post_decomp), "rows | per-path stocks:", nrow(tree_stock), "rows\n")

# ---- 4.5 Fixed component (identical in every realization) ------------------
# Trees without identity uncertainty (single path, no posterior, fallback):
# exported reconstruction for every census and interval.
fixed_flux_q <- map_tree_change[
    !treeID %in% sampled_trees,
    lapply(.SD, fsum),
    .SDcols = flux_cols, by = .(quadrat, CensusID_from, CensusID_to)
]
fixed_stock_q <- tree_census[
    !treeID %in% sampled_trees,
    .(TotalBA_m2 = fsum(TotalBA_m2)),
    by = .(quadrat, CensusID)
]

# ---- 4.6 Realizations ------------------------------------------------------
# One path per sampled tree and realization, drawn by weight (inverse-CDF on
# the cumulative weights of the tree's paths).
wtab <- multi[treeID %in% sampled_trees, .(treeID, path_idx, w)]
setorder(wtab, treeID, path_idx)
wtab[, cw := cumsum(w), by = treeID]
wtab[, cw := cw / cw[.N], by = treeID]
wtab[, cw_prev := shift(cw, fill = 0), by = treeID]
tree_index <- wtab[, .GRP, by = treeID]
wtab[tree_index, on = "treeID", tidx := i.GRP]
n_sampled <- nrow(tree_index)

add_fixed <- function(var_dt, fixed_dt, by_cols, value_cols) {
    out <- rbindlist(list(fixed_dt, var_dt), use.names = TRUE, fill = TRUE)
    out[, lapply(.SD, fsum), .SDcols = value_cols, by = by_cols]
}

realization_dir <- file.path(out_dir, "ba_mc_realizations_treeID")
if (write_tree_realizations && !dir.exists(realization_dir)) dir.create(realization_dir, recursive = TRUE)

set.seed(mc_seed)
all_quadrat_realizations <- vector("list", K_realizations)
all_stock_realizations <- vector("list", K_realizations)
for (k in seq_len(K_realizations)) {
    if (k == 1L || k %% 100L == 0L) cat(sprintf("  Realization %d / %d\n", k, K_realizations))
    u <- runif(n_sampled)
    ur <- u[wtab$tidx]
    sel <- wtab[ur > cw_prev & ur <= cw, .(treeID, path_idx)]
    k_flux <- post_decomp[sel, on = .(treeID, path_idx), nomatch = 0L]
    k_flux_q <- k_flux[, lapply(.SD, fsum), .SDcols = flux_cols, by = .(quadrat, CensusID_from, CensusID_to)]
    all_quadrat_realizations[[k]] <- add_fixed(k_flux_q, fixed_flux_q, c("quadrat", "CensusID_from", "CensusID_to"), flux_cols)[, realization := k]
    k_stock_q <- tree_stock[sel, on = .(treeID, path_idx), nomatch = 0L][, .(TotalBA_m2 = fsum(TotalBA_m2)), by = .(quadrat, CensusID)]
    all_stock_realizations[[k]] <- add_fixed(k_stock_q, fixed_stock_q, c("quadrat", "CensusID"), "TotalBA_m2")[, realization := k]
    if (write_tree_realizations) {
        k_tree <- rbindlist(list(
            map_tree_change[!treeID %in% sampled_trees, c("quadrat", "treeID", "CensusID_from", "CensusID_to", flux_cols), with = FALSE],
            k_flux[, c("quadrat", "treeID", "CensusID_from", "CensusID_to", flux_cols), with = FALSE]
        ))[, realization := k]
        write_feather(k_tree, file.path(realization_dir, sprintf("ba_mc_realization_treeID_%04d.feather", k)))
    }
}
all_quadrat_realizations <- rbindlist(all_quadrat_realizations)
all_stock_realizations <- rbindlist(all_stock_realizations)
all_quadrat_realizations[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
cat("[BA] realizations:", K_realizations, "| quadrat-flux rows:", nrow(all_quadrat_realizations), "| stock rows:", nrow(all_stock_realizations), "\n")

# ---- 4.7 Accounting and invariance checks ----------------------------------
plot_flux_real <- all_quadrat_realizations[, lapply(.SD, fsum), .SDcols = c(flux_cols, "DeltaBA_total"), by = .(CensusID_from, CensusID_to, realization)]
plot_stock_real <- all_stock_realizations[, .(TotalBA_m2 = fsum(TotalBA_m2)), by = .(CensusID, realization)]
acc <- plot_stock_real[, .(CensusID_to = CensusID, realization, s_to = TotalBA_m2)][
    plot_stock_real[, .(CensusID_from = CensusID, realization, s_from = TotalBA_m2)][
        plot_flux_real,
        on = .(CensusID_from, realization)
    ],
    on = .(CensusID_to, realization)
]
acc_err <- acc[, max(abs(DeltaBA_total - (s_to - s_from)))]
bio_check(
    acc_err < 1e-6,
    sprintf("In every realization, Growth + Loss + Gain equals the stock change (max error %.2e m2)", acc_err)
)
map_plot_stock <- map_quadrat_stock[, .(map = fsum(TotalBA_m2)), by = CensusID]
stock_dev <- plot_stock_real[map_plot_stock, on = "CensusID"][, max(abs(TotalBA_m2 - map))]
cat(sprintf("[BA] largest |MC plot stock - exported stock| across realizations: %.4f m2 (only alive-but-unmeasured stems can differ)\n", stock_dev))

# ---- 4.8 Summaries and outputs ----------------------------------------------
all_quadrat_summary <- summarise_flux(all_quadrat_realizations, c("quadrat", "CensusID_from", "CensusID_to"))
all_stock_summary <- all_stock_realizations[,
    {
        q <- fquantile(TotalBA_m2, c(0.025, 0.975))
        .(TotalBA_mean = fmean(TotalBA_m2), TotalBA_sd = fsd(TotalBA_m2), TotalBA_lwr = q[1L], TotalBA_upr = q[2L])
    },
    by = .(quadrat, CensusID)
]

write_feather(map_tree_change, file.path(out_dir, "ba_map_change_treeID.feather"))
write_feather(map_quadrat_change, file.path(out_dir, "ba_map_change_quadrat.feather"))
write_feather(map_quadrat_stock, file.path(out_dir, "ba_map_stock_quadrat.feather"))
if (write_quadrat_realizations) {
    write_feather(all_quadrat_realizations, file.path(out_dir, "ba_mc_realizations_quadrat.feather"))
    write_feather(all_stock_realizations, file.path(out_dir, "ba_mc_realizations_stock_quadrat.feather"))
}
write_feather(all_quadrat_summary, file.path(out_dir, "ba_mc_summary_quadrat.feather"))
write_feather(all_stock_summary, file.path(out_dir, "ba_mc_summary_stock_quadrat.feather"))
cat("[BA] outputs written to", out_dir, "\n")

# ---- 4.9 Diagnostics --------------------------------------------------------
# Monte Carlo error of the plot-level 95 % interval bounds: bootstrap over
# realizations.
per_ha <- 1 / 50 # plot total (m2) -> per hectare
mc_quantile_se <- function(x, B = 500L) {
    qs <- replicate(B, fquantile(sample(x, replace = TRUE), c(0.025, 0.975)))
    c(se_lwr = sd(qs[1, ]), se_upr = sd(qs[2, ]), width = unname(diff(fquantile(x, c(0.025, 0.975)))))
}
flux_long_real <- melt(plot_flux_real, id.vars = c("CensusID_from", "CensusID_to", "realization"), variable.name = "Component")
mc_err <- flux_long_real[, as.list(mc_quantile_se(value)), by = .(CensusID_from, Component)]
mc_err[, rel_se := fifelse(width > 0, pmax(se_lwr, se_upr) / width, 0)]

# Posterior probability of the exported partition (label-invariant): the
# weight of the paths that group the tree's observations exactly as the R
# tables do.
sig_path <- comp[order(treeID, path_idx, oid)][, canon := match(lab, unique(lab)), by = .(treeID, path_idx)][
    , .(sig = paste(canon, collapse = ",")),
    by = .(treeID, path_idx)
]
sig_exp <- U[treeID %in% sampled_trees][order(treeID, oid)][, canon := match(xstem, unique(xstem)), by = treeID][
    , .(sig_exp = paste(canon, collapse = ",")),
    by = treeID
]
sig_path <- multi[, .(treeID, path_idx, w)][sig_path, on = .(treeID, path_idx)]
exp_prob <- sig_exp[sig_path, on = "treeID"][, .(p_exported = sum(w[sig == sig_exp]), p_best = max(w)), by = treeID]

splice_methods <- NULL
if (report_splice_methods) {
    splice_methods <- s2[U[position %in% c("before window", "after window") & treeID %in% sampled_trees], on = .(treeID, StemPaths), nomatch = NA][
        , .N,
        by = .(position, method)
    ][order(position, -N)]
}
rm(s2)

diag_line("")
diag_line("# ba_mc_diagnostics.txt — identity-uncertainty Monte Carlo, generated ", format(Sys.time()))
diag_line("K_realizations = ", K_realizations, " | seed = ", mc_seed, " | anchor census = ", ANCHOR_START_CENSUS)
diag_line("")
diag_line("## Trees")
diag_line("trees with a posterior: ", uniqueN(post_full$treeID))
diag_line("  single path (no identity uncertainty): ", post_full[n_paths == 1L, uniqueN(treeID)])
diag_line(
    "  multi path, sampled: ", length(sampled_trees),
    sprintf(" (of which probabilistic engine: %d)", length(intersect(sampled_trees, prob_trees)))
)
diag_line("  multi path, fallback to exported reconstruction: ", length(fallback_trees))
diag_line("  multi path, not analysed (no quadrat / no measured observation): ", length(no_obs_trees))
no_post_trees <- setdiff(unique(rec$treeID), post_full$treeID)
diag_line("analysed trees without a posterior (exported reconstruction in every realization): ", length(no_post_trees))
diag_line("  single-stem tags: ", sum(no_post_trees %in% single_stem_trees))
diag_line(
    "  multi-stem trees (no sampled paths, e.g. measured in only one census up to the anchor): ",
    sum(!no_post_trees %in% single_stem_trees)
)
diag_line("")
diag_line("## Measured observations of multi-path trees (all censuses)")
for (r in seq_len(nrow(U[, .N, by = position]))) {
    p <- U[, .N, by = position][r]
    diag_line(sprintf("  %-28s %8d", p$position, p$N))
}
if (!is.null(splice_methods)) {
    diag_line("spliced observations by stage-2 method:")
    for (r in seq_len(nrow(splice_methods))) {
        diag_line(sprintf("  %-14s %-28s %6d", splice_methods$position[r], splice_methods$method[r], splice_methods$N[r]))
    }
}
diag_line("")
diag_line("## Exported reconstruction vs posterior (sampled trees, by engine)")
diag_line("Both engines are sampled the same way (path_count / sum(path_count)).")
diag_line("DP draws are exact posterior samples and repeat; probabilistic draws are")
diag_line("approximate and nearly all unique, so a 'most probable path' exists only for DP trees.")
exp_prob[, engine := fifelse(treeID %in% prob_trees, "probabilistic", "DP")]
path_stats <- multi[treeID %in% sampled_trees, .(n_paths = .N, once = sum(path_count == 1L)), by = treeID]
path_stats[, engine := fifelse(treeID %in% prob_trees, "probabilistic", "DP")]
for (e in c("DP", "probabilistic")) {
    ep <- exp_prob[engine == e]
    ps <- path_stats[engine == e]
    if (nrow(ep) == 0L) next
    diag_line(sprintf(
        "  %-13s trees %6d | median unique paths per tree %5.0f | paths drawn once %5.1f %% | mean posterior probability of the exported partition %.3f | exported partition among the sampled paths %.3f%s",
        e, nrow(ep), median(ps$n_paths), 100 * sum(ps$once) / sum(ps$n_paths), ep[, mean(p_exported)], ep[, mean(p_exported > 0)],
        if (e == "DP") sprintf(" | most probable path %.3f", ep[, mean(p_exported >= p_best - 1e-12)]) else ""
    ))
}
rm(path_stats)
diag_line("")
diag_line("## Invariance")
diag_line(sprintf("largest |MC plot stock - exported stock|: %.6f m2 (alive-but-unmeasured stems only)", stock_dev))
diag_line(sprintf("largest accounting error (Growth + Loss + Gain vs stock change): %.2e m2", acc_err))
diag_line("")
diag_line("## Monte Carlo error of the plot-level 95 % interval bounds (m2 per plot)")
for (r in seq_len(nrow(mc_err))) {
    e <- mc_err[r]
    diag_line(sprintf(
        "  C%d->C%d %-14s width %9.4f | SE(lwr) %8.4f | SE(upr) %8.4f | max SE / width %.3f",
        e$CensusID_from, e$CensusID_from + 1L, e$Component, e$width, e$se_lwr, e$se_upr, e$rel_se
    ))
}
cat("[BA] diagnostics written to", diag_file, "\n")

# ============================================================
# SECTION 5: Figures (identity uncertainty only)
# ============================================================
# Exported reconstruction – bold solid line.
# Identity MC            – 95 % interval ribbon across realizations + dashed
#                          center line (mean or median, per mc_center).
# Forest-level values are plot totals divided by the plot area (50 ha).
# Post-anchor censuses/intervals have zero identity uncertainty by construction.
# ============================================================

pal <- c(MAP = "#1b7a56", MC = "#c4520a")
flux_pal <- c(
    Growth_BA     = "#1b7a56",
    Loss_BA       = "#c4520a",
    Gain_BA       = "#5b57a8",
    DeltaBA_total = "#c4186a"
)
COL_OBS <- "#3d3d3a" # charcoal → exported reconstruction
COL_MOD <- "#2e9e75" # green    → posterior paths

theme_forest <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
        theme(
            panel.grid.major = element_line(colour = "#e8e5e0", linewidth = 0.35),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "#c8c4bc", fill = NA, linewidth = 0.5),
            panel.spacing = unit(0.8, "lines"),
            axis.title = element_text(size = rel(0.85), colour = "#555550"),
            axis.text = element_text(size = rel(0.78), colour = "#777770"),
            axis.ticks = element_line(colour = "#c8c4bc", linewidth = 0.3),
            strip.background = element_rect(fill = "#f2efe9", colour = "#c8c4bc", linewidth = 0.4),
            strip.text = element_text(
                size = rel(0.80), colour = "#444440",
                face = "bold", margin = margin(3, 6, 3, 6)
            ),
            legend.position = "top",
            legend.key.size = unit(0.85, "lines"),
            legend.text = element_text(size = rel(0.80), colour = "#555550"),
            legend.title = element_text(size = rel(0.82), colour = "#333330", face = "bold"),
            legend.background = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(
                size = rel(1.10), face = "bold", colour = "#222220",
                margin = margin(b = 3)
            ),
            plot.subtitle = element_text(
                size = rel(0.83), colour = "#777770",
                margin = margin(b = 8)
            ),
            plot.caption = element_text(
                size = rel(0.70), colour = "#aaaaaa",
                hjust = 1, margin = margin(t = 6)
            ),
            plot.margin = margin(10, 12, 10, 10)
        )
}

center_fun <- function(x) {
    if (mc_center == "median") fmedian(x, na.rm = TRUE) else fmean(x, na.rm = TRUE)
}

# ── Figure 1: BA stock per hectare ──────────────────────────────────────────────
stock_exp <- map_plot_stock[, .(CensusID, value = map * per_ha)][dates[, .(CensusID, Year)], on = "CensusID"]
stock_mc <- plot_stock_real[, .(value = TotalBA_m2 * per_ha), by = .(CensusID, realization)][,
    .(center = center_fun(value), lwr = fquantile(value, 0.025), upr = fquantile(value, 0.975)),
    by = CensusID
][dates[, .(CensusID, Year)], on = "CensusID"]
stock_exp <- stock_exp[CensusID >= first_plot_census]
stock_mc <- stock_mc[CensusID >= first_plot_census]

fig1 <- ggplot() +
    geom_ribbon(data = stock_mc, aes(Year, ymin = lwr, ymax = upr, fill = "MC"), alpha = 0.25) +
    geom_line(data = stock_mc, aes(Year, center, colour = "MC"), linewidth = 0.9, linetype = "dashed") +
    geom_line(data = stock_exp, aes(Year, value, colour = "MAP"), linewidth = 1.4) +
    scale_colour_manual(
        "Estimate",
        values = pal,
        labels = c(MAP = "Exported reconstruction", MC = "Identity MC (95 % interval)")
    ) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
        title = "Forest-level BA stock per hectare",
        subtitle = sprintf(
            "Identity uncertainty only (%d realizations) · stock can vary only through alive-but-unmeasured stems",
            K_realizations
        ),
        x = "Year",
        y = expression("BA (m"^2 ~ "ha"^
            {
                -1
            } * ")")
    ) +
    theme_forest()
print(fig1)

# ── Figure 2: annual BA flux components per hectare ────────────────────────────
flux_labels <- c(
    Growth_BA     = "Growth",
    Loss_BA       = "Loss (mortality)",
    Gain_BA       = "Gain (recruitment + ingrowth)",
    DeltaBA_total = "Net BA change"
)
flux_all <- c(flux_cols, "DeltaBA_total")

flux_exp <- map_quadrat_change[, lapply(.SD, fsum), .SDcols = flux_all, by = .(CensusID_from, CensusID_to)]
flux_exp <- census_pairs[, .(CensusID_from, Interval_yr, Year_mid)][flux_exp, on = "CensusID_from"]
flux_exp[, (flux_all) := lapply(.SD, function(x) x * per_ha / Interval_yr), .SDcols = flux_all]
flux_exp_long <- melt(flux_exp[CensusID_from >= first_plot_census],
    id.vars = c("CensusID_from", "Year_mid"), measure.vars = flux_all,
    variable.name = "Component", value.name = "value"
)

flux_mc <- census_pairs[, .(CensusID_from, Interval_yr, Year_mid)][plot_flux_real, on = "CensusID_from"]
flux_mc[, (flux_all) := lapply(.SD, function(x) x * per_ha / Interval_yr), .SDcols = flux_all]
flux_mc_ci <- melt(flux_mc[CensusID_from >= first_plot_census],
    id.vars = c("CensusID_from", "Year_mid", "realization"), measure.vars = flux_all,
    variable.name = "Component", value.name = "value"
)[, .(center = center_fun(value), lwr = fquantile(value, 0.025), upr = fquantile(value, 0.975)),
    by = .(CensusID_from, Year_mid, Component)
]

fig2 <- ggplot() +
    geom_ribbon(data = flux_mc_ci, aes(Year_mid, ymin = lwr, ymax = upr, fill = Component), alpha = 0.6) +
    geom_line(data = flux_mc_ci, aes(Year_mid, center, colour = Component), linewidth = 0.7, linetype = "dashed") +
    geom_line(data = flux_exp_long, aes(Year_mid, value), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "#bbbbaa", linewidth = 0.4) +
    facet_wrap(~Component, scales = "free_y", ncol = 1, labeller = as_labeller(flux_labels)) +
    scale_colour_manual(values = flux_pal, guide = "none") +
    scale_fill_manual(values = flux_pal, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
        title = "Forest-level annual BA fluxes per hectare",
        subtitle = sprintf(
            "Bold = exported reconstruction · dashed %s + ribbon = identity MC 95 %% interval (%d realizations)",
            mc_center, K_realizations
        ),
        x = "Midpoint year between censuses",
        y = expression("BA flux (m"^2 ~ "ha"^{
            -1
        } ~ "yr"^
            {
                -1
            } * ")")
    ) +
    theme_forest()
print(fig2)

# ── Figure 3: BA trajectories of focal trees — one page per tree ────────────────
# Each page shows the exported reconstruction (panel 0, charcoal) and the most
# probable posterior paths of one tree (completed with the splice).
build_tree_page <- function(tid, paths_dt, exported_dt, dates_dt, first_census, weights, n_show = 5L) {
    top <- weights[treeID == tid][order(-w)][seq_len(min(n_show, .N))]
    ap <- paths_dt[treeID == tid & path_idx %in% top$path_idx, .(path_idx, stemID = lab, CensusID, BA)]
    if (nrow(ap) == 0L) {
        message(sprintf("[Fig3] treeID %s not found among sampled paths — skipping.", tid))
        return(NULL)
    }
    obs <- exported_dt[treeID == tid, .(path_idx = 0L, stemID = as.character(stemID), CensusID, BA)]
    pp <- rbindlist(list(ap, obs), use.names = TRUE)
    pp <- dates_dt[, .(CensusID, Year)][pp, on = "CensusID"][CensusID >= first_census]
    pp[, path_type := fifelse(path_idx == 0L, "Observed", "Modelled")]
    lab_w <- setNames(sprintf("Path %d (p = %.2f)", top$path_idx, top$w), top$path_idx)
    path_labeller <- labeller(path_idx = function(x) ifelse(x == "0", "Exported reconstruction", lab_w[x]))
    # Lines only for stems seen in >= 2 censuses (single points need none).
    # (filtering on a count keeps every column even when no stem qualifies)
    pp[, n_obs := .N, by = .(path_idx, stemID)]
    pp_lines <- pp[n_obs > 1L]
    ggplot(pp, aes(x = Year, y = BA, group = stemID, colour = path_type)) +
        geom_line(data = pp_lines, linewidth = 0.8) +
        geom_point(size = 1.6, shape = 21, fill = "white", stroke = 0.7) +
        scale_colour_manual(values = c(Observed = COL_OBS, Modelled = COL_MOD), guide = "none") +
        scale_x_continuous(breaks = scales::pretty_breaks(4)) +
        scale_y_continuous(labels = scales::label_comma()) +
        facet_wrap(~path_idx, scales = "free_y", ncol = 2, labeller = path_labeller) +
        labs(
            title = paste0("BA trajectories · treeID ", tid),
            subtitle = "One line per stem · posterior paths ordered by probability",
            x = "Year", y = expression("BA (m"^2 * ")")
        ) +
        theme_forest()
}

path_counts_tree <- multi[treeID %in% sampled_trees, .N, by = treeID]
set.seed(mc_seed)
focal_trees <- sample(path_counts_tree[N >= 5L]$treeID, 5L)

fig3_path <- file.path(out_dir, "fig3_BA_trajectories.pdf")
pdf(fig3_path, width = 9, height = 8)
for (tid in focal_trees) {
    pg <- build_tree_page(
        tid = tid,
        paths_dt = comp_f,
        exported_dt = exp_stems,
        dates_dt = dates,
        first_census = first_plot_census,
        weights = multi[, .(treeID, path_idx, w)]
    )
    if (!is.null(pg)) print(pg)
}
dev.off()
cat("[Fig3] Written:", fig3_path, "\n")

ggsave(file.path(out_dir, "fig1_BA_stock.pdf"), fig1, width = 8, height = 4.5)
ggsave(file.path(out_dir, "fig2_BA_fluxes.pdf"), fig2, width = 8, height = 10)
cat("[Figs] fig1_BA_stock.pdf and fig2_BA_fluxes.pdf saved to", out_dir, "\n")
