# ==============================================================================
# BCI Stem Reconstruction — Basal Area Uncertainty Propagation
# ==============================================================================
#
# PURPOSE
# -------
# Propagates stem-identity uncertainty (encoded in posterior reconstruction
# paths produced by the dp_global engine) into forest-level basal area (BA)
# stocks and fluxes for the full BCI 50-ha plot across nine stem censuses
# (1985–2022/3).
#
# DESIGN OVERVIEW
# ---------------
# For each tree the engine outputs one or more *reconstruction paths*: ordered
# sequences of StemID assignments across censuses. Trees with only one path are
# deterministic (MAP-equivalent); trees with multiple paths have genuine
# identity ambiguity. Each Monte Carlo (MC) realization draws one path per
# ambiguous tree proportional to path probability, then aggregates BA across all
# trees and all quadrats to produce a single forest-level estimate. Repeating
# this K_realizations times yields an empirical posterior distribution of
# forest-level BA stocks and fluxes from which we read uncertainty (95 % CI).
#
# ANCHOR CENSUSES AND SCOPE
# -------------------------
# The DP engine reconstructs identities through the anchor census,
# ANCHOR_START_CENSUS = Census 7. Post-anchor censuses (C8+) have confirmed
# stemID and are appended as deterministic rows. MC uncertainty therefore
# applies only to pre-anchor intervals (C1–C6).
# See the ANCHOR_START_CENSUS block in Section 1 for full documentation.
#
# OUTPUTS (written to BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs/)
# --------
#   ba_map_change_treeID.feather        MAP tree-level flux per census pair
#   ba_map_change_quadrat.feather       MAP quadrat-level flux per census pair
#   ba_map_stock_quadrat.feather        MAP quadrat-level BA stock per census
#   ba_mc_realizations_quadrat.feather  K MC realization quadrat fluxes
#   ba_mc_realizations_stock_*.feather  K MC realization quadrat stocks
#   ba_mc_summary_quadrat.feather       Empirical 95 % CI of quadrat fluxes
#   ba_mc_summary_stock_quadrat.feather Empirical 95 % CI of quadrat stocks
#   ba_mc_realizations_treeID/          Per-realization tree-level feather files
#   fig1_stock.pdf                      Forest BA stock: MAP vs MC
#   fig2_fluxes.pdf                     Forest BA fluxes: MAP vs MC
#   fig3_trajectories.pdf               Individual-tree BA trajectories (6 trees)
#
# NOTES
# --------
# Additional methodological decisions and design choices are documented in
# biomass_stocks_fluxes.R and should be evaluated against individual study
# needs. In the present script, giant strangler ficus (> 500 mm DBH) are
# retained despite violating standard allometric assumptions and
# disproportionately inflating plot-level basal area. No correction is applied
# for the DBH bias introduced by measuring around buttresses in the first
# census, nor for buttressed trees in subsequent censuses. Bias-corrected
# productivity (G*) and mortality (M*) estimators from Kohyama et al. (2019) are
# not included. Moreover, we are not accounting for the lack of diameter growth
# in some palm species. For a valid output, these decisions should be evaluated.
# ==============================================================================

rm(list = ls())

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(collapse)
library(arrow)

# ============================================================
# SECTION 1: Configuration and data loading
# ============================================================
# Reads all 9 BCI stem census files, stacks them into a single
# long data.table (rec), and imputes missing measurement dates
# via a two-step median strategy: first within quadrat × census,
# then within the full census, so every row has an ExactDate.
# ============================================================

bci_stem_nums <- as.character(1:9)
census_list <- lapply(bci_stem_nums, function(num) {
    fp <- paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", num, ".Rdata")
    if (!file.exists(fp)) stop("Missing census file: ", fp)
    load(fp)
    get(paste0("bci.stem", num))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)
rec <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rec <- rec[!is.na(quadrat)]
rm(census_list, bci_stem_nums)

# Impute missing ExactDate: quadrat-census median, then plot-census median.
rec[, date_quad_census := median(ExactDate, na.rm = TRUE), by = .(quadrat, CensusID)]
rec[, date_plot_census := median(ExactDate, na.rm = TRUE), by = CensusID]
rec[, ExactDate := fifelse(is.na(ExactDate), date_quad_census, ExactDate)]
rec[is.na(ExactDate), ExactDate := date_plot_census]
rec[, c("date_quad_census", "date_plot_census") := NULL]

post_file <- "./BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds"

# ── Anchor census ──────────────────────────────────────────────────────────────
# ANCHOR_START_CENSUS is the first census where individual stem identities are
# known with certainty (from 2010 onward = Census 7). The DP reconstruction
# algorithm runs *backward* from this anchor: it finds the most probable
# stem-identity assignment for censuses 1 through ANCHOR_START_CENSUS, treating
# the anchor census as the fixed endpoint.
#
# Consequence for uncertainty propagation:
#   • Posterior paths only encode identity choices for StemPaths in C1–C7
#     (the DP segment). The 'recon' string in the paths file does NOT include
#     StemPaths for C8 or C9.
#   • Fluxes and stocks involving post-anchor censuses (C7→C8, C8→C9) are
#     deterministic: every MC realization must use the MAP (stemID-based)
#     values for those intervals.
#   • Comparison figures must be restricted to pre-anchor intervals; including
#     post-anchor intervals in the MC layers would show inflated mortality
#     (C7 stems appear dead because C8/C9 StemPaths are absent from
#     all_paths) and zero recruitment — artifacts, not uncertainty.
ANCHOR_START_CENSUS <- 7L

# Exclude the first census from all figures and plot summaries.
# The first census is usually omitted because buttressed trees were measured
# around the buttress at breast height, which introduces a strong DBH bias.
first_plot_census <- 2L

# Number of MC realizations to draw from the posterior distribution of paths.
K_realizations <- 200L

mc_center <- "mean" # choose "mean" or "median"
mc_center <- match.arg(mc_center, c("mean", "median"))

# ============================================================
# SECTION 2: Helper functions
# ============================================================
# ba_m2(dbh_mm)   – converts stem DBH (mm) to basal area in m²
#                   using the standard circle formula: BA = π/4 × (d/1000)².
#
# parse_recon(s)  – parses a posterior 'recon' string of the form
#                   "obs_row_id:ReconstructedStemID;..." into a two-column
#                   data.table (StemPaths, stemID). obs_row_id maps back to
#                   the original census row in rec via the StemPaths column.
#
# decompose_ba()  – given a merged from/to BA table (one row per stem per
#                   interval), classifies stems as survivor (observed in both
#                   censuses), death (present only in 'from'), or recruit
#                   (present only in 'to'), then returns the signed BA
#                   components: Growth (positive), Loss (negative), Gain
#                   (positive). Used identically for MAP and posterior loops.
#
# summarise_flux()– collapses the K-realization flux distribution to empirical
#                   mean/sd/95 % CI for each grouping level (quadrat or plot).
# ============================================================

ba_m2 <- function(dbh_mm) pi / 4 * (dbh_mm / 1000)^2

parse_recon <- function(recon_str) {
    # recon_str: a single character string e.g. "3:101;3:205;7:88"
    pairs <- strsplit(recon_str, ";", fixed = TRUE)[[1L]]
    parts <- strsplit(pairs, ":", fixed = TRUE)
    data.table(
        StemPaths = as.integer(vapply(parts, `[`, character(1L), 1L)),
        stemID    = as.integer(vapply(parts, `[`, character(1L), 2L))
    )
}

decompose_ba <- function(m, by_cols) {
    m[, status := fifelse(
        !is.na(BA_from) & !is.na(BA_to), "survivor",
        fifelse(!is.na(BA_from), "death", "recruit")
    )]
    m[, .(
        Growth_BA     = fsum((BA_to - BA_from) * (status == "survivor"), na.rm = TRUE),
        Loss_BA       = -fsum(BA_from * (status == "death"), na.rm = TRUE),
        Gain_BA       = fsum(BA_to * (status == "recruit"), na.rm = TRUE)
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
                Gain_lwr    = qa[1L],               Gain_upr   = qa[2L]
            )
        },
        by = by_cols
    ]
}

# ============================================================
# SECTION 3: MAP baseline BA decomposition
# ============================================================
# Uses observed (maximum a posteriori) DBH measurements — i.e. the single
# best-guess stem-identity assignment — to compute deterministic BA stocks
# and fluxes for every census and census interval. This is the "ground truth"
# reference against which MC uncertainty is assessed.
#
# Stems are matched across consecutive censuses on (quadrat, treeID, stemID).
# The decompose_ba() function then classifies each matched/unmatched stem as
# survivor / death / recruit and accumulates the signed BA components:
#   Growth_BA  – BA increment of stems alive in both censuses (positive)
#   Loss_BA    – BA of stems that died (negative)
#   Gain_BA    – BA of newly recruited stems (positive)
#
# Outputs:
#   map_tree_change    – tree-level flux per census pair
#   map_quadrat_change – quadrat-level flux (sum within each 20×20 m quadrat)
#   map_quadrat_stock  – quadrat-level BA stock (sum of stem BAs per census)
#
# These serve two roles downstream:
#   (a) Fixed component in every MC realization for trees with a single path.
#   (b) Reference lines in the comparison figures.
# ============================================================

tree_census <- rec[!is.na(dbh) & !is.na(stemID) & !is.na(treeID),
    .(TotalBA_m2 = fsum(ba_m2(dbh)), NumStems = .N),
    by = .(quadrat, treeID, CensusID)
]
setorder(tree_census, treeID, CensusID)
cat("[BA] treeID x census rows:", nrow(tree_census), "\n")

dates <- rec[!is.na(treeID), .(Date = median(ExactDate, na.rm = TRUE)), by = CensusID][order(CensusID)]
dates[, Year := as.integer(format(Date, "%Y"))]
cat("Need >= 2 censuses" = nrow(dates) >= 2L)

census_pairs <- data.table(
    CensusID_from = dates$CensusID[-nrow(dates)],
    CensusID_to   = dates$CensusID[-1L],
    Date_from     = dates$Date[-nrow(dates)],
    Date_to       = dates$Date[-1L]
)
census_pairs[, Date_mid := Date_from + (Date_to - Date_from) / 2]
census_pairs[, Year_mid := as.integer(format(Date_mid, "%Y"))]

stem_dt <- rec[
    !is.na(treeID) & !is.na(stemID),
    .(treeID, stemID = stemID, CensusID, BA = ba_m2(dbh), quadrat)
]

# MAP census-pair decomposition using shared decompose_ba().
map_change <- rbindlist(lapply(seq_len(nrow(census_pairs)), function(i) {
    cf <- census_pairs$CensusID_from[i]
    ct <- census_pairs$CensusID_to[i]
    sf <- stem_dt[CensusID == cf, .(quadrat, treeID, stemID, BA_from = BA)]
    st <- stem_dt[CensusID == ct, .(quadrat, treeID, stemID, BA_to = BA)]
    d <- decompose_ba(
        merge(sf, st, by = c("quadrat", "treeID", "stemID"), all = TRUE),
        by_cols = c("quadrat", "treeID")
    )
    d[, `:=`(
        CensusID_from = cf, CensusID_to = ct,
        Date_from = census_pairs$Date_from[i],
        Date_to = census_pairs$Date_to[i]
    )]
}))
cat(
    "[BA] MAP decomposition:", nrow(map_change), "intervals across",
    uniqueN(map_change$treeID), "treeIDs\n"
)

flux_cols <- c("Growth_BA", "Loss_BA", "Gain_BA")

# Quadrat-level MAP aggregations (deterministic; no CI needed here).
map_tree_change <- copy(map_change)

# Fluxes: 8 census intervals
map_quadrat_change <- map_tree_change[,
    lapply(.SD, fsum, na.rm = TRUE),
    .SDcols = flux_cols,
    by = .(quadrat, CensusID_from, CensusID_to)
]
map_quadrat_change[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]

# Stocks: 9 censuses
map_quadrat_stock <- tree_census[,
    .(TotalBA_m2 = fsum(TotalBA_m2, na.rm = TRUE)),
    by = .(quadrat, CensusID)
]
cat("[BA] MAP quadrat stock:", nrow(map_quadrat_stock), "quadrat×census rows\n")

# ============================================================
# SECTION 4: Posterior uncertainty propagation (MC realizations)
# ============================================================
# Reads the posterior reconstruction paths and propagates stem-identity
# uncertainty into forest-level BA stocks and fluxes.
#
# KEY CONCEPT — pre-anchor vs post-anchor censuses
# ─────────────────────────────────────────────────
# Posterior 'recon' strings only encode identity assignments for the
# pre-anchor DP segment (C1 through ANCHOR_START_CENSUS). Post-anchor
# rows (C8, C9, …) carry confirmed stemID and are NOT included in
# any posterior path. Therefore:
#   • all_paths is restricted to pre-anchor observations after the merge.
#   • post_decomp covers only census pairs where CensusID_to <=
#     ANCHOR_START_CENSUS (true sampling uncertainty).
#   • For post-anchor intervals the MAP flux (deterministic) is used
#     in every realization — post_decomp_postanchor handles this.
#
# TREE PARTITIONING
# ─────────────────
# Trees are split into two groups based on the number of posterior paths:
#   fixed component  – single-path trees. Their identity is unambiguous;
#                      MAP fluxes/stocks are contributed unchanged to every
#                      realization. (fixed_tree_change, fixed_stock)
#   variable component – multi-path trees. One path is drawn per realization
#                      proportionally to path_prob. (post_decomp, tree_stock)
#
# MC REALIZATION LOOP
# ────────────────────
# For each realization k = 1…K_realizations:
#   1. Sample one path per uncertain tree group (pre-drawn in sampled_paths).
#   2. Assemble per-tree flux table: fixed MAP trees + sampled uncertain trees
#      (pre-anchor from post_decomp) + deterministic post-anchor uncertain
#      trees (from post_decomp_postanchor).
#   3. Aggregate to quadrat level and store.
#   4. Repeat for stocks: fixed MAP + sampled pre-anchor + deterministic
#      post-anchor (from tree_census).
#
# Outputs:
#   all_quadrat_realizations – K × quadrat × interval flux table
#   all_stock_realizations   – K × quadrat × census stock table
#   (collapsed to empirical 95 % CI in all_quadrat_summary / all_stock_summary)
# ============================================================

post_full <- as.data.table(readRDS(post_file))
post_full[, group_id := treeID]
cat(
    "[BA] Posterior paths loaded:", nrow(post_full), "rows for",
    uniqueN(post_full$treeID), "treeIDs\n"
)

# Add path index if not present in the posterior file.
if (!"path_idx" %in% names(post_full)) {
    post_full[, path_idx := seq_len(.N), by = group_id]
}

# Guard against missing path_prob; normalise to sum = 1 within each group so
# that sample(prob = ...) is well-defined even when the engine returns raw
# log-probabilities or un-normalised scores.
post_full[, path_prob := fifelse(is.na(path_prob), 0, path_prob)]
post_full[, path_prob := {
    s <- fsum(path_prob)
    if (s > 0) path_prob / s else rep(1 / .N, .N)
}, by = group_id]

# ── Parse reconstruction strings → long stem table ────────────────────────────
# Each row in post_full holds a 'recon' string of the form
# "StemPaths:ReconstructedStemID;StemPaths:ReconstructedStemID;..." encoding the
# identity assignment for every stem observation in the pre-anchor DP segment.
# parse_recon() expands one string into one row per stem, keyed by StemPaths.
# The resulting all_paths table has one row per treeID × path × stem observation
# before the DBH merge.
all_paths <- rbindlist(mapply(
    function(grp, tid, pidx, pprob, recon_str) {
        dt <- parse_recon(recon_str)
        dt[, `:=`(group_id = grp, treeID = tid, path_idx = pidx, path_prob = pprob)]
        dt
    },
    post_full$group_id, post_full$treeID, post_full$path_idx,
    post_full$path_prob, post_full$recon,
    SIMPLIFY = FALSE
))

# Create compound stemID: "treeID_ReconstructedStemID". Used as a within-tree
# stem identifier when merging across censuses (two ReconstructedStemIDs from
# different trees could collide numerically; prefixing with treeID avoids this).
all_paths[, stemID := paste(treeID, stemID, sep = "_")]

# ── Merge with census observations to obtain DBH → BA per path × census ──────
# StemPaths in all_paths is the integer key used inside the 'recon' string. It maps
# to rec$StemPaths which indexes the original census row. Join key: (treeID,
# StemPaths). Do NOT join on stemID — all_paths$stemID is the compound stemID
# ("treeID_stemID"); they are in different namespaces.
obs_lookup <- rec[
    !is.na(treeID) & !is.na(StemPaths),
    .(quadrat, treeID, StemPaths, CensusID, dbh)
]
all_paths <- merge(all_paths, obs_lookup, by = c("treeID", "StemPaths"), all.x = TRUE)
all_paths <- all_paths[!is.na(dbh)] # drop path rows with no observed DBH
all_paths[, `:=`(BA = ba_m2(dbh), dbh = NULL)]
rm(obs_lookup)
gc()
setkey(all_paths, treeID)
cat(
    "[BA] Parsed path observations:", nrow(all_paths), "rows for",
    uniqueN(all_paths$treeID), "treeIDs\n"
)

# ── Posterior census-pair decomposition (PRE-ANCHOR INTERVALS ONLY) ───────────
# Posterior paths encode identity assignments only for censuses up to and
# including ANCHOR_START_CENSUS (the DP segment). Running decompose_ba() over
# post-anchor intervals would find zero rows in 'st' for C8/C9 (no obs_row_ids
# in all_paths for those censuses), classifying every C7 stem as dead and every
# C8 recruit from scratch — a pure artifact, not uncertainty.
# Solution: restrict to census pairs where CensusID_to <= ANCHOR_START_CENSUS.
pre_anchor_pairs <- census_pairs[CensusID_to <= ANCHOR_START_CENSUS]
post_decomp <- rbindlist(lapply(seq_len(nrow(pre_anchor_pairs)), function(i) {
    cf <- pre_anchor_pairs$CensusID_from[i]
    ct <- pre_anchor_pairs$CensusID_to[i]
    sf <- all_paths[CensusID == cf, .(group_id, quadrat, treeID, path_idx, stemID, BA_from = BA)]
    st <- all_paths[CensusID == ct, .(group_id, quadrat, treeID, path_idx, stemID, BA_to = BA)]
    d <- decompose_ba(
        merge(sf, st, by = c("group_id", "quadrat", "treeID", "path_idx", "stemID"), all = TRUE),
        by_cols = c("group_id", "quadrat", "treeID", "path_idx")
    )
    d[, `:=`(CensusID_from = cf, CensusID_to = ct)]
}))
post_decomp[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
cat(
    "[BA] Posterior decompositions:", nrow(post_decomp),
    "rows (pre-anchor intervals C1–C", ANCHOR_START_CENSUS, " only)\n"
)

# ── Per-path tree-level BA stocks (PRE-ANCHOR ONLY) ───────────────────────────
# Restrict to pre-anchor censuses. Post-anchor census stocks for uncertain trees
# are deterministic and will be pulled from tree_census in the realization loop.
tree_stock <- all_paths[CensusID <= ANCHOR_START_CENSUS,
    .(TotalBA_m2 = fsum(BA, na.rm = TRUE)),
    by = .(group_id, quadrat, treeID, path_idx, CensusID)
]
cat("[BA] Per-path tree stocks:", nrow(tree_stock), "rows (pre-anchor censuses)\n")

# ── Partition trees: single-path (fixed) vs multi-path (uncertain) ────────────
# fixed component : trees with one path only → MAP values used unchanged.
# variable component: trees with ≥ 2 paths → one path sampled per realization.
path_group_sizes <- post_full[, .N, by = group_id]
multi_path_groups <- path_group_sizes[N > 1L, group_id]
uncertain_treeIDs <- unique(post_full[group_id %in% multi_path_groups, treeID])
fixed_tree_change <- map_change[!treeID %in% uncertain_treeIDs]
fixed_tree_change[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
fixed_stock <- tree_census[!treeID %in% uncertain_treeIDs]

cat(
    "[BA] Uncertain groups:", length(multi_path_groups),
    "; uncertain treeIDs:", length(uncertain_treeIDs),
    "; fixed treeIDs:", uniqueN(fixed_tree_change$treeID), "\n"
)

# ── Post-anchor intervals for uncertain trees (deterministic) ─────────────────
# Censuses after ANCHOR_START_CENSUS have confirmed stemID — every path
# for these trees gives the same stem assignment. Their fluxes therefore equal
# MAP fluxes exactly. We pull MAP values for UNCERTAIN trees only from
# map_change (fixed trees are already fully covered by fixed_tree_change).
post_anchor_pairs <- census_pairs[CensusID_from >= ANCHOR_START_CENSUS]
if (nrow(post_anchor_pairs) > 0L) {
    post_decomp_postanchor <- map_change[
        treeID %in% uncertain_treeIDs &
            CensusID_from %in% post_anchor_pairs$CensusID_from
    ]
    post_decomp_postanchor[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
    cat(
        "[BA] Post-anchor intervals stored for", nrow(post_anchor_pairs),
        "census pair(s) (uncertain trees only, deterministic MAP values)\n"
    )
} else {
    post_decomp_postanchor <- data.table()
}

# Column sets for consistent subset-and-bind inside the realization loop.
tree_fixed_cols <- c("quadrat", "treeID", "CensusID_from", "CensusID_to", flux_cols)
stock_cols <- c("quadrat", "treeID", "CensusID", "TotalBA_m2")

out_dir <- "./BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs"
realization_dir <- file.path(out_dir, "ba_mc_realizations_treeID")
if (!dir.exists(realization_dir)) dir.create(realization_dir, recursive = TRUE)

all_quadrat_realizations <- vector("list", K_realizations)
all_stock_realizations <- vector("list", K_realizations)
k_trees_NtreeID <- vector("list", K_realizations)
k_stock_NtreeID <- vector("list", K_realizations)
# Pre-draw all path selections so each group gets exactly one path
# per realization, sampled proportionally to path_prob.
path_opts <- post_full[group_id %in% multi_path_groups, .(group_id, path_idx, path_prob)]
sampled_paths <- path_opts[,
    .(path_idx = sample(path_idx, K_realizations, replace = TRUE, prob = path_prob)),
    by = group_id
]
sampled_paths[, realization := seq_len(.N), by = group_id]
setkey(sampled_paths, group_id, realization)

for (k in seq_len(K_realizations)) {
    if (k == 1L || k %% 10L == 0L) cat(sprintf("  Realization %d / %d\n", k, K_realizations))
    # Select the sampled path for each uncertain group in this realization.
    sel <- sampled_paths[.(k), on = "realization", .(group_id, path_idx)]
    # ── Pre-anchor sampled fluxes and stocks ──────────────────────────────────
    k_var <- post_decomp[sel, on = .(group_id, path_idx), nomatch = 0L]
    k_var_stk <- tree_stock[sel, on = .(group_id, path_idx), nomatch = 0L]
    # ── Tree-level flux table for this realization ────────────────────────────
    # Combines three components:
    #   (1) fixed_tree_change : single-path trees — all census intervals, MAP.
    #   (2) k_var             : uncertain trees, pre-anchor intervals, sampled path.
    #   (3) post-anchor rows  : uncertain trees, post-anchor intervals, MAP
    #       (post_decomp_postanchor). These are deterministic because stemID
    #       is confirmed; every realization contributes the same MAP flux here.
    k_postanchor_flux <- if (nrow(post_decomp_postanchor) > 0L) {
        post_decomp_postanchor[, ..tree_fixed_cols]
    } else {
        NULL
    }
    k_tree <- rbindlist(
        c(
            list(
                fixed_tree_change[, ..tree_fixed_cols],
                k_var[, ..tree_fixed_cols]
            ),
            if (!is.null(k_postanchor_flux)) list(k_postanchor_flux) else list()
        ),
        use.names = TRUE, fill = TRUE
    )
    k_tree[, realization := k]
    write_feather(
        k_tree,
        file.path(realization_dir, sprintf("ba_mc_realization_treeID_%03d.feather", k))
    )
    # Quadrat-level: aggregate tree fluxes within each quadrat × interval.
    k_trees_NtreeID[[k]] <- k_tree[,
        .(N_unique_treeid = uniqueN(treeID)),
        by = .(CensusID_from, CensusID_to)
    ]
    all_quadrat_realizations[[k]] <- k_tree[,
        lapply(.SD, fsum, na.rm = TRUE),
        .SDcols = flux_cols,
        by = .(quadrat, CensusID_from, CensusID_to)
    ][, realization := k]
    rm(k_tree)
    gc()
    # ── Stock-level table for this realization ────────────────────────────────
    # (1) fixed_stock    : single-path trees, all censuses, MAP.
    # (2) k_var_stk      : uncertain trees, pre-anchor censuses, sampled path.
    # (3) post-anchor stocks: uncertain trees, post-anchor censuses (C8, C9, …).
    #     These come from tree_census (MAP) because stemID is confirmed;
    #     every realization gives the same BA stock for these censuses.
    k_postanchor_stk <- tree_census[
        treeID %in% uncertain_treeIDs & CensusID > ANCHOR_START_CENSUS,
        ..stock_cols
    ]
    k_stk <- rbindlist(list(
        fixed_stock[, ..stock_cols],
        k_var_stk[, ..stock_cols],
        k_postanchor_stk
    ), use.names = TRUE, fill = TRUE)
    k_stock_NtreeID[[k]] <- k_stk[, .(N_unique_treeid = uniqueN(treeID)), by = CensusID]
    all_stock_realizations[[k]] <- k_stk[,
        .(TotalBA_m2 = fsum(TotalBA_m2, na.rm = TRUE)),
        by = .(quadrat, CensusID)
    ][, realization := k]
    rm(k_stk, k_var, k_var_stk, k_postanchor_stk)
    gc()
}
all_quadrat_realizations <- rbindlist(all_quadrat_realizations)
all_stock_realizations <- rbindlist(all_stock_realizations)

all_quadrat_realizations[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
cat(
    "[BA] Realization loop complete:", K_realizations, "draws;",
    nrow(all_quadrat_realizations), "quadrat-flux rows;",
    nrow(all_stock_realizations), "stock rows\n"
)

# ── MC summary tables: empirical 95 % CI collapsed across realizations ────────
# all_quadrat_summary – quadrat × interval mean/sd/95 % CI for each flux component.
# all_stock_summary   – quadrat × census mean/sd/95 % CI for BA stock.
all_quadrat_summary <- summarise_flux(
    all_quadrat_realizations,
    c("quadrat", "CensusID_from", "CensusID_to")
)
all_stock_summary <- all_stock_realizations[,
    {
        q <- fquantile(TotalBA_m2, c(0.025, 0.975))
        .(
            TotalBA_mean = fmean(TotalBA_m2), TotalBA_sd = fsd(TotalBA_m2),
            TotalBA_lwr  = q[1L],             TotalBA_upr = q[2L]
        )
    },
    by = .(quadrat, CensusID)
]

# ---- Write outputs to disk --------------------------------------------------
out_dir <- "./BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_feather(map_tree_change, file.path(out_dir, "ba_map_change_treeID.feather"))
write_feather(map_quadrat_change, file.path(out_dir, "ba_map_change_quadrat.feather"))
write_feather(map_quadrat_stock, file.path(out_dir, "ba_map_stock_quadrat.feather"))
write_feather(all_quadrat_realizations, file.path(out_dir, "ba_mc_realizations_quadrat.feather"))
write_feather(all_stock_realizations, file.path(out_dir, "ba_mc_realizations_stock_quadrat.feather"))
write_feather(all_quadrat_summary, file.path(out_dir, "ba_mc_summary_quadrat.feather"))
write_feather(all_stock_summary, file.path(out_dir, "ba_mc_summary_stock_quadrat.feather"))
cat("[BA] Outputs written to", out_dir, "\n")

# ============================================================
# SECTION 5: Figures
# ============================================================
# Three figures compare MAP estimates against MC uncertainty using a consistent
# visual grammar throughout:
#   MAP  – bold solid line (+ bootstrap 95 % CI ribbon for stocks)
#   MC   – empirical 95 % CI ribbon + dashed center line
#            (mean or median, per mc_center)
#
# Figure 1 — Forest-level BA stock per hectare across all censuses.
#   MC layers are restricted to censuses <= ANCHOR_START_CENSUS because
#   post-anchor census stocks for uncertain trees are deterministic (no
#   identity uncertainty) and the spaghetti would collapse to a single line.
#
# Figure 2 — Forest-level BA flux components (Growth, Loss, Gain, Net ΔBA)
#   per hectare per year across census intervals, faceted by component.
#   MC and MAP layers are both restricted to pre-anchor intervals
#   (CensusID_from < ANCHOR_START_CENSUS). Post-anchor intervals have no path
#   uncertainty; including them would place MAP and MC on top of each other
#   trivially and inflate the apparent certainty of the pre-anchor estimates.
#
# Figure 3 — Individual-tree BA trajectories for five focal trees,
#   one tree per page in a single multi-page PDF. Each page overlays the
#   observed stem BA record (path 0, charcoal) against all posterior paths
#   (coloured lines), faceted by reconstruction path.
#
# KEY DESIGN RULES (applied uniformly):
# • MAP center = fmean() across quadrats. Using fmedian() would yield ~0 for
#   zero-inflated Loss/Gain distributions (most quadrats have no deaths or
#   recruits per interval), placing the MAP line far below the MC ribbon.
# • mc_center (mean or median) controls only the MC dashed summary line.
# • All MC layers respect ANCHOR_START_CENSUS filter.
# ============================================================

library(scales)

# ── Scaling & palettes ───────────────────────────────────────────────────────
scale_ha <- 10000 / (20 * 20) # 20 m × 20 m quadrats → per-hectare conversion

# Two-source palette (MAP vs MC)
pal <- c(MAP = "#1b7a56", MC = "#c4520a")

# Four flux-component palette
flux_pal <- c(
    Growth_BA     = "#1b7a56",
    Loss_BA       = "#c4520a",
    Gain_BA       = "#5b57a8",
    DeltaBA_total = "#c4186a"
)

# Individual-tree palette
COL_OBS <- "#3d3d3a" # charcoal → observed path (path 0)
COL_MOD <- "#2e9e75" # green    → imputed paths

# ── Shared theme ─────────────────────────────────────────────────────────────
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

# Shared guide: ribbon swatch + line inside the legend key
guide_two_source <- function() {
    guides(
        colour = guide_legend(override.aes = list(linewidth = 1.2, linetype = c("solid", "dashed"))),
        fill   = guide_legend(override.aes = list(alpha = 0.35))
    )
}

# ── center_fun: summarises the MC distribution across realizations ────────────
# Applied to the 250 per-realization forest-level means (one value per draw).
# NOTE: this is NOT used for the MAP — see the design rule in the section header.
center_fun <- function(x) {
    if (mc_center == "median") fmedian(x, na.rm = TRUE) else fmean(x, na.rm = TRUE)
}

# ── Figure 1: Forest-level BA stock per hectare ───────────────────────────────

# BUG FIX 1 & 2 — MAP stock bootstrap:
# Previously used center_fun() for both the center statistic and the bootstrap
# statistic. When mc_center == "median" this gave fmedian across raw per-quadrat
# BA values (one value per quadrat, ~1250 values in BCI). The MC center is
# center_fun() of 250 per-realization *means* — a distribution of means that
# clusters tightly near the overall forest mean by the CLT. The two "medians"
# measure completely different things and for skewed quadrat distributions they
# diverge badly, placing the MAP line far outside the MC ribbon.
# Fix: the MAP forest-level summary is always fmean() across quadrats, matching
# exactly what each MC realization computes (fmean across quadrats). The
# mc_center choice then only governs the MC posterior summary.
map_stock_boot <- map_quadrat_stock[
    , .(TotalBA_ha = TotalBA_m2 * scale_ha),
    by = .(quadrat, CensusID)
][,
    {
        v <- TotalBA_ha
        n <- length(v)
        B <- 1000L
        bm <- vapply(
            seq_len(B),
            function(i) fmean(v[sample.int(n, n, replace = TRUE)]), # FIX: was center_fun
            numeric(1L)
        )
        .(
            center = fmean(v), # FIX: was center_fun(v)
            lwr = fquantile(bm, 0.025),
            upr = fquantile(bm, 0.975)
        )
    },
    by = CensusID
]
map_stock_boot <- merge(map_stock_boot, dates[, .(CensusID, Year)], by = "CensusID")
map_stock_boot <- map_stock_boot[CensusID >= first_plot_census]

# MC: per-realization forest-level stock (mean per-ha across all quadrats).
# Restrict to pre-anchor censuses (<=ANCHOR_START_CENSUS): post-anchor stocks
# for uncertain trees are identical across paths (deterministic stemID),
# so the spaghetti would collapse to a single line there — not informative.
mc_stock_per_real <- all_stock_realizations[,
    .(TotalBA_ha = fmean(TotalBA_m2 * scale_ha)),
    by = .(CensusID, realization)
]
mc_stock_per_real <- merge(mc_stock_per_real, dates[, .(CensusID, Year)], by = "CensusID")
mc_stock_per_real <- mc_stock_per_real[CensusID >= first_plot_census & CensusID <= ANCHOR_START_CENSUS]

# MC: collapse K realization means to center + 95 % empirical CI.
mc_stock_ci <- mc_stock_per_real[,
    {
        q <- fquantile(TotalBA_ha, c(0.025, 0.975))
        .(center = center_fun(TotalBA_ha), lwr = q[1L], upr = q[2L])
    },
    by = CensusID
]
mc_stock_ci <- merge(mc_stock_ci, dates[, .(CensusID, Year)], by = "CensusID")
mc_stock_ci <- mc_stock_ci[CensusID >= first_plot_census & CensusID <= ANCHOR_START_CENSUS]

fig1 <- ggplot() +
    # # MC: spaghetti (one line per realization, pre-anchor censuses only)
    # geom_line(
    #     data = mc_stock_per_real,
    #     aes(Year, TotalBA_ha, group = realization),
    #     colour = pal["MC"], alpha = 0.1, linewidth = 0.25
    # ) +
    # MC: 95 % empirical CI ribbon (pre-anchor)
    geom_ribbon(
        data = mc_stock_ci,
        aes(Year, ymin = lwr, ymax = upr, fill = "MC"),
        alpha = 0.22
    ) +
    # MC: center line (dashed, pre-anchor)
    geom_line(
        data = mc_stock_ci,
        aes(Year, center, colour = "MC"),
        linewidth = 0.9, linetype = "dashed"
    ) +
    # MAP: bootstrap 95 % CI ribbon
    geom_ribbon(
        data = map_stock_boot,
        aes(Year, ymin = lwr, ymax = upr, fill = "MAP"),
        alpha = 0.22
    ) +
    # MAP: forest mean (solid, thicker — primary reference)
    geom_line(
        data = map_stock_boot,
        aes(Year, center, colour = "MAP"),
        linewidth = 1.4
    ) +
    scale_colour_manual(
        "Estimate",
        values = pal,
        labels = c(MAP = "MAP mean (bootstrap 95 % CI)", MC = "MC (empirical 95 % CI)")
    ) +
    scale_fill_manual(
        "Estimate",
        values = pal,
        labels = c(MAP = "MAP mean (bootstrap 95 % CI)", MC = "MC (empirical 95 % CI)")
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(labels = scales::label_comma()) +
    guide_two_source() +
    labs(
        title = "Forest-level BA stock per hectare",
        subtitle = sprintf(
            "MAP: solid mean + bootstrap 95 %% CI  ·  MC: dashed %s + empirical 95 %% CI",
            mc_center # mc_center governs only the MC dashed line, not the MAP
        ),
        x = "Year",
        y = expression("BA (m"^2 ~ "ha"^
            {
                -1
            } * ")")
    ) +
    theme_forest()

print(fig1)

# ── Figure 2: Forest-level annual BA flux components per hectare ────────────
# Flux components are annualised (divided by census interval length in years)
# and scaled to per-hectare units. All MC and MAP layers are restricted to
# pre-anchor intervals (CensusID_from < ANCHOR_START_CENSUS). Post-anchor
# intervals (C7→C8, C8→C9) are excluded: their MC paths are deterministic
# (zero identity uncertainty), so showing them would misleadingly suggest the
# MC ribbon collapses post-C7 due to higher certainty rather than absence of
# sampling. Forest-level values are always averaged with fmean() across
# quadrats (see MAP design rule in Section 5 header).

flux_labels <- c(
    Growth_BA     = "Growth",
    Loss_BA       = "Loss (mortality)",
    Gain_BA       = "Gain (recruitment + ingrowth)",
    DeltaBA_total = "Net \u0394BA"
)

census_pairs[, Interval_yr := as.numeric(difftime(Date_to, Date_from, units = "days")) / 365.25]

# MAP flux — restricted to pre-anchor intervals.
map_flux_center <- merge(
    map_quadrat_change, # [CensusID_from < ANCHOR_START_CENSUS],
    census_pairs[, .(CensusID_from, Interval_yr)],
    by = "CensusID_from",
    all.x = TRUE
)
map_flux_center[, (flux_cols) := lapply(.SD, function(x) x / Interval_yr), .SDcols = flux_cols]
map_flux_center <- map_flux_center[,
    lapply(.SD, fmean, na.rm = TRUE),
    .SDcols = flux_cols,
    by = CensusID_from
]
map_flux_center[, (flux_cols) := lapply(.SD, `*`, scale_ha), .SDcols = flux_cols]
map_flux_center[, DeltaBA_total := Growth_BA + Loss_BA + Gain_BA]
map_flux_long <- melt(
    map_flux_center,
    id.vars = "CensusID_from", variable.name = "Component", value.name = "value"
)
map_flux_long <- merge(map_flux_long, census_pairs[, .(CensusID_from, Year_mid)], by = "CensusID_from")
map_flux_long <- map_flux_long[CensusID_from >= first_plot_census]

# MC flux — all_quadrat_realizations already contains only pre-anchor intervals
# (post_decomp was restricted above). Annualise and scale.
flux_cols <- c("Growth_BA", "Loss_BA", "Gain_BA", "DeltaBA_total")
mc_flux_mean <- merge(
    all_quadrat_realizations[CensusID_from < ANCHOR_START_CENSUS],
    census_pairs[, .(CensusID_from, Interval_yr)],
    by = "CensusID_from",
    all.x = TRUE
)
mc_flux_mean[, (flux_cols) := lapply(.SD, function(x) x / Interval_yr), .SDcols = flux_cols]
mc_flux_mean <- mc_flux_mean[,
    lapply(.SD, fmean, na.rm = TRUE),
    .SDcols = flux_cols,
    by = .(CensusID_from, realization)
]
mc_flux_mean[, (flux_cols) := lapply(.SD, `*`, scale_ha), .SDcols = flux_cols]
mc_flux_long <- melt(
    mc_flux_mean,
    id.vars = c("CensusID_from", "realization"),
    variable.name = "Component", value.name = "value"
)
mc_flux_long <- merge(mc_flux_long, census_pairs[, .(CensusID_from, Year_mid)], by = "CensusID_from")
mc_flux_long <- mc_flux_long[CensusID_from >= first_plot_census]

# MC: collapse K realizations to center + 95 % empirical CI.
mc_flux_ci <- mc_flux_long[,
    {
        q <- fquantile(value, c(0.025, 0.975))
        .(center = center_fun(value), lwr = q[1L], upr = q[2L])
    },
    by = .(CensusID_from, Component)
]
mc_flux_ci <- merge(mc_flux_ci, census_pairs[, .(CensusID_from, Year_mid)], by = "CensusID_from")

fig2 <- ggplot() +
    # # MC: spaghetti (pre-anchor intervals only)
    # geom_line(
    #     data = mc_flux_long,
    #     aes(Year_mid, value, group = realization, colour = Component),
    #     alpha = 0.12, linewidth = 0.7
    # ) +
    # MC: 95 % ribbon
    geom_ribbon(
        data = mc_flux_ci,
        aes(Year_mid, ymin = lwr, ymax = upr, fill = Component),
        alpha = 0.7
    ) +
    # MC: center (dashed)
    geom_line(
        data = mc_flux_ci,
        aes(Year_mid, center, colour = Component),
        linewidth = 0.7, linetype = "dashed"
    ) +
    # MAP: bold solid (pre-anchor intervals; aligned with MC data range)
    geom_line(
        data = map_flux_long,
        aes(Year_mid, value),
        linewidth = 1
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "#bbbbaa", linewidth = 0.4) +
    facet_wrap(
        ~Component,
        scales = "free_y", ncol = 1,
        labeller = as_labeller(flux_labels)
    ) +
    scale_colour_manual(values = flux_pal, guide = "none") +
    scale_fill_manual(values = flux_pal, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
        title = "Forest-level annual BA fluxes per hectare",
        subtitle = sprintf(
            "Bold solid = MAP (mean)  \u00b7  dashed %s + ribbon = MC uncertainty",
            mc_center
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

# ── Figure 3: BA trajectories for five focal trees — one page per tree ─────────
# For each tree in focal_trees we overlay every posterior reconstruction path
# against the observed (MAP) stem record. Panels within each page are faceted
# by path index; path 0 is the observed record (charcoal); modelled paths are
# coloured green. All five pages are written to a single multi-page PDF via
# pdf() so the file can be scrolled in any PDF reader.
#
# The function build_tree_page() encapsulates the plot logic so the loop is
# kept clean. It returns NULL (with a message) if the treeID is not found in
# all_paths or stem_dt, so missing trees do not crash the loop.

build_tree_page <- function(tid, all_paths_dt, stem_dt_all, census_dates,
                            dates_dt, first_census) {
    # ── Subset posterior paths for this tree ──────────────────────────────────
    ap <- all_paths_dt[treeID == tid]
    if (nrow(ap) == 0L) {
        message(sprintf("[Fig3] treeID %d not found in all_paths — skipping.", tid))
        return(NULL)
    }
    ap[, group_id := NULL]
    ap[, StemPaths := NULL]

    # ── Add observed path (path_idx = 0) from stem_dt ────────────────────────
    obs <- stem_dt_all[treeID == tid]
    if (nrow(obs) == 0L) {
        message(sprintf("[Fig3] treeID %d not found in stem_dt — skipping.", tid))
        return(NULL)
    }
    obs[, path_idx := 0L]
    obs[, path_prob := 1.0]
    keep_cols <- intersect(names(ap), names(obs))
    obs <- obs[, ..keep_cols]

    pp <- rbindlist(list(ap[, ..keep_cols], obs), use.names = TRUE, fill = TRUE)
    pp[, path_type := fifelse(path_idx == 0L, "Observed", "Modelled")]
    pp <- merge(pp, dates_dt[, .(CensusID, Year)], by = "CensusID")
    pp <- pp[CensusID >= first_census]

    # ── Colour palette: charcoal for observed, gradient of greens for paths ───
    n_mod <- uniqueN(pp[path_idx != 0L, path_idx])
    mod_pal <- colorRampPalette(c("#a8d8c2", COL_MOD))(max(n_mod, 1L))
    path_ids <- sort(unique(pp$path_idx))
    path_col <- setNames(c(COL_OBS, mod_pal), c(0L, path_ids[path_ids != 0L]))

    # ── Facet labeller ────────────────────────────────────────────────────────
    path_labeller <- labeller(path_idx = function(x) {
        ifelse(x == "0", "Observed (path 0)", paste0("Path ", x))
    })

    ggplot(pp, aes(
        x      = Year,
        y      = BA,
        group  = interaction(stemID, path_idx),
        colour = factor(path_idx)
    )) +
        geom_line(
            data = pp[path_type == "Modelled"],
            linewidth = 0.55, alpha = 0.55
        ) +
        geom_line(
            data      = pp[path_type == "Observed"],
            linewidth = 1.1
        ) +
        geom_point(
            data = pp[path_type == "Observed"],
            size = 1.8, shape = 21, fill = "white", stroke = 0.8
        ) +
        scale_colour_manual(values = path_col, guide = "none") +
        scale_x_continuous(breaks = scales::pretty_breaks(4)) +
        scale_y_continuous(labels = scales::label_comma()) +
        facet_wrap(
            ~path_idx,
            scales    = "free_y",
            ncol      = 2,
            labeller  = path_labeller
        ) +
        labs(
            title = paste0("BA trajectories \u00b7 treeID ", tid),
            subtitle = paste0(
                "Each panel = one posterior reconstruction path  \u00b7  ",
                "Path 0 = observed record (charcoal)  \u00b7  ",
                n_mod, " modelled path(s)"
            ),
            x = "Year",
            y = expression("BA (m"^2 * ")"),
            caption = "Modelled paths in green; observed in charcoal"
        ) +
        theme_forest()
}

# ── Render: one page per focal tree in a single multi-page PDF ────────────────

# Five treeIDs are selected automatically for Figure 3 from those with
# 10 posterior paths, so the figure shows trees with identity ambiguity.
# Replace the selection logic below if you want a fixed manual set instead.
inc <- all_paths[, .(N_paths = uniqueN(path_idx)), by = treeID][N_paths == 10]$treeID

set.seed(42)
focal_trees <- sample(inc, 5L, replace = FALSE)

fig3_path <- file.path(out_dir, "fig3_BA_trajectories.pdf")
pdf(fig3_path, width = 9, height = 8)
for (tid in focal_trees) {
    pg <- build_tree_page(
        tid          = tid,
        all_paths_dt = all_paths,
        stem_dt_all  = stem_dt,
        census_dates = census_pairs,
        dates_dt     = dates,
        first_census = first_plot_census
    )
    if (!is.null(pg)) print(pg)
}
dev.off()
cat("[Fig3] Written:", fig3_path, "\n")

# ── Save Figures 1 and 2 ──────────────────────────────────────────────────────
ggsave(file.path(out_dir, "fig1_BA_stock.pdf"), fig1, width = 8, height = 4.5)
ggsave(file.path(out_dir, "fig2_BA_fluxes.pdf"), fig2, width = 8, height = 10)
cat("[Figs] fig1_stock.pdf and fig2_fluxes.pdf saved to", out_dir, "\n")
