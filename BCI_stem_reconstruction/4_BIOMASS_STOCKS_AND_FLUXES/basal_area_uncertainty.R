rm(list = ls())

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(collapse)
library(arrow)

# ============================================================
# SECTION 1: Load and prepare census data
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
# FIXME
rec <- rec[gx <= 200 & gy <= 200]

post_file <- "./BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds"

# ============================================================
# SECTION 2: Helper functions
# ============================================================
# ba_m2()         – converts stem DBH (mm) to basal area in m²
#                   using the standard circle formula.
#
# parse_recon()   – parses a single "StemPaths:StemID;..." string
#                   (one row from the posterior table) into a
#                   two-column data.table ready for merging with
#                   the census observations.
#
# decompose_ba()  – classifies stems as survivor / death / recruit
#                   by checking NA patterns in BA_from / BA_to,
#                   then sums Growth, Loss, Gain, and Net ΔBA.
#                   Shared by both the MAP loop (Section 3) and
#                   the posterior loop (Section 4) so the logic
#                   is defined exactly once.
#
# summarise_flux()– across-realization summary for the four flux
#                   columns: mean, sd, and empirical 95 % CI.
#                   Used to produce the MC summary tables written
#                   to disk; by_cols controls the grouping level
#                   (quadrat or plot).
# ============================================================

ba_m2 <- function(dbh_mm) pi / 4 * (dbh_mm / 1000)^2

parse_recon <- function(recon_str) {
    # recon_str: a single character string e.g. "3:101;3:205;7:88"
    pairs <- strsplit(recon_str, ";", fixed = TRUE)[[1L]]
    parts <- strsplit(pairs, ":", fixed = TRUE)
    data.table(
        StemPaths = as.integer(vapply(parts, `[`, character(1L), 1L)),
        StemID    = as.integer(vapply(parts, `[`, character(1L), 2L))
    )
}

# Classify stems and compute BA change components within each group.
# m must have columns BA_from and BA_to; by_cols defines the aggregation groups.
# Used identically for both MAP (map_change) and posterior (post_decomp) loops.
decompose_ba <- function(m, by_cols) {
    m[, status := fifelse(
        !is.na(BA_from) & !is.na(BA_to), "survivor",
        fifelse(!is.na(BA_from), "death", "recruit")
    )]
    m[, .(
        Growth_BA     = fsum((BA_to - BA_from) * (status == "survivor"), na.rm = TRUE),
        Loss_BA       = -fsum(BA_from * (status == "death"), na.rm = TRUE),
        Gain_BA       = fsum(BA_to * (status == "recruit"), na.rm = TRUE),
        DeltaBA_total = fsum(BA_to, na.rm = TRUE) - fsum(BA_from, na.rm = TRUE),
        NumSurvivors  = fsum(status == "survivor"),
        NumDeaths     = fsum(status == "death"),
        NumRecruits   = fsum(status == "recruit")
    ), by = by_cols]
}

# Summarise all four flux columns across realizations: mean, sd, empirical 95 % CI.
# Applied at both the tree and quadrat levels via the by_cols argument.
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

# ============================================================
# SECTION 3: MAP baseline BA decomposition
# ============================================================
# Uses the observed (maximum a posteriori) DBH measurements to
# compute deterministic BA stocks and fluxes. Outputs:
#   map_tree_change    – tree-level flux per census pair
#   map_quadrat_change – quadrat-level flux per census pair
#                        (sum of trees within each quadrat)
#   map_quadrat_stock  – quadrat-level BA stock per census
# These serve as (a) the fixed component in every MC realization
# for trees with a single reconstruction path, and (b) the bold
# reference line in the figures.
# ============================================================

tree_census <- rec[!is.na(dbh) & !is.na(stemID) & !is.na(treeID),
    .(TotalBA_m2 = fsum(ba_m2(dbh)), NumStems = .N),
    by = .(quadrat, treeID, CensusID)
]
setorder(tree_census, treeID, CensusID)
cat("[BA] treeID x census rows:", nrow(tree_census), "\n")

dates <- rec[!is.na(treeID), .(Date = median(ExactDate, na.rm = TRUE)), by = CensusID][order(CensusID)]
stopifnot("Need >= 2 censuses" = nrow(dates) >= 2L)
census_pairs <- data.table(
    CensusID_from = dates$CensusID[-nrow(dates)],
    CensusID_to   = dates$CensusID[-1L],
    Date_from     = dates$Date[-nrow(dates)],
    Date_to       = dates$Date[-1L]
)

stem_dt <- rec[
    !is.na(treeID) & !is.na(stemID),
    .(treeID, StemID = stemID, CensusID, BA = ba_m2(dbh), quadrat)
]

# MAP census-pair decomposition using shared decompose_ba().
map_change <- rbindlist(lapply(seq_len(nrow(census_pairs)), function(i) {
    cf <- census_pairs$CensusID_from[i]
    ct <- census_pairs$CensusID_to[i]
    sf <- stem_dt[CensusID == cf, .(quadrat, treeID, StemID, BA_from = BA)]
    st <- stem_dt[CensusID == ct, .(quadrat, treeID, StemID, BA_to = BA)]
    d <- decompose_ba(
        merge(sf, st, by = c("quadrat", "treeID", "StemID"), all = TRUE),
        by_cols = c("quadrat", "treeID")
    )
    d[, `:=`(
        CensusID_from = cf, CensusID_to = ct,
        Date_from = census_pairs$Date_from[i],
        Date_to = census_pairs$Date_to[i]
    )]
}))
map_change[, Interval_yr := as.numeric(difftime(Date_to, Date_from, units = "days")) / 365.25]
map_change[is.na(Interval_yr) | Interval_yr <= 0, Interval_yr := 5.0]
cat(
    "[BA] MAP decomposition:", nrow(map_change), "intervals across",
    uniqueN(map_change$treeID), "treeIDs\n"
)

flux_cols <- c("Growth_BA", "Loss_BA", "Gain_BA", "DeltaBA_total")

# Quadrat-level MAP aggregations (deterministic; no CI needed here).
map_tree_change <- copy(map_change)

# Fluxes: 8 census intervals
map_quadrat_change <- map_tree_change[,
    lapply(.SD, fsum, na.rm = TRUE),
    .SDcols = flux_cols,
    by = .(quadrat, CensusID_from, CensusID_to)
]

# Stocks: 9 censuses
map_quadrat_stock <- tree_census[,
    .(TotalBA_m2 = fsum(TotalBA_m2, na.rm = TRUE)),
    by = .(quadrat, CensusID)
]

# ============================================================
# SECTION 4: Posterior uncertainty propagation (MC realizations)
# ============================================================
# Reads the posterior reconstruction paths. Trees with a single
# path contribute their MAP fluxes/stocks to every realization
# (fixed component). Trees with multiple paths have one path
# drawn proportionally to path_prob in each realization (variable
# component). For K_realizations draws this produces:
#   all_quadrat_realizations – per-realization quadrat × interval
#                              flux table (all trees combined)
#   all_stock_realizations   – per-realization quadrat × census
#                              stock table
# These are then collapsed to empirical 95 % CI summary tables
# and written alongside the MAP outputs to disk.
# ============================================================

K_realizations <- 100L
set.seed(42L)
# mc_center <- "mean"    # choose "mean" or "median"
mc_center <- "median" # choose "mean" or "median"
mc_center <- match.arg(mc_center, c("mean", "median"))

post_full <- as.data.table(readRDS(post_file))
post_full[, group_id := treeID]

# Add path index if not present in the posterior file.
if (!"path_idx" %in% names(post_full)) {
    post_full[, path_idx := seq_len(.N), by = group_id]
}

# Guard against missing path_prob; normalise within each group.
post_full[, path_prob := fifelse(is.na(path_prob), 0, path_prob)]
post_full[, path_prob := {
    s <- fsum(path_prob)
    if (s > 0) path_prob / s else rep(1 / .N, .N)
}, by = group_id]
cat("[BA] Posterior paths:", nrow(post_full), "for", uniqueN(post_full$group_id), "groups\n")

# Parse reconstruction strings → long stem table.
# Each row in post_full has a "recon" string like "StemPaths:StemID;..."
# parse_recon() expands it to one row per stem, which is then merged
# with the observed DBH records to obtain BA values per path × census.
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

obs_lookup <- rec[
    !is.na(treeID) & !is.na(StemPaths),
    .(quadrat, treeID, StemPaths, CensusID, dbh)
]
all_paths <- merge(all_paths, obs_lookup, by = c("treeID", "StemPaths"), all.x = TRUE)
all_paths <- all_paths[!is.na(dbh)]
all_paths[, `:=`(BA = ba_m2(dbh), dbh = NULL)]
rm(obs_lookup)
gc()
setkey(all_paths, treeID)
cat("[BA] Parsed path observations:", nrow(all_paths), "\n")

# Posterior census-pair decomposition — reuses decompose_ba() identically to MAP.
# path_prob is dropped from sf/st to avoid a suffix collision in the full-outer merge.
post_decomp <- rbindlist(lapply(seq_len(nrow(census_pairs)), function(i) {
    cf <- census_pairs$CensusID_from[i]
    ct <- census_pairs$CensusID_to[i]
    sf <- all_paths[CensusID == cf, .(group_id, quadrat, treeID, path_idx, StemID, BA_from = BA)]
    st <- all_paths[CensusID == ct, .(group_id, quadrat, treeID, path_idx, StemID, BA_to = BA)]
    d <- decompose_ba(
        merge(sf, st, by = c("group_id", "quadrat", "treeID", "path_idx", "StemID"), all = TRUE),
        by_cols = c("group_id", "quadrat", "treeID", "path_idx")
    )
    d[, `:=`(CensusID_from = cf, CensusID_to = ct)]
}))
cat("[BA] Posterior decompositions:", nrow(post_decomp), "rows\n")

# Per-path tree-level BA stocks — used for individual tree trajectory plots.
tree_stock <- all_paths[,
    .(TotalBA_m2 = fsum(BA, na.rm = TRUE)),
    by = .(group_id, quadrat, treeID, path_idx, CensusID)
]

# Partition: uncertain = multi-path groups; fixed = single MAP path.
path_group_sizes <- post_full[, .N, by = group_id]
multi_path_groups <- path_group_sizes[N > 1L, group_id]
uncertain_treeIDs <- unique(post_full[group_id %in% multi_path_groups, treeID])
fixed_tree_change <- map_change[!treeID %in% uncertain_treeIDs]
fixed_stock <- tree_census[!treeID %in% uncertain_treeIDs]
cat(
    "[BA] Uncertain groups:", length(multi_path_groups),
    "; uncertain treeIDs:", length(uncertain_treeIDs), "\n"
)

# Column sets for consistent subset-and-bind inside the realization loop.
tree_fixed_cols <- c("quadrat", "treeID", "CensusID_from", "CensusID_to", flux_cols)
stock_cols <- c("quadrat", "treeID", "CensusID", "TotalBA_m2")

out_dir <- "./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/outputs"
realization_dir <- file.path(out_dir, "ba_mc_realizations_treeID")
if (!dir.exists(realization_dir)) dir.create(realization_dir, recursive = TRUE)

all_quadrat_realizations <- vector("list", K_realizations)
all_stock_realizations <- vector("list", K_realizations)

if (length(multi_path_groups) > 0L) {
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

        sel <- sampled_paths[.(k), on = "realization", .(group_id, path_idx)]
        k_var <- post_decomp[sel, on = .(group_id, path_idx), nomatch = 0L]
        k_var_stk <- tree_stock[sel, on = .(group_id, path_idx), nomatch = 0L]

        # Tree-level: fixed MAP trees + uncertain trees under this realization's paths.
        k_tree <- rbindlist(list(
            fixed_tree_change[, ..tree_fixed_cols],
            k_var[, ..tree_fixed_cols]
        ), use.names = TRUE)
        k_tree[, realization := k]

        write_feather(
            k_tree,
            file.path(realization_dir, sprintf("ba_mc_realization_treeID_%03d.feather", k))
        )

        # Quadrat-level: aggregate tree fluxes within each quadrat × interval.
        all_quadrat_realizations[[k]] <- k_tree[,
            lapply(.SD, fsum, na.rm = TRUE),
            .SDcols = flux_cols,
            by = .(quadrat, CensusID_from, CensusID_to)
        ][, realization := k]

        rm(k_tree)
        gc()

        # Stock-level: sum tree BAs within each quadrat × census.
        k_stk <- rbindlist(list(
            fixed_stock[, ..stock_cols],
            k_var_stk[, ..stock_cols]
        ), use.names = TRUE)
        all_stock_realizations[[k]] <- k_stk[,
            .(TotalBA_m2 = fsum(TotalBA_m2, na.rm = TRUE)),
            by = .(quadrat, CensusID)
        ][, realization := k]

        rm(k_stk, k_var, k_var_stk)
        gc()
    }

    all_quadrat_realizations <- rbindlist(all_quadrat_realizations)
    all_stock_realizations <- rbindlist(all_stock_realizations)
    gc()
    cat(
        "[BA] Tree realization files written:", K_realizations, ";",
        nrow(all_quadrat_realizations), "quadrat rows;",
        nrow(all_stock_realizations), "stock rows\n"
    )
} else {
    cat("[BA] No multi-path groups; realizations are identical to MAP.\n")
    all_quadrat_realizations <- copy(map_quadrat_change)[, realization := 1L]
    all_stock_realizations <- copy(map_quadrat_stock)[, realization := 1L]
    sampled_paths <- data.table(
        group_id = unique(tree_stock$group_id), path_idx = 1L, realization = 1L
    )
    setkey(sampled_paths, group_id, realization)
}

# ---- MC summary tables: empirical 95 % CI across realizations ---------------
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
out_dir <- "./BCI_stem_reconstruction/4_BIOMASS_STOCKS_AND_FLUXES/outputs"
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
# Four figures compare MAP estimates against MC uncertainty.
# A consistent visual grammar is applied throughout:
#   • MAP  – bold solid line (+ bootstrap 95 % CI ribbon for stocks)
#   • MC   – thin spaghetti lines + empirical 95 % CI ribbon
#            + dashed center line (mean or median, per mc_center)
#
# Figure 1: forest-level BA stock per hectare across all censuses
# Figure 2: forest-level BA flux components (Growth, Loss, Gain,
#            Net ΔBA) per hectare across census intervals, faceted
#            by component
# Figure 3a: per-path BA stock trajectories for one focal tree,
#             faceted by reconstruction path
# Figure 3b: annualised BA growth rate per census interval for the
#             same focal tree, modelled vs observed
#
# KEY DESIGN RULE (fixes Bugs 1 & 2): the MAP "center" is always
# fmean() across quadrats regardless of mc_center. mc_center
# controls only how the MC posterior distribution is summarised
# (mean vs median across the 250 realizations). Using fmedian()
# across raw per-quadrat values for the MAP is wrong because for
# zero-inflated fluxes (Loss_BA, Gain_BA) the median quadrat
# value is ~0 while the MC ribbon sits at the non-zero forest
# mean — placing the MAP line entirely outside the MC ribbon.
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

# MC: per-realization forest-level stock (mean per-ha across all quadrats).
mc_stock_per_real <- all_stock_realizations[,
    .(TotalBA_ha = fmean(TotalBA_m2 * scale_ha)),
    by = .(CensusID, realization)
]
# MC: collapse 250 realization means to center + 95 % empirical CI.
mc_stock_ci <- mc_stock_per_real[,
    {
        q <- fquantile(TotalBA_ha, c(0.025, 0.975))
        .(center = center_fun(TotalBA_ha), lwr = q[1L], upr = q[2L])
    },
    by = CensusID
]

fig1 <- ggplot() +
    # MC: spaghetti (one line per realization)
    geom_line(
        data = mc_stock_per_real[CensusID <= 7L],
        aes(CensusID, TotalBA_ha, group = realization),
        colour = pal["MC"], alpha = 0.12, linewidth = 0.4
    ) +
    # MC: 95 % empirical CI ribbon
    geom_ribbon(
        data = mc_stock_ci[CensusID <= 7L],
        aes(CensusID, ymin = lwr, ymax = upr, fill = "MC"),
        alpha = 0.22
    ) +
    # MC: center line (dashed)
    geom_line(
        data = mc_stock_ci[CensusID <= 7L],
        aes(CensusID, center, colour = "MC"),
        linewidth = 0.9, linetype = "dashed"
    ) +
    # MAP: bootstrap 95 % CI ribbon
    geom_ribbon(
        data = map_stock_boot,
        aes(CensusID, ymin = lwr, ymax = upr, fill = "MAP"),
        alpha = 0.22
    ) +
    # MAP: forest mean (solid, thicker — primary reference)
    geom_line(
        data = map_stock_boot,
        aes(CensusID, center, colour = "MAP"),
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
            "MAP: solid mean + bootstrap 95 %% CI  ·  MC: spaghetti + dashed %s + empirical 95 %% CI",
            mc_center # mc_center governs only the MC dashed line, not the MAP
        ),
        x = "Census",
        y = expression("BA (m"^2 ~ "ha"^
            {
                -1
            } * ")")
    ) +
    theme_forest()

print(fig1)

# ── Figure 2: Forest-level BA flux components per hectare ─────────────────────

flux_labels <- c(
    Growth_BA     = "Growth",
    Loss_BA       = "Loss (mortality)",
    Gain_BA       = "Gain (recruitment)",
    DeltaBA_total = "Net \u0394BA"
)

# BUG FIX 1 — MAP flux center:
# Previously used lapply(.SD, center_fun) which, when mc_center == "median",
# applied fmedian across all ~1250 per-quadrat flux values. For Loss_BA and
# Gain_BA most quadrats have zero deaths/recruits per interval, so the median
# across quadrats is ~0. The MC ribbon sits at the non-zero forest-level mean
# (because each MC realization first averages across quadrats, eliminating the
# zeros, then center_fun is applied across 250 near-identical means). The result
# is the MAP line stuck at 0 while the MC ribbon floats above it.
# Fix: always fmean() for the MAP forest-level flux, matching each MC realization.
map_flux_center <- map_quadrat_change[,
    lapply(.SD, fmean, na.rm = TRUE), # FIX: was center_fun (fmedian when mc_center=="median")
    .SDcols = flux_cols, by = CensusID_from
]
map_flux_center[, (flux_cols) := lapply(.SD, `*`, scale_ha), .SDcols = flux_cols]
map_flux_long <- melt(
    map_flux_center,
    id.vars = "CensusID_from", variable.name = "Component", value.name = "value"
)

# MC: per-realization forest-level fluxes (mean across quadrats, then scale).
mc_flux_mean <- all_quadrat_realizations[,
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

# MC: collapse to center + 95 % empirical CI across realizations.
mc_flux_ci <- mc_flux_long[,
    {
        q <- fquantile(value, c(0.025, 0.975))
        .(center = center_fun(value), lwr = q[1L], upr = q[2L])
    },
    by = .(CensusID_from, Component)
]

fig2 <- ggplot() +
    # MC: spaghetti
    geom_line(
        data = mc_flux_long[CensusID_from <= 7L],
        aes(CensusID_from, value, group = realization, colour = Component),
        alpha = 0.12, linewidth = 0.4
    ) +
    # MC: 95 % ribbon
    geom_ribbon(
        data = mc_flux_ci[CensusID_from <= 7L],
        aes(CensusID_from, ymin = lwr, ymax = upr, fill = Component),
        alpha = 0.22
    ) +
    # MC: center (dashed)
    geom_line(
        data = mc_flux_ci[CensusID_from <= 7L],
        aes(CensusID_from, center, colour = Component),
        linewidth = 0.8, linetype = "dashed"
    ) +
    # MAP: bold solid — BUG FIX 4: filter to same census range as MC layers
    # Previously map_flux_long had no filter, extending the MAP line one census
    # beyond the MC ribbon and making the last interval incomparable.
    geom_line(
        data = map_flux_long[CensusID_from <= 7L], # FIX: was map_flux_long (unfiltered)
        aes(CensusID_from, value, colour = Component),
        linewidth = 1.4
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
        title = "Forest-level BA fluxes per hectare",
        subtitle = sprintf(
            "Bold solid = MAP (mean)  \u00b7  Spaghetti + dashed %s + ribbon = MC uncertainty",
            mc_center
        ),
        x = "Census (from)",
        y = expression("BA flux (m"^2 ~ "ha"^
            {
                -1
            } * ")")
    ) +
    theme_forest()

print(fig2)

# ── Figure 3: Individual tree inspection ──────────────────────────────────────
# Shows one focal tree's BA trajectories across all posterior paths (3a) and
# its annualised growth rate per census interval, comparing the observed MAP
# record against the distribution across modelled paths (3b).

plot_trees <- 231607L

all_paths_plot <- all_paths[treeID %in% plot_trees]
all_paths_plot[, StemID := paste0(treeID, "_", StemID)]
all_paths_plot[, group_id := NULL]
all_paths_plot[, StemPaths := NULL]

stem_dt_plot <- stem_dt[treeID %in% all_paths_plot$treeID]
stem_dt_plot[, path_idx := 0L]
stem_dt_plot[, path_prob := 1L]

sort_names <- colnames(all_paths_plot)
stem_dt_plot <- stem_dt_plot[, ..sort_names]

plot_paths <- rbind(all_paths_plot, stem_dt_plot)
plot_paths[, path_type := fifelse(path_idx == 0L, "Observed", "Modelled")]

# Compute fluxes for the focal tree across all paths (path 0 = observed).
map_change_ind <- rbindlist(lapply(seq_len(nrow(census_pairs)), function(i) {
    cf <- census_pairs$CensusID_from[i]
    ct <- census_pairs$CensusID_to[i]
    sf <- plot_paths[CensusID == cf, .(quadrat, treeID, StemID, path_idx, BA_from = BA)]
    st <- plot_paths[CensusID == ct, .(quadrat, treeID, StemID, path_idx, BA_to = BA)]
    d <- decompose_ba(
        merge(sf, st, by = c("quadrat", "treeID", "StemID", "path_idx"), all = TRUE),
        by_cols = c("quadrat", "treeID", "StemID", "path_idx")
    )
    d[, `:=`(
        CensusID_from = cf, CensusID_to = ct,
        Date_from = census_pairs$Date_from[i],
        Date_to = census_pairs$Date_to[i]
    )]
}))
map_change_ind[, Interval_yr := as.numeric(difftime(Date_to, Date_from, units = "days")) / 365.25]
map_change_ind[is.na(Interval_yr) | Interval_yr <= 0, Interval_yr := 5.0]

# ── Figure 3a: BA stock trajectories faceted by path ─────────────────────────
n_mod <- uniqueN(plot_paths[path_idx != 0L, path_idx])
mod_pal <- colorRampPalette(c("#a8d8c2", COL_MOD))(max(n_mod, 1L))
path_ids <- sort(unique(plot_paths$path_idx))
path_cols <- setNames(
    c(COL_OBS, mod_pal),
    c(0L, path_ids[path_ids != 0L])
)

fig3a <- ggplot(
    plot_paths,
    aes(
        x      = CensusID,
        y      = BA,
        group  = interaction(StemID, path_idx),
        colour = factor(path_idx)
    )
) +
    geom_line(
        data = plot_paths[path_type == "Modelled"],
        linewidth = 0.55, alpha = 0.55
    ) +
    geom_line(
        data      = plot_paths[path_type == "Observed"],
        linewidth = 1.1
    ) +
    geom_point(
        data = plot_paths[path_type == "Observed"],
        size = 1.8, shape = 21, fill = "white", stroke = 0.8
    ) +
    scale_colour_manual(values = path_cols, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    scale_y_continuous(labels = scales::label_comma()) +
    facet_wrap(
        ~path_idx,
        scales = "free_y", ncol = 2,
        labeller = labeller(path_idx = function(x) {
            ifelse(x == "0", "Observed (path 0)", paste0("Path ", x))
        })
    ) +
    labs(
        title    = paste0("BA trajectories \u00b7 tree ", plot_trees),
        subtitle = "Each panel = one imputed path  \u00b7  Path 0 = observed record (charcoal)",
        x        = "Census",
        y        = expression("BA (m"^2 * ")"), # FIX 5: was cm²; ba_m2() returns m²
        caption  = "Modelled paths in green; observed in charcoal"
    ) +
    theme_forest()

print(fig3a)

# ── Figure 3b: Annualised BA growth rate by census interval ───────────────────
obs_data <- map_change_ind[path_idx == 0L]
mod_data <- map_change_ind[path_idx != 0L]

# BUG FIX 5 & 6 — annualise and correct units:
# Growth_BA was never divided by Interval_yr despite the subtitle/caption
# claiming "annual rate". Fixed by adding GrowthRate_m2yr = Growth_BA /
# Interval_yr. Also corrected the y-axis label from cm² to m² (ba_m2()
# returns m²; cm² appeared throughout Figs 3a and 3b).
obs_data[, GrowthRate_m2yr := Growth_BA / Interval_yr]
mod_data[, GrowthRate_m2yr := Growth_BA / Interval_yr]

mod_ribbon <- mod_data[, .(
    med = median(GrowthRate_m2yr, na.rm = TRUE),
    lo  = quantile(GrowthRate_m2yr, 0.10, na.rm = TRUE),
    hi  = quantile(GrowthRate_m2yr, 0.90, na.rm = TRUE)
), by = .(StemID, CensusID_to = CensusID_from + 1L)]

fig3b <- ggplot(
    obs_data,
    aes(x = CensusID_from + 1L, y = GrowthRate_m2yr, group = StemID)
) +
    # MC: 80 % interval ribbon
    geom_ribbon(
        data = mod_ribbon,
        aes(x = CensusID_to, ymin = lo, ymax = hi, group = StemID),
        fill = COL_MOD, alpha = 0.18,
        inherit.aes = FALSE
    ) +
    # MC: median (dashed)
    geom_line(
        data = mod_ribbon,
        aes(x = CensusID_to, y = med, group = StemID),
        colour = COL_MOD, linewidth = 0.9, linetype = "dashed",
        inherit.aes = FALSE
    ) +
    # Observed: solid + points
    geom_line(colour = COL_OBS, linewidth = 1.2) +
    geom_point(size = 2.0, shape = 21, fill = "white", colour = COL_OBS, stroke = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "#bbbbaa", linewidth = 0.4) +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
        title = paste0("Annual BA growth \u00b7 tree ", plot_trees),
        subtitle = "Charcoal solid = observed (MAP)  \u00b7  Green dashed + ribbon = modelled paths (median \u00b1 80 % interval)",
        x = "Census",
        y = expression("Annual BA growth (m"^2 ~ "yr"^
            {
                -1
            } * ")"), # FIX 5: was cm²
        caption = "Interval length standardised to annual rate"
    ) +
    theme_forest()

print(fig3b)

# ── Save all four figures ─────────────────────────────────────────────────────
ggsave(file.path(out_dir, "fig1_stock.pdf"), fig1, width = 8, height = 4.5)
ggsave(file.path(out_dir, "fig2_fluxes.pdf"), fig2, width = 8, height = 10)
ggsave(file.path(out_dir, "fig3a_trajectories.pdf"), fig3a, width = 9, height = 8)
ggsave(file.path(out_dir, "fig3b_growth.pdf"), fig3b, width = 8, height = 4.5)

# beepr::beep(2)
