rm(list = ls())

# ============================================================
# KEY PLOT INFORMATION: BCI 50-ha Plot — FULLY CORRECTED
# ============================================================
# BCI permanent plot: 1000 m (E–W) × 500 m (N–S) = 50 ha
# Base quadrats:      20 × 20 m (0.04 ha each); 1,250 total
#
# CORRECTIONS over original script:
#   (1) Per-hectare extrapolation via complete super-quadrat
#       grids of MULTIPLE sizes (20, 25, 50, 100, 125, 250 m),
#       each tiling the full 50-ha plot without gaps.
#       Bootstrap operates on per-super-quadrat values, not on
#       the raw total ÷ plot area — correctly propagating
#       quadrat-level spatial variation.
#   (2) Temporal per-ha trends (all censuses) use the same
#       bootstrapped quadrat approach instead of crude division.
#   (3) Complete grid zero-filling: empty super-quadrats are
#       included as zeros so the estimator is unbiased.
#   (4) Compositional change: linear trend per species across
#       all censuses + first-vs-last deltas for treeID counts
#       and basal area (winners / losers).
#   (5) New: demographic rates, species turnover (Jaccard),
#       diversity indices over time, DBH size-class shift,
#       spatial heatmaps (100 m sectors).
#   (6) Nine publication-ready plots saved to ./figures/.
# ============================================================

library(data.table)
library(boot)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(scales)

# ============================================================
# SECTION 1: CONSTANTS
# ============================================================

PLOT_X_M     <- 1000L                        # plot width  (m), E–W
PLOT_Y_M     <- 500L                         # plot height (m), N–S
PLOT_AREA_HA <- PLOT_X_M * PLOT_Y_M / 1e4   # = 50 ha
BASE_Q_M     <- 20L                          # base quadrat side (m)
N_CENSUSES   <- 9L
BOOT_R       <- 4999L
set.seed(42)

# Quadrat sizes (m, square) that tile 1000 × 500 exactly with NO gaps.
# Requirement: qs must divide both PLOT_X_M and PLOT_Y_M.
# GCD(1000, 500) = 500; valid sizes >= 20 m: 20, 25, 50, 100, 125, 250.
QUAD_SIZES_M <- c(20, 25, 50, 100, 125, 250)

# Approximate BCI census mid-years (1 through 9)
CENSUS_YEARS <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)

# Convenience: extract genus + species epithet from full Latin name
short_latin <- function(x) gsub("^(\\S+\\s+\\S+).*", "\\1", trimws(x))

# ============================================================
# SECTION 2: DATA LOADING
# ============================================================

message("Loading census data ...")
bci_nums    <- as.character(seq_len(N_CENSUSES))
census_list <- lapply(bci_nums, function(n) {
    fp <- sprintf("./data_paper_documentation/RTABLES/bci.stem%s.Rdata", n)
    if (!file.exists(fp)) stop("Missing file: ", fp)
    e <- new.env(parent = emptyenv())
    load(fp, envir = e)
    get(paste0("bci.stem", n), envir = e)
})
names(census_list) <- paste0("bci.stem", bci_nums)
rec <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rec <- rec[!is.na(quadrat)]
rm(census_list, bci_nums)
invisible(gc())

# ---- Species table --------------------------------------------------
load("./data_paper_documentation/RTABLES/bci.spptable.rdata")
bci.spptable <- as.data.table(bci.spptable)

# Build full Latin name (include infraspecific epithets when present)
if (all(c("InfraspecificRank", "InfraspecificEpithet") %in% names(bci.spptable))) {
    bci.spptable[, Latin := fifelse(
        !is.na(InfraspecificRank) &
            nchar(trimws(as.character(InfraspecificRank)))  > 0 &
            !is.na(InfraspecificEpithet) &
            nchar(trimws(as.character(InfraspecificEpithet))) > 0,
        paste(Latin, InfraspecificRank, InfraspecificEpithet),
        Latin
    )]
}
lf_cols  <- intersect(c("Lifeform_RFoster", "Lifeform_RPerez_SAguilar"), names(bci.spptable))
keep_spp <- c("sp", "Latin", lf_cols)
bci.spptable <- bci.spptable[, keep_spp, with = FALSE]
rec <- merge(rec, bci.spptable, by = "sp", all.x = TRUE)

# ---- Basal area (m²) per stem ---------------------------------------
rec[, ba_m2 := pi / 4 * (dbh / 1000)^2]

# ---- Ensure gx / gy are present ------------------------------------
# BCI stem tables carry gx/gy; derive from quadrat label as fallback.
if (!all(c("gx", "gy") %in% names(rec))) {
    message("  gx/gy not found — deriving centres from quadrat label (format 'XXYY')")
    rec[, gx := (as.integer(substr(quadrat, 1, 2)) + 0.5) * BASE_Q_M]
    rec[, gy := (as.integer(substr(quadrat, 3, 4)) + 0.5) * BASE_Q_M]
}

# ---- Census metadata lookup -----------------------------------------
census_meta <- data.table(
    censusID    = paste0("bci.stem", seq_len(N_CENSUSES)),
    census_num  = seq_len(N_CENSUSES),
    census_year = CENSUS_YEARS
)

# ============================================================
# SECTION 3: TREE-LEVEL AGGREGATION (alive stems only)
# ============================================================
# • Keep only Rstatus == "A" (alive resolved status).
# • Aggregate stems → tree (treeID) per census; sum BA, mean gx/gy.
# • gx/gy: for multi-stemmed trees take the mean position of alive
#   stems — sufficient for quadrat assignment.

rec_alive <- rec[Rstatus == "A" & !is.na(gx) & !is.na(gy)]

by_cols <- c("censusID", "CensusID", "treeID", "sp", "quadrat", "Latin", lf_cols)

rec_tree <- rec_alive[, .(
    ba_m2   = sum(ba_m2,  na.rm = TRUE),
    n_stems = .N,
    gx      = mean(gx,    na.rm = TRUE),
    gy      = mean(gy,    na.rm = TRUE)
), by = by_cols]

# Clip to strict plot interior (edge-case guard)
rec_tree[gx < 0,           gx := 0]
rec_tree[gx >= PLOT_X_M,   gx := PLOT_X_M - 1e-3]
rec_tree[gy < 0,           gy := 0]
rec_tree[gy >= PLOT_Y_M,   gy := PLOT_Y_M - 1e-3]

rec_tree <- merge(rec_tree, census_meta, by = "censusID", all.x = TRUE)
setorder(rec_tree, census_num, treeID)

# Validate census_num merge
if (any(is.na(rec_tree$census_num))) {
    warning("Some treeIDs have missing census_num; dropping affected rows.")
    rec_tree <- rec_tree[!is.na(census_num)]
}

all_cids   <- census_meta$censusID
FIRST_ID   <- all_cids[1]
LAST_ID    <- all_cids[N_CENSUSES]

last_trees <- rec_tree[censusID == LAST_ID]
last_stems <- rec_alive[censusID == LAST_ID]

total_trees <- nrow(last_trees)
total_stems <- nrow(last_stems)
total_ba    <- sum(last_trees$ba_m2, na.rm = TRUE)
n_species   <- uniqueN(last_trees$sp)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# ---- Bootstrap mean + 95% percentile CI ----
boot_mean_ci <- function(x, R = BOOT_R) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    if (n < 4) {
        se <- if (n > 1) sd(x) / sqrt(n) else 0
        return(c(mean = m, ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se))
    }
    b  <- boot::boot(x, function(d, i) mean(d[i]), R = R)
    ci <- tryCatch(
        boot::boot.ci(b, type = "perc")$percent[4:5],
        error = function(e) {
            se <- sd(x) / sqrt(n)
            c(m - 1.96 * se, m + 1.96 * se)
        }
    )
    c(mean = m, ci_lo = ci[1], ci_hi = ci[2])
}

# ---- Build complete super-quadrat grid and aggregate ----
# For quadrat size qs_m × qs_m:
#   • Trees assigned via floor(gx/qs_m), floor(gy/qs_m).
#   • COMPLETE grid created (n_x × n_y cells); empty cells = 0.
#   • Counts and BA scaled to per-ha.
# qs_m must divide both PLOT_X_M and PLOT_Y_M exactly.
make_superquad <- function(tree_dt, qs_m) {
    dt    <- copy(tree_dt)
    n_x   <- as.integer(PLOT_X_M / qs_m)
    n_y   <- as.integer(PLOT_Y_M / qs_m)
    qs_ha <- (qs_m * qs_m) / 1e4

    dt[, sxi  := pmin(as.integer(floor(gx / qs_m)), n_x - 1L)]
    dt[, syi  := pmin(as.integer(floor(gy / qs_m)), n_y - 1L)]
    dt[, sqid := sxi * n_y + syi]

    sq_obs <- dt[, .(
        n_trees = .N,
        ba_m2   = sum(ba_m2, na.rm = TRUE),
        n_sp    = uniqueN(sp)
    ), by = sqid]

    # Complete grid (zero-fill empties so the bootstrap is unbiased)
    sq_full <- merge(
        data.table(sqid = 0L:(n_x * n_y - 1L)),
        sq_obs, by = "sqid", all.x = TRUE
    )
    sq_full[is.na(n_trees), n_trees := 0L]
    sq_full[is.na(ba_m2),   ba_m2   := 0.0]
    sq_full[is.na(n_sp),    n_sp    := 0L]
    sq_full[, trees_ha := n_trees / qs_ha]
    sq_full[, ba_ha    := ba_m2   / qs_ha]
    sq_full
}

# ============================================================
# SECTION 4: OVERALL PLOT SUMMARY (last census)
# ============================================================

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  BCI 50-ha PLOT — Census %d (alive trees, DBH ≥ 10 mm)\n", N_CENSUSES))
cat(strrep("=", 64), "\n")
cat(sprintf("  Plot dimensions:                  %d × %d m = %g ha\n",
            PLOT_X_M, PLOT_Y_M, PLOT_AREA_HA))
cat(sprintf("  Total trees (treeIDs, alive):     %d\n",      total_trees))
cat(sprintf("  Total stems (stemIDs, alive):     %d\n",      total_stems))
cat(sprintf("  Species richness:                 %d\n",      n_species))
cat(sprintf("  Total basal area (m²):            %.2f\n",    total_ba))
cat(sprintf("  Crude trees ha⁻¹ (total/area):   %.1f\n",    total_trees / PLOT_AREA_HA))
cat(sprintf("  Crude BA ha⁻¹ (m²/ha):           %.2f\n",    total_ba    / PLOT_AREA_HA))

# ============================================================
# SECTION 5: SPECIES-LEVEL SUMMARY (last census)
# ============================================================

spp_sum <- last_trees[, .(
    n_trees     = .N,
    total_ba_m2 = sum(ba_m2, na.rm = TRUE),
    mean_ba_m2  = mean(ba_m2, na.rm = TRUE)
), by = .(sp, Latin)]
spp_sum[, rel_abund_pct := round(100 * n_trees     / total_trees, 3)]
spp_sum[, rel_ba_pct    := round(100 * total_ba_m2 / total_ba,    3)]
spp_sum[, IV            := rel_abund_pct + rel_ba_pct]   # Importance Value (max 200)

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY ABUNDANCE (Census 9)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[order(-n_trees),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct)], 10))

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY BASAL AREA (m²)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[order(-total_ba_m2),
    .(Latin, total_ba_m2, rel_ba_pct, n_trees, rel_abund_pct)], 10))

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY IMPORTANCE VALUE (abundance + BA, max = 200)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)], 10))

# ============================================================
# SECTION 6: ABUNDANCE / BA THRESHOLDS
# ============================================================

avg_n  <- mean(spp_sum$n_trees)
avg_ba <- mean(spp_sum$total_ba_m2)
q75_n  <- quantile(spp_sum$n_trees, 0.75)
q75_ba <- quantile(spp_sum$total_ba_m2, 0.75)

cat("\n", strrep("-", 64), "\n")
cat("  ABUNDANCE / BA THRESHOLDS\n")
cat(strrep("-", 64), "\n")
cat(sprintf("  Mean trees/species:       %6.1f  | 75th pct: %6.0f\n", avg_n, q75_n))
cat(sprintf("  Mean BA/species (m²):   %8.4f  | 75th pct: %8.4f\n", avg_ba, q75_ba))
cat("\n  >> Species above MEAN in both abundance & BA:\n")
print(spp_sum[n_trees > avg_n & total_ba_m2 > avg_ba][order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)])
cat("\n  >> Species above 75th pct in both abundance & BA:\n")
print(spp_sum[n_trees > q75_n & total_ba_m2 > q75_ba][order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)])

# ============================================================
# SECTION 7: LIFEFORM BREAKDOWN
# ============================================================

for (lfc in lf_cols) {
    lf <- last_trees[!is.na(get(lfc)), .(
        n_trees   = .N,
        n_species = uniqueN(sp),
        ba_m2     = sum(ba_m2, na.rm = TRUE)
    ), by = c(lfc)][order(-n_trees)]
    lf[, rel_n_pct  := round(100 * n_trees / total_trees, 1)]
    lf[, rel_ba_pct := round(100 * ba_m2   / total_ba,    1)]
    cat("\n", strrep("-", 64), "\n")
    cat(sprintf("  LIFEFORM BREAKDOWN — %s\n", lfc))
    cat(strrep("-", 64), "\n")
    print(lf)

    # Top 5 species per group (by abundance)
    cat(sprintf("\n  >> Top 5 species per group (%s):\n", lfc))
    for (g in sort(unique(na.omit(last_trees[[lfc]])))) {
        s <- last_trees[get(lfc) == g, .N, by = .(sp, Latin)][order(-N)]
        cat(sprintf("  [%s]\n", g))
        print(head(s[, .(Latin, N)], 5))
    }
}

# ============================================================
# SECTION 8: PER-HECTARE EXTRAPOLATION — MULTIPLE QUADRAT SIZES
# ============================================================
# KEY CORRECTION: For each quadrat size qs (last census):
#   1. Assign trees to qs×qs super-quadrats using (gx, gy).
#   2. Build the COMPLETE grid (n_x × n_y cells); fill empties = 0.
#   3. Convert each super-quadrat total to trees/ha and BA/ha.
#   4. Bootstrap the MEAN across super-quadrats → point estimate
#      + 95% CI that properly reflects spatial heterogeneity.
# This is correct because the unit of replication is the super-
# quadrat, not the whole plot, so quadrat-level variance is
# propagated into the CI rather than ignored.

message("Computing per-ha bootstrap by quadrat size ...")

perha_rows <- lapply(QUAD_SIZES_M, function(qs) {
    qs_ha <- (qs * qs) / 1e4
    message(sprintf("  %d × %d m (%.4f ha) ...", qs, qs, qs_ha))
    sq <- make_superquad(last_trees, qs)
    rt <- boot_mean_ci(sq$trees_ha)
    rb <- boot_mean_ci(sq$ba_ha)
    rs <- boot_mean_ci(sq$n_sp)
    data.table(
        quad_m      = qs,       quad_ha     = qs_ha,   n_quads     = nrow(sq),
        trees_ha    = rt["mean"], trees_ci_lo = rt["ci_lo"], trees_ci_hi = rt["ci_hi"],
        ba_ha       = rb["mean"], ba_ci_lo    = rb["ci_lo"], ba_ci_hi    = rb["ci_hi"],
        sp_per_quad = rs["mean"], sp_ci_lo    = rs["ci_lo"], sp_ci_hi    = rs["ci_hi"]
    )
})
perha_dt <- rbindlist(perha_rows)

cat("\n", strrep("=", 72), "\n")
cat(sprintf("  PER-HECTARE ESTIMATES BY QUADRAT SIZE  (boot 95%% CI, R = %d)\n", BOOT_R))
cat(strrep("=", 72), "\n")
print(perha_dt[, .(
    quad_m, quad_ha, n_quads,
    trees_ha    = round(trees_ha,    1),
    trees_ci_lo = round(trees_ci_lo, 1),
    trees_ci_hi = round(trees_ci_hi, 1),
    ba_ha       = round(ba_ha,       2),
    ba_ci_lo    = round(ba_ci_lo,    2),
    ba_ci_hi    = round(ba_ci_hi,    2)
)])

# ============================================================
# SECTION 9: DIVERSITY INDICES (last census)
# ============================================================

pi_vals    <- spp_sum$n_trees / total_trees
H_shannon  <- -sum(pi_vals * log(pi_vals))
D_conc     <- sum(pi_vals^2)
J_evenness <- H_shannon / log(n_species)

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  DIVERSITY INDICES — Census %d\n", N_CENSUSES))
cat(strrep("=", 64), "\n")
cat(sprintf("  Species richness (S):              %d\n",     n_species))
cat(sprintf("  Shannon-Wiener (H'):               %.4f\n",   H_shannon))
cat(sprintf("  Max H' [= ln(S)]:                  %.4f\n",   log(n_species)))
cat(sprintf("  Pielou's evenness (J = H'/ln S):   %.4f\n",   J_evenness))
cat(sprintf("  Simpson's concentration (D):       %.4f\n",   D_conc))
cat(sprintf("  Simpson's diversity (1–D):         %.4f\n",   1 - D_conc))
cat(sprintf("  Simpson's reciprocal (1/D):        %.2f\n",   1 / D_conc))

# ============================================================
# SECTION 10: TEMPORAL TRENDS — BOOTSTRAPPED PER-HA (all censuses)
# ============================================================
# KEY CORRECTION: use the same super-quadrat bootstrap approach
# (20 m base) for every census rather than crude total/plot area.

message("Computing bootstrapped per-ha temporal trends ...")

census_boot_list <- lapply(all_cids, function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0) return(NULL)
    sq   <- make_superquad(ct, BASE_Q_M)   # 20 m base quadrats, full grid
    rt   <- boot_mean_ci(sq$trees_ha)
    rb   <- boot_mean_ci(sq$ba_ha)
    meta <- census_meta[censusID == cid]
    data.table(
        censusID    = cid,
        census_num  = meta$census_num,
        census_year = meta$census_year,
        n_trees     = nrow(ct),
        n_species   = uniqueN(ct$sp),
        total_ba_m2 = round(sum(ct$ba_m2, na.rm = TRUE), 1),
        trees_ha    = rt["mean"], trees_ci_lo = rt["ci_lo"], trees_ci_hi = rt["ci_hi"],
        ba_ha       = rb["mean"], ba_ci_lo    = rb["ci_lo"], ba_ci_hi    = rb["ci_hi"]
    )
})
census_boot_dt <- rbindlist(census_boot_list)
setorder(census_boot_dt, census_num)

cat("\n", strrep("=", 72), "\n")
cat("  TEMPORAL TRENDS — BOOTSTRAPPED PER-HA (20 m quadrats, full grid)\n")
cat(strrep("=", 72), "\n")
print(census_boot_dt[, .(
    censusID, census_year, n_trees, n_species,
    trees_ha    = round(trees_ha,    1),
    trees_ci_lo = round(trees_ci_lo, 1),
    trees_ci_hi = round(trees_ci_hi, 1),
    ba_ha       = round(ba_ha,       2),
    ba_ci_lo    = round(ba_ci_lo,    2),
    ba_ci_hi    = round(ba_ci_hi,    2)
)])

# ============================================================
# SECTION 11: DIVERSITY INDICES ACROSS ALL CENSUSES
# ============================================================

census_div_list <- lapply(all_cids, function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0) return(NULL)
    pi_v <- ct[, .N, by = sp][, N / sum(N)]
    S    <- uniqueN(ct$sp)
    H    <- -sum(pi_v * log(pi_v))
    D_c  <- sum(pi_v^2)
    meta <- census_meta[censusID == cid]
    data.table(
        censusID    = cid,
        census_num  = meta$census_num,
        census_year = meta$census_year,
        S_richness  = S,
        H_shannon   = round(H,            4),
        J_evenness  = round(H / log(S),   4),
        D1mD        = round(1 - D_c,      4),
        D_recip     = round(1 / D_c,      2)
    )
})
census_div_dt <- rbindlist(census_div_list)
setorder(census_div_dt, census_num)

cat("\n", strrep("=", 64), "\n")
cat("  DIVERSITY INDICES ACROSS ALL CENSUSES\n")
cat(strrep("=", 64), "\n")
print(census_div_dt)

# ============================================================
# SECTION 12: COMPOSITIONAL CHANGE — WINNERS AND LOSERS
# ============================================================
# For each species:
#   (a) ΔTrees and ΔBA from Census 1 to Census 9.
#   (b) Linear regression of n_trees ~ census_num across all
#       censuses (species in ≥ 3 censuses); slope + p-value.
#   Classified as Increasing / Decreasing / Stable (p < 0.05).

message("Fitting species temporal trends ...")

# Per-species totals per census
spp_all_cens <- rec_tree[, .(
    n_trees = .N,
    ba_m2   = sum(ba_m2, na.rm = TRUE)
), by = .(sp, census_num, census_year)]
setorder(spp_all_cens, sp, census_num)

# First and last census counts
spp_c1 <- spp_all_cens[census_num == 1L,
    .(sp, n_c1 = n_trees, ba_c1 = ba_m2)]
spp_c9 <- spp_all_cens[census_num == N_CENSUSES,
    .(sp, n_c9 = n_trees, ba_c9 = ba_m2)]

# All species ever observed (some may be absent in C1 or C9)
all_sp <- unique(rec_tree[, .(sp, Latin)])
spp_change <- merge(
    merge(all_sp, spp_c1, by = "sp", all.x = TRUE),
    spp_c9, by = "sp", all.x = TRUE
)
spp_change[is.na(n_c1),  n_c1  := 0L]
spp_change[is.na(n_c9),  n_c9  := 0L]
spp_change[is.na(ba_c1), ba_c1 := 0.0]
spp_change[is.na(ba_c9), ba_c9 := 0.0]

spp_change[, delta_n       := n_c9  - n_c1]
spp_change[, delta_ba      := ba_c9 - ba_c1]
spp_change[, pct_change_n  := round(100 * delta_n  / pmax(n_c1,  1), 1)]
spp_change[, pct_change_ba := round(100 * delta_ba / pmax(ba_c1, 1e-9), 1)]

# Linear trend across all censuses (species with ≥ 3 appearances)
spp_ncens <- spp_all_cens[, .(n_cens = .N), by = sp]
spp_fit   <- spp_all_cens[sp %in% spp_ncens[n_cens >= 3L, sp]]

lm_trends <- spp_fit[, {
    mn <- tryCatch(lm(n_trees ~ census_num), error = function(e) NULL)
    mb <- tryCatch(lm(ba_m2   ~ census_num), error = function(e) NULL)
    list(
        slope_n  = if (!is.null(mn)) coef(mn)["census_num"]                   else NA_real_,
        p_n      = if (!is.null(mn)) summary(mn)$coefficients["census_num", 4] else NA_real_,
        slope_ba = if (!is.null(mb)) coef(mb)["census_num"]                   else NA_real_,
        p_ba     = if (!is.null(mb)) summary(mb)$coefficients["census_num", 4] else NA_real_
    )
}, by = sp]

spp_change <- merge(spp_change, lm_trends, by = "sp", all.x = TRUE)

spp_change[, trend_n := fcase(
    slope_n  >  0 & p_n  < 0.05, "Increasing",
    slope_n  <  0 & p_n  < 0.05, "Decreasing",
    !is.na(slope_n),               "Stable",
    default                      = "Insufficient data"
)]
spp_change[, trend_ba := fcase(
    slope_ba >  0 & p_ba < 0.05, "Increasing",
    slope_ba <  0 & p_ba < 0.05, "Decreasing",
    !is.na(slope_ba),              "Stable",
    default                      = "Insufficient data"
)]

n_new     <- nrow(spp_change[n_c1 == 0 & n_c9 > 0])
n_extinct <- nrow(spp_change[n_c1 > 0  & n_c9 == 0])

cat("\n", strrep("=", 64), "\n")
cat("  COMPOSITIONAL CHANGE: Census 1 → Census 9\n")
cat(strrep("=", 64), "\n")
cat(sprintf("  New arrivals  (absent C1, present C9): %d species\n", n_new))
cat(sprintf("  Local extinct (present C1, absent C9): %d species\n", n_extinct))
cat("\n  Linear trend summary — tree count:\n")
print(spp_change[, .N, by = trend_n][order(-N)])
cat("\n  Linear trend summary — basal area:\n")
print(spp_change[, .N, by = trend_ba][order(-N)])

N_SHOW <- 15

cat("\n  >> TOP 15 WINNERS — tree count (C9 − C1):\n")
print(head(spp_change[delta_n > 0][order(-delta_n),
    .(Latin, n_c1, n_c9, delta_n, pct_change_n, trend_n)], N_SHOW))
cat("\n  >> TOP 15 LOSERS — tree count (C9 − C1):\n")
print(head(spp_change[delta_n < 0][order(delta_n),
    .(Latin, n_c1, n_c9, delta_n, pct_change_n, trend_n)], N_SHOW))
cat("\n  >> TOP 15 WINNERS — basal area (m²):\n")
print(head(spp_change[delta_ba > 0][order(-delta_ba),
    .(Latin, ba_c1 = round(ba_c1, 2), ba_c9 = round(ba_c9, 2),
      delta_ba = round(delta_ba, 2), pct_change_ba, trend_ba)], N_SHOW))
cat("\n  >> TOP 15 LOSERS — basal area (m²):\n")
print(head(spp_change[delta_ba < 0][order(delta_ba),
    .(Latin, ba_c1 = round(ba_c1, 2), ba_c9 = round(ba_c9, 2),
      delta_ba = round(delta_ba, 2), pct_change_ba, trend_ba)], N_SHOW))

# ============================================================
# SECTION 13: SPECIES TURNOVER (Jaccard dissimilarity)
# ============================================================

spp_by_cens <- lapply(setNames(all_cids, all_cids),
    function(cid) unique(rec_tree[censusID == cid, sp]))

jaccard_d <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0L) return(NA_real_)
    1 - length(intersect(a, b)) / u
}

jacc_vs_c1 <- sapply(spp_by_cens, function(s) jaccard_d(spp_by_cens[[1]], s))
jacc_consec <- sapply(seq_len(N_CENSUSES - 1L),
    function(i) jaccard_d(spp_by_cens[[i]], spp_by_cens[[i + 1L]]))

turnover_dt <- data.table(
    censusID    = all_cids,
    census_year = CENSUS_YEARS,
    jacc_vs_C1  = round(jacc_vs_c1, 4),
    jacc_consec = c(NA_real_, round(jacc_consec, 4))
)

cat("\n", strrep("=", 64), "\n")
cat("  SPECIES TURNOVER (Jaccard dissimilarity; 0 = identical)\n")
cat(strrep("=", 64), "\n")
print(turnover_dt)

# ============================================================
# SECTION 14: DEMOGRAPHIC RATES (trees, consecutive censuses)
# ============================================================
# Annual mortality:  m = (1 − (N_surv/N_start)^(1/Δt)) × 100 %
# Annual recruitment: r = (N_recruits / (N_start × Δt)) × 100 %

demog_list <- lapply(seq_len(N_CENSUSES - 1L), function(i) {
    ids_i   <- rec_tree[censusID == all_cids[i],      unique(treeID)]
    ids_ip1 <- rec_tree[censusID == all_cids[i + 1L], unique(treeID)]
    n_start  <- length(ids_i)
    n_surv   <- sum(ids_i %in% ids_ip1)
    n_rec    <- length(ids_ip1) - n_surv
    dt_yr    <- CENSUS_YEARS[i + 1L] - CENSUS_YEARS[i]
    mort     <- (1 - (n_surv / n_start)^(1 / dt_yr)) * 100
    rec_rate <- (n_rec / (n_start * dt_yr)) * 100
    data.table(
        pair        = sprintf("C%d→C%d", i, i + 1L),
        year_start  = CENSUS_YEARS[i],
        year_end    = CENSUS_YEARS[i + 1L],
        dt_yr,
        n_start,
        n_survivors = n_surv,
        n_dead      = n_start - n_surv,
        n_recruits  = n_rec,
        mort_pct    = round(mort,     2),
        rec_pct     = round(rec_rate, 2)
    )
})
demog_dt <- rbindlist(demog_list)

cat("\n", strrep("=", 64), "\n")
cat("  DEMOGRAPHIC RATES (annual %, trees)\n")
cat(strrep("=", 64), "\n")
print(demog_dt)

# ============================================================
# SECTION 15: DBH SIZE CLASS DISTRIBUTION (Census 1 vs. Census 9)
# ============================================================

breaks_mm <- c(10, 20, 30, 50, 70, 100, 150, 200, 300, 500, Inf)
labels_mm <- c("10–20", "20–30", "30–50", "50–70", "70–100",
               "100–150", "150–200", "200–300", "300–500", "≥500")

stems_c1c9 <- rbind(
    rec_alive[censusID == FIRST_ID, .(censusID, dbh)],
    rec_alive[censusID == LAST_ID,  .(censusID, dbh)]
)
stems_c1c9 <- merge(stems_c1c9,
    census_meta[, .(censusID, census_year)], by = "censusID")
stems_c1c9[, dbh_class := cut(dbh, breaks = breaks_mm,
    labels = labels_mm, right = FALSE)]

dbh_dist <- stems_c1c9[!is.na(dbh_class), .N, by = .(census_year, dbh_class)]
dbh_dist[, tot := sum(N), by = census_year]
dbh_dist[, pct := round(100 * N / tot, 1)]

cat("\n", strrep("=", 64), "\n")
cat("  DBH SIZE CLASS — count (%) by census\n")
cat(strrep("=", 64), "\n")
wide_n <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "N",   fill = 0L)
wide_p <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "pct", fill = 0)
cat("  Counts:\n");      print(wide_n)
cat("  Proportions (%):\n"); print(wide_p)

# ============================================================
# SECTION 16: SPATIAL SUMMARY (100 m sectors)
# ============================================================

last_trees[, sec_x := floor(gx / 100) * 100]
last_trees[, sec_y := floor(gy / 100) * 100]
last_trees[sec_x >= PLOT_X_M, sec_x := PLOT_X_M - 100L]
last_trees[sec_y >= PLOT_Y_M, sec_y := PLOT_Y_M - 100L]

sector_obs <- last_trees[, .(
    n_trees = .N,
    ba_m2   = sum(ba_m2, na.rm = TRUE),
    n_sp    = uniqueN(sp)
), by = .(sec_x, sec_y)]

# Complete 10 × 5 grid of 1-ha sectors
full_sectors <- CJ(sec_x = seq(0L, 900L, 100L),
                   sec_y = seq(0L, 400L, 100L))
sector_sum <- merge(full_sectors, sector_obs, by = c("sec_x", "sec_y"), all.x = TRUE)
sector_sum[is.na(n_trees), n_trees := 0L]
sector_sum[is.na(ba_m2),   ba_m2   := 0.0]
sector_sum[is.na(n_sp),    n_sp    := 0L]
sector_sum[, trees_ha := n_trees / 1]   # 100 × 100 m = 1 ha
sector_sum[, ba_ha    := ba_m2   / 1]

cat("\n", strrep("=", 64), "\n")
cat("  SPATIAL HETEROGENEITY — 100 × 100 m sectors (1 ha each)\n")
cat(strrep("=", 64), "\n")
cat(sprintf("  Trees/ha:      mean = %6.1f, SD = %5.1f, CV = %4.1f%%\n",
    mean(sector_sum$trees_ha), sd(sector_sum$trees_ha),
    100 * sd(sector_sum$trees_ha) / mean(sector_sum$trees_ha)))
cat(sprintf("  BA/ha (m²):    mean = %6.2f, SD = %5.2f, CV = %4.1f%%\n",
    mean(sector_sum$ba_ha), sd(sector_sum$ba_ha),
    100 * sd(sector_sum$ba_ha) / mean(sector_sum$ba_ha)))
cat(sprintf("  Spp/sector:    mean = %6.1f, SD = %5.1f\n",
    mean(sector_sum$n_sp), sd(sector_sum$n_sp)))

# ============================================================
# SECTION 17: PLOTS (nine figures saved to ./figures/)
# ============================================================

message("Generating plots ...")
dir.create("figures", showWarnings = FALSE)

# Shared theme
theme_bci <- function(base = 11) {
    theme_classic(base_size = base) +
    theme(
        plot.title       = element_text(face = "bold",  size = base),
        plot.subtitle    = element_text(size = base - 1, color = "grey40"),
        axis.title       = element_text(size = base - 1),
        axis.text        = element_text(size = base - 2),
        legend.title     = element_text(size = base - 2),
        legend.text      = element_text(size = base - 2),
        strip.text       = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.35)
    )
}
YR_BRK <- CENSUS_YEARS

# -------------------------------------------------------------------
# PLOT 1: Temporal trends — trees/ha, BA/ha, richness, diversity
# -------------------------------------------------------------------

p1a <- ggplot(census_boot_dt, aes(x = census_year)) +
    geom_ribbon(aes(ymin = trees_ci_lo, ymax = trees_ci_hi),
                fill = "#1976D2", alpha = .20) +
    geom_line(aes(y = trees_ha),  color = "#1565C0", linewidth = .9) +
    geom_point(aes(y = trees_ha), color = "#1565C0", size = 2.2) +
    scale_x_continuous(breaks = YR_BRK) +
    labs(title = "Tree density", x = NULL,
         y = expression(Trees~ha^{-1})) + theme_bci()

p1b <- ggplot(census_boot_dt, aes(x = census_year)) +
    geom_ribbon(aes(ymin = ba_ci_lo, ymax = ba_ci_hi),
                fill = "#388E3C", alpha = .20) +
    geom_line(aes(y = ba_ha),  color = "#2E7D32", linewidth = .9) +
    geom_point(aes(y = ba_ha), color = "#2E7D32", size = 2.2) +
    scale_x_continuous(breaks = YR_BRK) +
    labs(title = "Basal area", x = NULL,
         y = expression(BA~(m^2~ha^{-1}))) + theme_bci()

p1c <- ggplot(census_boot_dt, aes(x = census_year, y = n_species)) +
    geom_line(color = "#7B1FA2", linewidth = .9) +
    geom_point(color = "#7B1FA2", size = 2.2) +
    scale_x_continuous(breaks = YR_BRK) +
    labs(title = "Species richness", x = "Census year",
         y = "No. species") + theme_bci()

div_long <- melt(census_div_dt, id.vars = "census_year",
    measure.vars = c("H_shannon", "J_evenness", "D1mD"),
    variable.name = "Index", value.name = "Value")
div_long[, Index := factor(Index,
    c("H_shannon", "J_evenness", "D1mD"),
    c("Shannon H'", "Pielou J",  "Simpson 1–D"))]

p1d <- ggplot(div_long, aes(x = census_year, y = Value, color = Index)) +
    geom_line(linewidth = .8) + geom_point(size = 1.8) +
    scale_color_brewer(palette = "Dark2") +
    scale_x_continuous(breaks = YR_BRK) +
    labs(title = "Diversity indices", x = "Census year",
         y = "Index value", color = NULL) +
    theme_bci() + theme(legend.position = "bottom")

plot1 <- (p1a | p1b) / (p1c | p1d) +
    plot_annotation(
        title = "BCI 50-ha Plot — Temporal Trends  (bootstrap 95% CI, 20 m quadrats)",
        theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
ggsave("figures/BCI_01_temporal_trends.pdf", plot1, width = 12, height = 9)
message("  Saved: BCI_01_temporal_trends.pdf")

# -------------------------------------------------------------------
# PLOT 2: Rank-abundance curve (last census, log-log)
# -------------------------------------------------------------------

rank_dt  <- copy(spp_sum)[order(-n_trees)][, rank := .I]
top20_lbl <- rank_dt[rank <= 20, .(rank, n_trees, lbl = short_latin(Latin))]

p2 <- ggplot(rank_dt, aes(x = rank, y = n_trees)) +
    geom_line(color = "#37474F", linewidth = .6) +
    geom_point(size = .7, color = "#37474F") +
    geom_text_repel(data = top20_lbl, aes(label = lbl),
        size = 2.6, max.overlaps = 25, segment.size = .25,
        color = "#B71C1C", fontface = "italic") +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    labs(
        title    = sprintf("Rank–abundance curve (Census %d)", N_CENSUSES),
        subtitle = "Top 20 species labelled",
        x        = "Species rank (log scale)",
        y        = "No. trees (log scale)"
    ) + theme_bci()
ggsave("figures/BCI_02_rank_abundance.pdf", p2, width = 9, height = 6)
message("  Saved: BCI_02_rank_abundance.pdf")

# -------------------------------------------------------------------
# PLOT 3: DBH size-class distribution — Census 1 vs. Census 9
# -------------------------------------------------------------------

p3 <- ggplot(dbh_dist,
    aes(x = dbh_class, y = N, fill = factor(census_year), group = census_year)) +
    geom_col(position = "dodge", color = "white", linewidth = .25) +
    scale_fill_manual(
        values = c("#81C784", "#1B5E20"),
        labels = c(paste0("Census 1 (", CENSUS_YEARS[1L],         ")"),
                   paste0("Census 9 (", CENSUS_YEARS[N_CENSUSES], ")"))
    ) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "DBH size-class distribution — Census 1 vs. Census 9",
         x = "DBH class (mm)", y = "No. stems", fill = NULL) +
    theme_bci() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("figures/BCI_03_dbh_distribution.pdf", p3, width = 9, height = 5)
message("  Saved: BCI_03_dbh_distribution.pdf")

# -------------------------------------------------------------------
# PLOT 4: Bootstrap CI by quadrat size (trees/ha and BA/ha)
# -------------------------------------------------------------------

p4a <- ggplot(perha_dt, aes(x = quad_ha, y = trees_ha)) +
    geom_ribbon(aes(ymin = trees_ci_lo, ymax = trees_ci_hi),
                fill = "#FF8F00", alpha = .30) +
    geom_line(color = "#E65100", linewidth = .9) +
    geom_point(aes(size = n_quads), color = "#E65100", show.legend = FALSE) +
    geom_text(aes(label = sprintf("%dm\n(n=%d)", quad_m, n_quads)),
              vjust = -0.6, size = 2.5, color = "#6D4C41") +
    scale_x_log10(labels = label_number(accuracy = .001)) +
    labs(
        title    = expression(Trees~ha^{-1}~"by quadrat size"),
        subtitle = "Mean ± 95% bootstrap CI  |  point size ∝ no. of quadrats",
        x        = "Quadrat area (ha, log scale)",
        y        = expression(Trees~ha^{-1})
    ) + theme_bci()

p4b <- ggplot(perha_dt, aes(x = quad_ha, y = ba_ha)) +
    geom_ribbon(aes(ymin = ba_ci_lo, ymax = ba_ci_hi),
                fill = "#6A1B9A", alpha = .30) +
    geom_line(color = "#4A148C", linewidth = .9) +
    geom_point(aes(size = n_quads), color = "#4A148C", show.legend = FALSE) +
    geom_text(aes(label = sprintf("%dm\n(n=%d)", quad_m, n_quads)),
              vjust = -0.6, size = 2.5, color = "#4A148C") +
    scale_x_log10(labels = label_number(accuracy = .001)) +
    labs(
        title    = expression(Basal~area~ha^{-1}~"by quadrat size"),
        subtitle = "Mean ± 95% bootstrap CI  |  point size ∝ no. of quadrats",
        x        = "Quadrat area (ha, log scale)",
        y        = expression(BA~(m^2~ha^{-1}))
    ) + theme_bci()

plot4 <- p4a | p4b
ggsave("figures/BCI_04_perha_by_quadsize.pdf", plot4, width = 12, height = 5.5)
message("  Saved: BCI_04_perha_by_quadsize.pdf")

# -------------------------------------------------------------------
# PLOT 5: Compositional winners and losers (tree count + basal area)
# -------------------------------------------------------------------

build_wl <- function(n_show = 15, var_delta, var_trend) {
    inc <- spp_change[get(var_trend) == "Increasing"]
    dec <- spp_change[get(var_trend) == "Decreasing"]
    setorderv(inc, var_delta, order = -1L)
    setorderv(dec, var_delta, order =  1L)
    inc <- inc[seq_len(min(.N, n_show))]
    dec <- dec[seq_len(min(.N, n_show))]
    dt  <- rbind(
        inc[, .(Latin, val = get(var_delta), dir = "Increasing")],
        dec[, .(Latin, val = get(var_delta), dir = "Decreasing")]
    )
    dt[, lbl := reorder(short_latin(Latin), val)]
    dt
}

wl_n  <- build_wl(var_delta = "delta_n",  var_trend = "trend_n")
wl_ba <- build_wl(var_delta = "delta_ba", var_trend = "trend_ba")

p5a <- ggplot(wl_n, aes(x = val, y = lbl, fill = dir)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = .4) +
    scale_fill_manual(values = c("Increasing" = "#43A047", "Decreasing" = "#E53935")) +
    labs(title    = "Abundance change — top 15 winners / losers",
         subtitle = "Species with significant linear trend (p < 0.05)",
         x        = "ΔTrees (Census 9 − Census 1)", y = NULL) +
    theme_bci() + theme(axis.text.y = element_text(face = "italic", size = 8))

p5b <- ggplot(wl_ba, aes(x = val, y = lbl, fill = dir)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = .4) +
    scale_fill_manual(values = c("Increasing" = "#43A047", "Decreasing" = "#E53935")) +
    labs(title    = "Basal area change — top 15 winners / losers",
         subtitle = "Species with significant linear trend (p < 0.05)",
         x        = expression(Delta~"BA (m², Census 9 − Census 1)"), y = NULL) +
    theme_bci() + theme(axis.text.y = element_text(face = "italic", size = 8))

plot5 <- p5a | p5b
ggsave("figures/BCI_05_winners_losers.pdf", plot5, width = 14, height = 9)
message("  Saved: BCI_05_winners_losers.pdf")

# -------------------------------------------------------------------
# PLOT 6: Spatial heatmaps — trees/ha, BA/ha, species richness
# -------------------------------------------------------------------

mk_tile <- function(fill_var, ttl, palette = "plasma", lbl = "") {
    ggplot(sector_sum, aes(x = sec_x, y = sec_y,
                           fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = .3) +
    scale_fill_viridis_c(option = palette, name = lbl) +
    coord_fixed(expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 900, 200)) +
    scale_y_continuous(breaks = seq(0, 400, 100)) +
    labs(title = ttl, x = "x (m)", y = "y (m)") + theme_bci()
}

p6a <- mk_tile("trees_ha", "Tree density",    "plasma",  "Trees\nha⁻¹")
p6b <- mk_tile("ba_ha",    "Basal area",       "magma",   expression(m^2~ha^{-1}))
p6c <- mk_tile("n_sp",     "Species richness", "viridis", "Species\nper ha")

plot6 <- p6a / p6b / p6c +
    plot_annotation(
        title = "BCI — Spatial distribution across 100 × 100 m sectors",
        theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
ggsave("figures/BCI_06_spatial_maps.pdf", plot6, width = 11, height = 14)
message("  Saved: BCI_06_spatial_maps.pdf")

# -------------------------------------------------------------------
# PLOT 7: Demographic rates (annual mortality and recruitment)
# -------------------------------------------------------------------

demog_long <- melt(demog_dt,
    id.vars      = c("pair", "year_start", "year_end"),
    measure.vars = c("mort_pct", "rec_pct"),
    variable.name = "rate_type", value.name = "rate")
demog_long[, year_mid  := (year_start + year_end) / 2]
demog_long[, rate_type := factor(rate_type,
    c("mort_pct", "rec_pct"), c("Mortality", "Recruitment"))]

p7 <- ggplot(demog_long, aes(x = year_mid, y = rate, color = rate_type)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_line(linewidth = .9) + geom_point(size = 2.2) +
    scale_color_manual(
        values = c("Mortality" = "#E53935", "Recruitment" = "#43A047")) +
    scale_x_continuous(breaks = YR_BRK) +
    labs(title = "Annual demographic rates — BCI trees",
         x = "Census mid-year", y = "Rate (% per year)", color = NULL) +
    theme_bci() + theme(legend.position = "bottom")
ggsave("figures/BCI_07_demographic_rates.pdf", p7, width = 8, height = 5)
message("  Saved: BCI_07_demographic_rates.pdf")

# -------------------------------------------------------------------
# PLOT 8: Lifeform composition over time
# -------------------------------------------------------------------

if ("Lifeform_RFoster" %in% names(rec_tree)) {
    lf_time <- rec_tree[!is.na(Lifeform_RFoster),
        .(n_trees = .N), by = .(census_year, Lifeform_RFoster)]

    p8 <- ggplot(lf_time,
        aes(x = factor(census_year), y = n_trees,
            fill = Lifeform_RFoster)) +
        geom_col(position = "fill", color = "white", linewidth = .2) +
        scale_y_continuous(labels = percent_format()) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = "Lifeform composition over time (Lifeform_RFoster)",
             x = "Census year", y = "Proportion", fill = "Lifeform") +
        theme_bci() + theme(legend.position = "right")
    ggsave("figures/BCI_08_lifeform_composition.pdf", p8, width = 9, height = 5)
    message("  Saved: BCI_08_lifeform_composition.pdf")
} else {
    message("  Plot 8 skipped — Lifeform_RFoster column absent.")
}

# -------------------------------------------------------------------
# PLOT 9: Per-species abundance trajectories — top 20 species
# -------------------------------------------------------------------

top20_sp <- spp_sum[order(-n_trees), sp][seq_len(min(20L, nrow(spp_sum)))]
traj_dt  <- spp_all_cens[sp %in% top20_sp]
traj_dt  <- merge(traj_dt, unique(rec_tree[, .(sp, Latin)]), by = "sp")
traj_dt[, short := short_latin(Latin)]

# Relative abundance vs. census total
cens_tot <- census_boot_dt[, .(census_year, cens_total = n_trees)]
traj_dt  <- merge(traj_dt, cens_tot, by = "census_year")
traj_dt[, rel_abund := n_trees / cens_total]

pal20 <- scales::hue_pal()(uniqueN(traj_dt$short))

p9a <- ggplot(traj_dt,
    aes(x = census_year, y = n_trees, color = short, group = short)) +
    geom_line(linewidth = .7) + geom_point(size = 1.5) +
    scale_x_continuous(breaks = YR_BRK) +
    scale_y_continuous(labels = label_comma()) +
    scale_color_manual(values = pal20) +
    labs(title = "Top 20 species — absolute abundance",
         x = "Census year", y = "No. trees", color = NULL) +
    theme_bci() + theme(legend.position  = "right",
                        legend.text      = element_text(size = 7, face = "italic"),
                        legend.key.height = unit(.55, "lines"))

p9b <- ggplot(traj_dt,
    aes(x = census_year, y = rel_abund, color = short, group = short)) +
    geom_line(linewidth = .7) + geom_point(size = 1.5) +
    scale_x_continuous(breaks = YR_BRK) +
    scale_y_continuous(labels = percent_format(accuracy = .1)) +
    scale_color_manual(values = pal20) +
    labs(title = "Top 20 species — relative abundance",
         x = "Census year", y = "Proportion of all trees", color = NULL) +
    theme_bci() + theme(legend.position  = "right",
                        legend.text      = element_text(size = 7, face = "italic"),
                        legend.key.height = unit(.55, "lines"))

plot9 <- p9a / p9b +
    plot_annotation(
        title = "BCI — Abundance trajectories of top 20 species",
        theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
ggsave("figures/BCI_09_species_trajectories.pdf", plot9, width = 12, height = 12)
message("  Saved: BCI_09_species_trajectories.pdf")

# ============================================================
# DONE
# ============================================================
cat("\n", strrep("=", 64), "\n")
cat("  ANALYSIS COMPLETE\n")
cat("  Plots saved to ./figures/  (BCI_01 – BCI_09)\n")
cat(strrep("=", 64), "\n")
