rm(list = ls())

# ============================================================
# GENERAL PLOT INFORMATION: BCI 50-ha Plot
# ============================================================
# BCI permanent plot: 1000 m (E–W) × 500 m (N–S) = 50 ha
# Base quadrats:      20 × 20 m (0.04 ha each); 1,250 total
#
# NOTE: These are basic descriptive summaries intended to give
# an overall view of forest structure and composition across
# the full observation period. Understanding temporal
# patterns and their ecological drivers requires more developed
# statistical analyses beyond what is shown here.
# ============================================================

library(data.table)
library(boot)

# ============================================================
# SECTION 1: CONSTANTS
# ============================================================

PLOT_X_M <- 1000L # plot width  (m), E–W
PLOT_Y_M <- 500L # plot height (m), N–S
PLOT_AREA_HA <- PLOT_X_M * PLOT_Y_M / 1e4 # = 50 ha
BASE_Q_M <- 20L # base quadrat side (m)
N_CENSUSES <- 9L
BOOT_R <- 4999L
set.seed(42)

# Quadrat sizes (m, square) that tile 1000 × 500 exactly with NO gaps.
# Requirement: qs must divide both PLOT_X_M and PLOT_Y_M.
# GCD(1000, 500) = 500; valid sizes >= 20 m: 20, 25, 50, 100, 125, 250.
QUAD_SIZES_M <- c(20, 25, 50, 100, 125, 250)

# Convenience: extract genus + species epithet from full Latin name
short_latin <- function(x) gsub("^(\\S+\\s+\\S+).*", "\\1", trimws(x))

# ============================================================
# SECTION 2: DATA LOADING
# ============================================================
message("Loading census data ...")
bci_nums <- as.character(seq_len(N_CENSUSES))
census_list <- lapply(bci_nums, function(n) {
    fp <- sprintf("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem%s.Rdata", n)
    if (!file.exists(fp)) stop("Missing file: ", fp)
    e <- new.env(parent = emptyenv())
    load(fp, envir = e)
    get(paste0("bci.stem", n), envir = e)
})
names(census_list) <- paste0("bci.stem", bci_nums)
rec <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rec <- rec[!is.na(quadrat)]
rm(census_list, bci_nums)
gc()

# ---- Species table --------------------------------------------------
load("./BCI_stem_reconstruction/DATA/RTABLES/bci.spptable.rdata")
bci.spptable <- as.data.table(bci.spptable)

# Build full Latin name (include infraspecific epithets when present)
if (all(c("InfraspecificRank", "InfraspecificEpithet") %in% names(bci.spptable))) {
    bci.spptable[, Latin := fifelse(
        !is.na(InfraspecificRank) &
            nchar(trimws(as.character(InfraspecificRank))) > 0 &
            !is.na(InfraspecificEpithet) &
            nchar(trimws(as.character(InfraspecificEpithet))) > 0,
        paste(Latin, InfraspecificRank, InfraspecificEpithet),
        Latin
    )]
}
lf_cols <- intersect(c("Lifeform_RFoster", "Lifeform_RPerez_SAguilar"), names(bci.spptable))
keep_spp <- c("sp", "Latin", lf_cols)
bci.spptable <- bci.spptable[, keep_spp, with = FALSE]
rec <- merge(rec, bci.spptable, by = "sp", all.x = TRUE)

# ---- Basal area (m²) per stem ---------------------------------------
rec[, ba_m2 := pi / 4 * (dbh / 1000)^2]

# ---- Census metadata lookup -----------------------------------------
get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

CENSUS_YEARS <- rec[, year := year(ExactDate)][
    , .(CensusID, year)
][
    , .(year = get_mode(year)),
    by = CensusID
][
    order(CensusID),
    year
]

census_meta <- data.table(
    censusID    = paste0("bci.stem", seq_len(N_CENSUSES)),
    census_num  = seq_len(N_CENSUSES),
    census_year = CENSUS_YEARS
)

# ============================================================
# SECTION 3: TREE-LEVEL AGGREGATION (alive stems only)
# ============================================================
# • Keep only Rstatus == "A" (alive resolved status).
# • One row per treeID per census: BA = sum of all stem BAs,
#   gx/gy = mean position of alive stems (used for super-quadrat
#   assignment via make_superquad.

rec_alive <- rec[Rstatus == "A" & !is.na(gx) & !is.na(gy)]

by_cols <- c("censusID", "CensusID", "treeID", "sp", "Latin", lf_cols)

rec_tree <- rec_alive[, .(
    ba_m2   = sum(ba_m2, na.rm = TRUE),
    n_stems = .N,
    gx      = mean(gx, na.rm = TRUE),
    gy      = mean(gy, na.rm = TRUE)
), by = by_cols]

# Clip to strict plot interior (edge-case guard)
rec_tree[gx < 0, gx := 0]
rec_tree[gx >= PLOT_X_M, gx := PLOT_X_M - 1e-3]
rec_tree[gy < 0, gy := 0]
rec_tree[gy >= PLOT_Y_M, gy := PLOT_Y_M - 1e-3]

rec_tree <- merge(rec_tree, census_meta, by = "censusID", all.x = TRUE)
setorder(rec_tree, census_num, treeID)

# Validate census_num merge
if (any(is.na(rec_tree$census_num))) {
    warning("Some treeIDs have missing census_num; dropping affected rows.")
    rec_tree <- rec_tree[!is.na(census_num)]
}

all_cids <- census_meta$censusID
# Exclude 1982 census
FIRST_ID <- all_cids[2]
LAST_ID <- all_cids[N_CENSUSES]

last_trees <- rec_tree[censusID == LAST_ID]
last_stems <- rec_alive[censusID == LAST_ID]

total_trees <- nrow(last_trees)
total_stems <- nrow(last_stems)
total_ba <- sum(last_trees$ba_m2, na.rm = TRUE)
n_species <- uniqueN(last_trees$sp)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# ---- Bootstrap mean + 95% HDI ----
boot_mean_ci <- function(x, R = BOOT_R, ci_level = 0.95) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    if (n < 4) {
        se <- if (n > 1) sd(x) / sqrt(n) else 0
        return(c(mean = m, ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se))
    }
    b <- boot::boot(x, function(d, i) mean(d[i]), R = R)
    boots <- as.numeric(b$t)
    ci <- HDInterval::hdi(boots, credMass = ci_level)
    if (!is.numeric(ci) || length(ci) != 2 || any(!is.finite(ci))) {
        se <- sd(x) / sqrt(n)
        ci <- c(m - 1.96 * se, m + 1.96 * se)
    }
    c(mean = m, ci_lo = ci[1], ci_hi = ci[2])
}

# ---- Build complete super-quadrat grid and aggregate ----
# For quadrat size qs_m × qs_m:
#   • Trees assigned via floor(gx/qs_m), floor(gy/qs_m).
#   • COMPLETE grid created (n_x × n_y cells); empty cells = 0.
#   • Counts and BA scaled to per-ha.
# qs_m must divide both PLOT_X_M and PLOT_Y_M exactly.
make_superquad <- function(tree_dt, qs_m) {
    dt <- copy(tree_dt)
    n_x <- as.integer(PLOT_X_M / qs_m)
    n_y <- as.integer(PLOT_Y_M / qs_m)
    qs_ha <- (qs_m * qs_m) / 1e4

    dt[, sxi := pmin(as.integer(floor(gx / qs_m)), n_x - 1L)]
    dt[, syi := pmin(as.integer(floor(gy / qs_m)), n_y - 1L)]
    dt[, sqid := sxi * n_y + syi]

    sq_obs <- dt[, .(
        n_trees = .N,
        ba_m2   = sum(ba_m2, na.rm = TRUE),
        n_sp    = uniqueN(sp)
    ), by = sqid]

    # Complete grid (zero-fill empties so the bootstrap is unbiased)
    sq_full <- merge(
        data.table(sqid = 0L:(n_x * n_y - 1L)),
        sq_obs,
        by = "sqid", all.x = TRUE
    )
    sq_full[is.na(n_trees), n_trees := 0L]
    sq_full[is.na(ba_m2), ba_m2 := 0.0]
    sq_full[is.na(n_sp), n_sp := 0L]
    sq_full[, trees_ha := n_trees / qs_ha]
    sq_full[, ba_ha := ba_m2 / qs_ha]
    sq_full
}

# ============================================================
# SECTION 4: OVERALL PLOT SUMMARY (last census)
# ============================================================

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  BCI 50-ha PLOT — Census %d (alive trees, DBH ≥ 10 mm)\n", N_CENSUSES))
cat(strrep("=", 64), "\n")
cat(sprintf(
    "  Plot dimensions:                  %d × %d m = %g ha\n",
    PLOT_X_M, PLOT_Y_M, PLOT_AREA_HA
))
cat(sprintf("  Total trees (treeIDs, alive):     %d\n", total_trees))
cat(sprintf("  Total stems (stemIDs, alive):     %d\n", total_stems))
cat(sprintf("  Species richness:                 %d\n", n_species))
cat(sprintf("  Total basal area (m²):            %.2f\n", total_ba))
cat(sprintf("  Crude trees ha⁻¹ (total/area):   %.1f\n", total_trees / PLOT_AREA_HA))
cat(sprintf("  Crude BA ha⁻¹ (m²/ha):           %.2f\n", total_ba / PLOT_AREA_HA))

# ============================================================
# SECTION 5: SPECIES-LEVEL SUMMARY (last census)
# ============================================================

spp_sum <- last_trees[, .(
    n_trees     = .N,
    total_ba_m2 = sum(ba_m2, na.rm = TRUE),
    mean_ba_m2  = mean(ba_m2, na.rm = TRUE)
), by = .(sp, Latin)]
spp_sum[, rel_abund_pct := round(100 * n_trees / total_trees, 3)]
spp_sum[, rel_ba_pct := round(100 * total_ba_m2 / total_ba, 3)]
spp_sum[, IV := rel_abund_pct + rel_ba_pct] # Importance Value (max 200)

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY ABUNDANCE (Census 9)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-n_trees),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct)
], 10))

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY BASAL AREA (m²)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-total_ba_m2),
    .(Latin, total_ba_m2, rel_ba_pct, n_trees, rel_abund_pct)
], 10))

cat("\n", strrep("-", 64), "\n")
cat("  TOP 10 SPECIES BY IMPORTANCE VALUE (abundance + BA, max = 200)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
], 10))

# ============================================================
# SECTION 6: ABUNDANCE / BA THRESHOLDS
# ============================================================

avg_n <- mean(spp_sum$n_trees)
avg_ba <- mean(spp_sum$total_ba_m2)
q75_n <- quantile(spp_sum$n_trees, 0.75)
q75_ba <- quantile(spp_sum$total_ba_m2, 0.75)

cat("\n", strrep("-", 64), "\n")
cat("  ABUNDANCE / BA THRESHOLDS\n")
cat(strrep("-", 64), "\n")
cat(sprintf("  Mean trees/species:       %6.1f  | 75th pct: %6.0f\n", avg_n, q75_n))
cat(sprintf("  Mean BA/species (m²):   %8.4f  | 75th pct: %8.4f\n", avg_ba, q75_ba))
cat("\n  >> Species above MEAN in both abundance & BA:\n")
print(spp_sum[n_trees > avg_n & total_ba_m2 > avg_ba][
    order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
])
cat("\n  >> Species above 75th pct in both abundance & BA:\n")
print(spp_sum[n_trees > q75_n & total_ba_m2 > q75_ba][
    order(-IV),
    .(Latin, n_trees, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
])

# ============================================================
# SECTION 7: LIFEFORM BREAKDOWN
# ============================================================

for (lfc in lf_cols) {
    lf <- last_trees[!is.na(get(lfc)), .(
        n_trees   = .N,
        n_species = uniqueN(sp),
        ba_m2     = sum(ba_m2, na.rm = TRUE)
    ), by = c(lfc)][order(-n_trees)]
    lf[, rel_n_pct := round(100 * n_trees / total_trees, 1)]
    lf[, rel_ba_pct := round(100 * ba_m2 / total_ba, 1)]
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
    # qs <- 20L
    qs_ha <- (qs * qs) / 1e4
    message(sprintf("  %d × %d m (%.4f ha) ...", qs, qs, qs_ha))
    sq <- make_superquad(last_trees, qs)
    rt <- boot_mean_ci(sq$trees_ha)
    rb <- boot_mean_ci(sq$ba_ha)
    rs <- boot_mean_ci(sq$n_sp)
    data.table(
        quad_m = qs, quad_ha = qs_ha, n_quads = nrow(sq),
        trees_ha = rt["mean"], trees_ci_lo.lower = rt["ci_lo.lower"], trees_ci_hi.upper = rt["ci_hi.upper"],
        ba_ha = rb["mean"], ba_ci_lo.lower = rb["ci_lo.lower"], ba_ci_hi.upper = rb["ci_hi.upper"],
        sp_per_quad = rs["mean"], sp_ci_lo.lower = rs["ci_lo.lower"], sp_ci_hi.upper = rs["ci_hi.upper"]
    )
})
perha_dt <- rbindlist(perha_rows)

cat("\n", strrep("=", 72), "\n")
cat(sprintf("  PER-HECTARE ESTIMATES BY QUADRAT SIZE  (boot 95%% CI, R = %d)\n", BOOT_R))
cat(strrep("=", 72), "\n")
print(perha_dt[, .(
    quad_m, quad_ha, n_quads,
    trees_ha = round(trees_ha, 1),
    trees_ci_lo = round(trees_ci_lo.lower, 1),
    trees_ci_hi = round(trees_ci_hi.upper, 1),
    ba_ha = round(ba_ha, 2),
    ba_ci_lo = round(ba_ci_lo.lower, 2),
    ba_ci_hi = round(ba_ci_hi.upper, 2)
)])

# ============================================================
# SECTION 9: DIVERSITY INDICES (last census)
# ============================================================

# pi_vals <- spp_sum$n_trees / total_trees
# H_shannon <- -sum(pi_vals * log(pi_vals))
# D_conc <- sum(pi_vals^2)
# J_evenness <- H_shannon / log(n_species)

pi_vals <- spp_sum$n_trees / sum(spp_sum$n_trees)
H_shannon <- -sum(pi_vals[pi_vals > 0] *
    log(pi_vals[pi_vals > 0]))
D_conc <- sum(pi_vals^2)
n_species <- sum(spp_sum$n_trees > 0)
J_evenness <- H_shannon / log(n_species)

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  DIVERSITY INDICES — Census %d\n", N_CENSUSES))
cat(strrep("=", 64), "\n")
cat(sprintf("  Species richness (S):              %d\n", n_species))
cat(sprintf("  Shannon-Wiener (H'):               %.4f\n", H_shannon))
cat(sprintf("  Max H' [= ln(S)]:                  %.4f\n", log(n_species)))
cat(sprintf("  Pielou's evenness (J = H'/ln S):   %.4f\n", J_evenness))
cat(sprintf("  Simpson's concentration (D):       %.4f\n", D_conc))
cat(sprintf("  Simpson's diversity (1–D):         %.4f\n", 1 - D_conc))
cat(sprintf("  Simpson's reciprocal (1/D):        %.2f\n", 1 / D_conc))

# ============================================================
# SECTION 10: TEMPORAL TRENDS — BOOTSTRAPPED PER-HA (all censuses)
# ============================================================
# KEY CORRECTION: use the same super-quadrat bootstrap approach
# (20 m base) for every census rather than crude total/plot area.

message("Computing bootstrapped per-ha temporal trends ...")

census_boot_list <- lapply(all_cids[2:9], function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0) {
        return(NULL)
    }
    sq <- make_superquad(ct, BASE_Q_M) # 20 m base quadrats, full grid
    rt <- boot_mean_ci(sq$trees_ha)
    rb <- boot_mean_ci(sq$ba_ha)
    meta <- census_meta[censusID == cid]
    data.table(
        censusID    = cid,
        census_num  = meta$census_num,
        census_year = meta$census_year,
        n_trees     = nrow(ct),
        n_species   = uniqueN(ct$sp),
        total_ba_m2 = round(sum(ct$ba_m2, na.rm = TRUE), 1),
        trees_ha    = rt["mean"], trees_ci_lo.lower = rt["ci_lo.lower"], trees_ci_hi.upper = rt["ci_hi.upper"],
        ba_ha       = rb["mean"], ba_ci_lo.lower    = rb["ci_lo.lower"], ba_ci_hi.upper    = rb["ci_hi.upper"]
    )
})
census_boot_dt <- rbindlist(census_boot_list)
setorder(census_boot_dt, census_num)

cat("\n", strrep("=", 72), "\n")
cat("  TEMPORAL TRENDS — BOOTSTRAPPED PER-HA (20 m quadrats, full grid)\n")
cat(strrep("=", 72), "\n")
print(census_boot_dt[, .(
    censusID, census_year, n_trees, n_species,
    trees_ha = round(trees_ha, 1),
    trees_ci_lo.lower = round(trees_ci_lo.lower, 1),
    trees_ci_hi.upper = round(trees_ci_hi.upper, 1),
    ba_ha = round(ba_ha, 2),
    ba_ci_lo.lower = round(ba_ci_lo.lower, 2),
    ba_ci_hi.upper = round(ba_ci_hi.upper, 2)
)])

# ============================================================
# SECTION 11: DIVERSITY INDICES ACROSS ALL CENSUSES
# ============================================================

census_div_list <- lapply(all_cids[2:9], function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0) {
        return(NULL)
    }
    pi_v <- ct[, .N, by = sp][, N / sum(N)]
    S <- uniqueN(ct$sp)
    H <- -sum(pi_v * log(pi_v))
    D_c <- sum(pi_v^2)
    meta <- census_meta[censusID == cid]
    data.table(
        censusID    = cid,
        census_num  = meta$census_num,
        census_year = meta$census_year,
        S_richness  = S,
        H_shannon   = round(H, 4),
        J_evenness  = round(H / log(S), 4),
        D1mD        = round(1 - D_c, 4),
        D_recip     = round(1 / D_c, 2)
    )
})
census_div_dt <- rbindlist(census_div_list)
setorder(census_div_dt, census_num)

cat("\n", strrep("=", 64), "\n")
cat("  DIVERSITY INDICES ACROSS ALL CENSUSES\n")
cat(strrep("=", 64), "\n")
print(census_div_dt)

# SECTION 13 has been removed; species turnover is not included in this version.

# ============================================================
# SECTION 13: DBH SIZE CLASS DISTRIBUTION (Census 1 vs. Census 9)
# ============================================================

breaks_mm <- c(10, 20, 30, 50, 70, 100, 150, 200, 300, 500, Inf)
labels_mm <- c(
    "10–20", "20–30", "30–50", "50–70", "70–100",
    "100–150", "150–200", "200–300", "300–500", "≥500"
)

stems_c2c9 <- rbind(
    rec_alive[censusID == FIRST_ID, .(censusID, dbh)],
    rec_alive[censusID == LAST_ID, .(censusID, dbh)]
)
stems_c2c9 <- merge(stems_c2c9,
    census_meta[, .(censusID, census_year)],
    by = "censusID"
)
stems_c2c9[, dbh_class := cut(dbh,
    breaks = breaks_mm,
    labels = labels_mm, right = FALSE
)]

dbh_dist <- stems_c2c9[!is.na(dbh_class), .N, by = .(census_year, dbh_class)]
dbh_dist[, tot := sum(N), by = census_year]
dbh_dist[, pct := round(100 * N / tot, 1)]

cat("\n", strrep("=", 64), "\n")
cat("  DBH SIZE CLASS — count (%) by census\n")
cat(strrep("=", 64), "\n")
wide_n <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "N", fill = 0L)
wide_p <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "pct", fill = 0)
cat("  Counts:\n")
print(wide_n)
cat("  Proportions (%):\n")
print(wide_p)
