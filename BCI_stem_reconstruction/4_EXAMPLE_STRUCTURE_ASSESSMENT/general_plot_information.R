rm(list = ls())
workspace_root <- getwd()

# ============================================================
# GENERAL PLOT INFORMATION: BCI 50-ha Plot
# ============================================================
# Script  : general_plot_information.R
# Purpose : Descriptive summaries of forest structure and
#           composition for the Barro Colorado Island (BCI)
#           50-ha permanent plot across all available censuses.
#
# Plot layout
#   Dimensions   : 1000 m (E–W) × 500 m (N–S) = 50 ha
#   Base quadrats: 20 × 20 m (0.04 ha each); 1,250 total
#   Censuses     : 9 total (Census 1 = 1982 excluded from
#                  temporal analyses; Censuses 2–9 retained)
#
# Key definitions used throughout
#   Individual   : one unique treeID. A tree may have multiple
#                  stems, but it is counted once as one individual.
#   Alive stem   : Rstatus == "A" with valid (gx, gy) coordinates.
#   Basal area   : π/4 × (DBH mm / 1000)² m² per stem;
#                  summed across all alive stems per treeID.
#   Per-hectare  : bootstrapped mean across super-quadrats
#
# Script organisation
#   Sec. 1  – Libraries & constants
#   Sec. 2  – Helper functions
#   Sec. 3  – Data loading
#   Sec. 4  – Tree-level aggregation (alive stems only)
#   Sec. 5  – Overall plot summary       (most recent census)
#   Sec. 6  – Species-level summaries    (most recent census)
#   Sec. 7  – Genus & family summaries   (most recent census)
#   Sec. 8  – Abundance / BA thresholds  (most recent census)
#   Sec. 9  – Lifeform breakdown         (most recent census)
#   Sec. 10 – Per-hectare estimates, multiple quadrat sizes
#   Sec. 11 – Diversity indices          (most recent census)
#   Sec. 12 – Temporal trends, bootstrapped per-ha (all censuses)
#   Sec. 13 – Diversity indices across all censuses
#   Sec. 14 – DBH size-class distribution (Census 2 vs. Census 9)
#
# Dependencies: data.table, boot, HDInterval
#
# NOTE: These are basic descriptive summaries intended to give
# an overall view of forest structure and composition. Understanding
# temporal patterns and their ecological drivers requires more
# developed statistical analyses beyond what is shown here.
# ============================================================

# ============================================================
# SECTION 1: LIBRARIES & CONSTANTS
# ============================================================

library(data.table)
library(boot)
# HDInterval is called via :: inside boot_mean_ci().
# Install once if needed: install.packages("HDInterval")

# ---- Plot geometry --------------------------------------------------
PLOT_X_M <- 1000L # plot width  (m), E–W
PLOT_Y_M <- 500L # plot height (m), N–S
PLOT_AREA_HA <- PLOT_X_M * PLOT_Y_M / 1e4 # = 50 ha

# ---- Quadrat settings -----------------------------------------------
# BASE_Q_M is the quadrat size used for all temporal bootstrap analyses.
# QUAD_SIZES_M lists every valid square size that tiles 1000 × 500 m
# with no gaps (must divide both dimensions exactly).
# GCD(1000, 500) = 500 → valid sizes ≥ 20 m: 20, 25, 50, 100, 125, 250.
BASE_Q_M <- 20L
QUAD_SIZES_M <- c(20, 25, 50, 100, 125, 250)

# ---- Census metadata ------------------------------------------------
N_CENSUSES <- 9L
# Census 1 (1982) is excluded from temporal trend analyses because it
# predates the standardised protocol used in Censuses 2–9.

# ---- Bootstrap settings ---------------------------------------------
BOOT_R <- 4999L # number of bootstrap replicates
set.seed(42) # for reproducibility

# ============================================================
# SECTION 2: HELPER FUNCTIONS
# ============================================================
# All functions are defined here so they are available
# throughout the remainder of the script.

# ---- Modal value ----------------------------------------------------
# Returns the most frequent element of x. Used to assign a single
# representative calendar year to each census (some censuses span
# two calendar years).
get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

# ---- Bootstrap mean with 95 % highest density interval --------------
# Bootstraps the mean of x (R replicates) and returns a 95 % HDI
# from the bootstrap distribution using HDInterval::hdi().
#
# Fallback to ±1.96 SE is applied when:
#   • n < 4 (too few observations to bootstrap reliably), or
#   • the HDI cannot be computed (non-finite values).
#
# Returns a named numeric vector with elements: mean, ci_lo, ci_hi.
# unname() is applied to the HDI bounds to prevent R from creating
# compound names (e.g. "ci_lo.lower") when the result is embedded
# in a larger named vector downstream.
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
    c(mean = m, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

# ---- Build a complete super-quadrat grid and aggregate --------------
# Assigns each tree to a qs_m × qs_m cell using its (gx, gy) position,
# constructs the COMPLETE grid (n_x × n_y cells), zero-fills empty
# cells, and returns per-hectare densities for trees, basal area, and
# species richness per cell.
#
# Why a complete grid?
#   Empty quadrats must contribute zeros to the bootstrap distribution
#   so that the CI reflects true spatial heterogeneity. Omitting them
#   would upward-bias both the mean and the CI width.
#
# Requirement: qs_m must divide PLOT_X_M and PLOT_Y_M exactly.
#              This is guaranteed for all values in QUAD_SIZES_M.
make_superquad <- function(tree_dt, qs_m) {
    dt <- copy(tree_dt)
    n_x <- as.integer(PLOT_X_M / qs_m)
    n_y <- as.integer(PLOT_Y_M / qs_m)
    qs_ha <- (qs_m * qs_m) / 1e4

    # 0-based integer cell indices; pmin() guards against floating-point
    # coordinates that land exactly on the right/top plot boundary.
    dt[, sxi := pmin(as.integer(floor(gx / qs_m)), n_x - 1L)]
    dt[, syi := pmin(as.integer(floor(gy / qs_m)), n_y - 1L)]
    dt[, sqid := sxi * n_y + syi]

    sq_obs <- dt[, .(
        n_trees = .N,
        ba_m2   = sum(ba_m2, na.rm = TRUE),
        n_sp    = uniqueN(sp)
    ), by = sqid]

    # Left-join to the full cell index to zero-fill unoccupied quadrats
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
# SECTION 3: DATA LOADING
# ============================================================
# Loads all 9 stem-census tables, stacks them into a single
# data.table (rec), attaches species taxonomy and lifeform codes,
# and computes stem-level basal area and census year metadata.

message("Loading census data ...")

bci_nums <- as.character(seq_len(N_CENSUSES))
census_list <- lapply(bci_nums, function(n) {
    fp <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RTABLES", paste0("bci.stem", n, ".Rdata"))
    if (!file.exists(fp)) stop("Missing file: ", fp)
    e <- new.env(parent = emptyenv())
    load(fp, envir = e)
    get(paste0("bci.stem", n), envir = e)
})
names(census_list) <- paste0("bci.stem", bci_nums)

# Stack all censuses; censusID column identifies each source table
rec <- rbindlist(census_list, fill = TRUE, idcol = "censusID")
rec <- rec[!is.na(quadrat)] # remove stems with no quadrat assignment
rm(census_list, bci_nums)
gc()

# ---- Species / taxonomy table ---------------------------------------
load(file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RTABLES", "bci.spptable.rdata"))
bci.spptable <- as.data.table(bci.spptable)

# Append infraspecific epithet to the Latin name when available
# (e.g. "Quercus robur subsp. robur")
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

# Retain only the columns needed downstream
lf_cols <- intersect(
    c("Lifeform_RFoster", "Lifeform_RPerez_SAguilar"),
    names(bci.spptable)
)
keep_spp <- c("sp", "Family", "Genus", "Latin", lf_cols)
bci.spptable <- bci.spptable[, keep_spp, with = FALSE]

# Join taxonomy onto the stem records
rec <- merge(rec, bci.spptable, by = "sp", all.x = TRUE)

# ---- Taper correction -------------------------------------
# source(file.path(workspace_root, "BCI_stem_reconstruction", "1_DATA_PREPARATION", "HELPER_FUNCTIONS", "taper_correction.R"))

# rec <- apply_taper_correction(rec,
#     dbh_col = "dbh", hom_col = "hom", wsg_col = NULL, output_col = "dbh_t",
#     taper_correction = TRUE, common_hom = 1.3, convert_units = TRUE,
#     verbose = TRUE, overwrite = TRUE
# )

# ---- Taper correction -------------------------------------
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

# with(rec[CensusID == 9], plot(dbh_raw, dbh))
# abline(a = 0, b = 1, col = "red")

# ---- Stem-level basal area (m²) -------------------------------------
# BA = π/4 × (DBH in m)²; DBH is stored in mm, hence ÷ 1000.
rec[, ba_m2 := pi / 4 * (dbh / 1000)^2]

# ---- Census year lookup ---------------------------------------------
# Some censuses span two calendar years; the modal year is used as
# the single representative year for each census.
CENSUS_YEARS <- rec[, year := year(ExactDate)][
    , .(CensusID, year)
][
    , .(year = get_mode(year)),
    by = CensusID
][
    order(CensusID), year
]

census_meta <- data.table(
    censusID    = paste0("bci.stem", seq_len(N_CENSUSES)),
    census_num  = seq_len(N_CENSUSES),
    census_year = CENSUS_YEARS
)

# ============================================================
# SECTION 4: TREE-LEVEL AGGREGATION (alive stems only)
# ============================================================
# Collapses the stem-level table to one row per individual
# (unique treeID) per census.
#
# Filtering rules applied before aggregation:
#   • Rstatus == "A"        : resolved alive status only
#   • !is.na(gx) & !is.na(gy) : valid spatial coordinates required
#
# Tree-level columns produced:
#   ba_m2   – sum of basal areas of all alive stems for that treeID
#   n_stems – count of alive stems (informational; not used as the
#             individual count, which is always the treeID count)
#   gx, gy  – mean coordinate of alive stems (for quadrat assignment)
#
# IMPORTANT: all n_individuals / n_trees counts in this script refer
# to unique treeIDs, not stems. last_trees has exactly one row per
# treeID in the most recent census, so .N on that table always gives
# the individual (treeID) count.

rec_alive <- rec[Rstatus == "A" & !is.na(gx) & !is.na(gy)]

by_cols <- c(
    "censusID", "CensusID", "treeID",
    "sp", "Family", "Genus", "Latin",
    lf_cols
)

rec_tree <- rec_alive[, .(
    ba_m2   = sum(ba_m2, na.rm = TRUE), # total BA of all alive stems
    n_stems = .N, # alive-stem count (≥1 per tree)
    gx      = mean(gx, na.rm = TRUE), # mean E–W position
    gy      = mean(gy, na.rm = TRUE) # mean N–S position
), by = by_cols]

# Clip coordinates to the strict plot interior (floating-point guard)
rec_tree[gx < 0, gx := 0]
rec_tree[gx >= PLOT_X_M, gx := PLOT_X_M - 1e-3]
rec_tree[gy < 0, gy := 0]
rec_tree[gy >= PLOT_Y_M, gy := PLOT_Y_M - 1e-3]

rec_tree <- merge(rec_tree, census_meta, by = "censusID", all.x = TRUE)
setorder(rec_tree, census_num, treeID)

# Validate that every tree received a census_num after the merge
if (any(is.na(rec_tree$census_num))) {
    warning("Some treeIDs have missing census_num after merge; dropping affected rows.")
    rec_tree <- rec_tree[!is.na(census_num)]
}

# ---- Census-level convenience subsets -------------------------------
all_cids <- census_meta$censusID

# Census 1 (1982) used a different sampling protocol and is excluded
# from temporal comparisons. Census 2 serves as the historical baseline.
FIRST_ID <- all_cids[2] # earliest comparable census
LAST_ID <- all_cids[N_CENSUSES] # most recent census

# One row per treeID in the most recent census
last_trees <- rec_tree[censusID == LAST_ID]
# One row per stemID in the most recent census (for stem-level analyses)
last_stems <- rec_alive[censusID == LAST_ID]

# ---- Plot-level scalars for the most recent census ------------------
total_trees <- nrow(last_trees) # unique individuals
total_stems <- nrow(last_stems) # total alive stems
total_ba <- sum(last_trees$ba_m2, na.rm = TRUE) # total BA (m²)
n_species <- uniqueN(last_trees$sp)
n_families <- uniqueN(last_trees$Family)
n_genera <- uniqueN(last_trees$Genus)

# ============================================================
# SECTION 5: OVERALL PLOT SUMMARY (most recent census)
# ============================================================

cat("\n", strrep("=", 64), "\n")
cat(sprintf(
    "  BCI 50-ha PLOT — Census %d (alive trees, DBH ≥ 10 mm)\n",
    N_CENSUSES
))
cat(strrep("=", 64), "\n")
cat(sprintf(
    "  Plot dimensions:                     %d × %d m = %g ha\n",
    PLOT_X_M, PLOT_Y_M, PLOT_AREA_HA
))
cat(sprintf("  Total individuals (treeIDs, alive):  %d\n", total_trees))
cat(sprintf("  Total stems       (stemIDs, alive):  %d\n", total_stems))
cat(sprintf("  Species richness:                    %d\n", n_species))
cat(sprintf("  Genera:                              %d\n", n_genera))
cat(sprintf("  Families:                            %d\n", n_families))
cat(sprintf("  Total basal area (m²):               %.2f\n", total_ba))
cat(sprintf("  Crude individuals ha⁻¹:              %.1f\n", total_trees / PLOT_AREA_HA))
cat(sprintf("  Crude BA ha⁻¹ (m²/ha):              %.2f\n", total_ba / PLOT_AREA_HA))

# ============================================================
# SECTION 6: SPECIES-LEVEL SUMMARIES (most recent census)
# ============================================================
# last_trees has one row per treeID, so .N = unique individuals.
# Relative abundance and relative BA are expressed as % of plot totals.
# Importance Value (IV) = rel. abundance % + rel. BA % (max = 200).

spp_sum <- last_trees[, .(
    n_individuals = .N, # unique treeIDs
    total_ba_m2   = sum(ba_m2, na.rm = TRUE),
    mean_ba_m2    = mean(ba_m2, na.rm = TRUE)
), by = .(sp, Latin)]
spp_sum[, rel_abund_pct := round(100 * n_individuals / total_trees, 3)]
spp_sum[, rel_ba_pct := round(100 * total_ba_m2 / total_ba, 3)]
spp_sum[, IV := rel_abund_pct + rel_ba_pct]

cat("\n", strrep("-", 64), "\n")
cat(sprintf("  TOP 10 SPECIES BY INDIVIDUALS — Census %d\n", N_CENSUSES))
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-n_individuals),
    .(Latin, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct)
], 10))

cat("\n", strrep("-", 64), "\n")
cat(sprintf("  TOP 10 SPECIES BY BASAL AREA (m²) — Census %d\n", N_CENSUSES))
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-total_ba_m2),
    .(Latin, total_ba_m2, rel_ba_pct, n_individuals, rel_abund_pct)
], 10))

cat("\n", strrep("-", 64), "\n")
cat(sprintf(
    "  TOP 10 SPECIES BY IMPORTANCE VALUE — Census %d\n",
    N_CENSUSES
))
cat("  (IV = rel. abundance %% + rel. BA %%; max = 200)\n")
cat(strrep("-", 64), "\n")
print(head(spp_sum[
    order(-IV),
    .(Latin, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
], 10))

# ============================================================
# SECTION 7: GENUS & FAMILY SUMMARIES (most recent census)
# ============================================================
# Individuals are counted as unique treeIDs. last_trees already has
# exactly one row per treeID, so .N gives the individual count
# directly — no additional deduplication is needed.
#
# Columns produced for both genus and family tables:
#   n_individuals : unique treeIDs (individuals, NOT stems)
#   n_species     : number of species in the taxon group
#   total_ba_m2   : summed basal area (m²) of all individuals
#   rel_abund_pct : % of plot-total individual count
#   rel_ba_pct    : % of plot-total basal area
#   IV            : Importance Value = rel_abund_pct + rel_ba_pct
#                   (max = 200; combines dominance in abundance & BA)

# ---- Genus-level summaries ------------------------------------------
genus_sum <- last_trees[!is.na(Genus), .(
    n_individuals = .N,
    n_species     = uniqueN(sp),
    total_ba_m2   = sum(ba_m2, na.rm = TRUE)
), by = Genus]
genus_sum[, rel_abund_pct := round(100 * n_individuals / total_trees, 2)]
genus_sum[, rel_ba_pct := round(100 * total_ba_m2 / total_ba, 2)]
genus_sum[, IV := rel_abund_pct + rel_ba_pct]

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  GENUS SUMMARIES — Census %d\n", N_CENSUSES))
cat(strrep("=", 64), "\n")

cat("\n  >> Top 15 genera by individuals (unique treeIDs):\n")
print(head(genus_sum[
    order(-n_individuals),
    .(Genus, n_species, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct)
], 15))

cat("\n  >> Top 15 genera by basal area (m²):\n")
print(head(genus_sum[
    order(-total_ba_m2),
    .(Genus, n_species, total_ba_m2, rel_ba_pct, n_individuals, rel_abund_pct)
], 15))

cat("\n  >> Top 15 genera by Importance Value (max = 200):\n")
print(head(genus_sum[
    order(-IV),
    .(Genus, n_species, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
], 15))

# ---- Family-level summaries -----------------------------------------
family_sum <- last_trees[!is.na(Family), .(
    n_individuals = .N,
    n_genera      = uniqueN(Genus),
    n_species     = uniqueN(sp),
    total_ba_m2   = sum(ba_m2, na.rm = TRUE)
), by = Family]
family_sum[, rel_abund_pct := round(100 * n_individuals / total_trees, 2)]
family_sum[, rel_ba_pct := round(100 * total_ba_m2 / total_ba, 2)]
family_sum[, IV := rel_abund_pct + rel_ba_pct]

cat("\n", strrep("=", 64), "\n")
cat(sprintf("  FAMILY SUMMARIES — Census %d\n", N_CENSUSES))
cat(strrep("=", 64), "\n")

cat("\n  >> Top 15 families by individuals (unique treeIDs):\n")
print(head(family_sum[
    order(-n_individuals),
    .(Family, n_genera, n_species, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct)
], 15))

cat("\n  >> Top 15 families by basal area (m²):\n")
print(head(family_sum[
    order(-total_ba_m2),
    .(Family, n_genera, n_species, total_ba_m2, rel_ba_pct, n_individuals, rel_abund_pct)
], 15))

cat("\n  >> Top 15 families by Importance Value (max = 200):\n")
print(head(family_sum[
    order(-IV),
    .(Family, n_genera, n_species, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
], 15))

# ============================================================
# SECTION 8: ABUNDANCE / BA THRESHOLDS (most recent census)
# ============================================================
# Identifies species simultaneously dominant in both abundance and
# basal area relative to the plot mean and the 75th percentile.
# These dual-dominance species tend to have the greatest structural
# and functional importance in the stand.

avg_n <- mean(spp_sum$n_individuals)
avg_ba <- mean(spp_sum$total_ba_m2)
q75_n <- quantile(spp_sum$n_individuals, 0.75)
q75_ba <- quantile(spp_sum$total_ba_m2, 0.75)

rare_species_half_thr <- avg_n / 2
n_species_below_half <- sum(spp_sum$n_individuals < rare_species_half_thr)
pct_species_below_half <- round(100 * n_species_below_half / n_species, 1)
n_individuals_below_half <- sum(spp_sum[n_individuals < rare_species_half_thr, n_individuals])
pct_individuals_below_half <- round(100 * n_individuals_below_half / total_trees, 1)

n_species_rare10 <- sum(spp_sum$n_individuals <= 10)
pct_species_rare10 <- round(100 * n_species_rare10 / n_species, 1)
n_individuals_rare10 <- sum(spp_sum[n_individuals <= 10, n_individuals])
n_singletons <- sum(spp_sum$n_individuals == 1)

cat("\n  RARITY SUMMARY\n")
cat(strrep("-", 64), "\n")
cat(sprintf("  Threshold: < %.1f individuals (half the species mean)\n", rare_species_half_thr))
cat(sprintf(
    "  Species below threshold: %d (%.1f%% of species), representing %d individuals (%.1f%% of total trees)\n",
    n_species_below_half, pct_species_below_half,
    n_individuals_below_half, pct_individuals_below_half
))
cat(sprintf(
    "  Species with ≤ 10 individuals: %d (%.1f%% of species), representing %d individuals\n",
    n_species_rare10, pct_species_rare10, n_individuals_rare10
))
cat(sprintf("  Singletons: %d species with exactly 1 individual\n", n_singletons))

cat("\n", strrep("-", 64), "\n")
cat("  ABUNDANCE / BA THRESHOLDS\n")
cat(strrep("-", 64), "\n")
cat(sprintf("  Mean individuals/species:   %6.1f  | 75th pct: %6.0f\n", avg_n, q75_n))
cat(sprintf("  Mean BA/species (m²):     %8.4f  | 75th pct: %8.4f\n", avg_ba, q75_ba))

cat("\n  >> Species above MEAN in both individuals & BA:\n")
print(spp_sum[n_individuals > avg_n & total_ba_m2 > avg_ba][
    order(-IV),
    .(Latin, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
])
cat("\n  >> Species above 75th percentile in both individuals & BA:\n")
print(spp_sum[n_individuals > q75_n & total_ba_m2 > q75_ba][
    order(-IV),
    .(Latin, n_individuals, rel_abund_pct, total_ba_m2, rel_ba_pct, IV)
])

# ============================================================
# SECTION 9: LIFEFORM BREAKDOWN (most recent census)
# ============================================================
# Summarises individuals, species richness, and basal area by lifeform
# category. Loops over all available lifeform classification columns
# found in the species table (lf_cols). Within each loop iteration,
# the top 5 species per group are listed by individual count.
# n_individuals counts unique treeIDs (not stems).

for (lfc in lf_cols) {
    lf <- last_trees[!is.na(get(lfc)), .(
        n_individuals = .N,
        n_species     = uniqueN(sp),
        ba_m2         = sum(ba_m2, na.rm = TRUE)
    ), by = c(lfc)][order(-n_individuals)]
    lf[, rel_n_pct := round(100 * n_individuals / total_trees, 1)]
    lf[, rel_ba_pct := round(100 * ba_m2 / total_ba, 1)]

    cat("\n", strrep("-", 64), "\n")
    cat(sprintf("  LIFEFORM BREAKDOWN — %s\n", lfc))
    cat(strrep("-", 64), "\n")
    print(lf)

    # Top 5 species per lifeform group by individual count (treeIDs)
    cat(sprintf("\n  >> Top 5 species per group (%s):\n", lfc))
    for (g in sort(unique(na.omit(last_trees[[lfc]])))) {
        s <- last_trees[get(lfc) == g, .N, by = .(sp, Latin)][order(-N)]
        cat(sprintf("  [%s]\n", g))
        print(head(s[, .(Latin, N)], 5))
    }
}

# ============================================================
# SECTION 10: PER-HECTARE ESTIMATES — MULTIPLE QUADRAT SIZES
#             (most recent census)
# ============================================================
# Spatial bootstrap methodology:
#   1. Assign each tree to a qs × qs super-quadrat via (gx, gy).
#   2. Build the COMPLETE grid (n_x × n_y cells); empty cells = 0.
#   3. Scale each cell total to trees ha⁻¹ and BA ha⁻¹.
#   4. Bootstrap the mean across cells → point estimate + 95 % HDI CI.
#
# The unit of replication is the super-quadrat, not the whole plot,
# so spatial heterogeneity is properly propagated into the CI.
# Omitting empty cells would upward-bias both the mean and the CI.

message("Computing per-ha bootstrap estimates by quadrat size ...")

perha_rows <- lapply(QUAD_SIZES_M, function(qs) {
    qs_ha <- (qs * qs) / 1e4
    message(sprintf("  %d × %d m (%.4f ha) ...", qs, qs, qs_ha))
    sq <- make_superquad(last_trees, qs)
    rt <- boot_mean_ci(sq$trees_ha)
    rb <- boot_mean_ci(sq$ba_ha)
    rs <- boot_mean_ci(sq$n_sp)
    data.table(
        quad_m      = qs,
        quad_ha     = qs_ha,
        n_quads     = nrow(sq),
        trees_ha    = rt["mean"], trees_ci_lo = rt["ci_lo"], trees_ci_hi = rt["ci_hi"],
        ba_ha       = rb["mean"], ba_ci_lo    = rb["ci_lo"], ba_ci_hi    = rb["ci_hi"],
        sp_per_quad = rs["mean"], sp_ci_lo    = rs["ci_lo"], sp_ci_hi    = rs["ci_hi"]
    )
})
perha_dt <- rbindlist(perha_rows)

cat("\n", strrep("=", 72), "\n")
cat(sprintf(
    "  PER-HECTARE ESTIMATES BY QUADRAT SIZE (bootstrap 95%% CI, R = %d)\n",
    BOOT_R
))
cat(strrep("=", 72), "\n")
print(perha_dt[, .(
    quad_m, quad_ha, n_quads,
    trees_ha = round(trees_ha, 1),
    trees_ci_lo = round(trees_ci_lo, 1),
    trees_ci_hi = round(trees_ci_hi, 1),
    ba_ha = round(ba_ha, 2),
    ba_ci_lo = round(ba_ci_lo, 2),
    ba_ci_hi = round(ba_ci_hi, 2)
)])

# ============================================================
# SECTION 11: DIVERSITY INDICES (most recent census)
# ============================================================
# Computed from the species abundance distribution in spp_sum.
# Abundance measure: unique treeIDs (individuals) per species,
# consistent with the individual-based convention used throughout.
#
# Indices reported:
#   S   – species richness (species with ≥1 individual)
#   H'  – Shannon-Wiener entropy (nats; natural log base)
#   J   – Pielou's evenness = H' / ln(S)  [0, 1]
#   D   – Simpson's concentration = Σ pᵢ²
#   1-D – Simpson's diversity
#   1/D – Simpson's reciprocal diversity

pi_vals <- spp_sum$n_individuals / sum(spp_sum$n_individuals)
H_shannon <- -sum(pi_vals[pi_vals > 0] * log(pi_vals[pi_vals > 0]))
D_conc <- sum(pi_vals^2)
n_species <- sum(spp_sum$n_individuals > 0) # species with ≥1 individual
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
# SECTION 12: TEMPORAL TRENDS — BOOTSTRAPPED PER-HA
#             (Censuses 2–9)
# ============================================================
# Applies the same spatial-bootstrap method as Section 10
# (20 m base quadrats, complete grid) to every census from 2 to 9.
# Census 1 (1982) is excluded due to protocol differences.
#
# n_trees in the output refers to unique treeIDs per census
# (i.e. individuals, not stems).

message("Computing bootstrapped per-ha temporal trends ...")

census_boot_list <- lapply(all_cids[-1], function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0L) {
        return(NULL)
    }
    sq <- make_superquad(ct, BASE_Q_M) # 20 m quadrats, complete grid
    rt <- boot_mean_ci(sq$trees_ha)
    rb <- boot_mean_ci(sq$ba_ha)
    meta <- census_meta[censusID == cid]
    data.table(
        censusID    = cid,
        census_num  = meta$census_num,
        census_year = meta$census_year,
        n_trees     = nrow(ct), # unique treeIDs
        n_species   = uniqueN(ct$sp),
        total_ba_m2 = round(sum(ct$ba_m2, na.rm = TRUE), 1),
        trees_ha    = rt["mean"], trees_ci_lo = rt["ci_lo"], trees_ci_hi = rt["ci_hi"],
        ba_ha       = rb["mean"], ba_ci_lo    = rb["ci_lo"], ba_ci_hi    = rb["ci_hi"]
    )
})
census_boot_dt <- rbindlist(census_boot_list)
setorder(census_boot_dt, census_num)

cat("\n", strrep("=", 72), "\n")
cat("  TEMPORAL TRENDS — BOOTSTRAPPED PER-HA (20 m quadrats, Censuses 2–9)\n")
cat(strrep("=", 72), "\n")
print(census_boot_dt[, .(
    censusID, census_year, n_trees, n_species,
    trees_ha = round(trees_ha, 1),
    trees_ci_lo = round(trees_ci_lo, 1),
    trees_ci_hi = round(trees_ci_hi, 1),
    ba_ha = round(ba_ha, 2),
    ba_ci_lo = round(ba_ci_lo, 2),
    ba_ci_hi = round(ba_ci_hi, 2)
)])

# ============================================================
# SECTION 13: DIVERSITY INDICES ACROSS ALL CENSUSES (2–9)
# ============================================================
# Recomputes Shannon-Wiener, Pielou's evenness, and Simpson's indices
# for each census using unique treeIDs as the abundance measure,
# consistent with the individual-based convention used throughout.

census_div_list <- lapply(all_cids[-1], function(cid) {
    ct <- rec_tree[censusID == cid]
    if (nrow(ct) == 0L) {
        return(NULL)
    }
    # Species proportional abundances based on individual (treeID) counts
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
cat("  DIVERSITY INDICES ACROSS ALL CENSUSES (2–9)\n")
cat(strrep("=", 64), "\n")
print(census_div_dt)

# ============================================================
# SECTION 14: DBH SIZE-CLASS DISTRIBUTION
#             (Census 2 vs. Census 9 — stem level)
# ============================================================
# DBH is a stem-level measurement, so this section reports STEM
# counts rather than individual (treeID) counts. This is the only
# section of the script in which stems, not individuals, are the
# counting unit (noted explicitly below).
#
# Comparison: Census 2 (earliest comparable) vs. Census 9 (most
# recent), bracketing the full standardised observation period.
#
# Size classes (mm), left-closed right-open intervals [lower, upper):
#   10–20, 20–30, 30–50, 50–70, 70–100,
#   100–150, 150–200, 200–300, 300–500, ≥500

breaks_mm <- c(10, 20, 30, 50, 70, 100, 150, 200, 300, 500, Inf)
labels_mm <- c(
    "10–20", "20–30", "30–50", "50–70", "70–100",
    "100–150", "150–200", "200–300", "300–500", "≥500"
)

stems_c2c9 <- rbind(
    rec_alive[censusID == FIRST_ID, .(censusID, dbh)],
    rec_alive[censusID == LAST_ID, .(censusID, dbh)]
)
stems_c2c9 <- merge(
    stems_c2c9,
    census_meta[, .(censusID, census_year)],
    by = "censusID"
)
stems_c2c9[, dbh_class := cut(
    dbh,
    breaks = breaks_mm,
    labels = labels_mm,
    right  = FALSE # left-closed [lower, upper)
)]

dbh_dist <- stems_c2c9[!is.na(dbh_class), .N, by = .(census_year, dbh_class)]
dbh_dist[, tot := sum(N), by = census_year]
dbh_dist[, pct := round(100 * N / tot, 1)]

cat("\n", strrep("=", 64), "\n")
cat("  DBH SIZE-CLASS DISTRIBUTION — stem counts by census\n")
cat("  (stem-level counts; left-closed intervals in mm)\n")
cat(strrep("=", 64), "\n")
wide_n <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "N", fill = 0L)
wide_p <- dcast(dbh_dist, dbh_class ~ census_year, value.var = "pct", fill = 0)
cat("  Counts:\n")
print(wide_n)
cat("  Proportions (%):\n")
print(wide_p)
