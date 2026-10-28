# =============================================================================
# 3_census-9-summary.R — the census-9 (2022–2023) summary of the Abstract and
# the Data Description of the Data Paper, from the files of this package: living trees,
# species, genera and families, density and basal area per hectare with 95%
# intervals, rarity, and the leading species and families.
#
#   - Living trees: trees with a living stem (Rstatus A); a tree's basal area
#     is the sum over its living stems.
#   - With use_taper_stocks = TRUE (below), diameters are corrected to 1.3 m
#     with Model 1 of Cushman et al. (2014) and the species terms of their
#     Table 3, palms (Arecaceae) keeping the recorded diameter; with FALSE,
#     every diameter is used as recorded, at its height of measurement.
#     Basal area = pi/4 * (DBH in m)^2.
#   - Per-hectare values: mean over the 1,250 quadrats of 20 x 20 m (trees
#     placed by the mean coordinates of their living stems, or at their
#     quadrat centre without coordinates), with the 95% highest-density
#     interval of 4,999 bootstrap means (seed 42).
#
# Run from the folder that holds data/ and code/ (the repository folder), in
# R or with:  Rscript code/3_census-9-summary.R
# Writes census-9_summary.csv, census-9_species.csv and census-9_families.csv
# to ../manuscript_outputs when that folder exists, else to outputs/.
# Needs data.table, boot and HDInterval.
# =============================================================================
rm(list = ls())

suppressMessages(library(data.table))
source(file.path("code", "load-bci-tables.R"))
DATA <- "data"
OUT <- if (dir.exists("../manuscript_outputs")) "../manuscript_outputs" else "outputs"
dir.create(OUT, showWarnings = FALSE)

# Basal area from the taper-corrected DBH (TRUE) or from the DBH as recorded,
# at its height of measurement (FALSE)
use_taper_stocks <- FALSE

bci <- load_bci_tables(DATA, censuses = 9)
spp <- as.data.table(bci$bci.spptable)
c9 <- as.data.table(bci$bci.stem9)
rm(bci)

# ---- Taper correction (Cushman et al. 2014, Model 1, Table 3) ------------------------
TAPER_SPECIES_2014 <- c(
    quaras = 0.4865, ceibpe = 0.0495, huracr = 0.0086, diptpa = 0.0150, alsebl = 0.1440,
    brosal = 0.3874, cavapl = 0.7982, anacex = 0.4371, tet2pa = -0.5141, tab2ar = 0.3576
)
taper_dbh <- function(dbh_mm, hom, sp) {
    hom[is.na(hom)] <- 1.3
    dbh_cm <- dbh_mm / 10
    dbh_cm[dbh_cm <= 0] <- NA_real_
    hom_log <- hom
    hom_log[hom_log <= 0] <- NA_real_
    term <- unname(TAPER_SPECIES_2014[sp])
    term[is.na(term)] <- 0
    b <- exp(-2.0205 + term - 0.5053 * log(dbh_cm) + 0.3748 * log(hom_log))
    out <- 10 * dbh_cm * exp(b * (hom - 1.3))
    out[!is.finite(out)] <- NA_real_
    out
}

# ---- Living stems and trees ----------------------------------------------------------------
c9 <- c9[Rstatus == "A" & !is.na(dbh) & !is.na(quadrat)]
c9 <- merge(c9, spp[, .(sp, Family, Genus, Latin)], by = "sp", all.x = TRUE)
if (use_taper_stocks) {
    c9[, dbh_ba := taper_dbh(dbh, hom, sp)]
    c9[Family %in% "Arecaceae" | is.na(dbh_ba), dbh_ba := dbh]
} else {
    c9[, dbh_ba := dbh]
}
c9[, ba_m2 := pi / 4 * (dbh_ba / 1000)^2]
c9[is.na(gx) | is.na(gy), `:=`(
    gx = as.integer(substr(quadrat, 1, 2)) * 20 + 10,
    gy = as.integer(substr(quadrat, 3, 4)) * 20 + 10
)]
trees <- c9[, .(ba_m2 = sum(ba_m2), gx = mean(gx), gy = mean(gy)), by = .(treeID, sp, Family, Genus, Latin)]
trees[gx < 0, gx := 0][gx >= 1000, gx := 1000 - 1e-3][gy < 0, gy := 0][gy >= 500, gy := 500 - 1e-3]
n_tot <- nrow(trees)
ba_tot <- sum(trees$ba_m2)

# ---- Per hectare, with bootstrap intervals over the 20 x 20 m quadrats ----------------------
q <- trees[, .(n = .N, ba = sum(ba_m2)), by = .(sqid = pmin(floor(gx / 20), 49) * 25 + pmin(floor(gy / 20), 24))]
q <- merge(data.table(sqid = 0:1249), q, by = "sqid", all.x = TRUE)
q[is.na(n), n := 0L][is.na(ba), ba := 0]
q[, `:=`(trees_ha = n / 0.04, ba_ha = ba / 0.04)]
set.seed(42)
boot_hdi <- function(v) {
    b <- boot::boot(v, function(d, i) mean(d[i]), R = 4999)
    HDInterval::hdi(as.numeric(b$t), credMass = 0.95)
}
ci_t <- boot_hdi(q$trees_ha)
ci_b <- boot_hdi(q$ba_ha)

# ---- Species  and families -------------------------------------------------------------------
sp9 <- trees[, .(individuals = .N, ba_m2 = sum(ba_m2)), by = .(sp, Latin, Family, Genus)]
sp9[, `:=`(pct_individuals = round(100 * individuals / n_tot, 2), pct_ba = round(100 * ba_m2 / ba_tot, 2))]
setorder(sp9, -individuals)
fam9 <- sp9[, .(genera = uniqueN(Genus), species = .N, individuals = sum(individuals), ba_m2 = sum(ba_m2)), by = Family]
fam9[, `:=`(pct_individuals = round(100 * individuals / n_tot, 2), pct_ba = round(100 * ba_m2 / ba_tot, 2))]
setorder(fam9, -individuals)
half <- mean(sp9$individuals) / 2

summary9 <- data.table(
    item = c(
        "living trees", "species", "genera", "families",
        "basal area, total (m2)", "density (trees/ha)", "density, 95% interval low", "density, 95% interval high",
        "basal area (m2/ha)", "basal area, 95% interval low", "basal area, 95% interval high",
        "half the mean abundance (individuals)", "species below half the mean abundance",
        "species below half the mean abundance (%)", "individuals of those species",
        "species with ten or fewer individuals", "individuals of those species", "singleton species"
    ),
    value = c(
        n_tot, uniqueN(trees$sp), uniqueN(trees$Genus), uniqueN(trees$Family),
        round(ba_tot, 2), round(mean(q$trees_ha)), round(ci_t[1]), round(ci_t[2]),
        round(mean(q$ba_ha), 2), round(ci_b[1], 2), round(ci_b[2], 2),
        round(half, 1), sp9[individuals < half, .N],
        round(100 * sp9[individuals < half, .N] / nrow(sp9), 1), sp9[individuals < half, sum(individuals)],
        sp9[individuals <= 10, .N], sp9[individuals <= 10, sum(individuals)], sp9[individuals == 1, .N]
    )
)
fwrite(summary9, file.path(OUT, "census-9_summary.csv"))
fwrite(sp9, file.path(OUT, "census-9_species.csv"))
fwrite(fam9, file.path(OUT, "census-9_families.csv"))

cat("Basal area from the", if (use_taper_stocks) "taper-corrected DBH" else "DBH as recorded", "\n")
print(summary9, nrows = 50)
cat("\nLeading species by individuals:\n")
print(sp9[1:3, .(Latin, individuals, pct_individuals)])
cat("Leading species by basal area:\n")
print(sp9[order(-ba_m2)][1:3, .(Latin, ba_m2 = round(ba_m2, 2), pct_ba)])
cat("Leading families by individuals:\n")
print(fam9[1:3, .(Family, genera, species, individuals, pct_individuals)])
cat("Leading families by basal area:\n")
print(fam9[order(-ba_m2)][1:3, .(Family, ba_m2 = round(ba_m2, 2), pct_ba)])
cat("Written to", normalizePath(OUT), "\n")
