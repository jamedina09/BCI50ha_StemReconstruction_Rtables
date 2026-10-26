# =============================================================================
# make-manuscript-tables.R — computes, from the files of this data package,
# every table and number of the Data Paper: Table 1 (censuses), the status
# counts of Table 4, the counts of the Data Description and the Limitations,
# and the census-9 summary of the Abstract.
#
# Run from the folder that holds data/ and code/:
#   Rscript code/make-manuscript-tables.R
# Needs data.table; boot and HDInterval for the census-9 intervals.
#
# Census-9 summary (Abstract):
#   - Living trees and stems: Rstatus A.
#   - Diameters are corrected to 1.3 m with Model 1 of Cushman et al. (2014)
#     and the species terms of their Table 3; palms (Arecaceae) keep the
#     recorded diameter. Basal area = pi/4 * (DBH in m)^2.
#   - Per-hectare density and basal area: mean over the 1,250 quadrats of
#     20 x 20 m (trees placed by the mean coordinates of their living stems, or
#     at their quadrat centre without coordinates), with a 95% highest-density
#     interval of 4,999 bootstrap means (seed 42).
# =============================================================================

suppressMessages(library(data.table))
source(file.path("code", "load-bci-tables.R"))
DATA <- "data"

say <- function(...) cat(sprintf(...), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

bci <- load_bci_tables(DATA)
spp <- as.data.table(bci$bci.spptable)
x <- rbindlist(lapply(bci[paste0("bci.stem", 1:9)], as.data.table))
x[, CensusID := as.integer(CensusID)]
rm(bci)

# ---- Table 1: censuses ----------------------------------------------------------
tab1 <- x[, .(
    first_date = min(ExactDate, na.rm = TRUE),
    last_date = max(ExactDate, na.rm = TRUE),
    records = sum(!is.na(DFstatus)),
    measurements = sum(!is.na(dbh)),
    living_stems = sum(Rstatus == "A"),
    living_trees = uniqueN(treeID[Rstatus == "A"])
), by = CensusID][order(CensusID)]
say("\nTable 1. Censuses")
print(tab1)

# ---- Table 4: status counts over the nine tables ----------------------------------
say("\nTable 4. Rstatus counts over the nine tables")
print(x[, .N, by = Rstatus][order(-N)])

# ---- Data Description ------------------------------------------------------------------
taxa <- spp[sp != "uniden"]
say("\nData Description")
say("Trees: %s; stems: %s; rows per census table: %s", fmt(uniqueN(x$treeID)), fmt(x[CensusID == 1, .N]), fmt(x[CensusID == 1, .N]))
say("Field records: %s; with a diameter: %s", fmt(sum(!is.na(x$DFstatus))), fmt(sum(!is.na(x$dbh))))
say("Trees with more than one stem: %s", fmt(x[CensusID == 1, .N, by = treeID][N > 1, .N]))
say("Taxa (species codes without uniden): %d in %d genera and %d families", nrow(taxa), uniqueN(taxa$Genus), uniqueN(taxa$Family))
say("Trees never identified (uniden): %s", fmt(x[sp == "uniden", uniqueN(treeID)]))
say("Trees with a reconstruction key (StemPaths): %s; without: %s",
    fmt(x[!is.na(StemPaths), uniqueN(treeID)]), fmt(uniqueN(x$treeID) - x[!is.na(StemPaths), uniqueN(treeID)]))
post <- load_bci_posterior(DATA)
say("Posterior samples: %s rows for %s trees", fmt(nrow(post)), fmt(uniqueN(post$treeID)))
say("Species codes whose name changed since Condit et al. (2019) (NameChange): %d", spp[!is.na(NameChange), .N])
say("Field-status records: %s", paste(x[!is.na(DFstatus), .N, by = DFstatus][order(-N), sprintf("%s %s", DFstatus, fmt(N))], collapse = "; "))

# Interval between consecutive censuses of a tree (mean date of its living stems)
td <- x[Rstatus == "A", .(d = mean(date)), by = .(treeID, CensusID)]
setorder(td, treeID, CensusID)
td[, `:=`(c0 = shift(CensusID), d0 = shift(d)), by = treeID]
iv <- td[!is.na(c0) & c0 == CensusID - 1L, .(median_years = round(median((d - d0) / 365.25), 2)), by = .(from = c0, to = CensusID)]
say("\nMedian interval between consecutive censuses of a tree (years)")
print(iv)
say("Census duration (years): %s", paste(sprintf("%d: %.2f", tab1$CensusID, as.numeric(tab1$last_date - tab1$first_date) / 365.25), collapse = "; "))

# ---- Limitations -------------------------------------------------------------------------
alive_no_dbh <- x[Rstatus == "A" & is.na(dbh)]
say("\nLimitations")
say("Alive cells without a diameter: %s (stumps, broken below: %s; missed measurements: %s)",
    fmt(nrow(alive_no_dbh)), fmt(alive_no_dbh[DFstatus %in% "broken below", .N]), fmt(alive_no_dbh[!DFstatus %in% "broken below", .N]))
say("Measurements taken away from 1.3 m: %.1f%%", 100 * x[!is.na(dbh) & !is.na(hom) & abs(hom - 1.3) > 0.01, .N] / x[!is.na(dbh), .N])
say("Trees without coordinates: %s", fmt(x[is.na(gx) | is.na(gy), uniqueN(treeID)]))

# ---- Census 9 (Abstract) ---------------------------------------------------------------------
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

c9 <- x[CensusID == 9 & Rstatus == "A" & !is.na(dbh) & !is.na(quadrat)]
c9 <- merge(c9, spp[, .(sp, Family, Genus, Latin)], by = "sp", all.x = TRUE)
c9[, dbh_t := taper_dbh(dbh, hom, sp)]
c9[Family %in% "Arecaceae" | is.na(dbh_t), dbh_t := dbh]
c9[, ba_m2 := pi / 4 * (dbh_t / 1000)^2]
c9[is.na(gx) | is.na(gy), `:=`(gx = as.integer(substr(quadrat, 1, 2)) * 20 + 10, gy = as.integer(substr(quadrat, 3, 4)) * 20 + 10)]
trees9 <- c9[, .(ba_m2 = sum(ba_m2), gx = mean(gx), gy = mean(gy)), by = .(treeID, sp, Family, Genus, Latin)]
trees9[gx < 0, gx := 0][gx >= 1000, gx := 1000 - 1e-3][gy < 0, gy := 0][gy >= 500, gy := 500 - 1e-3]

say("\nCensus 9 (2022-2023)")
say("Living trees: %s; living stems: %s; species: %d; genera: %d; families: %d",
    fmt(nrow(trees9)), fmt(nrow(c9)), uniqueN(trees9$sp), uniqueN(trees9$Genus), uniqueN(trees9$Family))
say("Basal area: %.2f m2 in total, %.2f m2/ha", sum(trees9$ba_m2), sum(trees9$ba_m2) / 50)

q <- trees9[, .(n = .N, ba = sum(ba_m2)), by = .(sqid = pmin(floor(gx / 20), 49) * 25 + pmin(floor(gy / 20), 24))]
q <- merge(data.table(sqid = 0:1249), q, by = "sqid", all.x = TRUE)
q[is.na(n), n := 0L][is.na(ba), ba := 0]
q[, `:=`(trees_ha = n / 0.04, ba_ha = ba / 0.04)]
if (requireNamespace("boot", quietly = TRUE) && requireNamespace("HDInterval", quietly = TRUE)) {
    set.seed(42)
    boot_hdi <- function(v) {
        b <- boot::boot(v, function(d, i) mean(d[i]), R = 4999)
        HDInterval::hdi(as.numeric(b$t), credMass = 0.95)
    }
    ci_t <- boot_hdi(q$trees_ha)
    ci_b <- boot_hdi(q$ba_ha)
    say("Density: %.0f trees/ha (95%% interval %.0f-%.0f); basal area: %.2f m2/ha (95%% interval %.2f-%.2f)",
        mean(q$trees_ha), ci_t[1], ci_t[2], mean(q$ba_ha), ci_b[1], ci_b[2])
} else {
    say("Install boot and HDInterval for the bootstrap intervals; density %.0f trees/ha, basal area %.2f m2/ha", mean(q$trees_ha), mean(q$ba_ha))
}

sp9 <- trees9[, .(individuals = .N, ba_m2 = sum(ba_m2)), by = .(sp, Latin, Family, Genus)]
n_tot <- nrow(trees9)
ba_tot <- sum(trees9$ba_m2)
half <- mean(sp9$individuals) / 2
say("Species below half the mean abundance (%.1f): %d (%.1f%%), with %s individuals (%.1f%%)",
    half, sp9[individuals < half, .N], 100 * sp9[individuals < half, .N] / nrow(sp9),
    fmt(sp9[individuals < half, sum(individuals)]), 100 * sp9[individuals < half, sum(individuals)] / n_tot)
say("Species with ten or fewer individuals: %d (%s individuals); singletons: %d",
    sp9[individuals <= 10, .N], fmt(sp9[individuals <= 10, sum(individuals)]), sp9[individuals == 1, .N])
say("Top species by individuals: %s", paste(sp9[order(-individuals)][1:3, sprintf("%s %s (%.2f%%)", Latin, fmt(individuals), 100 * individuals / n_tot)], collapse = "; "))
say("Top species by basal area: %s", paste(sp9[order(-ba_m2)][1:3, sprintf("%s %.2f m2 (%.2f%%)", Latin, ba_m2, 100 * ba_m2 / ba_tot)], collapse = "; "))
fam9 <- sp9[, .(species = .N, genera = uniqueN(Genus), individuals = sum(individuals), ba_m2 = sum(ba_m2)), by = Family]
say("Top families by individuals: %s", paste(fam9[order(-individuals)][1:3, sprintf("%s %s (%.2f%%)", Family, fmt(individuals), 100 * individuals / n_tot)], collapse = "; "))
say("Top families by basal area: %s", paste(fam9[order(-ba_m2)][1:3, sprintf("%s %.2f m2 (%.2f%%)", Family, ba_m2, 100 * ba_m2 / ba_tot)], collapse = "; "))
