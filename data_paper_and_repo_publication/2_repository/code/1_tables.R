# =============================================================================
# 1_tables.R — the figure and the tables of the Data Paper, from the files of
# this package. Each output holds everything the paper shows:
#   figure-1_tree-dynamics.png / .pdf   living trees through the censuses (above
#                                       zero: trees of census 1 still alive and
#                                       trees new since then; below zero: the same
#                                       two groups recorded alive earlier and not
#                                       alive now) and the period of each census;
#                                       figure-1_tree-dynamics.csv holds its
#                                       numbers and census dates
#   table_columns-census-tables.csv     the columns of the census tables
#   table_columns-species-table.csv     the columns of the species table
#   table_status-codes.csv              the analytical status codes (Rstatus)
#
# Run from the folder that holds data/ and code/ (the repository folder), in
# R or with:  Rscript code/1_tables.R
# Writes to ../manuscript_outputs when that folder exists, else to outputs/.
# Needs data.table.
# =============================================================================
rm(list = ls())

suppressMessages(library(data.table))
source(file.path("code", "load-bci-tables.R"))
DATA <- "data"
OUT <- if (dir.exists("../manuscript_outputs")) "../manuscript_outputs" else "outputs"
dir.create(OUT, showWarnings = FALSE)

bci <- load_bci_tables(DATA)
x <- rbindlist(lapply(bci[paste0("bci.stem", 1:9)], as.data.table))
census_year <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2023)

# ---- Figure 1: dynamics of the trees through the censuses ------------------------
# Above zero, the living trees (Rstatus A) of each census, split by when the tree
# was first alive: the trees of census 1 that are still alive (the cohort of 1982)
# and the trees that appeared after it. Below zero, the trees recorded alive in an
# earlier census that are not alive in this one, split the same way (cumulative).
# Every tree recorded alive up to a census is therefore above or below zero. A tree
# is one individual (treeID), whatever its number of stems. The percentage of the
# census-1 trees that are not alive is written on the cohort.
#
# Census period: the 2.5% and 97.5% quantiles of the record dates (95% of the
# records of a census were taken between them) and the median date, rather than
# the first and last dates, which a few early or late records can stretch
# (census 8 starts in September 2013 although most of it was measured in 2015).
# type = 1 returns dates on which records were taken. Census k is placed at its
# median date.
date_q <- function(d, p) as.Date(quantile(as.numeric(d), p, na.rm = TRUE, type = 1, names = FALSE), origin = "1970-01-01")
dec_year <- function(d) as.numeric(format(d, "%Y")) + (as.numeric(format(d, "%j")) - 1) / 365.25
period <- x[, .(
    date_q025 = date_q(ExactDate, 0.025),
    median_date = date_q(ExactDate, 0.5),
    date_q975 = date_q(ExactDate, 0.975)
), by = CensusID][order(CensusID)]

alive <- unique(x[Rstatus == "A", .(treeID, CensusID)])
first_alive <- alive[, .(first = min(CensusID)), by = treeID]
dyn <- first_alive[alive, on = "treeID"][, .(
    living_trees = .N,
    from_census_1 = sum(first == 1L),
    new_since_census_1 = sum(first > 1L)
), by = CensusID][order(CensusID)]
stopifnot(identical(dyn$living_trees, x[Rstatus == "A", uniqueN(treeID), by = CensusID][order(CensusID)]$V1))
# Recorded alive up to census k, and not alive in census k
dyn[, recorded_census_1 := sum(first_alive$first == 1L)]
dyn[, recorded_new := vapply(CensusID, function(k) sum(first_alive$first > 1L & first_alive$first <= k), numeric(1))]
dyn[, not_alive_census_1 := recorded_census_1 - from_census_1]
dyn[, not_alive_new := recorded_new - new_since_census_1]
dyn[, recorded_census_1 := NULL]
dyn[, recorded_new := NULL]
dyn[, pct_census_1_not_alive := round(100 * not_alive_census_1 / (from_census_1 + not_alive_census_1), 1)]
dyn[, pct_new := round(100 * new_since_census_1 / living_trees, 1)]
dyn <- merge(period, dyn, by = "CensusID")
dyn[, year := census_year[CensusID]]
fwrite(dyn, file.path(OUT, "figure-1_tree-dynamics.csv"))
print(dyn)

draw_figure_1 <- function() {
    col_cohort <- "#0072B2"
    col_new <- "#D55E00"
    ink <- "grey15"
    k <- 1e3
    ymax <- 285
    ymin <- -285
    t_med <- dec_year(dyn$median_date)
    n <- nrow(dyn)
    op <- par(mar = c(4.6, 4.6, 1, 1), family = "sans", las = 1, col.axis = ink, col.lab = ink)
    on.exit(par(op))
    plot(NA, xlim = c(1980.9, 2025.3), ylim = c(ymin, ymax), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    abline(h = seq(-250, 250, 50), col = "grey88", lwd = 0.8)
    # living trees, above zero
    polygon(c(t_med, rev(t_med)), c(dyn$from_census_1 / k, rep(0, n)), col = adjustcolor(col_cohort, 0.85), border = NA)
    polygon(c(t_med, rev(t_med)), c(dyn$living_trees / k, rev(dyn$from_census_1 / k)), col = adjustcolor(col_new, 0.85), border = NA)
    # trees recorded alive before and not alive now, below zero (cumulative)
    polygon(c(t_med, rev(t_med)), c(-dyn$not_alive_census_1 / k, rep(0, n)), col = adjustcolor(col_cohort, 0.35), border = NA)
    polygon(c(t_med, rev(t_med)), c(-(dyn$not_alive_census_1 + dyn$not_alive_new) / k, rev(-dyn$not_alive_census_1 / k)), col = adjustcolor(col_new, 0.35), border = NA)
    # period of each census: 95% of its records
    rect(dec_year(dyn$date_q025), ymin, dec_year(dyn$date_q975), ymax, col = adjustcolor("black", 0.14), border = NA)
    abline(h = 0, col = ink, lwd = 1)
    lines(t_med, dyn$living_trees / k, lwd = 2.2, col = ink)
    points(t_med, dyn$living_trees / k, pch = 21, bg = "white", col = ink, cex = 1.1, lwd = 1.6)
    x_total <- t_med
    x_total[1] <- x_total[1] + 0.5 # the first label would run over the left edge
    x_total[2] <- x_total[2] + 0.6 # and the second one over the first
    text(x_total, dyn$living_trees / k, formatC(dyn$living_trees, format = "d", big.mark = ","), pos = 3, offset = 0.8, col = ink, cex = 0.72)
    # one call per census: axis() drops labels that it finds overlapping
    for (i in seq_len(n)) {
        axis(1, at = t_med[i], labels = paste0(dyn$CensusID[i], "\n", dyn$year[i]), tck = -0.02, padj = 0.6, col = "grey60", col.ticks = "grey60")
    }
    axis(2, at = seq(-250, 250, 50), labels = abs(seq(-250, 250, 50)), tck = -0.02, col = "grey60", col.ticks = "grey60")
    mtext("Census and year", 1, 3.6, col = ink)
    mtext("Trees (thousands)", 2, 3.2, las = 0, col = ink)
    mtext("not alive", 2, 3.2, at = -140, las = 0, col = "grey40", cex = 0.8)
    mtext("alive", 2, 3.2, at = 140, las = 0, col = "grey40", cex = 0.8)
    # percentage of the census-1 trees that are not alive, on the cohort
    text(t_med[2:(n - 1)], dyn$from_census_1[2:(n - 1)] / k, sprintf("%.1f%%", dyn$pct_census_1_not_alive[2:(n - 1)]), pos = 1, offset = 0.55, col = "white", cex = 0.8, font = 2)
    text(t_med[n] - 0.3, dyn$from_census_1[n] / k, sprintf("%.1f%%", dyn$pct_census_1_not_alive[n]), adj = c(1, 1.6), col = "white", cex = 0.8, font = 2)
    text(2013, 262, "All living trees", col = ink, cex = 0.95, font = 2)
    text(1982.2, 66, "Trees of census 1\nstill alive", adj = c(0, 0.5), col = "white", cex = 0.9, font = 2)
    i <- n - 1L
    text(t_med[i], (dyn$living_trees[i] + dyn$from_census_1[i]) / 2 / k, "Trees new since 1982\n(still alive)", col = "white", cex = 0.9, font = 2)
    text(1982.2, 25, "White percentages: share of the census-1 trees that are not alive", adj = c(0, 0.5), col = "white", cex = 0.72)
    text(1982.2, 10, "Shaded bands: period in which 95% of the records of each census were taken", adj = c(0, 0.5), col = "white", cex = 0.72)
    # below zero
    text(t_med[n - 3], -dyn$not_alive_census_1[n - 3] / k / 2, "Trees of census 1\nnot alive", col = ink, cex = 0.9, font = 2)
    j <- n - 1L
    text(t_med[j], -(dyn$not_alive_census_1[j] + dyn$not_alive_new[j] / 2) / k, "Trees new since 1982\nnot alive", col = ink, cex = 0.9, font = 2)
    text(t_med[n] - 0.4, -dyn$not_alive_census_1[n] / k / 2, formatC(dyn$not_alive_census_1[n], format = "d", big.mark = ","), adj = c(1, 0.5), col = ink, cex = 0.8, font = 2)
    text(t_med[n] - 0.4, -(dyn$not_alive_census_1[n] + dyn$not_alive_new[n] / 2) / k, formatC(dyn$not_alive_new[n], format = "d", big.mark = ","), adj = c(1, 0.5), col = ink, cex = 0.8, font = 2)
}
png(file.path(OUT, "figure-1_tree-dynamics.png"), width = 7.5, height = 6.6, units = "in", res = 300)
draw_figure_1()
dev.off()
pdf(file.path(OUT, "figure-1_tree-dynamics.pdf"), width = 7.5, height = 6.6)
draw_figure_1()
dev.off()

# ---- Columns of the census tables -----------------------------------------------------
# Type: from the data. Content: what the column holds.
census_columns <- c(
    CensusID = "Census, 1–9",
    treeID = "Tree; the same in every census",
    stemID = "Stem number within its tree (1, 2, …); treeID and stemID together identify a stem",
    tag = "Field tag of the tree, with leading zeros and sign kept (e.g. 000006, -05599)",
    StemTag = "Field tag of the stem, used since 2010; earlier records may show the tag the stem received later",
    sp = "Species code; joins the species table",
    quadrat = "20 × 20 m quadrat, XXYY: column (gx/20) and row (gy/20)",
    gx = "Location of the tree, in m east of the western edge (0–1000); one value per tree",
    gy = "Location of the tree, in m north of the southern edge (0–500); one value per tree",
    dbh = "Diameter as recorded, in mm, at height hom; empty when not measured",
    hom = "Height of measurement as recorded, in m",
    ExactDate = "Date of the field record; empty without a record",
    date = "Date of every row, in days since 1 January 1960: the recorded date, else the most common date of the tree, the quadrat or the census",
    codes = "Field codes as recorded, separated by ;",
    DFstatus = "Field status as recorded: alive, broken below, dead, stem dead, missing; empty without a record",
    Rstatus = "Analytical status: P, A, G or D (status codes)",
    StemPaths = "Key of the record in the reconstruction, unique within its tree; the posterior samples list the measured records up to the tree's anchor census (2010 for most trees); empty for trees not reconstructed"
)
stopifnot(identical(names(x), names(census_columns)))
tab_census <- data.table(
    Column = names(x),
    Type = vapply(x, function(v) class(v)[1], character(1)),
    Content = unname(census_columns)
)
fwrite(tab_census, file.path(OUT, "table_columns-census-tables.csv"))
print(tab_census)

# ---- Columns of the species table -------------------------------------------------------
spp <- as.data.table(bci$bci.spptable)
species_columns <- c(
    sp = "Species code; joins the census tables",
    Order = "Order",
    Family = "Family",
    Genus = "Genus",
    SpeciesName = "Species epithet; for a morphospecies, its label",
    Latin = "Scientific name: genus and epithet, then the infraspecific name when there is one",
    InfraspecificRank = "Rank below species (var. or subsp.); empty when none",
    InfraspecificEpithet = "Epithet below species; empty when none",
    Authority = "Authors of the name",
    Synonyms = "Synonyms and earlier names, separated by commas",
    Lifeform_RFoster = "Growth form after R. B. Foster",
    Lifeform_RPerez_SAguilar = "Growth form and height class after R. Pérez and S. Aguilar (in Spanish)",
    CommonName = "Common names",
    Herbarium = "Herbarium collections",
    NameCondit2019 = "Name of the code in Condit et al. (2019a, 2019b), only when it differs from Latin",
    NameChange = "Kind of change from that name",
    Notes = "What no other column says, e.g. that uniden is not a taxon"
)
stopifnot(identical(names(spp), names(species_columns)))
tab_species <- data.table(
    Column = names(spp),
    Type = vapply(spp, function(v) class(v)[1], character(1)),
    Content = unname(species_columns)
)
fwrite(tab_species, file.path(OUT, "table_columns-species-table.csv"))
print(tab_species)

# ---- Analytical status codes (Rstatus) ---------------------------------------------------
# Evidence of life is an alive record, or a diameter on a broken-below, missing
# or status-less record. A dead record is dead or stem dead, or broken below
# without a diameter.
tab_status <- data.table(
    Code = c("P", "A", "G", "D"),
    Meaning = c("prior", "alive", "dead stem", "dead tree"),
    Rule = c(
        "Before the stem's first record. A stem with only missing records is P in every census",
        "From the first record to the last evidence of life, also in censuses without a record or with a dead record in between. An alive cell without a diameter is a missed measurement or a stump",
        "Dead, while another stem of its tree is A in that census or later",
        "Dead, and every stem of its tree dead"
    )
)
stopifnot(setequal(tab_status$Code, unique(x$Rstatus)))
fwrite(tab_status, file.path(OUT, "table_status-codes.csv"))
print(tab_status)

cat("Written to", normalizePath(OUT), "\n")
