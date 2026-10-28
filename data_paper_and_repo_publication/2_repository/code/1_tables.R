# =============================================================================
# 1_tables.R — the figure and the tables of the Data Paper, from the files of
# this package. Each output holds everything the paper shows:
#   figure-1_tree-dynamics.png / .pdf   trees through the censuses: above zero the
#                                       trees with a living stem at 1.3 m (of
#                                       census 1 and new since then), below zero
#                                       the trees alive without one (broken below)
#                                       and the dead; the period of each census
#   figure-1_tree-dynamics.csv          the numbers and census dates of the figure
#   figure-1_caption.txt                its caption, written from those numbers
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
# A tree is one individual (treeID), whatever its number of stems. In each census,
# every tree recorded alive up to then is in exactly one of four groups:
#   living stem      at least one stem alive at 1.3 m: Rstatus A on a record that is
#                    not a stump. These are the trees counted in a census (above zero)
#   broken below     no living stem at 1.3 m, but a stump record: DFstatus "broken
#                    below" without a diameter. The field crew recorded the tree alive,
#                    broken below 1.3 m, so it has no stem at breast height
#   alive later      neither in this census (missing, or recorded dead), but a living
#                    stem or a stump record in a later census
#   dead             no living stem and no stump record in this census or in any later
#                    one (cumulative)
# "Recorded alive" means a living stem at 1.3 m or a stump record. Trees of census 1:
# recorded alive in census 1; new trees: first recorded alive after census 1.
#
# The census tables code a stump record by the stem at 1.3 m, which is gone: G when
# another stem of the tree is alive then or later, D when none is (A where the same
# stem is alive later). The figure counts these trees as alive, broken below; the
# codes of their stump records are written to the figure data.
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

st <- x[, .(treeID, CensusID, Rstatus, stump = DFstatus %in% "broken below" & is.na(dbh), DFstatus)]
st[, stem13 := Rstatus == "A" & !stump]
tc <- st[, .(
    stem13 = any(stem13), stump = any(stump), anyA = any(Rstatus == "A"),
    stump_A = any(stump & Rstatus == "A"), stump_G = any(stump & Rstatus == "G"), stump_D = any(stump & Rstatus == "D")
), by = .(treeID, CensusID)]
life <- tc[stem13 | stump, .(first = min(CensusID), last = max(CensusID)), by = treeID]
n_trees_all <- uniqueN(x$treeID)
n_no_sign <- n_trees_all - nrow(life)
tc <- life[tc, on = "treeID", nomatch = NULL][first <= CensusID]
tc[, cohort := fifelse(first == 1L, "census_1", "new")]
tc[, group := fcase(stem13, "living_stem", stump, "broken_below", last > CensusID, "alive_later", default = "dead")]
grid <- CJ(CensusID = 1:9, cohort = c("census_1", "new"), group = c("living_stem", "broken_below", "alive_later", "dead"))
cnt <- tc[, .N, by = .(CensusID, cohort, group)][grid, on = c("CensusID", "cohort", "group")]
cnt[is.na(N), N := 0L]
dyn <- dcast(cnt, CensusID ~ cohort + group, value.var = "N")
setnames(dyn, sub("^census_1_", "", sub("^(census_1|new)_(.*)$", "\\2_\\1", names(dyn))))
dyn <- merge(dyn, tc[, .(
    with_A_stem = sum(anyA),
    broken_below_coded_A = sum(group == "broken_below" & stump_A),
    broken_below_coded_G = sum(group == "broken_below" & stump_G & !stump_A),
    broken_below_coded_D = sum(group == "broken_below" & stump_D & !stump_A & !stump_G)
), by = CensusID], by = "CensusID")
dyn[, `:=`(
    living_stem = living_stem_census_1 + living_stem_new,
    broken_below = broken_below_census_1 + broken_below_new,
    alive_later = alive_later_census_1 + alive_later_new,
    recorded_census_1 = living_stem_census_1 + broken_below_census_1 + alive_later_census_1 + dead_census_1,
    recorded_new = living_stem_new + broken_below_new + alive_later_new + dead_new
)]
setorder(dyn, CensusID)
# Checks: every tree is in one group; the groups agree with the census tables
dead_tc <- life[last < 9L, .(CensusID = (last + 1L):9L), by = treeID]
dead_st <- st[dead_tc, on = c("treeID", "CensusID"), nomatch = NULL]
stopifnot(
    all(dyn$recorded_census_1 == dyn$living_stem_census_1[1]), # the trees of census 1, in every census
    all(dyn$recorded_census_1 + dyn$recorded_new == vapply(1:9, function(k) life[first <= k, .N], integer(1))),
    nrow(life) + n_no_sign == n_trees_all,
    dyn$alive_later[9] == 0, # no later census
    !is.unsorted(dyn$dead_census_1), !is.unsorted(dyn$dead_new), # the dead only accumulate
    nrow(dead_tc) == dyn[, sum(dead_census_1 + dead_new)],
    all(dead_st$Rstatus == "D"), !any(dead_st$stump), # a dead tree: all its stems D, no stump record
    all(dyn$broken_below == dyn$broken_below_coded_A + dyn$broken_below_coded_G + dyn$broken_below_coded_D),
    all(dyn$with_A_stem - dyn$living_stem == dyn$broken_below_coded_A) # trees whose only A stems are stumps
)
n_census_1 <- dyn$recorded_census_1[1]
dyn[, `:=`(
    pct_census_1_living_stem = round(100 * living_stem_census_1 / recorded_census_1, 1),
    pct_census_1_broken_below = round(100 * broken_below_census_1 / recorded_census_1, 1),
    pct_census_1_dead = round(100 * dead_census_1 / recorded_census_1, 1),
    pct_new_living_stem = fifelse(recorded_new > 0, round(100 * living_stem_new / recorded_new, 1), NA_real_),
    pct_new_broken_below = fifelse(recorded_new > 0, round(100 * broken_below_new / recorded_new, 1), NA_real_),
    pct_new_dead = fifelse(recorded_new > 0, round(100 * dead_new / recorded_new, 1), NA_real_)
)]
dyn <- merge(period, dyn, by = "CensusID")
dyn[, year := census_year[CensusID]]
fwrite(dyn, file.path(OUT, "figure-1_tree-dynamics.csv"))
print(dyn[, .(CensusID, living_stem, broken_below, alive_later, dead_census_1, dead_new, with_A_stem,
              broken_below_coded_A, broken_below_coded_G, broken_below_coded_D)])

draw_figure_1 <- function() {
    col_cohort <- "#0072B2"
    col_new <- "#D55E00"
    col_band <- "grey45"
    ink <- "grey15"
    k <- 1e3
    ymax <- 285
    ymin <- -285
    t_med <- dec_year(dyn$median_date)
    n <- nrow(dyn)
    big <- function(v) formatC(v, format = "d", big.mark = ",")
    pct <- function(v) sprintf("%.1f%%", v)
    op <- par(mar = c(4.6, 4.6, 1, 1), family = "sans", las = 1, col.axis = ink, col.lab = ink)
    on.exit(par(op))
    plot(NA, xlim = c(1980.9, 2025.3), ylim = c(ymin, ymax), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    abline(h = seq(-250, 250, 50), col = "grey88", lwd = 0.8)
    area <- function(lower, upper, col) polygon(c(t_med, rev(t_med)), c(rep(upper, length.out = n) / k, rev(rep(lower, length.out = n) / k)), col = col, border = NA)
    # above zero: trees with a living stem at 1.3 m
    area(0, dyn$living_stem_census_1, adjustcolor(col_cohort, 0.85))
    area(dyn$living_stem_census_1, dyn$living_stem, adjustcolor(col_new, 0.85))
    # below zero, from zero down: alive without a living stem at 1.3 m (broken below,
    # or alive later), then the dead (cumulative)
    b1 <- -(dyn$broken_below + dyn$alive_later)
    b2 <- b1 - dyn$dead_census_1
    b3 <- b2 - dyn$dead_new
    area(b1, 0, adjustcolor(col_band, 0.9))
    area(b2, b1, adjustcolor(col_cohort, 0.35))
    area(b3, b2, adjustcolor(col_new, 0.35))
    # period of each census: 95% of its records
    rect(dec_year(dyn$date_q025), ymin, dec_year(dyn$date_q975), ymax, col = adjustcolor("black", 0.14), border = NA)
    abline(h = 0, col = ink, lwd = 1)
    lines(t_med, dyn$living_stem / k, lwd = 2.2, col = ink)
    points(t_med, dyn$living_stem / k, pch = 21, bg = "white", col = ink, cex = 1.1, lwd = 1.6)
    x_total <- t_med
    x_total[1] <- x_total[1] + 0.5 # the first label would run over the left edge
    x_total[2] <- x_total[2] + 0.6 # and the second one over the first
    text(x_total, dyn$living_stem / k, big(dyn$living_stem), pos = 3, offset = 0.8, col = ink, cex = 0.72)
    # one call per census: axis() drops labels that it finds overlapping
    for (i in seq_len(n)) {
        axis(1, at = t_med[i], labels = paste0(dyn$CensusID[i], "\n", dyn$year[i]), tck = -0.02, padj = 0.6, col = "grey60", col.ticks = "grey60")
    }
    axis(2, at = seq(-250, 250, 50), labels = abs(seq(-250, 250, 50)), tck = -0.02, col = "grey60", col.ticks = "grey60")
    mtext("Census and year", 1, 3.6, col = ink)
    mtext("Trees (thousands)", 2, 3.2, las = 0, col = ink)
    mtext("living stem at 1.3 m", 2, 3.2, at = 175, las = 0, col = "grey40", cex = 0.75)
    mtext("no living stem at 1.3 m", 2, 3.2, at = -185, las = 0, col = "grey40", cex = 0.75)
    # Numbers: each coloured area carries the percentage of its own group (trees of
    # census 1; new trees) with a living stem at 1.3 m (above zero) or dead (below);
    # the grey band carries its number of trees. The last census also carries the
    # number of trees of each area.
    mid_c1 <- dyn$living_stem_census_1 / 2 / k
    mid_new <- (dyn$living_stem_census_1 + dyn$living_stem) / 2 / k
    mid_band <- b1 / 2 / k
    mid_d1 <- (b1 + b2) / 2 / k
    mid_d2 <- (b2 + b3) / 2 / k
    inner <- 2:(n - 1)
    inner_d2 <- 5:(n - 1) # the dead new trees are too few before census 5 to hold a label
    inner_band <- 3:(n - 1) # the band is too thin in census 2 to hold a label
    text(t_med[inner], mid_c1[inner], pct(dyn$pct_census_1_living_stem[inner]), col = "white", cex = 0.8, font = 2)
    text(t_med[inner], mid_new[inner], pct(dyn$pct_new_living_stem[inner]), col = "white", cex = 0.8, font = 2)
    text(t_med[inner_band], mid_band[inner_band], big(-b1[inner_band]), col = "white", cex = 0.62, font = 2)
    text(t_med[inner], mid_d1[inner], pct(dyn$pct_census_1_dead[inner]), col = ink, cex = 0.8, font = 2)
    text(t_med[inner_d2], mid_d2[inner_d2], pct(dyn$pct_new_dead[inner_d2]), col = ink, cex = 0.8, font = 2)
    last <- function(y, count, p, col) text(t_med[n] - 0.4, y, sprintf("%s\n(%s)", big(count), pct(p)), adj = c(1, 0.5), col = col, cex = 0.8, font = 2)
    last(mid_c1[n], dyn$living_stem_census_1[n], dyn$pct_census_1_living_stem[n], "white")
    last(mid_new[n], dyn$living_stem_new[n], dyn$pct_new_living_stem[n], "white")
    text(t_med[n] - 0.4, mid_band[n], big(-b1[n]), adj = c(1, 0.5), col = "white", cex = 0.62, font = 2)
    last(mid_d1[n], dyn$dead_census_1[n], dyn$pct_census_1_dead[n], ink)
    last(mid_d2[n], dyn$dead_new[n], dyn$pct_new_dead[n], ink)
    # key, in the free lower-left corner
    key <- function(y, col, label) {
        rect(1982, y - 5.5, 1983.6, y + 5.5, col = col, border = NA)
        text(1984.2, y, label, adj = c(0, 0.5), col = ink, cex = 0.78)
    }
    key(-150, adjustcolor(col_cohort, 0.85), "Trees of census 1: living stem at 1.3 m")
    key(-167, adjustcolor(col_new, 0.85), "Trees new since 1982: living stem at 1.3 m")
    key(-184, adjustcolor(col_band, 0.9), "Alive, no living stem at 1.3 m (broken below)")
    key(-201, adjustcolor(col_cohort, 0.35), "Trees of census 1: dead")
    key(-218, adjustcolor(col_new, 0.35), "Trees new since 1982: dead")
}
png(file.path(OUT, "figure-1_tree-dynamics.png"), width = 7.5, height = 7.4, units = "in", res = 300)
draw_figure_1()
dev.off()
pdf(file.path(OUT, "figure-1_tree-dynamics.pdf"), width = 7.5, height = 7.4)
draw_figure_1()
dev.off()

# Caption of Figure 1, from the numbers of the figure
cap_n <- function(v) formatC(v, format = "d", big.mark = ",")
l9 <- dyn[CensusID == 9]
band <- dyn[CensusID %in% 2:8, broken_below + alive_later]
bb_share <- dyn[CensusID %in% 2:8, 100 * broken_below / (broken_below + alive_later)]
caption <- paste0(
    "Figure 1. Trees of the Barro Colorado Island 50-ha plot through nine censuses, 1982–2023. ",
    "A tree is one individual, whatever its number of stems. In each census, every tree recorded alive up to then is in one of three groups, from the census tables: ",
    "with a living stem at 1.3 m, that is, a stem with Rstatus A on a record that is not a stump (above zero; total written above the line; these are the trees counted in a census); ",
    "alive without a living stem at 1.3 m (grey, next to zero): broken below 1.3 m, with a record whose DFstatus is broken below and that has no diameter (a stump, which the field crew recorded alive), ",
    "or, for a few trees, missing or recorded dead in that census and alive in a later one; ",
    "and dead (pale colours, cumulative): no living stem at 1.3 m and no broken-below record in that census or in any later one. ",
    "Blue: the ", cap_n(n_census_1), " trees of census 1; orange: trees first recorded alive after census 1. ",
    "Percentages: share of the trees of the same group, recorded alive up to that census, that have a living stem at 1.3 m (above zero) or are dead (below zero); the rest of each group is in the grey band, whose number of trees is written in it; percentages are rounded to 0.1%, so the three shares of a group add up to 100% within 0.1%. ",
    "In 2022–2023, ", cap_n(l9$living_stem), " trees had a living stem at 1.3 m (", cap_n(l9$living_stem_census_1), " of census 1, ", l9$pct_census_1_living_stem, "%; ",
    cap_n(l9$living_stem_new), " new, ", l9$pct_new_living_stem, "%), ", cap_n(l9$broken_below), " were broken below (", cap_n(l9$broken_below_census_1), " of census 1, ",
    cap_n(l9$broken_below_new), " new) and ", cap_n(l9$dead_census_1 + l9$dead_new), " were dead (", cap_n(l9$dead_census_1), " of census 1, ", l9$pct_census_1_dead, "%; ",
    cap_n(l9$dead_new), " new, ", l9$pct_new_dead, "%). ",
    "From census 2 to 8, the grey band holds ", cap_n(min(band)), "–", cap_n(max(band)), " trees, ", round(min(bb_share)), "–", round(max(bb_share)), "% of them broken below. ",
    "The census tables code a stump record by the stem at 1.3 m, which is gone: G when another stem of the tree is alive then or later, D when none is; ",
    "in 2022–2023 the ", cap_n(l9$broken_below), " broken-below trees are coded D. ",
    "Shaded vertical bands: period in which 95% of the records of each census were taken (2.5% to 97.5% quantiles of the record dates); each census is placed at its median record date. ",
    cap_n(n_no_sign), " trees with only dead or missing records are not shown."
)
writeLines(strwrap(caption, width = 100), file.path(OUT, "figure-1_caption.txt"))
cat("\n", caption, "\n", sep = "")

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
