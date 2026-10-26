# =============================================================================
# 2_data-description.R — the numbers given in the Data Description, Summary
# of Methods and Limitations and Challenges of the Data Paper, from the files
# of this package.
#
# Run from the folder that holds data/ and code/ (the repository folder), in
# R or with:  Rscript code/2_data-description.R
# Writes data-description_numbers.csv to ../manuscript_outputs when that
# folder exists, else to outputs/. Needs data.table.
# =============================================================================

suppressMessages(library(data.table))
source(file.path("code", "load-bci-tables.R"))
DATA <- "data"
OUT <- if (dir.exists("../manuscript_outputs")) "../manuscript_outputs" else "outputs"
dir.create(OUT, showWarnings = FALSE)

bci <- load_bci_tables(DATA)
spp <- as.data.table(bci$bci.spptable)
x <- rbindlist(lapply(bci[paste0("bci.stem", 1:9)], as.data.table))
post <- as.data.table(load_bci_posterior(DATA))
rm(bci)

out <- list()
add <- function(section, item, value) out[[length(out) + 1L]] <<- data.table(section, item, value = as.character(value))

# ---- Data Description ------------------------------------------------------------
taxa <- spp[sp != "uniden"]
add("Data Description", "trees", uniqueN(x$treeID))
add("Data Description", "stems (rows per census table)", x[CensusID == 1, .N])
add("Data Description", "trees with more than one stem", x[CensusID == 1, .N, by = treeID][N > 1, .N])
add("Data Description", "field records (DFstatus not NA)", sum(!is.na(x$DFstatus)))
add("Data Description", "measurements (dbh not NA)", sum(!is.na(x$dbh)))
add("Data Description", "taxa (species codes without uniden)", nrow(taxa))
add("Data Description", "genera", uniqueN(taxa$Genus))
add("Data Description", "families", uniqueN(taxa$Family))
add("Data Description", "trees never identified (uniden)", x[sp == "uniden", uniqueN(treeID)])
add("Data Description", "species codes renamed since 2019 (NameChange)", spp[!is.na(NameChange), .N])
add("Data Description", "trees reconstructed (StemPaths not NA)", x[!is.na(StemPaths), uniqueN(treeID)])
add("Data Description", "trees not reconstructed", uniqueN(x$treeID) - x[!is.na(StemPaths), uniqueN(treeID)])
add("Data Description", "trees with posterior samples", uniqueN(post$treeID))
add("Data Description", "rows of the posterior samples", nrow(post))
for (s in x[!is.na(DFstatus), .N, by = DFstatus][order(-N), DFstatus]) {
    add("Data Description", paste0("DFstatus records: ", s), x[DFstatus %in% s, .N])
}

# Census duration and the interval between consecutive censuses of a tree
# (mean date of its living stems in each census)
for (k in 1:9) {
    d <- x[CensusID == k, range(ExactDate, na.rm = TRUE)]
    add("Data Description", sprintf("census %d duration (years)", k), round(as.numeric(diff(d)) / 365.25, 2))
}
td <- x[Rstatus == "A", .(d = mean(date)), by = .(treeID, CensusID)]
setorder(td, treeID, CensusID)
td[, `:=`(c0 = shift(CensusID), d0 = shift(d)), by = treeID]
iv <- td[!is.na(c0) & c0 == CensusID - 1L, .(years = median((d - d0) / 365.25)), by = .(c0, CensusID)][order(c0)]
for (i in seq_len(nrow(iv))) {
    add("Data Description", sprintf("median interval, census %d to %d (years)", iv$c0[i], iv$CensusID[i]), round(iv$years[i], 2))
}

# ---- Limitations and Challenges ------------------------------------------------------------
a <- x[Rstatus == "A" & is.na(dbh)]
add("Limitations", "alive cells without a diameter", nrow(a))
add("Limitations", "of which stumps (DFstatus broken below)", a[DFstatus %in% "broken below", .N])
add("Limitations", "of which missed measurements", a[!DFstatus %in% "broken below", .N])
add("Limitations", "measurements away from 1.3 m (%)",
    round(100 * x[!is.na(dbh) & !is.na(hom) & abs(hom - 1.3) > 0.01, .N] / x[!is.na(dbh), .N], 1))
add("Limitations", "trees without coordinates", x[is.na(gx) | is.na(gy), uniqueN(treeID)])
add("Limitations", "trees without quadrat", x[is.na(quadrat), uniqueN(treeID)])

res <- rbindlist(out)
fwrite(res, file.path(OUT, "data-description_numbers.csv"))
print(res, nrows = 200)
cat("Written to", normalizePath(file.path(OUT, "data-description_numbers.csv")), "\n")
