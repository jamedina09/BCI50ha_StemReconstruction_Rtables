# =============================================================================
# 1_tables.R — Tables 1–4 of the Data Paper, from the files of this package.
#   Table 1: the censuses (field dates, records, measurements, living stems
#            and trees)
#   Table 2: the data files (objects, rows, columns, size)
#   Table 3: the columns of the census tables (type, missing values)
#   Table 4: the counts of each analytical status (Rstatus)
#
# Run from the folder that holds data/ and code/ (the repository folder), in
# R or with:  Rscript code/1_tables.R
# Writes CSV files to ../manuscript_outputs when that folder exists, else to
# outputs/. Needs data.table.
# =============================================================================

suppressMessages(library(data.table))
source(file.path("code", "load-bci-tables.R"))
DATA <- "data"
OUT <- if (dir.exists("../manuscript_outputs")) "../manuscript_outputs" else "outputs"
dir.create(OUT, showWarnings = FALSE)

bci <- load_bci_tables(DATA)
x <- rbindlist(lapply(bci[paste0("bci.stem", 1:9)], as.data.table))

# ---- Table 1: censuses -----------------------------------------------------------
# Records: rows with a field record (DFstatus). Measurements: records with a
# diameter. Living stems: Rstatus A; living trees: trees with a living stem.
tab1 <- x[, .(
    first_date = min(ExactDate, na.rm = TRUE),
    last_date = max(ExactDate, na.rm = TRUE),
    records = sum(!is.na(DFstatus)),
    measurements = sum(!is.na(dbh)),
    living_stems = sum(Rstatus == "A"),
    living_trees = uniqueN(treeID[Rstatus == "A"])
), by = CensusID][order(CensusID)]
tab1[, year := c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2022)[CensusID]]
setcolorder(tab1, c("CensusID", "year"))
fwrite(tab1, file.path(OUT, "table-1_censuses.csv"))
print(tab1)

# ---- Table 2: data files -----------------------------------------------------------
files <- c(paste0("bci.stem", 1:9, ".Rdata"), "bci.spptable.rdata", "posterior_sampled_paths.rds")
objects <- c(bci[paste0("bci.stem", 1:9)], list(bci$bci.spptable, load_bci_posterior(DATA)))
tab2 <- data.table(
    file = files,
    rows = vapply(objects, nrow, integer(1)),
    columns = vapply(objects, ncol, integer(1)),
    size_MB = round(file.size(file.path(DATA, files)) / 1e6, 2)
)
fwrite(tab2, file.path(OUT, "table-2_files.csv"))
print(tab2)
cat(sprintf("Total: %.1f MB\n", sum(tab2$size_MB)))

# ---- Table 3: columns of the census tables -----------------------------------------------
tab3 <- data.table(
    column = names(x),
    type = vapply(x, function(v) class(v)[1], character(1)),
    missing_values = vapply(x, function(v) sum(is.na(v)), integer(1))
)
fwrite(tab3, file.path(OUT, "table-3_columns.csv"))
print(tab3)

# ---- Table 4: analytical status ------------------------------------------------------------
tab4 <- x[, .(cells = .N), by = Rstatus][order(match(Rstatus, c("P", "A", "G", "D")))]
fwrite(tab4, file.path(OUT, "table-4_status-counts.csv"))
print(tab4)

cat("Tables written to", normalizePath(OUT), "\n")
