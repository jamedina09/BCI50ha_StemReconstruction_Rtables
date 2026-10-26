# =============================================================================
# 1_tables.R — Tables 1–3 of the Data Paper, from the files of this package.
#   Table 1: the censuses (field dates, measurements, living trees)
#   Table 2: the data files (rows, columns, size)
#   Table 3: the columns of the census tables and their types
# Table 4 (the analytical status codes) holds definitions only.
#
# Run from the folder that holds data/ and code/ (the repository folder), in
# R or with:  Rscript code/1_tables.R
# Writes CSV files to ../manuscript_outputs when that folder exists, else to
# outputs/. Needs data.table.
# =============================================================================
rm(list = ls())

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
# Field dates: the 2.5% and 97.5% quantiles of the record dates (95% of the
# records were taken between them) and the median date, rather than the first
# and last dates, which a few early or late records can stretch (census 8
# starts in September 2013 although most of it was measured in 2015). type = 1
# returns dates on which records were taken.
date_q <- function(d, p) as.Date(quantile(as.numeric(d), p, na.rm = TRUE, type = 1, names = FALSE), origin = "1970-01-01")
tab1 <- x[, .(
    date_q025 = date_q(ExactDate, 0.025),
    median_date = date_q(ExactDate, 0.5),
    date_q975 = date_q(ExactDate, 0.975),
    # records = sum(!is.na(DFstatus)),
    # measurements = sum(!is.na(dbh)),
    # living_stems = sum(Rstatus == "A"),
    living_trees = uniqueN(treeID[Rstatus == "A"])
), by = CensusID][order(CensusID)]
tab1[, year := c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2023)[CensusID]]
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
    type = vapply(x, function(v) class(v)[1], character(1))
)
fwrite(tab3, file.path(OUT, "table-3_columns.csv"))
print(tab3)

cat("Tables written to", normalizePath(OUT), "\n")
