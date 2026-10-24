# ========================================================================
# SCRIPT: compare_viewfulltable_bci_old_and_oct13_2026.R
# PURPOSE: Reassess the three BCI `ViewFullTable` files kept in DATA/RAW:
#          what each one contains and how they differ, comparing every pair
#          of tables with the same method.
#
# Tables (name used in this script, file, what it is supposed to contain):
#   extract_2223  ViewFullTable_census2022-23.txt
#                 only the 2022-23 census
#   old           ViewFiles_bci_allcensuses/ViewFullTable_bci.csv
#                 all censuses
#   new           oct_13_2026/ViewFullTable_bci.csv
#                 all censuses, most up to date (received 2026-10-13)
#
# Table of Contents:
#   0. Setup
#   1. Configuration (paths of the tables, output folder, `check` and `expect`)
#   2. Load input files (everything as character, so values are compared
#      exactly as written in the files)
#   3. What each table contains (censuses, dates, rows, trees, stems, DBH)
#   4. Comparison helpers (rows keyed by DBHID)
#   5. Pairwise comparisons: extract_2223 vs old, old vs new, extract_2223 vs new
#   6. old vs new in detail: every changed cell and the species that were renamed
#   7. old vs new: rows that exist only in one version
#   8. Taxonomy names against ViewTaxonomy
#   9. Final check: compare DBH distributions using density (plots + KS test)
#  10. Verdicts and saved outputs (DATA/CHECKS)
#
# Main results (objects to look at):
#   - verdicts        : every expectation about the tables, TRUE / FALSE + detail
#   - rows_per_census : rows of each census in each of the three tables
#   - totals          : rows, trees, stems and sum of DBH of each table
#   - pairs_summary   : one row per pair of tables (rows shared, columns changed)
#   - changed_columns : old vs new, names of the columns whose values changed
#   - comparison      : old vs new, one row per column with the changed cells
#   - changes         : old vs new, one row per changed cell
#   - species_renamed : old vs new, species whose names changed
#   - rows_only_old / rows_only_new : records present in only one version
#   - report_densities_*  : DBH density comparison (L1/L2 distances, KS test)
#
# Notes:
#  - `DBHID` is the measurement key of ViewFullTable; rows are matched on it,
#    so the comparison does not depend on row order.
#  - NA vs NA counts as equal; NA vs value counts as a change.
#  - `check()` stops the script (broken input); `expect()` only reports: a ❌
#    from `expect()` is a difference between tables, which is what this script
#    is looking for.
# ========================================================================

# ---- 0. Setup ----
rm(list = ls())

library(data.table)
library(here)

setDTthreads(16) # tune this value for your machine

# ---- 1. Configuration ----
site <- "bci"
key_col <- "DBHID"

# TRUE: drop the 2020 census (PlotCensusNumber 9), which only the Oct 13 table
# has, so `new` and `old` hold the same censuses and only real differences
# remain. FALSE: compare the tables as delivered.
drop_census_9 <- TRUE

# Columns holding taxonomic names, and the codes that identify a species
taxonomy_cols <- c("Family", "Genus", "SpeciesName", "Subspecies")
species_key <- c("Mnemonic", "SpeciesID", "SubspeciesID")

file_extract <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "ViewFullTable_census2022-23.txt"
)
file_old <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses",
  paste0("ViewFullTable_", site, ".csv")
)
file_new <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026",
  paste0("ViewFullTable_", site, ".csv")
)
file_taxonomy <- here::here(
  "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses",
  paste0("ViewTaxonomy_", site, ".csv")
)

CHECKS_folder <- here::here("BCI_stem_reconstruction", "DATA", "CHECKS")
dir.create(CHECKS_folder, showWarnings = FALSE, recursive = TRUE)
out_prefix <- file.path(CHECKS_folder, "compare_viewfulltable_old_vs_oct13_2026") # old vs new
out_prefix_all <- file.path(CHECKS_folder, "compare_viewfulltable_all_tables") # the three tables
log_file <- paste0(out_prefix, "_log.txt")
cat("Comparison run:", format(Sys.time()), "\n", file = log_file)

# ---- Helper: check ----
# Loud check (the script is run line by line): prints ❌, raises an immediate
# warning, writes to the log, and only then stops.
check <- function(ok, msg) {
  if (isTRUE(all(ok))) {
    cat("✓", msg, "\n")
    return(invisible(TRUE))
  }
  full_msg <- paste0("CHECK FAILED: ", msg)
  cat("❌", full_msg, "\n")
  warning(full_msg, call. = FALSE, immediate. = TRUE)
  cat(full_msg, "\n", file = log_file, append = TRUE)
  stop(full_msg, call. = FALSE)
}

# ---- Helper: expect ----
# Records whether an expectation about the tables holds, in `verdicts`.
# Unlike `check()` it does not stop: a ❌ here is a difference between tables
# (a result), not a broken input. It is just as loud: ❌ line, immediate
# warning, and a line in the log.
verdicts <- data.table(expectation = character(), holds = logical(), detail = character())
expect <- function(ok, msg, detail = "") {
  ok <- isTRUE(all(ok))
  verdicts <<- rbind(verdicts[expectation != msg], data.table(expectation = msg, holds = ok, detail = detail))
  detail_txt <- if (nzchar(detail)) paste0(" [", detail, "]") else ""
  if (ok) {
    cat("✓", paste0(msg, detail_txt), "\n")
  } else {
    full_msg <- paste0("NOT TRUE: ", msg, detail_txt)
    cat("❌", full_msg, "\n")
    warning(full_msg, call. = FALSE, immediate. = TRUE)
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }
  invisible(ok)
}

# Number with thousands separator, for messages
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ---- 2. Load input files ----
check(file.exists(file_extract), paste("2022-23 extract exists:", file_extract))
check(file.exists(file_old), paste("old file exists:", file_old))
check(file.exists(file_new), paste("new file exists:", file_new))
check(file.exists(file_taxonomy), paste("taxonomy file exists:", file_taxonomy))

read_view <- function(file) {
  fread(
    file = file,
    sep = "\t",
    quote = "",
    colClasses = "character",
    na.strings = c("NA", "NULL", "")
  )
}

extract_2223 <- read_view(file_extract)
old <- read_view(file_old)
new <- read_view(file_new)
if (drop_census_9) new <- new[PlotCensusNumber != 9]
unique(old[, .(CensusID, PlotCensusNumber)])
unique(new[, .(CensusID, PlotCensusNumber)])

oldtreeid <- unique(old$TreeID)
newtreeid <- unique(new$TreeID)

length(oldtreeid)
length(newtreeid)

oldstemid <- unique(paste(old$TreeID, old$StemID, sep = "_"))
newstemid <- unique(paste(new$TreeID, new$StemID, sep = "_"))

length(oldstemid)
length(newstemid)

sum(as.numeric(old$DBH), na.rm = TRUE)
sum(as.numeric(new$DBH), na.rm = TRUE)

check(
  key_col %in% names(extract_2223) && key_col %in% names(old) && key_col %in% names(new),
  paste(key_col, "present in the three tables")
)
check(
  !anyNA(extract_2223[[key_col]]) && !anyNA(old[[key_col]]) && !anyNA(new[[key_col]]),
  paste(key_col, "has no missing values")
)
check(
  !anyDuplicated(extract_2223[[key_col]]) && !anyDuplicated(old[[key_col]]) && !anyDuplicated(new[[key_col]]),
  paste(key_col, "is unique within each table")
)

# ---- 3. What each table contains ----
# One row per census of each table: size and dates. This answers "which
# censuses are really in each file" before any comparison.
first_value <- function(v) if (all(is.na(v))) NA_character_ else min(v, na.rm = TRUE)
last_value <- function(v) if (all(is.na(v))) NA_character_ else max(v, na.rm = TRUE)

census_inventory <- function(x, label) {
  x[, .(
    dataset = label,
    rows = .N,
    trees = uniqueN(TreeID),
    stems = uniqueN(StemID),
    rows_with_DBH = sum(!is.na(DBH)),
    rows_without_date = sum(is.na(ExactDate)),
    first_date = first_value(ExactDate),
    last_date = last_value(ExactDate)
  ), keyby = .(PlotCensusNumber = as.integer(PlotCensusNumber), CensusID)]
}

census_inventory_all <- rbindlist(list(
  census_inventory(extract_2223, "extract_2223"),
  census_inventory(old, "old"),
  census_inventory(new, "new")
))
census_inventory_all[, dataset := factor(dataset, levels = c("extract_2223", "old", "new"))]
census_inventory_all

# Rows of each census in each table (NA = that census is not in the table)
rows_per_census <- dcast(census_inventory_all, PlotCensusNumber + CensusID ~ dataset, value.var = "rows")
rows_per_census <- merge(
  census_inventory_all[, .(first_date = min(first_date), last_date = max(last_date)), by = .(PlotCensusNumber, CensusID)],
  rows_per_census,
  by = c("PlotCensusNumber", "CensusID")
)
rows_per_census

expect(
  identical(names(extract_2223), names(old)) && identical(names(old), names(new)),
  "the three tables have the same columns in the same order",
  paste(ncol(new), "columns")
)

census_extract <- unique(extract_2223$CensusID)
years_extract <- unique(substr(extract_2223$ExactDate, 1, 4))
expect(
  length(census_extract) == 1 && all(years_extract %in% c("2022", "2023")),
  "extract_2223 contains only the 2022-23 census",
  paste0("CensusID ", paste(census_extract, collapse = ", "), "; years ", paste(sort(years_extract, na.last = TRUE), collapse = ", "))
)

census_only_old <- setdiff(old$CensusID, new$CensusID)
census_only_new <- setdiff(new$CensusID, old$CensusID)
expect(
  length(census_only_old) == 0 && length(census_only_new) == 0,
  "old and new contain the same censuses",
  paste0(
    "censuses only in old: ", if (length(census_only_old)) paste(census_only_old, collapse = ", ") else "none",
    "; only in new: ",
    if (length(census_only_new)) {
      census_inventory_all[dataset == "new" & CensusID %in% census_only_new, paste(sprintf(
        "PlotCensusNumber %s (CensusID %s, %s to %s, %s rows)", PlotCensusNumber, CensusID, first_date, last_date, fmt(rows)
      ), collapse = "; ")]
    } else {
      "none"
    }
  )
)

# Totals per table: the quick figures of Section 2 (trees, stems, sum of DBH)
# for the three tables, plus `new` restricted to the censuses that `old`
# contains. If the only difference in records is an extra census in `new`,
# that last row equals the row of `old`.
table_totals <- function(x, label) {
  data.table(
    dataset = label,
    rows = nrow(x),
    trees = uniqueN(x$TreeID),
    stems = uniqueN(paste(x$TreeID, x$StemID, sep = "_")),
    rows_with_DBH = sum(!is.na(x$DBH)),
    sum_DBH = sum(as.numeric(x$DBH), na.rm = TRUE)
  )
}
totals <- rbindlist(list(
  table_totals(extract_2223, "extract_2223"),
  table_totals(old, "old"),
  table_totals(new, "new"),
  table_totals(new[CensusID %in% old$CensusID], "new (only the censuses of old)")
))
totals

expect(
  isTRUE(all.equal(
    unlist(totals[dataset == "old", -"dataset"]),
    unlist(totals[dataset == "new (only the censuses of old)", -"dataset"])
  )),
  "new restricted to the censuses of old has the same rows, trees, stems and sum of DBH as old",
  totals[dataset %in% c("old", "new"), paste0(dataset, ": ", fmt(rows), " rows, sum of DBH ", fmt(sum_DBH), collapse = "; ")]
)

# ---- 4. Comparison helpers ----
# TRUE where a cell differs (NA vs NA is equal, NA vs value is a change)
cell_differs <- function(x, y) {
  (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}

# TRUE where both values are numbers that are equal (e.g. "1.30" vs "1.3")
same_number <- function(x, y) {
  xn <- suppressWarnings(as.numeric(x))
  yn <- suppressWarnings(as.numeric(y))
  !is.na(xn) & !is.na(yn) & xn == yn
}

# ---- Helper: compare_tables ----
# Purpose: Compare two ViewFullTable versions. Rows are matched on `key_col`,
#          so the result does not depend on row order or on the number of rows.
# Args: x, y - tables read with `read_view()`; label_x, label_y - their names
#       in the outputs
# Returns: a list with
#   structure       - rows and columns shared / present in only one table
#   comparison      - one row per shared column: cells that differ on the
#                     rows shared by both tables
#   changed_columns - columns with at least one differing cell
#   ix, iy          - positions of the shared records in x and in y (aligned)
#   ids_only_x, ids_only_y - keys of the records present in only one table
compare_tables <- function(x, y, label_x, label_y) {
  cols_only_x <- setdiff(names(x), names(y))
  cols_only_y <- setdiff(names(y), names(x))
  cols_shared <- setdiff(intersect(names(x), names(y)), key_col)

  ids_only_x <- setdiff(x[[key_col]], y[[key_col]])
  ids_only_y <- setdiff(y[[key_col]], x[[key_col]])
  ids_shared <- intersect(x[[key_col]], y[[key_col]])
  # Positions that put the shared records of x and of y in the same order
  ix <- match(ids_shared, x[[key_col]])
  iy <- match(ids_shared, y[[key_col]])

  structure_summary <- data.table(
    item = c(
      paste("rows", label_x), paste("rows", label_y), "rows shared",
      paste("rows only in", label_x), paste("rows only in", label_y),
      paste("columns only in", label_x), paste("columns only in", label_y)
    ),
    value = c(
      nrow(x), nrow(y), length(ids_shared), length(ids_only_x), length(ids_only_y),
      paste(cols_only_x, collapse = ", "), paste(cols_only_y, collapse = ", ")
    )
  )

  comparison <- rbindlist(lapply(cols_shared, function(col) {
    xv <- x[[col]][ix]
    yv <- y[[col]][iy]
    d <- cell_differs(xv, yv)
    data.table(
      column = col,
      n_changed = sum(d),
      pct_changed = round(100 * mean(d), 4),
      n_value_to_NA = sum(d & !is.na(xv) & is.na(yv)),
      n_NA_to_value = sum(d & is.na(xv) & !is.na(yv)),
      n_value_to_value = sum(d & !is.na(xv) & !is.na(yv)),
      n_format_only = sum(same_number(xv[d], yv[d])) # same number, written differently
    )
  }))
  comparison[, result := fifelse(n_changed == 0, "Identical", "Different")]
  setorder(comparison, -n_changed)

  list(
    label_x = label_x, label_y = label_y,
    structure = structure_summary,
    comparison = comparison,
    changed_columns = comparison[result == "Different", column],
    ix = ix, iy = iy,
    ids_only_x = ids_only_x, ids_only_y = ids_only_y
  )
}

# Rows per census among the records `ids` of table x (used for the records
# that exist in only one table of a pair)
rows_by_census <- function(x, ids) {
  x[get(key_col) %in% ids, .(rows = .N), keyby = .(PlotCensusNumber = as.integer(PlotCensusNumber), CensusID)]
}

# "Family (27,707 cells), Genus (...)" for the columns `cols` of a comparison
changed_string <- function(cmp, cols = cmp$changed_columns) {
  if (length(cols) == 0) {
    return("no column differs")
  }
  paste0(cols, " (", fmt(cmp$comparison[match(cols, column), n_changed]), " cells)", collapse = ", ")
}

# ---- Helper: expect_pair ----
# The same four expectations for any pair of tables (results go to `verdicts`):
#   1. same records in the censuses that both tables contain
#   2. every column agrees on the shared records
#   3. every column other than the taxonomy names agrees
#   4. DBH and HOM agree
expect_pair <- function(cmp, x, y) {
  pair <- paste(cmp$label_x, "vs", cmp$label_y)
  census_shared <- intersect(x$CensusID, y$CensusID)
  ids_x <- x[[key_col]][x$CensusID %in% census_shared]
  ids_y <- y[[key_col]][y$CensusID %in% census_shared]
  expect(
    setequal(ids_x, ids_y),
    paste0(pair, ": same records in the censuses that both tables contain"),
    sprintf(
      "%s census(es); %s records in %s, %s in %s",
      length(census_shared), fmt(length(ids_x)), cmp$label_x, fmt(length(ids_y)), cmp$label_y
    )
  )
  expect(
    length(cmp$changed_columns) == 0,
    paste0(pair, ": every column agrees on the shared records"),
    changed_string(cmp)
  )
  not_taxonomy <- setdiff(cmp$changed_columns, taxonomy_cols)
  expect(
    length(not_taxonomy) == 0,
    paste0(pair, ": every column other than the taxonomy names agrees on the shared records"),
    changed_string(cmp, not_taxonomy)
  )
  dbh_changed <- intersect(cmp$changed_columns, c("DBH", "HOM"))
  expect(
    length(dbh_changed) == 0,
    paste0(pair, ": DBH and HOM agree on the shared records"),
    changed_string(cmp, dbh_changed)
  )
}

# ---- 5. Pairwise comparisons ----

# ---- 5a. extract_2223 vs old ----
# The comparison of compare_tables1_and_tables2.R, redone on the key: is the
# 2022-23 extract the same as the 2022-23 census inside the all-census table?
cmp_extract_old <- compare_tables(extract_2223, old, "extract_2223", "old")
cmp_extract_old$structure
cmp_extract_old$comparison[result == "Different"] # empty = no column differs
rows_by_census(extract_2223, cmp_extract_old$ids_only_x) # records only in the extract
rows_by_census(old, cmp_extract_old$ids_only_y) # records only in old (its other censuses)
expect_pair(cmp_extract_old, extract_2223, old)

# ---- 5b. old vs new ----
# Structure: columns and rows
cmp_old_new <- compare_tables(old, new, "old", "new")
structure_summary <- cmp_old_new$structure
structure_summary

ids_only_old <- cmp_old_new$ids_only_x
ids_only_new <- cmp_old_new$ids_only_y
rows_by_census(old, ids_only_old) # records only in old
rows_by_census(new, ids_only_new) # records only in new

# Column-by-column comparison (shared rows)
# Align both versions on the key so row i is the same measurement in both.
old_s <- old[cmp_old_new$ix]
new_s <- new[cmp_old_new$iy]
check(identical(old_s[[key_col]], new_s[[key_col]]), "shared rows aligned on the key")
cols_shared <- setdiff(intersect(names(old), names(new)), key_col)

comparison <- cmp_old_new$comparison

# >>> RESULT: columns that have changed <<<
changed_columns <- cmp_old_new$changed_columns
changed_columns
comparison[result == "Different"]
expect_pair(cmp_old_new, old, new)

# ---- 5c. extract_2223 vs new ----
cmp_extract_new <- compare_tables(extract_2223, new, "extract_2223", "new")
cmp_extract_new$structure
cmp_extract_new$comparison[result == "Different"]
rows_by_census(extract_2223, cmp_extract_new$ids_only_x) # records only in the extract
rows_by_census(new, cmp_extract_new$ids_only_y) # records only in new (its other censuses)
expect_pair(cmp_extract_new, extract_2223, new)

# ---- 5d. The three pairs side by side ----
pairs <- list(cmp_extract_old, cmp_old_new, cmp_extract_new)
pairs_summary <- rbindlist(lapply(pairs, function(p) {
  data.table(
    pair = paste(p$label_x, "vs", p$label_y),
    rows_shared = length(p$ix),
    rows_only_first = length(p$ids_only_x),
    rows_only_second = length(p$ids_only_y),
    n_columns_changed = length(p$changed_columns),
    columns_changed = paste(p$changed_columns, collapse = ", ")
  )
}))
pairs_summary

pairs_columns <- rbindlist(lapply(pairs, function(p) {
  cbind(pair = paste(p$label_x, "vs", p$label_y), p$comparison)
}))
pairs_columns[result == "Different"]

# ---- 6. old vs new in detail: every changed cell ----
id_cols <- intersect(c("TreeID", "StemID", "Tag", "StemTag", "CensusID", "PlotCensusNumber"), names(new_s))

changes <- rbindlist(lapply(changed_columns, function(col) {
  d <- which(cell_differs(old_s[[col]], new_s[[col]]))
  cbind(
    new_s[d, c(key_col, id_cols), with = FALSE],
    data.table(column = col, old_value = old_s[[col]][d], new_value = new_s[[col]][d])
  )
}))
changes

# Most frequent old -> new changes per column (quick view of what happened)
changes_top <- changes[, .N, by = .(column, old_value, new_value)][order(column, -N)]
changes_top[, head(.SD, 10), by = column]

# Changed cells per column and census
if (nrow(changes) > 0) dcast(changes, column ~ PlotCensusNumber, value.var = "DBHID", fun.aggregate = length)

# Species that were renamed: the taxonomy names of every species code in old
# and in new, side by side
species_names <- function(x) unique(x[, c(species_key, taxonomy_cols), with = FALSE])
species_old <- species_names(old)
species_new <- species_names(new)
check(
  !anyDuplicated(species_old, by = species_key) && !anyDuplicated(species_new, by = species_key),
  "each species code has a single set of names within old and within new"
)

species_compared <- merge(species_old, species_new, by = species_key, all = TRUE, suffixes = c("_old", "_new"))
species_compared <- merge(species_compared, old[, .(n_rows_old = .N), by = species_key], by = species_key, all.x = TRUE)
species_compared <- merge(species_compared, new[, .(n_rows_new = .N), by = species_key], by = species_key, all.x = TRUE)
species_compared[, renamed := !is.na(n_rows_old) & !is.na(n_rows_new) & (
  cell_differs(Family_old, Family_new) | cell_differs(Genus_old, Genus_new) |
    cell_differs(SpeciesName_old, SpeciesName_new) | cell_differs(Subspecies_old, Subspecies_new)
)]

species_only_one <- species_compared[is.na(n_rows_old) | is.na(n_rows_new)] # codes present in only one table
species_only_one

species_renamed <- species_compared[renamed == TRUE][order(-n_rows_old)]
species_renamed[, name_old := trimws(paste(Genus_old, SpeciesName_old, fifelse(is.na(Subspecies_old), "", Subspecies_old)))]
species_renamed[, name_new := trimws(paste(Genus_new, SpeciesName_new, fifelse(is.na(Subspecies_new), "", Subspecies_new)))]
species_renamed[, .(Mnemonic, name_old, name_new, Family_old, Family_new, n_rows_old)]

expect(
  nrow(species_only_one) == 0,
  "old vs new: same species codes (Mnemonic, SpeciesID, SubspeciesID)",
  sprintf("%s codes; %s of them with different names", nrow(species_compared), nrow(species_renamed))
)

# ---- 7. old vs new: rows that exist only in one version ----
rows_only_old <- old[get(key_col) %in% ids_only_old]
rows_only_new <- new[get(key_col) %in% ids_only_new]

rows_only_old[, .N, keyby = .(CensusID, PlotCensusNumber, Status)]
rows_only_new[, .N, keyby = .(CensusID, PlotCensusNumber, Status)]

# Where and when were the rows only in new recorded, and do their trees and
# stems already exist in old?
if (nrow(rows_only_new) > 0) {
  rows_only_new[, .(
    rows = .N,
    trees = uniqueN(TreeID),
    stems = uniqueN(StemID),
    trees_not_in_old = uniqueN(TreeID[!TreeID %in% old$TreeID]),
    stems_not_in_old = uniqueN(StemID[!StemID %in% old$StemID]),
    quadrats = uniqueN(QuadratName),
    first_date = first_value(ExactDate),
    last_date = last_value(ExactDate),
    PX_min = min(as.numeric(PX), na.rm = TRUE),
    PX_max = max(as.numeric(PX), na.rm = TRUE),
    PY_min = min(as.numeric(PY), na.rm = TRUE),
    PY_max = max(as.numeric(PY), na.rm = TRUE)
  ), keyby = .(PlotCensusNumber, CensusID)]
}
uniqueN(new[!is.na(QuadratName), QuadratName]) # quadrats in the whole plot, for reference

# ---- 8. Taxonomy names against ViewTaxonomy ----
# Which of the three tables carries the names of ViewTaxonomy? Species are
# matched on SpeciesID and subspecies on SubspeciesID (the IDs did not change).
taxonomy <- read_view(file_taxonomy)
taxonomy_species <- unique(taxonomy[, .(SpeciesID, Family_tax = Family, Genus_tax = Genus, SpeciesName_tax = SpeciesName, in_tax = TRUE)])
taxonomy_subspecies <- unique(taxonomy[!is.na(SubspeciesID), .(SubspeciesID, Subspecies_tax = Subspecies, in_tax = TRUE)])
check(!anyDuplicated(taxonomy_species$SpeciesID), "ViewTaxonomy has one Family/Genus/SpeciesName per SpeciesID")
check(!anyDuplicated(taxonomy_subspecies$SubspeciesID), "ViewTaxonomy has one Subspecies name per SubspeciesID")

compare_with_taxonomy <- function(x, label) {
  sp <- merge(unique(x[, .(SpeciesID, Family, Genus, SpeciesName)]), taxonomy_species, by = "SpeciesID", all.x = TRUE)
  sub <- merge(unique(x[!is.na(SubspeciesID), .(SubspeciesID, Subspecies)]), taxonomy_subspecies, by = "SubspeciesID", all.x = TRUE)
  found <- !is.na(sp$in_tax)
  data.table(
    dataset = label,
    species = nrow(sp),
    species_not_in_ViewTaxonomy = sum(!found),
    Family_differs = sum(found & cell_differs(sp$Family, sp$Family_tax)),
    Genus_differs = sum(found & cell_differs(sp$Genus, sp$Genus_tax)),
    SpeciesName_differs = sum(found & cell_differs(sp$SpeciesName, sp$SpeciesName_tax)),
    species_with_any_difference = sum(found & (
      cell_differs(sp$Family, sp$Family_tax) | cell_differs(sp$Genus, sp$Genus_tax) |
        cell_differs(sp$SpeciesName, sp$SpeciesName_tax)
    )),
    subspecies = nrow(sub),
    subspecies_not_in_ViewTaxonomy = sum(is.na(sub$in_tax)),
    Subspecies_differs = sum(!is.na(sub$in_tax) & cell_differs(sub$Subspecies, sub$Subspecies_tax))
  )
}

taxonomy_check <- rbindlist(list(
  compare_with_taxonomy(extract_2223, "extract_2223"),
  compare_with_taxonomy(old, "old"),
  compare_with_taxonomy(new, "new")
))
taxonomy_check

for (i in seq_len(nrow(taxonomy_check))) {
  expect(
    taxonomy_check[i, species_not_in_ViewTaxonomy + subspecies_not_in_ViewTaxonomy + species_with_any_difference + Subspecies_differs] == 0,
    paste0(taxonomy_check$dataset[i], ": taxonomy names are those of ViewTaxonomy_", site, ".csv"),
    taxonomy_check[i, sprintf(
      "%s of %s species and %s of %s subspecies named differently",
      species_with_any_difference, species, Subspecies_differs, subspecies
    )]
  )
}

# ---- 9. Final check: compare DBH distributions using density ----
# DBH (and HOM) are part of the column-by-column comparison in Section 5.
# Here the DBH distributions of both versions are compared visually and with
# L1/L2 distances plus a Kolmogorov-Smirnov test.
dbh_cols <- intersect(c("DBH", "HOM"), cols_shared)
comparison[column %in% dbh_cols]

# DBH as numeric (the tables were loaded as character)
dbh_old_shared <- as.numeric(old_s$DBH) # old version, rows shared by both
dbh_new_shared <- as.numeric(new_s$DBH) # new version, rows shared by both
dbh_old_all <- as.numeric(old$DBH) # old version, all rows
dbh_new_all <- as.numeric(new$DBH) # new version, all rows (includes rows only in new)
dbh_only_new <- as.numeric(rows_only_new$DBH) # rows that exist only in the new version
# No rows only in new (e.g. `drop_census_9` is TRUE): nothing to compare, and
# `density()` cannot run on an empty vector
has_dbh_only_new <- any(!is.na(dbh_only_new))
# Median / max that give NA instead of a warning when there is no DBH
dbh_stat <- function(f, x) if (any(!is.na(x))) f(x, na.rm = TRUE) else NA_real_

# Measured DBH per version (rows with a DBH value) and largest cell difference
data.table(
  set = c("old (shared rows)", "new (shared rows)", "old (all rows)", "new (all rows)", "rows only in new"),
  n_rows = c(length(dbh_old_shared), length(dbh_new_shared), length(dbh_old_all), length(dbh_new_all), length(dbh_only_new)),
  n_with_DBH = c(sum(!is.na(dbh_old_shared)), sum(!is.na(dbh_new_shared)), sum(!is.na(dbh_old_all)), sum(!is.na(dbh_new_all)), sum(!is.na(dbh_only_new))),
  median_DBH = c(median(dbh_old_shared, na.rm = TRUE), median(dbh_new_shared, na.rm = TRUE), median(dbh_old_all, na.rm = TRUE), median(dbh_new_all, na.rm = TRUE), dbh_stat(median, dbh_only_new)),
  max_DBH = c(max(dbh_old_shared, na.rm = TRUE), max(dbh_new_shared, na.rm = TRUE), max(dbh_old_all, na.rm = TRUE), max(dbh_new_all, na.rm = TRUE), dbh_stat(max, dbh_only_new))
)
max(abs(dbh_old_shared - dbh_new_shared), na.rm = TRUE) # 0 = DBH identical on shared rows

# ---- Helper: compare_densities ----
# Purpose: Compute and plot density estimates for two numeric vectors and
#          report L1/L2 differences plus a Kolmogorov-Smirnov test.
# Args:
#   x1, x2 - numeric vectors (NA will be removed)
#   n_grid - grid size for density interpolation (default 2048)
#   plot_title - main title for the plots
# Returns: a list with L1/L2 distances, KS test object, density grids and
#          difference curve and summaries.
compare_densities <- function(x1, x2, n_grid = 2048, plot_title = "Density Comparison") {
  # Remove NAs
  x1 <- na.omit(x1)
  x2 <- na.omit(x2)
  # Compute densities
  d1 <- density(x1, n = n_grid)
  d2 <- density(x2, n = n_grid)
  # Align densities on same grid
  x_grid <- seq(min(c(d1$x, d2$x)), max(c(d1$x, d2$x)), length.out = n_grid)
  d1y <- approx(d1$x, d1$y, xout = x_grid, rule = 2)$y
  d2y <- approx(d2$x, d2$y, xout = x_grid, rule = 2)$y
  # Compute difference curves
  diff_curve <- d1y - d2y
  # Compute L1 and L2 differences
  library(pracma)
  L1_diff <- trapz(x_grid, abs(diff_curve))
  L2_diff <- sqrt(trapz(x_grid, diff_curve^2))
  # Kolmogorov-Smirnov test
  ks <- ks.test(x1, x2)
  # Plot densities and difference
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  plot(x_grid, d1y,
    type = "l", col = "blue", lwd = 2,
    main = paste(plot_title, "- Densities"), xlab = "Value", ylab = "Density"
  )
  lines(x_grid, d2y, col = "red", lwd = 2)
  legend("topright", legend = c("x1", "x2"), col = c("blue", "red"), lwd = 2)
  plot(x_grid, diff_curve,
    type = "l", col = "darkgreen", lwd = 2,
    main = "Difference Curve (x1 - x2)", xlab = "Value", ylab = "Difference"
  )
  abline(h = 0, col = "gray", lty = 2)
  # Return full report as a list
  report <- list(
    L1_difference = L1_diff,
    L2_difference = L2_diff,
    KS_test = ks,
    density_grid = x_grid,
    density_x1 = d1y,
    density_x2 = d2y,
    difference_curve = diff_curve,
    summary_x1 = summary(x1),
    summary_x2 = summary(x2)
  )
  return(report)
}

# ---- Helper: save_density_plot ----
# Purpose: Redraw the two panels of `compare_densities()` from its report and
#          save them as a PNG (no recomputation of the densities).
save_density_plot <- function(report, file, plot_title = "Density Comparison") {
  png(file, width = 1600, height = 1600, res = 200)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  plot(report$density_grid, report$density_x1,
    type = "l", col = "blue", lwd = 2,
    main = paste(plot_title, "- Densities"), xlab = "Value", ylab = "Density"
  )
  lines(report$density_grid, report$density_x2, col = "red", lwd = 2)
  legend("topright", legend = c("x1", "x2"), col = c("blue", "red"), lwd = 2)
  plot(report$density_grid, report$difference_curve,
    type = "l", col = "darkgreen", lwd = 2,
    main = "Difference Curve (x1 - x2)", xlab = "Value", ylab = "Difference"
  )
  abline(h = 0, col = "gray", lty = 2)
  dev.off()
  invisible(file)
}

# Run density comparison and print summary statistics
# (x1 = old version, x2 = new version; each plot is also saved to DATA/CHECKS)

# a) Shared rows: old vs new DBH of the same measurements
report_densities_shared <- compare_densities(dbh_old_shared, dbh_new_shared, plot_title = "DBH old vs new (shared rows)")
save_density_plot(report_densities_shared, paste0(out_prefix, "_dbh_density_shared_rows.png"), plot_title = "DBH old vs new (shared rows)")
report_densities_shared$L1_difference
report_densities_shared$L2_difference
report_densities_shared$KS_test

# b) All rows: full old table vs full new table (new includes the added rows)
report_densities_all <- compare_densities(dbh_old_all, dbh_new_all, plot_title = "DBH old vs new (all rows)")
save_density_plot(report_densities_all, paste0(out_prefix, "_dbh_density_all_rows.png"), plot_title = "DBH old vs new (all rows)")
report_densities_all$L1_difference
report_densities_all$L2_difference
report_densities_all$KS_test

# c) Rows only in new vs the old table (are the added records comparable?)
#    Skipped when there are no rows only in new (`drop_census_9` is TRUE)
if (has_dbh_only_new) {
  report_densities_only_new <- compare_densities(dbh_old_all, dbh_only_new, plot_title = "DBH old vs rows only in new")
  save_density_plot(report_densities_only_new, paste0(out_prefix, "_dbh_density_rows_only_new.png"), plot_title = "DBH old vs rows only in new")
  report_densities_only_new$L1_difference
  report_densities_only_new$L2_difference
  report_densities_only_new$KS_test
} else {
  cat("No DBH in rows only in new: density comparison c) skipped\n")
}

# d) Same three comparisons on log10(DBH): DBH is strongly right-skewed, so on
#    the raw scale the curves collapse near zero; the log scale shows the shape.
report_log_shared <- compare_densities(log10(dbh_old_shared), log10(dbh_new_shared), plot_title = "log10 DBH old vs new (shared rows)")
save_density_plot(report_log_shared, paste0(out_prefix, "_dbh_density_log10_shared_rows.png"), plot_title = "log10 DBH old vs new (shared rows)")

report_log_all <- compare_densities(log10(dbh_old_all), log10(dbh_new_all), plot_title = "log10 DBH old vs new (all rows)")
save_density_plot(report_log_all, paste0(out_prefix, "_dbh_density_log10_all_rows.png"), plot_title = "log10 DBH old vs new (all rows)")

if (has_dbh_only_new) {
  report_log_only_new <- compare_densities(log10(dbh_old_all), log10(dbh_only_new), plot_title = "log10 DBH old vs rows only in new")
  save_density_plot(report_log_only_new, paste0(out_prefix, "_dbh_density_log10_rows_only_new.png"), plot_title = "log10 DBH old vs rows only in new")
}

# e) The check of compare_tables1_and_tables2.R: the 2022-23 extract against the
#    same census inside old, and inside new (x1 = extract, x2 = full table)
dbh_extract <- as.numeric(extract_2223$DBH)
dbh_old_census_extract <- as.numeric(old[CensusID %in% census_extract, DBH])
dbh_new_census_extract <- as.numeric(new[CensusID %in% census_extract, DBH])

report_densities_extract_old <- compare_densities(dbh_extract, dbh_old_census_extract, plot_title = "DBH extract 2022-23 vs old (same census)")
save_density_plot(report_densities_extract_old, paste0(out_prefix_all, "_dbh_density_extract_2223_vs_old.png"), plot_title = "DBH extract 2022-23 vs old (same census)")
report_densities_extract_old$L1_difference
report_densities_extract_old$L2_difference
report_densities_extract_old$KS_test

report_densities_extract_new <- compare_densities(dbh_extract, dbh_new_census_extract, plot_title = "DBH extract 2022-23 vs new (same census)")
save_density_plot(report_densities_extract_new, paste0(out_prefix_all, "_dbh_density_extract_2223_vs_new.png"), plot_title = "DBH extract 2022-23 vs new (same census)")
report_densities_extract_new$L1_difference
report_densities_extract_new$L2_difference
report_densities_extract_new$KS_test

# ---- 10. Verdicts and saved outputs ----
# >>> RESULT: what holds and what does not, for the three tables <<<
verdicts

# old vs new
fwrite(structure_summary, paste0(out_prefix, "_structure.csv"))
fwrite(comparison, paste0(out_prefix, "_columns.csv"))
fwrite(changes, paste0(out_prefix, "_changed_cells.csv"))
fwrite(changes_top, paste0(out_prefix, "_changed_values_top.csv"))
fwrite(species_renamed, paste0(out_prefix, "_species_renamed.csv"))
fwrite(rows_only_old, paste0(out_prefix, "_rows_only_old.csv"))
fwrite(rows_only_new, paste0(out_prefix, "_rows_only_new.csv"))

# the three tables
fwrite(verdicts, paste0(out_prefix_all, "_verdicts.csv"))
fwrite(census_inventory_all, paste0(out_prefix_all, "_census_inventory.csv"))
fwrite(rows_per_census, paste0(out_prefix_all, "_rows_per_census.csv"))
fwrite(totals, paste0(out_prefix_all, "_totals.csv"))
fwrite(pairs_summary, paste0(out_prefix_all, "_pairs_summary.csv"))
fwrite(pairs_columns, paste0(out_prefix_all, "_pairs_columns.csv"))
fwrite(taxonomy_check, paste0(out_prefix_all, "_vs_viewtaxonomy.csv"))

cat(
  "Columns changed:", if (length(changed_columns)) paste(changed_columns, collapse = ", ") else "none", "\n",
  "Rows only in old:", length(ids_only_old), "| rows only in new:", length(ids_only_new), "\n",
  "Expectations that hold:", sum(verdicts$holds), "of", nrow(verdicts), "\n",
  file = log_file, append = TRUE
)
cat("✓ Outputs written to", CHECKS_folder, "\n")
