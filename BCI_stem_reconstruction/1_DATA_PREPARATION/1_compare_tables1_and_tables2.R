# ========================================================================
# SCRIPT: 1_compare_tables1_and_tables2.R
# PURPOSE: Reconcile two BCI `ViewFullTable` datasets and produce a single
#          validated, analysis-ready dataset with complete Tag x Census mapping.
#
# Detailed process (Table of Contents):
#   0. Setup (deterministic environment)
#      - Remove objects from previous sessions,
#      - Load required packages (`data.table`, `here`),
#      - Set number of threads for deterministic performance.
#   1. Configuration (environment-specific inputs)
#      - Define site code,
#      - Set paths for two source dataset versions,
#      - Define output/check directories and create them if missing.
#   2. Load input files (exact schema expectation)
#      - Load two versions as tab-separated with explicit colClasses,
#      - Enforce existence of the source files, and standardize missradianing values.
#   3. Quick QA (sanity checks)
#      - Quick `table()` counts on `ListOfTSM`/`Status`,
#      - Ensure `ExactDate` coverage includes 2022 and 2023.
#   4. Compare datasets (2022-2023 subset + data equivalence checks)
#      - Subset both tables to 2022/2023,
#      - Compare each column for identity,
#      - Compute DBH density and KS tests to confirm distributions match.
#   5. Prepare main table (canonicalization and normalization)
#      - Choose authoritative dataset (`ViewFullTable_1` unless differences found),
#      - Preserve original `CensusID` in `CensusID_raw`,
#      - Normalize factor IDs for joins and consistency.
#   6. Fill missing rows (Tag x Census completeness)
#      - Generate full tag-by-census grid,
#      - Add intentionally missing rows (tag/census combos) for run-level completeness.
#   7. Save outputs (persistent results and diagnostics)
#      - Write processed dataset to `DATA/PROCESSED`,
#      - Write QA outputs to `DATA/CHECKS`.
#
# Notes:
#  - Uses `data.table` for speed and memory efficiency.
#  - Input expectations: tab-delimited `ViewFullTable_*` files with columns
#    `ExactDate`, `Tag`, `StemTag`, `TreeID`, `StemID`, `CensusID`, `DBH`.
# ========================================================================

# ---- 0. Setup ----
# Clean workspace to avoid conflicts from previous sessions
rm(list = ls())

# Load required libraries
library(data.table) # fast data manipulation
library(here) # convenient file paths
# References: https://arelbundock.com/posts/dt_tb_df/index.html
#             https://rdatatable.gitlab.io/data.table/

# Set data.table to use multiple threads for performance (adjust as needed)
setDTthreads(15) # tune this value for your machine

# ---- 1. Configuration ----
# User-editable variables and file locations. Ensure the paths below point to
# the expected tab-delimited `ViewFullTable` files for the site.

site <- "bci" # Site code for BCI

# Input folders for the two datasets to compare
INPUT_folder_1 <- here::here("DATA", "RAW", "ViewFiles_bci_allcensuses")
INPUT_folder_2 <- here::here("DATA", "RAW")
# Expected files:
#   - file.path(INPUT_folder_1, paste0("ViewFullTable_", site, ".csv"))
#   - file.path(INPUT_folder_2, "ViewFullTable_census2022-23.txt")

# Output and diagnostics folders
OUTPUT_folder <- here::here("DATA", "PROCESSED")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
}

CHECK_folder <- here::here("DATA", "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
}

# End of configuration (adjust above values as needed for local runs)

# ---- 2. Load and verify input data ----
# Read the two versions of ViewFullTable. We expect tab-delimited tables with
# columns such as `ExactDate`, `Tag`, `StemTag`, `TreeID`, `StemID`, `DBH`.

file1 <- file.path(INPUT_folder_1, paste0("ViewFullTable_", site, ".csv"))
file2 <- file.path(INPUT_folder_2, "ViewFullTable_census2022-23.txt")
stopifnot(file.exists(file1))
stopifnot(file.exists(file2))

ViewFullTable_1 <- fread(
  file = file1,
  sep = "\t",
  colClasses = c(
    PlotName = "factor",
    PlotID = "factor",
    CensusID = "factor",
    QuadratName = "character",
    ExactDate = "Date"
  ),
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NULL", "")
)

# Quick QA: small frequency table to spot weird Status/ListOfTSM combos
table(ViewFullTable_1[, c("ListOfTSM", "Status")], useNA = "ifany")

# Load second version
ViewFullTable_2 <- fread(
  file = file2,
  sep = "\t",
  colClasses = c(
    PlotName = "factor",
    PlotID = "factor",
    CensusID = "factor",
    QuadratName = "character",
    ExactDate = "Date"
  ),
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NULL", "")
)

# ---- 3. Compare datasets (focus on 2022 & 2023) ----

# Quick QA for second file
table(ViewFullTable_2[, c("ListOfTSM", "Status")], useNA = "ifany")

# If files load without errors, move on to subset and compare below.

# ---- 4. Compare datasets (focus on 2022 & 2023) ----
# Subset both datasets to 2022 and 2023 and perform column-wise and
# distributional comparisons (DBH) to detect any substantive differences.
# Summarize available years in each dataset
dates_1 <- as.data.table(table(ViewFullTable_1[!is.na(ExactDate), .(ExactDate)]))[, year := year(ExactDate)]
dates_2 <- as.data.table(table(ViewFullTable_2[!is.na(ExactDate), .(ExactDate)]))[, year := year(ExactDate)]

# Confirm both datasets contain 2022 and 2023
data.table(table(dates_1$year))
data.table(table(dates_2$year))

# Add year column to both datasets
ViewFullTable_1[, year := year(ExactDate)]
ViewFullTable_2[, year := year(ExactDate)]

# Subset to only 2022 and 2023, and order for comparison
df1 <- copy(ViewFullTable_1)[year %in% c(2022, 2023)][order(TreeID, StemID)]
df2 <- copy(ViewFullTable_2)[year %in% c(2022, 2023)][order(TreeID, StemID)]

### Check for differences column-by-column
cols_to_compare <- sort(unique(c(names(df1), names(df2))))

# For each column, align factor levels and compare values
comparison <- rbindlist(lapply(cols_to_compare, function(col) {
  old_col <- addNA(as.factor(df1[[col]]))
  new_col <- addNA(as.factor(df2[[col]]))
  all_levels <- union(levels(old_col), levels(new_col))
  old_f <- factor(old_col, levels = all_levels)
  new_f <- factor(new_col, levels = all_levels)
  same_levels <- identical(levels(old_f), levels(new_f))
  same_values <- identical(as.integer(old_f), as.integer(new_f))
  data.table(
    column = col,
    result = if (same_levels && same_values) "Similar" else "Different"
  )
}))

# Print any columns that differ (should be empty if identical)
comparison[result == "Different"]
# No differences found

# ---- Helper: show_mismatches ----
# Purpose: Return rows where df1 and df2 differ for a given column.
# Args: col (string) - column name to compare
# Returns: data.table with mismatch row indices and differing values, or a
#          short message if no differences are found.
show_mismatches <- function(col) {
  old_col <- addNA(as.factor(df1[[col]]))
  new_col <- addNA(as.factor(df2[[col]]))
  all_levels <- union(levels(old_col), levels(new_col))
  old_f <- factor(old_col, levels = all_levels)
  new_f <- factor(new_col, levels = all_levels)
  mismatch_idx <- which(as.integer(old_f) != as.integer(new_f))
  if (length(mismatch_idx) == 0) {
    return(data.table(
      column = col,
      message = "No mismatches — columns are identical"
    ))
  }
  data.table(
    row = mismatch_idx,
    old_value = as.character(old_f)[mismatch_idx],
    new_value = as.character(new_f)[mismatch_idx]
  )
}

# Run mismatch check for all columns (should be empty if identical)
mismatches_list <- lapply(cols_to_compare, show_mismatches)
mismatches_list

# Final check: compare DBH distributions using density
d1 <- density(df1$DBH, na.rm = TRUE)
d2 <- density(df2$DBH, na.rm = TRUE)

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

# Run density comparison and print summary statistics
report_densities <- compare_densities(df1$DBH, df2$DBH)
report_densities$L1_difference
report_densities$L2_difference
report_densities$KS_test

# NOTE: Both datasets for years 2022 and 2023 are identical
# Proceed with the main dataset only

# ---- 5. Prepare main table and normalize IDs ----
# Having compared both versions and found no substantial differences for the
# years of interest, proceed using `ViewFullTable_1` as the canonical table.
ViewFullTable <- copy(ViewFullTable_1)[, year := NULL]

# Rename CensusID to CensusID_raw to preserve original labels, then assign
# standardized sequential CensusID values below.
setnames(ViewFullTable, "CensusID", "CensusID_raw")
# Ensure important identifier columns are factors for efficient joins and memory
ViewFullTable[, Tag := as.factor(Tag)]
ViewFullTable[, StemTag := as.factor(StemTag)]
ViewFullTable[, TreeID := as.factor(TreeID)]
ViewFullTable[, StemID := as.factor(StemID)]

# Summarize census dates for each CensusID_raw
censusid_dates <- ViewFullTable[, .(
  n_dates = uniqueN(ExactDate),
  min_date = min(ExactDate, na.rm = TRUE),
  mean_date = mean(ExactDate, na.rm = TRUE),
  max_date = max(ExactDate, na.rm = TRUE)
), by = CensusID_raw][order(mean_date)]

# Assign sequential CensusID values
censusid_dates[, CensusID := seq_len(.N)]

# Set keys for efficient joining
setkey(ViewFullTable, CensusID_raw)
setkey(censusid_dates, CensusID_raw)

# Join new CensusID into main table
ViewFullTable <-
  censusid_dates[, .(CensusID_raw, CensusID)][
    ViewFullTable,
    on = "CensusID_raw"
  ]

# Remove old CensusID_raw column (now replaced)
ViewFullTable[, CensusID_raw := NULL]

# ---- 6. Fill missing rows (Tag x CensusID completeness) ----
# Goal: create a complete grid of all Tag x CensusID combinations and add rows
# with NA measurements where these combinations are missing (so downstream
# longitudinal analyses have a full panel for every Tag).

sort(unique(round(ViewFullTable[!is.na(DBH)]$DBH, 1)))

# NOTE: 12 INDIVIDUALS WITH DBH < 10 CM — LIKELY TAPERING ERRORS
# NOTE: no individual with zero DBH
table(ViewFullTable[DBH < 10]$DBH)
range(ViewFullTable$DBH, na.rm = TRUE)

# sort print(inspectdf::inspect_na(ViewFullTable), n = 50) by tag and censusID
## NOTE: Tag and CensusID have no NAs, so we can use them to create a complete grid

tag_census_unique <- unique(copy(ViewFullTable[, .(Tag, CensusID)]))

# Create complete grid of all Tag-CensusID combinations
tag_ranges <- tag_census_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]

## Get all combinations tags and census IDs
## NOTE: LINE BELOW TAKES 34 MINUTES TO RUN (ON MY COMPUTER) FULL DATASET
complete_grid <- tag_ranges[, .(CensusID = seq.int(min_c, max_c)), by = Tag]

## Export the complete grid for future use (so we don't have to re-run above line)
saveRDS(complete_grid, file = here::here("DATA", "CHECKS", "complete_tag_census_grid.rds"))

# complete_grid <- readRDS(file = here::here("DATA", "CHECKS", "complete_tag_census_grid.rds"))

# Find which combinations are missing from the original data
missing_combinations <- complete_grid[!ViewFullTable,
  on = .(Tag, CensusID)
]

unique(ViewFullTable[, .(Tag, TreeID)])

# NOTE: 141 MISSING TAG-CENSUSID COMBINATIONS TO ADD
nrow(missing_combinations)

### LETS GET THE METADATA FOR THE MISSING INFORMATION
# Get metadata columns (excluding measurements)

diff_names <- setdiff(names(ViewFullTable), names(missing_combinations))

for (col in diff_names) {
  col_class <- class(ViewFullTable[[col]])
  if (any(col_class %in% c("factor", "character"))) {
    missing_combinations[, (col) := as.character(NA)]
  } else if (any(col_class %in% c("integer", "numeric"))) {
    missing_combinations[, (col) := as.numeric(NA)]
  } else if (inherits(ViewFullTable[[col]], "Date")) {
    missing_combinations[, (col) := as.Date(NA)]
  } else {
    missing_combinations[, (col) := NA]
  }
}

## sort columns to match ViewFullTable
setcolorder(missing_combinations, names(ViewFullTable))

# Scalar values (same for all rows)
missing_combinations[, PlotName := unique(ViewFullTable$PlotName)]
missing_combinations[, PlotID := unique(ViewFullTable$PlotID)]

# Tag-level attributes (one join instead of 13 separate joins)
tag_cols <- c(
  "Family", "Genus", "SpeciesName", "Mnemonic", "Subspecies",
  "SpeciesID", "SubspeciesID", "QuadratName", "QuadratID",
  "PX", "PY", "QX", "QY", "TreeID"
)

tag_lookup <- unique(ViewFullTable[, c("Tag", tag_cols), with = FALSE], by = "Tag")

missing_combinations[tag_lookup,
  (tag_cols) := mget(paste0("i.", tag_cols)),
  on = "Tag"
]

# CensusID-level attributes
missing_combinations[
  unique(ViewFullTable[, .(CensusID, PlotCensusNumber)]),
  PlotCensusNumber := i.PlotCensusNumber,
  on = "CensusID"
]

unique(missing_combinations$DBH)

# QuadratID + QuadratName + CensusID level attributes (one join instead of two)
census_quadrat_lookup <- unique(ViewFullTable[, .(QuadratID, QuadratName, CensusID, ExactDate, Date)])

missing_combinations[census_quadrat_lookup,
  `:=`(ExactDate = i.ExactDate, Date = i.Date),
  on = .(QuadratID, QuadratName, CensusID)
]

# Verify columns propagated
columns_to_propagate <- c(
  "PlotName", "PlotID", "Family", "Genus", "SpeciesName", "Mnemonic",
  "Subspecies", "SpeciesID", "SubspeciesID", "QuadratName", "QuadratID",
  "PX", "PY", "QX", "QY", "TreeID", "PlotCensusNumber", "ExactDate", "Date"
)

# Check for any NAs in propagated columns
sapply(missing_combinations[, columns_to_propagate, with = FALSE], function(x) sum(is.na(x)))

# print(inspectdf::inspect_na(missing_combinations), n = 50)

unique(names(ViewFullTable) == names(missing_combinations))

# Add the missing rows to the main dataset
ViewFullTable_no_missing_tags_census <- rbindlist(
  list(ViewFullTable[, missing_tag_census := FALSE], missing_combinations[, missing_tag_census := TRUE]),
  fill = TRUE,
  use.names = TRUE
)

# print(inspectdf::inspect_na(ViewFullTable_no_missing_tags_census), n = 50)

(nrow(ViewFullTable) + nrow(missing_combinations)) == nrow(ViewFullTable_no_missing_tags_census)

# # Sort by Tag and CensusID for cleaner viewing
setorder(ViewFullTable_no_missing_tags_census, Tag, CensusID)

## CHECK AGAIN IF MISSING TAG CENSUS COMBINATIONS

missing_combinations_final <- complete_grid[!ViewFullTable_no_missing_tags_census,
  on = .(Tag, CensusID)
]
# good

# ---- 7. Save outputs and diagnostics ----
# Save final validated table for downstream analysis and future reproducibility
saveRDS(ViewFullTable_no_missing_tags_census, here::here("DATA", "PROCESSED", "1_ViewFullTable_no_missing_tags_census.rds"))

## get mean of minimum census per tag (example)
mean(ViewFullTable_no_missing_tags_census[, .(min_c = min(CensusID)), by = Tag]$min_c)
## get mean of maximum census per tag (example)
mean(ViewFullTable_no_missing_tags_census[, .(max_c = max(CensusID)), by = Tag]$max_c)
## get the median number of censuses per tag
median(ViewFullTable_no_missing_tags_census[, .N, by = Tag]$N)
## get the mode number of censuses per tag
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(ViewFullTable_no_missing_tags_census[, .N, by = Tag]$N)

## why 9 repeats more?
tags_per_census <- unique(ViewFullTable_no_missing_tags_census[, .(Tag, CensusID)])

pdf("./DATA/CHECKS/ntags_per_censusID_plot.pdf",
  width = 8,
  height = 4
)
with(tags_per_census[, .N, by = CensusID][order(CensusID)], plot(CensusID, N,
  type = "b",
  xlab = "CensusID", ylab = "Number of unique Tags",
  main = "Number of unique Tags per CensusID"
))
dev.off()
