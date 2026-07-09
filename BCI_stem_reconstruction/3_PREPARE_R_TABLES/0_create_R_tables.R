################################################################################
# FORESTGEO STEM TABLE CREATION SCRIPT
################################################################################
# Purpose: Build ForestGEO Rtables and a species table from a reconstructed
# census dataset plus taxonomy input.
#
# Inputs:
#   - DATA/
#
# Outputs:
#   - DATA/RTABLES/[site].stem[n].Rdata
#   - DATA/RTABLES/[site].spptable.rdata
#   - DATA/CHECKS/ diagnostic CSV files
#
# Main pipeline:
#    1. Setup and load input data                          (Sections 1–2)
#    2. Standardize stems across all censuses                (Section 3)
#    3. Map raw Status values to canonical codes A/D/G/P/N  (Section 4)
#    4. Validate encounter histories before propagation      (Section 5)
#    5. Propagate no-data states into complete histories     (Section 6)
#    6. Correct PD/PG and resurrection anomalies            (Sections 7–8)
#    7. Safety-net: fix D/G stranded between two A’s        (Section 9)
#    8. Derive tree-level histories; apply per-census D←4G remap (Section 10)
#    9. Assess biology across all history versions          (Section 11)
#   10. Data-quality diagnostics and status × DBH audit    (Sections 12–13)
#   11. Export per-census R tables and species table        (Sections 14–15)
################################################################################
# Variable summary:
#   Status               raw field record, later mapped to A/D/G/P/N
#   original_status      stem encounter history before propagation
#   new_status           propagated stem history before D/G remap
#   corrected_new_status final stem history after tree-level D/G adjustment
#   tree_histories       tree-level history aggregated from stem histories
#   DBHs                 numeric matrix of stem DBH measurements by census
#   Rstatus              exported per-census corrected stem status
#   DFstatus             legacy ForestGEO status field for prior stems
################################################################################

# ========================================================================
# SECTION 1: SETUP AND INITIALIZATION
# ========================================================================
# Load libraries, configure paths, and define output column mappings.
# This section also ensures required folders exist before the pipeline runs.

# Remove all objects from workspace to start fresh
rm(list = ls())

# Load libraries ####
# Load data.table for fast data manipulation and better console output
library(data.table)
# References 1: https://arelbundock.com/posts/dt_tb_df/index.html
# References 2: https://rdatatable.gitlab.io/data.table/
# library(ggplot2)

getDTthreads() # Show number of threads data.table will use
# data.table is using 8 threads by default
# To change, use setDTthreads(<n>) where <n> is number of threads
# setDTthreads(1)

# ========================================================================
# USER CONFIGURATION - MODIFY THESE PARAMETERS
# ========================================================================
# ⚠️ REQUIRED: Set these variables before running the script

# 1. Set working directory path
main_path <- getwd()
cat("Working directory:", main_path, "\n")

# 2. Site identifier (must match file naming convention)
site <- "Huai Kha Khaeng" # ⚠️ CHANGE THIS for your site (e.g., "bci", "scbi", "hkk")
cat("Processing site:", site, "\n")

# 3. Input folder containing ViewFullTable and ViewTaxonomy CSV files
INPUT_folder <- file.path(main_path, "DATA")
# ℹ️ Files required:
#    - ViewFullTable_[site].csv (tab-delimited)
#    - ViewTaxonomy_[site].csv (tab-delimited)

# 4. Output folder for final .Rdata files
OUTPUT_folder <- file.path(main_path, "RTABLES")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
  cat("✓ Created OUTPUT folder:", OUTPUT_folder, "\n")
}

# 5. Diagnostics folder for QA/QC reports
CHECK_folder <- file.path(main_path, "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
  cat("✓ Created CHECKS folder:", CHECK_folder, "\n")
}

# 6. Export log messages to an external file showing the mistakes
log_file <- file.path(CHECK_folder, paste0("check_", site, ".txt"))
# create log file and populate with messsages
create_log_file <- function(log_file) {
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  file.create(log_file)
}

create_log_file(log_file)

# add message about ViewFullTable[Tag %in% "5319323"] to log file
print_to_log <- function(message, log_file, new_message = TRUE) {
  if (new_message) {
    cat(strrep("-", 60), file = log_file, append = TRUE, sep = "\n")
  }
  cat(message, file = log_file, append = TRUE, sep = "\n")
}

# ========================================================================

# ── Output column definitions ────────────────────────────────────────────────
# Two parallel vectors (same length, same order) that define:
#   ViewFullTable_columns_to_keep — original ForestGEO DB column names to keep
#   new_names_columns_to_keep     — standardized ForestGEO R-table names to rename them to
# Both are used later when subsetting and renaming columns for each census export.
#
# Column notes:
#   "new_status"                    — computed by this script; does not exist in raw data
#   "obs_row_id"                    — row identifier from the DP reconstruction pipeline
#   "DP_PosteriorReconstructedProb" — DP model posterior probability for ReconstructedStemID
#   "sp"                            — species mnemonic (short ForestGEO species code)
#   "gx", "gy"                      — stem X/Y coordinates in meters within the plot
#   "hom"                           — height of measurement (cm above ground where DBH taken)
#   "DFstatus"                      — legacy ForestGEO column name for the raw field status
#   "Rstatus"                       — script-computed corrected status (from "new_status")

# Columns to retain from ViewFullTable (original DB column names)
ViewFullTable_columns_to_keep <- c(
  "TreeID", "StemID", "Tag", "StemTag", "Mnemonic", # stem / tree identifiers
  "QuadratName", "PX", "PY", # spatial location in plot (meters)
  "DBHID", "CensusID", "DBH", "HOM", # measurement record
  "ExactDate", "Status", "ListOfTSM", # raw field status + measurement codes
  "Date", "new_status" # date + script-computed corrected status
)

# Standardized ForestGEO R-table column names (must match ViewFullTable_columns_to_keep
# one-to-one: same length and same order)
new_names_columns_to_keep <- c(
  "treeID", "stemID", "tag", "StemTag", "sp", # sp = species mnemonic
  "quadrat", "gx", "gy", # gx/gy = plot coordinates in meters
  "MeasureID", "CensusID", "dbh", "hom", # hom = height of measurement
  "ExactDate", "DFstatus", "codes", # DFstatus = raw field status (legacy name)
  "date", "Rstatus" # date + script-computed corrected status
)

# ---- Hard guard: the two column vectors must be 1-to-1 -------------------
# A typo or accidental edit to either vector would silently mis-rename
# columns at export time. Catch it now, before the pipeline starts.
stopifnot(
  length(ViewFullTable_columns_to_keep) == length(new_names_columns_to_keep),
  !anyDuplicated(ViewFullTable_columns_to_keep),
  !anyDuplicated(new_names_columns_to_keep)
)

# ========================================================================
# SECTION 2: LOAD INPUT DATA
# ========================================================================
# Load the reconstructed census table and taxonomy table, coerce key types,
# and write diagnostic reports for raw value distributions.

################################################################################
# HELPER FUNCTION: show_levels
################################################################################
# PURPOSE:
#   Display unique values in categorical columns for data validation
#
# PARAMETERS:
#   df         - Data frame or data.table to examine
#   n_to_print - Maximum unique values to display (default: 10; Inf = all)
#   output     - Output format: "print" (console) or "df" (data frame)
#
# RETURNS:
#   - If output="print": Prints formatted text showing unique values
#   - If output="df": Returns data frame with levels as rows
#
# EXAMPLE:
#   show_levels(ViewFullTable, n_to_print = Inf, output = "print")
################################################################################
show_levels <- function(df, n_to_print = 10, output = "print") {
  cat_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
  if (length(cat_cols) == 0) {
    if (output == "print") {
      cat("No character or factor columns found.\n")
      return(invisible(NULL))
    } else {
      return(data.frame())
    }
  }
  uniques_list <- lapply(cat_cols, function(col) {
    vals <- if (is.factor(df[[col]])) {
      levels(df[[col]])
    } else {
      unique(df[[col]])
    }
    vals <- vals[!is.na(vals)]
    if (!is.infinite(n_to_print)) {
      vals <- head(vals, n_to_print)
    }
    vals
  })
  names(uniques_list) <- cat_cols
  if (output == "print") {
    for (col in cat_cols) {
      vals <- uniques_list[[col]]
      total_vals <- length(if (is.factor(df[[col]])) levels(df[[col]]) else unique(df[[col]])) - sum(is.na(df[[col]]))
      cat(paste0(
        "Levels of ", col,
        " (showing ",
        ifelse(is.infinite(n_to_print), "all", paste0(length(vals), " of ", total_vals)),
        "): ",
        paste(vals, collapse = ", "),
        if (!is.infinite(n_to_print) && total_vals > n_to_print) " ..." else "",
        "\n\n\n"
      ))
    }
    return(invisible(NULL))
  } else if (output == "df") {
    max_len <- max(sapply(uniques_list, length))
    padded_list <- lapply(uniques_list, function(vals) {
      c(vals, rep(NA, max_len - length(vals)))
    })
    return(as.data.frame(padded_list, stringsAsFactors = FALSE))
  } else {
    stop("Invalid output option. Use 'print' or 'df'.")
  }
}

# ========================================================================
# LOAD RAW DATA FROM CSV FILES
# ========================================================================

cat("\n", strrep("=", 70), "\n")
cat("LOADING DATA FOR SITE:", toupper(site), "\n")
cat(strrep("=", 70), "\n\n")

# Load ViewFullTable_RAW (census data) ####
cat("📂 Loading ViewFullTable_RAW...\n")

ViewFullTable_RAW <- as.data.table(fread(file.path(
  INPUT_folder,
  "ViewFullTable_hkk.csv"
)))

# make all NULL to NA
ViewFullTable_RAW[ViewFullTable_RAW == "NULL"] <- NA

ViewFullTable_RAW[, TreeID := paste(TreeID, Tag, sep = "_")]

# ── Standardize column types ──────────────────────────────────────────────────
# The RDS file may carry incorrect types from the database export pipeline.
# Each conversion is explained below:
#   PlotName, PlotID         → factor    : categorical identifiers, not quantities
#   CensusID                 → numeric   : stored as factor in DB; must go through
#                                          as.character() first to avoid silent factor-
#                                          to-integer coercion (which would return
#                                          level indices, not the actual census numbers)
#   ExactDate                → Date      : DB stores as character string "YYYY-MM-DD"
#   Tag, TreeID, StemID, StemTag,
#   QuadratName              → character : codes like "0001" must stay as strings;
#                                          integer coercion silently drops leading zeros
ViewFullTable_RAW <- transform(ViewFullTable_RAW,
  PlotName = as.factor(PlotName),
  PlotID = as.factor(PlotID),
  CensusID = as.numeric(as.character(CensusID)), # factor → character → numeric
  ExactDate = as.Date(ExactDate),
  Tag = as.character(Tag),
  TreeID = as.character(TreeID),
  StemID = as.character(StemID),
  StemTag = as.character(StemTag),
  QuadratName = as.character(QuadratName) # preserve leading zeros
)

# 1. Create a copy and apply the function
ViewFullTable <- copy(ViewFullTable_RAW)

## DIAGNOSTIC: Verify Tag <-> TreeID 1:1 mapping
# Check 1 to 1 match tag and treeid
setkey(ViewFullTable, TreeID)
setkey(ViewFullTable, Tag)
ViewFullTable[, .(n_tags = uniqueN(Tag)), by = TreeID][order(-n_tags)]
ViewFullTable[, .(n_treeid = uniqueN(TreeID)), by = Tag][order(-n_treeid)]

# 1 tag with 2 treeids
ViewFullTable[Tag %in% "5319323", .(CensusID, Mnemonic, TreeID, Tag, StemID, StemTag, ListOfTSM, Status, DBH, ExactDate)]

print_to_log("Tag 5319323 has 2 TreeIDs:", log_file, new_message = TRUE)
print_to_log(
  capture.output(
    print(ViewFullTable[
      Tag %in% "5319323",
      .(CensusID, Mnemonic, TreeID, Tag, StemID, StemTag, ListOfTSM, Status, DBH, ExactDate)
    ])
  ),
  log_file,
  new_message = FALSE
)

ViewFullTable[is.na(Tag)]
ViewFullTable[is.na(TreeID)]

# ── DIAGNOSTIC: TreeID × CensusID completeness ──────────────────────────────────
# Purpose: identify TreeID that are missing one or more intermediate census records.
# "Complete" means every census between the TreeID's first and last observed census
#   is present. A gap means a row is missing from the raw ViewFullTable for that
#   TreeID × census combination.
#   only MISSING INTERMEDIATE censuses count (e.g., present in census 2 and 4
#   but absent in census 3 would be a gap of 1).
xraw_unique <- unique(ViewFullTable[, .(TreeID, CensusID)])
# Get the range per TreeID
TreeID_ranges <- xraw_unique[, .(min_c = min(as.character(as.numeric(CensusID))), max_c = max(as.character(as.numeric(CensusID)))), by = TreeID]
# Add expected count (how many censuses should exist)
TreeID_ranges[, expected_count := as.numeric(max_c) - as.numeric(min_c) + 1L]
# Get actual count per TreeID
actual_counts <- xraw_unique[, .(actual_count = .N), by = TreeID]
# Merge and compare
TreeID_check <- TreeID_ranges[actual_counts, on = "TreeID"]
TreeID_check[, complete := actual_count == expected_count]
# Summary
TreeID_check[, .N, by = complete]
# See TreeIDs with missing censuses
missing_TreeIDs <- TreeID_check[complete == FALSE]
missing_TreeIDs[, gap := expected_count - actual_count]
# Summary stats
cat("Total TreeIDs:", nrow(TreeID_check), "\n")
cat("Complete TreeIDs:", TreeID_check[complete == TRUE, .N], "\n")
cat("TreeIDs with gaps:", TreeID_check[complete == FALSE, .N], "\n")
if (nrow(missing_TreeIDs) > 0) {
  cat("Total missing observations:", sum(missing_TreeIDs$gap), "\n")
}

# ── DIAGNOSTIC: StemID × CensusID completeness ───────────────────────────────
# Purpose: same completeness check as above, but at the reconstructed stem level.
# Each unique StemID ("TreeID_ReconstructedStemID") should appear in every census
#   between its first and last observed census. Gaps here would indicate that
#   the StemID overwrite collapsed or lost rows for some stems.
xraw_unique <- unique(ViewFullTable[, .(StemID, CensusID)])
# Get the range per StemID
StemID_ranges <- xraw_unique[, .(min_c = min(as.character(as.numeric(CensusID))), max_c = max(as.character(as.numeric(CensusID)))), by = StemID]
# Add expected count (how many censuses should exist)
StemID_ranges[, expected_count := as.numeric(max_c) - as.numeric(min_c) + 1L]
# Get actual count per StemID
actual_counts <- xraw_unique[, .(actual_count = .N), by = StemID]
# Merge and compare
StemID_check <- StemID_ranges[actual_counts, on = "StemID"]
StemID_check[, complete := actual_count == expected_count]
# Summary
StemID_check[, .N, by = complete]
# See StemIDs with missing censuses
missing_StemIDs <- StemID_check[complete == FALSE]
missing_StemIDs[, gap := expected_count - actual_count]
# Summary stats
cat("Total StemIDs:", nrow(StemID_check), "\n")
cat("Complete StemIDs:", StemID_check[complete == TRUE, .N], "\n")
cat("StemIDs with gaps:", StemID_check[complete == FALSE, .N], "\n")
if (nrow(missing_StemIDs) > 0) {
  cat("Total missing observations:", sum(missing_StemIDs$gap), "\n")
}

# --------------------------------------------------------------------
# Verify data quality
# --------------------------------------------------------------------
cat("Checking Status vs ListOfTSM combinations:\n")
table(ViewFullTable[, c("ListOfTSM", "Status")], useNA = "ifany")

# Export unique values for QA/QC
# show_levels(ViewFullTable, n_to_print = Inf, output = "print")
fwrite(show_levels(ViewFullTable, n_to_print = Inf, output = "df"),
  file = file.path(CHECK_folder, paste0("ViewFullTable_levels_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

# ========================================================================
# SECTION 3: STANDARDIZE DATA FORMAT ACROSS CENSUSES
# ========================================================================
# Build a master stem list and force every census table to include the same
# StemIDs in the same row order. Missing stems are retained as NA rows so that
# encounter histories can be assembled consistently across all censuses.

# --------------------------------------------------------------------
# Create master stem list with fixed attributes
# --------------------------------------------------------------------

# Apply this immediately after the StemID overwrite block
# ViewFullTable <- ViewFullTable[!grepl("_NA$", StemID)]

# fixed_columns: the minimum identifier set that does not change across
# censuses. Reconstruction is keyed on ReconstructedStemID (legacy database
# StemIDs were rewritten in step 2); StemID below is the composite
# TreeID_ReconstructedStemID created earlier.
fixed_columns <- c(
  "PlotName", "PlotID",
  "Mnemonic",
  "QuadratName", "QuadratID", "PX", "PY", "QX", "QY",
  "TreeID", "Tag",
  "StemID", "StemNumber", "StemTag"
)

# Convert to data.table if not already (defensive coding)
setDT(ViewFullTable)

# Create master stem list: One row per stem, fixed attributes only
# unique() removes duplicate rows (same stem in multiple censuses)
unique_StemID <- unique(ViewFullTable[, ..fixed_columns])

## return the duplicated rows in unique_StemID using StemID as the identifier
duplicated_stems <- unique_StemID[duplicated(unique_StemID$StemID) | duplicated(unique_StemID$StemID, fromLast = TRUE)]

# VERIFY: If duplicates exist, there's a data quality issue (same stem with
# different fixed attributes)
if (anyDuplicated(unique_StemID$StemID)) {
  warning("Duplicate StemIDs found in unique_StemID! Check your data.")
} else {
  cat("✓ All StemIDs are unique\n")
}

# Show all unique values
# show_levels(unique_StemID, n_to_print = Inf, output = "print")
fwrite(show_levels(unique_StemID, n_to_print = Inf, output = "df"),
  file = file.path(CHECK_folder, paste0("unique_StemID_levels_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

## Check date information per census ####
# VERIFY: Census 1 should contain the earliest date
ViewFullTable[, .(
  n_dates = uniqueN(ExactDate),
  min_date = min(ExactDate, na.rm = TRUE),
  max_date = max(ExactDate, na.rm = TRUE)
), by = CensusID][order(CensusID)]

## Split ViewFullTable into separate census tables ####
# Creates a list where each element is one census
# Split by both PlotID and CensusID (in case multiple plots in dataset)
# drop=TRUE removes empty combinations
ViewFullTable_split_unbalanced <- split(ViewFullTable, by = c("PlotID", "CensusID"), drop = TRUE)
cat("Split data into", length(ViewFullTable_split_unbalanced), "Plot-Census groups.\n")

current_names <- names(ViewFullTable_split_unbalanced)
ord <- order(as.numeric(sub(".*\\.", "", current_names)))
cat("Current order:", paste(current_names, collapse = ", "), "\n")
cat("Order indices:", paste(ord, collapse = ", "), "\n")
cat("New order:", paste(current_names[ord], collapse = ", "), "\n")

ViewFullTable_split_unbalanced <- ViewFullTable_split_unbalanced[ord]
cat("Reordered by CensusID (after decimal).\n")

## Make all censuses the same format (one row per stemID, ordered in the same way across censuses) ####
# BEFORE: Each census has different stems (recruits, deaths cause row count differences)
# AFTER: All censuses have ALL stems in the SAME order (missing stems filled with NA)

# Show current dimensions - should be DIFFERENT across censuses
cat("\n📊 Current census dimensions (BEFORE standardization):\n")
lapply(ViewFullTable_split_unbalanced, dim)

print_to_log("Current census dimensions (nrows + columns) BEFORE standardization:", log_file, new_message = TRUE)
print_to_log(
  capture.output(
    lapply(ViewFullTable_split_unbalanced, dim)
  ),
  log_file,
  new_message = FALSE
)

## Standardize each census to have ALL stems in the SAME order ####
# This is the CORE operation of Section 3
# For each census:
#   1. Match its stems to the master list (unique_StemID)
#   2. Reorder rows to match master list order
#   3. Fill missing stems with NA rows
#   4. Fill in fixed attributes from master list
#   5. Fill in census identifiers (CensusID, PlotCensusNumber)
ViewFullTable_split <- lapply(ViewFullTable_split_unbalanced, function(X) {
  # Ensure it's a data.table
  if (!is.data.table(X)) setDT(X)
  # SAFETY CHECK: match() silently uses the first row when a StemID appears more
  # than once in X. Detect this before it happens so the root cause can be fixed.
  dup_ids <- X[duplicated(StemID), unique(StemID)]
  if (length(dup_ids) > 0L) {
    census_id <- unique(na.omit(X$CensusID))
    warning(sprintf(
      "Census %s: %d StemID(s) appear more than once — match() will silently keep only the first row. Duplicates: %s",
      paste(census_id, collapse = "/"),
      length(dup_ids),
      paste(head(dup_ids, 10), collapse = ", ")
    ))
  }
  print_to_log(
    sprintf(
      "Census %s: %d StemID(s) appear more than once — match() will silently keep only the first row. Duplicates: %s",
      paste(unique(na.omit(X$CensusID)), collapse = "/"),
      length(dup_ids),
      paste(head(dup_ids, 10), collapse = ", ")
    ),
    log_file,
    new_message = TRUE
  )
  # REORDER + FILL MISSING: match() returns indices or NA for missing stems
  # This aligns X's stems to match unique_StemID's order
  # NA indices create new rows filled with NA
  # NOTE: This relies on StemID being unique in unique_StemID (verified earlier)
  idx <- match(unique_StemID$StemID, X$StemID)
  X <- X[idx]
  # FILL FIXED ATTRIBUTES: Copy all fixed columns from master list
  # This ensures taxonomy, coordinates, etc. are consistent across censuses
  X[, (fixed_columns) := unique_StemID]
  # FILL CENSUS IDENTIFIERS: Propagate census info to all rows
  # unique(na.omit()) extracts the single non-NA value for this census
  if (length(unique(na.omit(X$CensusID))) > 1) {
    bad_vals <- unique(na.omit(X$CensusID))
    cat("Data quality issue: Multiple CensusID values in census ",
      unique(na.omit((X$PlotCensusNumber))), " (", paste(bad_vals, collapse = ", "), ")\n",
      sep = ""
    )
  }
  X[, CensusID := unique(na.omit(CensusID))]
  X[, PlotCensusNumber := unique(na.omit(PlotCensusNumber))]
  return(X)
})

print_to_log("First set of repeated StemID per census:", log_file, new_message = TRUE)
ViewFullTable_split_unbalanced[[2]][StemID %in% c("80482", "80483", "80522"), .(
  StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
)]
ViewFullTable_split[[2]][StemID %in% c("80482", "80483", "80522"), .(
  StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
)]

print_to_log(
  capture.output(
    ViewFullTable_split_unbalanced[[2]][StemID %in% c("80482", "80483", "80522"), .(
      StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
    )]
  ),
  log_file,
  new_message = FALSE
)

print_to_log("Second set of repeated StemID per census:", log_file, new_message = TRUE)
ViewFullTable_split_unbalanced[[6]][StemID %in% c("395212", "387760"), .(
  StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
)]
ViewFullTable_split[[6]][StemID %in% c("395212", "387760"), .(
  StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
)]

print_to_log(
  capture.output(
    ViewFullTable_split_unbalanced[[6]][StemID %in% c("395212", "387760"), .(
      StemID, CensusID, PlotCensusNumber, Tag, TreeID, Mnemonic, QuadratName, PX, PY, DBH, HOM, ExactDate, ListOfTSM, Status
    )]
  ),
  log_file,
  new_message = FALSE
)

# NOTE: match() keeps the first observation per StemID. When duplicate rows
# contain identical information this is the desired behaviour.

## VERIFICATION: Check that standardization worked ####
# All censuses should now have IDENTICAL dimensions (same # rows, same # columns)
cat("\n📊 Census dimensions (AFTER standardization):\n")
lapply(ViewFullTable_split, dim)

print_to_log("Census dimensions (nrows + columns) AFTER standardization:", log_file, new_message = TRUE)
print_to_log(
  capture.output(
    lapply(ViewFullTable_split, dim)
  ),
  log_file,
  new_message = FALSE
)

# Show first few rows of each census for visual inspection
cat("\n📋 First few rows of each census:\n")
lapply(ViewFullTable_split, head, 4)

## VERIFY: check wether length StemID is same as nrow of each census and print that it does
for (i in seq_along(ViewFullTable_split)) {
  n_stems <- nrow(unique_StemID)
  n_rows <- nrow(ViewFullTable_split[[i]])
  if (n_stems != n_rows) {
    warning(sprintf("Census %d: Number of stems (%d) does not match number of rows (%d)!", i, n_stems, n_rows))
  } else {
    cat(sprintf("✓ Census %d: Number of stems matches number of rows (%d)\n", i, n_stems))
  }
}

# Show first few rows of each census for visual inspection
cat("\n📋 First few rows of each census:\n")
lapply(ViewFullTable_split, head, 4)

# Explore NAs
lapply(ViewFullTable_split, inspectdf::inspect_na)

## VERIFY: check fixed columns are IDENTICAL across all censuses ####
# The fixed attributes should be the same for each stem in every census
# If not, there's a data quality issue
# Compare all pairs of censuses: Are fixed columns identical?
# combn() generates all pairs, then checks if fixed columns match
all_equal <- all(
  combn(seq_along(ViewFullTable_split), 2, simplify = TRUE, FUN = function(i) {
    identical(
      ViewFullTable_split[[i[1]]][, ..fixed_columns],
      ViewFullTable_split[[i[2]]][, ..fixed_columns]
    )
  })
)

# Report results
if (all_equal) {
  cat("✓ All fixed columns are IDENTICAL across censuses (standardization successful!)\n\n")
} else {
  cat("⚠ WARNING: Fixed columns differ between censuses. Investigating...\n\n")
}

# If differences found, show details
diffs <- combn(ViewFullTable_split, 2, simplify = FALSE, FUN = function(pair) {
  a <- pair[[1]][, ..fixed_columns]
  b <- pair[[2]][, ..fixed_columns]
  if (!identical(a, b)) {
    list(
      census1 = names(pair[1]),
      census2 = names(pair[2]),
      diff = a[apply(a != b, 1, any), , drop = FALSE] # Rows that differ
    )
  } else {
    NULL
  }
})

# Remove NULL entries (pairs with no differences)
diffs <- Filter(Negate(is.null), diffs)

if (length(diffs) > 0) {
  cat("Differences between censuses:\n")
  print(diffs)
} else {
  cat("No differences found in fixed columns\n")
}

# ========================================================================
# SECTION 4: STATUS CODE TRANSFORMATION
# ========================================================================
# Convert raw Status terms to the canonical codes A/D/G/P/N, resolve "broken
# below" using DBH, and assemble per-stem encounter history strings for later
# propagation and correction.

# --------------------------------------------------------------------
# Transform English status codes to single letters
# --------------------------------------------------------------------
# Gather Status AND DBH from all censuses into one long data.table.
# DBH is needed because the "broken below" status is interpreted in a
# DBH-aware way below: a broken-below stem with a measured DBH is
# treated as alive (A); without DBH it is treated as dead (D). See the
# "broken below" rule further down for the full rationale.
DT_Status <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(StemID, Status, DBH, census = i)]
}), use.names = TRUE, fill = TRUE)

check <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]]
}), use.names = TRUE, fill = TRUE)

if (identical(seq_along(ViewFullTable_split), unique(DT_Status$census))) {
  cat("✓ Census numbering in DT_Status is correct\n")
} else {
  warning("⚠ WARNING: Census numbering in DT_Status is incorrect!\n")
}

## Replace english words by corresponding codes ####
# ------------------------------------------------------------------------
# Resolve "broken below" using DBH BEFORE the wide pivot
# ------------------------------------------------------------------------
# RULE (site-specific field code — applied before the wide pivot):
#   "broken below" + DBH non-NA  → "alive"  (the stem WAS measured this
#                                  census, so the field crew recorded a
#                                  diameter; treat as alive even though
#                                  the field code says "broken below")
#   "broken below" + DBH NA      → "dead"   (no measurement recorded;
#                                  the stem is treated as dead)
#
# WHY HERE (and not later as a string gsub on the encounter histories):
#   The disambiguation needs the per-census DBH value.  Doing it on the
#   long DT_Status table keeps the rule local to a single observation
#   (StemID × census) and avoids any matrix indexing later.  After this
#   step the Status field no longer contains "broken below".
#
# Rule: "dead" / "stem dead" are taken as terminal here (mapped to "D" below).
# fix_resurrections() in Section 8 will later backfill any "D" that is
# contradicted by a later "A" for the same stem.
# ------------------------------------------------------------------------

n_broken_total <- DT_Status[Status == "broken below", .N]
n_broken_with_dbh <- DT_Status[Status == "broken below" & !is.na(DBH), .N]
n_broken_no_dbh <- DT_Status[Status == "broken below" & is.na(DBH), .N]
DT_Status[Status == "broken below" & !is.na(DBH), Status := "alive"]
DT_Status[Status == "broken below" & is.na(DBH), Status := "stem dead"]
cat(sprintf(
  "\n🔧 'broken below' resolved by DBH: %d total → %d alive (DBH present), %d dead (no DBH)\n",
  n_broken_total, n_broken_with_dbh, n_broken_no_dbh
))

# Pivot long to wide: one row per StemID, one column per census
original_status_wide <- dcast(DT_Status,
  formula = StemID ~ census,
  value.var = "Status"
) # fill cells with Status values

# Find reordering indices: match master stem order (unique_StemID) to current rows
idx <- match(
  unique_StemID$StemID, # desired order (from master list)
  original_status_wide$StemID
) # current StemID order in wide table

# Reorder rows to align with unique_StemID
original_status_wide <- original_status_wide[idx] # subset/reorder using idx
# head(original_status_wide)
original_status_full <- as.matrix(original_status_wide[, -1])
original_status <- original_status_full # Create working copy (preserve original)
head(original_status)
# Then proceed using original_status_full_dt as the status matrix

# Show frequency of status values before transformation
cat("\n📊 Full Status values before transformation:\n")
print(table(c(original_status), useNA = "ifany"))

# Show frequency of status by census values before transformation
sum_status_pre <- DT_Status[, .(nobs = .N), by = census:Status][order(Status, census)]
cat("\n📊 Full Status by Census before transformation:\n")
sum_status_pre[, .N, by = .(Status, census)][N > 1]

## Replace english words by corresponding codes ####
# Transform verbose English status terms into single-letter codes
# This standardizes the status vocabulary across different ForestGEO sites

# Step 1: "alive" → "A"
original_status <- gsub("alive", "A", original_status)
# original_status <- gsub("alive-not measured", "A", original_status)

# Step 2: NA (not measured) → "N" (temporary placeholder)
# "N" represents either "Prior" (census not started) or real missing data
# We'll distinguish between these in Section 6
original_status[is.na(original_status)] <- "N"

# Step 2b: Raw "missing" values should be treated as no-data, not dead.
# A raw "missing" observation does not imply the stem was confirmed dead.
original_status[original_status == "missing"] <- "N"

# The site-specific "broken below" field code was already resolved before
# the wide pivot (see above), so this gsub is a defensive no-op kept only
# for documentation.
original_status <- gsub("broken below", "G", original_status)

table(original_status, useNA = "ifany")

# Step 4: Everything else → "D"
# This catches: "dead", "stem dead", and any remaining non-A/N/G statuses.
# After Step 3, "broken below" has already been split into "alive"/"dead"
# using DBH, so the only stems that fall into "D" here are genuine deads
# or terminal statuses that were not explicitly marked as missing/no-data.
# The "D" label is treated as REAL unless the same stem reappears as
# "A" in a later census, in which case fix_resurrections() in Section
# 9 backfills it (the "death" was a temporary measurement gap).
unique(as.vector(original_status))
original_status[] <- ifelse(!original_status %in% c("A", "N", "G"), "D", original_status)

# Show frequency of status codes AFTER transformation
cat("\n📊 Status codes after transformation:\n")
print(table(c(original_status), useNA = "ifany"))

original_status_codes_summary <- data.table(original_status)
original_status_codes_summary <- melt(
  original_status_codes_summary,
  measure.vars = 1:ncol(original_status_codes_summary), # all columns
  variable.name = "census",
  value.name = "status_code"
)
original_status_codes_summary <- original_status_codes_summary[, .(nobs = .N), by = census:status_code][order(status_code, census)]
cat("\n📊 Full Status by Census after transformation:\n")
dcast(original_status_codes_summary,
  formula = status_code ~ census,
  value.var = "nobs"
)[order(status_code)]

print_to_log("Status codes to correct:", log_file, new_message = TRUE)
print_to_log(
  capture.output(
    dcast(original_status_codes_summary,
      formula = status_code ~ census,
      value.var = "nobs"
    )[order(status_code)]
  ),
  log_file,
  new_message = FALSE
)

# --------------------------------------------------------------------
# Create encounter history strings
# --------------------------------------------------------------------
# Convert matrix to concatenated strings (e.g., "AAAD" = alive 3x, dead once)
original_status <- do.call(paste0, as.data.frame(original_status, stringsAsFactors = FALSE))
new_status <- original_status # Working copy for propagation
head(new_status)

# --------------------------------------------------------------------
# Helper function for status pattern analysis
# --------------------------------------------------------------------
sort_table_status <- function(x, sort_by = "Freq", decreasing = TRUE) {
  tbl <- data.frame(table(x))
  if (!sort_by %in% names(tbl)) {
    stop(paste0("Column '", sort_by, "' not found. Available: ", paste(names(tbl), collapse = ", ")))
  }
  tbl[order(tbl[[sort_by]], decreasing = decreasing), ]
}

# --------------------------------------------------------------------
# Display and export encounter histories (before propagation)
# --------------------------------------------------------------------
cat("\n📋 Encounter history patterns (alphabetically):\n")
print(sort_table_status(new_status, sort_by = "x", decreasing = FALSE))

cat("\n📋 Encounter history patterns (by Frequency):\n")
print(sort_table_status(new_status, sort_by = "Freq", decreasing = TRUE))

tbl_sorted_before_propagation <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_before_propagation)
tbl_sorted_before_propagation[, `:=`(
  code_before_propagation = x,
  x = NULL,
  Freq_before_propagation = Freq,
  Freq = NULL,
  rowid = .I
)]

## Export final encounter history patterns before propagation ####
fwrite(tbl_sorted_before_propagation,
  file = file.path(CHECK_folder, paste0("encounter_history_patterns_before_propagation_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

# ========================================================================
# SECTION 5: VALIDATE ENCOUNTER HISTORIES
# ========================================================================
# Check raw stem history strings for impossible biological transitions before
# propagation. This flags histories with invalid codes, illegal state changes,
# or resurrection-like patterns in the raw data.

# Convert to character vector
histories_before_propagation <- as.vector(tbl_sorted_before_propagation$code_before_propagation)
head(histories_before_propagation)
# --------
# Define validation rules
# --------

# Allowed codes
allowed_codes <- c("A", "D", "N", "G")
all_combinations <- expand.grid(first = allowed_codes, second = allowed_codes)
# Get all possible combinations
expand.grid(first = allowed_codes, second = allowed_codes)[all_combinations$first != all_combinations$second, ]

# Define transitions that are biologically impossible
# e.g., Dead → Alive (D→A), Gone → Alive (G→A), etc.
invalid_transitions <- list(
  c("D", "A"), # D is dead, cannot become A
  c("G", "A"), # G is gone, cannot become A
  c("N", "D"), # N is no data, cannot confirm D directly
  c("G", "D"), # G is gone, cannot become D (already gone)
  c("A", "N"), # A is alive, cannot become N (should be measured)
  c("D", "N"), # D is dead, cannot become N (should remain D)
  c("G", "N"), # G is gone, cannot become N (should remain G)
  c("D", "G"), # D is dead, cannot become G (already dead)
  c("N", "G") # N is no data, cannot become G directly
)

# ----------------------------
# 3. Loop through each encounter history
# ----------------------------
check_histories <- function(histories, allowed_codes, invalid_transitions, nchars = length(ViewFullTable_split)) {
  # ----------------------------
  # Function to check stem encounter histories
  # ----------------------------
  # Initialize a data frame to store issues
  issues <- data.frame(
    History = character(),
    Issue = character(),
    stringsAsFactors = FALSE
  )
  # Loop through each history
  for (h in histories) {
    h <- trimws(h)
    chars <- unlist(strsplit(h, ""))
    # Collect issues for this history
    history_issues <- character()
    # ---- Rule 1: Length check ----
    if (length(chars) != nchars) {
      history_issues <- c(history_issues, paste0("Invalid length (should be ", nchars, " characters)"))
    }
    # ---- Rule 2: Allowed characters ----
    if (!all(chars %in% allowed_codes)) {
      invalid_chars <- chars[!chars %in% allowed_codes]
      history_issues <- c(
        history_issues,
        paste0("Contains invalid character(s): ", paste(invalid_chars, collapse = ","))
      )
    }
    # ---- Rule 3: First census logic ----
    if (length(chars) >= 1 && chars[1] %in% c("D", "G")) {
      history_issues <- c(history_issues, "Starts with D or G (cannot start dead or gone)")
    }
    # ---- Rule 4: Invalid direct transitions ----
    if (length(chars) >= 2) {
      for (i in 2:length(chars)) {
        pair <- c(chars[i - 1], chars[i])
        if (any(sapply(invalid_transitions, function(x) all(x == pair)))) {
          history_issues <- c(
            history_issues,
            sprintf("Invalid transition %s→%s at position %d", pair[1], pair[2], i)
          )
        }
      }
    }
    # ---- Rule 5: Irreversibility after death/gone ----
    first_DG_index <- which(chars %in% c("D", "G"))
    if (length(first_DG_index) > 0) {
      first_DG_index <- first_DG_index[1]
      if (first_DG_index < length(chars)) {
        later_states <- chars[(first_DG_index + 1):length(chars)]
        A_positions <- which(later_states == "A") + first_DG_index
        if (length(A_positions) > 0) {
          for (pos in A_positions) {
            history_issues <- c(
              history_issues,
              sprintf("Reappears alive after D/G at position %d", pos)
            )
          }
        }
      }
    }
    # ---- Add all issues for this history to the main data frame ----
    if (length(history_issues) > 0) {
      issues <- rbind(issues, data.frame(
        History = rep(h, length(history_issues)),
        Issue = history_issues,
        stringsAsFactors = FALSE
      ))
    }
  }
  # Return the issues data frame
  return(issues)
}

issues <- check_histories(
  histories = histories_before_propagation,
  allowed_codes = allowed_codes,
  invalid_transitions = invalid_transitions,
  nchars = length(ViewFullTable_split)
)

# ----------------------------
# 4. Print all detected issues
# ----------------------------
if (nrow(issues) == 0) {
  cat("No issues detected in any histories.\n")
} else {
  cat("Detected issues:\n")
  issues <- unique(issues)
  data.frame(sort(unique(issues$Issue)))
}

print_to_log("Detected issues in encounter histories:", log_file, new_message = TRUE)
print_to_log(
  capture.output(
    if (nrow(issues) == 0) {
      cat("No issues detected in any histories.\n")
    } else {
      cat("Detected issues:\n")
      issues <- unique(issues)
      data.frame(sort(unique(issues$Issue)))
    }
  ),
  log_file,
  new_message = FALSE
)

# ========================================================================
# SECTION 6: STATUS PROPAGATION RULES
# ========================================================================
# Resolve placeholder "N" values so each stem history is fully defined.
#   - ^N → P
#   - PN → PP
#   - AN → AD
#   - GN → GG
#   - DN → DD
# This turns partially observed histories into terminal or prior patterns
# suitable for the downstream D/G correction logic.

# --------------------------------------------------------------------
## RULE 1: First census N → P (Prior = not yet recruited)
# --------------------------------------------------------------------
# Example: "NAAA" → "PAAA" (recruited in census 2)
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("^N", tbl_sorted$x), ])
cat("\n Rule 1: First census N → P...")
new_status <- gsub("^N", "P", new_status)
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("^N", tbl_sorted$x), ])
print(tbl_sorted[grepl("^P", tbl_sorted$x), ])

# --------------------------------------------------------------------
## RULE 2: Propagate status 'P' forward until stem is first censused
# --------------------------------------------------------------------
# Pattern: PN → PP (not yet recruited)
print(tbl_sorted[grepl("PN", tbl_sorted$x), ])
cat("\n Applying Rule 2 (PN→PP)...\n")
while (any(grep("PN", new_status))) {
  new_status <- gsub("PN", "PP", new_status)
}
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("PN", tbl_sorted$x), ])
print(tbl_sorted[grepl("P", tbl_sorted$x), ])

# --------------------------------------------------------------------
## RULE 3: Propagate Alive forward (AN → AD, interim assignment)
# --------------------------------------------------------------------
print(tbl_sorted[grepl("AN", tbl_sorted$x), ])
cat("\n Applying Rule 3 (AN→AD)...\n")
while (any(grep("AN", new_status))) {
  new_status <- gsub("AN", "AD", new_status)
}
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("AN", tbl_sorted$x), ])
print(tbl_sorted[grepl("AD", tbl_sorted$x), ])

# --------------------------------------------------------------------
## RULE 4: Propagate status 'G' forward (GN → GG) ####
# --------------------------------------------------------------------
# Once a stem is gone (G), it stays gone unless explicitly recorded otherwise
# Pattern: "AAGN" → "AAGG"
# Example: Stem broken in census 3, still broken in census 4 (no data)
print(tbl_sorted[grepl("GN", tbl_sorted$x), ])
cat("\n Applying Rule 4 (GN→GG)...\n")
while (any(grep("GN", new_status))) {
  new_status <- gsub("GN", "GG", new_status)
}
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("GN", tbl_sorted$x), ])
print(tbl_sorted[grepl("GG", tbl_sorted$x), ])

# --------------------------------------------------------------------
## RULE 5: Propagate status 'D' forward (DN → DD) ####
# --------------------------------------------------------------------
# Once a stem/tree is dead (D), it stays dead unless explicitly recorded otherwise
# Pattern: "AADN" → "AADD"
# Example: Stem died in census 3, still dead in census 4 (no data)
print(tbl_sorted[grepl("DN", tbl_sorted$x), ])
cat("\n Applying Rule 5 (DN→DD)...\n")
while (any(grep("DN", new_status))) {
  new_status <- gsub("DN", "DD", new_status)
}
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("DN", tbl_sorted$x), ])
print(tbl_sorted[grepl("DD", tbl_sorted$x), ])

# --------------------------------------------------------------------
# Final status propagation summary
# --------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("✓✓✓ STATUS PROPAGATION COMPLETE ✓✓✓\n")
cat(strrep("=", 70), "\n\n")

# Print summary statistics
cat("Summary of status codes after propagation:\n")
status_summary <- table(unlist(strsplit(new_status, "")))
print(status_summary)
cat("\nTotal stems processed:", length(new_status), "\n")
cat("Unique encounter patterns:", length(unique(new_status)), "\n\n")
tbl_sorted_after_propagation <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_after_propagation)
tbl_sorted_after_propagation[, `:=`(
  code_after_propagation = x,
  x = NULL,
  Freq_after_propagation = Freq,
  Freq = NULL,
  rowid = .I
)]

## Export final encounter history patterns after propagation ####
fwrite(tbl_sorted_after_propagation,
  file = file.path(CHECK_folder, paste0("encounter_history_patterns_after_propagation_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

## COMPARE BEFORE AND AFTER PROPAGATION
comparison_tbl <- merge(
  tbl_sorted_before_propagation,
  tbl_sorted_after_propagation,
  by = "rowid",
  all = TRUE
)

cat("\n📋 Comparison of encounter history patterns BEFORE and AFTER propagation:\n")
comparison_tbl[Freq_before_propagation == Freq_after_propagation, ]

cat("\n📋 Patterns that CHANGED frequency due to propagation:\n")
comparison_tbl[Freq_before_propagation != Freq_after_propagation, ]

issues <- check_histories(
  histories = unique(as.vector(new_status)),
  allowed_codes = allowed_codes,
  invalid_transitions = invalid_transitions,
  nchars = length(ViewFullTable_split)
)

###############################################################
# VALIDATION: Encounter History Script (With Propagation)
# -------------------------------------------------------------
# Codes:
#   A = Alive
#   D = Dead
#   N = No data
#   G = Gone
#   P / PP / PPP / ... = Prior placeholders (before first measurement)
#
# Biological rules:
#   - Dead (D) or Gone (G) is terminal; cannot become Alive (A) or Prior (P)
#   - Prior (P) can propagate forward until first real measurement
#   - Alive (A) cannot appear after Dead (D) or Gone (G)
#   - Gone (G) propagates forward
#   - Dead (D) propagates forward
###############################################################

# ----------------------------
# 1. Define histories
# ----------------------------

histories_after_propagation <- as.vector(tbl_sorted_after_propagation$code_after_propagation)
unique(unlist(strsplit(histories_after_propagation, "")))

# ----------------------------
# 2. Define allowed codes
# ----------------------------
allowed_codes <- c("A", "D", "G", "P") # Extend as needed

all_combinations <- as.data.table(
  expand.grid(
    first = allowed_codes,
    second = allowed_codes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
)

# -----------------------------------------------------------------------------
# 2. Mark biologically valid transitions
# -----------------------------------------------------------------------------
all_combinations[, possible := fifelse(
  # Alive can stay alive, die, or disappear
  (first == "A" & second %in% c("A", "D", "G")) |
    # Dead stays dead (cannot return or disappear)
    (first == "D" & second == "D") |
    # Gone stays gone
    (first == "G" & second == "G") |
    # Prior can stay prior or transition to first observed alive
    (first == "P" & second %in% c("P", "A")),
  "Possible", "Impossible"
)]

# ----------------------------
# 3. Define invalid transitions
# ----------------------------
# invalid_transitions is a list of biologically impossible transitions between
# stem states:
#   - Dead cannot become alive:       D → A
#   - Gone cannot become alive:       G → A
#   - Gone cannot become dead:        G → D
#   - Prior cannot directly become dead:  P → D
#   - Dead cannot become gone:        D → G
#   - Prior cannot directly become gone:  P → G
#   - Alive cannot revert to prior:   A → P
#   - Dead cannot revert to prior:    D → P
#   - Gone cannot become prior:       G → P

all_combinations[possible == "Possible", .(first, second)]

invalid_transitions <- all_combinations[possible == "Impossible", .(first, second)]
# make invalid_transitions a list of vectors by row
invalid_transitions <- lapply(seq_len(nrow(invalid_transitions)), function(i) {
  as.vector(unlist(invalid_transitions[i, ]))
})

issues <- check_histories(histories_after_propagation, allowed_codes, invalid_transitions, nchars = length(ViewFullTable_split))

# ----------------------------
# 4. Print all detected issues
# ----------------------------
if (nrow(issues) == 0) {
  cat("No issues detected in any histories.\n")
} else {
  cat("Detected issues:\n")
  issues <- unique(issues)
  data.frame(sort(unique(issues$Issue)))
}

# ========================================================================
# SECTION 7: CORRECT PRIOR-TO-DEAD/GONE INCONSISTENCIES
# ========================================================================
# Detect PD/PG transitions and correct them using DBH evidence. If the D/G cell
# has no DBH, rewrite it to P; if the D/G cell has DBH, preserve the measured
# observation and flag the case for review.

# ------------------------------------------------------------------------
# Build DBH matrix aligned with new_status / unique_StemID order
# ------------------------------------------------------------------------
# Rows = stems (in unique_StemID order), Columns = censuses (1..N),
# Values = DBH (numeric, NA when not measured).
DT_DBH <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(StemID, DBH, census = i)]
}), use.names = TRUE, fill = TRUE)
DBHs_dt <- dcast(DT_DBH, formula = StemID ~ census, value.var = "DBH")
DBHs_dt <- DBHs_dt[match(unique_StemID$StemID, DBHs_dt$StemID)]
DBHs <- as.matrix(DBHs_dt[, -1])

# ------------------------------------------------------------------------
# HELPER FUNCTION: fix_PD_PG_inconsistencies
# ------------------------------------------------------------------------
# PURPOSE
#   Detect every (stem, census) cell where the propagated status is D or G
#   AND the previous census' status is P. Optionally use DBH as evidence
#   to decide what to do with each cell.
#
# PARAMETERS
#   status_vec : character vector of encounter histories (one element per
#                stem; each string of length = n_censuses, characters in
#                {A, D, G, P, N}).
#   dbh_matrix : numeric matrix, rows aligned with status_vec, columns =
#                censuses. NA = not measured.
#   stem_ids   : optional character vector of StemIDs (same length as
#                status_vec) used only for the log.
#   dbh_aware  : logical.
#                  * FALSE — naive mode: every PD / PG → PP regardless of DBH
#                    (matches the legacy behavior).
#                  * TRUE  — DBH-aware mode: cells where the D/G has NO DBH
#                    are rewritten to P (PD/PG → PP). Cells where the D/G
#                    HAS a DBH are LEFT UNCHANGED and FLAGGED.
#   log_file   : path of the audit log written for both modes (always
#                lists every PD/PG cell found, with its disposition).
#   verbose    : print summary to console.
#
# RETURNS
#   list(
#     status_vec = corrected character vector,
#     n_changed  = number of cells rewritten to "P",
#     n_flagged  = number of cells flagged (DBH-aware mode only),
#     flagged    = data.table of every PD/PG cell with disposition
#   )
# ------------------------------------------------------------------------
fix_PD_PG_inconsistencies <- function(status_vec,
                                      dbh_matrix,
                                      stem_ids = NULL,
                                      dbh_aware = TRUE,
                                      flag_dbh_action = c("promote_to_A", "keep"),
                                      log_file = file.path(CHECK_folder, "log_PD_PG.txt"),
                                      verbose = TRUE) {
  # ----------------------------------------------------------------------
  # ARGUMENTS
  #   status_vec      : char vector of encounter histories (one per stem)
  #   dbh_matrix      : numeric matrix [stem x census]; NA = not measured
  #   stem_ids        : optional StemIDs for the log
  #   dbh_aware       : TRUE  -> DBH evidence drives the decision per cell
  #                     FALSE -> legacy gsub: every PD/PG -> PP
  #   flag_dbh_action : what to do with PD/PG cells that DO have a DBH
  #                     (only used when dbh_aware = TRUE):
  #                       "promote_to_A" -> rewrite cell to "A" (the stem
  #                                          was alive: there's a DBH)
  #                       "keep"         -> leave cell alone, just log it
  #   log_file        : audit log (one row per unique PD/PG cell)
  #
  # DISPOSITION CODES IN THE LOG
  #   rewrite_to_P : PD/PG without DBH -> P
  #   promote_to_A : PD/PG  with DBH -> A  (dbh_aware + "promote_to_A")
  #   flagged_keep : PD/PG  with DBH left unchanged (dbh_aware + "keep")
  # ----------------------------------------------------------------------
  flag_dbh_action <- match.arg(flag_dbh_action)
  stopifnot(is.character(status_vec))
  status_lengths <- nchar(status_vec)
  if (length(unique(status_lengths)) != 1L) {
    stop("fix_PD_PG_inconsistencies: all encounter histories must have the same length.")
  }
  n_censuses <- status_lengths[1]
  n_stems <- length(status_vec)
  if (!is.matrix(dbh_matrix) || nrow(dbh_matrix) != n_stems || ncol(dbh_matrix) != n_censuses) {
    stop(sprintf(
      "fix_PD_PG_inconsistencies: dbh_matrix must be %d x %d (got %d x %d).",
      n_stems, n_censuses, nrow(dbh_matrix), ncol(dbh_matrix)
    ))
  }
  if (is.null(stem_ids)) stem_ids <- as.character(seq_len(n_stems))

  # IMPLEMENTATION (single matrix scan + bounded vectorized iteration):
  #   - Each pass is one full vectorized matrix scan (O(n_stems * n_censuses)
  #     in C) plus matrix-index assignment of the changed cells.
  #   - Iteration is needed because rewrite_to_P creates new PD/PG cells
  #     one column to the right (PDDD -> PPDD -> PPPD -> PPPP). The loop
  #     is bounded by n_censuses-1 and converges quickly because each pass
  #     touches an ever-smaller subset of cells.
  #   - Per-pass index records are appended to plain integer/numeric lists
  #     and the `flagged` data.table is built ONCE at the end (avoids
  #     per-iteration data.table allocation).
  #
  # Status -> character matrix for cell-level edits (rows=stems, cols=censuses)
  smat <- do.call(rbind, strsplit(status_vec, "", fixed = TRUE))

  # Per-pass record buffers (plain vectors, cheap to grow via list()).
  rec_rows <- vector("list")
  rec_cols <- vector("list")
  rec_codes <- vector("list")
  rec_dbh <- vector("list")
  rec_action <- vector("list")

  n_rewrite_P <- 0L
  n_promote_A <- 0L
  n_flag_keep <- 0L
  iter <- 0L

  repeat {
    iter <- iter + 1L
    prev <- smat[, 1:(n_censuses - 1), drop = FALSE]
    curr <- smat[, 2:n_censuses, drop = FALSE]

    hits <- which(prev == "P" & (curr == "D" | curr == "G"), arr.ind = TRUE)
    if (nrow(hits) == 0L) break

    cell_rows <- hits[, 1]
    cell_cols <- hits[, 2] + 1L
    cell_codes <- smat[cbind(cell_rows, cell_cols)] # "D" or "G"
    cell_dbh <- dbh_matrix[cbind(cell_rows, cell_cols)]
    has_dbh <- !is.na(cell_dbh)

    if (!dbh_aware) {
      action <- rep("rewrite_to_P", length(cell_rows))
    } else {
      action <- ifelse(has_dbh,
        if (flag_dbh_action == "promote_to_A") "promote_to_A" else "flagged_keep",
        "rewrite_to_P"
      )
    }

    # On flagged_keep mode, cells with DBH stay PD/PG forever -> the loop
    # would never end. Detect: if every hit is flagged_keep, we're done
    # after recording them once.
    is_P <- action == "rewrite_to_P"
    is_A <- action == "promote_to_A"
    is_K <- action == "flagged_keep"

    if (any(is_P)) smat[cbind(cell_rows[is_P], cell_cols[is_P])] <- "P"
    if (any(is_A)) smat[cbind(cell_rows[is_A], cell_cols[is_A])] <- "A"

    n_rewrite_P <- n_rewrite_P + sum(is_P)
    n_promote_A <- n_promote_A + sum(is_A)

    # Record only cells that were actually edited this pass (rewrite/promote);
    # flagged_keep cells are recorded once on the final pass below.
    keep_now <- is_P | is_A
    if (any(keep_now)) {
      k <- length(rec_rows) + 1L
      rec_rows[[k]] <- cell_rows[keep_now]
      rec_cols[[k]] <- cell_cols[keep_now]
      rec_codes[[k]] <- cell_codes[keep_now]
      rec_dbh[[k]] <- cell_dbh[keep_now]
      rec_action[[k]] <- action[keep_now]
    }

    # If this pass produced no edits (all hits were flagged_keep), record
    # those once and exit.
    if (!any(is_P) && !any(is_A)) {
      if (any(is_K)) {
        n_flag_keep <- n_flag_keep + sum(is_K)
        k <- length(rec_rows) + 1L
        rec_rows[[k]] <- cell_rows[is_K]
        rec_cols[[k]] <- cell_cols[is_K]
        rec_codes[[k]] <- cell_codes[is_K]
        rec_dbh[[k]] <- cell_dbh[is_K]
        rec_action[[k]] <- action[is_K]
      }
      break
    }
  }

  # Build flagged data.table ONCE from the accumulated index buffers.
  flagged <- if (length(rec_rows) > 0L) {
    all_rows <- unlist(rec_rows, use.names = FALSE)
    data.table(
      stem_idx = all_rows,
      StemID   = stem_ids[all_rows],
      census   = unlist(rec_cols, use.names = FALSE),
      pattern  = paste0("P", unlist(rec_codes, use.names = FALSE)),
      DBH      = unlist(rec_dbh, use.names = FALSE),
      action   = unlist(rec_action, use.names = FALSE)
    )
  } else {
    data.table(
      stem_idx = integer(), StemID = character(),
      census = integer(), pattern = character(), DBH = numeric(),
      action = character()
    )
  }

  # Reassemble status vector (vectorized; ~50x faster than apply)
  new_vec <- do.call(paste0, lapply(seq_len(n_censuses), function(j) smat[, j]))

  # Write log
  if (!is.null(log_file)) {
    con <- file(log_file, open = "w")
    on.exit(close(con), add = TRUE)
    writeLines(c(
      paste0("# log_PD_PG.txt  - generated ", format(Sys.time())),
      paste0("# mode            : ", if (dbh_aware) "DBH-aware" else "naive (gsub-equivalent)"),
      paste0("# flag_dbh_action : ", flag_dbh_action),
      paste0("# n_stems         : ", n_stems),
      paste0("# n_censuses      : ", n_censuses),
      paste0("# iterations      : ", iter),
      paste0("# n_rewrite_to_P  : ", n_rewrite_P, "  (PD/PG without DBH -> P)"),
      paste0("# n_promote_to_A  : ", n_promote_A, "  (PD/PG  with DBH -> A)"),
      paste0("# n_flagged_keep  : ", n_flag_keep, "  (PD/PG  with DBH left unchanged)"),
      "#"
    ), con)
    if (nrow(flagged) > 0L) {
      writeLines("# per-cell disposition (one row per unique PD/PG cell):", con)
      write.table(flagged, con,
        sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE
      )
    } else {
      writeLines("# no PD/PG cells found.", con)
    }
  }

  if (verbose) {
    cat(sprintf(
      "fix_PD_PG_inconsistencies: mode=%s flag_dbh_action=%s | rewrite_to_P=%d | promote_to_A=%d | flagged_keep=%d | iters=%d\n",
      if (dbh_aware) "DBH-aware" else "naive",
      flag_dbh_action, n_rewrite_P, n_promote_A, n_flag_keep, iter
    ))
    if (n_flag_keep > 0L) {
      cat(sprintf(
        "  ⚠ %d PD/PG cell(s) with DBH left unchanged. See: %s\n",
        n_flag_keep, log_file
      ))
    }
    if (n_promote_A > 0L) {
      cat(sprintf(
        "  ℹ %d PD/PG cell(s) with DBH promoted to A. See: %s\n",
        n_promote_A, log_file
      ))
    }
  }

  list(
    status_vec = new_vec,
    n_rewrite_P = n_rewrite_P,
    n_promote_A = n_promote_A,
    n_flag_keep = n_flag_keep,
    flagged = flagged
  )
}

# ------------------------------------------------------------------------
# Apply the correction
# ------------------------------------------------------------------------
cat("\n🔍 PD / PG patterns BEFORE correction:\n")
tbl_sorted_stem <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted_stem[grepl("PD|PG", tbl_sorted_stem$x), ])

# DBH-aware behaviour:
#   PD / PG cell with no DBH        → rewritten to P (rewrite_to_P)
#                                     → the "D" or "G" was a propagated
#                                       N→D/G that has no measurement to
#                                       support it; demote back to Prior.
#   PD / PG cell WITH a DBH         → promoted to A (promote_to_A)
#                                     → a measurement exists this census,
#                                       so the stem was alive (typically a
#                                       "broken below" record that already
#                                       went through the Section-4 rule).
#   flag_dbh_action = "keep"        → disabled here; would log the cell
#                                     and leave it untouched ("flagged_keep").
# Toggle dbh_aware = FALSE to reproduce the legacy gsub("PD|PG","PP")
# behaviour (no DBH inspection, all PD/PG → PP).
PD_PG_fix <- fix_PD_PG_inconsistencies(
  status_vec = new_status,
  dbh_matrix = DBHs,
  stem_ids = unique_StemID$StemID,
  dbh_aware = TRUE,
  log_file = file.path(CHECK_folder, "log_PD_PG.txt"),
  verbose = TRUE,
  flag_dbh_action = "promote_to_A"
)

new_status <- PD_PG_fix$status_vec

cat("\n🔍 PD / PG patterns AFTER correction (should be empty in naive mode;\n")
cat("   may still contain a few flagged cells in DBH-aware mode):\n")
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("PD|PG", tbl_sorted$x), ])

# ========================================================================
# SECTION 8: RESOLVE RESURRECTIONS (D→A, G→A)
# ========================================================================
# Correct impossible resurrection patterns using DBH evidence:
#   - A with DBH after D/G means the earlier D/G was wrong and the stem is
#     backfilled to A.
#   - A without DBH after D/G means the alive record is suspect and is
#     demoted to D or G.
# Every case is audited in log_resurrections.txt.

# ------------------------------------------------------------------------
# HELPER FUNCTION: fix_resurrections
# ------------------------------------------------------------------------
fix_resurrections <- function(status_vec,
                              dbh_matrix,
                              stem_ids = NULL,
                              dbh_aware = TRUE,
                              log_file = file.path(CHECK_folder, "log_resurrections.txt"),
                              verbose = TRUE) {
  # PARAMETERS
  #   status_vec : char vector of encounter histories (one per stem)
  #   dbh_matrix : numeric matrix [stem x census], NA = not measured
  #   stem_ids   : optional StemIDs for the log
  #   dbh_aware  : TRUE  -> per-cell decision uses DBH evidence
  #                FALSE -> legacy behavior: every DA/GA -> AA (backfill all)
  #   log_file   : audit log path
  #
  # RETURNS list with corrected status_vec + counts + flagged data.table.
  #
  # IMPLEMENTATION (single matrix pass + vectorized string ops):
  #   1. Build status matrix once.
  #   2. Find ALL (i, j) cells where smat[i, j-1] in {D,G} and smat[i, j] = A.
  #      Each such cell is an immediate D->A / G->A transition.
  #   3. Classify each cell by DBH at j:
  #        has DBH  -> backfill_to_A   (the A is real; earlier D/G is wrong)
  #        no  DBH  -> demote_to_DG    (the A is a zombie)
  #   4. Apply demotes by direct matrix-index assignment (single vectorized op).
  #   5. Apply backfills via gsub("DA|GA","AA",...) loop on the string vector
  #      (only on rows that actually need backfilling; vectorized C code).
  #   6. Write a single log of all classified cells (no per-pass overhead).
  stopifnot(is.character(status_vec))
  status_lengths <- nchar(status_vec)
  if (length(unique(status_lengths)) != 1L) {
    stop("fix_resurrections: all encounter histories must have the same length.")
  }
  n_censuses <- status_lengths[1]
  n_stems <- length(status_vec)
  if (!is.matrix(dbh_matrix) || nrow(dbh_matrix) != n_stems || ncol(dbh_matrix) != n_censuses) {
    stop(sprintf(
      "fix_resurrections: dbh_matrix must be %d x %d (got %d x %d).",
      n_stems, n_censuses, nrow(dbh_matrix), ncol(dbh_matrix)
    ))
  }
  if (is.null(stem_ids)) stem_ids <- as.character(seq_len(n_stems))

  # ---- single matrix scan -----------------------------------------------
  smat <- do.call(rbind, strsplit(status_vec, "", fixed = TRUE))

  prev <- smat[, 1:(n_censuses - 1), drop = FALSE]
  curr <- smat[, 2:n_censuses, drop = FALSE]
  hits <- which((prev == "D" | prev == "G") & curr == "A", arr.ind = TRUE)

  if (nrow(hits) == 0L) {
    if (verbose) cat("fix_resurrections: no D->A / G->A cells found.\n")
    return(list(
      status_vec = status_vec, n_backfill = 0L, n_demote = 0L,
      flagged = data.table()
    ))
  }

  cell_rows <- hits[, 1]
  cell_cols <- hits[, 2] + 1L # the "A" cell column
  prev_cols <- hits[, 2] # the "D"/"G" cell column
  prev_codes <- smat[cbind(cell_rows, prev_cols)]
  cell_dbh <- dbh_matrix[cbind(cell_rows, cell_cols)]
  has_dbh <- !is.na(cell_dbh)

  if (!dbh_aware) {
    action <- rep("backfill_to_A", length(cell_rows))
  } else {
    action <- ifelse(has_dbh, "backfill_to_A", "demote_to_DG")
  }
  is_back <- action == "backfill_to_A"
  is_dem <- action == "demote_to_DG"
  n_backfill <- sum(is_back)
  n_demote <- sum(is_dem)

  # ---- apply demotes (one vectorized matrix assignment) -----------------
  if (n_demote > 0L) {
    smat[cbind(cell_rows[is_dem], cell_cols[is_dem])] <- prev_codes[is_dem]
  }

  # ---- rebuild string vector (vectorized; ~50x faster than apply) -------
  new_vec <- do.call(paste0, lapply(seq_len(n_censuses), function(j) smat[, j]))

  # ---- apply backfills via gsub on the rows that need it ----------------
  # After demotes, rows with backfill cells still contain "DA" or "GA".
  # gsub propagates DDDA -> DDAA -> DAAA -> AAAA in a small loop on a
  # subset of rows only.
  if (n_backfill > 0L) {
    back_rows <- unique(cell_rows[is_back])
    sub <- new_vec[back_rows]
    while (any(grepl("DA|GA", sub, fixed = FALSE))) {
      sub <- gsub("DA|GA", "AA", sub)
    }
    new_vec[back_rows] <- sub
  }

  # ---- build log --------------------------------------------------------
  flagged <- data.table(
    stem_idx = cell_rows,
    StemID   = stem_ids[cell_rows],
    census   = cell_cols,
    pattern  = paste0(prev_codes, "A"),
    DBH      = cell_dbh,
    action   = action
  )

  if (!is.null(log_file)) {
    con <- file(log_file, open = "w")
    on.exit(close(con), add = TRUE)
    writeLines(c(
      paste0("# log_resurrections.txt - generated ", format(Sys.time())),
      paste0("# mode            : ", if (dbh_aware) "DBH-aware" else "naive (gsub-equivalent)"),
      paste0("# n_stems         : ", n_stems),
      paste0("# n_censuses      : ", n_censuses),
      paste0("# n_backfill_to_A : ", n_backfill, "  (D/G rewritten to A; cell had DBH)"),
      paste0("# n_demote_to_DG  : ", n_demote, "  (A rewritten to D/G; cell had no DBH)"),
      "#"
    ), con)
    if (nrow(flagged) > 0L) {
      writeLines("# per-cell disposition (one row per D->A or G->A cell):", con)
      write.table(flagged, con,
        sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE
      )
    } else {
      writeLines("# no resurrection cells found.", con)
    }
  }

  if (verbose) {
    cat(sprintf(
      "fix_resurrections: mode=%s | backfill_to_A=%d | demote_to_DG=%d\n",
      if (dbh_aware) "DBH-aware" else "naive",
      n_backfill, n_demote
    ))
    if (n_demote > 0L) {
      cat(sprintf(
        "  ⚠ %d zombie A cell(s) demoted (no DBH after D/G). See: %s\n",
        n_demote, log_file
      ))
    }
  }

  list(
    status_vec = new_vec,
    n_backfill = n_backfill,
    n_demote = n_demote,
    flagged = flagged
  )
}

# ------------------------------------------------------------------------
# Apply the correction
# ------------------------------------------------------------------------
cat("\n🔍 D→A / G→A patterns BEFORE correction:\n")
tbl_sorted_stem <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted_stem[grepl("DA|GA", tbl_sorted_stem$x), ])

# DBHs matrix was built earlier in Section 7 and is still aligned with
# unique_StemID / new_status.
RES_fix <- fix_resurrections(
  status_vec = new_status,
  dbh_matrix = DBHs,
  stem_ids   = unique_StemID$StemID,
  dbh_aware  = TRUE,
  log_file   = file.path(CHECK_folder, "log_resurrections.txt"),
  verbose    = TRUE
)
new_status <- RES_fix$status_vec

cat("\n🔍 D→A / G→A patterns AFTER correction (should be empty; if not, next section):\n")
tbl_sorted_stem <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted_stem[grepl("DA|GA", tbl_sorted_stem$x), ])

# ========================================================================
# SECTION 9: SAFETY NET — D/G BETWEEN TWO A’s
# ========================================================================
# Rewrite any D or G strictly between two alive observations to A. This
# final pass removes residual stranded dead/gone codes that are impossible
# given the stem was alive before and after the gap.

fix_DG_between_A <- function(status_vec,
                             dbh_matrix,
                             stem_ids = NULL,
                             log_file = file.path(
                               CHECK_folder,
                               "log_DG_between_A.txt"
                             ),
                             verbose = TRUE) {
  stopifnot(is.character(status_vec))
  status_lengths <- nchar(status_vec)
  if (length(unique(status_lengths)) != 1L) {
    stop("fix_DG_between_A: histories must all have the same length.")
  }
  n_censuses <- status_lengths[1]
  n_stems <- length(status_vec)
  if (!is.matrix(dbh_matrix) || nrow(dbh_matrix) != n_stems ||
    ncol(dbh_matrix) != n_censuses) {
    stop("fix_DG_between_A: dbh_matrix shape mismatch.")
  }
  if (is.null(stem_ids)) stem_ids <- as.character(seq_len(n_stems))

  smat <- do.call(rbind, strsplit(status_vec, "", fixed = TRUE))

  # First and last column-index of "A" in each row. Rows with < 2 A's
  # have no "between" region and are skipped.
  is_A <- smat == "A"
  any_A <- rowSums(is_A) >= 2L
  if (!any(any_A)) {
    if (verbose) cat("fix_DG_between_A: no stem has >= 2 A cells.\n")
    return(list(
      status_vec = status_vec, n_changed = 0L,
      flagged = data.table()
    ))
  }

  first_A <- max.col(is_A, ties.method = "first")
  last_A <- max.col(is_A, ties.method = "last")

  # Build a mask of "strictly between first and last A" for rows with >= 2 A.
  col_idx <- matrix(seq_len(n_censuses),
    nrow = n_stems,
    ncol = n_censuses, byrow = TRUE
  )
  between_msk <- any_A & (col_idx > first_A) & (col_idx < last_A)

  # Cells to rewrite: those that are D or G inside the between region.
  fix_msk <- between_msk & (smat == "D" | smat == "G")
  hits <- which(fix_msk, arr.ind = TRUE)
  if (nrow(hits) == 0L) {
    if (verbose) cat("fix_DG_between_A: no D/G stranded between A's.\n")
    return(list(
      status_vec = status_vec, n_changed = 0L,
      flagged = data.table()
    ))
  }

  cell_rows <- hits[, 1]
  cell_cols <- hits[, 2]
  prev_codes <- smat[cbind(cell_rows, cell_cols)] # "D" or "G"
  cell_dbh <- dbh_matrix[cbind(cell_rows, cell_cols)]

  # Apply the rewrite in one vectorized matrix assignment.
  smat[cbind(cell_rows, cell_cols)] <- "A"

  # Reassemble the string vector.
  new_vec <- do.call(
    paste0,
    lapply(seq_len(n_censuses), function(j) smat[, j])
  )

  flagged <- data.table(
    stem_idx = cell_rows,
    StemID = stem_ids[cell_rows],
    census = cell_cols,
    was = prev_codes,
    DBH = cell_dbh,
    has_DBH = !is.na(cell_dbh),
    first_A_at = first_A[cell_rows],
    last_A_at = last_A[cell_rows]
  )

  n_changed_D <- sum(prev_codes == "D")
  n_changed_G <- sum(prev_codes == "G")

  if (!is.null(log_file)) {
    con <- file(log_file, open = "w")
    on.exit(close(con), add = TRUE)
    writeLines(c(
      paste0("# log_DG_between_A.txt  - generated ", format(Sys.time())),
      paste0("# rule           : any D/G strictly between two A cells in the"),
      paste0("#                  same stem history is rewritten to A"),
      paste0("# n_stems        : ", n_stems),
      paste0("# n_censuses     : ", n_censuses),
      paste0("# n_D_to_A       : ", n_changed_D),
      paste0("# n_G_to_A       : ", n_changed_G),
      paste0("# n_total_cells  : ", nrow(flagged)),
      paste0(
        "# n_with_no_DBH  : ", sum(!flagged$has_DBH),
        "  (measurement gap; downstream user must interpolate)"
      ),
      paste0(
        "# n_with_DBH     : ", sum(flagged$has_DBH),
        "  (DBH was actually recorded; status code was wrong)"
      ),
      "#"
    ), con)
    write.table(flagged, con,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = TRUE
    )
  }

  if (verbose) {
    cat(sprintf(
      "fix_DG_between_A: D->A=%d, G->A=%d (%d cells; %d had DBH, %d were measurement gaps).\n",
      n_changed_D, n_changed_G, nrow(flagged),
      sum(flagged$has_DBH), sum(!flagged$has_DBH)
    ))
  }

  list(
    status_vec = new_vec,
    n_changed = nrow(flagged),
    flagged = flagged
  )
}

DGBA_fix <- fix_DG_between_A(
  status_vec = new_status,
  dbh_matrix = DBHs,
  stem_ids   = unique_StemID$StemID,
  log_file   = file.path(CHECK_folder, "log_DG_between_A.txt"),
  verbose    = TRUE
)
new_status <- DGBA_fix$status_vec

cat("\n🔍 D / G between A patterns AFTER safety-net (should be empty):\n")
tbl_sorted_stem <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted_stem[grepl("A[DG]+A", tbl_sorted_stem$x), ])

cat("\n✓✓✓ STATUS CORRECTION AFTER PROPAGATION COMPLETE ✓✓✓\n\n")

tbl_sorted_correction_after_propagation <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_correction_after_propagation)
tbl_sorted_correction_after_propagation[, `:=`(
  code_correction_after_propagation = x,
  x = NULL,
  Freq_correction_after_propagation = Freq,
  Freq = NULL,
  rowid = .I
)]

## Export final encounter history patterns after propagation ####
fwrite(tbl_sorted_correction_after_propagation,
  file = file.path(CHECK_folder, paste0("encounter_history_patterns_correction_after_propagation_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

###############################################################
# VALIDATION: Encounter History Script (With correction after Propagation)
# -------------------------------------------------------------
# Codes:
#   A = Alive
#   D = Dead
#   N = No data
#   G = Gone
#   P / PP / PPP / ... = Prior placeholders (before first measurement)
#
# Biological rules:
#   - Dead (D) or Gone (G) is terminal; cannot become Alive (A) or Prior (P)
#   - Prior (P) can propagate forward until first real measurement
#   - Alive (A) cannot appear after Dead (D) or Gone (G)
#   - Gone (G) propagates forward
#   - Dead (D) propagates forward
###############################################################

# ----------------------------
# 1. Define histories
# ----------------------------

histories_correction_after_propagation <- as.vector(tbl_sorted_correction_after_propagation$code_correction_after_propagation)

# ----------------------------
# 2. Define allowed codes
# ----------------------------
allowed_codes <- c("A", "D", "G", "P") # Extend as needed

all_combinations <- as.data.table(
  expand.grid(
    first = allowed_codes,
    second = allowed_codes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
)

# -----------------------------------------------------------------------------
# 2. Mark biologically valid transitions (post-correction)
# -----------------------------------------------------------------------------
# After Section 10.4, GD and DG are additional valid transitions:
#   GD : gone → dead  (tree fully dies at the next census)
#   DG : dead → gone  (new recruit makes tree alive, reclassifying stem)
all_combinations[, possible := fifelse(
  # Alive can stay alive, die, or disappear
  (first == "A" & second %in% c("A", "D", "G")) |
    # Dead stays dead OR becomes G (tree comes back alive via new recruit)
    (first == "D" & second %in% c("D", "G")) |
    # Gone stays gone OR becomes D (tree fully dies)
    (first == "G" & second %in% c("G", "D")) |
    # Prior can stay prior or transition to first observed alive
    (first == "P" & second %in% c("P", "A")),
  "Possible", "Impossible"
)]

# ----------------------------
# 3. Define transitions
# ----------------------------
# Extract valid transitions as two-character strings
valid_trans <- all_combinations[possible == "Possible", paste0(first, second)]
invalid_transitions <- all_combinations[possible == "Impossible", .(first, second)]
invalid_transitions <- split(invalid_transitions, seq(nrow(invalid_transitions)))

issues <- check_histories(
  histories = histories_correction_after_propagation,
  allowed_codes = allowed_codes,
  invalid_transitions = invalid_transitions,
  nchars = length(ViewFullTable_split)
)

# ----------------------------
# 4. Print all detected issues
# ----------------------------
if (nrow(issues) == 0) {
  cat("No issues detected in any histories.\n")
} else {
  cat("Detected issues:\n")
  issues <- unique(issues)
  data.frame(sort(unique(issues$Issue)))
}

# -----------------------------
# Diagram of allowed stem transitions
# -----------------------------

library(igraph)

# make valid_trans edges vector
valid_edges <- unlist(strsplit(valid_trans, ""))
invalid_edges <- unlist(invalid_transitions)

# Create graph
g_valid <- make_graph(edges = valid_edges, directed = TRUE)
g_invalid <- make_graph(edges = invalid_edges, directed = TRUE)

# Plot settings
png(
  filename = file.path(CHECK_folder, paste0("stem_transition_diagrams_", site, ".png")),
  width = 12,
  height = 6,
  units = "in",
  res = 100
)
par(mfrow = c(1, 2)) # side-by-side plots
plot(
  g_valid,
  vertex.size = 40,
  vertex.label.cex = 1.5,
  vertex.color = "green",
  edge.arrow.size = 0.8,
  main = "Allowed Stem Transitions (A, D, G, P)"
)

plot(
  g_invalid,
  vertex.size = 40,
  vertex.label.cex = 1.5,
  vertex.color = "red",
  edge.arrow.size = 0.8,
  main = "Not-allowed Stem Transitions (A, D, G, P)"
)
dev.off()

# ========================================================================
# SECTION 10: TREE-LEVEL STATUS CALCULATION AND D/G CORRECTION
# ========================================================================
# Aggregate stem histories by tree and enforce the D/G semantics
# at EACH CENSUS independently:
#   - single-stem trees       : G -> D (gone ≡ dead)
#   - multi-stem, tree NOT fully dead (has A or P at census j)
#                             : D -> G (stem lost, tree still alive)
#   - multi-stem, tree fully dead (all stems D or G, no A/P at census j)
#                             : G -> D (permanent death)
# A and P cells are never modified. This produces GD and DG as
# valid transitions in corrected_new_status (see Section 10.4).

#--------------------------------------------------------------
# 10.1  Tag <-> TreeID consistency check
#--------------------------------------------------------------
# Each Tag should map to exactly one TreeID (and vice versa). Anything
# else means the master stem table has duplicate identifiers and tree
# grouping below would be wrong.
tag_treeid_dt <- unique(unique_StemID[, .(Tag, TreeID)])
tags_per_treeid <- tag_treeid_dt[, .N, by = TreeID][N > 1L]
treeids_per_tag <- tag_treeid_dt[, .N, by = Tag][N > 1L]

# NOTE: i fixed this issue by combining treeid_tag in treeid earlier in the script

if (nrow(tags_per_treeid) == 0L && nrow(treeids_per_tag) == 0L) {
  cat("✓ Tag <-> TreeID mapping is consistent (1:1).\n")
} else {
  cat("⚠️  Inconsistent Tag <-> TreeID mapping detected.\n")
  if (nrow(treeids_per_tag) > 0L) {
    cat(sprintf("   %d Tag(s) map to multiple TreeIDs.\n", nrow(treeids_per_tag)))
  }
  if (nrow(tags_per_treeid) > 0L) {
    cat(sprintf("   %d TreeID(s) map to multiple Tags.\n", nrow(tags_per_treeid)))
  }
  inconsistent_tags <- treeids_per_tag$Tag
  inconsistent_treeids <- tags_per_treeid$TreeID
  check_inconsistent_tree_codes <- unique_StemID[
    Tag %in% inconsistent_tags | TreeID %in% inconsistent_treeids,
    .(Tag, TreeID, StemID)
  ]
}

#--------------------------------------------------------------
# 10.2  Group stem-level new_status by TreeID (one tree = one element)
#--------------------------------------------------------------
DT_ns <- data.table(TreeID = unique_StemID$TreeID, new_status = new_status)
new_status_split <- DT_ns[, .(new_status_list = list(new_status)), by = TreeID]
new_status_split <- setNames(
  new_status_split[["new_status_list"]],
  new_status_split[["TreeID"]]
)

# ========================================================================
# Exploration of all possible tree life-history sequences
# ========================================================================

# ========================================================================
# TREE-LEVEL ANALYSIS: HELPER FUNCTIONS AND VALIDATION
# ========================================================================
# This section defines functions to compute tree-level status from stem data
# and validates biological plausibility of life-history sequences.

# --------
# Define valid status transitions
# --------
allowed_codes <- c("A", "D", "G", "P")
all_combinations <- as.data.table(expand.grid(
  first = allowed_codes,
  second = allowed_codes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
))

# Mark biologically valid transitions
# GD (gone→dead) and DG (dead→gone) are valid after the Section 10.4
# per-census D↔G correction.
all_combinations[, possible := fifelse(
  (first == "A" & second %in% c("A", "D", "G")) |
    (first == "D" & second %in% c("D", "G")) |
    (first == "G" & second %in% c("G", "D")) |
    (first == "P" & second %in% c("P", "A")),
  "Possible", "Impossible"
)]

valid_trans <- all_combinations[possible == "Possible", paste0(first, second)]

# --------
# Function: valid_seq - Check if stem sequence is biologically valid
# --------
valid_seq <- function(x, valid_trans) {
  if (!("A" %in% x)) {
    return(FALSE)
  } # Must have at least one "A"
  pairs <- paste0(x[-length(x)], x[-1])
  all(pairs %in% valid_trans)
}

# --------
# Function: tree_state - Compute tree-level status from stem matrix
# --------
# KEY DECISION POINTS:
#   1. If ANY stem is "A" → tree is "A"
#   2. If all stems are "P" → tree is "P"
#   3. If mix of P and D/G:
#      - Check future censuses for any "A"
#      - If future A exists → tree is "A" (was alive, just not all stems recruited)
#      - If no future A → tree is "D" (never became alive)
#   4. If only D/G (no P, no A) → tree is "D"
# ---------------------------------------------------------------------
tree_state <- function(stem_states) {
  n_censuses <- nrow(stem_states)
  tree_seq <- character(n_censuses)
  for (t in 1:n_censuses) {
    row <- stem_states[t, ]
    if ("A" %in% row) {
      # Any alive stem → tree alive
      tree_seq[t] <- "A"
    } else if (all(row == "P")) {
      # All stems unobserved → tree unobserved
      tree_seq[t] <- "P"
    } else if (any(row == "P") && any(row %in% c("D", "G"))) {
      # P + D/G combination → tree exists only if future A exists
      if (t < n_censuses && any(stem_states[(t + 1):n_censuses, ] == "A")) {
        tree_seq[t] <- "A" # tree exists due to future alive stem
      } else {
        tree_seq[t] <- "D" # no future alive → tree dead
      }
    } else if (any(row %in% c("D", "G"))) {
      # Only dead/gone → tree dead
      tree_seq[t] <- "D"
    } else {
      # Fallback → treat as P
      tree_seq[t] <- "P"
    }
  }
  paste0(tree_seq, collapse = "")
}

# ---------------------------------------------------------------------
# Function: tree_exists_check
# ---------------------------------------------------------------------
# Validates whether a tree's life history is biologically plausible.
#
# VALIDATION RULES:
#   - tree with any "A" (alive) stem is valid at that census
#   - tree with only "P" stems is valid (not yet recruited)
#   - tree with P+D/G combination is ONLY valid if future A exists
#     (i.e., the tree eventually became alive, so P+D/G just means
#     partial recruitment - some stems recruited, others didn't yet)
#   - tree with only D/G stems is valid (all dead/gone)
#
# PARAMETERS:
#   stem_states: Matrix where rows = censuses, columns = stems
#
# RETURNS:
#   TRUE if tree life history is biologically valid
#   FALSE if any census violates biological rules
#
# WHY THIS CHECK:
#   Prevents impossible scenarios like:
#   - A stem marked "P" (will appear later) but "D" (already dead)
#     with no future "A" to prove the tree ever existed
# ---------------------------------------------------------------------
tree_exists_check <- function(stem_states) {
  n_censuses <- nrow(stem_states)
  # Logical vector to store validity per census
  valid <- logical(n_censuses)
  for (t in 1:n_censuses) {
    row <- stem_states[t, ]
    if ("A" %in% row) {
      # Any alive stem → census valid
      valid[t] <- TRUE
    } else if (any(row == "P")) {
      # If P exists, check for future alive stems
      if (t < n_censuses && any(stem_states[(t + 1):n_censuses, ] == "A")) {
        valid[t] <- TRUE # tree exists due to future alive stem
      } else if (any(row %in% c("D", "G"))) {
        # P + D/G and no future A → invalid.
        # NOTE (known limitation): this flags trees ending in P+D/G with
        # no future census as biologically impossible. Some are real
        # (the tree was never recorded alive within the observation
        # window) so this check is intentionally NOT enforced inside
        # compute_tree_for_row(); see RUN_TREE_EXISTS_DIAGNOSTIC below.
        valid[t] <- FALSE
      } else {
        # Only P → valid
        valid[t] <- TRUE
      }
    } else if (any(row %in% c("D", "G"))) {
      # Dead/gone stems only → valid
      valid[t] <- TRUE
    } else {
      # Fallback
      valid[t] <- TRUE
    }
  }
  # Return TRUE only if all censuses are valid
  all(valid)
}

# ---------------------------------------------------------------------
# Function: compute_tree_for_row
# ---------------------------------------------------------------------
# Computes the overall tree sequence from a vector of stem histories.
# Steps:
# 1. Converts each stem string to a character matrix
#    (rows = censuses, cols = stems)
# 2. Checks if the tree is biologically plausible
# 3. If valid, computes the tree sequence using tree_state()
# ---------------------------------------------------------------------
compute_tree_for_row <- function(stem_strings) {
  # Convert vector of stem sequences into a matrix
  mat <- do.call(
    cbind,
    lapply(stem_strings, function(s) strsplit(as.character(s), "")[[1]])
  )
  # NOTE: tree_exists_check() is intentionally DISABLED here.
  # When enabled it returns NA for any tree ending in P + D/G with no
  # future A, which is too strict (some such trees are real — a stem
  # was tagged but never measured alive within the observation window).
  # The strict variant is kept as compute_tree_for_row_checking() below
  # and is exercised by the optional RUN_TREE_EXISTS_DIAGNOSTIC block in
  # Section 10.3 to quantify how many trees would be flagged.
  # if (!tree_exists_check(mat)) {
  #   return(NA_character_) # return NA if tree invalid
  # }
  # Compute tree sequence
  tree_state(mat)
}

compute_tree_for_row_checking <- function(stem_strings) {
  # Strict variant of compute_tree_for_row() that DOES enforce the
  # tree_exists_check() guard. Used only by the optional diagnostic in
  # Section 10.3 (RUN_TREE_EXISTS_DIAGNOSTIC). See note in the standard
  # variant above for why it is not the default.
  mat <- do.call(
    cbind,
    lapply(stem_strings, function(s) strsplit(as.character(s), "")[[1]])
  )
  if (!tree_exists_check(mat)) {
    return(NA_character_) # return NA if tree invalid
  }
  # Compute tree sequence
  tree_state(mat)
}

# ========================================================================
# 10.3  APPLY TREE-LEVEL CALCULATION TO ACTUAL DATA
# ========================================================================

# Compute tree-level encounter history for every TreeID (one pass).
tree_histories_list <- lapply(new_status_split, compute_tree_for_row)

# Optional diagnostic: re-run with the strict tree_exists_check() guard
# and report which trees would have been flagged invalid. OFF by default
# because the strict check returns NA for trees ending in P+D/G (known
# false positive).
RUN_TREE_EXISTS_DIAGNOSTIC <- FALSE
if (RUN_TREE_EXISTS_DIAGNOSTIC) {
  tree_histories_list_checking <- lapply(new_status_split, compute_tree_for_row_checking)
  diff_indices <- which(!mapply(identical, tree_histories_list, tree_histories_list_checking))
  cat(sprintf(
    "  diagnostic: tree_exists_check would change %d / %d tree histories.\n",
    length(diff_indices), length(tree_histories_list)
  ))
  rm(tree_histories_list_checking)
}

# Validate results: every tree history must obey the allowed transitions.
is_valid <- vapply(
  tree_histories_list,
  function(s) valid_seq(strsplit(as.character(s), "", fixed = TRUE)[[1]], valid_trans),
  logical(1)
)
cat("Valid tree histories:", sum(is_valid), "out of", length(is_valid), "\n")
if (any(!is_valid)) {
  cat("Invalid values are:\n")
  print(unique(unlist(tree_histories_list[!is_valid])))
}

# Match tree status to each stem's TreeID
tree_histories <- tree_histories_list[as.character(unique_StemID$TreeID)]
# Convert list to matrix (rows = stems, columns = censuses)
tree_histories <- do.call(rbind, tree_histories)
cat("✓ Tree status matched to stem level\n\n")

# ========================================================================
# 10.4  BUILD STATUS MATRIX AND APPLY PER-CENSUS D/G CORRECTION
# ========================================================================
# Convert the per-stem encounter history strings into a character matrix
# (rows = stems, columns = censuses), then apply the D↔4G correction
# vectorised over all stems at once.

# Split each string in 'new_status' into individual characters
split_chars_new_status <- strsplit(new_status, "")
# Combine into matrix while keeping original names (if any)
new_status_matrix <- do.call(rbind, split_chars_new_status)

# Restore rownames if available
if (!is.null(names(split_chars_new_status))) {
  rownames(new_status_matrix) <- names(split_chars_new_status)
} else {
  # fallback: use sequence numbers
  rownames(new_status_matrix) <- seq_along(split_chars_new_status)
}

# ------------------------------------------------------------------
# 10.5  APPLY THE D/G CORRECTION (per-census, vectorized)
# ------------------------------------------------------------------
# RULES (applied independently at each census j):
#   Single-stem trees (all censuses):
#     G -> D  (gone = dead for a single-stem tree)
#
#   Multi-stem trees (per census j):
#     D -> G  if tree is NOT fully dead at j  (has A or P stem)
#             → the stem is lost/broken, not permanently dead;
#               G must come before D in the history
#     G -> D  if tree IS fully dead at j  (all stems are D or G,
#             no A and no P) → gone becomes permanent death
#
# Resulting valid transitions in corrected_new_status:
#   same as new_status PLUS GD (gone→dead) and DG (dead→gone)
# A and P cells are NEVER touched. Enforced by the stopifnot below.
# ------------------------------------------------------------------
stopifnot(nrow(new_status_matrix) == nrow(unique_StemID))

n_cens <- ncol(new_status_matrix)
TreeID_vec <- unique_StemID$TreeID

# Per-stem row info: single-stem vs multi-stem trees.
stem_info_dt <- data.table(
  row_idx = seq_len(nrow(new_status_matrix)),
  TreeID  = TreeID_vec
)
stem_info_dt[, n_stems_in_tree := .N, by = TreeID]
setorder(stem_info_dt, row_idx)

single_stem_row <- stem_info_dt$n_stems_in_tree == 1L
multi_stem_row <- stem_info_dt$n_stems_in_tree > 1L

corrected_new_status_matrix <- new_status_matrix

# ---- (a) Single-stem trees: G → D across ALL censuses --------------------
# For a single-stem tree "gone" is biologically equivalent to "dead".
corrected_new_status_matrix[(new_status_matrix == "G") & single_stem_row] <- "D"

# ---- (b) Multi-stem trees: per-census D ↔ G correction ------------------
# At EACH census independently:
#   D → G : stem is dead but another stem in the same tree is alive (A) →
#           the tree is still alive; the stem is dead but the tree is not permanently dead.
#   G → D : no stem in the tree is alive at this census →
#           the whole tree is dead; "gone" becomes permanent death.
# Order within each census: D→G first, then G→D, so a D that was just
# promoted to G is immediately re-evaluated.
n_D_to_G <- 0L
n_G_to_D <- 0L

for (j in seq_len(n_cens)) {
  col_j <- corrected_new_status_matrix[, j]
  # Per-tree: are ALL stems definitively terminal (D or G only)?
  # A stem with A or P blocks both conversions:
  #   - A means the tree is alive right now.
  #   - P means a future recruit may make the tree alive.
  # Both cases mean the tree is not yet fully dead.
  tree_all_DG <- as.logical(tapply(col_j %in% c("D", "G"), TreeID_vec, all)[TreeID_vec])
  # D → G: tree is NOT fully dead (has A or P) → a dead stem is merely
  # lost/broken; it should be G, not D.
  # This also covers the case where only P stems remain (no A yet): the stem
  # must not carry D before G in the history.
  D_to_G_j <- (col_j == "D") & multi_stem_row & !tree_all_DG
  corrected_new_status_matrix[D_to_G_j, j] <- "G"
  n_D_to_G <- n_D_to_G + sum(D_to_G_j)
  # Re-read column after D→G; recompute tree_all_DG
  col_j <- corrected_new_status_matrix[, j]
  tree_all_DG <- as.logical(tapply(col_j %in% c("D", "G"), TreeID_vec, all)[TreeID_vec])
  # G → D: tree IS fully dead (no A, no P) → gone becomes permanently dead.
  G_to_D_j <- (col_j == "G") & multi_stem_row & tree_all_DG
  corrected_new_status_matrix[G_to_D_j, j] <- "D"
  n_G_to_D <- n_G_to_D + sum(G_to_D_j)
}

# Hard guard: A and P cells must be byte-identical between input and output.
AP_in <- new_status_matrix %in% c("A", "P")
AP_out <- corrected_new_status_matrix %in% c("A", "P")
stopifnot(identical(AP_in, AP_out))
stopifnot(identical(
  new_status_matrix[AP_in],
  corrected_new_status_matrix[AP_in]
))

# Reassemble per-stem string vector (vectorized; ~50x faster than apply).
corrected_new_status <- do.call(
  paste0,
  lapply(seq_len(n_cens), function(j) corrected_new_status_matrix[, j])
)
names(corrected_new_status) <- rownames(new_status_matrix)

cat(sprintf(
  "\u2713 Section 10 D/G correction (per-census): D->G=%d cells, G->D=%d cells (A/P preserved).\n",
  n_D_to_G, n_G_to_D
))

# ---- Per-census D/G consistency validation --------------------------------
# Rule 1: single-stem trees must have no G remaining.
# Rule 2: multi-stem trees — no D may coexist with an A in the same tree
#         at the same census.
# Rule 3: multi-stem trees — no G may remain when no stem in the tree is A
#         at the same census.
cat("\n\U0001F50D Validating per-census D/G consistency...\n")
validation_errors <- 0L
for (j in seq_len(n_cens)) {
  col_j <- corrected_new_status_matrix[, j]
  tree_all_DG <- as.logical(tapply(col_j %in% c("D", "G"), TreeID_vec, all)[TreeID_vec])
  n_single_G <- sum((col_j == "G") & single_stem_row)
  n_D_nonDG <- sum((col_j == "D") & multi_stem_row & !tree_all_DG)
  n_G_all_DG <- sum((col_j == "G") & multi_stem_row & tree_all_DG)
  if (n_single_G > 0L) {
    cat(sprintf("  \u274C Census %d: %d single-stem G cells remain\n", j, n_single_G))
  }
  if (n_D_nonDG > 0L) {
    cat(sprintf("  \u274C Census %d: %d multi-stem D cells in a tree with A or P\n", j, n_D_nonDG))
  }
  if (n_G_all_DG > 0L) {
    cat(sprintf("  \u274C Census %d: %d multi-stem G cells where all stems are D/G\n", j, n_G_all_DG))
  }
  validation_errors <- validation_errors + n_single_G + n_D_nonDG + n_G_all_DG
}

if (validation_errors == 0L) {
  cat("  \u2713 Per-census D/G consistency check passed (0 violations).\n")
} else {
  warning(sprintf("Per-census D/G consistency check FAILED: %d violation(s) found.", validation_errors))
}

# ========================================================================
# 10.6  CHECK STEMS WITH SAME TERMINAL CODE AT FIRST AND LAST CENSUS
# ========================================================================
# Identify stems whose corrected history begins and ends with the same
# terminal code (G or D). These stems are terminal throughout the full
# observation window; their first-census raw metadata is exported for review.
first_code <- corrected_new_status_matrix[, 1]
last_code <- corrected_new_status_matrix[, n_cens]
same_first_last_G <- first_code == "G" & last_code == "G"
same_first_last_D <- first_code == "D" & last_code == "D"
same_first_last <- same_first_last_G | same_first_last_D

cat(sprintf(
  "✓ Section 11.7 check: %d stems start and end with the same terminal code.\n",
  sum(same_first_last)
))
cat(sprintf(
  "   all-G first/last: %d, all-D first/last: %d\n",
  sum(same_first_last_G), sum(same_first_last_D)
))

if (sum(same_first_last) > 0L) {
  first_census_dt <- ViewFullTable_split[[1]][
    , .(StemID, TreeID,
      raw_first_status = Status, raw_first_dbh = DBH,
      raw_first_ListOfTSM = ListOfTSM, raw_first_HOM = HOM
    )
  ]
  idx_first <- match(unique_StemID$StemID, first_census_dt$StemID)

  same_first_last_dt <- data.table(
    TreeID = unique_StemID$TreeID,
    StemID = unique_StemID$StemID,
    first_code = first_code,
    last_code = last_code,
    raw_first_status = first_census_dt$raw_first_status[idx_first],
    raw_first_dbh = first_census_dt$raw_first_dbh[idx_first],
    raw_first_ListOfTSM = first_census_dt$raw_first_ListOfTSM[idx_first],
    raw_first_HOM = first_census_dt$raw_first_HOM[idx_first],
    original_status = original_status,
    corrected_history = corrected_new_status
  )[same_first_last]

  fwrite(same_first_last_dt,
    file.path(CHECK_folder, "section11_7_first_last_terminal_status.csv"),
    sep = ",",
    na = "",
    row.names = FALSE
  )
  cat(sprintf(
    "   details written to: %s\n",
    file.path(CHECK_folder, "section11_7_first_last_terminal_status.csv")
  ))
}

rm(
  stem_info_dt, AP_in, AP_out, single_stem_row, multi_stem_row,
  n_D_to_G, n_G_to_D, TreeID_vec
)

tbl_sorted_new_status <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_new_status)

tbl_sorted_corrected_new_status <- sort_table_status(corrected_new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_corrected_new_status)

tbl_sorted_tree_histories <- sort_table_status(as.vector(tree_histories), sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_tree_histories)

# ========================================================================
# SECTION 11: BIOLOGY ASSESSMENT
# ========================================================================
# Compare original_status, new_status, corrected_new_status, and
# tree_histories for illegal transitions and stranded D/G patterns.
# Write a single QA report summarising pipeline behaviour across all
# four history versions.

assess_biology <- function(status_vec, label, valid_trans) {
  n <- length(status_vec)
  n_cens <- nchar(status_vec[1])
  pairs <- substring(
    rep(status_vec, each = n_cens - 1L),
    rep(seq_len(n_cens - 1L), n),
    rep(seq_len(n_cens - 1L), n) + 1L
  )
  pair_tab <- sort(table(pairs), decreasing = TRUE)
  illegal <- setdiff(names(pair_tab), valid_trans)
  illegal_tab <- pair_tab[illegal]

  # Rows with at least one illegal transition.
  any_illegal <- colSums(matrix(!(pairs %in% valid_trans),
    nrow = n_cens - 1L, ncol = n
  )) > 0L
  n_bad_rows <- sum(any_illegal)

  # "Stranded D/G between two A" pattern.
  n_stranded <- sum(grepl("A[DG]+A", status_vec))

  list(
    label        = label,
    n_stems      = n,
    n_bad_rows   = n_bad_rows,
    n_stranded   = n_stranded,
    illegal_tab  = illegal_tab,
    legal_tab    = pair_tab[intersect(names(pair_tab), valid_trans)]
  )
}

versions <- list(
  list(label = "1. original_status      (raw, pre-propagation)", vec = original_status),
  list(label = "2. new_status           (after Sections 7-10b)", vec = new_status),
  list(label = "3. corrected_new_status (after Section 11 D<->G)", vec = corrected_new_status),
  list(label = "4. tree_histories       (tree-level)", vec = as.vector(tree_histories))
)

biology_log <- file.path(CHECK_folder, "biology_assessment.txt")
con <- file(biology_log, open = "w")
writeLines(c(
  paste0("# biology_assessment.txt  - generated ", format(Sys.time())),
  "# allowed transitions: AA AD AG DD DG GD GG PP PA",
  "# GD (gone->dead) and DG (dead->gone) are valid after Section 10.4 correction",
  "# illegal: every other 2-letter substring",
  ""
), con)

for (v in versions) {
  rep_v <- assess_biology(v$vec, v$label, valid_trans)
  writeLines(c(
    strrep("-", 72),
    rep_v$label,
    sprintf("  n_stems              : %d", rep_v$n_stems),
    sprintf("  n_rows_with_illegal  : %d", rep_v$n_bad_rows),
    sprintf(
      "  n_rows_with_A[DG]+A  : %d  (D/G stranded between two A)",
      rep_v$n_stranded
    ),
    "  illegal transition counts:"
  ), con)
  if (length(rep_v$illegal_tab) > 0L) {
    for (k in seq_along(rep_v$illegal_tab)) {
      writeLines(sprintf(
        "      %s : %d",
        names(rep_v$illegal_tab)[k],
        rep_v$illegal_tab[k]
      ), con)
    }
  } else {
    writeLines("      (none)", con)
  }
  writeLines("  legal transition counts:", con)
  for (k in seq_along(rep_v$legal_tab)) {
    writeLines(sprintf(
      "      %s : %d",
      names(rep_v$legal_tab)[k],
      rep_v$legal_tab[k]
    ), con)
  }
  writeLines("", con)
}

writeLines(c(
  strrep("=", 72),
  "EXPECTED RESULTS",
  "  1. original_status   : MAY contain any illegal transitions (raw).",
  "  2. new_status        : DA/GA = 0, PD/PG = 0, A[DG]+A = 0.",
  "  3. corrected_new_st. : DA/GA = 0, PD/PG = 0, A[DG]+A = 0; GD and DG now valid (per-census correction, Section 10.4).",
  "  4. tree_histories    : DA/GA = 0, PD/PG = 0, A[DG]+A = 0.",
  ""
), con)
close(con)

cat(sprintf("\n📄 Biology assessment written to: %s\n", biology_log))

# Print the most important lines to stdout so the user sees pass/fail
# without opening the file.
summary_msg <- vapply(versions, function(v) {
  r <- assess_biology(v$vec, v$label, valid_trans)
  sprintf(
    "  %s\n      illegal rows = %d, A[DG]+A rows = %d",
    v$label, r$n_bad_rows, r$n_stranded
  )
}, character(1))
cat("\nBiology assessment summary:\n")
cat(paste(summary_msg, collapse = "\n"), "\n")
rm(versions, summary_msg)

# ========================================================================
# SECTION 12: DATA QUALITY REPORTS AND DIAGNOSTICS
# ========================================================================
# Export CSV diagnostics that document key edge cases and status changes.
# Covers resurrection patterns, never-alive stems, gone-then-alive stems,
# and a summary of how raw statuses were transformed through the pipeline.

cat("\n📋 Creating diagnostic reports...\n")

# Optional view for manual inspection (commented out for batch processing)

## DIAGNOSTIC 1: Export resurrected trees ####
# Identify and export cases where trees were resurrected (dead → alive)
# This helps QA/QC the data to find potential measurement errors

# Find all stems belonging to trees that showed D→G conversion
# (stems where tree was alive but stem was marked dead)

# Build wide DBH/Status/ListOfTSM via rbindlist + dcast:
DT_long <- rbindlist(
  lapply(seq_along(ViewFullTable_split), function(i) {
    ViewFullTable_split[[i]][, .(TreeID, StemID, DBH, Status, ListOfTSM, census = i)]
  }),
  use.names = TRUE, fill = TRUE
)

DBH_wide <- dcast(DT_long,
  formula = TreeID + StemID ~ census,
  value.var = "DBH"
)

Status_wide <- dcast(DT_long,
  formula = TreeID + StemID ~ census,
  value.var = "Status"
)

TSM_wide <- dcast(DT_long,
  formula = TreeID + StemID ~ census,
  value.var = "ListOfTSM"
)

problem_dt <- DBH_wide[Status_wide, on = .(TreeID, StemID)][TSM_wide, on = .(TreeID, StemID)]

# Find reordering indices: match master stem order (unique_StemID) to current rows
idx <- match(
  unique_StemID$StemID, # desired order (from master list)
  problem_dt$StemID
) # current StemID order in wide table

# Reorder rows to align with unique_StemID
problem_dt <- problem_dt[idx] # subset/reorder using idx
problem_dt[, `:=`(
  original_status = original_status,
  new_status = new_status,
  corrected_new_status = corrected_new_status,
  tree_histories = tree_histories
)]

problem_df <- problem_dt

# Filter to only resurrected stems
problem.ID <- unique_StemID$StemID[grepl("DA", original_status)]
problem_df <- problem_df[StemID %in% problem.ID, ]

# Create descriptive column names (DBH_1, DBH_2, Status_1, Status_2, etc.)
names(problem_df)[-c(1, 2, ncol(problem_df) - 3, ncol(problem_df) - 2, ncol(problem_df) - 1, ncol(problem_df))] <-
  paste(rep(c("DBH", "Status", "ListOfTSM"), each = length(ViewFullTable_split)), seq_along(ViewFullTable_split), sep = "_")

# Export to CSV for review
fwrite(problem_df, file = file.path(CHECK_folder, "stem_status_resurected.csv"))
cat("  ✓ Saved: stem_status_resurected.csv\n")

## DIAGNOSTIC 2: Export status transformation summary ####
# Create a summary table showing all unique status transformation patterns
# Useful for QA/QC and understanding what the script changed

# Combine all status versions into unique rows
X_dt <- data.table(original_status, new_status, corrected_new_status, tree_histories)
# rename VA with tree_histories
# set names explicity old to new
setnames(X_dt, c("original_status", "new_status", "corrected_new_status", "V1"), c("original_status", "new_status", "corrected_new_status", "tree_histories"))

# check <- X_dt[, count := .N, by = .(tree_histories)]
# setorder(check, tree_histories)
# head(check, 100)

X <- unique(X_dt)
setorder(X, original_status)
fwrite(X, file = file.path(CHECK_folder, "status_changed.csv"))
cat("  ✓ Saved: status_changed.csv\n")

## DIAGNOSTIC 3: Flag stems that never showed alive status ####
# Find cases where P goes directly to G or D, without first going to A
# These are stems that existed before monitoring but never showed alive status
# Could indicate: 1) stems that died before first census, 2) data quality issues

problem.ID <- unique_StemID$StemID[
  grepl("PG|PD|PM|P$", new_status)
]

cat(sprintf("  Stems never alive (P→P/D/G): %d\n", length(problem.ID)))

fwrite(as.data.table(ViewFullTable)[StemID %in% problem.ID],
  file = file.path(CHECK_folder, "subset_stems_never_alive.csv")
)
cat("  ✓ Saved: subset_stems_never_alive.csv\n")

## DIAGNOSTIC 4: Flag stems that reappeared after being gone ####
# Find cases where G goes to A (gone then alive)
# These are stems that were gone/lost but then reappeared as alive
# Important to review as they may indicate:
#   - Tagging errors (different stem given same TreeID)
#   - Data entry errors
#   - Misidentification

problem.ID <- unique_StemID$StemID[grepl("GA", original_status)]

cat(sprintf("  Stems with GA (gone→alive): %d\n", length(problem.ID)))

# Filter to only stems with GA pattern
problem_df <- problem_dt[StemID %in% problem.ID, ]

# Create descriptive column names
names(problem_df)[-c(1, 2, ncol(problem_df) - 3, ncol(problem_df) - 2, ncol(problem_df) - 1, ncol(problem_df))] <-
  paste(rep(c("DBH", "Status", "ListOfTSM"), each = length(ViewFullTable_split)), seq_along(ViewFullTable_split), sep = "_")

fwrite(problem_df, file = file.path(CHECK_folder, "stems_G_then_A.csv"))
cat("  ✓ Saved: stems_G_then_A.csv\n")

## DIAGNOSTIC 5: Export unique status transformation examples ####
# Create a reference file showing ONE EXAMPLE of each unique status transformation pattern
# Includes count of how many stems followed each pattern
# Very useful for understanding and documenting the status cleaning logic

cat("\n📊 Creating status transformation reference file...\n")

# Keep only one row per unique combination of status transformations
# This creates a reference showing each possible transformation pattern
problem_df <- problem_dt[!duplicated(problem_dt[, c("original_status", "new_status", "corrected_new_status")]), ]

# Create descriptive column names
names(problem_df)[-c(1, 2, ncol(problem_df) - 3, ncol(problem_df) - 2, ncol(problem_df) - 1, ncol(problem_df))] <-
  paste(rep(c("DBH", "Status", "ListOfTSM"), each = length(ViewFullTable_split)), seq_along(ViewFullTable_split), sep = "_")

# Add column showing how many stems have each status pattern
# This gives context about which patterns are common vs rare

DT_counts <- data.table(original_status, new_status)[
  , .N,
  by = .(original_status, new_status, corrected_new_status)
]

problem_df <- merge(problem_df[, tree_histories := NULL], DT_counts,
  by = c("original_status", "new_status", "corrected_new_status"),
  all.x = TRUE
)
setnames(problem_df, "N", "number_of_cases")

fwrite(problem_df, file = file.path(CHECK_folder, paste0(site, "_examples_of_what_happens_for_status.csv")))
cat(sprintf("  ✓ Saved: %s_examples_of_what_happens_for_status.csv\n", site))

## DIAGNOSTIC 6: Export all combinations

problem_df <- problem_dt
# Create descriptive column names (DBH_1, DBH_2, Status_1, Status_2, etc.)
names(problem_df)[-c(1, 2, ncol(problem_df) - 3, ncol(problem_df) - 2, ncol(problem_df) - 1, ncol(problem_df))] <-
  paste(rep(c("DBH", "Status", "ListOfTSM"), each = length(ViewFullTable_split)), seq_along(ViewFullTable_split), sep = "_")

# if DBH_1 is not NA, then, 1
problem_df[, dbh1 := ifelse(!is.na(get("DBH_1")), 1, 0)]
problem_df[, dbh2 := ifelse(!is.na(get("DBH_2")), 1, 0)]
problem_df[, dbh3 := ifelse(!is.na(get("DBH_3")), 1, 0)]
problem_df[, dbh4 := ifelse(!is.na(get("DBH_4")), 1, 0)]
problem_df[, dbh5 := ifelse(!is.na(get("DBH_5")), 1, 0)]
problem_df[, dbh6 := ifelse(!is.na(get("DBH_6")), 1, 0)]
problem_df[, dbh7 := ifelse(!is.na(get("DBH_7")), 1, 0)]
# problem_df[, dbh8 := ifelse(!is.na(get("DBH_8")), 1, 0)]
# problem_df[, dbh9 := ifelse(!is.na(get("DBH_9")), 1, 0)]

write.csv(unique(problem_df[, .(dbh1, dbh2, dbh3, dbh4, original_status, new_status, corrected_new_status, tree_histories)]),
  file = file.path(CHECK_folder, paste0(site, "_all_combinations_of_dbh_and_status.csv")),
  row.names = FALSE
)

# ========================================================================
# SECTION 13: STATUS × DBH SUPPORT SUMMARY
# ========================================================================
# Audit the final corrected_new_status_matrix before export: verify that
# no illegal transitions remain, that A/P cells are preserved, and
# summarise status codes by whether a DBH measurement exists.

stopifnot(
  is.matrix(corrected_new_status_matrix),
  nrow(corrected_new_status_matrix) == nrow(DBHs),
  ncol(corrected_new_status_matrix) == ncol(DBHs)
)

# ---- (a) and (b) final biology + A/P-preservation audit ------------------
final_pairs <- substring(
  rep(corrected_new_status, each = ncol(corrected_new_status_matrix) - 1L),
  rep(seq_len(ncol(corrected_new_status_matrix) - 1L), length(corrected_new_status)),
  rep(seq_len(ncol(corrected_new_status_matrix) - 1L), length(corrected_new_status)) + 1L
)
illegal <- setdiff(unique(final_pairs), valid_trans)
if (length(illegal) > 0L) {
  warning(sprintf(
    "corrected_new_status_matrix still contains illegal transitions: %s. See biology_assessment.txt for details.",
    paste(illegal, collapse = ", ")
  ))
} else {
  cat("✓ corrected_new_status_matrix passes the biology check (no illegal transitions).\n")
}

stranded <- sum(grepl("A[DG]+A", corrected_new_status))
if (stranded > 0L) {
  warning(sprintf(
    "corrected_new_status still has %d stem(s) with D/G stranded between two A's.",
    stranded
  ))
} else {
  cat("✓ No D/G stranded between two A's in corrected_new_status_matrix.\n")
}

# A/P cells must match between new_status_matrix and corrected_new_status_matrix.
AP_in <- new_status_matrix %in% c("A", "P")
AP_out <- corrected_new_status_matrix %in% c("A", "P")
stopifnot(identical(AP_in, AP_out))
stopifnot(identical(
  new_status_matrix[AP_in],
  corrected_new_status_matrix[AP_in]
))
cat("✓ A and P cells were preserved between new_status_matrix and corrected_new_status_matrix.\n")

rm(final_pairs, illegal, stranded, AP_in, AP_out)

# ---- (c) status × DBH × census summary -----------------------------------
n_cens <- ncol(corrected_new_status_matrix)

# Long-format cell-level table (one row per stem×census).
cells_dt <- data.table(
  status  = as.vector(corrected_new_status_matrix),
  has_DBH = !is.na(as.vector(DBHs)),
  census  = rep(seq_len(n_cens), each = nrow(corrected_new_status_matrix))
)

# Per-census × status × has_DBH counts.
summary_long <- cells_dt[
  , .(n_cells = .N),
  by = .(census, status, has_DBH)
][order(census, status, has_DBH)]

# Overall (all censuses pooled).
summary_overall <- cells_dt[
  , .(n_cells = .N),
  by = .(status, has_DBH)
][order(status, has_DBH)]
summary_overall[, pct := round(100 * n_cells / sum(n_cells), 3)]

# Wide pivot for compact viewing: rows = status, cols = (has_DBH, census).
summary_wide <- dcast(
  summary_long,
  status + has_DBH ~ census,
  value.var = "n_cells",
  fill = 0L
)

# Save CSVs.
fwrite(summary_long,
  file = file.path(CHECK_folder, "status_x_dbh_summary_long.csv")
)
fwrite(summary_overall,
  file = file.path(CHECK_folder, "status_x_dbh_summary_overall.csv")
)
fwrite(summary_wide,
  file = file.path(CHECK_folder, "status_x_dbh_summary_wide.csv")
)

# Human-readable plain-text version.
txt_path <- file.path(CHECK_folder, "status_x_dbh_summary.txt")
con <- file(txt_path, open = "w")
writeLines(c(
  paste0("# status_x_dbh_summary.txt - generated ", format(Sys.time())),
  "# canonical post-pipeline summary of corrected_new_status_matrix",
  sprintf(
    "# total cells = %d stems × %d censuses = %d",
    nrow(corrected_new_status_matrix), n_cens,
    nrow(corrected_new_status_matrix) * n_cens
  ),
  "",
  "## OVERALL (all censuses pooled)",
  "## status  has_DBH       n_cells       pct"
), con)
for (i in seq_len(nrow(summary_overall))) {
  writeLines(sprintf(
    "   %-6s %-5s %12d  %8.3f%%",
    summary_overall$status[i],
    as.character(summary_overall$has_DBH[i]),
    summary_overall$n_cells[i],
    summary_overall$pct[i]
  ), con)
}
writeLines(c(
  "",
  "## EXPECTATIONS (red flags if violated)",
  "   - status = 'P' & has_DBH = TRUE   should be 0  (P means not yet recruited)",
  "   - status = 'A' & has_DBH = FALSE  is OK but downstream MUST interpolate DBH",
  "   - status = 'D' & has_DBH = TRUE   are 'broken-below-with-measurement' or last-alive measurements",
  "   - status = 'G' & has_DBH = TRUE   same biological meaning, multi-stem",
  "",
  "## PER-CENSUS (wide pivot)"
), con)
write.table(summary_wide, con,
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

# Red-flag summary block.
n_P_with_DBH <- summary_overall[
  status == "P" & has_DBH == TRUE,
  sum(n_cells)
]
n_A_no_DBH <- summary_overall[
  status == "A" & has_DBH == FALSE,
  sum(n_cells)
]
n_D_with_DBH <- summary_overall[
  status == "D" & has_DBH == TRUE,
  sum(n_cells)
]
n_G_with_DBH <- summary_overall[
  status == "G" & has_DBH == TRUE,
  sum(n_cells)
]
n_other <- summary_overall[
  !(status %in% c("A", "D", "G", "P")),
  sum(n_cells)
]

writeLines(c(
  "",
  "## RED-FLAG COUNTS",
  sprintf(
    "   P with DBH (must be 0)            : %d",
    if (length(n_P_with_DBH)) n_P_with_DBH else 0L
  ),
  sprintf(
    "   A without DBH (interpolate)       : %d",
    if (length(n_A_no_DBH)) n_A_no_DBH else 0L
  ),
  sprintf(
    "   D with DBH (broken/last-measure)  : %d",
    if (length(n_D_with_DBH)) n_D_with_DBH else 0L
  ),
  sprintf(
    "   G with DBH (broken/last-measure)  : %d",
    if (length(n_G_with_DBH)) n_G_with_DBH else 0L
  ),
  sprintf(
    "   cells with status not in {A,D,G,P}: %d",
    if (length(n_other)) n_other else 0L
  )
), con)
close(con)

cat(sprintf(
  "\n📄 status x DBH summary written to:\n  %s\n  %s\n  %s\n  %s\n",
  txt_path,
  file.path(CHECK_folder, "status_x_dbh_summary_long.csv"),
  file.path(CHECK_folder, "status_x_dbh_summary_overall.csv"),
  file.path(CHECK_folder, "status_x_dbh_summary_wide.csv")
))

# Echo the headline numbers to the console.
cat("\nStatus × DBH overall summary:\n")
print(summary_overall)
if (length(n_P_with_DBH) && n_P_with_DBH > 0L) {
  warning(sprintf(
    "❌ %d 'P' cells have a DBH measurement. This is biologically impossible.",
    n_P_with_DBH
  ))
}

rm(
  cells_dt, summary_long, summary_wide,
  n_P_with_DBH, n_A_no_DBH, n_D_with_DBH, n_G_with_DBH, n_other,
  con, txt_path
)

# ========================================================================
# SECTION 14: EXPORT CENSUS TABLES
# ========================================================================
# Rename and subset each census table to ForestGEO R-table format, set
# the legacy DFstatus field for prior stems, then save each census as
# both a .Rdata object and a .csv file.

for (census in seq_along(ViewFullTable_split)) {
  ViewFullTable_split[[census]]$new_status <- corrected_new_status_matrix[, census]
}

# ========================================================================
# COORDINATE IMPUTATION: fill NA PX / PY across censuses
# ========================================================================
# Because the unique-cell selection can leave PX / PY as NA in some
# censuses for a given tree, we impute missing coordinates using
# observed values of the same tree (matched by TreeID + Tag).
#
# IMPUTATION RULE:
#   - Exactly 1 unique non-NA value for a tree → fill all NAs for that tree.
#   - More than 1 distinct non-NA value       → leave NAs untouched.
#     (Trees with conflicting coordinates are left for manual review.)
#   - No non-NA value at all                  → remains NA.
#
# The `multi_val` argument (mean / median / mode) is available for
# numeric columns ONLY when `strict = FALSE`.  The default `strict = TRUE`
# follows the simpler rule above and works for both numeric and character
# columns (e.g. QuadratName).
#
# ARGUMENTS:
#   split_list  — list of data.tables, one per census (ViewFullTable_split)
#   id_cols     — columns that uniquely identify a tree
#   coord_cols  — location columns to impute (numeric or character)
#   strict      — if TRUE (default): fill only when a single unique value
#                 exists; if FALSE: aggregate multi-value trees with multi_val
#   multi_val   — aggregation rule used only when strict = FALSE:
#                 "mean" | "median" | "mode" (numeric columns only)
#
# RETURNS:
#   The same list with NAs in coord_cols filled in-place.
# ========================================================================

coord_mode <- function(x) {
  # Statistical mode for a vector of any type (NA-aware).
  # For ties, returns the first value in order of first appearance.
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

impute_tree_coords <- function(split_list,
                               id_cols = c("TreeID", "Tag"),
                               coord_cols = c("PX", "PY", "QuadratName"),
                               strict = TRUE,
                               multi_val = c("mode", "median", "mean")) {
  # multi_val controls aggregation when strict = FALSE AND a tree has more
  # than one distinct non-NA value:
  #   numeric columns  → mean / median / mode (as chosen)
  #   character columns → always mode (only sensible aggregate for strings)
  multi_val <- match.arg(multi_val)
  num_agg_fun <- switch(multi_val,
    mean   = function(x) mean(as.numeric(x), na.rm = TRUE),
    median = function(x) median(as.numeric(x), na.rm = TRUE),
    mode   = coord_mode
  )
  # 1. Pool all location observations across every census
  pool <- rbindlist(
    lapply(split_list, function(dt) dt[, c(id_cols, coord_cols), with = FALSE]),
    use.names = TRUE, fill = TRUE
  )
  # 2. Compute one reference value per tree × column.
  #    Character columns always fall back to mode when strict = FALSE.
  ref <- pool[, lapply(.SD, function(col) {
    vals <- col[!is.na(col)]
    n_uniq <- length(unique(vals))
    if (n_uniq == 0L) {
      return(if (is.character(col) || is.factor(col)) NA_character_ else NA_real_)
    } # no data at all → stay NA
    if (n_uniq == 1L) {
      return(vals[[1L]])
    } # single value → use it
    # More than 1 distinct value:
    if (strict) {
      return(if (is.character(col) || is.factor(col)) NA_character_ else NA_real_)
    } # strict mode → leave NAs alone
    # Non-strict: aggregate by column type
    if (is.character(col) || is.factor(col)) coord_mode(vals) else num_agg_fun(vals)
  }), .SDcols = coord_cols, by = id_cols]
  # Diagnostic: summarise how many trees have ambiguous coordinates
  pool_uniq <- pool[, lapply(.SD, function(col) {
    length(unique(col[!is.na(col)]))
  }), .SDcols = coord_cols, by = id_cols]
  for (cc in coord_cols) {
    n0 <- pool_uniq[get(cc) == 0L, .N]
    n1 <- pool_uniq[get(cc) == 1L, .N]
    nm <- pool_uniq[get(cc) > 1L, .N]
    cat(sprintf(
      "  [impute_tree_coords] %-15s | no data: %d | single value (fill): %d | multiple values (%s): %d\n",
      cc, n0, n1,
      if (strict) "leave NA" else if (is.character(pool[[cc]]) || is.factor(pool[[cc]])) "mode" else multi_val,
      nm
    ))
  }
  # 3. Fill NAs in each census using the reference table
  split_list <- lapply(seq_along(split_list), function(i) {
    dt <- copy(split_list[[i]])
    orig_cols <- names(dt)
    n_na_before <- sum(vapply(coord_cols, function(cc) sum(is.na(dt[[cc]])), integer(1L)))
    merged <- merge(dt, ref, by = id_cols, suffixes = c("", ".ref"), all.x = TRUE)
    for (cc in coord_cols) {
      ref_col <- paste0(cc, ".ref")
      if (ref_col %in% names(merged)) {
        merged[is.na(get(cc)), (cc) := get(ref_col)]
        set(merged, j = ref_col, value = NULL)
      }
    }
    setcolorder(merged, orig_cols)
    n_na_after <- sum(vapply(coord_cols, function(cc) sum(is.na(merged[[cc]])), integer(1L)))
    cat(sprintf(
      "  [impute_tree_coords] census %d: filled %d NA location cell(s); %d remain (ambiguous trees)\n",
      i, n_na_before - n_na_after, n_na_after
    ))
    merged
  })
  split_list
}

cat("🗺️  Imputing missing PX / PY / QuadratName across censuses...\n")

ViewFullTable_split <- impute_tree_coords(
  split_list = ViewFullTable_split,
  id_cols = c("TreeID", "Tag"),
  coord_cols = c("PX", "PY", "QuadratName"),
  strict = FALSE,
  multi_val = "mode"
)
cat("✓ Location imputation complete.\n\n")

#   [impute_tree_coords] PX              | no data: 0 | single value (fill): 260286 | multiple values (mode): 3808
#   [impute_tree_coords] PY              | no data: 0 | single value (fill): 260900 | multiple values (mode): 3194
#   [impute_tree_coords] QuadratName     | no data: 0 | single value (fill): 263946 | multiple values (mode): 148
#   [impute_tree_coords] census 1: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 2: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 3: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 4: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 5: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 6: filled 0 NA location cell(s); 0 remain (ambiguous trees)
#   [impute_tree_coords] census 7: filled 0 NA location cell(s); 0 remain (ambiguous trees)
# ✓ Location imputation complete.

cat("💾 Exporting census tables to .Rdata files...\n")

## Format and export each census as ForestGEO R table ####
# Loop through each census and export in standardized format
check_data <- rbindlist(lapply(ViewFullTable_split, function(dt) dt[, ..ViewFullTable_columns_to_keep]))
# check unique px per TreeID and unique py per TreeID
check_coordinates <- check_data[, .(TreeID, PX, PY, QuadratName)]
check_coordinates[, c("n_px", "n_py", "n_quadrat") :=
  .(uniqueN(PX), uniqueN(PY), uniqueN(QuadratName)),
by = TreeID
]
check_coordinates <- unique(check_coordinates[
  n_px > 1L |
    n_py > 1L |
    n_quadrat > 1L,
  .(TreeID, n_px, n_py, n_quadrat)
])

check_coordinates

fwrite(check_coordinates, file.path(CHECK_folder, "repeated_coordinates.csv"))

impute_tree_dates <- function(split_list,
                              id_cols = c("TreeID", "Tag"),
                              date_col = "ExactDate",
                              quadrat_col = "QuadratName",
                              strict = TRUE) {
  # Process each census independently
  split_list <- lapply(seq_along(split_list), function(i) {
    dt <- copy(split_list[[i]])
    n_na_before <- sum(is.na(dt[[date_col]]))

    if (n_na_before == 0L) {
      cat(sprintf("  [impute_tree_dates] census %d: no missing dates\n", i))
      return(dt)
    }

    # Step 1: Within-tree imputation (multiple observations per tree in same census)
    tree_ref <- dt[!is.na(get(date_col)),
      .(ref_date = {
        vals <- get(date_col)
        n_uniq <- length(unique(vals))
        if (n_uniq == 1L) vals[1L] else date_mode(vals)
      }),
      by = id_cols
    ]

    dt <- merge(dt, tree_ref, by = id_cols, all.x = TRUE, suffixes = c("", ".tree_ref"))
    dt[is.na(get(date_col)), (date_col) := ref_date]
    dt[, ref_date := NULL]

    n_na_after_tree <- sum(is.na(dt[[date_col]]))
    n_filled_tree <- n_na_before - n_na_after_tree

    if (!strict && n_na_after_tree > 0L) {
      # Step 2: Quadrat-level imputation
      quadrat_ref <- dt[!is.na(get(date_col)),
        .(ref_date = {
          vals <- get(date_col)
          n_uniq <- length(unique(vals))
          if (n_uniq == 1L) vals[1L] else date_mode(vals)
        }),
        by = quadrat_col
      ]

      dt <- merge(dt, quadrat_ref,
        by = quadrat_col, all.x = TRUE,
        suffixes = c("", ".quad_ref")
      )
      dt[is.na(get(date_col)), (date_col) := ref_date]
      dt[, ref_date := NULL]

      n_na_final <- sum(is.na(dt[[date_col]]))
      n_filled_quadrat <- n_na_after_tree - n_na_final

      cat(sprintf(
        "  [impute_tree_dates] census %d: filled %d from tree-level, %d from quadrat-level; %d remain NA\n",
        i, n_filled_tree, n_filled_quadrat, n_na_final
      ))
    } else {
      cat(sprintf(
        "  [impute_tree_dates] census %d: filled %d from tree-level%s; %d remain NA\n",
        i, n_filled_tree,
        if (strict) " (strict mode: no quadrat imputation)" else "",
        n_na_after_tree
      ))
    }

    dt
  })

  split_list
}

# Mode function for dates
date_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ViewFullTable_split <- impute_tree_dates(
  split_list = ViewFullTable_split,
  id_cols = c("TreeID", "Tag"),
  date_col = "ExactDate",
  quadrat_col = "QuadratName",
  strict = FALSE
)

# [impute_tree_dates] census 1: filled 18551 from tree-level, 228199 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 2: filled 22315 from tree-level, 203986 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 3: filled 35348 from tree-level, 177833 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 4: filled 27171 from tree-level, 152580 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 5: filled 12161 from tree-level, 74307 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 6: filled 34293 from tree-level, 115633 from quadrat-level; 0 remain NA
# [impute_tree_dates] census 7: filled 31451 from tree-level, 101777 from quadrat-level; 0 remain NA

cat("✓ Dates imputation complete.\n\n")

cat("💾 Exporting census tables to .Rdata files...\n")

## Format and export each census as ForestGEO R table ####
# Loop through each census and export in standardized format
check_data <- rbindlist(lapply(ViewFullTable_split, function(dt) dt[, ..ViewFullTable_columns_to_keep]))
setorder(check_data, TreeID, StemID, CensusID)
# check unique date per TreeID
check_dates <- unique(check_data[, .(TreeID, ExactDate, CensusID)])

# get nunique exactdate per treeid and census
check_dates[, c("n_dates") :=
  .(uniqueN(ExactDate)),
by = .(CensusID, TreeID)
]

fwrite(check_dates[n_dates > 1L, .(TreeID, CensusID, n_dates)], file.path(CHECK_folder, "repeated_dates.csv"))

site <- "hkk"

for (census in seq_along(ViewFullTable_split)) {
  cat(sprintf("  Processing census %d...\n", census))
  # Extract current census data and ensure it's a data.table
  X <- as.data.table(ViewFullTable_split[[census]])
  # NOTE: SPLIT TREEID
  X[, TreeID := stringr::str_split_fixed(TreeID, "_", 2)[, 1]]
  setorder(X, Tag, StemID, CensusID)
  # Validate that every required column is present BEFORE subsetting,
  # so a missing column produces a clear error rather than a cryptic
  # data.table failure deep inside the export loop.
  missing_cols <- setdiff(ViewFullTable_columns_to_keep, names(X))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Census %d is missing required columns: %s",
      census, paste(missing_cols, collapse = ", ")
    ))
  }
  # Select only the columns needed for ForestGEO R format
  # ..ViewFullTable_columns_to_keep references the variable defined in Section 1
  X <- X[, ..ViewFullTable_columns_to_keep]
  # Rename columns from database format to ForestGEO R format
  # Example: TreeID → treeID, Mnemonic → sp, PX → gx, etc.
  setnames(X, old = ViewFullTable_columns_to_keep, new = new_names_columns_to_keep)
  # Set DFstatus field for stems with "prior" status
  # DFstatus is a legacy field used in older ForestGEO analyses
  # Only stems with status='P' get DFstatus="prior"
  X[Rstatus == "P", DFstatus := "prior"]
  # Convert to data.frame for compatibility with legacy R code
  # Many ForestGEO functions expect data.frame, not data.table
  fwrite(X, file = file.path(OUTPUT_folder, sprintf("%s.stem%d.csv", site, census)))
  X <- as.data.frame(X)
  print(head(X))
  # Create R object with standardized name: [site].stem[census#]
  # Example: "hkk.stem1", "hkk.stem2", etc.
  assign(paste0(site, ".stem", census), X)
  # Export as .Rdata file to OUTPUT_folder
  # This creates files like: hkk.stem1.Rdata, hkk.stem2.Rdata, etc.
  save(
    list = paste0(site, ".stem", census),
    file = file.path(OUTPUT_folder, paste0(site, ".stem", census, ".Rdata"))
  )
  cat(sprintf("    ✓ Saved: %s.stem%d.Rdata (%d rows)\n", site, census, nrow(X)))
}

cat(sprintf("\n✓✓ All %d census tables exported successfully\n\n", length(ViewFullTable_split)))

# ========================================================================
# SECTION 15: EXPORT SPECIES TABLE
# ========================================================================
# Build a ForestGEO species taxonomy table from ViewTaxonomy, keep only
# species observed in this plot, clean NULL values, and export it as
# .rdata and .csv.

# --------------------------------------------------------------------
# --------------------------------------------------------------------
ViewTaxonomy <- fread(
  file = file.path(main_path, "DATA", "ViewTaxonomy_hkk.csv"),
  sep = "\t",
  stringsAsFactors = FALSE, # Keep text as character strings
  na.strings = c("NA", "NULL", "")
)

# Export unique values for QA/QC
show_levels(ViewTaxonomy, n_to_print = Inf, output = "print")
fwrite(show_levels(ViewTaxonomy, n_to_print = Inf, output = "df"),
  file = file.path(CHECK_folder, paste0("ViewTaxonomy_levels_", site, ".csv")),
  sep = ",",
  na = "",
  row.names = FALSE
)

cat("📚 Creating and exporting species table...\n")

## Transform ViewTaxonomy into ForestGEO species table format ####

# Create independent copy to avoid modifying original ViewTaxonomy
sptable <- copy(ViewTaxonomy)

# Sort alphabetically by species code (Mnemonic)
# This makes the table easier to browse and reference
setorder(sptable, Mnemonic)

# Create full Latin name by combining Genus and SpeciesName
# Example: Genus="Acer" + SpeciesName="rubrum" → Latin="Acer rubrum"
sptable[, Latin := paste(Genus, SpeciesName)]

# Select columns needed for ForestGEO format
cols_to_keep <- c(
  "Mnemonic", # Species code (e.g., "ACERUB")
  "Latin", # Full scientific name
  "Genus", # Genus name
  "SpeciesName", # Species epithet
  "Family", # Taxonomic family
  "SpeciesID", # Numeric species identifier
  "Authority", # Taxonomic authority
  "IDLevel", # Identification confidence level
  "ListOfOldNames", # Historical synonyms
  "Subspecies", # Subspecies information
  "wsg"
)
sptable <- sptable[, ..cols_to_keep]

# Rename columns to ForestGEO standard names
# Database names → ForestGEO R names
setnames(sptable,
  old = c(
    "Mnemonic", "Latin", "Genus", "SpeciesName", "Family", "SpeciesID",
    "Authority", "IDLevel", "ListOfOldNames", "Subspecies"
  ),
  new = c(
    "sp", # Species code
    "Latin", # Full Latin name (unchanged)
    "Genus", # Genus name (unchanged)
    "Species", # Species epithet (was SpeciesName)
    "Family", # Family (unchanged)
    "SpeciesID", # Species ID (unchanged)
    "Authority", # Authority (unchanged)
    "IDLevel", # ID level (unchanged)
    "syn", # Synonyms (was ListOfOldNames)
    "subsp" # Subspecies (was Subspecies)
  )
)

# Filter to only species that actually appear in the census data
# This removes species from the taxonomy list that aren't in this plot
species_in_census <- unique(ViewFullTable$Mnemonic)
sptable <- sptable[sp %in% species_in_census]

cat(sprintf("  Species in this plot: %d\n", nrow(sptable)))

# Clean up "NULL" and empty strings by converting to NA
# Some databases export NULL values as the string "NULL"
for (col in names(sptable)) {
  sptable[get(col) == "NULL", (col) := NA]
  sptable[get(col) == "", (col) := NA]
}

# Convert to data.frame for compatibility with legacy ForestGEO code
sptable <- as.data.frame(sptable)

# Export species table as .Rdata file
obj_name <- paste0(site, ".spptable")
assign(obj_name, sptable)
save(list = obj_name, file = file.path(OUTPUT_folder, paste0(obj_name, ".rdata")))

fwrite(sptable,
  file = file.path(OUTPUT_folder, sprintf("%s.spptable.csv", site))
)

cat(sprintf("  ✓ Saved: %s.spptable.rdata\n\n", site))
