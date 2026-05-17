################################################################################
# FORESTGEO STEM TABLE CREATION SCRIPT
################################################################################
# Purpose: Transform ViewFullTable_site.csv and ViewTaxonomy_site.csv into R tables
#          ([site].stem[census#].rdata and [site].spptable.rdata)
#
# Authors:
#   - Original: Valentine Herrmann (HerrmannV@si.edu)
#   - Updated: Jose Medina (MedinaJA@si.edu)
#
# Environment:
#   - R version 4.5.2 (2025-10-31)
#   - data.table version 1.17.8
#
# Last modified: [November 17th, 2025]
################################################################################
#
# TABLE OF CONTENTS
# =================
# Section  1: Setup and Initialization
# Section  2: Load Input Data
# Section  3: Standardize Data Format Across Censuses
# Section  4: Status Code Transformation (English → Single Letter) + Encounter Histories
#            Includes the DBH-aware "broken below" rule (see Section 4 header).
# Section  5: [VALIDATION] Check Encounter History Patterns (pre-propagation QA)
# Section  6: [reserved]
# Section  7: Status Propagation Rules (N → P/A/D/G)
# Section  8: [CORRECTION] Fix Prior-to-Dead/Gone (DBH-aware)
# Section  9-10: [CORRECTION] Resolve Resurrections D→A / G→A (DBH-aware merged)
# Section 10b: [SAFETY NET] Resolve D / G stranded between two A's
# Section 11: Tree-Level Status Calculation & D/G Correction (vectorized)
# Section 12: Apply Corrected Status to Census Data
# Section 12b: [VALIDATION] Unified Biology Assessment of all 4 history versions
# Section 13: Generate Data Quality Reports
# Section 13b: [VALIDATION] Status × DBH summary + final audit before export
# Section 14: Export Census Tables (.Rdata)
# Section 15: Export Species Table (.Rdata)
#
# VARIABLE NAMING CONVENTIONS (single source of truth)
# =====================================================
#   Status               raw English status from ViewFullTable ("alive",
#                        "dead", "broken below", "missing", …). Mutated
#                        in Section 4 by the broken-below DBH-aware rule
#                        BEFORE the wide pivot.
#   original_status      per-stem encounter history string built from the
#                        post-broken-below `Status` (e.g. "NNAAAD"); never
#                        mutated again — used as the pre-propagation
#                        baseline.
#   new_status           working per-stem encounter history. Updated by:
#                          Section 7  N→P/A/D/G propagation
#                          Section 8  fix_PD_PG_inconsistencies (DBH-aware)
#                          Section 9-10 fix_resurrections (DBH-aware)
#                          Section 10b fix_DG_between_A safety net
#   corrected_new_status per-stem history after the Section 11 D↔G remap
#                        (single-stem trees: G→D; multi-stem with any
#                        live stem at last census: D→G; A and P preserved).
#                        This is the canonical exported per-stem code.
#   tree_histories       per-tree encounter history aggregated from all
#                        stems of a TreeID via tree_state() (Section 11).
#   DBHs                 numeric matrix [n_stems x n_censuses] of DBH,
#                        built once in Section 8 and reused throughout.
#   Rstatus              column name in the EXPORTED tables holding the
#                        per-census slice of corrected_new_status
#                        (renamed from "new_status" in Section 14).
#   DFstatus             legacy ForestGEO column kept for backward
#                        compatibility; carries the raw `Status` string
#                        (post broken-below resolution). Set to "prior"
#                        on rows where Rstatus == "P".
################################################################################
#
# SCRIPT OVERVIEW:
# This script processes ForestGEO census data to create standardized R data tables.
# It handles status propagation across censuses, identifies resurrected stems,
# and exports cleaned data in .Rdata format for analysis.
#
# MAIN FEATURES:
# - Uses data.table for fast processing and clean console output
# - Verbose diagnostic messages showing status transformations
# - Generates QA/QC reports for data quality review
# - Handles multi-stem trees correctly (e.g., D vs G status codes)
# - Detects and corrects invalid stem resurrection cases
#
# INPUT FILES:
# - ViewFullTable_[site].csv: Tab-delimited census data (one row per stem per census)
# - ViewTaxonomy_[site].csv: Tab-delimited species taxonomy data
#
# OUTPUT FILES:
# - [site].stem[N].Rdata: One file per census with standardized format
# - [site].spptable.rdata: Species taxonomy table
# - Multiple CSV diagnostic files in ./DATA/CHECKS/
#
# WORKFLOW (HIGH-LEVEL PIPELINE):
# 1. Load raw census (ViewFullTable) and taxonomy (ViewTaxonomy) data
# 2. Standardize census format (same stem set & order across all censuses)
# 3. Transform raw English status values to single-letter codes (A/D/G/P/N)
# 4. Validate encounter histories for biological consistency (pre-propagation QA)
# 5. Apply status propagation rules (resolve N → P / D / G and forward fill terminal states)
# 6. Remove invalid resurrection patterns (alive after dead/gone with no DBH)
# 7. Backfill validated resurrections (dead → alive with DBH evidence) as continuously alive
# 8. Correct Prior→Dead/Gone inconsistencies (replace PD / PG with PP)
# 9. Compute tree-level life histories and apply final biologically consistent D vs G mapping
#    (single-stem: D = dead tree; multi-stem: G = gone stem while tree persists)
# 10. Apply corrected stem statuses to census tables and export per-census stem files
# 11. Export species taxonomy table (spptable)
# 12. Generate diagnostic QA/QC reports (status changes, edge cases, resurrection review)
#
# STATUS CODE DEFINITIONS:
#   A = Alive (stem measured and alive)
#   D = Dead (entire tree dead – single-stem trees, or all stems dead)
#   G = Gone (stem dead/lost but tree still alive – multi-stem trees only)
#   P = Prior (prior to first observation)
#   N = No data (missing measurement - converted to A/D/G/P during processing)
######################################################

# ========================================================================
# SECTION 1: SETUP AND INITIALIZATION
# ========================================================================
# PURPOSE:
#   Clean workspace, load required libraries, and define user parameters
#
# WHAT THIS SECTION DOES:
#   1. Removes all objects from memory to start fresh
#   2. Loads data.table library for fast data processing
#   3. Sets input/output folder paths
#   4. Creates output directories if they don't exist
#   5. Defines site name and processing parameters
#   6. Defines column names for data formatting
#
# USER PARAMETERS TO CONFIGURE:
#   - INPUT_folder: Where to find ViewFullTable and ViewTaxonomy CSV files
#   - OUTPUT_folder: Where to save .Rdata output files
#   - CHECK_folder: Where to save diagnostic CSV reports
#   - site: Site code (e.g., "hkk", "bci", "scbi")
# ========================================================================

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
site <- "bci" # ⚠️ CHANGE THIS for your site (e.g., "bci", "scbi", "hkk")
cat("Processing site:", site, "\n")

# 3. Input folder containing ViewFullTable and ViewTaxonomy CSV files
INPUT_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA", "PROCESSED")
# ℹ️ Files required:
#    - ViewFullTable_[site].csv (tab-delimited)
#    - ViewTaxonomy_[site].csv (tab-delimited)

# 4. Output folder for final .Rdata files
OUTPUT_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA", "RTABLES")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
  cat("✓ Created OUTPUT folder:", OUTPUT_folder, "\n")
}

# 5. Diagnostics folder for QA/QC reports
CHECK_folder <- file.path(main_path, "BCI_stem_reconstruction", "3_PREPARE_R_TABLES", "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
  cat("✓ Created CHECKS folder:", CHECK_folder, "\n")
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
  "Date", "new_status", # date + script-computed corrected status
  "obs_row_id", "DP_PosteriorReconstructedProb" # DP reconstruction metadata
)

# Standardized ForestGEO R-table column names (must match ViewFullTable_columns_to_keep
# one-to-one: same length and same order)
new_names_columns_to_keep <- c(
  "treeID", "stemID", "tag", "StemTag", "sp", # sp = species mnemonic
  "quadrat", "gx", "gy", # gx/gy = plot coordinates in meters
  "MeasureID", "CensusID", "dbh", "hom", # hom = height of measurement
  "ExactDate", "DFstatus", "codes", # DFstatus = raw field status (legacy name)
  "date", "Rstatus", # date + script-computed corrected status
  "obs_row_id", "PosteriorProb" # DP reconstruction metadata
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
# PURPOSE:
#   Load census and taxonomy data from tab-delimited CSV files
#
# WHAT THIS SECTION DOES:
#   1. Loads ViewFullTable (census data) using fread() for speed
#   2. Displays data structure and unique values for verification
#   3. Applies special filtering if needed
#   4. Loads ViewTaxonomy (species data) using fread()
#   5. Displays taxonomy structure for verification
#
# MAIN FUNCTIONS USED:
#   - fread(): Fast data loading from data.table package
#   - show_levels(): Custom function to display unique values in each column
#
# DATA VALIDATION:
#   The show_levels() function helps verify:
#   - Status values (alive, dead, etc.)
#   - No unexpected values in key columns
#
# IMPORTANT NOTE:
#   The script treats "NA", "NULL", and "" as missing values
#   Set stringsAsFactors = FALSE to keep text as characters
# ========================================================================

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
  # PURPOSE:
  #   Show all unique values in categorical columns for data validation
  # PARAMETERS:
  #   df: Data frame or data.table to examine
  #   n_to_print: Max number of unique values to display (Inf = show all)
  #   output: "print" to print values, "df" to return data frame
  # OUTPUT:
  #   If output="print": Prints formatted text showing unique values
  #   If output="df": Returns data frame with levels as rows, columns as vars
  # -----------------------------------------------------------------------
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

ViewFullTable_RAW <- as.data.table(readRDS(file.path(
  INPUT_folder,
  "complete_dataset_with_reconstructed_stemids.rds"
)))

# ── Standardize column types ──────────────────────────────────────────────────
# The RDS file may carry incorrect types from the database export pipeline.
# Each conversion is explained below:
#   PlotName, PlotID         → factor    : categorical identifiers, not quantities
#   CensusID                 → numeric   : stored as factor in DB; must go through
#                                          as.character() first to avoid silent factor-
#                                          to-integer coercion (which would return
#                                          level indices, not the actual census numbers)
#   ExactDate                → Date      : DB stores as character string "YYYY-MM-DD"
#   Tag, StemID, StemTag,
#   TrueStemID               → character : DB exports as integer or factor; character
#                                          preserves leading zeros and avoids silent
#                                          numeric coercions during joins
#   ReconstructedStemID      → character : DP assigns integer IDs; character is required
#                                          because we later paste() these into a new
#                                          StemID string (NA values must flow as "NA")
#   QuadratName              → character : codes like "0001" must stay as strings;
#                                          integer coercion silently drops leading zeros
ViewFullTable_RAW <- transform(ViewFullTable_RAW,
  PlotName = as.factor(PlotName),
  PlotID = as.factor(PlotID),
  CensusID = as.numeric(as.character(CensusID)), # factor → character → numeric
  ExactDate = as.Date(ExactDate),
  Tag = as.character(Tag),
  TreeID = as.character(TreeID),
  ReconstructedStemID = as.character(ReconstructedStemID),
  StemID = as.character(StemID),
  StemTag = as.character(StemTag),
  TrueStemID = as.character(TrueStemID),
  QuadratName = as.character(QuadratName) # preserve leading zeros
)

# ── DIAGNOSTIC: Verify ReconstructedStemID coverage relative to StemID ───────
# Goal: confirm that DP assigned a ReconstructedStemID to EVERY row with a valid
#   DBH measurement. Dead/broken rows (no DBH) may have NA ReconstructedStemID —
#   that is expected and handled by the fill steps later.
ViewFullTable_RAW[!is.na(StemID) & single_stem_tags == FALSE, .(
  n_stems = .N,
  na_levels_StemID = sum(is.na(StemID)),
  na_levels_ReconstructedStemID = sum(is.na(ReconstructedStemID)),
  n_stems_with_StemID = sum(!is.na(StemID)),
  n_stems_with_ReconstructedStemID = sum(!is.na(ReconstructedStemID)),
  n_stems_with_both_IDs = sum(!is.na(StemID) & !is.na(ReconstructedStemID)),
  n_stems_with_neither_ID = sum(is.na(StemID) & is.na(ReconstructedStemID))
), by = CensusID][order(CensusID)]

# 1. Create a copy and apply the function
ViewFullTable <- copy(ViewFullTable_RAW)

# ── FINAL STEP: Overwrite StemID with stable TreeID_ReconstructedStemID ──────────
# WHY: The original DB StemID is NOT stable across censuses. ForestGEO can assign
#   different integer IDs to the same physical stem in different censuses (e.g.,
#   when a stem resprout-codes and is re-entered with a new StemID). The
#   ReconstructedStemID assigned by DP IS stable: it tracks the same biological
#   stem across all censuses regardless of DB StemID changes.
# CONSTRUCTION: paste(TreeID, ReconstructedStemID, sep="_") produces a globally
#   unique key combining tree identity (TreeID) with stem identity (ReconstructedStemID).
# Raw_StemID: the original DB StemID is preserved in this new column for
#   traceability and cross-checking against the raw database export.
# IMPORTANT: Rows with NA ReconstructedStemID become StemID = "TreeID_NA".
#   These are dead/broken stems for which no DP identity could be established
#   (see analysis above). They carry no DBH and do not affect size or growth analyses.
ViewFullTable[, Raw_StemID := StemID] # preserve original DB StemID

setkey(ViewFullTable, TreeID)
setkey(ViewFullTable, Tag)

## Check wether tag and treeid can be interchanged
ViewFullTable[, .(n_tags = uniqueN(Tag)), by = TreeID][order(-n_tags)]
ViewFullTable[, .(n_treeid = uniqueN(TreeID)), by = Tag][order(-n_treeid)]

ViewFullTable[is.na(Tag)]
ViewFullTable[is.na(TreeID)]

ViewFullTable[, StemID := paste(TreeID, ReconstructedStemID, sep = "_")] # new stable StemID

# ── DIAGNOSTIC: Report "TreeID_NA" StemIDs ───────────────────────────────────────
# Any row still with NA ReconstructedStemID becomes StemID = "TreeID_NA" after the
# paste() above. This is expected for dead/broken stems that DP never processed
# . We report count, status
# distribution, and census distribution so the analyst can verify the scale.
stemid_problem <- ViewFullTable[grepl("_NA$", StemID)]
cat("\nRows with 'TreeID_NA' StemID (dead/broken stems without DP identity):", nrow(stemid_problem), "\n")
if (nrow(stemid_problem) > 0L) {
  cat("\nBy Status:\n")
  print(stemid_problem[, .N, by = Status][order(-N)])
  cat("\nBy CensusID:\n")
  print(stemid_problem[, .N, by = CensusID][order(CensusID)])
}

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
# PURPOSE:
#   Create a standardized format where ALL censuses have the SAME stems
#   in the SAME order, even if some stems weren't measured in some censuses
#
# THE PROBLEM:
#   - Not all stems appear in all censuses (recruits, dead stems, missing data)
#   - Different censuses have different numbers of rows
#   - Difficult to compare status across censuses
#
# THE SOLUTION:
#   - Create a master list of ALL stems (unique_StemID)
#   - Force each census to have ALL stems in the SAME order
#   - Fill in missing stems with NA for variable data
#   - Fill in fixed attributes (e.g., ID's, taxonomy, coordinates) for all stems
#
# WHAT THIS SECTION DOES:
#   1. Identifies columns that don't change across censuses (fixed_columns)
#   2. Creates unique_StemID with one row per stem and fixed attributes
# .  3. Idendity census date ranges for verification
#   4. Splits ViewFullTable into separate data frames (one per census)
#   5. Standardizes each census to have all stems in the same order
#   6. Fills in fixed attributes and census identifiers
#   7. Verifies all censuses now have identical dimensions
#
# VERIFICATION CHECKS:
#   - All StemIDs are unique (no duplicates)
#   - All censuses have the same number of rows
# .  - Census dates are in expected order
#   - Fixed columns are identical across all censuses
#
# WHY THIS MATTERS:
#   This makes it possible to create "encounter history" strings where
#   position i in the string corresponds to census i for ALL stems
#   Example: "AAAD" = alive in censuses 1-3, dead in census 4
# ========================================================================

# --------------------------------------------------------------------
# Create master stem list with fixed attributes
# --------------------------------------------------------------------

# Apply this immediately after the StemID overwrite block
# ViewFullTable <- ViewFullTable[!grepl("_NA$", StemID)]

# To select the fixed columns, this is an iterative process. For other forestgeo
# sites, and for original bci data, the stemid is different However, the
# reconstructed algorithm change the wrong used stemids So, we need to base the
# R table reconstruction using the reconstructedstemids Iterative because we
# need to find the minimum identifiers that are fixed across censuses, and that
# can be used to reconstruct the stemid

fixed_columns <- c(
  "PlotName", "PlotID",
  "Mnemonic",
  "QuadratName", "QuadratID", # "PX", "PY", # "QX", "QY",
  "TreeID",
  "Tag",
  "StemID" # ,
  # "StemTag"
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

## VERIFICATION: Check that standardization worked ####
# All censuses should now have IDENTICAL dimensions (same # rows, same # columns)
cat("\n📊 Census dimensions (AFTER standardization):\n")
lapply(ViewFullTable_split, dim)

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
# SECTION 4: STATUS CODE TRANSFORMATION - FOCUS AT THE STEM LEVEL, THEN TREES
# ========================================================================
# PURPOSE:
#   Convert English status terms to standardized single-letter codes
#   and create "encounter history" strings for each stem across censuses
#
# WHAT THIS SECTION DOES:
#   1. Transforms English status ("alive", "dead", etc.) to codes (A, D, G, P)
#   2. Handles missing data by assigning "N" (No data) placeholder
#   3. Extracts status from each census into separate columns
#   4. Creates encounter history strings (e.g., "AAAD" = alive 3x, then dead)
#
# STATUS CODE SYSTEM:
#   A = Alive (stem was measured and alive)
#   D = Dead (stem confirmed dead)
#   G = Gone/Lost (stem missing/dead in multi-stem tree with other living stems)
#   P = Prior (stem exists in later census, but this census hadn't started yet)
#   N = No data (stem not measured, status unknown - temporary placeholder)
#
# WHY ENCOUNTER HISTORIES MATTER:
#   Pattern "AAAD" tells us: stem alive in censuses 1-3, died before census 4
#   Pattern "NAAA" tells us: stem recruited after census 1
#   Pattern "AADN" tells us: potential data quality issue (dead then missing data)
#   These patterns are used in Sections 5 and 7 for status propagation rules
#
# MAIN FUNCTIONS USED:
#   - case_when(): Transform English to codes with pattern matching
#   - lapply() + do.call(cbind): Extract status from each census into columns
#   - paste0(): Concatenate status codes into encounter history string
# ========================================================================

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
# RULE (BCI-specific, requested by data owner):
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
# NOTE on rule for plain "dead" / "stem dead":
#   These are taken at face value here (mapped to "D" further below).
#   The pipeline still treats a "D" as real UNLESS the same stem appears
#   alive ("A") in a later census, in which case fix_resurrections() in
#   Section 9 backfills the spurious deads.  No change is needed here.
# ------------------------------------------------------------------------
n_broken_total <- DT_Status[Status == "broken below", .N]
n_broken_with_dbh <- DT_Status[Status == "broken below" & !is.na(DBH), .N]
n_broken_no_dbh <- DT_Status[Status == "broken below" & is.na(DBH), .N]
DT_Status[Status == "broken below" & !is.na(DBH), Status := "alive"]
DT_Status[Status == "broken below" & is.na(DBH), Status := "dead"]
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

# NOTE: For BCI data, broken below in some cases have DBH and in other cases it doesnt.
# Step 3: "broken below" → already resolved above by DBH-aware rule
# (rows have been rewritten to "alive" or "dead" before the wide pivot,
# so this gsub is a defensive no-op — kept for documentation).
original_status <- gsub("broken below", "G", original_status)

# Step 4: Everything else → "D"
# This catches: "dead", "stem dead", "missing", and any remaining
# non-A/N/G statuses.  After Step 3, "broken below" has already been
# split into "alive"/"dead" using DBH, so the only stems that fall
# into "D" here are genuine deads (no DBH measurement implied).
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
# SECTION 5: VALIDATE ENCOUNTER HISTORIES FOR BIOLOGICAL CONSISTENCY
# ========================================================================
# PURPOSE:
#   Check encounter histories for biological impossibilities before propagation
#
# VALIDATION RULES:
#   - Length must match number of censuses
#   - Only allowed codes: A (Alive), D (Dead), N (No data), G (Gone)
#   - Cannot start with D or G (stems must be alive when first observed)
#   - Terminal states: D→D only, G→G only (no resurrection)
#   - Invalid transitions: D→A, G→A, D→N, G→N, etc.
#
# OUTPUT:
#   - Reports violations to console
#   - Exports encounter history patterns to CSV for review
# ========================================================================

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

# ========================================================================
# SECTION 7: STATUS PROPAGATION RULES
# ========================================================================
# PURPOSE:
#   Propagate status codes forward to resolve no-data (N) values
#
# PROPAGATION RULES:
#   1. ^N → P: First census no-data = Prior (not yet recruited)
#   2. PN → PP: Prior propagates until stem first measured
#   3. AN → AD: Alive + no-data (interim; final D/G resolved at tree level)
#   4. GN → GG: Gone propagates (terminal state)
#   5. DN → DD: Dead propagates (terminal state)
#
# while() loops used because patterns may need multiple passes
#   Example: "ANNN" requires 3 iterations → ADNN → ADDN → ADDD
# ========================================================================

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

# length(as.vector(new_status))

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
# SECTION 8: CORRECT PRIOR-TO-DEAD/GONE INCONSISTENCIES
# ========================================================================
# PURPOSE:
#   Detect and correct cases where stems marked "P" (Prior - not yet recruited)
#   transition directly to "D" (Dead) or "G" (Gone) without ever being alive.
#
# WHY THIS IS A VIOLATION:
#   "P" means the stem had not yet been recruited / first observed. A stem
#   cannot be "Dead" or "Gone" before it has ever been "Alive". So any
#   adjacent pair PD or PG inside an encounter history is biologically
#   inconsistent and must be resolved.
#
# TWO POSSIBLE RESOLUTIONS PER PD / PG CELL:
#   (a) The D/G cell has NO DBH measurement. The most parsimonious fix is
#       PD → PP / PG → PP: the stem was simply not yet recruited, so push
#       Prior forward by one census.
#   (b) The D/G cell DOES have a DBH measurement. This means the stem was
#       actually measured at that census, so it cannot have been Prior the
#       step before. The fix here is NOT PD → PP (that would silently drop
#       a real measurement). These cases are biologically suspicious and
#       should be flagged for manual review (and probably handled by a
#       dedicated rule like "promote the surrounding P's to A").
#
# WHAT THIS SECTION DOES:
#   1. Builds a DBH matrix aligned with the propagated status matrix.
#   2. Calls fix_PD_PG_inconsistencies() which (optionally) uses the DBH
#      matrix as evidence and writes a per-cell flag log to log_PD_PG.txt.
#   3. Reports before/after pattern frequencies.
# ========================================================================

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
# SECTION 9-10: RESOLVE RESURRECTIONS (D→A, G→A)  — DBH-AWARE
# ========================================================================
# PURPOSE
#   Resolve every D→A or G→A transition in the propagated encounter
#   histories. After Section 7's propagation, any DA or GA substring is by
#   definition biologically impossible (D and G are terminal). It must be
#   one of two things:
#
#     (1) The "A" cell is real (has a DBH measurement) → the earlier D/G
#         codes were wrong. Backfill alive ("data-error" interpretation).
#         DDAA → AAAA  (was alive all along, just unobserved)
#
#     (2) The "A" cell is a "zombie" (no DBH) → the "A" itself is wrong.
#         Demote it to D (or G if the previous code was G).
#         DDDA → DDDD  (still dead)
#
# WHY MERGE THE TWO OLD SECTIONS
#   The legacy Section 9 only handled D→A (not G→A) and only when the raw
#   English status was literally "alive", missing cases that propagated
#   through AN→AD chains. Section 10 then ran an unconditional gsub
#   "DA|GA"→"AA" that would silently backfill any zombies Section 9 missed.
#   The merged function below evaluates each cell once with full DBH
#   evidence and writes a per-cell audit log.
#
# AUDIT
#   Every D→A / G→A cell is logged to log_resurrections.txt with action:
#     backfill_to_A : cell has DBH; backfill earlier D/G to A
#     demote_to_DG  : cell has no DBH; rewrite A to D (or G)
# ========================================================================

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

# DBHs matrix was built earlier in Section 8 and is still aligned with
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

cat("\n🔍 D→A / G→A patterns AFTER correction (should be empty):\n")
tbl_sorted_stem <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted_stem[grepl("DA|GA", tbl_sorted_stem$x), ])

# ========================================================================
# SECTION 10b: SAFETY NET — D / G STRANDED BETWEEN TWO A's
# ========================================================================
# BIOLOGY
#   A stem that was observed alive (A) at census t1 and alive again at
#   census t2 (t2 > t1) MUST have been continuously alive from t1 to t2.
#   Any D or G in between is, by definition, an observation/recording
#   error. The cell may have a missing DBH — that is a measurement gap,
#   not a death — and the downstream user is expected to interpolate it.
#
# WHY fix_resurrections IS NOT ENOUGH
#   fix_resurrections() only triggers a backfill when an "A" cell HAS a
#   DBH (the cell itself proves the stem is alive). It then propagates
#   the backfill through contiguous DA/GA runs. It will NOT cross a real
#   "A" in the middle of the chain. Section 7's N→D and N→G propagation
#   can also leave a row like  AADGGA  where neither the central D nor G
#   was ever a real measurement, but they are flanked by alive cells on
#   both sides.
#
# WHAT THIS PASS DOES (single vectorized matrix scan)
#   For each stem row, find the first A position and the last A position.
#   Every cell strictly between them whose status is D or G is rewritten
#   to A and logged with its DBH (NA = no measurement, present = real).
#   A and P cells are NEVER touched. The A→D→A and A→G→A patterns are
#   logged separately.
# ========================================================================

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
# SECTION 11: TREE-LEVEL STATUS CALCULATION AND FINAL D/G CORRECTION
# ========================================================================
# PURPOSE:
#   Calculate tree-level status from stem-level data and apply the final
#   correction to distinguish between D (dead tree) and G (gone stem) based
#   on whether the tree has multiple stems
#
# THE CORE LOGIC:
#   - Single-stem trees: Dead stems are marked "D" (whole tree dead)
#   - Multi-stem trees:
#        - Dead stems are marked "G" when at least one stem is alive
#        - Dead stems are marked "D" when all stems are dead
#
# WHAT THIS SECTION DOES:
#   1. Validates Tag-TreeID consistency across the dataset
#   2. Groups stems by tree (using Tag or TreeID as tree identifier)
#   3. Simulates all possible tree life-history sequences to define valid patterns
#   4. Calculates tree-level status considering:
#      - If ANY stem is alive (A) → tree is alive
#      - If ANY stem is gone (G) → tree is alive (some stems lost)
#      - If ALL stems are prior (P) → tree is prior
#      - If ALL stems are dead (D) → tree is dead
#   5. Handles "P" (Prior) status with future-looking logic:
#      - If P + D/G combination exists, check if future censuses show alive stems
#      - If future A exists → tree was alive (backfill tree status to A)
#      - If no future A → tree is dead
#   6. Applies final D/G correction based on stem count:
#      - Multi-stem trees: Convert all D → G
#      - Single-stem trees: Convert all G → D
#
# KEY FUNCTIONS:
#   - tree_state(): Computes tree-level status from stem matrix
#   - tree_exists_check(): Validates biological plausibility
#   - compute_tree_for_row(): Wrapper to compute tree sequence
#   - correct_stem_status(): Applies final D/G correction
#
# OUTPUT:
#   - tree_histories: Tree-level encounter histories
#   - corrected_new_status: Final stem-level status codes with D/G corrected
#
# WHY THIS MATTERS:
#   This is the final step that ensures biological consistency:
#   - D means "whole tree dead" (only for single-stem trees)
#   - G means "stem gone but tree alive" (only for multi-stem trees)
#   This distinction is CRITICAL for accurate mortality calculations.
# ========================================================================

#--------------------------------------------------------------
# 11.1  Tag <-> TreeID consistency check
#--------------------------------------------------------------
# Each Tag should map to exactly one TreeID (and vice versa). Anything
# else means the master stem table has duplicate identifiers and tree
# grouping below would be wrong.
tag_treeid_dt <- unique(unique_StemID[, .(Tag, TreeID)])
tags_per_treeid <- tag_treeid_dt[, .N, by = TreeID][N > 1L]
treeids_per_tag <- tag_treeid_dt[, .N, by = Tag][N > 1L]

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
# 11.2  Group stem-level new_status by TreeID (one tree = one element)
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
all_combinations[, possible := fifelse(
  (first == "A" & second %in% c("A", "D", "G")) |
    (first == "D" & second == "D") |
    (first == "G" & second == "G") |
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
  # Section 11.4 to quantify how many trees would be flagged.
  # if (!tree_exists_check(mat)) {
  #   return(NA_character_) # return NA if tree invalid
  # }
  # Compute tree sequence
  tree_state(mat)
}

compute_tree_for_row_checking <- function(stem_strings) {
  # Strict variant of compute_tree_for_row() that DOES enforce the
  # tree_exists_check() guard. Used only by the optional diagnostic in
  # Section 11.4 (RUN_TREE_EXISTS_DIAGNOSTIC). See note in the standard
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
# 11.4  APPLY TREE-LEVEL CALCULATION TO ACTUAL DATA
# ========================================================================

# Compute tree-level encounter history for every TreeID (one pass).
tree_histories_list <- lapply(new_status_split, compute_tree_for_row)

# Optional diagnostic: re-run with the strict tree_exists_check() guard
# and report which trees would have been flagged invalid. This is OFF by
# default because the strict check returns NA for trees ending in P+D/G,
# which is a known false positive (see FIXME in tree_exists_check()).
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
# SECTION 12: APPLY CORRECTED STATUS TO CENSUS DATA FRAMES
# ========================================================================
# PURPOSE:
#   Apply the corrected status codes (new_status_G) back to the census data frames
#   so they can be exported with the cleaned status information
#
# WHAT THIS SECTION DOES:
#   1. Split new_status_G string back into individual characters (one per census)
#   2. For each census data frame, add a new column "new_status"
#   3. Populate new_status with the corrected status code for that census
#
# WHY WE CREATE "new_status" COLUMN:
#   - Preserves original "Status" column for comparison/verification
#   - Allows users to see what was changed
#
# OUTPUT:
#   Each ViewFullTable_split[[i]] now has an additional column: new_status
#   This column contains the corrected single-letter status code (A/D/G/P)
#
# MAIN FUNCTIONS USED:
#   - strsplit(): Split concatenated string into character vector
#   - sapply(): Extract character at position i for census i
#   - :=: data.table assignment operator for adding column by reference
# ========================================================================

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

correct_stem_status <- function(stem_mat) {
  #------------------------------------------------------------
  # Purpose:
  #   Correct stem status codes based on the number of stems
  #   and whether the tree is completely dead or not.
  #
  # Rules:
  #   - Valid codes: "A", "D", "G", "P"
  #   - Single-stem trees:
  #        * "G" → "D" (gone means dead)
  #        * "D" stays "D"
  #
  #   - Multi-stem trees:
  #        * If all stems are dead/gone by last census → "G" → "D"
  #        * Otherwise (some still alive) → "D" → "G"
  #
  # Checks:
  #   - Input must be a matrix or data.frame
  #   - No empty input
  #   - Only valid codes (A, D, G, P)
  #   - Warn if any unexpected or missing values appear
  #------------------------------------------------------------
  #---------------------------
  # 1. Input validation
  #---------------------------
  if (is.null(stem_mat) || length(stem_mat) == 0) {
    stop("Input stem_mat is empty or NULL.")
  }
  # Convert to matrix if needed
  if (is.data.frame(stem_mat)) {
    stem_mat <- as.matrix(stem_mat)
  } else if (!is.matrix(stem_mat)) {
    stop("Input must be a matrix or data.frame.")
  }
  # Ensure character mode
  mode(stem_mat) <- "character"
  #---------------------------
  # 2. Check for invalid or unexpected codes
  #---------------------------
  allowed_codes <- c("A", "D", "G", "P")
  # Detect invalid entries
  invalid_entries <- unique(stem_mat[!stem_mat %in% allowed_codes])
  if (length(invalid_entries) > 0) {
    stop(
      "Invalid or unexpected codes detected in stem_mat: ",
      paste(invalid_entries, collapse = ", ")
    )
  }
  #---------------------------
  # 3. Basic dimensions
  #---------------------------
  n_stems <- nrow(stem_mat)
  n_cens <- ncol(stem_mat)
  if (n_stems == 0 || n_cens == 0) {
    stop("Matrix has no rows or columns.")
  }
  #---------------------------
  # 4. Single-stem trees
  #---------------------------
  if (n_stems == 1) {
    # Convert "G" → "D" (since gone = dead for single stems)
    stem_mat[stem_mat == "G"] <- "D"
    return(stem_mat)
  }
  #---------------------------
  # 5. Multi-stem trees
  #---------------------------
  is_D <- (stem_mat == "D")
  is_G <- (stem_mat == "G")
  is_P <- (stem_mat == "P")
  # is_A <- (stem_mat == "A")
  # Determine whether all stems are non-alive (D/G/P) for each census
  all_dead_by_census <- (colSums(is_D | is_G | is_P) == n_stems)
  #---------------------------
  # 6. Apply correction rules
  #---------------------------
  if (all_dead_by_census[n_cens]) {
    # If all stems are dead/gone by last census → "G" → "D"
    stem_mat[stem_mat == "G"] <- "D"
  } else {
    # Otherwise (some alive remain) → "D" → "G"
    stem_mat[stem_mat == "D"] <- "G"
  }
  #---------------------------
  # 7. Return corrected matrix
  #---------------------------
  return(stem_mat)
}

# ------------------------------------------------------------------
# 11.6  APPLY THE D/G CORRECTION (vectorized over all stems at once)
# ------------------------------------------------------------------
# The function `correct_stem_status()` defined above operates on one tree
# at a time. Splitting/applying it across ~330k trees is slow and the
# logic is simple enough to express as a few matrix operations. The block
# below produces an identical `corrected_new_status` in one vectorized
# pass.
#
# RULES (single source of truth, matching correct_stem_status):
#   - single-stem tree (n_stems == 1)        : G -> D
#   - multi-stem  tree, no live stem at last : G -> D     (tree fully dead)
#   - multi-stem  tree, any live stem at last: D -> G     (some stems lost)
# A and P cells are NEVER touched here. This is enforced by an assert at
# the end of the block.
# ------------------------------------------------------------------
stopifnot(nrow(new_status_matrix) == nrow(unique_StemID))

n_cens <- ncol(new_status_matrix)
last_col <- new_status_matrix[, n_cens]
TreeID_vec <- unique_StemID$TreeID

# Per-tree summaries broadcast back to per-stem rows.
tree_dt <- data.table(
  row_idx = seq_len(nrow(new_status_matrix)),
  TreeID = TreeID_vec,
  is_alive_last = last_col == "A"
)
tree_dt[, n_stems_in_tree := .N, by = TreeID]
tree_dt[, any_alive_at_last := any(is_alive_last), by = TreeID]
setorder(tree_dt, row_idx)

# Per-row classification of which remap (if any) applies.
single_stem_row <- tree_dt$n_stems_in_tree == 1L
all_dead_row <- (tree_dt$n_stems_in_tree > 1L) & !tree_dt$any_alive_at_last
some_alive_row <- (tree_dt$n_stems_in_tree > 1L) & tree_dt$any_alive_at_last

# Build masks the same shape as the matrix, then apply both remaps.
G_to_D_mask <- (new_status_matrix == "G") & (single_stem_row | all_dead_row)
D_to_G_mask <- (new_status_matrix == "D") & some_alive_row

corrected_new_status_matrix <- new_status_matrix
corrected_new_status_matrix[G_to_D_mask] <- "D"
corrected_new_status_matrix[D_to_G_mask] <- "G"

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
  "✓ Section 11 D/G correction: G->D=%d cells, D->G=%d cells (A/P preserved).\n",
  sum(G_to_D_mask), sum(D_to_G_mask)
))

rm(
  tree_dt, AP_in, AP_out, single_stem_row, all_dead_row, some_alive_row,
  G_to_D_mask, D_to_G_mask, last_col, TreeID_vec
)

tbl_sorted_new_status <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_new_status)

tbl_sorted_corrected_new_status <- sort_table_status(corrected_new_status, sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_corrected_new_status)

tbl_sorted_tree_histories <- sort_table_status(as.vector(tree_histories), sort_by = "x", decreasing = FALSE)
setDT(tbl_sorted_tree_histories)

# ========================================================================
# SECTION 12b: UNIFIED BIOLOGY ASSESSMENT (CUMULATIVE PIPELINE REPORT)
# ========================================================================
# PURPOSE
#   Produce a single tab-separated report that quantifies, for the four
#   stem-history versions in the pipeline, how many illegal biological
#   transitions remain. This is the canonical post-pipeline QA file:
#
#       original_status         : as decoded from the raw ViewFullTable
#                                 (before any propagation / fix)
#       new_status              : after Section 7 propagation + Sections
#                                 8, 9-10, 10b corrections
#       corrected_new_status    : after Section 11 stem-level D<->G remap
#       tree_histories          : tree-level encounter history
#
# BIOLOGY RULES (single source of truth)
#   Allowed two-cell transitions:
#       A->A  A->D  A->G                (alive -> alive / death / gone)
#       D->D                             (dead is terminal)
#       G->G                             (gone is terminal)
#       P->P  P->A                       (prior -> prior / first alive)
#   Everything else is biologically impossible.
#
# WHAT GETS REPORTED
#   - Total number of stems with at least one illegal transition.
#   - Per-version counts of every illegal 2-letter substring (DA, GA,
#     PD, PG, DG, GD, AP, DP, GP).
#   - "Stranded D/G" pattern A[DG]+A : the row has at least one D or G
#     surrounded by alive cells. After Section 10b this should be 0.
#   - One-line summary of expected behaviour for each version.
#
# OUTPUT
#   CHECK_folder / biology_assessment.txt
# ========================================================================

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
  "# allowed transitions: AA AD AG DD GG PP PA",
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
  "  3. corrected_new_st. : same as (2) (Section 11 only swaps D<->G).",
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
# SECTION 13: DATA QUALITY REPORTS AND DIAGNOSTICS
# ========================================================================
# PURPOSE:
#   Generate diagnostic CSV files to help identify data quality issues
#   and understand status transformations applied by the script
#
# OUTPUT FILES (saved to CHECK_folder):
#   1. stem_status_resurected.csv - stem with resurrection patterns
#   2. status_changed.csv - Summary of all status transformations
#   3. subset_stems_never_alive.csv - Stems with P→D/G (never alive after prior)
#   4. subset_stems_gone_to_alive.csv - Stems with GA (gone then reappeared)
#
# WHY THIS MATTERS:
#   These reports help users:
#   - Verify the script made correct decisions
#   - Identify potential data entry errors in the raw data
#   - Understand biological patterns (true resurrections, stem losses)
#   - Make informed decisions about data quality before analysis
#
# MAIN PATTERNS TO REVIEW:
#   - Resurrections (DA): Could be data errors or real sprouting
#   - Never alive (P→D/G): Stems that where never there?
#   - Reappearances (GA): Likely tagging errors or misidentifications
# ========================================================================

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
problem_df[, dbh8 := ifelse(!is.na(get("DBH_8")), 1, 0)]
problem_df[, dbh9 := ifelse(!is.na(get("DBH_9")), 1, 0)]

write.csv(unique(problem_df[, .(dbh1, dbh2, dbh3, dbh4, original_status, new_status, corrected_new_status, tree_histories)]),
  file = file.path(CHECK_folder, paste0(site, "_all_combinations_of_dbh_and_status.csv")),
  row.names = FALSE
)

# ========================================================================
# SECTION 13b: STATUS × DBH SUPPORT SUMMARY  +  FINAL AUDIT OF
#              corrected_new_status_matrix (the object that gets exported)
# ========================================================================
# PURPOSE
#   `corrected_new_status_matrix` is the canonical, fully-cleaned object
#   that Section 14 writes to disk. Before exporting, we:
#     (a) verify it contains no illegal biological transitions;
#     (b) verify A and P cells were never altered by Sections 8–11;
#     (c) summarize every cell by status code AND by whether a DBH
#         measurement is present.
#
# WHY THE DBH BREAKDOWN MATTERS
#     status="A" with no DBH  : alive but unmeasured -> needs interpolation
#     status="A" with    DBH  : healthy, real measurement
#     status="D" with    DBH  : dead but a DBH was recorded that census
#                               (broken-below code or last living measure)
#     status="D" with no DBH  : standard "dead, no measurement" cell
#     status="G" with    DBH  : multi-stem, broken-below-but-measured
#     status="G" with no DBH  : multi-stem, gone, not measured
#     status="P" with    DBH  : SHOULD BE ZERO. Any non-zero value is a
#                               pipeline bug (P means "not yet recruited").
#
# OUTPUT FILES
#     CHECK_folder/status_x_dbh_summary.csv          (overall + per census)
#     CHECK_folder/status_x_dbh_summary.txt          (human-readable)
# ========================================================================

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
# PURPOSE:
#   Format each census data frame according to ForestGEO standards and
#   export as .Rdata files for use in downstream analyses
#
# WHAT THIS SECTION DOES:
#   1. For each census, select only the required columns
#   2. Rename columns to ForestGEO standard names
#   3. Set DFstatus field (legacy field for "prior" status)
#   4. Export as .Rdata file with standardized name
#
# OUTPUT FILES:
#   [site].stem1.Rdata, [site].stem2.Rdata, ..., [site].stem[n].Rdata
#   Saved to: OUTPUT_folder
#
# FORESTGEO R TABLE FORMAT:
#   Required columns (renamed from database format):
#     - treeID (from TreeID)
#     - stemID (from StemID)
#     - TreeID (from TreeID)
#     - StemTag (from StemTag)
#     - sp (from Mnemonic - species code)
#     - quadrat (from QuadratName)
#     - gx, gy (from PX, PY - tree coordinates in meters)
#     - hom (from HOM - height of measurement)
#     - dbh (from DBH - diameter at breast height)
#     - codes (from ListOfTSM - tree status modifiers)
#     - status (from new_status - corrected status code)
#     - date (from ExactDate - census date)
#     - DFstatus (legacy field, set to "prior" when status is "P")
#
# WHY THIS FORMAT:
#   - Standardized across all ForestGEO sites
#   - Compatible with CTFSRPackage and other ForestGEO R tools
#   - Enables cross-site comparisons and analyses
#
# MAIN FUNCTIONS USED:
#   - setnames(): Rename columns in data.table
#   - :=: data.table assignment by reference
#   - assign() + save(): Create and export R objects
# ========================================================================

for (census in seq_along(ViewFullTable_split)) {
  ViewFullTable_split[[census]]$new_status <- corrected_new_status_matrix[, census]
}

# table(ViewFullTable_split[[9]][, ..ViewFullTable_columns_to_keep]$new_status, useNA = "ifany")

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
# ViewFullTable_split[[1]][, ..ViewFullTable_columns_to_keep]
# inspectdf::inspect_na(ViewFullTable_split[[1]][, ..ViewFullTable_columns_to_keep])

# FIXME: fill in dates

## create sample dataset for impute_tree_coords testing
test_dt <- list(
  data.table(
    TreeID = c(1, 1, 1, 2, 2, 2, 2),
    Tag = c("A", "A", "A", "B", "B", "B", "B"),
    PX = c(NA, NA, 2, 5, 5, 4, NA),
    PY = c(NA, NA, 3, 10, 10, 9, NA),
    QuadratName = c("Q", NA, NA, "Q1", "Q1", NA, "Q3")
  )
)

test_dt

test_dt_corr <- impute_tree_coords(
  split_list = test_dt,
  id_cols = c("TreeID", "Tag"),
  coord_cols = c("PX", "PY", "QuadratName"),
  strict = FALSE,
  multi_val = "mode"
)

test_dt_corr

ViewFullTable_split <- impute_tree_coords(
  split_list = ViewFullTable_split,
  id_cols = c("TreeID", "Tag"),
  coord_cols = c("PX", "PY", "QuadratName"),
  strict = FALSE,
  multi_val = "mode"
)
cat("✓ Location imputation complete.\n\n")

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

if (dir.exists("./BCI_stem_reconstruction/3_PREPARE_R_TABLES/VIEWFULLTABLE_CHECKS")) {
  cat("✓ CHECK folder exists: ./BCI_stem_reconstruction/3_PREPARE_R_TABLES/VIEWFULLTABLE_CHECKS\n")
} else {
  dir.create("./BCI_stem_reconstruction/3_PREPARE_R_TABLES/VIEWFULLTABLE_CHECKS", recursive = TRUE)
  cat("✓ Created CHECK folder: ./BCI_stem_reconstruction/3_PREPARE_R_TABLES/VIEWFULLTABLE_CHECKS\n")
}

fwrite(check_coordinates, "./BCI_stem_reconstruction/3_PREPARE_R_TABLES/VIEWFULLTABLE_CHECKS/repeated_coordinates.csv")




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

# ============================================================================
# TEST DATA
# ============================================================================

set.seed(123)

# Census 1
census1_test <- data.table(
  TreeID = c(1, 1, 1, 1, 2, 2, 3, 4, 5, 6),
  Tag = c(101, 101, 101, 101, 102, 102, 103, 104, 105, 106),
  QuadratName = c("Q1", "Q1", "Q1", "Q1", "Q1", "Q1", "Q2", "Q2", "Q2", "Q3"),
  ExactDate = as.Date(c(
    "1981-06-12", "1981-06-12", "1981-06-15", NA, # Tree 1: multiple obs, mode = 1981-06-12
    "1981-06-12", NA, # Tree 2: one value, one NA
    NA, # Tree 3: NA (should get Q2 mode)
    "1981-06-20", # Tree 4: has value
    NA, # Tree 5: NA (should get Q2 mode)
    NA # Tree 6: NA (Q3 has no other dates)
  )),
  PX = c(10, 10, 10, 10, 20, 20, 30, 40, 50, 60),
  PY = c(15, 15, 15, 15, 25, 25, 35, 45, 55, 65)
)

# Census 2
census2_test <- data.table(
  TreeID = c(1, 2, 2, 2, 3, 4, 5, 6, 7),
  Tag = c(101, 102, 102, 102, 103, 104, 105, 106, 107),
  QuadratName = c("Q1", "Q1", "Q1", "Q1", "Q2", "Q2", "Q2", "Q3", "Q3"),
  ExactDate = as.Date(c(
    "1985-07-10", # Tree 1: different date than census 1 (expected!)
    "1985-07-10", "1985-07-10", "1985-07-12", # Tree 2: mode = 1985-07-10
    NA, # Tree 3: NA (should get Q2 value)
    "1985-07-15", # Tree 4: has value
    "1985-07-15", # Tree 5: has value
    NA, # Tree 6: NA
    "1985-08-01" # Tree 7: has value
  )),
  PX = c(10, 20, 20, 20, 30, 40, 50, 60, 70),
  PY = c(15, 25, 25, 25, 35, 45, 55, 65, 75)
)

test_data <- list(census1_test, census2_test)

# ============================================================================
# RUN TESTS
# ============================================================================

cat("\n===== STRICT MODE (only within-tree imputation) =====\n")
result_strict <- impute_tree_dates(test_data, strict = TRUE)

cat("\n\nCensus 1 results (strict):\n")
print(result_strict[[1]][order(TreeID)])

cat("\n\nCensus 2 results (strict):\n")
print(result_strict[[2]][order(TreeID)])

cat("\n\n===== NON-STRICT MODE (tree + quadrat imputation) =====\n")
result_nonstrict <- impute_tree_dates(test_data, strict = FALSE)

cat("\n\nCensus 1 results (non-strict):\n")
print(result_nonstrict[[1]][order(TreeID)])

cat("\n\nCensus 2 results (non-strict):\n")
print(result_nonstrict[[2]][order(TreeID)])

# ============================================================================
# CHECK TESTS
# ============================================================================

result_strict[[1]] <- result_strict[[1]][, .(TreeID, QuadratName, ExactDate)]
result_strict[[2]] <- result_strict[[2]][, .(TreeID, QuadratName, ExactDate)]

result_nonstrict[[1]] <- result_nonstrict[[1]][, .(TreeID, QuadratName, ExactDate)]
result_nonstrict[[2]] <- result_nonstrict[[2]][, .(TreeID, QuadratName, ExactDate)]

list(
  census1_test[, .(TreeID, QuadratName, ExactDate)],
  census2_test[, .(TreeID, QuadratName, ExactDate)]
)
result_strict
result_nonstrict

ViewFullTable_split <- impute_tree_dates(
  split_list = ViewFullTable_split,
  id_cols = c("TreeID", "Tag"),
  date_col = "ExactDate",
  quadrat_col = "QuadratName",
  strict = FALSE
)
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

inc <- unique(check_dates[n_dates > 1L]$TreeID)

for (census in seq_along(ViewFullTable_split)) {
  cat(sprintf("  Processing census %d...\n", census))
  # Extract current census data and ensure it's a data.table
  X <- as.data.table(ViewFullTable_split[[census]])
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

# sort(unique(check_data$Mnemonic))

# ========================================================================
# SECTION 15: EXPORT SPECIES TABLE
# ========================================================================
# PURPOSE:
#   Create and export the species taxonomy table in ForestGEO format
#   Maps species codes (Mnemonic) to full taxonomic information
#
# WHAT THIS SECTION DOES:
#   1. Copy ViewTaxonomy data (to avoid modifying original)
#   2. Sort species alphabetically by Mnemonic code
#   3. Create Latin name by combining Genus + SpeciesName
#   4. Select and reorder columns for ForestGEO format
#   5. Rename columns to ForestGEO standard names
#   6. Export as [site].spptable.rdata
#
# OUTPUT FILE:
#   [site].spptable.rdata (e.g., hkk.spptable.rdata)
#   Saved to: OUTPUT_folder
#
# FORESTGEO SPECIES TABLE FORMAT:
#   Required columns:
#     - sp (species code - the Mnemonic like "ACERUB")
#     - Latin (full scientific name like "Acer rubrum")
#     - Genus (genus name)
#     - SpeciesName (species epithet)
#     - Family (taxonomic family)
#     - speciesID (numeric species identifier)
#     - authority (taxonomic authority who named the species)
#     - IDLevel (identification confidence level)
#     - oldnames (historical names/synonyms)
#     - subsp (subspecies information)
#
# WHY THIS MATTERS:
#   - Links species codes in census tables to full taxonomy
#   - Enables filtering by family, genus, or species
#   - Documents taxonomic authorities and synonyms
#   - Essential for cross-site comparisons with standardized names
#
# MAIN FUNCTIONS USED:
#   - copy(): Create independent copy of data.table
#   - setorder(): Sort data.table by reference
#   - paste(): Concatenate Genus and SpeciesName
#   - setnames(): Rename columns to ForestGEO standards
# ========================================================================

# --------------------------------------------------------------------
# --------------------------------------------------------------------
ViewTaxonomy <- fread(
  file = file.path("./BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.txt"),
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
  "Family", # Taxonomic family
  "Latin", # Full scientific name
  "Genus", # Genus name
  "SpeciesName", # Species epithet
  "Subspecies", # Subspecies information
  # "SpeciesID", # Numeric species identifier
  "Authority", # Taxonomic authority
  # "IDLevel", # Identification confidence level
  "Synonyms", # Historical synonyms
  "Lifeform_RFoster",
  "Lifeform_RPerez_SAguilar",
  "CommonName", # Common name (if available)
  "Herbarium" # Herbarium voucher information (if available)
)

sptable <- sptable[, ..cols_to_keep]

# Rename columns to ForestGEO standard names
# Database names → ForestGEO R names
setnames(sptable,
  old = c(
    "Mnemonic"
  ),
  new = c(
    "sp"
  )
)

# sptable[!is.na(Subspecies)]

# Filter to only species that actually appear in the census data
# This removes species from the taxonomy list that aren't in this plot
species_in_census <- sort(unique(ViewFullTable$Mnemonic))
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
