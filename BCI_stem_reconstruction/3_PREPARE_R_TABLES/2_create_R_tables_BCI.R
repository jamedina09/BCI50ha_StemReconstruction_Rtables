################################################################################
# FORESTGEO STEM TABLE CREATION SCRIPT
################################################################################
# Purpose: Build the ForestGEO R tables (one per census) from the
# reconstructed census dataset of stage 2.
#
# Run from the project root (paths are built from the working directory),
# after 2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R.
#
# Inputs:
#   - BCI_stem_reconstruction/DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds
#
# Outputs:
#   - BCI_stem_reconstruction/DATA/RTABLES/[site].stem[n].Rdata (and [site].stem[n].csv)
#     (the species table is built from these tables by 3_create_species_table_BCI.R)
#   - BCI_stem_reconstruction/DATA/CHECKS/, four files:
#       check_[site].txt                    every check (✓ / ❌) and the summary
#                                           tables, in run order: read this first
#       status_patterns_[site].csv          how the raw records become Rstatus:
#                                           one row per encounter history
#       field_records_to_review_[site].csv  raw records that need a look (data
#                                           manager): one row per case, `reason`
#       reconstruction_to_review_[site].csv stems whose reconstruction needs a
#                                           look: one row per case, `reason`
#
# Packages: data.table, inspectdf.
#
# Main pipeline:
#    1. Setup and load input data                          (Sections 1–2)
#    2. Standardize stems across all censuses                (Section 3)
#    3. Map raw Status values to raw codes A/D/N           (Section 4)
#    4. Validate encounter histories before propagation      (Section 5)
#    5. Propagate no-data states into complete histories     (Section 6)
#    6. Check first records; fill each stem's lifespan       (Sections 7–8)
#    7. Assert no D/G is stranded between two A's            (Section 9)
#    8. Derive tree-level histories; apply tree-aware D/G    (Section 10)
#    9. Assess biology across all history versions          (Section 11)
#   10. Status patterns and review lists; status × DBH audit (Sections 12–13)
#   11. One location per tree; record dates; export per-census R tables;
#       write the review lists (Section 14)
# The status and DBH rules live in rstatus_functions.R (same folder), which
# tests/test_Rstatus.R also uses.
################################################################################
# Variable summary:
#   Status               raw field record, mapped to A/D/N (Section 4)
#   original_status      stem encounter history before propagation
#   new_status           propagated stem history before D/G remap
#   corrected_new_status final stem history after tree-level D/G adjustment
#   tree_histories       tree-level history aggregated from stem histories
#   DBHs                 numeric matrix of raw stem DBH measurements by census
#                        (the exported dbh, unchanged)
#   Rstatus              exported per-census corrected stem status
#   DFstatus             raw field status (ForestGEO column name), never modified
################################################################################
# BIOLOGICAL CONTRACT (checked with bio_check() before export)
#   Stem identity  : taken from the DP reconstruction (stage 2) and never
#                    changed here. Suspected identity breaks are only reported
#                    (CHECKS/reconstruction_to_review_[site].csv).
#   Evidence of life: a raw "alive" record, or a DBH on a "broken below",
#                    "missing" or status-less record.
#   Dead record    : "dead" / "stem dead" (with or without DBH) and "broken
#                    below" without DBH. A DBH on a dead record is not
#                    evidence of life (the stem is A there only if it is
#                    alive later); the DBH itself is kept.
#   P (prior)      : only before a stem's first record.
#   A (alive)      : from a stem's first record to its last evidence of life.
#                    Missed censuses and false deaths in between are A:
#                    alive later => never dead.
#   Stump          : an A without DBH whose DFstatus is "broken below": the
#                    stem lives later (e.g. a field-tagged stump measured in a
#                    later census), so it is alive, with no measurable DBH at
#                    1.3 m. Stage 4 keeps stump censuses out of the stock.
#   Dead (G or D)  : after the last evidence of life, from the first census
#                    without one; a stem never alive is dead from its first
#                    record.
#   Tree alive at census j: some stem of the tree is A at j OR at any later
#                    census. A tree that loses all its registered stems and
#                    later has a living stem was alive the whole time; it was
#                    just not registered.
#   G (dead stem)  : dead stem while its tree is alive at j.
#   D (dead tree)  : dead stem and no stem of the tree is A at j or later.
#                    D is absorbing: nothing in the tree lives after.
#   DBH            : exported as recorded on every record, A, G or D (G/D
#                    cells with a DBH: reason dbh_on_dead_record in
#                    CHECKS/field_records_to_review_[site].csv); a P cell
#                    never has one; never imputed.
#   ExactDate      : the field date exactly as recorded (NA where there is no
#                    record); never imputed, like dbh and DFstatus.
#   date           : days since 1960-01-01 (the ForestGEO convention) on every
#                    row: the recorded date, else the modal field date of the
#                    tree, quadrat or census in that census (Section 14).
#   Location       : one (gx, gy, quadrat) per tree, the same in every census
#                    and for every stem (Section 14).
#   stemID         : the reconstructed stem number within its tree (1, 2,
#                    ...): treeID + stemID identify a stem. Internally the
#                    script keys stems on "TreeID_ReconstructedStemID".
#   Legal stem transitions: PP PA PG PD AA AG AD GG GD DD
#   Legal tree transitions: PP PA PD AA AD DD
#   Everything else is illegal (e.g. DG, DA, GA, AP).
################################################################################

# ========================================================================
# SECTION 1: SETUP AND INITIALIZATION
# ========================================================================
# Load libraries, configure paths, and define output column mappings.
# This section also ensures required folders exist before the pipeline runs.

# Remove all objects from workspace to start fresh
rm(list = ls())

# Load libraries ####
library(data.table)
# References 1: https://arelbundock.com/posts/dt_tb_df/index.html
# References 2: https://rdatatable.gitlab.io/data.table/
# library(ggplot2)

getDTthreads() # Show number of threads data.table will use
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

# 3. Input folder containing the DP-reconstructed census table
INPUT_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA", "PROCESSED")
# ℹ️ File required:
#    - complete_dataset_final_with_reconstructed_stemids.rds (stage 2 output)

# 4. Output folder for final .Rdata files
OUTPUT_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA", "RTABLES")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
  cat("✓ Created OUTPUT folder:", OUTPUT_folder, "\n")
}

# 5. Diagnostics folder for QA/QC reports
CHECK_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA", "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
  cat("✓ Created CHECKS folder:", CHECK_folder, "\n")
}

# 6. Log file for the messages of the checks
log_file <- file.path(CHECK_folder, paste0("check_", site, ".txt"))
# Create an empty log file (an existing one is replaced)
create_log_file <- function(log_file) {
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  file.create(log_file)
}

create_log_file(log_file)

# Append a message to the log file (optionally preceded by a separator line)
print_to_log <- function(message, log_file, new_message = TRUE) {
  if (new_message) {
    cat(strrep("-", 60), file = log_file, append = TRUE, sep = "\n")
  }
  cat(message, file = log_file, append = TRUE, sep = "\n")
}

# A table as text lines for the log, printed at a fixed width so that the log
# does not depend on the width of the console
log_table <- function(x) {
  op <- options(width = 250L)
  on.exit(options(op))
  capture.output(print(x))
}

print_to_log(c(
  sprintf("check_%s.txt: checks and summaries of 2_create_R_tables_BCI.R, run %s", site, format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Every check is listed (✓ passed, ❌ failed: the script stopped there). The review lists are",
  sprintf("field_records_to_review_%s.csv and reconstruction_to_review_%s.csv (counts by reason at the end).", site, site)
), log_file, new_message = FALSE)

# ── Hard check that stays visible in interactive (line-by-line) runs ────────
# When this script is run line by line in the R console, a bare stop() is
# easy to miss. On failure bio_check() prints a ❌ line (count + examples),
# raises an immediate warning, writes the same text to the log file, and only
# then stops. On success it prints a ✓ line and writes it to the log, so the
# log lists every check.
#   ok       : TRUE when the check passes (NA counts as failure)
#   msg      : the requirement being checked, phrased as the expected state
#   examples : optional IDs of offending stems/trees (first 10 are shown)
#   n_bad    : optional number of offending cases
bio_check <- function(ok, msg, examples = NULL, n_bad = NULL) {
  if (isTRUE(all(ok))) {
    cat("✓", msg, "\n")
    print_to_log(paste("✓", msg), log_file, new_message = FALSE)
    return(invisible(TRUE))
  }
  n_txt <- if (!is.null(n_bad)) sprintf(" [%d case(s)]", n_bad) else ""
  ex_txt <- if (length(examples) > 0L) {
    paste0(" | examples: ", paste(head(unique(examples), 10), collapse = ", "))
  } else {
    ""
  }
  full_msg <- paste0("CHECK FAILED: ", msg, n_txt, ex_txt)
  cat("❌", full_msg, "\n")
  warning(full_msg, call. = FALSE, immediate. = TRUE)
  print_to_log(full_msg, log_file, new_message = TRUE)
  stop(full_msg, call. = FALSE)
}

# ── Legal status transitions (biological contract, see header) ──────────────
# Single source of truth for every validator in Sections 5–13.
status_codes <- c("A", "D", "G", "P")
valid_stem_trans <- c("PP", "PA", "PG", "PD", "AA", "AG", "AD", "GG", "GD", "DD")
valid_tree_trans <- c("PP", "PA", "PD", "AA", "AD", "DD")

# Status and DBH rules (shared with tests/test_Rstatus.R)
source(file.path(main_path, "BCI_stem_reconstruction", "3_PREPARE_R_TABLES", "rstatus_functions.R"))

# Convert a set of legal 2-letter transitions into the list of illegal
# (first, second) pairs expected by check_histories().
make_invalid_transitions <- function(valid_trans, codes = status_codes) {
  pairs <- expand.grid(first = codes, second = codes, stringsAsFactors = FALSE)
  pairs <- pairs[!paste0(pairs$first, pairs$second) %in% valid_trans, ]
  lapply(seq_len(nrow(pairs)), function(i) c(pairs$first[i], pairs$second[i]))
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
#   "obs_row_id"                    — row number of the record within its tag in stage 2
#                                     (exported as "StemPaths", the key of the posterior paths)
#   "sp"                            — species mnemonic (short ForestGEO species code)
#   "gx", "gy"                      — tree X/Y coordinates in meters within the plot
#   "hom"                           — height of measurement (m above ground where DBH taken)
#   "DFstatus"                      — ForestGEO column name for the raw field status
#   "Rstatus"                       — script-computed corrected status (from "new_status")

# # Columns to retain from ViewFullTable (original DB column names)
# ViewFullTable_columns_to_keep <- c(
#   "TreeID", "StemID", "Tag", "StemTag", "Mnemonic", # stem / tree identifiers
#   "QuadratName", "PX", "PY", # spatial location in plot (meters)
#   "DBHID", "CensusID", "DBH", "HOM", # measurement record
#   "ExactDate", "Status", "ListOfTSM", # raw field status + measurement codes
#   "Date", "new_status", # date + script-computed corrected status
#   "obs_row_id" # DP reconstruction metadata
# )

# # Standardized ForestGEO R-table column names (must match ViewFullTable_columns_to_keep
# # one-to-one: same length and same order)
# new_names_columns_to_keep <- c(
#   "treeID", "stemID", "tag", "StemTag", "sp", # sp = species mnemonic
#   "quadrat", "gx", "gy", # gx/gy = plot coordinates in meters
#   "MeasureID", "CensusID", "dbh", "hom", # hom = height of measurement
#   "ExactDate", "DFstatus", "codes", # DFstatus = raw field status (legacy name)
#   "date", "Rstatus", # date + script-computed corrected status
#   "StemPaths" # DP reconstruction metadata
# )

# Columns to retain from ViewFullTable (original DB column names)
ViewFullTable_columns_to_keep <- c(
  "CensusID", "TreeID", "StemID", "Tag", "StemTag", "Mnemonic", # stem / tree identifiers
  "QuadratName", "PX", "PY", # spatial location in plot (meters)
  "DBH", "HOM",
  "ExactDate", "Date", # raw field date; Date = filled date (days since 1960-01-01)
  "ListOfTSM", # field codes
  "Status", # raw field status
  "new_status", # script-computed corrected status
  "obs_row_id" # row number within the tag in stage 2
)

# Standardized ForestGEO R-table column names (must match ViewFullTable_columns_to_keep
# one-to-one: same length and same order)
new_names_columns_to_keep <- c(
  "CensusID", "treeID", "stemID", "tag", "StemTag", "sp", # sp = species mnemonic
  "quadrat", "gx", "gy", # gx/gy = plot coordinates in meters
  "dbh", "hom",
  "ExactDate", "date", # ExactDate raw; date = days since 1960-01-01, every row
  "codes",
  "DFstatus", # DFstatus = raw field status
  "Rstatus", # script-computed corrected status
  "StemPaths" # obs_row_id of stage 2: the key of the posterior path files
)

# ---- Hard guard: the two column vectors must be 1-to-1 -------------------
# A typo or accidental edit to either vector would silently mis-rename
# columns at export time. Catch it now, before the pipeline starts.
bio_check(
  length(ViewFullTable_columns_to_keep) == length(new_names_columns_to_keep) &&
    !anyDuplicated(ViewFullTable_columns_to_keep) &&
    !anyDuplicated(new_names_columns_to_keep),
  "Export column vectors map 1-to-1 with no duplicates"
)

# ========================================================================
# SECTION 2: LOAD INPUT DATA
# ========================================================================
# Load the reconstructed census table and coerce key types.

# ========================================================================
# LOAD RECONSTRUCTED CENSUS DATA (RDS FROM STAGE 2)
# ========================================================================

cat("\n", strrep("=", 70), "\n")
cat("LOADING DATA FOR SITE:", toupper(site), "\n")
cat(strrep("=", 70), "\n\n")

# Load ViewFullTable_RAW (census data) ####
cat("📂 Loading ViewFullTable_RAW...\n")

ViewFullTable_RAW <- as.data.table(readRDS(file.path(
  INPUT_folder,
  "complete_dataset_final_with_reconstructed_stemids.rds"
)))
# Stem identity comes from the DP and must survive the numeric normalisation
# unchanged: a non-numeric ID would silently become NA and be merged into the
# "TreeID_NA" pseudo-stem, overriding the DP.
raw_rsid <- as.character(ViewFullTable_RAW$ReconstructedStemID)
bad_rsid <- unique(raw_rsid[!is.na(raw_rsid) & is.na(suppressWarnings(as.numeric(raw_rsid)))])
bio_check(
  length(bad_rsid) == 0L,
  "Every non-NA ReconstructedStemID is numeric (no DP identity lost on conversion)",
  examples = bad_rsid,
  n_bad = length(bad_rsid)
)
rm(raw_rsid, bad_rsid)
ViewFullTable_RAW[, ReconstructedStemID := as.numeric(as.character(ReconstructedStemID))] # normalise IDs (e.g. "01" -> 1); made character below

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
# Goal: see how many rows with a database StemID have a ReconstructedStemID,
#   per census (multiple-stem tags). That every row with a DBH has one is
#   checked with bio_check() further below.
ViewFullTable_RAW[!is.na(StemID) & single_stem_tags == FALSE, .(
  n_stems = .N,
  na_levels_StemID = sum(is.na(StemID)),
  na_levels_ReconstructedStemID = sum(is.na(ReconstructedStemID)),
  n_stems_with_StemID = sum(!is.na(StemID)),
  n_stems_with_ReconstructedStemID = sum(!is.na(ReconstructedStemID)),
  n_stems_with_both_IDs = sum(!is.na(StemID) & !is.na(ReconstructedStemID)),
  n_stems_with_neither_ID = sum(is.na(StemID) & is.na(ReconstructedStemID))
), by = CensusID][order(CensusID)]

# Working copy
ViewFullTable <- copy(ViewFullTable_RAW)

# ── Overwrite StemID with the stable key TreeID_ReconstructedStemID ─────────────
# WHY: The database StemID is not stable across censuses: the same physical
#   stem can carry different StemIDs in different censuses. The
#   ReconstructedStemID of stage 2 follows the same stem across all censuses.
# CONSTRUCTION: paste(TreeID, ReconstructedStemID, sep="_") produces a globally
#   unique key combining tree identity (TreeID) with stem identity (ReconstructedStemID).
# Raw_StemID: the original DB StemID is preserved in this new column for
#   traceability and cross-checking against the raw database export.
# IMPORTANT: Rows with NA ReconstructedStemID become StemID = "TreeID_NA".
#   These are rows without a stem identity from stage 2 (in the current data:
#   placeholder rows of a census in which the tree has no record, with no
#   status and no DBH). That they carry no DBH is checked below; they are
#   dropped at the start of Section 3.
ViewFullTable[, Raw_StemID := StemID] # preserve original DB StemID

# unique(ViewFullTable$StemID)
# unique(ViewFullTable$Raw_StemID)
# unique(ViewFullTable$ReconstructedStemID)

setkey(ViewFullTable, Tag)

## DIAGNOSTIC: Verify Tag <-> TreeID 1:1 mapping
ViewFullTable[, .(n_tags = uniqueN(Tag)), by = TreeID][order(-n_tags)]
ViewFullTable[, .(n_treeid = uniqueN(TreeID)), by = Tag][order(-n_treeid)]

ViewFullTable[is.na(Tag)]
ViewFullTable[is.na(TreeID)]

# sort(table(ViewFullTable$ReconstructedStemID))
ViewFullTable[, StemID := paste(TreeID, ReconstructedStemID, sep = "_")] # new stable StemID

# ── DIAGNOSTIC: Report "TreeID_NA" StemIDs ───────────────────────────────────────
# Any row with NA ReconstructedStemID becomes StemID = "TreeID_NA" after the
# paste() above. Their count, status distribution and census distribution are
# reported so the scale can be verified.
stemid_problem <- ViewFullTable[grepl("_NA$", StemID)]
cat("\nRows with 'TreeID_NA' StemID (dead/broken stems without DP identity):", nrow(stemid_problem), "\n")
# A measured stem must keep its DP identity: if a row with a DBH had no
# ReconstructedStemID it would be pooled into the "TreeID_NA" pseudo-stem.
bio_check(
  stemid_problem[!is.na(DBH), .N] == 0L,
  "Rows without a DP stem identity (StemID 'TreeID_NA') carry no DBH",
  examples = stemid_problem[!is.na(DBH), StemID],
  n_bad = stemid_problem[!is.na(DBH), .N]
)
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
TreeID_ranges <- xraw_unique[, .(min_c = min(as.numeric(CensusID)), max_c = max(as.numeric(CensusID))), by = TreeID]
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
StemID_ranges <- xraw_unique[, .(min_c = min(as.numeric(CensusID)), max_c = max(as.numeric(CensusID))), by = StemID]
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

# ========================================================================
# SECTION 3: STANDARDIZE DATA FORMAT ACROSS CENSUSES
# ========================================================================
# Build a master stem list and force every census table to include the same
# StemIDs in the same row order. Missing stems are retained as NA rows so that
# encounter histories can be assembled consistently across all censuses.

# --------------------------------------------------------------------
# Create master stem list with fixed attributes
# --------------------------------------------------------------------

# Placeholder rows without a stem identity are dropped. Stage 1 adds a
# placeholder row (no StemID, no status, no DBH, no codes) for a census in
# which the tree has no record between its first and last record; stage 2
# gives it a stem identity only when one stem can own it. Kept, the others
# would become one "TreeID_NA" stem per tree, P in every census, a stem that
# does not exist. The tree's real stems keep their rows in those censuses.
no_identity <- ViewFullTable[, grepl("_NA$", StemID)]
bio_check(
  ViewFullTable[no_identity, all(is.na(ReconstructedStemID) & is.na(Raw_StemID) & is.na(Status) & is.na(DBH) & is.na(ListOfTSM))],
  "Rows without a stem identity are empty placeholder rows (no stem identity, StemID, status, DBH or codes)",
  examples = ViewFullTable[no_identity & !(is.na(Raw_StemID) & is.na(Status) & is.na(DBH) & is.na(ListOfTSM)), unique(TreeID)],
  n_bad = ViewFullTable[no_identity & !(is.na(Raw_StemID) & is.na(Status) & is.na(DBH) & is.na(ListOfTSM)), .N]
)
cat(sprintf(
  "Placeholder rows without a stem identity dropped: %d rows in %d trees\n",
  sum(no_identity), ViewFullTable[no_identity, uniqueN(TreeID)]
))
ViewFullTable <- ViewFullTable[!no_identity]
bio_check(!any(grepl("_NA$", ViewFullTable$StemID)), "Every remaining row has a stem identity (no 'TreeID_NA' stem)")
rm(no_identity)

# fixed_columns: the minimum identifier set that does not change across
# censuses. StemID below is the composite TreeID_ReconstructedStemID created
# in Section 2 (the database StemID is kept in Raw_StemID).
fixed_columns <- c(
  "PlotName", "PlotID",
  "Mnemonic",
  "QuadratName", "QuadratID", # "PX", "PY", # "QX", "QY",
  "TreeID",
  "Tag",
  "StemID" # StemTag is not fixed: before 2010 the database has no stem tags
)

# Convert to data.table if not already (defensive coding)
setDT(ViewFullTable)

# Create master stem list: One row per stem, fixed attributes only
# unique() removes duplicate rows (same stem in multiple censuses)
unique_StemID <- unique(ViewFullTable[, ..fixed_columns])

## return the duplicated rows in unique_StemID using StemID as the identifier
duplicated_stems <- unique_StemID[duplicated(unique_StemID$StemID) | duplicated(unique_StemID$StemID, fromLast = TRUE)]

# VERIFY: If duplicates exist, there's a data quality issue (same stem with
# different fixed attributes). Every later step assumes one row per StemID.
bio_check(
  anyDuplicated(unique_StemID$StemID) == 0L,
  "Master stem list has one row per StemID (fixed attributes constant across censuses)",
  examples = duplicated_stems$StemID,
  n_bad = uniqueN(duplicated_stems$StemID)
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
rows_recorded <- vapply(ViewFullTable_split_unbalanced, nrow, integer(1))

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
  # than once in X. Stop before it happens so the root cause can be fixed.
  dup_ids <- X[duplicated(StemID), unique(StemID)]
  bio_check(
    length(dup_ids) == 0L,
    sprintf(
      "Census %s: each StemID appears once (match() would otherwise keep only the first row)",
      paste(unique(na.omit(X$CensusID)), collapse = "/")
    ),
    examples = dup_ids,
    n_bad = length(dup_ids)
  )
  # REORDER + FILL MISSING: match() returns indices or NA for missing stems
  # This aligns X's stems to match unique_StemID's order
  # NA indices create new rows filled with NA
  # NOTE: This relies on StemID being unique in unique_StemID (verified earlier)
  idx <- match(unique_StemID$StemID, X$StemID)
  X <- X[idx]
  # FILL FIXED ATTRIBUTES: Copy all fixed columns from master list
  # This keeps species code, quadrat and IDs identical across censuses
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

# NOTE: match() would keep only the first row of a StemID that appears twice
# in a census; the check inside the loop stops the script before that happens.

## VERIFICATION: Check that standardization worked ####
# All censuses should now have IDENTICAL dimensions (same # rows, same # columns)
cat("\n📊 Census dimensions (AFTER standardization):\n")
lapply(ViewFullTable_split, dim)

print_to_log(c(
  "Rows per census: recorded, and after standardization (every stem in every census):",
  log_table(data.table(
    census = seq_along(ViewFullTable_split), rows_recorded = rows_recorded,
    rows_standardized = vapply(ViewFullTable_split, nrow, integer(1))
  ))
), log_file, new_message = TRUE)
rm(rows_recorded)

# Show first few rows of each census for visual inspection
cat("\n📋 First few rows of each census:\n")
lapply(ViewFullTable_split, head, 4)

## VERIFY: check wether length StemID is same as nrow of each census and print that it does
for (i in seq_along(ViewFullTable_split)) {
  n_stems <- nrow(unique_StemID)
  n_rows <- nrow(ViewFullTable_split[[i]])
  bio_check(
    n_stems == n_rows && identical(ViewFullTable_split[[i]]$StemID, unique_StemID$StemID),
    sprintf("Census %d: rows (%d) match the master stem list (%d) in the same order", i, n_rows, n_stems)
  )
}

# Show first few rows of each census for visual inspection
cat("\n📋 First few rows of each census:\n")
lapply(ViewFullTable_split, head, 4)

# Explore NAs
lapply(ViewFullTable_split, inspectdf::inspect_na)
# (The fixed columns are copied from the master stem list into every census
# above, so they are the same in every census.)

# ========================================================================
# SECTION 4: STATUS CODE TRANSFORMATION
# ========================================================================
# Convert raw Status terms to the raw codes A (alive), D (dead) and N (no
# data), resolve "broken below" using DBH, and assemble per-stem encounter
# history strings for later propagation and correction. P is assigned in
# Section 6 and G in Section 10; neither comes from raw data.

# --------------------------------------------------------------------
# Raw records -> codes A / D / N (rstatus_raw_code(), rstatus_functions.R)
# --------------------------------------------------------------------
# Gather Status AND DBH from all censuses into one long data.table. The DBH
# matters because a DBH on a "broken below", "missing" or status-less record
# is evidence of life: the stem (or its resprout) was measured.
DT_Status <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(StemID, Status, DBH, census = i)]
}), use.names = TRUE, fill = TRUE)

bio_check(
  identical(seq_along(ViewFullTable_split), unique(DT_Status$census)),
  "Census numbering in DT_Status is 1..n in order"
)

# ------------------------------------------------------------------------
# Rules, applied to each stem x census record by rstatus_raw_code():
#   A  "alive" (with or without DBH), and "broken below" / "missing" / no
#      status WITH a DBH (measured, so alive).
#   D  "dead" and "stem dead" (with or without DBH), and "broken below"
#      without DBH. A D followed later by an A is a false death: the
#      lifespan rule (Section 8) turns it into A (a broken-below record then
#      becomes a stump: alive, no measurable DBH; counted in the log, Section
#      13); when the tree lives on through other stems, Section 10 turns a
#      real death into G.
#   N  "missing" or no status without DBH, and no row (no record).
# A DBH on a "dead" / "stem dead" record stays on a D here. If the stem is
# alive later (false death), Section 8 makes the cell A; otherwise the stem is
# dead and the cell becomes G or D. In both cases the DBH is exported as
# recorded (Section 13 lists the DBHs that sit on G/D cells).
# P is assigned in Section 6 and G in Section 10; neither comes from raw data.
# ------------------------------------------------------------------------
table(DT_Status$Status, useNA = "ifany")

# Every raw status must be one this script knows how to interpret. Unknown
# values would get no defined meaning (rstatus_raw_code() also refuses them),
# so stop and ask for the new value to be mapped.
known_raw_status <- c("alive", "dead", "stem dead", "broken below", "missing")
unknown_status <- setdiff(unique(na.omit(DT_Status$Status)), known_raw_status)
bio_check(
  length(unknown_status) == 0L,
  paste0("Every raw Status is one of: ", paste(known_raw_status, collapse = ", "), " (or NA)"),
  examples = unknown_status,
  n_bad = length(unknown_status)
)

# A DBH can be evidence of life, so every DBH must be a real measurement: a
# 0 or negative "no value" placeholder would create life where there is none.
bad_dbh <- DT_Status[!is.na(DBH) & DBH <= 0]
bio_check(
  nrow(bad_dbh) == 0L,
  "Every recorded DBH is a positive measurement (no 0 / negative placeholders)",
  examples = bad_dbh$StemID,
  n_bad = nrow(bad_dbh)
)
rm(bad_dbh)

# Log the records the DBH decides: non-"alive" records that are alive because
# of their DBH, and dead records that carry a DBH.
dbh_alive_tab <- DT_Status[
  !is.na(DBH) & (is.na(Status) | Status %in% c("broken below", "missing")),
  .N,
  by = .(raw_status = Status)
][order(-N)]
dead_dbh_tab <- DT_Status[!is.na(DBH) & Status %in% c("dead", "stem dead"), .N, by = .(raw_status = Status)][order(-N)]
cat("\n🔧 Records alive because of their DBH (by raw status):\n")
print(dbh_alive_tab)
cat("🔧 Dead records with a DBH (dead records; the DBH is kept; the stem is A there only if alive later):\n")
print(dead_dbh_tab)
print_to_log("Records alive because of their DBH (by raw status):", log_file, new_message = TRUE)
print_to_log(log_table(dbh_alive_tab), log_file, new_message = FALSE)
print_to_log("Dead records with a DBH (dead records; the DBH is kept; the stem is A there only if alive later):", log_file, new_message = TRUE)
print_to_log(log_table(dead_dbh_tab), log_file, new_message = FALSE)
cat(sprintf("🔧 'broken below' without DBH coded as dead: %d records\n", DT_Status[Status %in% "broken below" & is.na(DBH), .N]))

DT_Status[, code := rstatus_raw_code(Status, DBH)]

# The pivot below must place exactly one value in each StemID × census cell;
# with duplicates dcast() would silently count values instead.
bio_check(
  DT_Status[, .N, by = .(StemID, census)][, all(N == 1L)],
  "One Status value per StemID × census before the pivot"
)

# Show frequency of raw status values by census
sum_status_pre <- DT_Status[, .(nobs = .N), by = census:Status][order(Status, census)]
cat("\n📊 Full Status by Census before transformation:\n")
sum_status_pre[, .N, by = .(Status, census)][N > 1]

# Pivot long to wide: one row per StemID, one column per census
original_status_wide <- dcast(DT_Status,
  formula = StemID ~ census,
  value.var = "code"
) # fill cells with the codes A / D / N

# Find reordering indices: match master stem order (unique_StemID) to current rows
idx <- match(
  unique_StemID$StemID, # desired order (from master list)
  original_status_wide$StemID
) # current StemID order in wide table

# Reorder rows to align with unique_StemID
original_status_wide <- original_status_wide[idx] # subset/reorder using idx
original_status_full <- as.matrix(original_status_wide[, -1])
original_status <- original_status_full # Create working copy (preserve original)
head(original_status_wide)
tail(original_status_wide)

bio_check(
  identical(original_status_wide$StemID, unique_StemID$StemID) && all(original_status %in% c("A", "D", "N")),
  "Raw code matrix is in master stem order and holds only A, D and N"
)

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
  log_table(
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
# Frequency table of the values of x (a data.frame with columns x and Freq),
# sorted by `sort_by` ("x" or "Freq"). Stops on any other column name.
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

# Allowed codes in RAW histories: A (alive), D (dead), N (no data).
# G never comes from raw data: it is derived from tree context in Section 10
# (see the biological contract in the header), so it is not a raw code.
allowed_codes <- c("A", "D", "N")

# Raw anomalies to be resolved by later sections (informational here; the
# contract is enforced with bio_check() after Sections 7-10):
invalid_transitions <- list(
  c("D", "A"), # dead then alive: resolved by the lifespan rule (Section 8)
  c("A", "N"), # alive then no record: death inferred, AN -> AD (Section 6)
  c("D", "N"), # dead then no record: stays dead, DN -> DD (Section 6)
  c("N", "D") # no record, then dead: dead from the first census without evidence of life (Sections 6, 10)
)

# ----------------------------
# check_histories(): list the problems of each encounter history
# ----------------------------
#   histories           : character vector, one string per history
#   allowed_codes       : codes a history may contain
#   invalid_transitions : list of c(first, second) pairs that may not follow
#                         each other
#   nchars              : expected length (number of censuses)
# Rules: (1) length, (2) allowed codes, (3) first census is not D or G,
# (4) no invalid transition, (5) no A after a D or G.
# Returns a data.frame(History, Issue), one row per problem (no rows = clean).
check_histories <- function(histories, allowed_codes, invalid_transitions, nchars = length(ViewFullTable_split)) {
  # Data frame to store issues
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
      history_issues <- c(history_issues, "Starts with D or G (dead in the first census)")
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
# Print all detected issues
# ----------------------------
if (nrow(issues) == 0) {
  cat("No issues detected in any histories.\n")
} else {
  cat("Detected issues:\n")
  issues <- unique(issues)
  data.frame(sort(unique(issues$Issue)))
}

# (The counts of these raw transitions are in the biology assessment of the
# log, Section 11.)

# ========================================================================
# SECTION 6: STATUS PROPAGATION RULES
# ========================================================================
# Resolve placeholder "N" values so each stem history is fully defined
# (rstatus_propagate(), rstatus_functions.R):
#   - ^N → P, PN → PP  P until the stem's first record (not yet recruited),
#                      e.g. "NAAA" → "PAAA" (recruited in census 2)
#   - AN → AD          no record after alive: dead from the first census
#                      without a record (undone in Section 8 if the stem is
#                      recorded alive later), e.g. "AANN" → "AADD"
#   - DN → DD          no record after a dead record: stays dead
# No G exists yet (G is derived from tree context in Section 10), so there is
# no G propagation rule here.
cat("\n Patterns to propagate (before):\n")
tbl_sorted <- sort_table_status(new_status, sort_by = "x", decreasing = FALSE)
print(tbl_sorted[grepl("^N|PN|AN|DN", tbl_sorted$x), ])
new_status <- rstatus_propagate(new_status)
bio_check(
  !any(grepl("N", new_status, fixed = TRUE)),
  "No N (no-record placeholder) left after propagation"
)

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

## COMPARE BEFORE AND AFTER PROPAGATION
# Per stem: which raw history became which propagated history.
comparison_tbl <- data.table(
  before_propagation = original_status,
  after_propagation = new_status
)[, .(n_stems = .N), by = .(before_propagation, after_propagation)][order(-n_stems)]

cat(sprintf(
  "\n📋 Stems whose history changed during propagation: %d of %d\n",
  comparison_tbl[before_propagation != after_propagation, sum(n_stems)],
  length(new_status)
))
cat("\n📋 Most common changes (raw → propagated):\n")
print(head(comparison_tbl[before_propagation != after_propagation], 20))

###############################################################
# VALIDATION: stem histories after propagation (informational)
# -------------------------------------------------------------
# Checked against the legal stem transitions of the biological contract
# (valid_stem_trans, Section 1). DA is EXPECTED here (false deaths): it is
# resolved by Section 8 and enforced with bio_check() afterwards.
###############################################################
histories_after_propagation <- as.vector(tbl_sorted_after_propagation$code_after_propagation)
unique(unlist(strsplit(histories_after_propagation, "")))

issues <- check_histories(
  histories = histories_after_propagation,
  allowed_codes = status_codes,
  invalid_transitions = make_invalid_transitions(valid_stem_trans),
  nchars = length(ViewFullTable_split)
)

if (nrow(issues) == 0) {
  cat("No issues detected in any histories.\n")
} else {
  cat("Detected issues (expected before Section 8):\n")
  issues <- unique(issues)
  data.frame(sort(unique(issues$Issue)))
}

# ========================================================================
# SECTION 7: FIRST RECORDS — P ONLY BEFORE A STEM'S FIRST RECORD
# ========================================================================
# P means "not yet recorded": after Section 6 it only fills the censuses
# before a stem's first record (asserted below). A stem whose first record
# is dead is in the population from that record:
#   - alive later: the dead record was a false death, and the lifespan rule
#     (Section 8) makes it A;
#   - never alive: dead from its first record, G or D in Section 10.
# The DBH matrix built here is used by Sections 8, 10 and 13.

# ------------------------------------------------------------------------
# Build DBH matrix aligned with new_status / unique_StemID order
# ------------------------------------------------------------------------
# Rows = stems (in unique_StemID order), Columns = censuses (1..N),
# Values = DBH (numeric, NA when not measured).
DT_DBH <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(StemID, DBH, census = i)]
}), use.names = TRUE, fill = TRUE)
bio_check(
  DT_DBH[, .N, by = .(StemID, census)][, all(N == 1L)],
  "One DBH value per StemID × census before the pivot"
)
DBHs_dt <- dcast(DT_DBH, formula = StemID ~ census, value.var = "DBH")
DBHs_dt <- DBHs_dt[match(unique_StemID$StemID, DBHs_dt$StemID)]
DBHs <- as.matrix(DBHs_dt[, -1])

first_rec_dead <- grepl("^P*D", new_status)
has_A <- grepl("A", new_status, fixed = TRUE)
cat(sprintf(
  "🔎 Stems whose first record is dead: %d (alive later → A from that record: %d; never alive → G/D from that record: %d)\n",
  sum(first_rec_dead), sum(first_rec_dead & has_A), sum(first_rec_dead & !has_A)
))
cat(sprintf("🔎 Stems never recorded (P in every census): %d\n", sum(grepl("^P+$", new_status))))
p_after_record <- grep("[ADG]P", new_status)
bio_check(
  length(p_after_record) == 0L,
  "P only before a stem's first record (no P after A or D)",
  examples = unique_StemID$StemID[p_after_record],
  n_bad = length(p_after_record)
)
rm(first_rec_dead, has_A, p_after_record)

# ========================================================================
# SECTION 8: LIFESPAN RULE — RESOLVE RESURRECTIONS (D→A, G→A)
# ========================================================================
# Biological contract: a stem is alive from its first record to its last
# evidence of life (Section 4). Every A cell at this point is evidence of
# life, so any D/G followed later by an A was not a real death (a false
# death, including a dead first record of a stem alive later): it is
# backfilled to A (alive later => never dead), however long the gap. A DBH
# on such a false-death record is a real measurement and is kept.
# fix_resurrections() is therefore called with dbh_aware = FALSE (backfill
# all). Its DBH-aware mode, which turns an A without DBH that follows a D/G
# into D/G, is not used: it would contradict the contract.
# The dead records of stems alive later are listed in
# field_records_to_review_[site].csv (Section 12).

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
  dbh_aware  = FALSE, # lifespan rule: backfill every D/G followed by an A
  log_file   = NULL, # the cases go to the review list (Section 12)
  verbose    = TRUE
)
new_status <- RES_fix$status_vec

# ========================================================================
# SECTION 9: ASSERT — NO RESURRECTION, NO D/G BETWEEN TWO A's
# ========================================================================
# The lifespan rule in Section 8 leaves no D/G between two A cells; both
# conditions are asserted here.
resurrect_idx <- grep("DA|GA", new_status)
bio_check(
  length(resurrect_idx) == 0L,
  "No stem is alive after being dead (no DA/GA) after the lifespan rule",
  examples = unique_StemID$StemID[resurrect_idx],
  n_bad = length(resurrect_idx)
)
stranded_idx <- grep("A[DG]+A", new_status)
bio_check(
  length(stranded_idx) == 0L,
  "No D/G stranded between two A's after the lifespan rule",
  examples = unique_StemID$StemID[stranded_idx],
  n_bad = length(stranded_idx)
)

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

###############################################################
# VALIDATION: stem histories after Sections 7-9
# -------------------------------------------------------------
# Every history must now use only legal stem transitions
# (valid_stem_trans, Section 1). G does not exist yet: it is derived from
# tree context in Section 10.
###############################################################
histories_correction_after_propagation <- as.vector(tbl_sorted_correction_after_propagation$code_correction_after_propagation)

issues <- check_histories(
  histories = histories_correction_after_propagation,
  allowed_codes = status_codes,
  invalid_transitions = make_invalid_transitions(valid_stem_trans),
  nchars = length(ViewFullTable_split)
)
bio_check(
  nrow(issues) == 0L,
  "Stem histories after Sections 7-9 use only legal stem transitions",
  examples = unique(issues$History),
  n_bad = length(unique(issues$History))
)

# ========================================================================
# SECTION 10: TREE-LEVEL STATUS CALCULATION AND D/G CORRECTION
# ========================================================================
# Aggregate stem histories by tree and assign D/G from the tree's life:
#   tree_last_A = last census in which ANY stem of the tree is A.
#   A dead cell (D or G) at census j becomes
#     G  if j <= tree_last_A  (the tree is alive now or later: stem dead,
#                              tree alive — possibly unregistered)
#     D  if j >  tree_last_A  (no stem of the tree lives again: tree dead)
# Consequences (biological contract, see header):
#   - GD is a valid transition (the tree dies after the stem did);
#   - DG is NOT valid: D means the whole tree is dead, and a dead tree
#     cannot come back to life;
#   - a stem never alive (dead from its first record) has no A, so it never
#     keeps its tree alive;
#   - single-stem trees get D directly (tree_last_A = the stem's own).
# A and P cells are never modified.

#--------------------------------------------------------------
# 10.1  Tag <-> TreeID consistency check
#--------------------------------------------------------------
# Each Tag should map to exactly one TreeID (and vice versa). Anything
# else means the master stem table has duplicate identifiers and tree
# grouping below would be wrong.
tag_treeid_dt <- unique(unique_StemID[, .(Tag, TreeID)])
tags_per_treeid <- tag_treeid_dt[, .N, by = TreeID][N > 1L]
treeids_per_tag <- tag_treeid_dt[, .N, by = Tag][N > 1L]

# Tree-level D/G below groups stems by TreeID, so the grouping must be sound.
bio_check(
  nrow(tags_per_treeid) == 0L && nrow(treeids_per_tag) == 0L,
  "Tag <-> TreeID mapping is 1:1",
  examples = c(tags_per_treeid$TreeID, treeids_per_tag$Tag),
  n_bad = nrow(tags_per_treeid) + nrow(treeids_per_tag)
)

# A tree is one individual, so all of its stems must be the same species;
# otherwise stems of different plants were grouped into one tree.
multi_species <- unique_StemID[, .(n_sp = uniqueN(Mnemonic)), by = TreeID][n_sp > 1L]
bio_check(
  nrow(multi_species) == 0L,
  "Every tree has a single species across all its stems",
  examples = multi_species$TreeID,
  n_bad = nrow(multi_species)
)
rm(multi_species)

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
# TREE-LEVEL ANALYSIS: HELPER FUNCTIONS AND VALIDATION
# ========================================================================
# This section defines functions to compute tree-level status from stem data
# and validates biological plausibility of life-history sequences.

# Legal transitions: valid_stem_trans / valid_tree_trans (Section 1).

# --------
# Function: tree_state - Compute tree-level status from stem matrix
# --------
# stem_states: character matrix, rows = censuses, columns = stems of one tree.
# Returns the tree history as one string (one character per census).
# KEY DECISION POINTS:
#   1. If ANY stem is "A" → tree is "A"
#   2. If all stems are "P" → tree is "P"
#   3. If mix of P and D/G:
#      - Check future censuses for any "A"
#      - If future A exists → tree is "A" (alive but unregistered: a stem
#        of the tree is recorded alive later)
#      - If no future A → tree is "D" (no stem of the tree lives again;
#        P stems here are recorded later only as dead, or never)
#   4. If only D/G (no P, no A) → tree is "D"
# This per-tree loop is an independent implementation of the contract; it is
# compared cell by cell with the vectorized Section 10.5 result (check I4).
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
        # P + D/G and no future A → flagged. Under the biological contract
        # this is simply a DEAD tree whose P stems are never alive (recorded
        # later only as dead, or never recorded), so the check is NOT
        # enforced inside compute_tree_for_row(); it only feeds the
        # informational RUN_TREE_EXISTS_DIAGNOSTIC count below.
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
  # When enabled it returns NA for any tree with P + D/G and no future A.
  # Under the contract those are simply dead trees holding P stems that are
  # never alive, so tree_state() codes them D, which is correct.
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

# Optional diagnostic: re-run with the strict tree_exists_check() guard and
# count the trees it would reject. Those are trees with P + D stems and no
# future A: under the contract they are simply DEAD trees whose P stems are
# never alive, so the count is informational (dead trees that hold stems
# recorded later only as dead, or never), not an error. Set to FALSE to skip
# (it re-runs the tree loop).
RUN_TREE_EXISTS_DIAGNOSTIC <- TRUE
if (RUN_TREE_EXISTS_DIAGNOSTIC) {
  tree_histories_list_checking <- lapply(new_status_split, compute_tree_for_row_checking)
  diff_indices <- which(!mapply(identical, tree_histories_list, tree_histories_list_checking))
  cat(sprintf(
    "  diagnostic: tree_exists_check would change %d / %d tree histories.\n",
    length(diff_indices), length(tree_histories_list)
  ))
  rm(tree_histories_list_checking)
}

# Validate results: every tree history must use only legal tree transitions.
# A tree never recorded stays all-P; a tree recorded only as dead is P until
# its first record and D from then on (PD is legal).
tree_hist_vec <- unlist(tree_histories_list, use.names = TRUE)
tree_pairs_ok <- Reduce(`&`, lapply(
  seq_len(length(ViewFullTable_split) - 1L),
  function(k) substr(tree_hist_vec, k, k + 1L) %in% valid_tree_trans
))
bio_check(
  all(tree_pairs_ok),
  sprintf("Every tree history uses only legal tree transitions (%s)", paste(valid_tree_trans, collapse = " ")),
  examples = names(tree_hist_vec)[!tree_pairs_ok],
  n_bad = sum(!tree_pairs_ok)
)
cat(sprintf(
  "   trees never alive: %d of %d (never recorded, all P: %d)\n",
  sum(!grepl("A", tree_hist_vec, fixed = TRUE)), length(tree_hist_vec), sum(grepl("^P+$", tree_hist_vec))
))

# Match tree status to each stem's TreeID
tree_histories <- tree_histories_list[as.character(unique_StemID$TreeID)]
# Convert list to matrix (rows = stems, columns = censuses)
tree_histories <- do.call(rbind, tree_histories)
cat("✓ Tree status matched to stem level\n\n")

# ========================================================================
# 10.4  BUILD STATUS MATRIX AND APPLY THE TREE-AWARE D/G RULE
# ========================================================================
# Convert the per-stem encounter history strings into a character matrix
# (rows = stems, columns = censuses), then apply the tree-aware D/G rule
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
# 10.5  APPLY THE TREE-AWARE D/G RULE (vectorized)
# ------------------------------------------------------------------
# For every stem : stem_last_A = last census in which the stem is A
#                  (0 if the stem is never alive)
# For every tree : tree_last_A = max(stem_last_A) over its stems
# Every dead cell (D or G) at census j becomes
#   G  if j <= tree_last_A   the tree is alive now or later (possibly
#                            unregistered): stem dead, tree alive
#   D  if j >  tree_last_A   no stem of the tree lives again: tree dead
# A stem never alive has stem_last_A = 0, so it can never keep a dead tree
# alive. Single-stem trees get D directly. A and P cells are NEVER touched
# (checked below). The rule itself is rstatus_tree_dg() (rstatus_functions.R);
# stem_info_dt below recomputes the tree life span for the checks.
# ------------------------------------------------------------------
bio_check(
  nrow(new_status_matrix) == nrow(unique_StemID),
  "Status matrix has one row per master stem"
)

n_cens <- ncol(new_status_matrix)
n_stems <- nrow(new_status_matrix)
TreeID_vec <- unique_StemID$TreeID

# Per-stem row info (rows stay in unique_StemID order).
stem_info_dt <- data.table(row_idx = seq_len(n_stems), TreeID = TreeID_vec)
stem_info_dt[, n_stems_in_tree := .N, by = TreeID]
single_stem_row <- stem_info_dt$n_stems_in_tree == 1L

# First / last census in which each stem is A. A stem never alive gets
# first = n_cens + 1 ("never") and last = 0, so both stay integers.
never_cens <- n_cens + 1L
is_A_num <- (new_status_matrix == "A") * 1L
has_A <- rowSums(is_A_num) > 0L
stem_info_dt[, `:=`(
  stem_first_A = ifelse(has_A, max.col(is_A_num, ties.method = "first"), never_cens),
  stem_last_A = ifelse(has_A, max.col(is_A_num, ties.method = "last"), 0L)
)]
# Tree life span: first and last census in which ANY stem of the tree is A
# (tree_first_A = n_cens + 1 and tree_last_A = 0 for trees never alive).
stem_info_dt[, `:=`(
  tree_first_A = min(stem_first_A),
  tree_last_A = max(stem_last_A)
), by = TreeID]

col_idx <- matrix(seq_len(n_cens), nrow = n_stems, ncol = n_cens, byrow = TRUE)
tree_last_A_mat <- matrix(stem_info_dt$tree_last_A, nrow = n_stems, ncol = n_cens)
tree_alive_cell <- col_idx <= tree_last_A_mat
dead_cell <- new_status_matrix == "D" | new_status_matrix == "G"

corrected_new_status_matrix <- rstatus_tree_dg(new_status_matrix, TreeID_vec)
bio_check(
  all(corrected_new_status_matrix[dead_cell & tree_alive_cell] == "G") &&
    all(corrected_new_status_matrix[dead_cell & !tree_alive_cell] == "D"),
  "rstatus_tree_dg() gives G when the tree is alive then or later and D otherwise"
)
n_to_G <- sum(dead_cell & tree_alive_cell & new_status_matrix != "G")
n_to_D <- sum(dead_cell & !tree_alive_cell & new_status_matrix != "D")

# A and P cells must be identical between input and output.
AP_in <- new_status_matrix == "A" | new_status_matrix == "P"
AP_out <- corrected_new_status_matrix == "A" | corrected_new_status_matrix == "P"
bio_check(
  identical(AP_in, AP_out) &&
    identical(new_status_matrix[AP_in], corrected_new_status_matrix[AP_in]),
  "A and P cells are unchanged by the tree-aware D/G rule"
)

# Reassemble per-stem string vector.
corrected_new_status <- do.call(
  paste0,
  lapply(seq_len(n_cens), function(j) corrected_new_status_matrix[, j])
)
names(corrected_new_status) <- rownames(new_status_matrix)

cat(sprintf(
  "✓ Section 10 tree-aware D/G: %d cells changed D->G (tree alive), %d changed G->D (tree dead); A/P preserved.\n",
  n_to_G, n_to_D
))

# ---- Tree-level invariants of the biological contract -------------------
cat("\n\U0001F50D Validating tree-level biology...\n")
is_D_num <- (corrected_new_status_matrix == "D") * 1L
has_D <- rowSums(is_D_num) > 0L
stem_info_dt[, stem_first_D := ifelse(has_D, max.col(is_D_num, ties.method = "first"), never_cens)]
stem_info_dt[, tree_first_D := min(stem_first_D), by = TreeID]

# I1: no stem of a tree is A at or after the first census where the tree is D.
i1_bad <- stem_info_dt[tree_first_D < never_cens & tree_last_A >= tree_first_D, unique(TreeID)]
bio_check(
  length(i1_bad) == 0L,
  "I1: no stem is alive at or after its tree is dead (no tree-level zombies)",
  examples = i1_bad,
  n_bad = length(i1_bad)
)

# I2: every G cell belongs to a tree with an A at the same or a later census.
i2_bad <- which(rowSums(corrected_new_status_matrix == "G" & !tree_alive_cell) > 0L)
bio_check(
  length(i2_bad) == 0L,
  "I2: every G (dead stem) belongs to a tree alive at that census or later",
  examples = unique_StemID$StemID[i2_bad],
  n_bad = length(i2_bad)
)

# I3: single-stem trees have no G (a single dead stem means a dead tree).
i3_bad <- which(single_stem_row & rowSums(corrected_new_status_matrix == "G") > 0L)
bio_check(
  length(i3_bad) == 0L,
  "I3: single-stem trees have no G",
  examples = unique_StemID$StemID[i3_bad],
  n_bad = length(i3_bad)
)

# I4: the tree status implied by the corrected stems (P before the tree's
# first record, A from then to its last A, D after) must equal
# tree_histories, computed independently by tree_state() in Section 10.3.
is_rec_num <- (corrected_new_status_matrix != "P") * 1L
stem_info_dt[, stem_first_rec := ifelse(rowSums(is_rec_num) > 0L, max.col(is_rec_num, ties.method = "first"), never_cens)]
stem_info_dt[, tree_first_rec := min(stem_first_rec), by = TreeID]
tree_first_rec_mat <- matrix(stem_info_dt$tree_first_rec, nrow = n_stems, ncol = n_cens)
tree_derived <- matrix("P", nrow = n_stems, ncol = n_cens)
tree_derived[col_idx >= tree_first_rec_mat & col_idx <= tree_last_A_mat] <- "A"
tree_derived[col_idx >= tree_first_rec_mat & col_idx > tree_last_A_mat] <- "D"
tree_derived_vec <- do.call(paste0, lapply(seq_len(n_cens), function(j) tree_derived[, j]))
i4_bad <- which(tree_derived_vec != as.vector(tree_histories))
bio_check(
  length(i4_bad) == 0L,
  "I4: tree status implied by the stem codes equals tree_histories (tree_state)",
  examples = unique_StemID$TreeID[i4_bad],
  n_bad = uniqueN(unique_StemID$TreeID[i4_bad])
)

# I5: a P cell (not yet in the population) never has a DBH.
i5_bad <- which(rowSums(corrected_new_status_matrix == "P" & !is.na(DBHs)) > 0L)
bio_check(
  length(i5_bad) == 0L,
  "I5: P (prior) cells never carry a DBH",
  examples = unique_StemID$StemID[i5_bad],
  n_bad = length(i5_bad)
)

rm(
  stem_info_dt, AP_in, AP_out, single_stem_row, TreeID_vec,
  is_A_num, has_A, is_D_num, has_D, col_idx, tree_first_rec_mat, is_rec_num,
  tree_last_A_mat, tree_alive_cell, dead_cell, tree_derived,
  tree_derived_vec, n_to_G, n_to_D, i1_bad, i2_bad, i3_bad, i4_bad, i5_bad,
  never_cens
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

# assess_biology(): transition counts of a vector of histories.
#   status_vec  : character vector, one history per stem (or tree)
#   label       : name of this version, copied to the result
#   valid_trans : legal 2-letter transitions
# Returns list(label, n_stems, n_bad_rows = histories with an illegal
# transition, n_stranded = histories with D/G between two A, illegal_tab and
# legal_tab = counts per transition).
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
  list(label = "1. original_status      (raw, pre-propagation)", vec = original_status, trans = valid_stem_trans, enforce = FALSE),
  list(label = "2. new_status           (after Sections 7-9)", vec = new_status, trans = valid_stem_trans, enforce = TRUE),
  list(label = "3. corrected_new_status (after Section 10.5 tree-aware D/G)", vec = corrected_new_status, trans = valid_stem_trans, enforce = TRUE),
  list(label = "4. tree_histories       (tree-level, one per stem)", vec = as.vector(tree_histories), trans = valid_tree_trans, enforce = TRUE)
)

# Tree × census composition of the exported stem codes (reused in Section 12).
tree_cens <- data.table(
  TreeID = rep(unique_StemID$TreeID, times = ncol(corrected_new_status_matrix)),
  census = rep(seq_len(ncol(corrected_new_status_matrix)), each = nrow(corrected_new_status_matrix)),
  s      = as.vector(corrected_new_status_matrix)
)[, .(nA = sum(s == "A"), nG = sum(s == "G"), nD = sum(s == "D"), nP = sum(s == "P")), by = .(TreeID, census)]
n_unregistered <- tree_cens[nG > 0L & nA == 0L, .N]
n_tree_zombie_cells <- tree_cens[nD > 0L & (nA > 0L | nG > 0L), .N]

# The assessment goes to the log (check_[site].txt).
fmt_counts <- function(tab) {
  if (length(tab) == 0L) "none" else paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", ")
}
bio_lines <- c(
  "Biology assessment (transitions between consecutive censuses in each version of the histories):",
  paste0("  legal stem transitions: ", paste(valid_stem_trans, collapse = " "), " | legal tree transitions: ", paste(valid_tree_trans, collapse = " ")),
  "  G = stem dead, tree alive then or later; D = tree dead (absorbing)"
)
for (v in versions) {
  rep_v <- assess_biology(v$vec, v$label, v$trans)
  bio_lines <- c(
    bio_lines,
    sprintf(
      "  %s: %d histories | with an illegal transition: %d | with D/G between two A: %d",
      trimws(rep_v$label), rep_v$n_stems, rep_v$n_bad_rows, rep_v$n_stranded
    ),
    sprintf("      illegal: %s", fmt_counts(rep_v$illegal_tab)),
    sprintf("      legal  : %s", fmt_counts(rep_v$legal_tab))
  )
}
bio_lines <- c(
  bio_lines,
  sprintf("  tree-census cells with a D stem next to an A or G stem (must be 0): %d", n_tree_zombie_cells),
  sprintf("  tree-census cells with G stems and no A stem (tree alive but unregistered: a stem lives later): %d", n_unregistered),
  "  Expected: version 1 (raw) may hold illegal transitions; versions 2-4 hold none."
)
print_to_log(bio_lines, log_file, new_message = TRUE)
rm(bio_lines, fmt_counts)

# Print the most important lines to stdout so the user sees pass/fail
# without opening the file, then enforce the contract.
cat("\nBiology assessment summary:\n")
for (v in versions) {
  r <- assess_biology(v$vec, v$label, v$trans)
  cat(sprintf("  %s\n      illegal rows = %d, A[DG]+A rows = %d\n", v$label, r$n_bad_rows, r$n_stranded))
  if (v$enforce) {
    bio_check(
      r$n_bad_rows == 0L && r$n_stranded == 0L,
      sprintf("%s has no illegal transitions", trimws(v$label)),
      examples = names(r$illegal_tab),
      n_bad = r$n_bad_rows
    )
  }
}
bio_check(
  n_tree_zombie_cells == 0L,
  "No tree-census has a D stem next to an A or G stem",
  n_bad = n_tree_zombie_cells
)
cat(sprintf("   tree-census cells with unregistered survival (accepted): %d\n", n_unregistered))
rm(versions)

# ========================================================================
# SECTION 12: STATUS PATTERNS AND REVIEW LISTS
# ========================================================================
# status_patterns_[site].csv: one row per encounter history (raw codes, DBH
# recorded or not, exported Rstatus, tree status) with its number of stems
# and one example stem.
# Two review lists, one row per case, written at the end of Section 14:
#   field_records_to_review_[site].csv   raw records that need a look (for the
#                                        data manager)
#   reconstruction_to_review_[site].csv  stems whose reconstruction needs a
#                                        look (stem identity, stage 2)
# Columns: reason, treeID, tag, stemID (the exported stem number; NA for a
# tree-level case), census, dbh, DFstatus and Rstatus (of that stem in that
# census) and details. Nothing is changed here; the counts per reason are at
# the end of the log.

cat("\n📋 Status patterns and review lists...\n")

# ---- Status patterns ---------------------------------------------------------
# raw_history: A evidence of life, D dead record, N no record (Section 4);
# dbh_recorded: 1 where a DBH was recorded in that census.
dbh_recorded <- do.call(paste0, lapply(seq_len(ncol(DBHs)), function(j) as.integer(!is.na(DBHs[, j]))))
status_patterns <- data.table(
  raw_history = original_status,
  dbh_recorded = dbh_recorded,
  Rstatus_history = corrected_new_status,
  tree_history = as.vector(tree_histories),
  treeID = unique_StemID$TreeID,
  stemID = as.integer(sub(".*_", "", unique_StemID$StemID))
)[, .(n_stems = .N, example_treeID = treeID[1L], example_stemID = stemID[1L]),
  by = .(raw_history, dbh_recorded, Rstatus_history, tree_history)
]
setorder(status_patterns, -n_stems, raw_history)
fwrite(status_patterns, file.path(CHECK_folder, paste0("status_patterns_", site, ".csv")))
cat(sprintf("  ✓ Saved: status_patterns_%s.csv (%d patterns)\n", site, nrow(status_patterns)))
rm(dbh_recorded, status_patterns)

# ---- Helpers for the review lists ---------------------------------------------
# Raw Status per stem x census, rows in master stem order (as DBHs).
Status_mat <- do.call(cbind, lapply(ViewFullTable_split, function(x) as.character(x$Status)))
bio_check(
  identical(dim(Status_mat), dim(corrected_new_status_matrix)),
  "Raw Status matrix has one row per master stem and one column per census"
)
is_status <- function(x) !is.na(Status_mat) & Status_mat == x
# Column of a CensusID in the census tables
census_index <- function(census_id) {
  match(as.numeric(census_id), vapply(ViewFullTable_split, function(x) as.numeric(x$CensusID[1L]), numeric(1)))
}
stem_number <- function(stem_key) sub(".*_", "", stem_key)
# Rows of a review list: stems (keys "TreeID_<stem number>", with the census
# of each case) or, with `tree`, whole trees. NULL when there is no case.
review_rows <- function(reason, stem_key = NULL, census = NA_integer_, details = "", tree = NULL) {
  n <- if (is.null(tree)) length(stem_key) else length(tree)
  if (n == 0L) {
    return(NULL)
  }
  census <- rep_len(as.integer(census), n)
  if (is.null(tree)) {
    i <- match(stem_key, unique_StemID$StemID)
    ij <- cbind(i, fcoalesce(census, 1L))
    has_c <- !is.na(census)
    data.table(
      reason = reason, treeID = unique_StemID$TreeID[i], tag = unique_StemID$Tag[i],
      stemID = as.integer(stem_number(stem_key)), census = census,
      dbh = fifelse(has_c, as.numeric(DBHs[ij]), NA_real_),
      DFstatus = fifelse(has_c, Status_mat[ij], NA_character_),
      Rstatus = fifelse(has_c, corrected_new_status_matrix[ij], NA_character_),
      details = rep_len(details, n)
    )
  } else {
    data.table(
      reason = reason, treeID = tree, tag = unique_StemID$Tag[match(tree, unique_StemID$TreeID)],
      stemID = NA_integer_, census = census,
      dbh = NA_real_, DFstatus = NA_character_, Rstatus = NA_character_,
      details = rep_len(details, n)
    )
  }
}
field_review <- list() # raw records (data manager)
recon_review <- list() # stem reconstruction

# ---- Dead records and stumps of stems alive later ------------------------------
# A "dead" / "stem dead" record of a stem recorded alive later is a false
# death: the lifespan rule (Section 8) exports it as A. A "broken below"
# record without DBH of a stem alive later is a stump: A, with no measurable
# DBH at 1.3 m (stage 4 keeps it out of the stock); it is A because stage 2
# put it on a stem that lives later (e.g. a field-tagged stump).
false_dead <- which((is_status("dead") | is_status("stem dead")) & corrected_new_status_matrix == "A", arr.ind = TRUE)
field_review[["dead_record_of_stem_alive_later"]] <- review_rows(
  "dead_record_of_stem_alive_later", unique_StemID$StemID[false_dead[, 1]], false_dead[, 2],
  paste("raw codes", original_status[false_dead[, 1]])
)
stumps <- which(is_status("broken below") & is.na(DBHs) & corrected_new_status_matrix == "A", arr.ind = TRUE)
recon_review[["stump_of_stem_alive_later"]] <- review_rows(
  "stump_of_stem_alive_later", unique_StemID$StemID[stumps[, 1]], stumps[, 2],
  paste("raw codes", original_status[stumps[, 1]])
)
cat(sprintf(
  "  Dead records of stems alive later: %d | stumps of stems alive later: %d\n",
  nrow(false_dead), nrow(stumps)
))

# ---- Stems never alive and stems never recorded --------------------------------
# Never alive: recorded only as dead (e.g. a broken-below record without DBH,
# then dead); P until the first record, then G or D; they never keep a tree
# alive. Never recorded: P in every census (only "missing" records without
# DBH). One row per stem, at its first record.
never_alive <- which(!grepl("A", corrected_new_status, fixed = TRUE) & grepl("[GD]", corrected_new_status))
recon_review[["stem_never_alive"]] <- review_rows(
  "stem_never_alive", unique_StemID$StemID[never_alive], as.integer(regexpr("[GD]", corrected_new_status[never_alive])),
  paste("raw codes", original_status[never_alive])
)
never_recorded <- which(grepl("^P+$", corrected_new_status))
recon_review[["stem_never_recorded"]] <- review_rows(
  "stem_never_recorded", unique_StemID$StemID[never_recorded], NA_integer_,
  paste("raw codes", original_status[never_recorded])
)
cat(sprintf("  Stems never alive: %d | never recorded (P in every census): %d\n", length(never_alive), length(never_recorded)))

# ------------------------------------------------------------------------
# Long stem × census table of exported codes and raw evidence (deaths and
# new stems). Raw Status in ViewFullTable_split is the field record (NA = no
# row in that census); rows are in unique_StemID order.
# ------------------------------------------------------------------------
diag_long <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(
    TreeID, StemID, Raw_StemID,
    census = i,
    raw_status = Status,
    DBH,
    ReconstructionMethod,
    DP_PosteriorReconstructedProb,
    Rstatus = corrected_new_status_matrix[, i]
  )]
}))
setkey(diag_long, StemID, census)
diag_long[, `:=`(
  prev_Rstatus = shift(Rstatus),
  prev_DBH = shift(DBH),
  prev_Raw_StemID = shift(Raw_StemID)
), by = StemID]
diag_long[, n_stems_in_tree := uniqueN(StemID), by = TreeID]

## Evidence behind each stem death, per census (to the log) ####
# First dead cell (G or D) after A, cross-tabulated with the raw field
# status of that census. "<no record>" = the death is inferred from the
# absence of any row (older censuses did not record dead secondary stems).
deaths_dt <- diag_long[Rstatus %in% c("G", "D") & prev_Rstatus == "A"]
death_evidence <- dcast(
  deaths_dt[, .N, by = .(raw_status = fifelse(is.na(raw_status), "<no record>", raw_status), census)],
  raw_status ~ census,
  value.var = "N", fill = 0L
)
cat("\n📋 Evidence behind each stem death (raw status in the death census):\n")
print(death_evidence)
print_to_log("Evidence behind each stem death (raw status in the death census):", log_file, new_message = TRUE)
print_to_log(log_table(death_evidence), log_file, new_message = FALSE)

## Trees alive but unregistered (count only) ####
# Tree × census cells where the tree has dead stems (G) but no living stem
# (A): a stem of the tree is recorded alive in a LATER census, so by the
# contract the tree was alive but not registered. Accepted biology.
unreg_txt <- sprintf(
  "Trees alive but unregistered (G stems, no A stem; a stem lives later): %d tree-census cells in %d trees",
  tree_cens[nG > 0L & nA == 0L, .N], tree_cens[nG > 0L & nA == 0L, uniqueN(TreeID)]
)
cat(" ", unreg_txt, "\n")
print_to_log(unreg_txt, log_file, new_message = TRUE)

## New stems: identity-break candidates and new stems of 50 mm or more ####
# New stem = P -> A. On a tree recorded in an earlier census it is ingrowth,
# otherwise a new tree (recruitment).
tree_first <- tree_cens[nA + nG + nD > 0L, .(tree_first = min(census)), by = TreeID]
new_stems <- tree_first[diag_long[Rstatus == "A" & prev_Rstatus == "P"], on = "TreeID"]
new_stems[, kind := fifelse(census > tree_first, "ingrowth (tree recorded before)", "new tree (recruitment)")]
# Candidate identity breaks: same tree and census, exactly one stem dies and
# exactly one new stem appears at a similar size (new stem >= 20 mm, 0.8-1.5 x
# the dead stem's last DBH). If both are the same physical stem, the tables
# hold a false death plus a false new stem. "same database StemID": the
# database kept one StemID across the two (strong evidence of a split);
# "different": the database also has two stems (resprout after breakage, or
# renumbering).
deaths_1 <- deaths_dt[, .(
  n_deaths = .N,
  dead_StemID = StemID[1L],
  dead_last_DBH = prev_DBH[1L],
  dead_last_Raw_StemID = prev_Raw_StemID[1L]
), by = .(TreeID, census)]
new_1 <- new_stems[, .(
  n_new = .N,
  new_StemID = StemID[1L],
  new_DBH = DBH[1L],
  new_Raw_StemID = Raw_StemID[1L],
  new_method = ReconstructionMethod[1L],
  new_posterior_prob = DP_PosteriorReconstructedProb[1L]
), by = .(TreeID, census)]
id_breaks <- merge(deaths_1, new_1, by = c("TreeID", "census"))[n_deaths == 1L & n_new == 1L]
id_breaks[, size_ratio := new_DBH / dead_last_DBH]
id_breaks <- id_breaks[new_DBH >= 20 & size_ratio >= 0.8 & size_ratio <= 1.5]
id_breaks[, same_db := !is.na(new_Raw_StemID) & !is.na(dead_last_Raw_StemID) & new_Raw_StemID == dead_last_Raw_StemID]
recon_review[["identity_break_candidate"]] <- review_rows(
  "identity_break_candidate", id_breaks$new_StemID, id_breaks$census,
  sprintf(
    "new stem the size of stem %s, dead in this census (last DBH %s mm, ratio %.2f, %s database StemID); method %s, posterior %.2f",
    stem_number(id_breaks$dead_StemID), id_breaks$dead_last_DBH, id_breaks$size_ratio,
    fifelse(id_breaks$same_db, "same", "different"), id_breaks$new_method, id_breaks$new_posterior_prob
  )
)
# New stems of 50 mm or more at their first measurement (a new stem should be
# close to the 10 mm threshold): often identity breaks or resprouts.
large_new <- new_stems[DBH >= 50]
recon_review[["new_stem_50mm_or_more"]] <- review_rows(
  "new_stem_50mm_or_more", large_new$StemID, large_new$census,
  sprintf(
    "%s; first DBH %s mm; tree with %d stems; method %s, posterior %.2f",
    large_new$kind, large_new$DBH, large_new$n_stems_in_tree, large_new$ReconstructionMethod, large_new$DP_PosteriorReconstructedProb
  )
)
cat(sprintf(
  "  Candidate identity breaks: %d (same database StemID: %d) | new stems >= 50 mm: %d (>= 100 mm: %d)\n",
  nrow(id_breaks), id_breaks[same_db == TRUE, .N], nrow(large_new), large_new[DBH >= 100, .N]
))

## Raw "stem dead" exported as D ####
# "stem dead" implies the tree was still alive, but no other stem of the
# tree is alive in this or any later census, so the tree is exported as
# dead (a field-code conflict).
n_in_tree <- unique_StemID[, .N, by = TreeID][unique_StemID, on = "TreeID", N]
sd_D <- which(is_status("stem dead") & corrected_new_status_matrix == "D", arr.ind = TRUE)
field_review[["stem_dead_but_tree_dead"]] <- review_rows(
  "stem_dead_but_tree_dead", unique_StemID$StemID[sd_D[, 1]], sd_D[, 2],
  sprintf("raw 'stem dead', but no stem of the tree is alive then or later; tree with %d stems", n_in_tree[sd_D[, 1]])
)
cat(sprintf("  Raw 'stem dead' exported as D (tree dead): %d cells\n", nrow(sd_D)))

## Likely duplicate measurements ####
# When the point of measurement (POM) of a trunk was raised, the crew could
# measure the trunk at the old and at the new height in the same census, and
# the database kept the two records as two stems. Signature: two stems of one
# tree, both >= 10 cm, measured on the same date of one census at heights of
# measurement at least 0.5 m apart, with taper-corrected DBHs within 20 % of
# each other. Two similar stems measured at the same height cannot be told
# from a duplicate, and 1982 has no HOM records, so neither is listed. One row
# per pair, on the stem measured lower.
dup_obs <- ViewFullTable[!is.na(DBH) & DBH >= 100 & !is.na(ExactDate), .(
  TreeID, Tag,
  census = CensusID, ExactDate, StemID, Raw_StemID, DBH,
  HOM = fcoalesce(as.numeric(HOM), 1.3), tc = dbh_with_best_candidate_taper_corrected
)]
dup_span <- ViewFullTable[!is.na(DBH), .(first_census = min(CensusID), last_census = max(CensusID)), by = StemID]
dup_tree_first <- ViewFullTable[!is.na(DBH), .(tree_first = min(CensusID)), by = TreeID]
dup <- dup_obs[dup_obs, on = .(TreeID, census, ExactDate), allow.cartesian = TRUE, nomatch = 0L][
  i.HOM - HOM >= 0.5 & pmin(tc, i.tc) / pmax(tc, i.tc) >= 0.8
]
dup <- dup[, .(
  TreeID, census,
  StemID_low = StemID, HOM_low = HOM,
  StemID_high = i.StemID, DBH_high = i.DBH, HOM_high = i.HOM,
  DBH_low = DBH, dbh_ratio = round(pmin(tc, i.tc) / pmax(tc, i.tc), 3)
)]
dup <- dup_span[, .(StemID_low = StemID, first_low = first_census, last_low = last_census)][dup, on = "StemID_low"]
dup <- dup_span[, .(StemID_high = StemID, first_high = first_census, last_high = last_census)][dup, on = "StemID_high"]
dup <- dup_tree_first[dup, on = "TreeID"]
dup[, one_series_starts_or_ends_here := (first_low == census) != (first_high == census) | (last_low == census) != (last_high == census)]
dup[, impossible_new_stem := census > tree_first & ((first_low == census & DBH_low >= 260) | (first_high == census & DBH_high >= 260))]
field_review[["duplicate_measurement"]] <- review_rows(
  "duplicate_measurement", dup$StemID_low, census_index(dup$census),
  sprintf(
    "same date as stem %s (HOM %.2f m, DBH %s mm); this stem HOM %.2f m; taper-corrected DBH ratio %.2f%s%s",
    stem_number(dup$StemID_high), dup$HOM_high, dup$DBH_high, dup$HOM_low, dup$dbh_ratio,
    fifelse(dup$one_series_starts_or_ends_here, "; one of the two series starts or ends here", ""),
    fifelse(dup$impossible_new_stem, "; one of the two would be a new stem >= 260 mm on a tree measured before", "")
  )
)
cat(sprintf(
  "  Likely duplicate measurements: %d pairs in %d trees\n",
  nrow(dup), uniqueN(dup$TreeID)
))

rm(
  false_dead, stumps, never_alive, never_recorded, diag_long, deaths_dt, death_evidence, unreg_txt,
  tree_first, new_stems, deaths_1, new_1, id_breaks, large_new, n_in_tree, sd_D,
  dup_obs, dup_span, dup_tree_first, dup
)

# ========================================================================
# SECTION 13: STATUS × DBH SUPPORT SUMMARY
# ========================================================================
# Audit the final corrected_new_status_matrix before export: verify that
# no illegal transitions remain and that A/P cells are preserved; list the
# DBHs recorded on dead (G/D) records (field review list); summarise status
# codes by whether a DBH was recorded (log). The exported dbh is the raw DBH
# (DBHs), unchanged.

bio_check(
  is.matrix(corrected_new_status_matrix) &&
    nrow(corrected_new_status_matrix) == nrow(DBHs) &&
    ncol(corrected_new_status_matrix) == ncol(DBHs),
  "Final status matrix and DBH matrix have the same shape"
)

# ---- (a) and (b) final biology + A/P-preservation audit ------------------
final_pairs <- substring(
  rep(corrected_new_status, each = ncol(corrected_new_status_matrix) - 1L),
  rep(seq_len(ncol(corrected_new_status_matrix) - 1L), length(corrected_new_status)),
  rep(seq_len(ncol(corrected_new_status_matrix) - 1L), length(corrected_new_status)) + 1L
)
illegal <- setdiff(unique(final_pairs), valid_stem_trans)
bio_check(
  length(illegal) == 0L,
  "Final stem histories contain only legal stem transitions (counts in the biology assessment of the log)",
  examples = illegal,
  n_bad = length(illegal)
)

stranded <- which(grepl("A[DG]+A", corrected_new_status))
bio_check(
  length(stranded) == 0L,
  "No D/G stranded between two A's in the final stem histories",
  examples = unique_StemID$StemID[stranded],
  n_bad = length(stranded)
)

# A/P cells must match between new_status_matrix and corrected_new_status_matrix.
AP_in <- new_status_matrix == "A" | new_status_matrix == "P"
AP_out <- corrected_new_status_matrix == "A" | corrected_new_status_matrix == "P"
bio_check(
  identical(AP_in, AP_out) &&
    identical(new_status_matrix[AP_in], corrected_new_status_matrix[AP_in]),
  "A and P cells preserved between new_status_matrix and corrected_new_status_matrix"
)

rm(final_pairs, illegal, stranded, AP_in, AP_out)

# ---- DBH: exported exactly as recorded --------------------------------------
# No DBH is removed or imputed: the exported dbh is the field record. A DBH
# sits on a G/D cell only when it was recorded on a "dead" / "stem dead"
# record of a stem that is never alive later; those records go to the field
# review list (reason dbh_on_dead_record). The DBH of a false-death record
# (stem alive later) sits on an A cell. A P cell has no record, so no DBH.
dd_idx <- which((corrected_new_status_matrix == "G" | corrected_new_status_matrix == "D") & !is.na(DBHs), arr.ind = TRUE)
dbh_dead_dt <- data.table(
  row = dd_idx[, 1], census = dd_idx[, 2],
  TreeID = unique_StemID$TreeID[dd_idx[, 1]],
  StemID = unique_StemID$StemID[dd_idx[, 1]],
  DBH = DBHs[dd_idx],
  Rstatus = corrected_new_status_matrix[dd_idx]
)
dbh_dead_dt[, `:=`(
  raw_status = mapply(function(i, j) ViewFullTable_split[[j]]$Status[i], row, census),
  codes = mapply(function(i, j) ViewFullTable_split[[j]]$ListOfTSM[i], row, census),
  original_status = original_status[row],
  corrected_history = corrected_new_status[row]
)]
field_review[["dbh_on_dead_record"]] <- review_rows(
  "dbh_on_dead_record", dbh_dead_dt$StemID, dbh_dead_dt$census,
  paste0("codes ", fcoalesce(dbh_dead_dt$codes, ""), " | raw codes ", dbh_dead_dt$original_status)
)
dead_dbh_cells <- DT_Status[!is.na(DBH) & Status %in% c("dead", "stem dead"), .(row = match(StemID, unique_StemID$StemID), census)]
dead_dbh_R <- corrected_new_status_matrix[cbind(dead_dbh_cells$row, dead_dbh_cells$census)]
cat(sprintf(
  "🔎 DBH on dead / stem dead records: %d, all kept | on G/D (stem never alive later): %d | on A (false death, stem alive later): %d\n",
  nrow(dead_dbh_cells), sum(dead_dbh_R != "A"), sum(dead_dbh_R == "A")
))
bio_check(
  nrow(dbh_dead_dt) == sum(dead_dbh_R != "A") && all(dbh_dead_dt$raw_status %in% c("dead", "stem dead")),
  "A DBH sits on a G/D cell only on a dead / stem dead record of a stem never alive later",
  examples = dbh_dead_dt[!raw_status %in% c("dead", "stem dead"), StemID],
  n_bad = dbh_dead_dt[!raw_status %in% c("dead", "stem dead"), .N]
)
bio_check(
  !any(corrected_new_status_matrix == "P" & !is.na(DBHs)),
  "No P cell carries a DBH (P is only before a stem's first record)",
  n_bad = sum(corrected_new_status_matrix == "P" & !is.na(DBHs))
)
rm(dd_idx, dead_dbh_cells, dead_dbh_R)

# ---- (c) status × DBH × census summary (to the log) ---------------------
# An A without DBH is a missed measurement (stage 4 interpolates its DBH) or
# a stump (DFstatus "broken below": stage 4 keeps it out of the stock).
n_cens <- ncol(corrected_new_status_matrix)
cells_dt <- data.table(
  status  = as.vector(corrected_new_status_matrix),
  has_DBH = !is.na(as.vector(DBHs)), # exported dbh (raw)
  census  = rep(seq_len(n_cens), each = nrow(corrected_new_status_matrix))
)
summary_long <- cells_dt[, .(n_cells = .N), by = .(census, status, has_DBH)][order(census, status, has_DBH)]
summary_overall <- cells_dt[, .(n_cells = .N), by = .(status, has_DBH)][order(status, has_DBH)]
summary_overall[, pct := round(100 * n_cells / sum(n_cells), 3)]
summary_wide <- dcast(summary_long, status + has_DBH ~ census, value.var = "n_cells", fill = 0L)
n_P_with_DBH <- summary_overall[status == "P" & has_DBH == TRUE, sum(n_cells)]
n_A_no_DBH <- summary_overall[status == "A" & has_DBH == FALSE, sum(n_cells)]
n_other <- summary_overall[!(status %in% c("A", "D", "G", "P")), sum(n_cells)]
n_stump <- sum(is_status("broken below") & is.na(DBHs) & corrected_new_status_matrix == "A")
stump_txt <- sprintf(
  "A without DBH: %d = %d stumps (DFstatus 'broken below'; stage 4 keeps them out of the stock) + %d missed measurements (stage 4 interpolates their DBH)",
  n_A_no_DBH, n_stump, n_A_no_DBH - n_stump
)
print_to_log(c(
  sprintf(
    "Status x DBH of the exported tables (%d stems x %d censuses = %d cells), all censuses:",
    nrow(corrected_new_status_matrix), n_cens, nrow(cells_dt)
  ),
  log_table(summary_overall),
  paste0("  ", stump_txt),
  "  G / D with DBH: recorded on a dead record and kept (reason dbh_on_dead_record); P with DBH must be 0",
  "Per census:",
  log_table(summary_wide)
), log_file, new_message = TRUE)
cat("\nStatus × DBH overall summary:\n")
print(summary_overall)
cat(" ", stump_txt, "\n")
bio_check(
  sum(n_P_with_DBH) == 0L,
  "No P cell carries a DBH in the exported dbh",
  n_bad = sum(n_P_with_DBH)
)
bio_check(
  sum(n_other) == 0L,
  "Every exported status is one of A, D, G, P",
  n_bad = sum(n_other)
)

rm(
  cells_dt, summary_long, summary_overall, summary_wide,
  n_P_with_DBH, n_A_no_DBH, n_other, n_stump, stump_txt
)

# ========================================================================
# SECTION 14: EXPORT CENSUS TABLES
# ========================================================================
# Write Rstatus into each census table, assign one location per tree, fill
# the date column, rename and subset each table to the ForestGEO R-table
# format, sort the rows by treeID and stemID (as numbers), then save each
# census as both a .Rdata object and a .csv file. dbh, DFstatus and ExactDate
# are exported exactly as recorded; stemID is the reconstructed stem number
# within its tree. At the end, write the two review lists (Sections 12-14)
# and their counts to the log.

# Rows of every census table are in unique_StemID order (checked in
# Section 3), so the status and DBH matrices align column by column.
for (census in seq_along(ViewFullTable_split)) {
  bio_check(
    identical(ViewFullTable_split[[census]]$StemID, unique_StemID$StemID),
    sprintf("Census %d rows are in master stem order before Rstatus is written", census)
  )
  bio_check(
    identical(as.numeric(ViewFullTable_split[[census]]$DBH), as.numeric(DBHs[, census])),
    sprintf("Census %d: exported dbh is the raw DBH, unchanged", census)
  )
  ViewFullTable_split[[census]]$new_status <- corrected_new_status_matrix[, census]
}

# ========================================================================
# LOCATION: one (PX, PY) pair per tree for every stem and census
# ========================================================================
# A tree does not move. Its raw coordinates can differ between censuses
# (re-measurement, rounding, entry errors), and PX / PY are NA in censuses
# where a stem has no row. Assuming an error is a one-off, the position
# that most censuses agree on is the correct one:
#   (a) each census casts ONE vote per tree: the modal (PX, PY) pair among
#       that tree's stems in that census (ties → smallest x, then y), so a
#       census with many stems does not outvote the others;
#   (b) the tree position is the modal pair across census votes
#       (ties → the most recent census);
#   (c) that pair is written to every stem of the tree in every census;
#   (d) QuadratName is derived from the pair (20 m quadrats named "XXYY"),
#       with the modal raw quadrat as fallback for trees without coordinates;
#   (e) trees whose raw pairs spread more than conflict_threshold_m, trees
#       with no coordinates, and raw quadrats that disagree with the pair
#       go to the field review list.
# PX and PY are always chosen TOGETHER, so the pair was actually recorded.
# ========================================================================
# Arguments: split_list (list of census tables in the same row order) and the
# names of its tree, x, y and quadrat columns; quadrat_size and the plot
# dimensions plot_x, plot_y (m); conflict_threshold_m.
# Returns list(split_list = the tables with the tree position written in,
# conflicts, missing_trees, quadrat_mismatch).
impute_tree_location <- function(split_list,
                                 tree_col = "TreeID",
                                 x_col = "PX",
                                 y_col = "PY",
                                 quadrat_col = "QuadratName",
                                 quadrat_size = 20,
                                 plot_x = 1000,
                                 plot_y = 500,
                                 conflict_threshold_m = 1) {
  # Pool every recorded pair with its census.
  pool <- rbindlist(lapply(seq_along(split_list), function(i) {
    dt <- split_list[[i]]
    keep <- !is.na(dt[[x_col]]) & !is.na(dt[[y_col]])
    data.table(
      tree = dt[[tree_col]][keep], census = i,
      x = dt[[x_col]][keep], y = dt[[y_col]][keep]
    )
  }))

  # (a) One vote per tree × census.
  census_vote <- pool[, .N, by = .(tree, census, x, y)]
  setorder(census_vote, tree, census, -N, x, y)
  census_vote <- census_vote[census_vote[, .I[1L], by = .(tree, census)]$V1]

  # (b) Tree position: most census votes; ties → most recent census.
  tree_vote <- census_vote[, .(n_votes = .N, last_census = max(census)), by = .(tree, x, y)]
  setorder(tree_vote, tree, -n_votes, -last_census)
  ref <- tree_vote[tree_vote[, .I[1L], by = tree]$V1, .(tree, x_ref = x, y_ref = y, n_votes)]
  n_censuses_xy <- census_vote[, .(n_censuses_with_xy = .N), by = tree]
  ref <- n_censuses_xy[ref, on = "tree"]

  # (d) Quadrat derived from the pair; modal raw quadrat as fallback.
  qx <- pmin(floor(ref$x_ref / quadrat_size), plot_x / quadrat_size - 1)
  qy <- pmin(floor(ref$y_ref / quadrat_size), plot_y / quadrat_size - 1)
  ref[, quadrat_ref := sprintf("%02d%02d", as.integer(qx), as.integer(qy))]

  quad_all <- rbindlist(lapply(split_list, function(dt) {
    keep <- !is.na(dt[[quadrat_col]])
    data.table(tree = dt[[tree_col]][keep], q = dt[[quadrat_col]][keep])
  }))
  quad_pool <- unique(quad_all)
  quad_counts <- quad_all[, .N, by = .(tree, q)]
  setorder(quad_counts, tree, -N, q)
  quad_mode <- quad_counts[quad_counts[, .I[1L], by = tree]$V1, .(tree, q)]

  # (e) Review tables.
  spread <- unique(pool[, .(tree, x, y)])[, .(
    n_distinct_pairs = .N,
    spread_m = sqrt(diff(range(x))^2 + diff(range(y))^2)
  ), by = tree]
  conflicts <- ref[spread[spread_m > conflict_threshold_m], on = "tree"]
  setorder(conflicts, -spread_m)
  all_trees <- unique(unlist(lapply(split_list, function(dt) dt[[tree_col]])))
  missing_trees <- setdiff(all_trees, ref$tree)
  quadrat_mismatch <- quad_pool[ref[, .(tree, quadrat_ref)], on = "tree", nomatch = 0L][q != quadrat_ref]
  setnames(quadrat_mismatch, "q", "raw_quadrat")

  cat(sprintf(
    "  [impute_tree_location] trees with a position: %d | no coordinates: %d | raw spread > %g m: %d | raw quadrat != derived: %d trees\n",
    nrow(ref), length(missing_trees), conflict_threshold_m, nrow(conflicts), uniqueN(quadrat_mismatch$tree)
  ))

  # (c) Write the tree position (and quadrat) to every stem and census.
  split_list <- lapply(split_list, function(dt) {
    dt <- copy(dt)
    tr <- dt[[tree_col]]
    i_ref <- match(tr, ref$tree)
    q_new <- ref$quadrat_ref[i_ref]
    no_xy <- is.na(q_new)
    q_new[no_xy] <- quad_mode$q[match(tr[no_xy], quad_mode$tree)]
    set(dt, j = x_col, value = ref$x_ref[i_ref])
    set(dt, j = y_col, value = ref$y_ref[i_ref])
    set(dt, j = quadrat_col, value = q_new)
    dt
  })

  setnames(conflicts, "tree", tree_col)
  setnames(quadrat_mismatch, "tree", tree_col)
  list(
    split_list = split_list,
    conflicts = conflicts,
    missing_trees = missing_trees,
    quadrat_mismatch = quadrat_mismatch
  )
}

cat("🗺️  Assigning one location per tree (modal PX/PY pair across censuses)...\n")
location_fix <- impute_tree_location(ViewFullTable_split)
ViewFullTable_split <- location_fix$split_list
lc <- location_fix$conflicts
field_review[["location_spread_over_1m"]] <- review_rows(
  "location_spread_over_1m",
  tree = lc$TreeID,
  details = sprintf(
    "%d different (gx, gy) pairs, %.1f m apart; position used: %.1f, %.1f (the pair most censuses agree on)",
    lc$n_distinct_pairs, lc$spread_m, lc$x_ref, lc$y_ref
  )
)
field_review[["no_coordinates"]] <- review_rows(
  "no_coordinates",
  tree = location_fix$missing_trees,
  details = "no (gx, gy) in any census; quadrat as recorded"
)
qm <- location_fix$quadrat_mismatch
field_review[["quadrat_differs_from_coordinates"]] <- review_rows(
  "quadrat_differs_from_coordinates",
  tree = qm$TreeID,
  details = sprintf("recorded quadrat %s; quadrat of the coordinates %s (used)", qm$raw_quadrat, qm$quadrat_ref)
)
rm(lc, qm)

# A tree has one location: one (PX, PY, QuadratName) across all its stems
# and censuses (trees without coordinates keep NA PX/PY and one quadrat).
loc_check <- unique(rbindlist(lapply(ViewFullTable_split, function(dt) {
  dt[, .(TreeID, PX, PY, QuadratName)]
})))[, .N, by = TreeID][N > 1L]
bio_check(
  nrow(loc_check) == 0L,
  "Every tree has exactly one (PX, PY, QuadratName) across all stems and censuses",
  examples = loc_check$TreeID,
  n_bad = nrow(loc_check)
)
rm(location_fix, loc_check)
cat("✓ Location assignment complete.\n\n")

# ========================================================================
# DATES: ExactDate exported raw; `date` filled on every row
# ========================================================================
# ExactDate is the field date of a record and is exported exactly as
# recorded (NA where there is no record), like dbh and DFstatus.
# `date` is the ForestGEO R-table date: days since 1960-01-01, the convention
# of the ForestGEO database (its Date field, equal to ExactDate on every
# record) and of the CTFS / fgeo R functions, which compute census intervals
# from it (a recruit's interval starts at its P row, a death's ends at its
# first G/D row, so every row needs one). It is the recorded date where there
# is one; otherwise the MODE of the recorded dates of, in order,
#   (1) the same tree, same census   : stems of one tree are measured together
#   (2) the same quadrat, same census: crews census a quadrat within days
#   (3) the whole census             : last resort (e.g. trees without a quadrat)
# The mode (not the median or mean) is always a day on which the field crew
# was recording, and the tree's own date comes first because the records of
# one quadrat can span months (1982). Votes come only from recorded dates,
# and ties go to the earliest tied date so every run gives the same answer.
# ========================================================================

# Modal date of a vector; ties → earliest date.
date_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(as.Date(NA))
  }
  cnt <- data.table(d = x)[, .N, by = d]
  setorder(cnt, -N, d)
  cnt$d[1L]
}

# Modal date per group of a data.table; ties → earliest date.
date_mode_by <- function(dt, by_col, date_col) {
  cnt <- dt[, .N, by = c(by_col, date_col)]
  setorderv(cnt, c(by_col, "N", date_col), c(1L, -1L, 1L))
  cnt[cnt[, .I[1L], by = by_col]$V1, c(by_col, date_col), with = FALSE]
}

# Filled date of every row of every census (a list of Date vectors, one per
# census); ExactDate itself is not changed.
fill_tree_dates <- function(split_list,
                            tree_col = "TreeID",
                            date_col = "ExactDate",
                            quadrat_col = "QuadratName") {
  lapply(seq_along(split_list), function(i) {
    dt <- split_list[[i]]
    dates <- dt[[date_col]]
    n_recorded <- sum(!is.na(dates))
    # Only dates recorded in the data vote.
    recorded <- dt[!is.na(dates)]

    # (1) same tree, same census
    ref_tree <- date_mode_by(recorded, tree_col, date_col)
    need <- is.na(dates)
    dates[need] <- ref_tree[[date_col]][match(dt[[tree_col]][need], ref_tree[[tree_col]])]
    n_tree <- sum(need & !is.na(dates))

    # (2) same quadrat, same census
    ref_quad <- date_mode_by(recorded[!is.na(recorded[[quadrat_col]])], quadrat_col, date_col)
    need <- is.na(dates)
    dates[need] <- ref_quad[[date_col]][match(dt[[quadrat_col]][need], ref_quad[[quadrat_col]])]
    n_quad <- sum(need & !is.na(dates))

    # (3) whole census
    need <- is.na(dates)
    dates[need] <- date_mode(recorded[[date_col]])
    n_census <- sum(need & !is.na(dates))

    cat(sprintf(
      "  [fill_tree_dates] census %d: %d recorded | filled %d from tree, %d from quadrat, %d from census mode | %d without a date\n",
      i, n_recorded, n_tree, n_quad, n_census, sum(is.na(dates))
    ))
    dates
  })
}

# Recorded-date window of each census (before imputation), for the checks.
census_date_window <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  d <- ViewFullTable_split[[i]]$ExactDate
  data.table(census = i, lo = min(d, na.rm = TRUE), hi = max(d, na.rm = TRUE))
}))

recorded_dates <- lapply(ViewFullTable_split, function(dt) dt$ExactDate)
raw_date_field <- lapply(ViewFullTable_split, function(dt) dt$Date) # database Date (days since 1960-01-01)

cat("📅 Filling the date column (days since 1960-01-01; modal field date of tree, quadrat, census)...\n")
filled_dates <- fill_tree_dates(ViewFullTable_split)
date_origin <- as.Date("1960-01-01")
for (i in seq_along(ViewFullTable_split)) {
  set(ViewFullTable_split[[i]], j = "Date", value = as.numeric(filled_dates[[i]] - date_origin))
}

dates_long <- rbindlist(lapply(seq_along(ViewFullTable_split), function(i) {
  ViewFullTable_split[[i]][, .(
    StemID,
    census = i, ExactDate, recorded = recorded_dates[[i]],
    date = Date, raw_date = raw_date_field[[i]]
  )]
}))
changed <- dates_long[!((is.na(ExactDate) & is.na(recorded)) | (!is.na(ExactDate) & !is.na(recorded) & ExactDate == recorded))]
bio_check(
  nrow(changed) == 0L,
  "ExactDate is exported exactly as recorded (no date imputed, removed or replaced)",
  examples = changed$StemID,
  n_bad = nrow(changed)
)
na_dates <- dates_long[is.na(date)]
bio_check(
  nrow(na_dates) == 0L,
  "Every row has a date",
  examples = na_dates$StemID,
  n_bad = nrow(na_dates)
)
off_record <- dates_long[!is.na(recorded) & date != as.numeric(recorded - date_origin)]
bio_check(
  nrow(off_record) == 0L,
  "date equals ExactDate (days since 1960-01-01) wherever a date was recorded",
  examples = off_record$StemID,
  n_bad = nrow(off_record)
)
off_db <- dates_long[!is.na(raw_date) & date != raw_date]
bio_check(
  nrow(off_db) == 0L,
  "date equals the database's Date field wherever it was recorded",
  examples = off_db$StemID,
  n_bad = nrow(off_db)
)
cat(sprintf(
  "📅 ExactDate: %d recorded (exported as is) | date: %d recorded + %d filled = every row\n",
  dates_long[!is.na(recorded), .N], dates_long[!is.na(recorded), .N], dates_long[is.na(recorded), .N]
))
dates_long[, date_d := date_origin + date]
dates_long <- census_date_window[dates_long, on = "census"]
out_window <- dates_long[date_d < lo | date_d > hi]
bio_check(
  nrow(out_window) == 0L,
  "Every date lies inside its census's recorded date window",
  examples = out_window$StemID,
  n_bad = nrow(out_window)
)
setkey(dates_long, StemID, census)
dates_long[, prev_date := shift(date), by = StemID]
non_monotonic <- dates_long[!is.na(prev_date) & date <= prev_date]
bio_check(
  nrow(non_monotonic) == 0L,
  "Dates strictly increase from census to census for every stem",
  examples = non_monotonic$StemID,
  n_bad = nrow(non_monotonic)
)
rm(
  dates_long, changed, na_dates, off_record, off_db, out_window, non_monotonic, census_date_window,
  recorded_dates, raw_date_field, filled_dates
)

cat("✓ Dates complete.\n\n")

cat("💾 Exporting census tables to .Rdata files...\n")

## Format and export each census as ForestGEO R table ####
# Loop through each census and export in standardized format
check_data <- rbindlist(lapply(ViewFullTable_split, function(dt) dt[, ..ViewFullTable_columns_to_keep]))
setorder(check_data, TreeID, StemID, CensusID)
# check unique date per TreeID (dated rows)
check_dates <- unique(check_data[!is.na(ExactDate), .(TreeID, ExactDate, CensusID)])

# get nunique exactdate per treeid and census
check_dates[, c("n_dates") :=
  .(uniqueN(ExactDate)),
by = .(CensusID, TreeID)
]

rd <- unique(check_dates[n_dates > 1L, .(TreeID, CensusID, n_dates)])
field_review[["two_dates_in_census"]] <- review_rows(
  "two_dates_in_census",
  tree = rd$TreeID, census = census_index(rd$CensusID),
  details = sprintf("%d different field dates on the stems of this tree", rd$n_dates)
)
rm(rd)

# ---- Row order and column types ---------------------------------------------
# Rows are sorted by treeID, then by the stem number within the tree, both
# compared as numbers: as text, tree "100" would come before tree "20" and
# stem 10 before stem 2. Sorting changes no column: treeID stays text (as in
# the stage-2 table), and tag, StemTag and quadrat stay text because their
# leading zeros and signs are part of the code (tags "000006", "-05599",
# quadrat "0001").
bio_check(
  all(grepl("^[0-9]+$", unique_StemID$TreeID)),
  "Every TreeID is a whole number written in digits (it sorts as a number without loss)",
  examples = unique_StemID$TreeID[!grepl("^[0-9]+$", unique_StemID$TreeID)]
)
bio_check(
  all(grepl("_[0-9]+$", unique_StemID$StemID)),
  "Every internal stem key ends in its stem number (TreeID_<number>)",
  examples = unique_StemID$StemID[!grepl("_[0-9]+$", unique_StemID$StemID)]
)
export_types <- c(
  CensusID = "numeric", treeID = "character", stemID = "integer", tag = "character", StemTag = "character",
  sp = "character", quadrat = "character", gx = "numeric", gy = "numeric", dbh = "numeric", hom = "numeric",
  ExactDate = "Date", date = "numeric", codes = "character", DFstatus = "character", Rstatus = "character",
  StemPaths = "integer"
)
bio_check(
  identical(names(export_types), new_names_columns_to_keep),
  "The export types list every exported column, in export order"
)

export_stem_order <- vector("list", length(ViewFullTable_split))
for (census in seq_along(ViewFullTable_split)) {
  cat(sprintf("  Processing census %d...\n", census))
  # Extract current census data and ensure it's a data.table
  X <- as.data.table(ViewFullTable_split[[census]])
  # Sort by treeID and stem number as numbers (see above)
  X[, `:=`(.tree_num = as.numeric(TreeID), .stem_num = as.numeric(sub(".*_", "", StemID)))]
  setorder(X, .tree_num, .stem_num)
  X[, c(".tree_num", ".stem_num") := NULL]
  # Validate that every required column is present BEFORE subsetting,
  # so a missing column produces a clear error rather than a cryptic
  # data.table failure deep inside the export loop.
  missing_cols <- setdiff(ViewFullTable_columns_to_keep, names(X))
  bio_check(
    length(missing_cols) == 0L,
    sprintf("Census %d has every column required for export", census),
    examples = missing_cols
  )
  # Select only the columns needed for ForestGEO R format
  # ..ViewFullTable_columns_to_keep references the variable defined in Section 1
  X <- X[, ..ViewFullTable_columns_to_keep]
  # Rename columns from database format to ForestGEO R format
  # Example: TreeID → treeID, Mnemonic → sp, PX → gx, etc.
  setnames(X, old = ViewFullTable_columns_to_keep, new = new_names_columns_to_keep)
  # DFstatus is the raw field status of the dataset (ForestGEO column name),
  # exported exactly as recorded, NA where the dataset has no record. It is
  # never modified (fixed rule); Rstatus holds the corrected status.
  raw_status_census <- ViewFullTable_split[[census]]$Status[match(X$stemID, ViewFullTable_split[[census]]$StemID)]
  bio_check(
    identical(as.character(X$DFstatus), as.character(raw_status_census)),
    sprintf("Census %d: DFstatus is the raw Status, unchanged", census)
  )
  bio_check(
    all(X$Rstatus %in% status_codes),
    sprintf("Census %d: Rstatus only takes values in {A, D, G, P}", census),
    examples = unique(X$Rstatus[!X$Rstatus %in% status_codes])
  )
  export_stem_order[[census]] <- X$stemID
  # stemID: the reconstructed stem number within its tree (1, 2, ...). The
  # internal key is "TreeID_ReconstructedStemID"; its second part is
  # exported, so treeID + stemID identify a stem. (Rows without a DP
  # identity were dropped in Section 3, so no exported stemID is NA.)
  rs_part <- substring(X$stemID, nchar(X$treeID) + 2L)
  bad_key <- paste(X$treeID, rs_part, sep = "_") != X$stemID | !(rs_part == "NA" | grepl("^[0-9]+$", rs_part))
  bio_check(
    !any(bad_key),
    sprintf("Census %d: every internal stem key is TreeID_<reconstructed stem number>", census),
    examples = X$stemID[bad_key],
    n_bad = sum(bad_key)
  )
  X[, stemID := as.integer(fifelse(rs_part == "NA", NA_character_, rs_part))]
  bio_check(
    !anyDuplicated(X[, .(treeID, stemID)]),
    sprintf("Census %d: treeID + stemID identify exactly one row", census),
    n_bad = sum(duplicated(X[, .(treeID, stemID)]))
  )
  rm(rs_part, bad_key)
  bio_check(
    identical(order(as.numeric(X$treeID), X$stemID), seq_len(nrow(X))),
    sprintf("Census %d: rows are sorted by treeID, then stemID (as numbers)", census)
  )
  col_types <- vapply(X, function(v) class(v)[1L], "")
  bio_check(
    identical(names(col_types), names(export_types)) && all(col_types == export_types),
    sprintf("Census %d: every column has its type (IDs and codes text where they carry leading zeros, measurements numeric)", census),
    examples = names(col_types)[col_types != export_types[names(col_types)]]
  )
  rm(col_types)
  # Write the CSV, then convert to data.frame: many ForestGEO R functions
  # expect a data.frame, not a data.table
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

# ForestGEO census functions compare census tables row by row, so every
# table must list the same stems in the same order.
bio_check(
  all(vapply(export_stem_order, identical, logical(1), export_stem_order[[1]])),
  "All exported census tables list the same stems in the same order"
)
rm(export_stem_order, export_types)

cat(sprintf("\n✓✓ All %d census tables exported successfully\n\n", length(ViewFullTable_split)))

# ========================================================================
# REVIEW LISTS (built in Sections 12-14)
# ========================================================================
# One row per case; `reason` says what to look at (README of this folder).
# Nothing in them changes the R tables.
empty_review <- review_rows("none", tree = "none")[0L]
field_review_dt <- rbindlist(c(list(empty_review), field_review), use.names = TRUE)
recon_review_dt <- rbindlist(c(list(empty_review), recon_review), use.names = TRUE)
setorder(field_review_dt, reason, treeID, stemID, census)
setorder(recon_review_dt, reason, treeID, stemID, census)
fwrite(field_review_dt, file.path(CHECK_folder, paste0("field_records_to_review_", site, ".csv")))
fwrite(recon_review_dt, file.path(CHECK_folder, paste0("reconstruction_to_review_", site, ".csv")))
review_counts <- function(x) {
  if (nrow(x) == 0L) {
    return("  (no case)")
  }
  x[, .(cases = .N, trees = uniqueN(treeID)), keyby = reason][, sprintf("  %-34s %8d cases in %7d trees", reason, cases, trees)]
}
n_checks <- sum(startsWith(readLines(log_file, encoding = "UTF-8"), "✓"))
review_txt <- c(
  sprintf("field_records_to_review_%s.csv (raw records that need a look, for the data manager):", site),
  review_counts(field_review_dt),
  sprintf("reconstruction_to_review_%s.csv (stems whose reconstruction needs a look):", site),
  review_counts(recon_review_dt),
  "",
  sprintf("All %d checks passed; %d R tables exported.", n_checks, length(ViewFullTable_split))
)
print_to_log(review_txt, log_file, new_message = TRUE)
cat(paste(review_txt, collapse = "\n"), "\n")
rm(empty_review, review_counts, review_txt, n_checks)
