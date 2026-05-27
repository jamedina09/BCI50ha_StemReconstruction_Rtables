# ========================================================================
# SCRIPT: compare_tables1_and_tables2.R
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

# ---- 1. Configuration ----
# User-editable variables and file locations. Ensure the paths below point to
# the expected tab-delimited `ViewFullTable` files for the site.

site <- "bci" # Site code for BCI

# Input folders for the two datasets to compare
INPUT_folder_1 <- here::here("BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses")

# Output and diagnostics folders
OUTPUT_folder <- here::here("BCI_stem_reconstruction", "DATA", "PROCESSED")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
}

CHECK_folder <- here::here("BCI_stem_reconstruction", "DATA", "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
}

# End of configuration (adjust above values as needed for local runs)

# ---- 2. Load and verify input data ----
# Read the two versions of ViewFullTable. We expect tab-delimited tables with
# columns such as `ExactDate`, `Tag`, `StemTag`, `TreeID`, `StemID`, `DBH`.

file1 <- file.path(INPUT_folder_1, paste0("ViewFullTable_", site, ".csv"))
stopifnot(file.exists(file1))

ViewFullTable <- fread(
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

# ! FIXME:
# Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
spp_new[especie %in% "rohrii"]

# this needs to be replaced with the correct code
unique(sp_bci_raw_input[SpeciesName == "officinalis", .(Mnemonic, Family, Genus, SpeciesName)])

# Correct
sp_bci_raw_input[Mnemonic == "pterro"]
# Replace with pterro
sp_bci_raw_input[Mnemonic == "pterof"]

# FIXME: pterof code in BCI data needs to be replaced by pterro.
# !

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

# NOTE: 13 INDIVIDUALS WITH DBH < 10 CM — LIKELY TAPERING ERRORS
# NOTE: no individual with zero DBH
table(ViewFullTable[DBH < 10]$DBH)
range(ViewFullTable$DBH, na.rm = TRUE)

# sort print(inspectdf::inspect_na(ViewFullTable), n = 50) by tag and censusID
## NOTE: Tag and CensusID have no NAs, so we can use them to create a complete grid

tag_census_unique <- unique(copy(ViewFullTable[, .(Tag, CensusID)]))

# Create complete grid of all Tag-CensusID combinations
tag_ranges <- tag_census_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]

## Get all combinations tags and census IDs
complete_grid <- tag_ranges[, .(CensusID = seq.int(min_c, max_c)), by = Tag]

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

unique(names(ViewFullTable) == names(missing_combinations))

# Add the missing rows to the main dataset
ViewFullTable_no_missing_tags_census <- rbindlist(
  list(ViewFullTable[, missing_tag_census := FALSE], missing_combinations[, missing_tag_census := TRUE]),
  fill = TRUE,
  use.names = TRUE
)

(nrow(ViewFullTable) + nrow(missing_combinations)) == nrow(ViewFullTable_no_missing_tags_census)

# Sort by Tag and CensusID for cleaner viewing
setorder(ViewFullTable_no_missing_tags_census, Tag, CensusID)

## CHECK AGAIN IF MISSING TAG CENSUS COMBINATIONS

missing_combinations_final <- complete_grid[!ViewFullTable_no_missing_tags_census,
  on = .(Tag, CensusID)
]
# good

# ========================================================================
# SCRIPT: 2_check_revised_viewfulltable.R
# PURPOSE: Identify and correct DBH / HOM measurement issues before stem
#          identification. This script produces diagnostic files and a
#          taper-corrected ViewFullTable for downstream steps.
#
# QUICK TOC:
#   0. Setup & packages
#   1. Load input data
#   2. Measurement error checks (HOM & DBH overview)
#   3. Fix 1: select max HOM per census
#   4. Recompute QC using corrected HOM
#   5. Check missing Tag × CensusID combinations
#   6. Detect potential DBH measurement errors
#     6.1 DBH outlier detection (log-difference method)
#     6.2 Stem-level summary & scoring
#   7. Export: Suspicious measurements for manual QC
#   8. Prepare & taper correction
#     8.1 Pre-taper checks and DBH candidate preparation
#     8.2 Apply taper correction (DBH → DBHC)
#   9. Finalize complete table (add missing rows)
# ========================================================================

# # ---- 0. Setup & packages ----
# rm(list = ls())

# # Core packages
# library(data.table) # fast table manipulation
library(ggplot2) # diagnostics/plots
# library(here) # project-relative paths
# library(patchwork)

# # Set data.table threads to use available physical cores (tune if needed)
# getDTthreads() # current threads
# setDTthreads(parallel::detectCores(logical = FALSE))
# getDTthreads() # show updated threads
# setDTthreads(1)
# # Notes:
# # - Keep operations in data.table style for memory and speed.
# # - Where possible, caching (saveRDS/readRDS) is used to avoid expensive recomputations.

# ---- 1. Load input data ----
# Read previously validated ViewFullTable and Condit-derived datasets used for
# cross-checks and reference comparisons.
ViewFullTable <- ViewFullTable_no_missing_tags_census # as.data.table(readRDS(here::here("DATA", "PROCESSED", "1_ViewFullTable_no_missing_tags_census.rds")))
ViewFullTable <- ViewFullTable[, missing_tag_census := NULL] # remove auxiliary column
## add rowid for tracking
ViewFullTable[, RowID := .I]

# Check wether all tags have the whole observations between minimum and maximum census
ViewFullTable[, .(CensusID = list(sort(unique(CensusID)))), by = Tag][
  , expected := .(list(seq(min(unlist(CensusID)), max(unlist(CensusID))))),
  by = Tag
][
  , complete := identical(CensusID[[1]], expected[[1]])
][
  , .N,
  by = complete
]
# NOTE: Yes, they do. All tags have complete observations between their minimum and maximum census.

# ---- 2. Measurement error checks ----
# Run several QC steps to identify likely measurement errors (DBH/HOM).
# Key steps: choose one measurement per Tag x Stem x Census (highest HOM),
# compute temporal log-differences to flag outliers, propose candidate values,
# and run taper correction.

# ---- Stem identification: notes & assumptions ----
# Variables and assumptions used in later stem matching and QC steps:
#   - Tag: plot + visible tag — primary grouping identifier across censuses
#   - CensusID: integer temporal index (1,2,3…) used to order measurements
#   - DBH: measured diameter at the recorded HOM (may be NA if not measured)
#   - HOM: height of measurement (metres) — used to select one measurement per census and for taper correction in buttressed trees
#   - Species / taxonomy: helps limit candidate matches when matching stems
# Historical notes: stem tags were assigned later; early censuses may lack StemTag and duplicate measurements may exist for buttressed trees.

# The number of tree IDs and Tags is the same
length(unique(ViewFullTable$Tag))
length(unique(ViewFullTable$TreeID))

# compare ntags per census with ntreeids per census
unique(ViewFullTable[, .N, by = .(CensusID, Tag)][, .N, by = CensusID] == ViewFullTable[, .N, by = .(CensusID, TreeID)][, .N, by = CensusID])

# Check whether Tag and TreeID have the same number of observations
inc_tag <- ViewFullTable[, .N, by = .(Tag)]
inc_treeid <- ViewFullTable[, .N, by = .(TreeID)]
# Yes, they do
all(inc_tag$N == inc_treeid$N)

# Note: Suzanne's main comments
# In the case of BCI, this is complicated by the fact that stem tags were not assigned to the
# multiple stems of a tree until around 2020. Following are some  observations:

# At BCI, the stems of a tree were not given tags in the first censuses. The stems of
# trees with multiple stems were tagged starting in 2010.

# Because there was no reliable and unequivocal way to match stems from a census to
# those in the previous ones for multiple-stemmed trees, new stemids were assigned
# to each of these multiple stems in each census  in the database.

# For the R stem tables, Rick Condit developed an algorithm to try to match stems
# across censuses when it was clear which stem was which throughout the censuses.

# When the measurement dieerences between stems were distinct enough throughout
# the censuses —like in tags 001112 and 003036— he could confidently match stems
# from one census to another. Then, the latest stemID was applied retroactively to
# earlier censuses for consistency.

# Check example tags where stem IDs were reassigned
inc <- c("001112", "003036")

ViewFullTable[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM)
][!is.na(DBH)]

## Function to plot stems
plot_stem <- function(data, tag) {
  p <- ggplot(
    data[Tag %in% tag],
    aes(
      x = factor(CensusID),
      y = DBH,
      color = factor(StemID),
      group = factor(StemID)
    )
  ) +
    geom_line() +
    geom_point() +
    labs(
      title = paste("DBH over time for Tag", tag),
      x = "CensusID",
      y = "DBH"
    ) +
    theme_minimal()
  if (length(tag) > 1) {
    p <- p + facet_wrap(~Tag,
      scales = "free_x",
      nrow = ceiling(length(tag) / 2)
    )
  }
  return(p)
}

plot_stem(ViewFullTable, inc)

# NOTE: Measurement-selection and taper-correction rules used downstream
#   - When a (Tag, StemID, Census) has multiple observations (e.g. tag 003036,
#     census 3), select the row with the highest HOM (smallest DBH due to taper).
#   - StemIDs distinguish stems within a census, enabling the HOM selection step.
#   - Taper correction is needed because stem matching across censuses requires
#     all measurements standardized at HOM = 1.3 m (DBH). It is applied to:
#       (i) HOM > 1.3 m with code "B" (buttressed),
#       (ii) HOM < 1.3 m (correct upward to 1.3 m),
#       (iii) HOM > 1.3 m without "B" (e.g. swelling/bifurcation at 1.3 m).
#   - Correction is used only for stem matching; original DBH values in the R
#     tables remain unchanged.
#   - Absence of any stem tag after 2010 implies a single-stemmed tree (at
#     least across the last three censuses).

# When the measurements between stems were not clearly distinguishable throughout
# In some censuses the original stemIDs were retained, as in the case of tag 151991.

inc <- c("151991")

plot_stem(ViewFullTable, inc)

ViewFullTable[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM, Mnemonic)
][!is.na(DBH)]

bci_dryad[Tag %in% inc][
  order(Tag, StemTag, TreeID, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, codes, sp)
][!is.na(DBH)]

# Note: please make sure to use one measurement for the B (buttressed) trees where
# more than one measurement per census may occur. The one where the hom is
# highest should be used for the R tables.

# ---- 3. Fix 1: select observation with highest HOM per group ----
# Problem: multiple DBH/HOM records per Tag x Stem x Census. Rule: choose the
# observation with the highest HOM (prefer higher measurement point for buttressed trees).
# Quick check of HOM distribution
table(ViewFullTable$HOM, useNA = "ifany")

# --
# Step 1: Identify groups with multiple observations (same Tag, TreeID, StemTag, StemID, CensusID)
# --
multi_obs_groups <- ViewFullTable[
  , .N,
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
][N > 1, .(Tag, TreeID, StemTag, StemID, CensusID)]

inc <- multi_obs_groups$Tag[1:6]

plot_stem(ViewFullTable, inc)

# --
# Step 2: For multi-row groups, compute quality control (QC) flags and HOM_max
# --
multi_qc <- ViewFullTable[multi_obs_groups, on = .(Tag, TreeID, StemTag, StemID, CensusID)][
  , c(
    "QC_has_any_zero_HOM_per_group", # Is there any HOM==0 in group?
    "QC_has_all_zero_HOM_per_group", # Are all HOM==0 in group?
    "QC_has_any_NA_HOM_per_group", # Is there any HOM==NA in group?
    "HOM_max", # Maximum HOM in group
    "QC_has_valid_HOM", # Is there any valid HOM?
    "QC_is_max_HOM", # Is this row the max HOM?
    "QC_is_tied_max_HOM", # Is max HOM tied (multiple rows)?
    "QC_all_HOM_NA" # Are all HOM NA?
  ) := {
    HOM_clean <- HOM
    HOM_clean[is.na(HOM_clean)] <- -Inf # Treat NA as -Inf for max
    m <- max(HOM_clean) # Find max HOM
    has_any_zero <- any(HOM == 0, na.rm = TRUE)
    has_all_zero <- all(HOM == 0, na.rm = TRUE)
    has_any_NA <- any(is.na(HOM))
    has_valid_HOM <- m > -Inf
    all_HOM_NA <- m == -Inf
    is_max <- HOM == m & has_valid_HOM
    n_max <- sum(is_max, na.rm = TRUE)
    is_tied_max <- is_max & n_max > 1
    list(
      has_any_zero,
      has_all_zero,
      has_any_NA,
      ifelse(has_valid_HOM, m, NA_real_), # HOM_max
      has_valid_HOM,
      is_max,
      is_tied_max,
      all_HOM_NA
    )
  },
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
]

print(inspectdf::inspect_na(multi_qc), n = 70)

# --
# Step 3: For single-row groups, assign default QC flags
# --
single_obs <- ViewFullTable[!multi_obs_groups, on = .(Tag, TreeID, StemTag, StemID, CensusID)]
single_obs[
  , c(
    "QC_has_any_zero_HOM_per_group",
    "QC_has_all_zero_HOM_per_group",
    "QC_has_any_NA_HOM_per_group",
    "HOM_max",
    "QC_has_valid_HOM",
    "QC_is_max_HOM",
    "QC_is_tied_max_HOM",
    "QC_all_HOM_NA"
  ) := list(
    HOM == 0, # any zero?
    HOM == 0, # all zero? (single row)
    is.na(HOM), # any NA?
    HOM, # HOM_max = HOM
    !is.na(HOM), # has valid HOM
    TRUE, # single row is max by definition
    FALSE, # cannot be tied
    is.na(HOM) # all NA?
  )
]

# --
# Step 4: Combine multi- and single-row groups for HOM assessment
# --
ViewFullTable_hom_assessment <- rbindlist(list(multi_qc, single_obs), use.names = TRUE)
setorder(ViewFullTable_hom_assessment, RowID)

nrow(ViewFullTable) == nrow(ViewFullTable_hom_assessment) # should be TRUE

# Create HOM_link column to store corrected HOM values (start with original HOM)
ViewFullTable_hom_assessment[, HOM_link := HOM]

# If all HOM in group are zero, set HOM_link to 1.3 (default correction)
ViewFullTable_hom_assessment[
  QC_has_all_zero_HOM_per_group == TRUE,
  HOM_link := 1.3
]

# Groups with `QC_has_any_zero_HOM_per_group == TRUE` are the same as those with `QC_has_all_zero_HOM_per_group == TRUE`.
ViewFullTable_hom_assessment[, .N, by = .(QC_has_any_zero_HOM_per_group, QC_has_all_zero_HOM_per_group)]
ViewFullTable_hom_assessment[, QC_has_any_zero_HOM_per_group := NULL]
ViewFullTable_hom_assessment[, QC_has_all_zero_HOM_per_group := NULL]

# If all HOM in group are NA, set HOM_link to 1.3
ViewFullTable_hom_assessment[
  QC_all_HOM_NA == TRUE,
  HOM_link := 1.3
]

ViewFullTable_hom_assessment[, QC_all_HOM_NA := NULL]

# Inspect HOM values for groups with `QC_has_any_NA_HOM_per_group == TRUE`
table(ViewFullTable_hom_assessment[
  QC_has_any_NA_HOM_per_group == TRUE
]$HOM, useNA = "ifany")

# All those groups have HOM values of 1.3; set any NA `HOM_link` to 1.3
table(ViewFullTable_hom_assessment[
  QC_has_any_NA_HOM_per_group == TRUE
]$HOM_link, useNA = "ifany")

ViewFullTable_hom_assessment[
  QC_has_any_NA_HOM_per_group == TRUE & is.na(HOM_link),
  HOM_link := 1.3
]
ViewFullTable_hom_assessment[, QC_has_any_NA_HOM_per_group := NULL]

# Remove QC_has_valid_HOM column after fixing invalid HOM
table(ViewFullTable_hom_assessment[QC_has_valid_HOM == FALSE]$HOM_link, useNA = "ifany")
# all non-valid HOM's (i.e., NAs, ZEROS) are fixed
ViewFullTable_hom_assessment[, QC_has_valid_HOM := NULL]

# For groups with tied max HOM, select the row with largest DBH
ViewFullTable_hom_assessment[QC_is_tied_max_HOM == TRUE,
  QC_is_max_DBH_among_tied_HOM := DBH == max(DBH, na.rm = TRUE),
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
]
table(ViewFullTable_hom_assessment$QC_is_max_DBH_among_tied_HOM, useNA = "ifany")

ViewFullTable_hom_assessment <- ViewFullTable_hom_assessment[QC_is_max_DBH_among_tied_HOM == TRUE | is.na(QC_is_max_DBH_among_tied_HOM)]
ViewFullTable_hom_assessment[, QC_is_tied_max_HOM := NULL]
ViewFullTable_hom_assessment[, QC_is_max_DBH_among_tied_HOM := NULL]

# Check for any remaining NA or infinite HOM_link values
table(ViewFullTable_hom_assessment$HOM_link, useNA = "ifany")
table(is.na(ViewFullTable_hom_assessment$HOM_link))
table(is.infinite(ViewFullTable_hom_assessment$HOM_link))

# Check number of rows again
# Those groups with tied HOM were resolved by selecting the row with the largest DBH
# so, there were
nrow(ViewFullTable) - nrow(ViewFullTable_hom_assessment) # number removed

# Get rows that differ between `ViewFullTable_hom_assessment` and `ViewFullTable`
diff_rows <- fsetdiff(
  ViewFullTable,
  ViewFullTable_hom_assessment[, names(ViewFullTable), with = FALSE]
)

ViewFullTable[Tag %in% diff_rows$Tag][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM)
][order(Tag, StemID, CensusID)]
# Six rows were removed due to tied HOM and selection of the largest DBH

# ---- 4. Recompute QC using corrected HOM ----
# Re-run the HOM QC using `HOM_link` (corrected HOM values) to ensure a single
# representative observation per (Tag, Stem, CensusID) remains.

# Identify multi-observation groups again (after HOM correction)
multi_obs_groups <- ViewFullTable_hom_assessment[
  , .N,
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
][N > 1, .(Tag, TreeID, StemTag, StemID, CensusID)]

# For multi-row groups, compute QC flags and HOM_max using HOM_link
multi_qc <- ViewFullTable_hom_assessment[multi_obs_groups, on = .(Tag, TreeID, StemTag, StemID, CensusID)][
  , c(
    "HOM_max",
    "QC_has_valid_HOM",
    "QC_is_max_HOM",
    "QC_is_tied_max_HOM"
  ) := {
    HOM_clean <- HOM_link
    HOM_clean[is.na(HOM_clean)] <- -Inf
    m <- max(HOM_clean)
    has_valid_HOM <- m > -Inf
    is_max <- HOM_link == m & has_valid_HOM
    n_max <- sum(is_max, na.rm = TRUE)
    is_tied_max <- is_max & n_max > 1
    list(
      ifelse(has_valid_HOM, m, NA_real_), # HOM_max
      has_valid_HOM,
      is_max,
      is_tied_max
    )
  },
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
]

# For single-row groups, assign default QC flags
single_obs <- ViewFullTable_hom_assessment[!multi_obs_groups, on = .(Tag, TreeID, StemTag, StemID, CensusID)]
single_obs[
  , c(
    "HOM_max",
    "QC_has_valid_HOM",
    "QC_is_max_HOM",
    "QC_is_tied_max_HOM"
  ) := list(
    HOM_link, # HOM_max = HOM_link
    !is.na(HOM_link), # has valid HOM
    TRUE, # single row is max by definition
    FALSE # cannot be tied
  )
]

# Combine multi- and single-row groups for final corrected HOM table
ViewFullTable_hom_corrected <- rbindlist(list(multi_qc, single_obs), use.names = TRUE)

nrow(ViewFullTable_hom_corrected) == nrow(ViewFullTable_hom_assessment) # should be TRUE

# Remove QC_has_valid_HOM column (all should be TRUE now)
table(ViewFullTable_hom_corrected$QC_has_valid_HOM, useNA = "ifany")
ViewFullTable_hom_corrected[, QC_has_valid_HOM := NULL]

# For any remaining tied max HOM, select row with largest DBH
ViewFullTable_hom_corrected[QC_is_tied_max_HOM == TRUE,
  QC_is_max_DBH_among_tied_HOM := DBH == max(DBH, na.rm = TRUE),
  by = .(Tag, TreeID, StemTag, StemID, CensusID)
]
ViewFullTable_hom_corrected <- ViewFullTable_hom_corrected[QC_is_max_DBH_among_tied_HOM == TRUE | is.na(QC_is_max_DBH_among_tied_HOM)]
ViewFullTable_hom_corrected[, QC_is_tied_max_HOM := NULL]
ViewFullTable_hom_corrected[, QC_is_max_DBH_among_tied_HOM := NULL]

# Final check: all rows should have HOM_max and QC_is_max_HOM
sort(table(ViewFullTable_hom_corrected$HOM_max, useNA = "ifany"))
table(ViewFullTable_hom_corrected$QC_is_max_HOM, useNA = "ifany")

# HOM_link
# Inspect examples where `QC_is_max_HOM` is FALSE
ViewFullTable_hom_corrected[Tag %in% ViewFullTable_hom_corrected[QC_is_max_HOM == FALSE]$Tag[1]]
ViewFullTable_hom_corrected[Tag %in% ViewFullTable_hom_corrected[QC_is_max_HOM == FALSE]$Tag[2]]
ViewFullTable_hom_corrected[Tag %in% ViewFullTable_hom_corrected[QC_is_max_HOM == FALSE]$Tag[3]]

ViewFullTable_hom_corrected_clean <- ViewFullTable_hom_corrected[QC_is_max_HOM == TRUE]
ViewFullTable_hom_corrected_clean[, HOM_link := NULL]
table(ViewFullTable_hom_corrected_clean$QC_is_max_HOM, useNA = "ifany")
ViewFullTable_hom_corrected_clean[, QC_is_max_HOM := NULL]
ViewFullTable_hom_corrected_clean[, HOM_max := NULL]

# set order and check number of rows again
setorder(ViewFullTable_hom_corrected_clean, RowID)

# count observations per group after correction
table(ViewFullTable_hom_corrected_clean[, .N, by = .(Tag, TreeID, StemTag, StemID, CensusID)]$N)

# Print memory size in MB (optional)
# print(object.size(ViewFullTable_hom_corrected_clean), units = "MB")

ViewFullTable_hom_corrected_clean
ViewFullTable_hom_assessment

# Remove intermediate data.tables to free memory
# rm(ViewFullTable)
rm(ViewFullTable_hom_assessment)
rm(multi_qc)
rm(single_obs)
rm(ViewFullTable_hom_corrected)
gc()

# Include rowid for tracking
ViewFullTable_hom_corrected_clean[, RowIDN1 := .I]

# ---- 5. Check missing Tag × CensusID combinations ----
# Compute expected vs actual counts per Tag and identify gaps (tags with missing censuses).
xraw_unique <- unique(ViewFullTable_hom_corrected_clean[, .(Tag, CensusID)])
# Get the range per tag
tag_ranges <- xraw_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]
# Add expected count (how many censuses should exist)
tag_ranges[, expected_count := max_c - min_c + 1L]
# Get actual count per tag
actual_counts <- xraw_unique[, .(actual_count = .N), by = Tag]
# Merge and compare
tag_check <- tag_ranges[actual_counts, on = "Tag"]
tag_check[, complete := actual_count == expected_count]
# Summary
tag_check[, .N, by = complete]

# See tags with missing censuses
missing_tags <- tag_check[complete == FALSE]
missing_tags[, gap := expected_count - actual_count]

# Summary stats
cat("Total tags:", nrow(tag_check), "\n")
cat("Complete tags:", tag_check[complete == TRUE, .N], "\n")
cat("Tags with gaps:", tag_check[complete == FALSE, .N], "\n")
if (nrow(missing_tags) > 0) {
  cat("Total missing observations:", sum(missing_tags$gap), "\n")
}

# ---- 6. Detect potential DBH measurement errors ----
# Method: compute temporal log-differences within stems (d_prev, d_next, d_span)
# and flag rows where DBH deviates sharply from both previous and next measurements
# but the span is small (likely data-entry typo that corrects in subsequent census).
ViewFullTable_hom_corrected_clean

# Check example tags where stem IDs were reassigned
inc <- c("001112")

# NOTE: log-DBH-difference checks for data-entry errors are reliable only
# for stems with a consistent StemID across censuses. Stems whose StemIDs
# change cannot be tracked confidently and are excluded from that check.
# This correction feeds the matching algorithm; user-facing DBH values in
# the R tables are not modified.

ViewFullTable_hom_corrected_clean[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM)
][!is.na(DBH)]

# ---- 6.1 DBH introduction errors (log-difference method) ----
# We compute log-transformed differences (d_prev, d_next, d_span) for each
# stem and flag entries satisfying multiple criteria (1.5x, 3x or data-driven
# thresholds). Rows flagged as `entry_error_any == TRUE` are suspicious and will
# receive `dbh_candidate` as the geometric mean of neighbors where available.
#
# Step 0: Separate valid DBH and NA DBH
# --
valid_DBH <- ViewFullTable_hom_corrected_clean[!is.na(DBH)]
NA_DBH <- ViewFullTable_hom_corrected_clean[is.na(DBH)]

nrow(valid_DBH) + nrow(NA_DBH) == nrow(ViewFullTable_hom_corrected_clean) # should be TRUE

# --
# 1) Sort & index
# --
setkey(valid_DBH, Tag, StemTag, TreeID, StemID, CensusID)

# --
# 2) Compute log directly in shift (no temp col)
# --
valid_DBH[, log_DBH := log(DBH)]

# Rationale: shifting log(DBH) is necessary to identify multiplicative data-entry errors
# Once corrected, DBH values can be tracked without log transformation
valid_DBH[, `:=`(
  log_prev = shift(log_DBH, type = "lag"),
  log_next = shift(log_DBH, type = "lead")
), by = .(Tag, StemTag, TreeID, StemID)]

# --
# 3) Cleanup
# --
# Compute differences using cached shifts
valid_DBH[, `:=`(
  d_prev = log_DBH - log_prev,
  d_next = log_next - log_DBH,
  d_span = log_next - log_prev
)]

# Verify computed differences for example Tag
valid_DBH[Tag == "000006", .(CensusID, Tag, StemTag, TreeID, StemID, DBH, log_DBH, log_prev, log_next, d_prev, d_next, d_span)][order(CensusID)]
valid_DBH[Tag == "001112", .(CensusID, Tag, StemTag, TreeID, StemID, DBH, log_DBH, log_prev, log_next, d_prev, d_next, d_span)][order(CensusID)]
# The example shows the differences are computed correctly

# --
# Step 3: Compute data-driven threshold
# --
thr_data <- quantile(abs(valid_DBH$d_prev), 0.99, na.rm = TRUE)
# --
# Step 4: Row-level error flag
# --
threshold_log1p5 <- log(1.5)
threshold_log3 <- log(3)

# Main idea:
# Normal growth: d_prev and d_next are similar, small positive values
# Data error: d_prev is huge jump up, d_next is huge jump down (or vice versa), but d_span is normal

# Flags a measurement as an error if **ALL** of these are true:
# 1. Has both previous and next measurements (not NA)
# 2. **Big jump in** (large `d_prev`) AND **big jump out** (large `d_next`)
# 3. BUT the **span is normal** (small `d_span`)
valid_DBH[, entry_error_any := !is.na(d_prev) & !is.na(d_next) & (
  (abs(d_prev) > threshold_log1p5 & abs(d_next) > threshold_log1p5 & abs(d_span) <= threshold_log1p5) |
    (abs(d_prev) > threshold_log3 & abs(d_next) > threshold_log3 & abs(d_span) <= threshold_log3) |
    (abs(d_prev) > thr_data & abs(d_next) > thr_data & abs(d_span) <= thr_data)
)]

# Why log transformation?
# In log-space, multiplicative errors are symmetric.
# Without logs, dividing errors look smaller than multiplying errors, making detection harder.

# --
# Step 4a: Numeric error score
# --
valid_DBH[, error_score := fifelse(is.na(d_prev) & is.na(d_next), NA_real_, pmax(abs(d_prev), abs(d_next), na.rm = TRUE))]
# --
# Step 5: Remove intermediate columns
# --
valid_DBH[, c("log_DBH", "d_prev", "d_next", "d_span") := NULL]
# --
# Step 6: Stem-level summary with fix for all-NA error_score
# --
# stem_summary <- valid_DBH[
#   , .(
#     stem_has_any_error = any(entry_error_any, na.rm = TRUE),
#     stem_max_error_score = max(error_score, na.rm = TRUE)
#   ),
#   by = .(Tag, StemTag, TreeID, StemID)
# ]
# stem_summary[is.infinite(stem_max_error_score), stem_max_error_score := NA_real_]

# ---- 6.2 Stem-level summary & scoring ----
# Separate stems with single vs multiple measurements and summarize
# `entry_error_any` and `stem_max_error_score` to prioritize manual review.
row_counts <- valid_DBH[, .N, by = .(Tag, StemTag, TreeID, StemID)]

# Split into single-row and multi-row groups
single_row_data <- valid_DBH[row_counts[N == 1], on = .(Tag, StemTag, TreeID, StemID)]
multi_row_data <- valid_DBH[row_counts[N > 1], on = .(Tag, StemTag, TreeID, StemID)]

# Summarize each separately
single_summary <- copy(single_row_data)[
  , .(Tag, StemTag, TreeID, StemID)
]
single_summary[, stem_has_any_error := FALSE] # single-row stems cannot have errors
single_summary[, stem_max_error_score := NA_real_] # set to NA for single-row stems

setkey(multi_row_data, Tag, StemTag, TreeID, StemID)

multi_summary <- multi_row_data[
  , .(
    stem_has_any_error = any(entry_error_any, na.rm = TRUE),
    stem_max_error_score = max(error_score, na.rm = TRUE)
  ),
  by = key(multi_row_data) # groups by key
]

# Replace Inf with NA
single_summary[is.infinite(stem_max_error_score), stem_max_error_score := NA_real_]
multi_summary[is.infinite(stem_max_error_score), stem_max_error_score := NA_real_]

# Optional: combine both summaries
stem_summary <- rbind(single_summary, multi_summary)

# Assign in-place
setkey(valid_DBH, Tag, StemTag, TreeID, StemID)
valid_DBH[stem_summary, `:=`(
  stem_has_any_error = i.stem_has_any_error,
  stem_max_error_score = i.stem_max_error_score
)]
# --
# Step 7: Handle NA DBH rows explicitly with NA_real_
# --
NA_DBH[, `:=`(
  entry_error_any = NA,
  error_score = NA_real_,
  stem_has_any_error = NA,
  stem_max_error_score = NA_real_
)]
# --
# Step 8: Recombine
# --

setdiff(names(valid_DBH), names(NA_DBH))

# NA_DBH[, log_DBH := NA_real_]
NA_DBH[, log_prev := NA_real_]
NA_DBH[, log_next := NA_real_]
# NA_DBH[, d_prev := NA_real_]
# NA_DBH[, d_next := NA_real_]
# NA_DBH[, d_span := NA_real_]

ViewFullTable_measurement_error_indication <- rbindlist(list(valid_DBH, NA_DBH), use.names = TRUE)
# setorder(ViewFullTable_measurement_error_indication, Tag, StemTag, TreeID, StemID, CensusID)
setorder(ViewFullTable_measurement_error_indication, RowIDN1)

# how many errors detected?
table(ViewFullTable_measurement_error_indication$entry_error_any, useNA = "ifany")

# Only stems needing review
unique(ViewFullTable_measurement_error_indication[stem_has_any_error == TRUE, .(Tag, StemTag, TreeID, StemID)])

# Add candidate DBH values using geometric mean of neighbors
# Compute geometric mean of previous and next valid DBH when available
ViewFullTable_measurement_error_indication[, dbh_candidate := fifelse(
  entry_error_any & !is.na(log_prev) & !is.na(log_next),
  exp((log_prev + log_next) / 2),
  NA_real_
)]

# ---- 8. Prepare DBH candidates (for taper correction) ----
# Prepare data to correct for possible measurement errors prior to taper
# correction. For rows flagged `entry_error_any == TRUE`, use `dbh_candidate`
# when available, otherwise keep original DBH. This produces
# `dbh_with_best_candidate` used for taper correction.
ViewFullTable_measurement_error_indication[, dbh_with_best_candidate := fifelse(
  entry_error_any & !is.na(dbh_candidate),
  dbh_candidate,
  DBH
)]

ViewFullTable_measurement_error_indication[stem_has_any_error == TRUE][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, dbh_candidate, dbh_with_best_candidate)
][order(Tag, StemTag, TreeID, StemID, CensusID)]

# ---- 8.1 Pre-taper check: missing Tag × CensusID combinations ----
# Quick pre-taper check: compute expected vs actual counts per Tag and identify
# any gaps before taper correction is applied.
xraw_unique <- unique(ViewFullTable_measurement_error_indication[, .(Tag, CensusID)])
# Get the range per tag
tag_ranges <- xraw_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]
# Add expected count (how many censuses should exist)
tag_ranges[, expected_count := max_c - min_c + 1L]
# Get actual count per tag
actual_counts <- xraw_unique[, .(actual_count = .N), by = Tag]
# Merge and compare
tag_check <- tag_ranges[actual_counts, on = "Tag"]
tag_check[, complete := actual_count == expected_count]
# Summary
tag_check[, .N, by = complete]

# See tags with missing censuses
missing_tags <- tag_check[complete == FALSE]
missing_tags[, gap := expected_count - actual_count]

# Summary stats
cat("Total tags:", nrow(tag_check), "\n")
cat("Complete tags:", tag_check[complete == TRUE, .N], "\n")
cat("Tags with gaps:", tag_check[complete == FALSE, .N], "\n")
if (nrow(missing_tags) > 0) {
  cat("Total missing observations:", sum(missing_tags$gap), "\n")
}

# Note: review `missing_tags` to decide whether to add NA rows prior to taper correction

# ---- 8.2 Apply taper correction (DBH → DBHC) ----
# Prepare `HOM_for_taper_correction` and run taper correction to compute DBH at 1.3 m (DBHC).
# Inputs: `dbh_with_best_candidate` (numeric) and `HOM_for_taper_correction`.
# The function `apply_taper_correction()` (in `taper_correction.R`) returns a column
# named by `output_col` (here: `dbh_with_best_candidate_taper_corrected`).
# Create `HOM_for_taper_correction` column for taper correction
ViewFullTable_measurement_error_indication[, HOM_for_taper_correction := HOM]

# Load taper utilities (provides `apply_taper_correction()` and `taper()`)
source(here::here("BCI_stem_reconstruction", "1_DATA_PREPARATION", "HELPER_FUNCTIONS", "taper_correction.R"))

# Inspect HOM_for_taper_correction values
range(ViewFullTable_measurement_error_indication$HOM_for_taper_correction, na.rm = TRUE)
unique(ViewFullTable_measurement_error_indication$HOM_for_taper_correction)
table(sort(ViewFullTable_measurement_error_indication$HOM_for_taper_correction), useNA = "ifany")

# Round HOM to 3 decimals
ViewFullTable_measurement_error_indication[, HOM_for_taper_correction := round(HOM_for_taper_correction, 3)]

# After rounding, inspect again
range(ViewFullTable_measurement_error_indication$HOM_for_taper_correction, na.rm = TRUE)
unique(ViewFullTable_measurement_error_indication$HOM_for_taper_correction)
table(sort(ViewFullTable_measurement_error_indication$HOM_for_taper_correction), useNA = "ifany")

# Check rows where HOM_for_taper_correction == 0 and DBH availability
ViewFullTable_measurement_error_indication[!is.na(dbh_with_best_candidate) & HOM_for_taper_correction == 0]
# All rows with HOM_for_taper_correction == 0 have NA DBH, so set those HOM_for_taper_correction values to NA
table(ViewFullTable_measurement_error_indication[HOM_for_taper_correction == 0]$DBH, useNA = "ifany")
ViewFullTable_measurement_error_indication[HOM_for_taper_correction == 0, HOM_for_taper_correction := NA_real_]

# Re-check range after conversion
range(ViewFullTable_measurement_error_indication$HOM_for_taper_correction, na.rm = TRUE)

# (quantile(ViewFullTable_measurement_error_indication$dbh_with_best_candidate, na.rm = TRUE) / 10) / 100

## lOAD SPECIES DATA
load("./BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.RData")

# ## Load growth forms
bci.spptable[, Lifeform := tolower(Lifeform_RPerez_SAguilar)]
growth_forms <- bci.spptable[, .(Mnemonic, Lifeform)]

species_to_use_tapper_corrected_dbh <- growth_forms[
  is.na(Lifeform) | grepl(pattern = "árbol|arbusto", x = Lifeform)
]$Mnemonic

ViewFullTable_taper_corrected <- apply_taper_correction(ViewFullTable_measurement_error_indication,
  # dbh input in mm
  dbh_col = "dbh_with_best_candidate",
  hom_col = "HOM_for_taper_correction",
  wsg_col = NULL,
  output_col = "dbh_with_best_candidate_taper_corrected_raw",
  taper_correction = TRUE,
  common_hom = 1.3,
  convert_units = FALSE, # Disable unit conversion (BCI hom units appear fine)
  verbose = TRUE
)

ViewFullTable_taper_corrected[, dbh_with_best_candidate_taper_corrected := fifelse(
  Mnemonic %in% species_to_use_tapper_corrected_dbh,
  dbh_with_best_candidate_taper_corrected_raw,
  dbh_with_best_candidate
)][, dbh_with_best_candidate_taper_corrected_raw := NULL]

# check whether number stems that where corrected
ViewFullTable_taper_corrected[HOM_for_taper_correction == 1.3 &
  dbh_with_best_candidate_taper_corrected != dbh_with_best_candidate]
# Only HOM != 1.3 were corrected
# Those 83 rows differ only by tiny decimals; rounding fixes most cases

# Check example tags where stem IDs were reassigned
inc <- c("001112")

ViewFullTable_taper_corrected[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(
    Tag, StemTag, TreeID, StemID, CensusID, DBH,
    dbh_with_best_candidate, dbh_with_best_candidate_taper_corrected,
    HOM_for_taper_correction, HOM, ListOfTSM
  )
][!is.na(dbh_with_best_candidate)]

# Order `ViewFullTable_taper_corrected` for inspection
setorder(ViewFullTable_taper_corrected, RowIDN1)

# Show a small sample for quick inspection
ViewFullTable_taper_corrected[1:100, .(
  Tag, StemTag, TreeID, StemID, CensusID,
  dbh_with_best_candidate, dbh_with_best_candidate_taper_corrected,
  HOM_for_taper_correction, ListOfTSM
)]

# ---- 9. Finalize complete table ----

select_cols <- c(
  intersect(
    names(ViewFullTable_taper_corrected),
    names(ViewFullTable)
  ),
  # "log_prev",
  # "log_next",
  # "stem_has_any_error",
  "dbh_with_best_candidate_taper_corrected"
)

ViewFullTable_taper_corrected <- ViewFullTable_taper_corrected[, ..select_cols]

# make new rowid
ViewFullTable_taper_corrected[, RowID := .I]

#### CHECK SINGLE VS MULTIPLE-STEMS TAGS ----
## ---- Mark single stem tags ------------------------------------------------------
id_single_stem_tags <- ViewFullTable_taper_corrected[
  , .(
    all_stemtag_na = all(is.na(StemTag)),
    one_stemid = uniqueN(StemID[!is.na(StemID)]) == 1
  ),
  by = Tag
]

table(id_single_stem_tags$all_stemtag_na)
table(id_single_stem_tags$one_stemid)
table(id_single_stem_tags$all_stemtag_na & id_single_stem_tags$one_stemid)

# # Do the counts match the number of unique tags?
# sum(table(id_single_stem_tags$all_stemtag_na)) == length(unique(ViewFullTable_taper_corrected$Tag))
# sum(table(id_single_stem_tags$one_stemid)) == length(unique(ViewFullTable_taper_corrected$Tag))
# sum(table(id_single_stem_tags$all_stemtag_na & id_single_stem_tags$one_stemid)) == length(unique(ViewFullTable_taper_corrected$Tag))

# get percentage of tags with one stem or multiple stems
tags_info <- round(table(id_single_stem_tags$all_stemtag_na & id_single_stem_tags$one_stemid) / nrow(id_single_stem_tags) * 100, 1)
# FALSE  TRUE
#  28.7  71.3

## Export percentage of single vs multiple stem tags
single_stem_percentage <- tags_info["TRUE"]
multiple_stem_percentage <- tags_info["FALSE"]

# select tags with all StemTag NA but only one unique StemID
# These have been already identified as single individual
tags_with_one_stemid_no_stemtag <- id_single_stem_tags[all_stemtag_na == TRUE & one_stemid == TRUE, Tag]
length(tags_with_one_stemid_no_stemtag)

# Two key point to consider:
# 1) For these Tags, do all NA in stemid have DBH?
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]
# There are 81 cases where StemID is NA and DBH is NA too. This is fine

# 2) count census per stemid for these tags to confirm they are
# single-stemmed across all censuses. If any have multiple StemIDs across
# censuses, that would be unexpected and worth investigating.
unique(ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag, .N, by = .(Tag, StemID, CensusID)]$N)
# All N are 1, which is consistent with single-stemmed tags having one StemID across all censuses.

# 3) Do all DBH have stemid?
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & !is.na(DBH) & is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]
# aLL DBH have StemID except for 81 cases where both DBH and StemID are NA.

# 4) What about the NA-DBH?
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & is.na(DBH) & !is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]
# 240595 rows with NA-DBH with StemID.
# There are several cases with NA DBH values and valid StemIDs. The final R tables will
# correct this, since if there is no DBH recorded in a census (e.g., the last census),
# that census will be treated as the last observed one, and the individual may be considered dead.
# Ultimately, the last observation (i.e., death) is determined by the presence of a DBH value,
# not by the presence of a StemID when DBH is NA.

# Quickly inspect the key variables
inspectdf::inspect_na(ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag, .(Tag, StemTag, TreeID, StemID, CensusID, DBH)])

## ---- Mark multi-stemmed tags that need to run DP ------------------------------------------------------

single_stem_tags <- id_single_stem_tags[Tag %in% tags_with_one_stemid_no_stemtag]$Tag
multiple_stem_tags <- id_single_stem_tags[!Tag %in% tags_with_one_stemid_no_stemtag]$Tag

# ! SANITY CHECKS
# Are the counts (single + multiple) correct?
(length(single_stem_tags) + length(multiple_stem_tags)) ==
  length(unique(ViewFullTable_taper_corrected$Tag)) # should be TRUE

# Are the number of rows correct between filterres observations?
nrow(ViewFullTable_taper_corrected[Tag %in% single_stem_tags, .(Tag, CensusID, StemTag, StemID, DBH, dbh_with_best_candidate_taper_corrected)]) +
  nrow(ViewFullTable_taper_corrected[Tag %in% multiple_stem_tags, .(Tag, CensusID, StemTag, StemID, DBH, dbh_with_best_candidate_taper_corrected)]) ==
  nrow(ViewFullTable_taper_corrected)

# Are Tags repeated in single vs multiple stem groups?
intersect(as.character(single_stem_tags), as.character(multiple_stem_tags)) # should be character(0)
# No, there are no tags in common between the single and multiple stem groups,
# which is consistent with our classification.

## ---- Create column indicating single vs. multiple stem tags ------------------------------------------------------
ViewFullTable_single_vs_multiple_stem_tags <- copy(ViewFullTable_taper_corrected)
ViewFullTable_single_vs_multiple_stem_tags[, single_stem_tags := Tag %in% single_stem_tags]
setorder(ViewFullTable_single_vs_multiple_stem_tags, RowID)
rm(ViewFullTable_taper_corrected)
gc()

# ! SANITY CHECKS
length(unique(ViewFullTable_single_vs_multiple_stem_tags[single_stem_tags == TRUE]$Tag)) == length(single_stem_tags)
length(unique(ViewFullTable_single_vs_multiple_stem_tags[single_stem_tags == FALSE]$Tag)) == length(multiple_stem_tags)

table(ViewFullTable_single_vs_multiple_stem_tags$single_stem_tags, useNA = "ifany")

# # ------------------------------------------------------------------
# # Section: trajectory comparison against Dryad 'condit' stemIDs
# # ------------------------------------------------------------------
# # Summarize DBH by Tag and stem identifier in both datasets and check if
# # the sorted sums match.  If they differ, the reconstructed trajectory
# # deviates from the original assignment.
# # remove indices for faster processing
# setindex(ViewFullTable_single_vs_multiple_stem_tags, NULL) # removes all indices
# indat_compare <- copy(ViewFullTable_single_vs_multiple_stem_tags)
# table(indat_compare$single_stem_tags)
# indat_compare <- indat_compare[CensusID <= 8]
# indat_compare <- indat_compare[Tag %in% tags_with_one_stemid_no_stemtag]
# indat_compare <- indat_compare[, Tag := as.character(Tag)][, StemID := as.numeric(as.character(StemID))]

# # ! Check how condit did when codes are involved
# # add DBH per Tag and ReconstructedStemID
# quantile(indat_compare$DBH, na.rm = TRUE)
# quantile(bci_dryad$DBH, na.rm = TRUE)

# indat_compare_summary <- indat_compare[!is.na(DBH), sum(DBH, na.rm = TRUE), by = .(Tag, StemID)]
# indat_compare_summary[, TotalDBH := V1]
# indat_compare_summary[, V1 := NULL]

# dryad_summary <- bci_dryad[!is.na(DBH), sum(DBH, na.rm = TRUE), by = .(Tag, StemID)]
# dryad_summary[, TotalDBH_Condit := V1]
# dryad_summary[, V1 := NULL]

# # for each StemID in dryad summary, sort TotalDBH_Condit. do the same for each
# # ReconstructedStemID in indat summary. Then, compare the sorted TotalDBH_Condit
# # with the sorted TotalDBH for each tag. If they are the same, then the
# # trajectories are the same, if they are different, then the trajectories are
# # different.

# # Pre-split data (fast!)
# dryad_split <- split(dryad_summary[Tag %in% tags_with_one_stemid_no_stemtag], by = "Tag", keep.by = FALSE)
# indat_split <- split(indat_compare_summary[Tag %in% tags_with_one_stemid_no_stemtag], by = "Tag", keep.by = FALSE)

# length(dryad_split)
# length(indat_split)

# compare_differences <- function(tag) {
#   dryad_tag <- dryad_split[[as.character(tag)]]
#   indat_tag <- indat_split[[as.character(tag)]]
#   if (is.null(dryad_tag) || is.null(indat_tag)) {
#     return(NA)
#   }
#   setorder(dryad_tag, TotalDBH_Condit)
#   setorder(indat_tag, TotalDBH)
#   all(dryad_tag$TotalDBH_Condit == indat_tag$TotalDBH)
# }

# setkey(dryad_summary, Tag)
# setkey(indat_compare_summary, Tag)

# comparison_results <- dryad_summary[indat_compare_summary,
#   .(equal = all(sort(TotalDBH_Condit) == sort(TotalDBH))),
#   by = .EACHI,
#   on = "Tag"
# ]
# (table(comparison_results$equal, useNA = "ifany") / nrow(comparison_results)) * 100
# # TRUE
# #  100

# ## ALL SINGLE STEM TAGS ARE THE SAME AS CONDIT, WHICH IS EXPECTED SINCE THEY WERE NOT CHANGED IN THE RECONSTRUCTION.
# ## THIS PROCESS SUGGEST THAT DP RECONSTRCTION SHOULD FOCUS ON THE MULTIPLE-STEMMED TAGS.

# ViewFullTable_single_vs_multiple_stem_tags
# setorder(ViewFullTable_single_vs_multiple_stem_tags, RowID)

# # SAVE SINGLE VS MULTI-STEM COLUMN TABLE
# saveRDS(
#   ViewFullTable_single_vs_multiple_stem_tags,
#   here::here("DATA", "PROCESSED", "5_ViewFullTable_single_vs_multiple_stem_tags.rds")
# )

# ---------------------------------------------------------------------------
# merge the categorical labels back into the full observation table so we can
# quantify how many records fall into each bucket.  This also reveals any
# mnemonic mismatches that produced NA values.
ViewFullTable_single_vs_multiple_stem_tags <- merge(ViewFullTable_single_vs_multiple_stem_tags, unique(growth_forms[, .(Mnemonic, Lifeform)]), by = "Mnemonic", all.x = TRUE)
setorder(ViewFullTable_single_vs_multiple_stem_tags, RowID)

# Summaries
ViewFullTable_single_vs_multiple_stem_tags[, .N, by = Lifeform][order(-N)]

unique(ViewFullTable_single_vs_multiple_stem_tags[is.na(Lifeform), .(Mnemonic, Genus, SpeciesName)])

tags_per_growth_form <- unique(ViewFullTable_single_vs_multiple_stem_tags[, .(Tag, single_stem_tags, Lifeform)])
# check individuals epr growth form with multiple vs single stem
tags_per_growth_form[
  , .N,
  by = .(single_stem_tags, Lifeform)
][order(single_stem_tags, -N)]

# save the enriched observation table (typo in filename preserved for legacy)
saveRDS(ViewFullTable_single_vs_multiple_stem_tags, "./BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_taper_corrected_growth_forms.rds")

# make sure SpeciesName in ViewFullTable_single_vs_multiple_stem_tags matches your full "Genus species" format
nobs_growth_form <- ViewFullTable_single_vs_multiple_stem_tags[!is.na(DBH)][
  , .N,
  by = .(Tag, SpeciesName, Lifeform)
]

nobs_growth_form[
  , .N,
  by = Lifeform
][order(-N)]

# ---- 10. Check all tags x censusid are complete ----
# Compute expected vs actual counts per Tag and identify any gaps (tags with missing censuses)
# using the taper-corrected table. Add missing Tag × CensusID rows below if needed.
ViewFullTable_single_vs_multiple_stem_tags_unique <- unique(ViewFullTable_single_vs_multiple_stem_tags[, .(Tag, CensusID)])
# Get the range per tag
tag_ranges <- ViewFullTable_single_vs_multiple_stem_tags_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]
# Add expected count (how many censuses should exist)
tag_ranges[, expected_count := max_c - min_c + 1L]
# Get actual count per tag
actual_counts <- ViewFullTable_single_vs_multiple_stem_tags_unique[, .(actual_count = .N), by = Tag]
# Merge and compare
tag_check <- tag_ranges[actual_counts, on = "Tag"]
tag_check[, complete := actual_count == expected_count]
# Summary
tag_check[, .N, by = complete]

# See tags with missing censuses
missing_tags <- tag_check[complete == FALSE]
missing_tags[, gap := expected_count - actual_count]

# Summary stats
cat("Total tags:", nrow(tag_check), "\n")
cat("Complete tags:", tag_check[complete == TRUE, .N], "\n")
cat("Tags with gaps:", tag_check[complete == FALSE, .N], "\n")
if (nrow(missing_tags) > 0) {
  cat("Total missing observations:", sum(missing_tags$gap), "\n")
}
