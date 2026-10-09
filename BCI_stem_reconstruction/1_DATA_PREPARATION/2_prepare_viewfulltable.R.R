# =============================================================================
# 1_prepare_viewfulltable.R.R
#
# Purpose: Load and prepare the BCI ViewFullTable for stem identification.
#          Give every tag a row at every census between its first and last
#          record, keep one measurement per stem and census (highest HOM),
#          replace likely DBH entry errors, apply the taper correction, and
#          label single- vs multiple-stem tags.
#
# Run from the project root: the species table is loaded with a path relative
# to the working directory.
#
# Inputs:
#   BCI_stem_reconstruction/DATA/RAW/ViewFiles_bci_allcensuses/ViewFullTable_bci.csv
#     (tab-separated; DBH in mm, HOM in m)
#   BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.RData
#     (from 0_prepare_species_tables.R; growth forms)
# Output:
#   BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
#     the ViewFullTable columns (CensusID renumbered 1..n by mean census date,
#     DBH as recorded) plus
#       dbh_with_best_candidate_taper_corrected  DBH (mm) after the entry-error
#                         replacement and, for species with a tree or shrub
#                         growth form (or none), the taper correction to 1.3 m
#       single_stem_tags  TRUE for tags with one StemID and no StemTag
#       Lifeform          growth form of the species (lower case, as in the
#                         species list)
#       RowID             row number
# Packages: data.table, ggplot2, inspectdf.
# =============================================================================

# Clear all objects from the workspace to avoid accidental contamination
rm(list = ls())

library(data.table)
library(ggplot2) # diagnostic plots
# References: https://arelbundock.com/posts/dt_tb_df/index.html
#             https://rdatatable.gitlab.io/data.table/

# ---- 1. Configuration ----
# User-editable variables and file locations. Ensure the paths below point to
# the expected tab-delimited `ViewFullTable` file for the site.

site <- "bci" # Site code for BCI

# Input folder (raw ViewFullTable export)
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
  workspace_root <- dirname(workspace_root)
}

INPUT_folder_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses")

# Output and diagnostics folders
OUTPUT_folder <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "PROCESSED")
if (!dir.exists(OUTPUT_folder)) {
  dir.create(OUTPUT_folder, recursive = TRUE)
}

CHECK_folder <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS")
if (!dir.exists(CHECK_folder)) {
  dir.create(CHECK_folder, recursive = TRUE)
}

# End of configuration (adjust above values as needed for local runs)

# ---- 2. Load and verify the input data ----
# Read the BCI ViewFullTable from the raw data folder. The table must contain
# measurement columns like `ExactDate`, `Tag`, `StemTag`, `TreeID`, `StemID`, and `DBH`.

file1 <- file.path(INPUT_folder_1, paste0("ViewFullTable_", site, ".csv"))

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

# Check script BCI_stem_reconstruction/1_DATA_PREPARATION/0_prepare_species_tables.R:
# Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
unique(ViewFullTable[SpeciesName == "rohrii", .(Mnemonic, Family, Genus, SpeciesName)])
nrow(ViewFullTable[SpeciesName %in% "rohrii"])

# the rows with code pterof (36) need to be replaced with the rohrii code
unique(ViewFullTable[SpeciesName == "officinalis", .(Mnemonic, Family, Genus, SpeciesName)])
nrow(ViewFullTable[SpeciesName %in% "officinalis"])

# Correct
ViewFullTable[Mnemonic == "pterro"]
# Incorrect
ViewFullTable[Mnemonic == "pterof"]
# Replace
ViewFullTable[Mnemonic == "pterof", Mnemonic := "pterro"]
# Check
ViewFullTable[Mnemonic == "pterof"]

# length(unique(ViewFullTable$Mnemonic))

# Rename CensusID to CensusID_raw to preserve original labels, then assign
# standardized sequential CensusID values below.
setnames(ViewFullTable, "CensusID", "CensusID_raw")
# Identifier columns as factors
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

# Assign sequential CensusID values (censuses ordered by their mean date)
censusid_dates[, CensusID := seq_len(.N)]

setkey(ViewFullTable, CensusID_raw)
setkey(censusid_dates, CensusID_raw)

# Join new CensusID into main table
ViewFullTable <-
  censusid_dates[, .(CensusID_raw, CensusID)][
    ViewFullTable,
    on = "CensusID_raw"
  ]

ViewFullTable[, CensusID_raw := NULL]

# ---- 3. Fill missing Tag × CensusID rows ----
# For every tag, build the grid of all censuses between its first and its last
# record and insert a placeholder row (no stem, no DBH) for each census in
# which the tag has no row, so that every tag has a complete panel.

# DBH unit check first (the panel fill starts further below)

sort(unique(round(ViewFullTable[!is.na(DBH)]$DBH, 1)))

plot(density(ViewFullTable[!is.na(DBH)]$DBH))

# Explore DBH units by comparing quantiles under different assumptions.
# The values confirm that DBH is recorded in millimetres.
check_diameter_units <- sort(ViewFullTable[!is.na(DBH)]$DBH)

head(check_diameter_units, 100) # mm
tail(check_diameter_units, 100) ## check next lines

ViewFullTable[DBH == "8169"]

quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 10 # convert to cm from mm # GOOD
quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 10 / 100 # convert to m from mm # GOOD
quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 100 # convert to m from cm #! WRONG

# Tag and CensusID have no missing values, so they can safely define the complete panel.

# Tag and TreeID are redundant; use Tag as the primary identifier for data preparation.
ViewFullTable[, .(n_tags = uniqueN(Tag), n_treeids = uniqueN(TreeID))]

tag_census_unique <- unique(copy(ViewFullTable[, .(Tag, CensusID)]))

# Create complete grid of all Tag-CensusID combinations
tag_ranges <- tag_census_unique[, .(min_c = min(CensusID), max_c = max(CensusID)), by = Tag]

## Every census between the first and last record of each tag
complete_grid <- tag_ranges[, .(CensusID = seq.int(min_c, max_c)), by = Tag]

# Find which combinations are missing from the original data
missing_combinations <- complete_grid[!ViewFullTable,
  on = .(Tag, CensusID)
]

unique(ViewFullTable[, .(Tag, TreeID)])

# There are missing Tag × CensusID rows that must be added to complete the panel.
nrow(missing_combinations)

# Build placeholder rows for the missing combinations and fill metadata columns.

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
missing_combinations

# Tag-level attributes, taken from the first row of each tag
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

# ExactDate and Date of a placeholder row: a date recorded for its quadrat in
# that census (joined on QuadratID + QuadratName + CensusID)
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
cat("Number of missing Tag-CensusID combinations after adding missing rows:", nrow(missing_combinations_final), "\n")

# Carry forward the complete Tag × CensusID table; drop auxiliary column; assign row index.
ViewFullTable <- ViewFullTable_no_missing_tags_census
ViewFullTable <- ViewFullTable[, missing_tag_census := NULL]
ViewFullTable[, RowID := .I]

# Verify all tags have contiguous observations between their first and last census.
ViewFullTable[, .(CensusID = list(sort(unique(CensusID)))), by = Tag][
  , expected := .(list(seq(min(unlist(CensusID)), max(unlist(CensusID))))),
  by = Tag
][
  , complete := identical(CensusID[[1]], expected[[1]])
][
  , .N,
  by = complete
]
# Confirmed: all tags now have complete observations between their first and last census.

# ---- Stem identification: notes & assumptions ----
# Variables and assumptions used in later stem matching and QC steps:
#   - Tag: plot + visible tag — primary grouping identifier across censuses
#   - CensusID: integer temporal index (1,2,3…) used to order measurements
#   - DBH: measured diameter at the recorded HOM (may be NA if not measured)
#   - HOM: height of measurement (metres) — used to select one measurement per census and for taper correction in buttressed trees
#   - Species / taxonomy: helps limit candidate matches when matching stems

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

# ---- Stem ID history and matching context ----
# StemIDs change across censuses for some tags, so stage 2 re-identifies the stems from the taper-corrected DBH.

# Example tags with retroactive StemID reassignment:
inc <- c("001112", "003036")

ViewFullTable[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM)
][!is.na(DBH)]

## plot_stem(): DBH against census for the tags in `tag` (a character vector),
## one line per StemID and one panel per tag. Returns a ggplot object.
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

# Section 4 keeps the row with the highest HOM of each (Tag, StemID, CensusID) group.
# The taper correction (section 8.2) then standardizes DBH to HOM = 1.3 m for the matching in stage 2.

# Some tags kept original StemIDs because stems were not clearly distinguishable.
inc <- c("151991")
plot_stem(ViewFullTable, inc)

ViewFullTable[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM, Mnemonic)
][!is.na(DBH)]

# ---- 4. HOM correction: select max HOM per group ----
# For groups with multiple DBH/HOM records per (Tag × Stem × CensusID), keep
# the row with the highest HOM (preferred measurement point for buttressed trees).
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

# ---- 5. Recompute QC using corrected HOM ----
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

# One remaining tied max HOM, select row with largest DBH (129 vs 130 DBH, 1 mm difference)
ViewFullTable_hom_corrected[QC_is_tied_max_HOM == TRUE]
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

ViewFullTable_hom_corrected_clean
ViewFullTable_hom_assessment

# check differences in row berween `ViewFullTable_hom_assessment` and `ViewFullTable_hom_corrected_clean`
diff_rows <- fsetdiff(
  ViewFullTable_hom_assessment[, names(ViewFullTable), with = FALSE],
  ViewFullTable_hom_corrected_clean[, names(ViewFullTable), with = FALSE]
)

ViewFullTable_hom_corrected_clean[Tag %in% diff_rows$Tag[1]]
ViewFullTable_hom_assessment[Tag %in% diff_rows$Tag[1]]
# in census 4, we selected one with the highest HOM, which is the correct one for the R tables

# Remove intermediate data.tables to free memory
rm(ViewFullTable_hom_assessment)
rm(multi_qc)
rm(single_obs)
rm(ViewFullTable_hom_corrected)
gc()

# Include rowid for tracking
ViewFullTable_hom_corrected_clean[, RowIDN1 := .I]

# ---- 6. Check missing Tag × CensusID combinations ----
# Compute expected vs. actual counts per Tag and verify no gaps remain after HOM correction.
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

# ---- 7. Detect DBH measurement errors ----
# Method: compute temporal log-differences within stems (d_prev, d_next, d_span)
# and flag rows where DBH deviates sharply from both previous and next measurements
# but the span is small (likely a data-entry typo that self-corrects in the next census).
ViewFullTable_hom_corrected_clean

# Check example tags where stem IDs were reassigned
inc <- c("001112")

# The differences are computed within one (Tag, StemTag, TreeID, StemID)
# series, so a stem whose StemID changes between censuses is checked as
# separate series. The recorded `DBH` column is not altered: the replacement
# value goes to `dbh_with_best_candidate` (section 8), which feeds the taper
# correction.

ViewFullTable_hom_corrected_clean[Tag %in% inc][
  order(Tag, StemID, CensusID)
][
  , .(Tag, StemTag, TreeID, StemID, CensusID, DBH, HOM, ListOfTSM)
][!is.na(DBH)]

# ---- 7.1 DBH outlier detection (log-difference method) ----
# Compute log-transformed differences (d_prev, d_next, d_span) for each stem
# and flag entries that meet any of three criteria (a factor of 1.5, a factor
# of 3, or the data-driven threshold).
# Rows flagged as `entry_error_any == TRUE` receive `dbh_candidate`, the
# geometric mean of their two neighbors.
#
# Step 0: Separate valid DBH and NA DBH
# --
valid_DBH <- ViewFullTable_hom_corrected_clean[!is.na(DBH)]
NA_DBH <- ViewFullTable_hom_corrected_clean[is.na(DBH)]

nrow(valid_DBH) + nrow(NA_DBH) == nrow(ViewFullTable_hom_corrected_clean) # should be TRUE

# --
# Step 1: Sort & index
# --
setkey(valid_DBH, Tag, StemTag, TreeID, StemID, CensusID)

# --
# Step 2: log DBH and its previous / next value within each stem series
# --
valid_DBH[, log_DBH := log(DBH)]

# Rationale: on the log scale a multiplicative data-entry error (e.g. a
# misplaced decimal) has the same size upward and downward
valid_DBH[, `:=`(
  log_prev = shift(log_DBH, type = "lag"),
  log_next = shift(log_DBH, type = "lead")
), by = .(Tag, StemTag, TreeID, StemID)]

# Differences to the previous and to the next measurement, and between those two
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
# Step 3: Compute data-driven threshold (99th percentile of |d_prev|)
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

# Flags a measurement as an error if **ALL** of these are true for at least
# one of the three thresholds:
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
# Step 4a: Numeric error score (the larger of |d_prev| and |d_next|)
# --
valid_DBH[, error_score := fifelse(is.na(d_prev) & is.na(d_next), NA_real_, pmax(abs(d_prev), abs(d_next), na.rm = TRUE))]
# --
# Step 5: Remove intermediate columns
# --
valid_DBH[, c("log_DBH", "d_prev", "d_next", "d_span") := NULL]
# --
# Step 6: Stem-level summary (see 7.2)
# --
# stem_summary <- valid_DBH[
#   , .(
#     stem_has_any_error = any(entry_error_any, na.rm = TRUE),
#     stem_max_error_score = max(error_score, na.rm = TRUE)
#   ),
#   by = .(Tag, StemTag, TreeID, StemID)
# ]
# stem_summary[is.infinite(stem_max_error_score), stem_max_error_score := NA_real_]

# ---- 7.2 Stem-level summary & scoring ----
# Summarize error flags by stem for review and downstream processing.
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

# Combine both summaries
stem_summary <- rbind(single_summary, multi_summary)

setkey(valid_DBH, Tag, StemTag, TreeID, StemID)
valid_DBH[stem_summary, `:=`(
  stem_has_any_error = i.stem_has_any_error,
  stem_max_error_score = i.stem_max_error_score
)]
# --
# Step 7: Rows without a DBH get NA in the error columns
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

NA_DBH[, log_prev := NA_real_]
NA_DBH[, log_next := NA_real_]

ViewFullTable_measurement_error_indication <- rbindlist(list(valid_DBH, NA_DBH), use.names = TRUE)
setorder(ViewFullTable_measurement_error_indication, RowIDN1)

# how many errors detected?
table(ViewFullTable_measurement_error_indication$entry_error_any, useNA = "ifany")

# Only stems needing review
unique(ViewFullTable_measurement_error_indication[stem_has_any_error == TRUE, .(Tag, StemTag, TreeID, StemID)])

# Add candidate DBH values for likely entry errors using the geometric mean of adjacent valid measurements.
ViewFullTable_measurement_error_indication[, dbh_candidate := fifelse(
  entry_error_any & !is.na(log_prev) & !is.na(log_next),
  exp((log_prev + log_next) / 2),
  NA_real_
)]

# ---- 8. Prepare DBH candidates (for taper correction) ----
# Use the corrected candidate DBH for flagged entries, otherwise keep the original DBH.
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

# ---- 8.2 Apply taper correction ----
# Prepare `HOM_for_taper_correction` and run the taper correction to compute the
# DBH at 1.3 m (`dbh_with_best_candidate_taper_corrected`).
# Inputs: `dbh_with_best_candidate` (mm) and `HOM_for_taper_correction` (m).
ViewFullTable_measurement_error_indication[, HOM_for_taper_correction := HOM]

# # Load taper utilities (provides `apply_taper_correction()` and `taper()`)
# source(here("BCI_stem_reconstruction", "1_DATA_PREPARATION", "HELPER_FUNCTIONS", "taper_correction.R"))

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

# Load species table with growth-form classifications (object bci.spptable;
# path relative to the working directory)
load("./BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.RData")

# check database
setdiff(bci.spptable$Mnemonic, ViewFullTable_measurement_error_indication$Mnemonic)
setdiff(ViewFullTable_measurement_error_indication$Mnemonic, bci.spptable$Mnemonic)
# uniden is an unidentified species

# Growth forms (Spanish labels of the species list, lower case)
bci.spptable[, Lifeform := tolower(Lifeform_RPerez_SAguilar)]
growth_forms <- bci.spptable[, .(Mnemonic, Lifeform)]

# Species that get the taper-corrected DBH: those whose growth form contains
# "árbol" (tree, which includes "árbol estrangulador", strangler) or "arbusto"
# (shrub), and those without a growth form. Palms ("palma ...") and tree ferns
# ("helecho arbóreo") keep the uncorrected DBH.
species_to_use_tapper_corrected_dbh <- growth_forms[
  is.na(Lifeform) | grepl(pattern = "árbol|arbusto", x = Lifeform)
]$Mnemonic

# DBH is corrected for taper using Cushman et al. 2014.
# Taper adjusts DBH to what it would be at 1.3 m when measured higher (e.g., above
# buttresses). The corrected value is stored as
# `dbh_with_best_candidate_taper_corrected`.
#
# The equation is evaluated for every row; a DBH measured at 1.3 m is returned
# unchanged (the exponent hom - 1.3 is zero), so only rows with another HOM
# change. The corrected value is then kept for the species selected above.

# NOTE: dbh should be in cm for the equation.
# Sanity check
# Explore DBH units by comparing quantiles under different assumptions.
# The values confirm that DBH is recorded in millimetres.
check_diameter_units <- sort(ViewFullTable_measurement_error_indication[!is.na(dbh_with_best_candidate)]$dbh_with_best_candidate)

head(check_diameter_units, 100) # mm
tail(check_diameter_units, 100) ## check next lines

quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 10 # convert to cm from mm # GOOD
quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 10 / 100 # convert to m from mm # GOOD
quantile(check_diameter_units, probs = c(seq(0, 1, 0.25), 0.95, 0.99), na.rm = TRUE) / 100 # convert to m from cm #! WRONG

# taper_2014(): DBH at `common_hom` from a DBH measured at height `hom`
# (taper model of Cushman et al. 2014):
#   b      = exp(-2.0205 - 0.5053 * log(dbh_cm) + 0.3748 * log(hom))
#   dbh_at = dbh_cm / exp(-b * (hom - common_hom))
# dbh_mm : DBH in mm; hom : height of measurement in m (NA is read as
# common_hom); common_hom : target height in m.
# Returns the corrected DBH in mm; NA where the DBH or the HOM is not
# positive. Stops when dbh_mm and hom differ in length.
taper_2014 <- function(dbh_mm, hom, common_hom = 1.3) {
  if (length(dbh_mm) != length(hom)) {
    stop("'dbh_mm' and 'hom' must have the same length")
  }
  dbh_mm <- as.numeric(dbh_mm)
  hom <- as.numeric(hom)
  # Replace NA heights with 1.3 m (do not modify valid measured heights)
  hom_na <- is.na(hom)
  hom[hom_na] <- common_hom
  # convert dbh from mm to cm for the model
  dbh_cm <- dbh_mm / 10
  # Protect against log(0) or negative inputs by coercing non-positive values to NA
  dbh_cm[dbh_cm <= 0] <- NA_real_
  hom_for_log <- hom
  hom_for_log[hom_for_log <= 0] <- NA_real_
  b <- exp(-2.0205 - 0.5053 * log(dbh_cm) + 0.3748 * log(hom_for_log))
  out <- dbh_cm / (exp(-b * (hom - common_hom)))
  # convert back to mm and set invalid values to NA
  out_mm <- out * 10
  out_mm[is.na(out_mm) | is.infinite(out_mm)] <- NA_real_
  return(out_mm)
}

ViewFullTable_taper_corrected <- copy(ViewFullTable_measurement_error_indication)

# are there na homs with dbh candidates?
ViewFullTable_taper_corrected[!is.na(dbh_with_best_candidate) & is.na(HOM_for_taper_correction)]
# only 4; the function reads a missing HOM as 1.3 m

ViewFullTable_taper_corrected[
  ,
  dbh_with_best_candidate_taper_corrected_raw := taper_2014(
    dbh_mm = dbh_with_best_candidate, hom = HOM_for_taper_correction, common_hom = 1.3
  )
]

# check plot
with(
  ViewFullTable_taper_corrected[CensusID == 9],
  plot(dbh_with_best_candidate_taper_corrected_raw ~ dbh_with_best_candidate, pch = 16, cex = 0.5)
)
abline(a = 0, b = 1, col = "red")

# use the taper-corrected value for the species selected above, the uncorrected one for the others
ViewFullTable_taper_corrected[, dbh_with_best_candidate_taper_corrected := fifelse(
  Mnemonic %in% species_to_use_tapper_corrected_dbh,
  dbh_with_best_candidate_taper_corrected_raw,
  dbh_with_best_candidate
)][, dbh_with_best_candidate_taper_corrected_raw := NULL]

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

# ---- 10. Classify single- vs. multiple-stem tags ----
# A tag is single-stemmed if all its StemTags are NA and it has exactly one
# unique StemID across all censuses. All other tags go to the DP algorithm.
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

# get percentage of tags with one stem or multiple stems
tags_info <- round(table(id_single_stem_tags$all_stemtag_na & id_single_stem_tags$one_stemid) / nrow(id_single_stem_tags) * 100, 1)
# FALSE  TRUE
#  28.7  71.3

# select tags with all StemTag NA but only one unique StemID
# These have been already identified as single individual
tags_with_one_stemid_no_stemtag <- id_single_stem_tags[all_stemtag_na == TRUE & one_stemid == TRUE, Tag]
length(tags_with_one_stemid_no_stemtag)

# Check the single-stem tags with missing StemID or DBH values.
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]

# Confirm these tags use at most one StemID across censuses.
unique(ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag, .N, by = .(Tag, StemID, CensusID)]$N)

# Inspect any DBH values without StemID and any NA-DBH rows.
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & !is.na(DBH) & is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]
ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag & is.na(DBH) & !is.na(StemID), .(Tag, CensusID, StemTag, StemID, DBH)]

# Quickly inspect the key variables
inspectdf::inspect_na(ViewFullTable_taper_corrected[Tag %in% tags_with_one_stemid_no_stemtag, .(Tag, StemTag, TreeID, StemID, CensusID, DBH)])

## ---- Mark multi-stemmed tags that need to run DP ------------------------------------------------------

single_stem_tags <- id_single_stem_tags[Tag %in% tags_with_one_stemid_no_stemtag]$Tag
multiple_stem_tags <- id_single_stem_tags[!Tag %in% tags_with_one_stemid_no_stemtag]$Tag

# ! SANITY CHECKS
# Are the counts (single + multiple) correct?
(length(single_stem_tags) + length(multiple_stem_tags)) ==
  length(unique(ViewFullTable_taper_corrected$Tag)) # should be TRUE

# Do the rows of the two groups add up to the whole table?
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

# ---------------------------------------------------------------------------
# Add the growth form (`Lifeform`) of each species to the observation table.
# The summaries below count the records per growth form and show the
# mnemonics without one.
ViewFullTable_single_vs_multiple_stem_tags <- merge(ViewFullTable_single_vs_multiple_stem_tags, unique(growth_forms[, .(Mnemonic, Lifeform)]), by = "Mnemonic", all.x = TRUE)
setorder(ViewFullTable_single_vs_multiple_stem_tags, RowID)

# Summaries
ViewFullTable_single_vs_multiple_stem_tags[, .N, by = Lifeform][order(-N)]

unique(ViewFullTable_single_vs_multiple_stem_tags[is.na(Lifeform), .(Mnemonic, Genus, SpeciesName)])

tags_per_growth_form <- unique(ViewFullTable_single_vs_multiple_stem_tags[, .(Tag, single_stem_tags, Lifeform)])
# Check individuals per growth form by single vs. multiple stem tag
tags_per_growth_form[
  , .N,
  by = .(single_stem_tags, Lifeform)
][order(single_stem_tags, -N)]

# Save the stage-1 output: the prepared observation table with growth forms.
saveRDS(ViewFullTable_single_vs_multiple_stem_tags, file.path(OUTPUT_folder, "ViewFullTable_single_vs_multiple_stem_tags.rds"))

# Measured rows per tag, then tags with a measurement per growth form
nobs_growth_form <- ViewFullTable_single_vs_multiple_stem_tags[!is.na(DBH)][
  , .N,
  by = .(Tag, SpeciesName, Lifeform)
]

nobs_growth_form[
  , .N,
  by = Lifeform
][order(-N)]

# the only observation with NA lifeform is an unidentified species.
unique(ViewFullTable_single_vs_multiple_stem_tags[is.na(Lifeform), .(Mnemonic, Genus, SpeciesName)])

# ---- 11. Check all Tag × CensusID combinations are complete ----
# Final completeness check on the taper-corrected, growth-form-enriched table.
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
