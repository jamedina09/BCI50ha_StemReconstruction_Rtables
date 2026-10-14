# =============================================================================
# 0_extract_ViewTaxonomy.R
#
# Purpose: Extract the database taxonomy (one row per species code) from the
#          current ViewFullTable export, for 1_prepare_species_tables.R.
#
# Run from the project root (or from BCI_stem_reconstruction/).
#
# Input:   BCI_stem_reconstruction/DATA/RAW/oct_13_2026/ViewFullTable_bci.csv
# Outputs: BCI_stem_reconstruction/DATA/RAW/oct_13_2026/ViewTaxonomy_bci.csv
#          and ViewTaxonomy_bci.rds, with Mnemonic, Family, Genus,
#          SpeciesName, Subspecies, SpeciesID, SubspeciesID.
# =============================================================================

# Clear all objects from the workspace to avoid accidental contamination
rm(list = ls())

library(data.table)

site <- "bci" # Site code for BCI

# Input folder (raw ViewFullTable export)
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
  workspace_root <- dirname(workspace_root)
}

INPUT_folder_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026")

file1 <- file.path(INPUT_folder_1, paste0("ViewFullTable_", site, ".csv"))

# Loud check (the script is run line by line): prints ❌, raises an immediate
# warning, writes to the log, and only then stops.
log_file <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS", "extract_ViewTaxonomy_log.txt")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
check <- function(ok, msg) {
  if (isTRUE(all(ok))) {
    cat("✓", msg, "\n")
    return(invisible(TRUE))
  }
  full_msg <- paste0("CHECK FAILED: ", msg)
  cat("❌", full_msg, "\n")
  warning(full_msg, call. = FALSE, immediate. = TRUE)
  cat(format(Sys.time()), full_msg, "\n", file = log_file, append = TRUE)
  stop(full_msg, call. = FALSE)
}

# Only the taxonomy columns are read
taxonomy_cols <- c("Mnemonic", "Family", "Genus", "SpeciesName", "Subspecies", "SpeciesID", "SubspeciesID")
ViewTaxonomy_raw <- fread(
  file = file1,
  sep = "\t",
  select = taxonomy_cols,
  na.strings = c("NA", "NULL", "")
)

# One row per distinct combination of names and IDs
ViewTaxonomy <- unique(ViewTaxonomy_raw)
setorder(ViewTaxonomy, Mnemonic)

# Codes with more than one row: protte and quaras each have one row with a
# SubspeciesID but no subspecies name, and one with neither. The row with the
# SubspeciesID but no name is dropped (a rule, not row numbers, so it does not
# depend on the row order of the export).
dup_codes <- ViewTaxonomy[duplicated(Mnemonic), unique(Mnemonic)]
ViewTaxonomy[Mnemonic %in% dup_codes]

ViewTaxonomy <- ViewTaxonomy[!(
  Mnemonic %in% dup_codes & is.na(Subspecies) & !is.na(SubspeciesID)
)]

check(!anyDuplicated(ViewTaxonomy$Mnemonic), "one row per code in the database taxonomy")
check(
  setequal(ViewTaxonomy$Mnemonic, unique(ViewTaxonomy_raw$Mnemonic)),
  "every code of the export is kept"
)

fwrite(
  ViewTaxonomy,
  file = file.path(INPUT_folder_1, paste0("ViewTaxonomy_", site, ".csv"))
)

saveRDS(
  ViewTaxonomy,
  file = file.path(INPUT_folder_1, paste0("ViewTaxonomy_", site, ".rds"))
)

# Rows per code (one each)
ViewTaxonomy[, .N, by = Mnemonic][order(-N)]
