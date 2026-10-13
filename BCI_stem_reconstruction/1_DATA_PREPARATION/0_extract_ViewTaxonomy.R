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

ViewTaxonomy_raw <- ViewFullTable[, .(Mnemonic, Family, Genus, SpeciesName, Subspecies, SpeciesID, SubspeciesID)]
ViewTaxonomy <- unique(copy(ViewTaxonomy_raw))
setorder(ViewTaxonomy, Mnemonic)

## add row id
ViewTaxonomy[, RowID := .I]

inc <- ViewTaxonomy[is.na(Subspecies) & !is.na(SubspeciesID)]$Mnemonic

ViewTaxonomy[Mnemonic %in% inc]

ViewTaxonomy <- ViewTaxonomy[!RowID %in% c(236, 257)]

length(unique(ViewTaxonomy$Mnemonic)) == length(unique(ViewTaxonomy_raw$Mnemonic))

ViewTaxonomy[, RowID := NULL]

fwrite(
  ViewTaxonomy,
  file = file.path(INPUT_folder_1, paste0("ViewTaxonomy_", site, ".csv"))
)

saveRDS(
  ViewTaxonomy,
  file = file.path(INPUT_folder_1, paste0("ViewTaxonomy_", site, ".rds"))
)

# count N for Mnemonic in ViewTaxonomy
ViewTaxonomy[, .N, by = Mnemonic][order(-N)]

