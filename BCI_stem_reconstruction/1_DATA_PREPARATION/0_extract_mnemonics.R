# Clear all objects from the workspace to avoid accidental contamination
rm(list = ls())

library(data.table)


site <- "bci" # Site code for BCI

# Input folder (raw ViewFullTable export)
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
  workspace_root <- dirname(workspace_root)
}

INPUT_folder_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses")

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

mnemonics_raw <- data.table(sort(unique(ViewFullTable$Mnemonic)))
# rename V1 to Mnemonic
setnames(mnemonics_raw, "V1", "Mnemonic")

saveRDS(mnemonics_raw, file = file.path(workspace_root, "BCI_stem_reconstruction", "1_DATA_PREPARATION", "mnemonics_raw.rds"))
