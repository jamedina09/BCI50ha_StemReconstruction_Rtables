rm(list = ls())

workspace_root <- getwd()

library(data.table)

bci_stem_nums <- as.character(1:9)

census_list <- lapply(bci_stem_nums, function(num) {
    filepath <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RTABLES", paste0("bci.stem", num, ".Rdata"))
    if (!file.exists(filepath)) stop("Missing census file: ", filepath)
    load(filepath)
    get(paste0("bci.stem", num))
})
names(census_list) <- paste0("bci.stem", bci_stem_nums)

df_stem <- data.table::rbindlist(census_list, fill = TRUE, idcol = "censusID")
rm(census_list, bci_stem_nums)

# lapply nrow
lapply(census_list, nrow)
