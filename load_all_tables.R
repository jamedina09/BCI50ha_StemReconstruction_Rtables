rm(list = ls())

library(data.table)

# for (i in 1:9) {
#     load(paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", i, ".Rdata"))
#     ## Combine all the data.tables into one
#     if (i == 1) {
#         bci.stem <- get(paste0("bci.stem", i))
#     } else {
#         bci.stem <- rbind(bci.stem, get(paste0("bci.stem", i)))
#     }
#     ## clean up the environment
#     rm(list = paste0("bci.stem", i))
# }

# bci.stem <- as.data.table(bci.stem)


rol <- data.table(readxl::read_excel("./BCI_stem_reconstruction/DATA/RAW/sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx"))
rol <- janitor::clean_names(rol)
o13 <- fread("./BCI_stem_reconstruction/DATA/RAW/oct_13_2026/ViewTaxonomy_bci.csv")
o16 <- fread("./BCI_stem_reconstruction/DATA/RAW/oct_16_2026/ViewFullTable_bci.csv")

inc <- c("swars1", "swars2", "protte", "quaras", "appuse", "pipeco", "sapisp")

rol[codigo %in% inc, .(codigo, genero, especie)][order(codigo)]

o13[Mnemonic %in% inc, .(Mnemonic, Genus, SpeciesName, SpeciesName, Subspecies, SubspeciesID)][order(Mnemonic)]

unique(o16[Mnemonic %in% inc, .(Mnemonic, Genus, SpeciesName, SpeciesName, Subspecies, SubspeciesID)][order(Mnemonic)])
