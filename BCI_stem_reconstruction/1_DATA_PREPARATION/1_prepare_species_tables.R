# =============================================================================
# 1_prepare_species_tables.R
#
# Purpose: Prepare the BCI 50-ha plot species table, one row per species code
#          (mnemonic) of the plot. The database taxonomy is corrected with the
#          emails of R. Pérez and S. Lao and completed from the Flora de Panamá
#          list, which has priority. The names are then checked, step by step,
#          against the current World Flora Online (WFO) backbone (a local copy
#          of the June 2026 release) with the WorldFlora package. WFO only
#          suggests corrections, mainly of spelling and of the format of the
#          authority; each one is adopted by hand, nothing is changed
#          automatically.
#
# Run from the project root (or from BCI_stem_reconstruction/), line by line.
#
# Sections:
#    0. Setup: packages, check helper, input paths
#    1. Load the data
#    2. Compare the old and new database taxonomy (information only)
#    3. Fix 1: species list of R. Pérez and S. Aguilar (email 1)
#    4. Fix 2: name changes of S. Lao and the reply of R. Pérez (email 2)
#    5. Fix 3: Flora de Panamá list (has priority)
#    6. Codes replaced in the inventory; growth forms; English column names
#    7. Notes on known data problems
#    8. Check the names against the current WFO backbone
#    9. Review the WFO suggestions (step by step)
#   10. Complete and check the remaining fields
#   11. Assemble the species table and export
#
# Inputs (BCI_stem_reconstruction/DATA/RAW/):
#   oct_13_2026/ViewTaxonomy_bci.rds               database taxonomy (from
#                                                  0_extract_ViewTaxonomy.R)
#   oct_13_2026/ViewFullTable_bci.csv              current inventory (codes in the plot)
#   ViewFiles_bci_allcensuses/ViewFullTable_bci.csv  previous inventory (section 2 only)
#   sp_tables/email_1/Lista_bci_mnemonics_formadevida.xlsx  list of R. Pérez & S. Aguilar
#   sp_tables/email_2/bci_name_changes.xlsx        name changes of S. Lao
#   sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx
#   sp_tables/wfo_plant_list_2026_06/classification.csv   WFO backbone, June 2026
#                                                  (Zenodo doi:10.5281/zenodo.20782718)
# Outputs:
#   DATA/SPP_TABLE/bci_spptable.txt, bci_spptable.csv (both tab-separated) and
#   bci_spptable.RData (object bci.spptable): one row per code of the plot
#   except uniden, plus ficuau (kept without records, see 6a), with Mnemonic, Order, Family, Genus, SpeciesName,
#   InfraspecificRank, InfraspecificEpithet, Authority, Synonyms,
#   Lifeform_RFoster, Lifeform_RPerez_SAguilar, CommonName, Herbarium, Notes.
#   DATA/CHECKS/species_wfo_review.csv (each name next to its WFO record),
#   species_wfo_candidates.csv (names with several WFO records),
#   species_wfo_accepted_changes.csv (changes adopted in the review),
#   prepare_species_tables_log.txt (failed checks).
# Packages: data.table, WorldFlora, stringr, compareDF, janitor, readxl.
# =============================================================================

# Clear all objects from the workspace to avoid accidental contamination
rm(list = ls())

# data.table print options
options(
    datatable.print.class = FALSE,
    datatable.print.keys = TRUE,
    datatable.verbose = FALSE
)

# =============================================================================
# 0. SETUP
# =============================================================================

# --- 0a. Packages ------------------------------------------------------------
library(data.table) # Fast table manipulation
library(WorldFlora) # Matching plant names against the WFO backbone
library(stringr) # String helpers (squish, case)
library(compareDF) # Comparison of the old and new database taxonomy (section 2)
library(janitor) # clean_names() for the spreadsheet column names

# Number of threads data.table will use
data.table::getDTthreads()

# --- 0b. Paths ---------------------------------------------------------------
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
    workspace_root <- dirname(workspace_root)
}

# Database taxonomy: one row per code, extracted from the current inventory by
# 0_extract_ViewTaxonomy.R
ViewTaxonomy <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026", "ViewTaxonomy_bci.rds")
# Species list of R. Pérez and S. Aguilar (email 1)
TAXONOMY_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_1", "Lista_bci_mnemonics_formadevida.xlsx")
# Name changes of S. Lao (email 2)
TAXONOMY_2 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_2", "bci_name_changes.xlsx")
# Flora de Panamá list (R. Pérez, received 2026-10-09)
flora_panama_oct_9_2026 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "flora_panama_oct_9_2026", "Lista de especies de la Flora de Panamá_nov25 original.xlsx")
# Previous inventory export, used only for the comparison of section 2
ViewFullTable_old <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses", "ViewFullTable_bci.csv")
# Current inventory (the table 0_extract_ViewTaxonomy.R reads)
ViewFullTable_new <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026", "ViewFullTable_bci.csv")
# Current WFO backbone: World Flora Online Plant List June 2026 (v.2026-06,
# CC0), file _DwC_backbone_R.zip of doi:10.5281/zenodo.20782718, unzipped
WFO_BACKBONE <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "wfo_plant_list_2026_06", "classification.csv")
WFO_BACKBONE_VERSION <- "2026-06"

# --- 0c. Check helper ----------------------------------------------------------
# Loud check (the script is run line by line): prints ❌, raises an immediate
# warning, writes to the log, and only then stops.
log_file <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS", "prepare_species_tables_log.txt")
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

# Empty values in the raw tables: NA, "" or the text "NULL"
is_empty <- function(x) is.na(x) | x %in% c("", "NULL")

# =============================================================================
# 1. LOAD THE DATA
# =============================================================================
viewt <- as.data.table(readRDS(ViewTaxonomy))
tax1 <- as.data.table(readxl::read_xlsx(TAXONOMY_1))
tax2 <- as.data.table(readxl::read_xlsx(TAXONOMY_2))
flor_pan <- as.data.table(readxl::read_xlsx(flora_panama_oct_9_2026))

# Previous inventory: only the taxonomy columns of the database taxonomy
inc <- colnames(viewt)
viewftold <- fread(ViewFullTable_old, select = inc)

viewtold <- unique(viewftold[, ..inc])
setorder(viewtold, Mnemonic)

# =============================================================================
# 2. COMPARE THE OLD AND NEW DATABASE TAXONOMY (INFORMATION ONLY)
# =============================================================================
# Names of every code in the previous inventory export against those of the
# current one (the database renamed 49 species in October 2026). Nothing here
# changes the species table.

# One row per code in the current taxonomy
viewt[, .N, by = Mnemonic][order(-N)]

# The previous export has two rows for protte and quaras: one with an empty
# subspecies name but a SubspeciesID, one without either. The first is dropped
# (the same rule 0_extract_ViewTaxonomy.R applies to the current export).
viewtold[, .N, by = Mnemonic][order(-N)]
viewtold[Mnemonic %in% viewtold[duplicated(Mnemonic), Mnemonic]]
viewtold <- viewtold[!(
    Mnemonic %in% viewtold[duplicated(Mnemonic), Mnemonic] &
        is_empty(Subspecies) & !is_empty(SubspeciesID)
)]
check(!anyDuplicated(viewtold$Mnemonic), "one row per code in the previous database taxonomy")

viewt <- viewt[, .(Mnemonic, Family, Genus, SpeciesName, Subspecies)]
viewtold <- viewtold[, .(Mnemonic, Family, Genus, SpeciesName, Subspecies)]

# The text "NULL" and empty strings mean no value
viewt[viewt == "NULL"] <- NA
viewtold[viewtold == "NULL"] <- NA
viewt[viewt == ""] <- NA
viewtold[viewtold == ""] <- NA

check(identical(viewt$Mnemonic, viewtold$Mnemonic), "the old and new database taxonomy have the same codes, in the same order")

# Codes whose names changed ("+" new, "-" old)
ctable_mnemonic <- compare_df(viewt, viewtold, c("Mnemonic"))
ctable_mnemonic$comparison_df

# To view the comparison as a table, or save it as a spreadsheet:
# create_output_table(ctable_mnemonic)
# create_output_table(ctable_mnemonic, output_type = "xlsx", file_name = "./BCI_stem_reconstruction/1_DATA_PREPARATION/compare_2025_vs_oct2026_taxonomy_from_viewfulltables.xlsx")

rm(viewtold, viewftold)

# =============================================================================
# 3. FIX 1: SPECIES LIST OF R. PÉREZ AND S. AGUILAR (EMAIL 1)
# =============================================================================
# Email from Rolando Pérez (sent May 14, 2026, forwarded Sep 28, 2026), with
# the species list and known corrections, kept here for reference:

# Monday, September 28, 2026 at 9:07:43
# Subject: Lista de las 50 has
# Date: Thursday, May 14, 2026 at 09:25:16
# From: Perez, Rolando
# To: Lao, Suzanne, Medina-Vega, Jose A.
# CC: Aguilar, Salomon, Mitre, David S., Davies, Stuart J.
# Attachments: Lista_bci_mnemonics_formadevida.xlsx
# Buenos días a todos,
# Esta es la lista de especies actualizada y elaborada por R. Pérez y S. Aguilar registradas en todos
# los censos de la parcela de 50 hectáreas de BCI.
# Quiero aclarar algunas correcciones que realizamos:
# El nombre Appunia siebertiii es correcto y actualizado, Morinda siebertii es un sinónimo.
# Apeiba "hybrida" murió.
# Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una
# con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla.
# Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
# Beilschmiedia pendula y Quararibea asterolepis son especies que no existen para BCI. Las
# especies correctas son Beilschmiedia tovarensis y Quararibea stenophylla, ambas fueron
# confirmadas por los especialistas en nuestras muestras de herbario.
# Añadimos las formas de vida que utilizo Robin Foster con pequeñas modificaciones realizadas por
# R. Pérez y S. Aguilar. También está la forma de vida que Salomón y Yo le hemos asignado a las
# especies de la parcela en español.
# Los nombres comunes de las especies y nuestras colecciones de herbario.
# Tengan en cuenta esto para las publicaciones
# Esperamos sea de utilidad.
# Saludos
# Rolando

# The database has 329 codes and the list 327: codes not in the list
setdiff(viewt$Mnemonic, tax1$codigo)

# Database taxonomy next to the list (all codes of both)
viewt_tax1 <- merge(
    viewt,
    tax1,
    by.x = "Mnemonic",
    by.y = "codigo",
    all.x = TRUE,
    all.y = TRUE
)

# The corrections of the email, one by one, recorded in `email_1_comment`

# 1) El nombre Appunia siebertiii es correcto y actualizado, Morinda siebertii es un sinónimo.
viewt_tax1[genero == "Appunia", email_1_comment := "Appunia siebertii es correcto; Morinda siebertii es un sinónimo"]

# 2) Apeiba "hybrida" murió.
viewt_tax1[Mnemonic == "apeihy", email_1_comment := "Apeiba hybrida murió"]

# 3) Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una
# con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla.
viewt_tax1[Mnemonic == "nects1", email_1_comment := "Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla."]
viewt_tax1[Mnemonic == "nects3", email_1_comment := "Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla."]

# 4) Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
viewt_tax1[Genus == "Pterocarpus"]
viewt_tax1[Mnemonic == "pterof", email_1_comment := "Pterocarpus officinalis no existe en la parcela; corresponde a P. rohrii"]
viewt_tax1[Mnemonic == "pterro", email_1_comment := "Pterocarpus officinalis no existe en la parcela; corresponde a P. rohrii"]

# 5) Beilschmiedia pendula y Quararibea asterolepis son especies que no existen para BCI. Las
# especies correctas son Beilschmiedia tovarensis y Quararibea stenophylla, ambas fueron
# confirmadas por los especialistas en nuestras muestras de herbario.
viewt_tax1[Genus == "Beilschmiedia"]
viewt_tax1[Genus == "Beilschmiedia" & SpeciesName == "tovarensis", email_1_comment := "Beilschmiedia pendula no existe para BCI; la especie correcta es Beilschmiedia tovarensis"]
viewt_tax1[Genus == "Quararibea"]
viewt_tax1[Genus == "Quararibea" & SpeciesName == "stenophylla", email_1_comment := "Quararibea asterolepis no existe para BCI; la especie correcta es Quararibea stenophylla"]

viewt_tax1[!is.na(email_1_comment), .(Mnemonic, Family, Genus, SpeciesName, email_1_comment)]

# The growth forms (R. Foster; R. Pérez & S. Aguilar), common names and
# herbarium collections come with the list (its columns).

# Family, genus and species: database vs list. The list is kept; where it has
# no value (the codes it does not include), the database value is used.
viewt_tax1[(is.na(familia) & !is.na(Family)) | Family != familia, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(familia), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(familia), familia := Family]

viewt_tax1[(is.na(genero) & !is.na(Genus)) | Genus != genero, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(genero), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(genero), genero := Genus]

viewt_tax1[(is.na(especie) & !is.na(SpeciesName)) | SpeciesName != especie, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(especie), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(especie), especie := SpeciesName]

# Codes with a subspecies in the database (kept in `Subspecies` until section 5)
viewt_tax1[!is.na(Subspecies)]

# The database names are no longer needed
viewt_tax1[, Family := NULL]
viewt_tax1[, Genus := NULL]
viewt_tax1[, SpeciesName := NULL]

# Done: email 1 included

# =============================================================================
# 4. FIX 2: NAME CHANGES OF S. LAO AND THE REPLY OF R. PÉREZ (EMAIL 2)
# =============================================================================
# Email from Suzanne Lao (July 28, 2026) with the name changes, and the reply
# of Rolando Pérez with two more corrections, kept here for reference.

# First email (Suzanne Lao)

# From: Lao, Suzanne <LAOZ@si.edu>
# Sent: Tuesday, July 28, 2026 2:28 PM
# To: Medina-Vega, Jose A. <MedinaJA@si.edu>
# Cc: Davies, Stuart J. <DaviesS@si.edu>; Perez, Rolando <PerezRo@si.edu>
# Subject: BCI species name changes

# Hola José,

# Adjunto le envío el listado de los nombres de especies de la parcela de BCI. Las primeras
# columnas son los nombres nuevos y las columnas al final son los nombres viejos.

# Solamente un código cambió: todos los "ficutr" cambiaron a "ficuar".
# NOTE: "ficuar" is a typo for ficuau (Ficus aurea); there is no code ficuar.
# Me avisa si este listado es diferente al que le envió Rolando. Es posible que el listado de
# Rolando esté más actualizado.

# Gracias y saludos,
# Suzanne

# The table has one row per code; Suzanne's list has 330 rows, with two for
# protte
viewt_tax1[, .N, by = Mnemonic][order(-N)]
tax2[, .N, by = mnemonic][order(-N)]

# protte: one row with the SubspeciesID 15 but no subspecies name, one with
# neither; the first is dropped (the same rule as in section 2)
tax2[mnemonic %in% tax2[duplicated(mnemonic), mnemonic]]
tax2 <- tax2[!(
    mnemonic %in% tax2[duplicated(mnemonic), mnemonic] &
        is_empty(subspecies) & !is_empty(subspeciesid)
)]
check(!anyDuplicated(tax2$mnemonic), "one row per code in the list of S. Lao")

setorder(viewt_tax1, Mnemonic)
setorder(tax2, mnemonic)

# Codes in only one of the two: ficutr in the table, ficuau in Suzanne's list
setdiff(viewt_tax1$Mnemonic, tax2$mnemonic)
setdiff(tax2$mnemonic, viewt_tax1$Mnemonic)

viewt_tax1[Mnemonic == "ficutr"]
tax2[mnemonic == "ficuau"]

# Species and subspecies IDs of the list are not used
tax2[, speciesid := NULL]
tax2[, subspeciesid := NULL]
tax2[, speciesidold := NULL]
tax2[, subspeciesidold := NULL]

# Second email (Rolando Pérez, reply to Suzanne Lao)

# Subject: RE: BCI species name changes
# Date: Wednesday, August 5, 2026 at 06:07:09 Eastern Standard Time
# From: Perez, Rolando
# To: Lao, Suzanne, Medina-Vega, Jose A.
# CC: Davies, Stuart J

# Buenos días a todos,

# En el listado hay que actualizar dos nombres: (stemgr ahora es Tabernaemontana
# grandiflora Jacq.) y (pla1pi ahora es Platymiscium pinnatum (Jacq.) Dugand).
# Los nombres que aparecen en la lista no son sinónimos. Los sinónimos de ambas
# especies son: Stemmadenia grandiflora y Platymiscium polystachyum.

# Saludos
# Rolando

# The two corrections of Rolando: the table already has both names (from the
# list of email 1); only the synonym of pla1pi is missing
tax2[mnemonic %in% c("stemgr", "pla1pi")]
viewt_tax1[Mnemonic %in% c("stemgr", "pla1pi")]

viewt_tax1[Mnemonic %in% c("pla1pi"), sinonimos := "Platymiscium polystachyum"]
viewt_tax1[Mnemonic %in% c("stemgr", "pla1pi"), .(Mnemonic, genero, especie, autoridad, sinonimos)]

# Family, genus, species and subspecies: the table vs Suzanne's list. The
# table (the list of email 1) is kept: Suzanne says that Rolando's list may be
# more up to date. The rows printed are the differences, for information.
mefam_susfam <- merge(
    viewt_tax1[, .(Mnemonic, familia)],
    tax2[, .(mnemonic, family, familyold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
mefam_susfam[familia != family]

megen_susgen <- merge(
    viewt_tax1[, .(Mnemonic, genero)],
    tax2[, .(mnemonic, genus, genusold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
megen_susgen[genero != genus]

mesp_sussp <- merge(
    viewt_tax1[, .(Mnemonic, especie)],
    tax2[, .(mnemonic, speciesname, speciesnameold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
mesp_sussp[especie != speciesname]

mesubsp_sussubsp <- merge(
    viewt_tax1[, .(Mnemonic, Subspecies)],
    tax2[, .(mnemonic, subspecies, subspeciesold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
mesubsp_sussubsp[Subspecies != subspecies] # no differences

# ficutr and ficuau. In the database, ficutr carries the name Ficus aurea
# (there is no Ficus trigonata in the plot), while the list of email 1 gives
# ficutr as Ficus trigonata. Suzanne's change ficutr -> ficuau puts these trees
# under the code of Ficus aurea. The change is not decided yet, so both codes
# are kept: a row for ficuau is added here from Suzanne's list, next to the
# ficutr row (see 6a). The recode of the inventory is a switch in
# 2_prepare_viewfulltable.R.
viewt[Mnemonic == "ficutr"]
viewt[Genus == "Ficus"][order(SpeciesName)]

viewt_tax1_updated <- rbindlist(
    list(
        viewt_tax1,
        data.table(
            Mnemonic = "ficuau",
            familia = tax2[mnemonic == "ficuau"]$family,
            genero = tax2[mnemonic == "ficuau"]$genus,
            especie = tax2[mnemonic == "ficuau"]$speciesname,
            email_2_comment = "Solamente un codigo cambio: todos los 'ficutr' cambiaron a 'ficuau'"
        )
    ),
    use.names = TRUE,
    fill = TRUE
)

# Codes of the table not in the database taxonomy (ficuau only)
setdiff(viewt_tax1_updated$Mnemonic, viewt$Mnemonic)

# Done: email 1, email 2 and Rolando's reply included

# =============================================================================
# 5. FIX 3: FLORA DE PANAMÁ LIST (HAS PRIORITY)
# =============================================================================
# Email from Rolando Pérez (Oct 9, 2026) with the Flora de Panamá list, the
# most up-to-date list (its spreadsheet was last modified on 2026-10-07). Where
# it has a value, that value is kept; the other sources only fill its gaps.

# Hola José,

# Creo que esta lista te puede aclarar algunas interrogantes, incluye árboles, arbustos,
# palmas y lianas de la Flora de Panamá. Es lo más actualizado que tenemos hasta el
# momento. MNC indica la presencia en el Monumento de Barro Colorado, todavía sigo
# trabajando en esto.

# Tiene los códigos de 6 letras que Susy hizo, algunos son registros recientes y
# hay que agregarles un código, pocos casos, si tienes alguna pregunta con gusto
# hablamos.

# Incluye nombres comunes, sinónimos, forma de vida (todavía por trabajar más), si
# tenemos fotos, testigos de herbario y descripciones. Estamos trabajando en esto para
# antes de salir de STRI montar todo en la página web de ForestGEO, esa es nuestra meta,
# tener una flora digital de plantas leñosas de Panamá. Es como dejar nuestra herencia
# botánica para las futuras generaciones.

# Saludos y espero te sea útil.

# Rolando

viewt_tax1_updated

# Column names in lower case without spaces or accents; the Flora de Panamá
# columns get the prefix "pan_" so both tables can be merged
flor_pan <- clean_names(flor_pan)
viewt_tax1_updated <- clean_names(viewt_tax1_updated)
colnames(flor_pan) <- paste0("pan_", colnames(flor_pan))

# The text "NULL" means no value
viewt_tax1_updated <- viewt_tax1_updated[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]
flor_pan <- flor_pan[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]

# Each code must appear only once in the Flora de Panamá list, otherwise the
# merge would duplicate rows
check(
    !anyDuplicated(flor_pan[pan_codigo %in% viewt_tax1_updated$mnemonic, pan_codigo]),
    "each code of the table appears at most once in the Flora de Panamá list"
)

viewt_tax1_updated_flor_pan <- merge(
    viewt_tax1_updated,
    flor_pan,
    by.x = "mnemonic",
    by.y = "pan_codigo",
    all.x = TRUE,
    all.y = FALSE
)

# Columns that are not used: photos and texts (both lists), the presence in
# the Barro Colorado Nature Monument (MNC), and the growth form of the Flora de
# Panamá list (incomplete; the growth forms of the list of email 1 are used,
# see 6b)
viewt_tax1_updated_flor_pan[, texto := NULL][, fotos := NULL][, pan_texto := NULL][, pan_fotos := NULL][, pan_mnbc := NULL]
viewt_tax1_updated_flor_pan[, pan_forma_de_vida := NULL]

# Codes of the table vs the database taxonomy: only ficuau is new (it has no
# records in the inventory unless ficutr is recoded, see 6a)
setdiff(viewt_tax1_updated_flor_pan$mnemonic, viewt$Mnemonic)
setdiff(viewt$Mnemonic, viewt_tax1_updated_flor_pan$mnemonic)

# Each field: the Flora de Panamá value where it has one, otherwise the value
# of the table (`<field>_corrected`). The rows printed before and after are
# those where the two sources differ or one is missing.

# orden
viewt_tax1_updated_flor_pan[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden)]
viewt_tax1_updated_flor_pan[, orden_corrected := fcoalesce(pan_orden, orden)]
viewt_tax1_updated_flor_pan[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden, orden_corrected)]
viewt_tax1_updated_flor_pan[, orden := NULL][, pan_orden := NULL]

# familia
viewt_tax1_updated_flor_pan[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia)]
viewt_tax1_updated_flor_pan[, familia_corrected := fcoalesce(pan_familia, familia)]
viewt_tax1_updated_flor_pan[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia, familia_corrected)]
viewt_tax1_updated_flor_pan[, familia := NULL][, pan_familia := NULL]

# genero
viewt_tax1_updated_flor_pan[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero)]
viewt_tax1_updated_flor_pan[, genero_corrected := fcoalesce(pan_genero, genero)]
viewt_tax1_updated_flor_pan[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero, genero_corrected)]
viewt_tax1_updated_flor_pan[, genero := NULL][, pan_genero := NULL]

# especie (for the two varieties of Swartzia simplex the Flora de Panamá
# species includes the variety, e.g. "simplex var. grandiflora"; it is split in
# section 8)
viewt_tax1_updated_flor_pan[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie)]
viewt_tax1_updated_flor_pan[, especie_corrected := fcoalesce(pan_especie, especie)]
viewt_tax1_updated_flor_pan[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie, especie_corrected)]
viewt_tax1_updated_flor_pan[, especie := NULL][, pan_especie := NULL]

# The database subspecies is now part of especie_corrected
viewt_tax1_updated_flor_pan[!is.na(subspecies)]
viewt_tax1_updated_flor_pan[, subspecies := NULL]

# autoridad
viewt_tax1_updated_flor_pan[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad)]
viewt_tax1_updated_flor_pan[, autoridad_corrected := fcoalesce(pan_autoridad, autoridad)]
viewt_tax1_updated_flor_pan[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad, autoridad_corrected)]
viewt_tax1_updated_flor_pan[, autoridad := NULL][, pan_autoridad := NULL]

# sinonimos
viewt_tax1_updated_flor_pan[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos)]
viewt_tax1_updated_flor_pan[, sinonimos_corrected := fcoalesce(pan_sinonimos, sinonimos)]
viewt_tax1_updated_flor_pan[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos, sinonimos_corrected)]
viewt_tax1_updated_flor_pan[, sinonimos := NULL][, pan_sinonimos := NULL]

# nombre comun (the column is spelled "nombre_comum" in both lists)
viewt_tax1_updated_flor_pan[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum)]
viewt_tax1_updated_flor_pan[, nombre_comun_corrected := fcoalesce(pan_nombre_comum, nombre_comum)]
viewt_tax1_updated_flor_pan[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum, nombre_comun_corrected)]
viewt_tax1_updated_flor_pan[, nombre_comum := NULL][, pan_nombre_comum := NULL]

# herbario
viewt_tax1_updated_flor_pan[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario)]
viewt_tax1_updated_flor_pan[, herbario_corrected := fcoalesce(pan_herbario, herbario)]
viewt_tax1_updated_flor_pan[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario, herbario_corrected)]
viewt_tax1_updated_flor_pan[, herbario := NULL][, pan_herbario := NULL]

# Spot checks: rows discussed in the emails
viewt_tax1_updated_flor_pan[mnemonic == "appuse"]
viewt_tax1_updated_flor_pan[mnemonic == "ficutr"]
viewt_tax1_updated_flor_pan[mnemonic == "ficuau"]
viewt_tax1_updated_flor_pan[!is.na(email_1_comment)]

# =============================================================================
# 6. CODES REPLACED IN THE INVENTORY; GROWTH FORMS; ENGLISH COLUMN NAMES
# =============================================================================

# --- 6a. Codes replaced in the inventory ------------------------------------
# A code of the inventory replaced by another code. The same replacement is
# made in 2_prepare_viewfulltable.R, before the species table is merged:
#   pterof -> pterro  all Pterocarpus in the plot are P. rohrii (Rolando, email 1)
# After the replacement pterof has no records in the plot, so its row is
# removed from the species table here.
mnemonic_replace <- c(pterof = "pterro")

# ficutr -> ficuau (Suzanne, email 2: the trees coded ficutr are Ficus aurea,
# whose code is ficuau) is not decided yet. Both codes are kept in the species
# table: ficutr with the trees of the inventory, and ficuau without records
# until the inventory is recoded (switch REPLACE_FICUTR_WITH_FICUAU in
# 2_prepare_viewfulltable.R). Both have the same growth form (6b).
codes_without_records <- "ficuau"

# rows of the replaced codes and of the codes that replace them
viewt_tax1_updated_flor_pan[mnemonic %in% c(names(mnemonic_replace), mnemonic_replace)]

check(
    all(mnemonic_replace %in% viewt_tax1_updated_flor_pan$mnemonic),
    "the codes that replace other codes have a row in the species table"
)

# remove the replaced codes
viewt_tax1_updated_flor_pan <- viewt_tax1_updated_flor_pan[!mnemonic %in% names(mnemonic_replace)]

# The species table must now hold exactly the codes of the inventory after the
# replacement, plus the codes kept without records (ficuau). The codes come from the current inventory (oct_13_2026), the
# same table 0_extract_ViewTaxonomy.R reads; only the Mnemonic column is loaded.
inventory_mnemonics <- unique(fread(ViewFullTable_new, sep = "\t", select = "Mnemonic")$Mnemonic)
inventory_mnemonics <- unique(fifelse(
    inventory_mnemonics %in% names(mnemonic_replace),
    mnemonic_replace[inventory_mnemonics],
    inventory_mnemonics
))
setdiff(viewt_tax1_updated_flor_pan$mnemonic, inventory_mnemonics) # in the table, not in the plot (ficuau)
setdiff(inventory_mnemonics, viewt_tax1_updated_flor_pan$mnemonic) # in the plot, not in the table
check(
    setequal(viewt_tax1_updated_flor_pan$mnemonic, c(inventory_mnemonics, codes_without_records)),
    "the species table holds exactly the codes present in the plot after the replacement, plus ficuau"
)
check(
    !anyDuplicated(viewt_tax1_updated_flor_pan$mnemonic),
    "one row per mnemonic in the species table"
)
nrow(viewt_tax1_updated_flor_pan)

# --- 6b. Growth forms --------------------------------------------------------
# The growth forms are those of the list of email 1: R. Pérez & S. Aguilar
# (f_de_vida_r_perez_s_aguilar, the column 2_prepare_viewfulltable.R reads)
# and R. Foster (f_de_vida_r_foster: T canopy tree, M midstory tree,
# U understory treelet, S shrub, TE strangler). The growth forms of the Flora
# de Panamá list are not used: they are incomplete (decided with the
# botanists, 2026-10-14).
# Downstream only the class is used (tree, shrub, palm, tree fern, strangler),
# read from the Spanish label: 2_prepare_viewfulltable.R corrects the DBH for
# taper when the label contains "árbol" or "arbusto" (or is empty), and the
# stem identification picks the growth parameters of each class.

# Species without a growth form of R. Pérez & S. Aguilar
viewt_tax1_updated_flor_pan[
    is.na(f_de_vida_r_perez_s_aguilar),
    .(mnemonic, genero_corrected, especie_corrected, f_de_vida_r_foster, f_de_vida_r_perez_s_aguilar)
]

# Growth forms given by the botanists (2026-10-14):
#   ficutr, ficuau          the same trees (see 6a), stranglers (the list gave
#                           ficutr "árbol grande mayor a 20m" and has no ficuau)
#   apeihy, nects1, nects3  trees
viewt_tax1_updated_flor_pan[mnemonic %in% c("ficutr", "ficuau"), f_de_vida_r_perez_s_aguilar := "árbol estrangulador mayor a 20m"]
viewt_tax1_updated_flor_pan[mnemonic %in% c("apeihy", "nects1", "nects3"), f_de_vida_r_perez_s_aguilar := "árbol"]

# R. Foster's code is empty for ficuau, apeihy, nects1 and nects3 (they are not
# in his list under these codes), and is T (canopy tree) for ficutr; change it
# by hand if the botanists give another value, e.g.
# viewt_tax1_updated_flor_pan[mnemonic == "ficuau", f_de_vida_r_foster := "TE"]

# Every code has a growth form, except uniden (unidentified stems, left out of
# the species table in section 8)
viewt_tax1_updated_flor_pan[is.na(f_de_vida_r_perez_s_aguilar), .(mnemonic, genero_corrected, especie_corrected)]
check(
    viewt_tax1_updated_flor_pan[mnemonic != "uniden", !anyNA(f_de_vida_r_perez_s_aguilar)],
    "every species has a growth form of R. Pérez & S. Aguilar"
)

# Class of each growth form, with the rules of the stem identification
# (1_main_cpp_chunk_bci.R): "estrangulador" before "árbol"; empty = tree
lifeform_class <- function(x) {
    x <- tolower(x)
    fcase(
        grepl("estrangulador", x), "strangler",
        grepl("árbol", x), "tree",
        grepl("palma", x), "palm",
        grepl("helecho", x), "tree fern",
        grepl("arbusto", x), "shrub"
    )
}
check(
    viewt_tax1_updated_flor_pan[!is.na(f_de_vida_r_perez_s_aguilar), !anyNA(lifeform_class(f_de_vida_r_perez_s_aguilar))],
    "every growth form of R. Pérez & S. Aguilar is read as tree, shrub, palm, tree fern or strangler"
)

# Growth forms by class, and against R. Foster's code
dcast(
    viewt_tax1_updated_flor_pan[, .N, by = .(class = lifeform_class(f_de_vida_r_perez_s_aguilar), foster = f_de_vida_r_foster)],
    class ~ foster,
    value.var = "N",
    fill = 0
)

# Ficus and other stranglers
viewt_tax1_updated_flor_pan[
    genero_corrected == "Ficus" | lifeform_class(f_de_vida_r_perez_s_aguilar) %in% "strangler",
    .(mnemonic, genero_corrected, especie_corrected, f_de_vida_r_foster, f_de_vida_r_perez_s_aguilar)
]

# --- 6c. English column names -----------------------------------------------
# Column names: the suffix "_corrected" is dropped, and the Spanish names (from
# "orden" onwards) are translated to English, one to one
colnames(viewt_tax1_updated_flor_pan) <- gsub("_corrected$", "", colnames(viewt_tax1_updated_flor_pan))
colnames(viewt_tax1_updated_flor_pan)

names_es_en <- c(
    orden = "order",
    familia = "family",
    genero = "genus",
    especie = "species",
    autoridad = "authority",
    sinonimos = "synonyms",
    nombre_comun = "common_name",
    herbario = "herbarium"
)

check(
    all(names(names_es_en) %in% colnames(viewt_tax1_updated_flor_pan)),
    "all Spanish column names to translate are in the species table"
)
check(
    !any(unname(names_es_en) %in% colnames(viewt_tax1_updated_flor_pan)),
    "no English column name is already in the species table"
)

setnames(viewt_tax1_updated_flor_pan, old = names(names_es_en), new = unname(names_es_en))

colnames(viewt_tax1_updated_flor_pan)

viewt_tax1_updated_flor_pan

# =============================================================================
# 7. NOTES ON KNOWN DATA PROBLEMS
# =============================================================================
# Notes exported with the species table (column Notes)
viewt_tax1_updated_flor_pan[, notes := NA_character_]

# Apeiba "hybrida" (apeihy): one tree, alive in 1982 and 1985, recorded dead
# in 1990.
viewt_tax1_updated_flor_pan[mnemonic %in% "apeihy", notes := "No longer present: recorded alive until 1985, dead in 1990"]
viewt_tax1_updated_flor_pan[mnemonic %in% "apeihy"]

# Nectandra s1 (nects1) and s3 (nects3) are Lauraceae morphospecies. nects3
# (one tree) was alive in 1982 and 1985 and recorded dead in 1990; nects1 is
# still alive (one tree in 2022-23) and should be followed to collect it
# fertile and identify it.
viewt_tax1_updated_flor_pan[mnemonic %in% "nects1"]
viewt_tax1_updated_flor_pan[mnemonic %in% "nects3"]
viewt_tax1_updated_flor_pan[mnemonic %in% "nects3", notes := "No longer present: recorded alive until 1985, dead in 1990"]

# pterro includes the trees recorded as pterof (Pterocarpus officinalis does
# not occur in the plot; R. Pérez, email 1)
viewt_tax1_updated_flor_pan[mnemonic %in% "pterro", notes := "Previously identified as 'pterof' (Pterocarpus officinalis) are incorrect; those are 'pterro' (Pterocarpus rohrii)"]

# =============================================================================
# 8. CHECK THE NAMES AGAINST THE CURRENT WFO BACKBONE
# =============================================================================
# Split every name into its parts and check each one against the World Flora
# Online (WFO) backbone of June 2026 (WFO_BACKBONE, a local file). The names in
# this table are those of the Flora de Panamá list and stay as the reference:
# the WFO record of each name is added in its own columns (wfo_*) and is used
# only to suggest corrections, mainly of spelling, that are adopted one by one
# in section 9. Nothing is changed automatically.
#
# Columns created here:
#   species_panama   species column exactly as it came from section 5
#                    (e.g. "simplex var. grandiflora"), kept for traceability
#   species          specific epithet only ("simplex")
#   rank             infraspecific rank ("var.", "subsp.", "f."); NA for species
#   variety          infraspecific name ("grandiflora"); NA for species
#   name_type        "species", "infraspecific" or "morphospecies"
#   wfo_check        TRUE for the rows whose name is checked against WFO
#   scientific_name  full name without authority, built from the parts above
# uniden (unidentified stems) is not a species and is left out of the table.
spp_new <- copy(viewt_tax1_updated_flor_pan[mnemonic != "uniden"])

# Helpers used below.
# TRUE where two values differ (NA vs NA is equal, NA vs a value differs)
differs <- function(x, y) {
    (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}
# Empty text fields as NA (the WFO backbone stores missing values as "")
blank_to_na <- function(x) {
    if (!is.character(x)) {
        return(x)
    }
    fifelse(!is.na(x) & trimws(x) == "", NA_character_, x)
}

# --- 8a. Format the names and split the infraspecific ones -------------------
# Genus capitalised and species epithet in lower case, without repeated spaces
# ("simplex  var. grandiflora" has two spaces in the Flora de Panamá list).
# Rows changed by this formatting:
genus_before <- spp_new$genus
species_before <- spp_new$species
spp_new[, genus := str_to_sentence(str_squish(genus))]
spp_new[, species := str_to_lower(str_squish(species))]
spp_new[differs(genus, genus_before) | differs(species, species_before), .(mnemonic, genus, species)]
rm(genus_before, species_before)

# Names with an infraspecific rank come as one string in `species`
# ("<species> <rank> <name>", e.g. swars1 and swars2, the two varieties of
# Swartzia simplex). They are split so that the two codes are told apart by
# `rank` + `variety` and not merged into the same binomial.
spp_new[, species_panama := species]

infra_ranks <- c("var.", "subsp.", "ssp.", "f.")
infra_pattern <- "^([a-z-]+) (var\\.|subsp\\.|ssp\\.|f\\.) ([a-z-]+)$"

spp_new[, rank := NA_character_]
spp_new[, variety := NA_character_]
spp_new[grepl(infra_pattern, species_panama), `:=`(
    species = sub(infra_pattern, "\\1", species_panama),
    rank = sub(infra_pattern, "\\2", species_panama),
    variety = sub(infra_pattern, "\\3", species_panama)
)]
# one spelling per rank
spp_new[rank == "ssp.", rank := "subsp."]

# rows that were split
spp_new[!is.na(rank), .(mnemonic, genus, species_panama, species, rank, variety, authority)]

# --- 8b. Classify the kind of name -------------------------------------------
# morphospecies: field names that are not published names ("sp.", "sp.4_...",
#                or a label such as "hybrid" / "fuzzy"). They are recognised by
#                the "sp." prefix or by having no authority, and are not
#                checked against WFO.
# infraspecific: species + rank + infraspecific name
# species:       genus + a single specific epithet
spp_new[, name_type := fcase(
    grepl("^sp\\.", species) | is.na(authority), "morphospecies",
    !is.na(rank), "infraspecific",
    grepl("^[a-z-]+$", species), "species"
)]
spp_new[, wfo_check := name_type %in% c("species", "infraspecific")]

spp_new[, .N, by = .(name_type, wfo_check)]
# morphospecies: set aside, not checked against WFO
spp_new[name_type == "morphospecies", .(mnemonic, family, genus, species, authority, notes)]

# Full name without authority, built from the separated parts. This is the
# Flora de Panamá name that is looked up in WFO.
spp_new[, scientific_name := fifelse(
    name_type %in% "infraspecific",
    paste(genus, species, rank, variety),
    paste(genus, species)
)]

# --- 8c. Checks on the separated names ---------------------------------------
check(
    !anyNA(spp_new$name_type),
    "every row is classified as species, infraspecific or morphospecies"
)
check(
    spp_new[wfo_check == TRUE, all(grepl("^[A-Z][a-z-]+$", genus) & grepl("^[a-z-]+$", species))],
    "names to check have a one-word genus and a one-word specific epithet"
)
check(
    spp_new[name_type == "infraspecific", all(rank %in% setdiff(infra_ranks, "ssp.") & grepl("^[a-z-]+$", variety))] &&
        spp_new[name_type != "infraspecific", all(is.na(rank) & is.na(variety))],
    "rank and variety are filled for infraspecific names and empty for all others"
)
check(
    !any(grepl("(^| )(var|subsp|ssp|f)\\.( |$)", spp_new$species)),
    "no infraspecific rank is left inside the species column"
)
check(
    !anyDuplicated(spp_new[wfo_check == TRUE, scientific_name]),
    "names to check are unique (no two codes share the same full name)"
)

# --- 8d. Authority -----------------------------------------------------------
# Every name to check needs its authority: it is compared with the WFO
# authority, and it decides between WFO records that share a name (8h).
check(
    spp_new[wfo_check == TRUE, !anyNA(authority)],
    "every name to check has an authority"
)
# Annotations that are not part of an authority (sensu, auct., nom. dub.,
# ined.) would distort that comparison. Rows that have one (should be empty;
# correct them by hand if not):
spp_new[
    grepl("sensu|auct\\.|nom\\. dub\\.|ined\\.", authority, ignore.case = TRUE),
    .(mnemonic, scientific_name, authority)
]

# The table, separated: reference names (Flora de Panamá) by part
spp_new[, .(mnemonic, name_type, wfo_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "species"]
spp_new[, .(mnemonic, name_type, wfo_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "morphospecies"]
spp_new[, .(mnemonic, name_type, wfo_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "infraspecific"]

# --- 8e. Load the current WFO backbone ---------------------------------------
# World Flora Online Plant List June 2026 (see WFO_BACKBONE in 0b). Read as
# WorldFlora expects it, with empty fields as "".
wfo_backbone <- fread(WFO_BACKBONE, encoding = "UTF-8")
check(
    all(c(
        "taxonID", "scientificName", "scientificNameAuthorship", "taxonRank",
        "family", "genus", "specificEpithet", "infraspecificEpithet",
        "taxonomicStatus", "acceptedNameUsageID"
    ) %in% names(wfo_backbone)),
    "the WFO backbone has the columns used below"
)
cat("WFO backbone", WFO_BACKBONE_VERSION, ":", format(nrow(wfo_backbone), big.mark = ","), "records\n")

# Accepted name of every record, to show the accepted name of a synonym (9h)
wfo_accepted_names <- wfo_backbone[, .(
    acceptedNameUsageID = taxonID,
    wfo_accepted_name = scientificName,
    wfo_accepted_author = scientificNameAuthorship
)]

# --- 8f. Names checked against WFO -------------------------------------------
# Only the names marked in 8b (wfo_check) are checked. The morphospecies are
# not published names; they are kept as they are (see 9i).
spp_new_to_check <- spp_new[wfo_check == TRUE]
spp_new_not_checked <- spp_new[wfo_check == FALSE]

cat(nrow(spp_new_not_checked), "rows are morphospecies and are not checked against WFO\n")
spp_new_not_checked[, .(mnemonic, family, genus, species)]
cat(nrow(spp_new_to_check), "rows are checked against WFO\n")

# Input of WFO.match(): the full name without authority, and the authority
wfo_input <- spp_new_to_check[, .(mnemonic, spec.name = scientific_name, Authorship = authority)]

# --- 8g. Match each name in the backbone (WorldFlora) ------------------------
# Kindt, R. (2020). WorldFlora: An R package for exact and fuzzy matching of
# plant names against the World Flora Online taxonomic backbone data.
# Applications in Plant Sciences 8(9): e11388.
# WFO.match() looks for the exact name first. When it is not found it looks for
# the closest spelling (fuzzy matching, Fuzzy = 0.1), and failing that for the
# first two words, or the genus only. It returns every record that carries the
# matched name (one row per record), with the distance between the Flora de
# Panamá authority and the authority of the record (Auth.dist).
# acceptedNameUsageID.match = FALSE keeps the record of the name itself, also
# when it is a synonym: the Flora de Panamá name has priority, and the accepted
# name of a synonym is only shown for information (9h).
wfo_candidates <- as.data.table(WFO.match(
    spec.data = as.data.frame(wfo_input),
    WFO.data = wfo_backbone,
    spec.name = "spec.name",
    Authorship = "Authorship",
    acceptedNameUsageID.match = FALSE,
    Fuzzy = 0.1,
    verbose = FALSE
))
wfo_candidates <- wfo_candidates[, lapply(.SD, blank_to_na)]

check(
    setequal(wfo_candidates$mnemonic, wfo_input$mnemonic),
    "WFO.match returned every name that was sent"
)

# Records per name: exact matches (Fuzzy FALSE) and spelling matches (Fuzzy TRUE)
wfo_candidates[, n_records := .N, by = mnemonic]
wfo_candidates[, .(names = uniqueN(mnemonic), records = .N), by = .(Matched, Fuzzy)]

# --- 8h. One record per name -------------------------------------------------
# Several WFO records can carry the same name (homonyms by other authors, or
# duplicates). The record kept is, in this order, the one whose authority is
# closest to the Flora de Panamá authority (smallest Auth.dist), an accepted
# one, and the one with the smallest taxonID (the rules of WorldFlora's
# WFO.one()).
wfo_candidates[, status_order := fcase(
    taxonomicStatus %in% "Accepted", 1L,
    taxonomicStatus %in% "Synonym", 2L,
    default = 3L
)]
setorder(wfo_candidates, mnemonic, Auth.dist, status_order, taxonID, na.last = TRUE)

wfo_choice <- wfo_candidates[, .(taxonID = taxonID[1]), by = mnemonic]

# Names with more than one record: every candidate, and the one kept
wfo_candidates[wfo_choice, on = "mnemonic", kept := !differs(taxonID, i.taxonID)]
wfo_candidates[
    n_records > 1,
    .(mnemonic, spec.name, Authorship, taxonID, scientificNameAuthorship, taxonomicStatus, Auth.dist, kept)
]

# To keep another record for a name, set its taxonID by hand, e.g.
# wfo_choice[mnemonic == "xxxxxx", taxonID := "wfo-0000000000"]
# Records chosen by hand (write the lines here):

# The record kept for each name (one of its candidates)
wfo_matched <- wfo_candidates[wfo_choice, on = .(mnemonic, taxonID)]
check(
    nrow(wfo_matched) == nrow(wfo_input) && !anyNA(wfo_matched$Matched),
    "one WFO record per name checked, chosen among its candidates"
)

# --- 8i. Compare each WFO record with the Flora de Panamá name ---------------
# One row per name checked: the name parts of the species table next to those
# of its WFO record, and four comparisons:
#   name_check    "same"; "spelling differs" (the name was found with another
#                 spelling of the genus, species or variety); or
#                 "not matched at its rank" (no record, or a record of another
#                 rank, e.g. only the genus)
#   author_check  "same"; "format only" (same text apart from spaces and
#                 capital letters, e.g. "C. DC." and "C.DC."); "different"; or
#                 "not compared" (the WFO record has no authority)
#   family_check  "same" or "differs"
#   wfo_status    status of the record in WFO (Accepted, Synonym, Unchecked)

# Rank of each name in the words of the WFO column taxonRank
rank_words <- c("var." = "variety", "subsp." = "subspecies", "f." = "form")

wfo_review <- merge(
    spp_new_to_check[, .(mnemonic, name_type, family, genus, species, rank, variety, authority, scientific_name)],
    wfo_matched[, .(
        mnemonic, Matched, Fuzzy, Fuzzy.dist, Auth.dist, n_records,
        wfo_id = taxonID, wfo_name = scientificName, wfo_author = scientificNameAuthorship,
        wfo_rank = taxonRank, wfo_family = family, wfo_genus = genus,
        wfo_species = specificEpithet, wfo_variety = infraspecificEpithet,
        wfo_status = taxonomicStatus, acceptedNameUsageID
    )],
    by = "mnemonic",
    all.x = TRUE
)
wfo_review <- merge(wfo_review, wfo_accepted_names, by = "acceptedNameUsageID", all.x = TRUE)
wfo_review[, acceptedNameUsageID := NULL]
setkey(wfo_review, mnemonic)

wfo_review[, submitted_rank := fifelse(is.na(rank), "species", unname(rank_words[rank]))]

wfo_review[, name_check := fcase(
    !Matched | is.na(wfo_rank) | differs(wfo_rank, submitted_rank),
    "not matched at its rank",
    differs(genus, wfo_genus) | differs(species, wfo_species) |
        (name_type == "infraspecific" & differs(variety, wfo_variety)),
    "spelling differs",
    default = "same"
)]

# Authority compared without spaces and capital letters for "format only"
author_key <- function(x) tolower(gsub(" ", "", x))
wfo_review[, author_check := fcase(
    is.na(wfo_author), "not compared",
    authority == wfo_author, "same",
    author_key(authority) == author_key(wfo_author), "format only",
    default = "different"
)]

wfo_review[, family_check := fifelse(differs(family, wfo_family), "differs", "same")]

# Version of the backbone, kept with the review
wfo_review[, wfo_backbone_version := WFO_BACKBONE_VERSION]

check(
    !anyNA(wfo_review$name_check) && !anyNA(wfo_review$author_check),
    "every name checked has a name and an author comparison"
)

# =============================================================================
# 9. REVIEW THE WFO SUGGESTIONS (STEP BY STEP)
# =============================================================================
# The Flora de Panamá names have priority, and nothing below changes them
# automatically. Each step prints what the WFO record suggests, and a
# suggestion is adopted only by copying its line (printed by
# print_accept_lines()) into the place indicated and running it.
# WFO is used mainly to correct spelling: the spelling of the name parts (9c)
# and the format of the authority (9d). Differences in content (9b, 9e, 9f)
# are checked by hand. Family and status in WFO are shown for information only
# (9g, 9h) and are not applied.
# If a record is chosen by hand in 8h, run 8i again before this section.

# Copy of the names before the review, to record what is changed (9j)
spp_panama <- copy(spp_new)

# print_accept_lines(): print, for each suggestion, the line that adopts it,
# ready to copy. Nothing is run. `suggestions` needs the columns mnemonic,
# column (column of spp_new to change) and suggested (new value).
print_accept_lines <- function(suggestions) {
    if (nrow(suggestions) == 0) {
        cat("No suggestions\n")
        return(invisible(NULL))
    }
    cat(sprintf(
        "spp_new[mnemonic == %s, %s := %s]\n",
        encodeString(suggestions$mnemonic, quote = '"'),
        suggestions$column,
        encodeString(suggestions$suggested, quote = '"')
    ), sep = "")
    invisible(NULL)
}

# --- 9a. Overview ------------------------------------------------------------
wfo_review[, .N, by = name_check]
wfo_review[, .N, by = author_check]
wfo_review[, .N, by = family_check]
wfo_review[, .N, by = wfo_status]
wfo_review[, .N, by = .(several_records = n_records > 1)] # record chosen in 8h

# --- 9b. Names not matched at their rank -------------------------------------
# No WFO record with this name, or only a record of another rank (e.g. the
# genus). Check these names in the Flora de Panamá list and the literature;
# keep them unless you find an error, and in that case correct the column by
# hand, e.g.
# spp_new[mnemonic == "xxxxxx", species := "..."]
wfo_review[
    name_check == "not matched at its rank",
    .(mnemonic, scientific_name, authority, wfo_name, wfo_rank, wfo_author)
]

# Corrections (write the lines here):

# --- 9c. Spelling of genus, species and variety ------------------------------
# The name was found in WFO with another spelling. Review each suggestion
# (distance: number of letters that differ in the full name; check the name in
# Tropicos or IPNI if unsure) and adopt only true spelling mistakes.
spelling_suggestions <- rbindlist(list(
    wfo_review[
        name_check == "spelling differs" & differs(genus, wfo_genus),
        .(mnemonic, column = "genus", current = genus, suggested = wfo_genus, distance = Fuzzy.dist)
    ],
    wfo_review[
        name_check == "spelling differs" & differs(species, wfo_species),
        .(mnemonic, column = "species", current = species, suggested = wfo_species, distance = Fuzzy.dist)
    ],
    wfo_review[
        name_check == "spelling differs" & name_type == "infraspecific" & differs(variety, wfo_variety),
        .(mnemonic, column = "variety", current = variety, suggested = wfo_variety, distance = Fuzzy.dist)
    ]
))
spelling_suggestions

# Lines that adopt each suggestion: copy the ones you adopt below and run them
print_accept_lines(spelling_suggestions)

# Adopted spelling corrections (paste the lines here):

# --- 9d. Authority: format only ----------------------------------------------
# Same authority written with other spaces or capital letters (WFO writes
# initials without spaces, e.g. "C.DC." for "C. DC."). Adopt the ones you want
# written as in WFO.
author_format <- wfo_review[
    author_check == "format only",
    .(mnemonic, scientific_name, column = "authority", current = authority, suggested = wfo_author)
]
author_format

# Lines that adopt each suggestion: copy the ones you adopt below and run them
print_accept_lines(author_format)

# Adopted authority formats (paste the lines here):
spp_new[mnemonic == "ade1tr", authority := "(Müll.Arg.) Hemsl."]
spp_new[mnemonic == "alchco", authority := "Pax & K.Hoffm."]
spp_new[mnemonic == "alibed", authority := "(Rich.) A.Rich. ex DC."]
spp_new[mnemonic == "amaico", authority := "(Lam. ex Poir.) Delprete & C.H.Perss."]
spp_new[mnemonic == "andiin", authority := "(W.Wright) Kunth ex DC."]
spp_new[mnemonic == "aspicr", authority := "Benth. ex Müll.Arg."]
spp_new[mnemonic == "ast1st", authority := "L.H.Bailey"]
spp_new[mnemonic == "bactba", authority := "L.H.Bailey"]
spp_new[mnemonic == "bactc1", authority := "L.H.Bailey"]
spp_new[mnemonic == "bactc2", authority := "L.H.Bailey"]
spp_new[mnemonic == "beilpe", authority := "(Klotzsch & H.Karst. ex Meisn.) Sachiko Nishida"]
spp_new[mnemonic == "boropa", authority := "Delprete & C.H.Perss."]
spp_new[mnemonic == "brosgu", authority := "(Aubl.) Huber ex Ducke"]
spp_new[mnemonic == "carica", authority := "(Jacq.) A.DC."]
spp_new[mnemonic == "cestme", authority := "G.Don"]
spp_new[mnemonic == "chr1ec", authority := "L.O.Williams"]
spp_new[mnemonic == "clidde", authority := "(D.Don) Michelang."]
spp_new[mnemonic == "cnempe", authority := "(Hook.) R.M.Tryon"]
spp_new[mnemonic == "cordbi", authority := "A.DC."]
spp_new[mnemonic == "cou2cu", authority := "(Duchass. ex Griseb.) C.M.Taylor"]
spp_new[mnemonic == "couthe", authority := "(Jacq.) K.Schum."]
spp_new[mnemonic == "crotbi", authority := "Müll.Arg."]
spp_new[mnemonic == "des2pa", authority := "(B.L.Rob.) Saff."]
spp_new[mnemonic == "drypst", authority := "G.L.Webster"]
spp_new[mnemonic == "entesc", authority := "(Benth.) A.R.M.Luz & E.R.Souza"]
spp_new[mnemonic == "eugega", authority := "(C.Wright ex Griseb.) Krug & Urb."]
spp_new[mnemonic == "eugeoe", authority := "O.Berg"]
spp_new[mnemonic == "faraoc", authority := "(L.) A.Rich."]
spp_new[mnemonic == "gustsu", authority := "(Kunth) O.Berg"]
spp_new[mnemonic == "hampap", authority := "(Donn.Sm.) Standl."]
spp_new[mnemonic == "herrpu", authority := "(Pittier) R.E.Schult."]
spp_new[mnemonic == "ingaru", authority := "G.Don"]
spp_new[mnemonic == "jac1co", authority := "(Aubl.) D.Don"]
spp_new[mnemonic == "jac2sp", authority := "(Aubl.) A.DC."]
spp_new[mnemonic == "koanwe", authority := "(B.L.Rob.) R.M.King & H.Rob."]
spp_new[mnemonic == "laciag", authority := "(P.J.Bergius) Rusby"]
spp_new[mnemonic == "lindla", authority := "C.Presl"]
spp_new[mnemonic == "lozapi", authority := "(S.F.Blake) L.B.Sm."]
spp_new[mnemonic == "micoho", authority := "Donn.Sm."]
spp_new[mnemonic == "micoim", authority := "(Sw.) D.Don"]
spp_new[mnemonic == "neeaam", authority := "Donn.Sm."]
spp_new[mnemonic == "oenoma", authority := "H.Karst."]
spp_new[mnemonic == "olmeas", authority := "(Pittier) C.C.Berg"]
spp_new[mnemonic == "phoeci", authority := "(Sw.) R.Rohde"]
spp_new[mnemonic == "pipeca", authority := "C.DC."]
spp_new[mnemonic == "pipeco", authority := "(Miq.) C.DC."]
spp_new[mnemonic == "pipecu", authority := "C.DC."]
spp_new[mnemonic == "pipeim", authority := "C.DC."]
spp_new[mnemonic == "pipelo", authority := "C.DC."]
spp_new[mnemonic == "pit1ma", authority := "(Pittier) M.V.B.Soares, M.P.Morim & Iganci"]
spp_new[mnemonic == "pochqu", authority := "(Seem.) W.S.Alverson & M.C.Duarte"]
spp_new[mnemonic == "protpa", authority := "(Rose) I.M.Johnst."]
spp_new[mnemonic == "psycli", authority := "K.Krause"]
spp_new[mnemonic == "psycra", authority := "(Aubl.) G.Nicholson"]
spp_new[mnemonic == "sch1zo", authority := "(Mutis ex L.f.) Wess.Boer"]
spp_new[mnemonic == "schipa", authority := "(Vell.) S.F.Blake"]
spp_new[mnemonic == "sennda", authority := "(Britton & Rose) H.S.Irwin & Barneby"]
spp_new[mnemonic == "sipacr", authority := "(Poepp. & Endl.) A.DC."]
spp_new[mnemonic == "sipapa", authority := "(Beurl.) A.DC."]
spp_new[mnemonic == "socrex", authority := "(Mart.) H.Wendl."]
spp_new[mnemonic == "sponra", authority := "Donn.Sm."]
spp_new[mnemonic == "sterap", authority := "(Jacq.) H.Karst."]
spp_new[mnemonic == "sympgl", authority := "L.f."]
spp_new[mnemonic == "tab1gu", authority := "(Seem.) S.O.Grose"]
spp_new[mnemonic == "tachve", authority := "van der Werff & N.Zamora"]
spp_new[mnemonic == "termam", authority := "(J.F.Gmel.) Exell"]
spp_new[mnemonic == "tripcu", authority := "Fisch. & C.A.Mey. ex C.A.Mey."]
spp_new[mnemonic == "turpoc", authority := "(Sw.) G.Don"]
spp_new[mnemonic == "virosp", authority := "D.Santam."]
spp_new[mnemonic == "xyl2ch", authority := "Donn.Sm."]
spp_new[mnemonic == "xyl2ol", authority := "Donn.Sm."]
spp_new[mnemonic == "zantp1", authority := "P.Wilson"]
spp_new[mnemonic == "zantse", authority := "P.Wilson"]
spp_new[mnemonic == "zuelgu", authority := "(A.Rich.) Warb."]

# --- 9e. Authority: different text -------------------------------------------
# The authority differs in content. `wfo` is the authority of the WFO record,
# not a recommendation: where WFO adds the author of the basionym in
# parentheses, or corrects a spelling, it is likely right; where it drops an
# author that the Flora de Panamá list has, keep the Flora de Panamá authority.
# Check in IPNI or Tropicos when unsure.
author_different <- wfo_review[
    author_check == "different",
    .(mnemonic, scientific_name, column = "authority", current = authority, wfo = wfo_author, wfo_status)
]
author_different

# Lines that copy the WFO authority as it is (only where WFO is right)
print_accept_lines(author_different[, .(mnemonic, column, suggested = wfo)])

# Adopted or hand-corrected authorities (paste or write the lines here):
# spp_new[mnemonic == "anacex", authority := "Skeels"]
spp_new[mnemonic == "ficuc2", authority := "(Liebm.) Miq."]
spp_new[mnemonic == "heisco", authority := "Standl."]
spp_new[mnemonic == "maclti", authority := "(L.) D.Don ex G.Don"]
spp_new[mnemonic == "ouralu", authority := "(Kunth) Engl."]
spp_new[mnemonic == "pse2sp", authority := "(Sw.) Griseb."]
spp_new[mnemonic == "tab1ro", authority := "(Bertol.) DC."]
# spp_new[mnemonic == "thevah", authority := "A.DC."]
spp_new[mnemonic == "thevah", authority := "(L.) Vahl"]
spp_new[mnemonic == "tri2tu", authority := "(Triana & Planch.) C.DC."]
spp_new[mnemonic == "tropra", authority := "(L.) Urb."]
spp_new[mnemonic == "zantpr", authority := "(Sw.) Sw."]

# --- 9f. Authority not compared ----------------------------------------------
# The WFO record has no authority. Check these authorities by eye and correct
# them by hand if needed, e.g.
# spp_new[mnemonic == "xxxxxx", authority := "..."]
wfo_review[
    author_check == "not compared",
    .(mnemonic, scientific_name, authority, wfo_id)
]

# Corrections (write the lines here):

# --- 9g. Family (information only) -------------------------------------------
# The family of the WFO record differs from the Flora de Panamá family. The
# Flora de Panamá family is kept.
wfo_review[
    family_check == "differs",
    .(mnemonic, scientific_name, family, wfo_family)
]

# --- 9h. Names that are not accepted in WFO (information only) ---------------
# Flora de Panamá has priority, so these names are kept. The status and the
# accepted name in WFO are shown for information (e.g. WFO treats the two
# varieties of Swartzia simplex as synonyms of the species).
wfo_review[
    differs(wfo_status, "Accepted"),
    .(mnemonic, scientific_name, wfo_status, wfo_accepted_name, wfo_accepted_author)
]

# --- 9i. Morphospecies (not checked against WFO) -----------------------------
# Kept as they are; their missing fields are completed in section 10.
spp_new_not_checked[, .(mnemonic, family, genus, species, notes)]

# --- 9j. Changes made to the Flora de Panamá names ---------------------------
# Rebuild the full name from the (possibly corrected) parts
spp_new[, scientific_name := fifelse(
    name_type %in% "infraspecific",
    paste(genus, species, rank, variety),
    paste(genus, species)
)]

# Every value that now differs from the Flora de Panamá name, one row per change
name_cols <- c("family", "genus", "species", "rank", "variety", "authority")
accepted_changes <- rbindlist(lapply(name_cols, function(col) {
    before <- spp_panama[, .(mnemonic, panama = get(col))]
    after <- spp_new[, .(mnemonic, now = get(col))]
    changed <- merge(before, after, by = "mnemonic")[differs(panama, now)]
    changed[, column := col]
    changed
}))
setcolorder(accepted_changes, c("mnemonic", "column", "panama", "now"))
accepted_changes

# Save the review, the candidate records and the changes made
checks_folder <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS")
fwrite(wfo_review, file.path(checks_folder, "species_wfo_review.csv"))
fwrite(
    wfo_candidates[n_records > 1, .(mnemonic, spec.name, Authorship, taxonID, scientificName, scientificNameAuthorship, taxonomicStatus, Auth.dist, kept)],
    file.path(checks_folder, "species_wfo_candidates.csv")
)
fwrite(accepted_changes, file.path(checks_folder, "species_wfo_accepted_changes.csv"))

# =============================================================================
# 10. COMPLETE AND CHECK THE REMAINING FIELDS
# =============================================================================

# --- 10a. Order --------------------------------------------------------------
# Rows without an order (morphospecies not in the Flora de Panamá list)
spp_new[is.na(order), .(mnemonic, family, genus, species)]

# The order is taken from the other species of the same family
family_orders <- spp_new[!is.na(order), .(n_orders = uniqueN(order), order = order[1]), by = family]
check(all(family_orders$n_orders == 1), "each family has a single order in the table")
family_orders[family %in% spp_new[is.na(order), family]]

order_malvaceae <- family_orders[family == "Malvaceae", order]
order_lauraceae <- family_orders[family == "Lauraceae", order]
spp_new[is.na(order) & family == "Malvaceae", order := order_malvaceae]
spp_new[is.na(order) & family == "Lauraceae", order := order_lauraceae]

spp_new[is.na(order), .(mnemonic, family, genus, species)] # should be empty
check(!anyNA(spp_new$order), "every species has an order")

# --- 10b. Growth form ---------------------------------------------------------
# Growth forms of R. Pérez & S. Aguilar, completed in 6b
spp_new[, .N, by = f_de_vida_r_perez_s_aguilar][order(-N)]
check(
    !anyNA(spp_new$f_de_vida_r_perez_s_aguilar),
    "every species of the table has a growth form of R. Pérez & S. Aguilar"
)

# --- 10c. Synonyms -----------------------------------------------------------
# A synonym must not repeat the name itself. For a variety, the old name of
# the same variety is a valid synonym (e.g. Swartzia simplex var. ochnacea for
# swars2, its name in the database until October 2026).
spp_new[
    !is.na(synonyms) &
        unname(mapply(
            function(name, syn) name %in% trimws(strsplit(syn, "[,;]")[[1]]),
            scientific_name,
            synonyms
        )),
    .(mnemonic, scientific_name, synonyms)
]

# Remove the repeated name by hand, e.g.
# spp_new[mnemonic == "xxxxxx", synonyms := NA_character_]

# =============================================================================
# 11. ASSEMBLE THE SPECIES TABLE AND EXPORT
# =============================================================================
# Columns of the species table (left) and their names in the exported table
# (right). The names are those used by the downstream scripts.
cols_export <- c(
    mnemonic = "Mnemonic", # BCI mnemonic code (primary key)
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "SpeciesName", # specific epithet
    rank = "InfraspecificRank", # e.g. "var."; NA for species
    variety = "InfraspecificEpithet", # NA for species
    authority = "Authority",
    synonyms = "Synonyms",
    f_de_vida_r_foster = "Lifeform_RFoster", # growth form, R. Foster
    f_de_vida_r_perez_s_aguilar = "Lifeform_RPerez_SAguilar", # growth form, R. Pérez & S. Aguilar
    common_name = "CommonName",
    herbarium = "Herbarium", # herbarium vouchers
    notes = "Notes"
)

full_out <- spp_new[, names(cols_export), with = FALSE]
setnames(full_out, old = names(cols_export), new = unname(cols_export))
setorder(full_out, Mnemonic)

names(full_out)

check(!anyDuplicated(full_out$Mnemonic), "one row per mnemonic in the exported table")
check(
    setequal(full_out$Mnemonic, setdiff(c(inventory_mnemonics, codes_without_records), "uniden")),
    "the exported table holds every code of the plot except uniden, plus ficuau"
)
check(
    full_out[, !anyNA(Order) && !anyNA(Family) && !anyNA(Genus) && !anyNA(SpeciesName)],
    "order, family, genus and species are filled for every row"
)
check(
    full_out[, all(is.na(InfraspecificRank) == is.na(InfraspecificEpithet))],
    "infraspecific rank and epithet are both filled or both empty"
)

# Rows per growth form (R. Pérez & S. Aguilar)
full_out[, .N, by = Lifeform_RPerez_SAguilar]

# --- 11a. Export -------------------------------------------------------------
bci.spptable <- full_out

output_path <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "SPP_TABLE")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Tab-delimited plain text — portable, version-control friendly
fwrite(
    bci.spptable,
    file.path(output_path, "bci_spptable.txt"),
    sep = "\t"
)
message("Exported: DATA/SPP_TABLE/bci_spptable.txt")

# R binary format — for direct use in downstream R scripts
save(
    bci.spptable,
    file = file.path(output_path, "bci_spptable.RData")
)
message("Exported: DATA/SPP_TABLE/bci_spptable.RData")

# Same table with a .csv extension (also tab-separated)
fwrite(
    bci.spptable,
    file.path(output_path, "bci_spptable.csv"),
    sep = "\t"
)
message("Exported: DATA/SPP_TABLE/bci_spptable.csv")
