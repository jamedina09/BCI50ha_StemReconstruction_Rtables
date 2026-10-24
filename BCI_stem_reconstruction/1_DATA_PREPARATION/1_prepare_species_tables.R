# =============================================================================
# 1_prepare_species_tables.R
#
# Purpose: Prepare the BCI 50-ha plot species table, one row per species code
#          (mnemonic) of the plot. The names come from the botanists: the
#          species list of R. Pérez and S. Aguilar (email 1), the name changes
#          of S. Lao and the reply of R. Pérez (email 2), and the Flora de
#          Panamá list, which has priority. The database export
#          (ViewFullTable) gives the species codes; its names are not up to
#          date and are used only for the three codes that no list names (the
#          morphospecies apeihy, nects1 and nects3), checked against the email
#          of R. Pérez. The names are then checked, step by step, against the
#          current World Flora Online (WFO) backbone (a local copy of the June
#          2026 release) with the WorldFlora package. WFO only suggests
#          corrections, mainly of spelling and of the format of the authority;
#          each one is adopted by hand, nothing is changed automatically.
#
# Run from the project root (or from BCI_stem_reconstruction/), line by line.
#
# Sections:
#    0. Setup: packages, check helper, input paths
#    1. Load the data
#    2. Species codes of the plot and their database names
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
#   oct_13_2026/ViewFullTable_bci.csv              current inventory export: its
#                                                  species codes, and the names
#                                                  of the codes no list names
#   sp_tables/email_1/Lista_bci_mnemonics_formadevida.xlsx  list of R. Pérez & S. Aguilar
#   sp_tables/email_2/bci_name_changes.xlsx        name changes of S. Lao
#   sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx
#   sp_tables/wfo_plant_list_2026_06/classification.csv   WFO backbone, June 2026
#                                                  (Zenodo doi:10.5281/zenodo.20782718)
# Outputs:
#   DATA/SPP_TABLE/bci_spptable.txt, bci_spptable.csv (both tab-separated) and
#   bci_spptable.RData (object bci.spptable): one row per code of the plot
#   except uniden, with Mnemonic, Order, Family, Genus, SpeciesName,
#   InfraspecificRank, InfraspecificEpithet, Authority, Synonyms,
#   Lifeform_RFoster, Lifeform_RPerez_SAguilar, CommonName, Herbarium, Notes.
#   DATA/CHECKS/species_wfo_review.csv (each name next to its WFO record),
#   species_wfo_candidates.csv (names with several WFO records),
#   species_wfo_accepted_changes.csv (changes adopted in the review),
#   species_wfo_vs_flora_panama.csv (for the botanists: names whose WFO record
#   differs from the Flora de Panamá in family, genus, species, infraspecific
#   name or authority, part by part),
#   species_table_vs_wfo.csv (names of the species table, after the review,
#   that still differ from their WFO record, or that WFO treats as synonyms),
#   species_synonyms_vs_wfo.csv (each synonym of the table and whether WFO
#   places it in the same taxon; information only),
#   prepare_species_tables_log.txt (failed checks).
# Packages: data.table, WorldFlora, stringr, janitor, readxl.
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
library(janitor) # clean_names() for the spreadsheet column names

# Number of threads data.table will use
data.table::getDTthreads()

# --- 0b. Paths ---------------------------------------------------------------
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
    workspace_root <- dirname(workspace_root)
}

# Current inventory export: its species codes (Mnemonic), and the names of the
# codes that no list names (section 3)
ViewFullTable_new <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026", "ViewFullTable_bci.csv")
# Species list of R. Pérez and S. Aguilar (email 1)
TAXONOMY_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_1", "Lista_bci_mnemonics_formadevida.xlsx")
# Name changes of S. Lao (email 2)
TAXONOMY_2 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_2", "bci_name_changes.xlsx")
# Flora de Panamá list (R. Pérez, received 2026-10-09)
flora_panama_oct_9_2026 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "flora_panama_oct_9_2026", "Lista de especies de la Flora de Panamá_nov25 original.xlsx")
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
check(
    file.exists(c(ViewFullTable_new, TAXONOMY_1, TAXONOMY_2, flora_panama_oct_9_2026, WFO_BACKBONE)),
    "every input file exists"
)
tax1 <- as.data.table(readxl::read_xlsx(TAXONOMY_1))
tax2 <- as.data.table(readxl::read_xlsx(TAXONOMY_2))
flor_pan <- as.data.table(readxl::read_xlsx(flora_panama_oct_9_2026))

# =============================================================================
# 2. SPECIES CODES OF THE PLOT AND THEIR DATABASE NAMES
# =============================================================================
# The six-letter codes recorded in the inventory export, with the names the
# database gives them. The database names (Family, Genus, SpeciesName) are not
# up to date: they are used only for the codes that no list names (section 3).
# Census 9 (2020, partial), which 2_prepare_viewfulltable.R removes, adds no
# code.
db_names <- unique(fread(
    ViewFullTable_new,
    sep = "\t",
    select = c("Mnemonic", "Family", "Genus", "SpeciesName"),
    na.strings = c("", "NA", "NULL")
))
plot_codes <- sort(unique(db_names$Mnemonic))
check(
    !anyNA(plot_codes) && all(nchar(plot_codes) == 6L),
    "every species code of the plot is a six-letter code"
)
length(plot_codes) # 329 codes, uniden (unidentified stems) included (export of 2026-10-13)

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

# The list has one row per code, and every code of the list is a code of the
# plot. Codes of the plot not in the list (329 against 327): pterof (named in
# the Flora de Panamá list, section 5, and replaced by pterro in section 6a)
# and uniden (unidentified stems, left out of the table in section 8).
check(!anyDuplicated(tax1$codigo), "one row per code in the list of email 1")
setdiff(tax1$codigo, plot_codes)
check(all(tax1$codigo %in% plot_codes), "every code of the list of email 1 is a code of the plot")
setdiff(plot_codes, tax1$codigo)

# One row per code of the plot, with the columns of the list
spp_email <- merge(
    data.table(Mnemonic = plot_codes),
    tax1,
    by.x = "Mnemonic",
    by.y = "codigo",
    all.x = TRUE
)

# The corrections of the email, one by one, recorded in `email_1_comment`

# 1) El nombre Appunia siebertiii es correcto y actualizado, Morinda siebertii es un sinónimo.
spp_email[genero == "Appunia", email_1_comment := "Appunia siebertii es correcto; Morinda siebertii es un sinónimo"]

# 2) Apeiba "hybrida" murió.
spp_email[Mnemonic == "apeihy", email_1_comment := "Apeiba hybrida murió"]

# 3) Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una
# con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla.
spp_email[Mnemonic == "nects1", email_1_comment := "Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla."]
spp_email[Mnemonic == "nects3", email_1_comment := "Nectandra s1 y Nectandra s3 son dos morfoespecies de Lauraceae, una murió. Solo queda una con vida y sugiero darle un seguimiento para colectarla fértil y poder identificarla."]

# 4) Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
spp_email[Mnemonic %in% c("pterof", "pterro"), .(Mnemonic, familia, genero, especie)]
spp_email[Mnemonic == "pterof", email_1_comment := "Pterocarpus officinalis no existe en la parcela; corresponde a P. rohrii"]
spp_email[Mnemonic == "pterro", email_1_comment := "Pterocarpus officinalis no existe en la parcela; corresponde a P. rohrii"]

# 5) Beilschmiedia pendula y Quararibea asterolepis son especies que no existen para BCI. Las
# especies correctas son Beilschmiedia tovarensis y Quararibea stenophylla, ambas fueron
# confirmadas por los especialistas en nuestras muestras de herbario.
spp_email[genero == "Beilschmiedia", .(Mnemonic, familia, genero, especie)]
spp_email[genero == "Beilschmiedia" & especie == "tovarensis", email_1_comment := "Beilschmiedia pendula no existe para BCI; la especie correcta es Beilschmiedia tovarensis"]
spp_email[genero == "Quararibea", .(Mnemonic, familia, genero, especie)]
spp_email[genero == "Quararibea" & especie == "stenophylla", email_1_comment := "Quararibea asterolepis no existe para BCI; la especie correcta es Quararibea stenophylla"]

spp_email[!is.na(email_1_comment), .(Mnemonic, familia, genero, especie, email_1_comment)]

# The growth forms (R. Foster; R. Pérez & S. Aguilar), common names and
# herbarium collections come with the list (its columns).

# Three codes of the list have no name at all (no family, genus or species):
# the morphospecies apeihy, nects1 and nects3, which are not in the Flora de
# Panamá list either. Their names are taken from the database export (section
# 2) and checked against the email (points 2 and 3 above): the genus must be
# that of the email, and the family that of the email (Lauraceae for the two
# Nectandra) and of the other species of the genus in the list. The labels
# differ and are not compared: the database has those of Condit et al. (2019),
# "hybrid", "fuzzy" and "sp.4_(tiny_leaf)", and the email "hybrida", s1 and
# s3. (The other codes of the plot without a name are not in the list: pterof,
# named by the Flora de Panamá list in section 5, and uniden.) The order of
# the three is completed in section 10.
codes_without_name <- spp_email[
    Mnemonic %in% tax1$codigo & is_empty(familia) & is_empty(genero) & is_empty(especie),
    Mnemonic
]
codes_without_name

# The three codes as the email names them
names_email_1 <- data.table(
    Mnemonic = c("apeihy", "nects1", "nects3"),
    name_email_1 = c('Apeiba "hybrida"', "Nectandra s1", "Nectandra s3"),
    genus_email_1 = c("Apeiba", "Nectandra", "Nectandra"),
    family_email_1 = c(NA, "Lauraceae", "Lauraceae")
)
check(
    setequal(codes_without_name, names_email_1$Mnemonic),
    "the codes without a name in the list are the morphospecies of the email (apeihy, nects1, nects3)"
)

# Their names in the database export
names_db <- db_names[Mnemonic %in% codes_without_name]
check(
    setequal(names_db$Mnemonic, codes_without_name) && !anyDuplicated(names_db$Mnemonic) &&
        !names_db[, anyNA(Family) || anyNA(Genus) || anyNA(SpeciesName)],
    "the database export gives one full name (family, genus, species) to each code without a name"
)

# Family of the other species of each genus in the list
genus_family <- unique(spp_email[
    genero %in% names_db$Genus & !Mnemonic %in% codes_without_name & !is_empty(familia),
    .(Genus = genero, family_list = familia)
])
genus_family # Apeiba: Malvaceae; Nectandra: Lauraceae
check(
    !anyDuplicated(genus_family$Genus) && setequal(genus_family$Genus, names_db$Genus),
    "each genus of the codes without a name has other species in the list, all of one family"
)

# The database names next to the email and the list
names_db <- merge(names_db, names_email_1, by = "Mnemonic")
names_db <- merge(names_db, genus_family, by = "Genus")
names_db[, .(Mnemonic, Family, Genus, SpeciesName, name_email_1, family_email_1, family_list)]
check(
    names_db[, all(Genus == genus_email_1)],
    "the database genus of each code without a name is the genus of the email"
)
check(
    names_db[, all((is.na(family_email_1) | Family == family_email_1) & Family == family_list)],
    "the database family of each code without a name is that of the email and of its genus in the list"
)

spp_email[names_db, on = "Mnemonic", `:=`(familia = i.Family, genero = i.Genus, especie = i.SpeciesName)]
spp_email[Mnemonic %in% codes_without_name, .(Mnemonic, familia, genero, especie, email_1_comment)]

# Every code of the list now has a family, a genus and a species
check(
    spp_email[Mnemonic %in% tax1$codigo, !any(is_empty(familia) | is_empty(genero) | is_empty(especie))],
    "every code of the list of email 1 has a family, a genus and a species"
)

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
# NOTE: this change is not applied. The new code is not a code of the plot
# NOTE: this change was flagged as incorrect by the team after discussion
# (the database has no such code and no tree carries it): the trees keep the
# code ficutr.
# Me avisa si este listado es diferente al que le envió Rolando. Es posible que el listado de
# Rolando esté más actualizado.

# Gracias y saludos,
# Suzanne

# The table has one row per code; Suzanne's list has 330 rows, with two for
# protte
spp_email[, .N, by = Mnemonic][order(-N)]
tax2[, .N, by = mnemonic][order(-N)]

# protte: one row with the SubspeciesID 15 but no subspecies name, one with
# neither; the first is dropped
tax2[mnemonic %in% tax2[duplicated(mnemonic), mnemonic]]
tax2 <- tax2[!(
    mnemonic %in% tax2[duplicated(mnemonic), mnemonic] &
        is_empty(subspecies) & !is_empty(subspeciesid)
)]
check(!anyDuplicated(tax2$mnemonic), "one row per code in the list of S. Lao")

setorder(spp_email, Mnemonic)
setorder(tax2, mnemonic)

# Codes in only one of the two: ficutr in the table and, in Suzanne's list, the
# new code she gives to ficutr (not a code of the plot, not used)
setdiff(spp_email$Mnemonic, tax2$mnemonic)
setdiff(tax2$mnemonic, spp_email$Mnemonic)

spp_email[Mnemonic == "ficutr"]
tax2[mnemonic %in% setdiff(tax2$mnemonic, spp_email$Mnemonic)]

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
spp_email[Mnemonic %in% c("stemgr", "pla1pi")]

spp_email[Mnemonic %in% c("pla1pi"), sinonimos := "Platymiscium polystachyum"]
spp_email[Mnemonic %in% c("stemgr", "pla1pi"), .(Mnemonic, genero, especie, autoridad, sinonimos)]

# Family, genus and species: the table vs Suzanne's list. The table (the list
# of email 1) is kept: Suzanne says that Rolando's list may be more up to
# date. The rows printed are the differences, for information.
mefam_susfam <- merge(
    spp_email[, .(Mnemonic, familia)],
    tax2[, .(mnemonic, family, familyold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
mefam_susfam[familia != family]

megen_susgen <- merge(
    spp_email[, .(Mnemonic, genero)],
    tax2[, .(mnemonic, genus, genusold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
megen_susgen[genero != genus]

mesp_sussp <- merge(
    spp_email[, .(Mnemonic, especie)],
    tax2[, .(mnemonic, speciesname, speciesnameold)],
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)
mesp_sussp[especie != speciesname]

# Suzanne's subspecies (information only: the infraspecific names of the table
# come from the Flora de Panamá list, section 5)
tax2[!is_empty(subspecies) | !is_empty(subspeciesold), .(mnemonic, genus, speciesname, subspecies, subspeciesold)]

# ficutr: the list of email 1 and the Flora de Panamá list (section 5) give
# Ficus trigonata, and the table keeps it; Suzanne's list renames the code
# (not applied, see above). The trees keep the code ficutr.
spp_email[Mnemonic == "ficutr", .(Mnemonic, familia, genero, especie, autoridad)]
tax2[mnemonic %in% c("ficutr", setdiff(tax2$mnemonic, spp_email$Mnemonic)), .(mnemonic, family, genus, speciesname, genusold, speciesnameold)]

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

spp_email

# Column names in lower case without spaces or accents; the Flora de Panamá
# columns get the prefix "pan_" so both tables can be merged
flor_pan <- clean_names(flor_pan)
spp_email <- clean_names(spp_email)
colnames(flor_pan) <- paste0("pan_", colnames(flor_pan))

# The text "NULL" means no value
spp_email <- spp_email[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]
flor_pan <- flor_pan[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]

# Each code must appear only once in the Flora de Panamá list, otherwise the
# merge would duplicate rows
check(
    !anyDuplicated(flor_pan[pan_codigo %in% spp_email$mnemonic, pan_codigo]),
    "each code of the table appears at most once in the Flora de Panamá list"
)

spp_flora <- merge(
    spp_email,
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
spp_flora[, texto := NULL][, fotos := NULL][, pan_texto := NULL][, pan_fotos := NULL][, pan_mnbc := NULL]
spp_flora[, pan_forma_de_vida := NULL]

# Codes of the plot without a row in the Flora de Panamá list: the
# morphospecies apeihy, nects1, nects3 and uniden
setdiff(spp_flora$mnemonic, flor_pan$pan_codigo)
check(setequal(spp_flora$mnemonic, plot_codes), "the table holds every code of the plot")

# Each field: the Flora de Panamá value where it has one, otherwise the value
# of the table (`<field>_corrected`). The rows printed before and after are
# those where the two sources differ or one is missing.

# orden
spp_flora[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden)]
spp_flora[, orden_corrected := fcoalesce(pan_orden, orden)]
spp_flora[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden, orden_corrected)]
spp_flora[, orden := NULL][, pan_orden := NULL]

# familia
spp_flora[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia)]
spp_flora[, familia_corrected := fcoalesce(pan_familia, familia)]
spp_flora[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia, familia_corrected)]
spp_flora[, familia := NULL][, pan_familia := NULL]

# genero
spp_flora[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero)]
spp_flora[, genero_corrected := fcoalesce(pan_genero, genero)]
spp_flora[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero, genero_corrected)]
spp_flora[, genero := NULL][, pan_genero := NULL]

# especie (for the two varieties of Swartzia simplex the Flora de Panamá
# species includes the variety, e.g. "simplex var. grandiflora"; it is split in
# section 8)
spp_flora[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie)]
spp_flora[, especie_corrected := fcoalesce(pan_especie, especie)]
spp_flora[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie, especie_corrected)]
spp_flora[, especie := NULL][, pan_especie := NULL]

# autoridad
spp_flora[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad)]
spp_flora[, autoridad_corrected := fcoalesce(pan_autoridad, autoridad)]
spp_flora[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad, autoridad_corrected)]
spp_flora[, autoridad := NULL][, pan_autoridad := NULL]

# sinonimos
spp_flora[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos)]
spp_flora[, sinonimos_corrected := fcoalesce(pan_sinonimos, sinonimos)]
spp_flora[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos, sinonimos_corrected)]
spp_flora[, sinonimos := NULL][, pan_sinonimos := NULL]

# nombre comun (the column is spelled "nombre_comum" in both lists)
spp_flora[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum)]
spp_flora[, nombre_comun_corrected := fcoalesce(pan_nombre_comum, nombre_comum)]
spp_flora[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum, nombre_comun_corrected)]
spp_flora[, nombre_comum := NULL][, pan_nombre_comum := NULL]

# herbario
spp_flora[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario)]
spp_flora[, herbario_corrected := fcoalesce(pan_herbario, herbario)]
spp_flora[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario, herbario_corrected)]
spp_flora[, herbario := NULL][, pan_herbario := NULL]

# Spot checks: rows discussed in the emails
spp_flora[mnemonic == "appuse"]
spp_flora[mnemonic == "ficutr"]
spp_flora[!is.na(email_1_comment)]

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

# rows of the replaced codes and of the codes that replace them
spp_flora[mnemonic %in% c(names(mnemonic_replace), mnemonic_replace)]

check(
    all(mnemonic_replace %in% spp_flora$mnemonic),
    "the codes that replace other codes have a row in the species table"
)

# remove the replaced codes
spp_flora <- spp_flora[!mnemonic %in% names(mnemonic_replace)]

# The species table must now hold exactly the codes of the plot (section 2)
# after the replacement
inventory_mnemonics <- unique(fifelse(
    plot_codes %in% names(mnemonic_replace),
    mnemonic_replace[plot_codes],
    plot_codes
))
setdiff(spp_flora$mnemonic, inventory_mnemonics) # in the table, not in the plot
setdiff(inventory_mnemonics, spp_flora$mnemonic) # in the plot, not in the table
check(
    setequal(spp_flora$mnemonic, inventory_mnemonics),
    "the species table holds exactly the codes present in the plot after the replacement"
)
check(
    !anyDuplicated(spp_flora$mnemonic),
    "one row per mnemonic in the species table"
)
nrow(spp_flora)

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
spp_flora[
    is.na(f_de_vida_r_perez_s_aguilar),
    .(mnemonic, genero_corrected, especie_corrected, f_de_vida_r_foster, f_de_vida_r_perez_s_aguilar)
]

# Growth forms given by the botanists (2026-10-14):
#   ficutr                  strangler (the list gave "árbol grande mayor a 20m")
#   apeihy, nects1, nects3  trees
spp_flora[mnemonic == "ficutr", f_de_vida_r_perez_s_aguilar := "árbol estrangulador mayor a 20m"]
spp_flora[mnemonic %in% c("apeihy", "nects1", "nects3"), f_de_vida_r_perez_s_aguilar := "árbol"]

# R. Foster's code is empty for apeihy, nects1 and nects3 (they are not in his
# list under these codes), and is T (canopy tree) for ficutr; change it by hand
# if the botanists give another value, e.g.
# spp_flora[mnemonic == "ficutr", f_de_vida_r_foster := "TE"]

# Every code has a growth form, except uniden (unidentified stems, left out of
# the species table in section 8)
spp_flora[is.na(f_de_vida_r_perez_s_aguilar), .(mnemonic, genero_corrected, especie_corrected)]
check(
    spp_flora[mnemonic != "uniden", !anyNA(f_de_vida_r_perez_s_aguilar)],
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
    spp_flora[!is.na(f_de_vida_r_perez_s_aguilar), !anyNA(lifeform_class(f_de_vida_r_perez_s_aguilar))],
    "every growth form of R. Pérez & S. Aguilar is read as tree, shrub, palm, tree fern or strangler"
)

# Growth forms by class, and against R. Foster's code
dcast(
    spp_flora[, .N, by = .(class = lifeform_class(f_de_vida_r_perez_s_aguilar), foster = f_de_vida_r_foster)],
    class ~ foster,
    value.var = "N",
    fill = 0
)

# Ficus and other stranglers
spp_flora[
    genero_corrected == "Ficus" | lifeform_class(f_de_vida_r_perez_s_aguilar) %in% "strangler",
    .(mnemonic, genero_corrected, especie_corrected, f_de_vida_r_foster, f_de_vida_r_perez_s_aguilar)
]

# --- 6c. English column names -----------------------------------------------
# Column names: the suffix "_corrected" is dropped, and the Spanish names (from
# "orden" onwards) are translated to English, one to one
colnames(spp_flora) <- gsub("_corrected$", "", colnames(spp_flora))
colnames(spp_flora)

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
    all(names(names_es_en) %in% colnames(spp_flora)),
    "all Spanish column names to translate are in the species table"
)
check(
    !any(unname(names_es_en) %in% colnames(spp_flora)),
    "no English column name is already in the species table"
)

setnames(spp_flora, old = names(names_es_en), new = unname(names_es_en))

colnames(spp_flora)

spp_flora

# =============================================================================
# 7. NOTES ON KNOWN DATA PROBLEMS
# =============================================================================
# Notes exported with the species table (column Notes)
spp_flora[, notes := NA_character_]

# Apeiba "hybrida" (apeihy): one tree, alive in 1982 and 1985, recorded dead
# in 1990.
spp_flora[mnemonic %in% "apeihy", notes := "No longer present: recorded alive until 1985, dead in 1990"]
spp_flora[mnemonic %in% "apeihy"]

# Nectandra s1 (nects1) and s3 (nects3) are Lauraceae morphospecies. nects3
# (one tree) was alive in 1982 and 1985 and recorded dead in 1990; nects1 is
# still alive (one tree in 2022-23) and should be followed to collect it
# fertile and identify it.
spp_flora[mnemonic %in% "nects1"]
spp_flora[mnemonic %in% "nects3"]
spp_flora[mnemonic %in% "nects3", notes := "No longer present: recorded alive until 1985, dead in 1990"]

# pterro includes the trees recorded as pterof (Pterocarpus officinalis does
# not occur in the plot; R. Pérez, email 1). The note gives only the reason:
# the code pterof and its trees are described by the stage-3 species table
# (3_create_species_table_BCI.R: the Notes of pterro, with the number of
# trees), so they are not repeated here.
spp_flora[mnemonic %in% "pterro", notes := "Pterocarpus officinalis does not occur in the plot (R. Pérez, 2026)"]

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
spp_new <- copy(spp_flora[mnemonic != "uniden"])

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
#   author_check  "same"; "format only" (the same authors written another
#                 way: spaces, periods, capital letters, accents, e.g.
#                 "C. DC." and "C.DC."); "different" (other authors); or
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

# Authority reduced to its authors: without spaces, periods, capital letters
# and accents, so that two ways of writing the same authors are equal (the
# rule of every authority comparison below, 8j and 9j, and of stage 3,
# authority_key() in 3_PREPARE_R_TABLES/species_table_functions.R)
authors_only <- function(x) {
    x <- tolower(gsub("[[:space:].]", "", x))
    chartr("áàâäãåéèêëíìîïóòôöõøúùûüñçý", "aaaaaaeeeeiiiioooooouuuuncy", x)
}
wfo_review[, author_check := fcase(
    is.na(wfo_author), "not compared",
    authority == wfo_author, "same",
    authors_only(authority) == authors_only(wfo_author), "format only",
    default = "different"
)]

wfo_review[, family_check := fifelse(differs(family, wfo_family), "differs", "same")]

# Version of the backbone, kept with the review
wfo_review[, wfo_backbone_version := WFO_BACKBONE_VERSION]

check(
    !anyNA(wfo_review$name_check) && !anyNA(wfo_review$author_check),
    "every name checked has a name and an author comparison"
)

# --- 8j. Parts of each name where WFO differs from the Flora de Panamá -------
# For the botanists: one row per name checked whose WFO record differs from the
# Flora de Panamá name in the family, genus, specific epithet, infraspecific
# name (rank and epithet) or authority, with both versions of each part and
# `differing_parts` naming the parts that differ. The authority is compared by
# its authors only: another way of writing the same authors (spaces, periods,
# capital letters, accents, e.g. "C. DC." and "C.DC.") is not a difference;
# where the WFO record has no authority, the authority is not compared. The
# Flora de Panamá values are those before the review of section 9. Written to
# DATA/CHECKS/species_wfo_vs_flora_panama.csv (9j). authors_only() is
# defined in 8i.
# Rank abbreviation of each WFO taxonRank ("variety" -> "var.")
rank_abbreviations <- setNames(names(rank_words), rank_words)

wfo_vs_flora <- wfo_review[, .(
    mnemonic,
    flora_name = paste(scientific_name, authority),
    wfo_name = fifelse(is.na(wfo_author), wfo_name, paste(wfo_name, wfo_author)),
    flora_family = family,
    wfo_family,
    flora_genus = genus,
    wfo_genus,
    flora_species = species,
    wfo_species,
    flora_infraspecific = fifelse(is.na(rank), NA_character_, paste(rank, variety)),
    wfo_infraspecific = fifelse(
        is.na(wfo_variety),
        NA_character_,
        paste(fcoalesce(unname(rank_abbreviations[wfo_rank]), wfo_rank), wfo_variety)
    ),
    flora_authority = authority,
    wfo_authority = wfo_author,
    name_check,
    wfo_status,
    wfo_id
)]
wfo_vs_flora[, differing_parts := gsub("^, |, $", "", gsub("(, )+", ", ", paste(
    fifelse(differs(flora_family, wfo_family), "family", ""),
    fifelse(differs(flora_genus, wfo_genus), "genus", ""),
    fifelse(differs(flora_species, wfo_species), "species", ""),
    fifelse(differs(flora_infraspecific, wfo_infraspecific), "infraspecific", ""),
    fifelse(
        !is.na(wfo_authority) & differs(authors_only(flora_authority), authors_only(wfo_authority)),
        "authority", ""
    ),
    sep = ", "
)))]
wfo_vs_flora[differing_parts == "", differing_parts := NA_character_]
setcolorder(wfo_vs_flora, c("mnemonic", "differing_parts"))
wfo_vs_flora <- wfo_vs_flora[!is.na(differing_parts)]

# Names whose WFO record differs, and the count per part
wfo_vs_flora[, .(mnemonic, differing_parts, flora_name, wfo_name)]
wfo_vs_flora[, .N, by = differing_parts][order(-N)]
check(
    setequal(
        wfo_vs_flora[grepl("family", differing_parts), mnemonic],
        wfo_review[family_check == "differs", mnemonic]
    ) &&
        all(wfo_vs_flora[grepl("authority", differing_parts), mnemonic] %in% wfo_review[author_check == "different", mnemonic]) &&
        all(wfo_review[name_check != "same", mnemonic] %in% wfo_vs_flora$mnemonic),
    "the parts that differ agree with the comparisons of 8i"
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
# The same authors written another way: spaces, periods, capital letters or
# accents (WFO writes initials without spaces, e.g. "C.DC." for "C. DC.").
# Adopt the ones you want written as in WFO.
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
spp_new[mnemonic == "maclti", authority := "(L.) D.Don ex G.Don"]
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
spp_new[mnemonic == "anacex", authority := "Skeels"]
spp_new[mnemonic == "ficuc2", authority := "(Liebm.) Miq."]
spp_new[mnemonic == "heisco", authority := "Standl."]
spp_new[mnemonic == "ouralu", authority := "(Kunth) Engl."]
spp_new[mnemonic == "pse2sp", authority := "(Sw.) Griseb."]
spp_new[mnemonic == "tab1ro", authority := "(Bertol.) DC."]
spp_new[mnemonic == "thevah", authority := "A.DC."]
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

# Where the species table (after the review) still differs from its WFO
# record: one row per name checked that differs in family, genus, species,
# infraspecific name or authority ("authority" when the authors differ,
# "authority format" when only the writing differs: spaces, periods, capital
# letters, accents), or whose WFO record is not an accepted name (status and
# accepted name shown). The morphospecies are not checked against WFO and are
# not listed. Written to DATA/CHECKS/species_table_vs_wfo.csv.
table_vs_wfo <- merge(
    spp_new[wfo_check == TRUE, .(
        mnemonic,
        table_name = paste(scientific_name, authority),
        table_family = family,
        table_genus = genus,
        table_species = species,
        table_infraspecific = fifelse(is.na(rank), NA_character_, paste(rank, variety)),
        table_authority = authority
    )],
    wfo_review[, .(
        mnemonic,
        wfo_name = fifelse(is.na(wfo_author), wfo_name, paste(wfo_name, wfo_author)),
        wfo_family,
        wfo_genus,
        wfo_species,
        wfo_infraspecific = fifelse(
            is.na(wfo_variety),
            NA_character_,
            paste(fcoalesce(unname(rank_abbreviations[wfo_rank]), wfo_rank), wfo_variety)
        ),
        wfo_authority = wfo_author,
        wfo_status,
        wfo_accepted_name = fifelse(
            is.na(wfo_accepted_author), wfo_accepted_name, paste(wfo_accepted_name, wfo_accepted_author)
        ),
        wfo_id
    )],
    by = "mnemonic"
)
check(
    nrow(table_vs_wfo) == nrow(wfo_review),
    "every name checked against WFO is compared with the species table"
)
table_vs_wfo[, differing_parts := gsub("^, |, $", "", gsub("(, )+", ", ", paste(
    fifelse(differs(table_family, wfo_family), "family", ""),
    fifelse(differs(table_genus, wfo_genus), "genus", ""),
    fifelse(differs(table_species, wfo_species), "species", ""),
    fifelse(differs(table_infraspecific, wfo_infraspecific), "infraspecific", ""),
    fcase(
        is.na(wfo_authority), "",
        differs(authors_only(table_authority), authors_only(wfo_authority)), "authority",
        differs(table_authority, wfo_authority), "authority format",
        default = ""
    ),
    sep = ", "
)))]
table_vs_wfo[differing_parts == "", differing_parts := NA_character_]
table_vs_wfo[wfo_status %in% "Accepted", wfo_accepted_name := NA_character_]
setcolorder(table_vs_wfo, c("mnemonic", "differing_parts", "wfo_status"))
table_vs_wfo <- table_vs_wfo[!is.na(differing_parts) | differs(wfo_status, "Accepted")]

# Names of the species table that differ from WFO (with the current inputs:
# none by any part; 4 names are synonyms in WFO: guargr, laetth, swars1, swars2)
table_vs_wfo[, .N, by = .(differing_parts, wfo_status)]
table_vs_wfo[, .(mnemonic, differing_parts, wfo_status, table_name, wfo_name, wfo_accepted_name)]

# Save the review, the candidate records and the changes made
checks_folder <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS")
fwrite(wfo_review, file.path(checks_folder, "species_wfo_review.csv"))
fwrite(
    wfo_candidates[n_records > 1, .(mnemonic, spec.name, Authorship, taxonID, scientificName, scientificNameAuthorship, taxonomicStatus, Auth.dist, kept)],
    file.path(checks_folder, "species_wfo_candidates.csv")
)
fwrite(accepted_changes, file.path(checks_folder, "species_wfo_accepted_changes.csv"))
# For the botanists (8j); with a byte order mark so that Excel shows the
# accents of the authorities (Müll.Arg.) correctly
fwrite(wfo_vs_flora, file.path(checks_folder, "species_wfo_vs_flora_panama.csv"), bom = TRUE)
fwrite(table_vs_wfo, file.path(checks_folder, "species_table_vs_wfo.csv"), bom = TRUE)

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

# --- 10d. Synonyms compared with WFO (information only) -----------------------
# Each synonym of the table (column synonyms, names separated by commas, without
# authority) looked up in the WFO backbone by its full name. A synonym should
# be a name of the same taxon as the code. synonym_check:
#   "same taxon in WFO"        the synonym and the code's name lead to the same
#                              accepted WFO record
#   "same species in WFO"      the synonym leads to a subspecies or variety of
#                              the species WFO accepts for the code (e.g.
#                              Bactris augustinea -> Bactris major var. major)
#   "another taxon in WFO"     every WFO record of the synonym leads to another
#                              accepted name: WFO places it in another species,
#                              or accepts it as a species of its own
#   "not found in WFO"         no WFO record carries this name (closest WFO name
#                              of the same genus in closest_wfo_name)
#   "repeats the name"         the synonym is the code's own name (10c)
#   "code not checked"         the code is a morphospecies (not checked, 8b)
# Nothing is changed. Written to DATA/CHECKS/species_synonyms_vs_wfo.csv.

# Accepted record of every WFO record (its own taxonID when accepted)
wfo_accepted_of <- wfo_backbone[, .(
    taxonID,
    accepted_id = fifelse(taxonomicStatus %in% "Accepted" | acceptedNameUsageID %in% c(NA, ""), taxonID, acceptedNameUsageID)
)]
wfo_accepted_of[wfo_backbone, on = .(accepted_id = taxonID), `:=`(
    accepted_name = i.scientificName,
    accepted_author = i.scientificNameAuthorship
)]

# Accepted record of each code (from the record kept in 8h), and its species
code_accepted <- merge(
    wfo_review[, .(mnemonic, wfo_id)],
    wfo_accepted_of[, .(wfo_id = taxonID, code_accepted_id = accepted_id, code_accepted_name = accepted_name)],
    by = "wfo_id",
    all.x = TRUE
)
code_accepted[, code_accepted_species := word(code_accepted_name, 1, 2)]

synonyms_long <- spp_new[!is.na(synonyms), .(synonym = str_squish(unlist(strsplit(synonyms, ",")))), by = .(mnemonic, scientific_name)]
synonyms_long <- synonyms_long[synonym != ""]
synonym_records <- merge(
    synonyms_long,
    wfo_backbone[, .(synonym = scientificName, taxonID, taxonomicStatus)],
    by = "synonym",
    all.x = TRUE,
    allow.cartesian = TRUE
)
synonym_records[wfo_accepted_of, on = "taxonID", `:=`(accepted_id = i.accepted_id, accepted_name = i.accepted_name, accepted_author = i.accepted_author)]
synonym_records[code_accepted, on = "mnemonic", `:=`(code_accepted_id = i.code_accepted_id, code_accepted_species = i.code_accepted_species)]

synonyms_vs_wfo <- synonym_records[, .(
    wfo_records = sum(!is.na(taxonID)),
    same_taxon = any(accepted_id == code_accepted_id, na.rm = TRUE),
    same_species = any(word(accepted_name, 1, 2) == code_accepted_species, na.rm = TRUE),
    checked_code = !is.na(code_accepted_id[1]),
    wfo_status = paste(sort(unique(na.omit(taxonomicStatus))), collapse = "/"),
    wfo_accepted_name = paste(unique(na.omit(str_squish(paste(accepted_name, fcoalesce(accepted_author, ""))))), collapse = "; ")
), by = .(mnemonic, scientific_name, synonym)]
synonyms_vs_wfo[, synonym_check := fcase(
    synonym == scientific_name, "repeats the name",
    !checked_code, "code not checked",
    wfo_records == 0L, "not found in WFO",
    same_taxon, "same taxon in WFO",
    same_species, "same species in WFO",
    default = "another taxon in WFO"
)]

# Names not found: the closest WFO name of the same genus (at most 2 letters
# apart), e.g. another gender ending of the epithet
not_found <- synonyms_vs_wfo[synonym_check == "not found in WFO", unique(synonym)]
closest <- vapply(not_found, function(s) {
    pool <- wfo_backbone[genus == word(s, 1) & taxonRank %in% "species", unique(scientificName)]
    if (!length(pool)) {
        return(NA_character_)
    }
    d <- drop(adist(s, pool))
    if (min(d) <= 2L) pool[which.min(d)] else NA_character_
}, character(1))
synonyms_vs_wfo[, closest_wfo_name := unname(closest[synonym])]
synonyms_vs_wfo[, c("same_taxon", "same_species", "checked_code") := NULL]
setcolorder(synonyms_vs_wfo, c("mnemonic", "scientific_name", "synonym", "synonym_check"))
setorder(synonyms_vs_wfo, mnemonic, synonym)
check(
    nrow(synonyms_vs_wfo) == nrow(synonyms_long) && !anyNA(synonyms_vs_wfo$synonym_check),
    "every synonym of the table is compared with WFO once"
)

# With the current inputs: 203 synonyms of 157 codes; the ones that are not the same taxon
# or species in WFO are printed for review
synonyms_vs_wfo[, .N, by = synonym_check][order(-N)]
synonyms_vs_wfo[
    !synonym_check %in% c("same taxon in WFO", "same species in WFO"),
    .(mnemonic, scientific_name, synonym, synonym_check, wfo_accepted_name, closest_wfo_name)
]
fwrite(synonyms_vs_wfo, file.path(checks_folder, "species_synonyms_vs_wfo.csv"), bom = TRUE)
rm(wfo_accepted_of, code_accepted, synonyms_long, synonym_records, not_found, closest)

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
    setequal(full_out$Mnemonic, setdiff(inventory_mnemonics, "uniden")),
    "the exported table holds every code of the plot except uniden"
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
