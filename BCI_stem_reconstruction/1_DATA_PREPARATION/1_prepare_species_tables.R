# =============================================================================
# 1_prepare_species_tables.R
#
# Purpose: Prepare the BCI 50-ha plot species table, one row per mnemonic of
#          the plot: the database taxonomy, corrected with the emails of
#          R. Pérez and S. Lao and completed from the Flora de Panamá list,
#          which has priority. The names are then checked against TNRS (World
#          Flora Online) step by step: TNRS only suggests corrections, mainly
#          spelling and authority format, and each one is reviewed and
#          accepted by hand; nothing is changed automatically.
#
# Run from the project root (or from BCI_stem_reconstruction/), line by line.
# The TNRS calls need an internet connection.
#
# Inputs (BCI_stem_reconstruction/DATA/RAW/):
#   oct_13_2026/ViewTaxonomy_bci.rds               database taxonomy (from
#                                                  0_extract_ViewTaxonomy.R)
#   oct_13_2026/ViewFullTable_bci.csv              current inventory (codes in the plot)
#   ViewFiles_bci_allcensuses/ViewFullTable_bci.csv  previous inventory (taxonomy comparison)
#   sp_tables/email_1/Lista_bci_mnemonics_formadevida.xlsx  list of R. Pérez & S. Aguilar
#   sp_tables/email_2/bci_name_changes.xlsx        name changes of S. Lao
#   sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx
# Outputs:
#   DATA/SPP_TABLE/bci_spptable.txt, bci_spptable.csv (both tab-separated) and
#   bci_spptable.RData (object bci.spptable): one row per code of the plot
#   except uniden, with Mnemonic, Order, Family, Genus, SpeciesName,
#   InfraspecificRank, InfraspecificEpithet, Authority, Synonyms,
#   Lifeform_RFoster, Lifeform_RPerez_SAguilar, Lifeform_FloraPanama,
#   CommonName, Herbarium, Notes.
#   DATA/CHECKS/species_tnrs_results_raw.rds (TNRS answer),
#   species_tnrs_review.csv, species_tnrs_accepted_changes.csv,
#   prepare_species_tables_log.txt (failed checks).
# Packages: data.table, TNRS, stringr, compareDF, janitor, readxl.
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
# LOAD REQUIRED PACKAGES
# =============================================================================
library(data.table) # Fast table manipulation
library(TNRS) # Taxonomic Name Resolution Service
library(stringr) # String helpers (trimming, case, etc.)
library(compareDF) # Compare data frames and generate difference tables
library(janitor) # Clean column names, check duplicates, etc.

# Number of threads data.table will use
data.table::getDTthreads()

# =============================================================================
# DEFINE INPUT FILE PATHS
# =============================================================================
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
    workspace_root <- dirname(workspace_root)
}

# load raw Mnemonics

ViewTaxonomy <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026", "ViewTaxonomy_bci.rds")
TAXONOMY_1 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_1", "Lista_bci_mnemonics_formadevida.xlsx")
TAXONOMY_2 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "email_2", "bci_name_changes.xlsx")
flora_panama_oct_9_2026 <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "flora_panama_oct_9_2026", "Lista de especies de la Flora de Panamá_nov25 original.xlsx")
ViewFullTable_old <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses", "ViewFullTable_bci.csv")
# current inventory (the table 0_extract_ViewTaxonomy.R reads)
ViewFullTable_new <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "oct_13_2026", "ViewFullTable_bci.csv")

# =============================================================================
# 1. LOAD ALL DATA
# =============================================================================

viewt <- as.data.table(readRDS(ViewTaxonomy))
tax1 <- as.data.table(readxl::read_xlsx(TAXONOMY_1))
tax2 <- as.data.table(readxl::read_xlsx(TAXONOMY_2))
flor_pan <- as.data.table(readxl::read_xlsx(flora_panama_oct_9_2026))
viewftold <- as.data.table(fread(ViewFullTable_old))

inc <- colnames(viewt)

viewtold <- unique(viewftold[, ..inc])
setorder(viewtold, Mnemonic)

# =============================================================================
# 2. COMPARE OLD vs. NEW taxonomy in the viewfulltable datasets
# =============================================================================

viewt[, .N, by = Mnemonic][order(-N)]

viewtold[, .N, by = Mnemonic][order(-N)]
viewtold[, RowID := .I]
viewtold[Mnemonic == "protte"]
viewtold[Mnemonic == "quaras"]
viewtold <- viewtold[!RowID %in% c(236, 257)]
viewtold[, RowID := NULL]

viewt[, .N, by = Mnemonic][order(-N)]
viewtold[, .N, by = Mnemonic][order(-N)]

viewt <- viewt[, .(Mnemonic, Family, Genus, SpeciesName, Subspecies)]
viewtold <- viewtold[, .(Mnemonic, Family, Genus, SpeciesName, Subspecies)]

# is NA or is NULL, make NA
viewt[viewt == "NULL"] <- NA
viewtold[viewtold == "NULL"] <- NA

# if character "", make NA
viewt[viewt == ""] <- NA
viewtold[viewtold == ""] <- NA

unique(viewt$Mnemonic == viewtold$Mnemonic)

ctable_mnemonic <- compare_df(viewt, viewtold, c("Mnemonic"))
ctable_mnemonic$comparison_df

# create_output_table(ctable_mnemonic)

# create_output_table(ctable_mnemonic, output_type = "xlsx", file_name = "./BCI_stem_reconstruction/1_DATA_PREPARATION/compare_2025_vs_oct2026_taxonomy_from_viewfulltables.xlsx")

# viewt[Mnemonic == "appuse"]
# viewt[Mnemonic == "clidse"]

rm(viewtold)

# =============================================================================
# 2. ASSESS FIX 1
# =============================================================================

# NOTE: The following email from Rolando Pérez (Sep 28, 2026) provides context for the species list and known corrections. It is included here for reference.

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

## So far, I have 329 mnemonic codes. Rolando provided 327. Lets see which ones i have extra
setdiff(viewt$Mnemonic, tax1$codigo)

## Lets merge my table with Rolando's table to proceed with the checking process
viewt_tax1 <- merge(
    viewt,
    tax1,
    by.x = "Mnemonic",
    by.y = "codigo",
    all.x = TRUE,
    all.y = TRUE
)

## Lets assess comments one by one. I will create a new column called "notes" to store the comments and observations.


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
viewt_tax1[Genus == "Beilschmiedia" & SpeciesName == "tovarensis", email_1_comment := "Beilschmiedia pendula no existe para BCI; la especie correcta es Beilschmiedia tovarensis"]
viewt_tax1[Genus == "Quararibea" & SpeciesName == "stenophylla", email_1_comment := "Quararibea asterolepis no existe para BCI; la especie correcta es Quararibea stenophylla"]

viewt_tax1[!is.na(email_1_comment), .(Mnemonic, Family, Genus, SpeciesName, email_1_comment)]

# Añadimos las formas de vida que utilizo Robin Foster con pequeñas modificaciones realizadas por
# R. Pérez y S. Aguilar. También está la forma de vida que Salomón y Yo le hemos asignado a las
# especies de la parcela en español.
# Los nombres comunes de las especies y nuestras colecciones de herbario.

## check differences in family, genus, species, subspecies, and order

viewt_tax1[(is.na(familia) & !is.na(Family)) | Family != familia, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(familia), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(familia), familia := Family]

viewt_tax1[(is.na(genero) & !is.na(Genus)) | Genus != genero, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(genero), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(genero), genero := Genus]

viewt_tax1[(is.na(especie) & !is.na(SpeciesName)) | SpeciesName != especie, .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(especie), .(Mnemonic, Family, familia, Genus, genero, SpeciesName, especie, email_1_comment)]
viewt_tax1[is.na(especie), especie := SpeciesName]

viewt_tax1[!is.na(Subspecies)]

viewt_tax1[, Family := NULL]
viewt_tax1[, Genus := NULL]
viewt_tax1[, SpeciesName := NULL]

# FINISHED: Rolando email 1 included

# =============================================================================
# 2. ASSESS FIX 2
# =============================================================================

# NOTE: The following emails from Suzanne Lao (July 28 , 2026) provides corrections.
# NOTE: Rolando Perez responded to the email making an additional correction. The email is included here for reference.

# First Suzanne email

# From: Lao, Suzanne <LAOZ@si.edu>
# Sent: Tuesday, July 28, 2026 2:28 PM
# To: Medina-Vega, Jose A. <MedinaJA@si.edu>
# Cc: Davies, Stuart J. <DaviesS@si.edu>; Perez, Rolando <PerezRo@si.edu>
# Subject: BCI species name changes

# Hola José,

# Adjunto le envío el listado de los nombres de especies de la parcela de BCI. Las primeras
# columnas son los nombres nuevos y las columnas al final son los nombres viejos.

# Solamente un código cambió: todos los "ficutr" cambiaron a "ficuar".
# NOTE: ficur is wrong, Suzanne ment ficuau, not ficuar. There is no ficuar.
# Me avisa si este listado es diferente al que le envió Rolando. Es posible que el listado de
# Rolando esté más actualizado.

# Gracias y saludos,
# Suzanne

## So far, I have 329 mnemonic codes. Suzanne provided 330, Let's investigate the difference
# i have one per mnemonic
viewt_tax1[, .N, by = Mnemonic][order(-N)]

# Suzanne has 2 per mnemonic - the problematic protte
tax2[, .N, by = mnemonic][order(-N)]
tax2[mnemonic == "protte"]
# lets exclude it to make a clean join
tax2[, RowID := .I]
tax2[mnemonic == "protte"]
tax2 <- tax2[!RowID %in% c(267)]
tax2[, RowID := NULL]

# lets count again
# i have 1 per mnemonic
setorder(viewt_tax1, Mnemonic)
viewt_tax1[, .N, by = Mnemonic][order(-N)]

# suzame has 1 per mnemonic
setorder(tax2, mnemonic)
tax2[, .N, by = mnemonic][order(-N)]

## lets find the differences
# i have "ficutr" and Suzanne has "ficuar"
setdiff(viewt_tax1$Mnemonic, tax2$mnemonic)
setdiff(tax2$mnemonic, viewt_tax1$Mnemonic)

viewt_tax1[Mnemonic == "ficutr"]
# viewt_tax1[Mnemonic == "ficuau"]
# tax2[mnemonic == "ficutr"]
tax2[mnemonic == "ficuau"]

tax2[, speciesid := NULL]
tax2[, subspeciesid := NULL]
tax2[, speciesidold := NULL]
tax2[, subspeciesidold := NULL]

# Second, Rolands email correcting Suzanne email

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

## Make corrections Rolando suggested

# first explore whats going on
tax2[mnemonic %in% c("stemgr", "pla1pi")]

# viewt_tax1 is already corrected, however, it only missed the sinonimo for pla1pi
viewt_tax1[Mnemonic %in% c("stemgr", "pla1pi")]

viewt_tax1[Mnemonic %in% c("pla1pi")]
viewt_tax1[Mnemonic %in% c("pla1pi"), sinonimos := "Platymiscium polystachyum"]

## check
viewt_tax1[Mnemonic %in% c("stemgr", "pla1pi"), .(Mnemonic, genero, especie, autoridad, sinonimos)]

# FINISHED: Rolando email 2 included

tax2
setorder(tax2, mnemonic)

viewt_tax1
setorder(viewt_tax1, Mnemonic)

susfam <- tax2[, .(mnemonic, family, familyold)]
mefam <- viewt_tax1[, .(Mnemonic, familia)]

mefam_susfam <- merge(
    mefam,
    susfam,
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)

mefam_susfam[familia != family]
## There are several inconsistencies in the family names between my table and Suzanne's table. I will keep my family names as they are correct.
## ones provided by Rolando. Also, Suzanne agrees in the email that Rolando's are more up to date.

susgen <- tax2[, .(mnemonic, genus, genusold)]
megen <- viewt_tax1[, .(Mnemonic, genero)]

megen_susgen <- merge(
    megen,
    susgen,
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)

megen_susgen[genero != genus]

## Same with genera

# Now with species
sussp <- tax2[, .(mnemonic, speciesname, speciesnameold)]
mesp <- viewt_tax1[, .(Mnemonic, especie)]

mesp_sussp <- merge(
    mesp,
    sussp,
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)

mesp_sussp[especie != speciesname]

## same with specues

# now with subspecies

sussubsp <- tax2[, .(mnemonic, subspecies, subspeciesold)]
mesubsp <- viewt_tax1[, .(Mnemonic, Subspecies)]

mesubsp_sussubsp <- merge(
    mesubsp,
    sussubsp,
    by.x = "Mnemonic",
    by.y = "mnemonic",
    all.x = TRUE,
    all.y = TRUE
)

mesubsp_sussubsp[Subspecies != subspecies]

# we are good here.

## so far, viewt_tax1 is correct

## now, whats missing regarding the mnemonic codes
setdiff(viewt_tax1$Mnemonic, tax2$mnemonic)
setdiff(tax2$mnemonic, viewt_tax1$Mnemonic)
##

## We are back to email 2, were Suzanne is telling us that all "ficutr" should be changed to "ficuar".
# NOTE: ficuar is wrong. It should be ficuau

# lets check again what does viewt_tax1 has for "ficutr" and "ficuau"
viewt_tax1[Mnemonic == "ficutr"] # it has ficus trigonata
viewt_tax1[Mnemonic == "ficuau"] # it does not have ficus aurea

# why viewt_tax1 has "ficutr" and not "ficuau"? Lets check the original viewt table
viewt[Mnemonic == "ficutr"] # there is ficus aurea in the plot, but the mnemonic is wrong. It should be ficuau, not ficutr
viewt[Mnemonic == "ficuau"]

# the original viewt has ficutr with name ficus aurea, which, logically, should be ficuau
# lets see wether it has ficus trigonata
viewt_tax1[Mnemonic == "ficutr"]

viewt[Genus == "Ficus"][order(SpeciesName)]
# no, the plot does not have ficus trigonata, so the mnemonic ficutr is wrong.
# It should be ficuar However, if I exclude it from viewt_tax1, I will have a
# missing mnemonic code. So, I will keep it in the table, but I will add a
# comment to the email_2_comment column.

setdiff(viewt_tax1$Mnemonic, tax2$mnemonic)
setdiff(tax2$mnemonic, viewt_tax1$Mnemonic)

tax2[mnemonic == "ficuau"]

# Include ficuar mnemonic to viewt_tax1
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

## lets check wether all species are included in the viwefulltables

setdiff(viewt_tax1_updated$Mnemonic, viewt$Mnemonic)

# FINISHED: Rolando email 1 and 2, Suzanne email 1 included

# NOTE: ficutr is replaced with ficuau in the BCI inventory in
# 2_prepare_viewfulltable.R; the ficutr row is removed from this table in section
# "4. MNEMONICS REPLACED IN THE INVENTORY".

# =============================================================================
# 3. ASSESS FIX 3
# =============================================================================

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

# clean columnames for flor_pan
flor_pan <- clean_names(flor_pan)
viewt_tax1_updated <- clean_names(viewt_tax1_updated)

# add the prefix "pan" to the column names of flor_pan to avoid name clashes
colnames(flor_pan) <- paste0("pan_", colnames(flor_pan))

# make all NA as NA in viewt_tax1_updated
viewt_tax1_updated <- viewt_tax1_updated[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]
flor_pan <- flor_pan[, lapply(.SD, function(x) ifelse(x == "NULL", NA, x))]

viewt_tax1_updated_flor_pan <- merge(
    viewt_tax1_updated,
    flor_pan,
    by.x = "mnemonic",
    by.y = "pan_codigo",
    all.x = TRUE,
    all.y = FALSE
)

viewt_tax1_updated_flor_pan[, texto := NULL][, fotos := NULL][, pan_texto := NULL][, pan_fotos := NULL][, pan_mnbc := NULL]

# check wether all mnemonic in viewt_tax1_updated_flor_pan are present in the viewfulltable

setdiff(viewt_tax1_updated_flor_pan$mnemonic, viewt$Mnemonic)
setdiff(viewt$Mnemonic, viewt_tax1_updated_flor_pan$mnemonic)

# yes, the new one doesnt have ficuau, its because I include it and this one
# needs to be replaced in the BCI inventory (raw data). I will do that in the
# next script.

## now, lets compare one by one and accept the most updated one

# orden
viewt_tax1_updated_flor_pan[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden)]
## create orden_corrected, when they are different or one is not available, I will use the pan_orden, otherwise I will use orden
viewt_tax1_updated_flor_pan[, orden_corrected := fcoalesce(pan_orden, orden)]
viewt_tax1_updated_flor_pan[(is.na(orden) | is.na(pan_orden)) | orden != pan_orden, .(mnemonic, orden, pan_orden, orden_corrected)]
viewt_tax1_updated_flor_pan[, orden := NULL][, pan_orden := NULL]

# familia
viewt_tax1_updated_flor_pan[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia)]
# create familia_corrected, when they are different or one is not available, I will use the pan_familia, otherwise I will use familia
viewt_tax1_updated_flor_pan[, familia_corrected := fcoalesce(pan_familia, familia)]
viewt_tax1_updated_flor_pan[(is.na(familia) | is.na(pan_familia)) | familia != pan_familia, .(mnemonic, familia, pan_familia, familia_corrected)]
viewt_tax1_updated_flor_pan[, familia := NULL][, pan_familia := NULL]

# genero
viewt_tax1_updated_flor_pan[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero)]
# create genero_corrected, when they are different or one is not available, I will use the pan_genero, otherwise I will use genero
viewt_tax1_updated_flor_pan[, genero_corrected := fcoalesce(pan_genero, genero)]
viewt_tax1_updated_flor_pan[(is.na(genero) | is.na(pan_genero)) | genero != pan_genero, .(mnemonic, genero, pan_genero, genero_corrected)]
viewt_tax1_updated_flor_pan[, genero := NULL][, pan_genero := NULL]

# especies
viewt_tax1_updated_flor_pan[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie)]
# create especie_corrected, when they are different or one is not available, I will use the pan_especie, otherwise I will use especie
viewt_tax1_updated_flor_pan[, especie_corrected := fcoalesce(pan_especie, especie)]
viewt_tax1_updated_flor_pan[(is.na(especie) | is.na(pan_especie)) | especie != pan_especie, .(mnemonic, especie, pan_especie, especie_corrected)]
viewt_tax1_updated_flor_pan[, especie := NULL][, pan_especie := NULL]

# subspecies is already included in especies_corrected. Lets remove it from the table
viewt_tax1_updated_flor_pan[!is.na(subspecies)]
viewt_tax1_updated_flor_pan[, subspecies := NULL]

# autoridad
viewt_tax1_updated_flor_pan[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad)]
# create autoridad_corrected, when they are different or one is not available, I will use the pan_autoridad, otherwise I will use autoridad
viewt_tax1_updated_flor_pan[, autoridad_corrected := fcoalesce(pan_autoridad, autoridad)]
viewt_tax1_updated_flor_pan[(is.na(autoridad) | is.na(pan_autoridad)) | autoridad != pan_autoridad, .(mnemonic, autoridad, pan_autoridad, autoridad_corrected)]
viewt_tax1_updated_flor_pan[, autoridad := NULL][, pan_autoridad := NULL]

# sinonimos
viewt_tax1_updated_flor_pan[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos)]
# create sinonimos_corrected, when they are different or one is not available, I will use the pan_sinonimos, otherwise I will use sinonimos
viewt_tax1_updated_flor_pan[, sinonimos_corrected := fcoalesce(pan_sinonimos, sinonimos)]
viewt_tax1_updated_flor_pan[(is.na(sinonimos) | is.na(pan_sinonimos)) | sinonimos != pan_sinonimos, .(mnemonic, sinonimos, pan_sinonimos, sinonimos_corrected)]
viewt_tax1_updated_flor_pan[, sinonimos := NULL][, pan_sinonimos := NULL]

# nombre comun
viewt_tax1_updated_flor_pan[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum)]
# create nombre_comun_corrected, when they are different or one is not available, I will use the pan_nombre_comum, otherwise I will use nombre_comum
viewt_tax1_updated_flor_pan[, nombre_comun_corrected := fcoalesce(pan_nombre_comum, nombre_comum)]
viewt_tax1_updated_flor_pan[(is.na(nombre_comum) | is.na(pan_nombre_comum)) | nombre_comum != pan_nombre_comum, .(mnemonic, nombre_comum, pan_nombre_comum, nombre_comun_corrected)]
viewt_tax1_updated_flor_pan[, nombre_comum := NULL][, pan_nombre_comum := NULL]

# herbario
viewt_tax1_updated_flor_pan[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario)]
# create herbario_corrected, when they are different or one is not available, I will use the pan_herbario, otherwise I will use herbario
viewt_tax1_updated_flor_pan[, herbario_corrected := fcoalesce(pan_herbario, herbario)]
viewt_tax1_updated_flor_pan[(is.na(herbario) | is.na(pan_herbario)) | herbario != pan_herbario, .(mnemonic, herbario, pan_herbario, herbario_corrected)]
viewt_tax1_updated_flor_pan[, herbario := NULL][, pan_herbario := NULL]

viewt_tax1_updated_flor_pan[mnemonic == "appuse"]

viewt_tax1_updated_flor_pan[mnemonic == "ficutr"]

viewt_tax1_updated_flor_pan[mnemonic == "ficuau"]

viewt_tax1_updated_flor_pan[!is.na(email_1_comment)]

# =============================================================================
# 4. MNEMONICS REPLACED IN THE INVENTORY
# =============================================================================
# Two codes of the inventory are replaced by another code. The same replacement
# is made in 2_prepare_viewfulltable.R, before the species table is merged:
#   ficutr -> ficuau  the trees recorded as ficutr are Ficus aurea, whose code
#                     is ficuau (Suzanne, email 2)
#   pterof -> pterro  all Pterocarpus in the plot are P. rohrii (Rolando, email 1)
# After the replacement ficutr and pterof have no records in the plot, so their
# rows are removed from the species table here.
mnemonic_replace <- c(ficutr = "ficuau", pterof = "pterro")

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

# rows of the replaced codes and of the codes that replace them
viewt_tax1_updated_flor_pan[mnemonic %in% c(names(mnemonic_replace), mnemonic_replace)]

check(
    all(mnemonic_replace %in% viewt_tax1_updated_flor_pan$mnemonic),
    "the codes that replace ficutr and pterof have a row in the species table"
)

# remove the replaced codes
viewt_tax1_updated_flor_pan <- viewt_tax1_updated_flor_pan[!mnemonic %in% names(mnemonic_replace)]

# the species table must now hold exactly the codes of the inventory after the
# replacement. The codes come from the current inventory (oct_13_2026), the
# same table 0_extract_ViewTaxonomy.R reads; only the Mnemonic column is loaded.
inventory_mnemonics <- unique(fread(ViewFullTable_new, sep = "\t", select = "Mnemonic")$Mnemonic)
inventory_mnemonics <- unique(fifelse(
    inventory_mnemonics %in% names(mnemonic_replace),
    mnemonic_replace[inventory_mnemonics],
    inventory_mnemonics
))
setdiff(viewt_tax1_updated_flor_pan$mnemonic, inventory_mnemonics) # in the table, not in the plot
setdiff(inventory_mnemonics, viewt_tax1_updated_flor_pan$mnemonic) # in the plot, not in the table
check(
    setequal(viewt_tax1_updated_flor_pan$mnemonic, inventory_mnemonics),
    "the species table holds exactly the codes present in the plot after the replacement"
)
check(
    !anyDuplicated(viewt_tax1_updated_flor_pan$mnemonic),
    "one row per mnemonic in the species table"
)
nrow(viewt_tax1_updated_flor_pan)

viewt_tax1_updated_flor_pan

## remove suffix "_corrected" from the column names
colnames(viewt_tax1_updated_flor_pan) <- gsub("_corrected$", "", colnames(viewt_tax1_updated_flor_pan))

colnames(viewt_tax1_updated_flor_pan)

## translate the Spanish column names (from "orden" onwards) to English, one to one
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
# 1.1. CHECK RAW TAXONOMY DATA
# =============================================================================
# Known data problems in the new taxonomy table: the lines below print the
# rows concerned and record a note (column `notes`) where one applies.
viewt_tax1_updated_flor_pan[, notes := NA_character_]

# Apeiba "hybrida" is no longer present in the plot.
viewt_tax1_updated_flor_pan[mnemonic %in% "apeihy", notes := "No longer present; last seen 1990"]
viewt_tax1_updated_flor_pan[mnemonic %in% "apeihy"]

# Nectandra s1 and Nectandra s3 are Lauraceae morphospecies; one of them died.
# Only one remains alive, and it should be monitored for fertility-based identification.
viewt_tax1_updated_flor_pan[mnemonic %in% "nects1"]
viewt_tax1_updated_flor_pan[mnemonic %in% "nects3"]

viewt_tax1_updated_flor_pan[mnemonic %in% "nects3", notes := "No longer present; last seen 1990"]

# Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
viewt_tax1_updated_flor_pan[mnemonic %in% "pterro", notes := "Previously identified as 'pterof' (Pterocarpus officinalis) are incorrect; those are 'pterro' (Pterocarpus rohrii)"]

# =============================================================================
# 2. PREPARE NEW TAXONOMY FOR TNRS VALIDATION
# =============================================================================
# Split every name into its parts and mark which rows hold a name that a name
# service can validate. The names in this table are those of the Flora de
# Panamá list and stay as the reference: TNRS / WorldFlora results are added
# later in their own columns and are used only to detect spelling mistakes.
#
# Columns created here:
#   species_panama   species column exactly as it came from the merge
#                    (e.g. "simplex var. grandiflora"), kept for traceability
#   species          specific epithet only ("simplex")
#   rank             infraspecific rank ("var.", "subsp.", "f."); NA for species
#   variety          infraspecific name ("grandiflora"); NA for species
#   name_type        "species", "infraspecific" or "morphospecies"
#   tnrs_check       TRUE for the rows whose name can be sent to a name service
#   scientific_name  full name without authority, built from the parts above
spp_new <- copy(viewt_tax1_updated_flor_pan[mnemonic != "uniden"])

# Genus must be capitalised for TNRS.
spp_new[, genus := str_to_sentence(str_squish(genus))]

# Species epithets should be lower-case, without repeated spaces
# ("simplex  var. grandiflora" has two spaces in the Flora de Panamá list).
spp_new[, species := str_to_lower(str_squish(species))]

# --- 2a. Separate species, infraspecific rank and infraspecific name ---------
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

# --- 2a.1 Classify the kind of name ------------------------------------------
# morphospecies: field names that are not published names ("sp.", "sp.4_...",
#                or a label such as "hybrid" / "fuzzy"). They are recognised by
#                the "sp." prefix or by having no authority, and cannot be
#                validated by a name service.
# infraspecific: species + rank + infraspecific name
# species:       genus + a single specific epithet
spp_new[, name_type := fcase(
    grepl("^sp\\.", species) | is.na(authority), "morphospecies",
    !is.na(rank), "infraspecific",
    grepl("^[a-z-]+$", species), "species"
)]
spp_new[, tnrs_check := name_type %in% c("species", "infraspecific")]

spp_new[, .N, by = .(name_type, tnrs_check)]
# morphospecies: set aside, not sent to the name service
spp_new[name_type == "morphospecies", .(mnemonic, family, genus, species, authority, notes)]

# Full name without authority, built from the separated parts. This is the
# Flora de Panamá name that the name-service result is compared with.
spp_new[, scientific_name := fifelse(
    name_type %in% "infraspecific",
    paste(genus, species, rank, variety),
    paste(genus, species)
)]

# --- 2a.2 Checks on the separated names --------------------------------------
check(
    !anyNA(spp_new$name_type),
    "every row is classified as species, infraspecific or morphospecies"
)
check(
    spp_new[tnrs_check == TRUE, all(grepl("^[A-Z][a-z-]+$", genus) & grepl("^[a-z-]+$", species))],
    "names to validate have a one-word genus and a one-word specific epithet"
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
    !anyDuplicated(spp_new[tnrs_check == TRUE, scientific_name]),
    "names to validate are unique (no two codes share the same full name)"
)

# --- 2a.3 Authority ----------------------------------------------------------
# Strip non-authority annotations (sensu, auct., nom. dub., ined.) from the
# authority field; they confuse the TNRS parser
spp_new[, authority := stringr::str_replace_all(
    authority,
    regex("sensu.*|auct\\.|nom\\. dub\\.|ined\\.", ignore_case = TRUE),
    ""
)]
spp_new[, authority := str_squish(authority)]
spp_new[authority == "", authority := NA_character_]
check(
    spp_new[tnrs_check == TRUE, !anyNA(authority)],
    "every name to validate still has an authority after removing annotations"
)

# TNRS only accepts families ending in -aceae, so other families are omitted
# from the generated name string.
cat("Families not ending in -aceae (will be excluded from name string):\n")
spp_new[
    !is.na(family) & !grepl("aceae$", family, ignore.case = TRUE),
    .(family)
]

# The table, separated: reference names (Flora de Panamá) by part
spp_new[, .(mnemonic, name_type, tnrs_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "species"]
spp_new[, .(mnemonic, name_type, tnrs_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "morphospecies"]
spp_new[, .(mnemonic, name_type, tnrs_check, family, genus, species, rank, variety, authority, scientific_name)][name_type == "infraspecific"]


# --- 2b. Build the TNRS name string -----------------------------------------
# build_tnrs_name(): query string for TNRS, "<Family> <Genus> <species>
# [<rank> <infraspecific name>] [<authority>]", vectorised over its arguments.
# Empty or NA parts are left out; the family is only included when it ends in
# -aceae; an infraspecific name without a rank gets "subsp.". Returns a
# character vector.

build_tnrs_name <- function(
  family,
  genus,
  species,
  infra_rank = NULL, # e.g. "subsp.", "var.", "f."
  infra_name = NULL,
  authority = NULL
) {
    name <- genus
    name <- ifelse(
        !is.na(species) & species != "",
        paste(name, species),
        name
    )
    # Infraspecific rank + epithet
    has_infra <- !is.na(infra_name) & infra_name != ""
    name <- ifelse(
        has_infra,
        paste(
            name,
            ifelse(is.na(infra_rank) | infra_rank == "", "subsp.", infra_rank),
            infra_name
        ),
        name
    )
    name <- ifelse(
        !is.na(authority) & authority != "",
        paste(name, authority),
        name
    )
    # Prepend family if valid
    name <- ifelse(
        !is.na(family) &
            family != "" &
            grepl("aceae$", family, ignore.case = TRUE),
        paste(family, name),
        name
    )
    return(name)
}

# Row ID: the key that joins the TNRS results back to the species table
spp_new[, ID := .I]

spp_new[, name_string := build_tnrs_name(
    family = family,
    genus = genus,
    species = species,
    infra_rank = rank,
    infra_name = variety,
    authority = authority
)]

# Spot-check: confirm name strings look correct before sending to TNRS
spp_new[, .(ID, mnemonic, name_type, name_string)]
spp_new[name_type == "infraspecific", .(ID, mnemonic, name_string)]

# --- 2c. Names sent to TNRS --------------------------------------------------
# Only the names marked in 2a.1 (tnrs_check) are sent. The morphospecies are
# not published names; they are kept as they are (see 3i).
spp_new_to_check <- spp_new[tnrs_check == TRUE]
spp_new_not_sent <- spp_new[tnrs_check == FALSE]

cat(nrow(spp_new_not_sent), "rows are morphospecies and are not sent to TNRS\n")
spp_new_not_sent[, .(ID, mnemonic, family, genus, species)]
cat(nrow(spp_new_to_check), "rows will be sent to TNRS\n")

# TNRS merges rows that send the same string, so each string must be unique
check(
    !anyDuplicated(spp_new_to_check$name_string),
    "every name string sent to TNRS is unique"
)

# Prepare the two-column input expected by TNRS(): ID + name_string
tnrs_input <- spp_new_to_check[, .(ID, name_string)]

# Helpers to compare the TNRS output with the species table.
# TNRS returns "" for empty fields; turn them into NA so they compare like the
# columns of the species table.
blank_to_na <- function(x) {
    if (is.logical(x)) x <- as.character(x)
    if (!is.character(x)) {
        return(x)
    }
    fifelse(!is.na(x) & trimws(x) == "", NA_character_, x)
}
# TRUE where two values differ (NA vs NA is equal, NA vs a value differs)
differs <- function(x, y) {
    (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}

# --- 2d. TNRS parse mode: check that each name is read as intended ----------
# Parse mode only splits each string into its parts (no matching). The parts
# must be the ones separated in 2a; a difference means the string would be
# misread in resolve mode. Run this before resolve mode.

# Boyle, B. L., Matasci, N., Mozzherin, D., Rees, T., Barbosa, G. C., Kumar
# Sajja, R., & Enquist, B. J. (2021). Taxonomic Name Resolution Service, version
# 5.1. In Botanical Information and Ecology Network. https://tnrs.biendata.org/
parsed <- TNRS(
    taxonomic_names = tnrs_input,
    mode = "parse"
)
parsed <- as.data.table(parsed)[, lapply(.SD, blank_to_na)]
parsed[, ID := as.integer(ID)]

check(
    setequal(parsed$ID, tnrs_input$ID),
    "TNRS parsed every name that was sent"
)

parse_review <- merge(
    spp_new_to_check[, .(ID, mnemonic, genus, species, rank, variety, authority)],
    parsed[, .(ID, Genus, Specific_epithet, Infraspecific_rank, Infraspecific_epithet, Author, Unmatched_terms)],
    by = "ID"
)

# Genus, species, rank or variety read differently by the parser (should be empty).
# The TNRS service now and then returns, for one ID, the parse of another name
# (seen 2026-10-13: guazul came back as "Annona acuminata Saff." with its own
# name in Unmatched_terms). That is a fault of the service, not of the name:
# run the parse call again until this table is empty.
parse_review[
    differs(genus, Genus) | differs(species, Specific_epithet) |
        differs(rank, Infraspecific_rank) | differs(variety, Infraspecific_epithet)
]
check(
    parse_review[, !any(
        differs(genus, Genus) | differs(species, Specific_epithet) |
            differs(rank, Infraspecific_rank) | differs(variety, Infraspecific_epithet)
    )],
    "TNRS parsed every name into the parts of section 2a (if not, run the parse call again)"
)

# Authorities read differently. The parser drops accents and commas (e.g.
# "Müll. Arg." is read as "Mull. Arg."), which is harmless; review the rest,
# and any text the parser could not place (Unmatched_terms).
parse_review[
    differs(authority, Author) | !is.na(Unmatched_terms),
    .(mnemonic, authority, Author, Unmatched_terms)
]

# --- 2e. TNRS resolve mode: match each name in the backbone ------------------
# Only run after confirming parse output looks correct.
# sources: World Flora Online (WFO); add "wcvp" (World Checklist of Vascular
# Plants, Kew) to compare with a second backbone.
# matches = "best": return only the single best match per name.
results <- TNRS(
    taxonomic_names = tnrs_input,
    sources = c("wfo"),
    # sources = c("wfo", "wcvp"),
    matches = "best"
)
results_dt <- as.data.table(results)

# Keep the raw answer: TNRS results change when the backbone is updated, and
# the review below can be redone from this file without querying TNRS again.
checks_folder <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "CHECKS")
tnrs_results_file <- file.path(checks_folder, "species_tnrs_results_raw.rds")
saveRDS(results_dt, tnrs_results_file)
# results_dt <- readRDS(tnrs_results_file)

# --- 2f. Check the IDs returned ----------------------------------------------
# TNRS merges rows that submit the same string into one result with a
# comma-separated ID (e.g. "12,13"). The strings are unique (2c), so there
# should be none.
results_dt <- results_dt[, lapply(.SD, blank_to_na)]
results_dt[grepl(",", ID)] # merged IDs (should be empty)
check(!any(grepl(",", results_dt$ID)), "no TNRS result merges two names")

results_dt[, ID := as.integer(ID)]
check(
    setequal(results_dt$ID, tnrs_input$ID) && !anyDuplicated(results_dt$ID),
    "one TNRS result for every name sent"
)
setorder(results_dt, ID)

# Each result must belong to the name sent with its ID: the service can return,
# for one ID, the answer of another name (see 2d). Rows where the genus or
# species TNRS received is not the one sent (should be empty):
id_mixups <- merge(
    spp_new_to_check[, .(ID, mnemonic, genus, species, name_string)],
    results_dt[, .(ID, Genus_submitted, Specific_epithet_submitted, Name_submitted)],
    by = "ID"
)[differs(genus, Genus_submitted) | differs(species, Specific_epithet_submitted)]
id_mixups
check(
    nrow(id_mixups) == 0,
    "every TNRS result belongs to the name sent with its ID (if not, run the resolve call again)"
)

# --- 2g. Compare each TNRS match with the Flora de Panamá name ---------------
# One row per name sent: the name parts of the species table next to the parts
# TNRS matched, and four comparisons:
#   name_check        "same"; "spelling differs" (TNRS matched the name with
#                     another spelling of the genus, species, rank or
#                     variety); or "not matched at its rank" (no match, or a
#                     match only at a higher rank, e.g. the genus)
#   author_check      "same"; "spacing only" (same text apart from spaces,
#                     e.g. "C. DC." and "C.DC."); "different"; or
#                     "not compared" (TNRS returned no author)
#   family_check      "same" or "differs" (family of the matched name in the
#                     backbone)
#   Taxonomic_status  status of the matched name in the backbone (Accepted,
#                     Synonym, ...)
# The author compared is the author of the matched name (Author_matched), not
# that of the accepted name: for a synonym the accepted name has another author.

# Rank of each name in the words TNRS uses in Name_matched_rank
rank_words <- c("var." = "variety", "subsp." = "subspecies", "f." = "forma")

tnrs_review <- merge(
    spp_new_to_check[, .(ID, mnemonic, name_type, family, genus, species, rank, variety, authority, scientific_name)],
    results_dt[, .(
        ID, Name_matched, Name_matched_rank, Overall_score,
        Genus_matched, Genus_score,
        Specific_epithet_matched, Specific_epithet_score,
        Infraspecific_rank_matched = Infraspecific_rank,
        Infraspecific_epithet_matched, Infraspecific_epithet_score,
        Author_matched, Name_matched_accepted_family,
        Taxonomic_status, Accepted_name, Accepted_name_author, Accepted_family,
        Unmatched_terms, Name_matched_url
    )],
    by = "ID",
    all.x = TRUE
)

tnrs_review[, submitted_rank := fifelse(is.na(rank), "species", unname(rank_words[rank]))]

tnrs_review[, name_check := fcase(
    is.na(Name_matched) | differs(Name_matched_rank, submitted_rank),
    "not matched at its rank",
    differs(genus, Genus_matched) | differs(species, Specific_epithet_matched) |
        (name_type == "infraspecific" &
            (differs(rank, Infraspecific_rank_matched) | differs(variety, Infraspecific_epithet_matched))),
    "spelling differs",
    default = "same"
)]

tnrs_review[, author_check := fcase(
    is.na(Author_matched), "not compared",
    authority == Author_matched, "same",
    gsub(" ", "", authority) == gsub(" ", "", Author_matched), "spacing only",
    default = "different"
)]

tnrs_review[, family_check := fifelse(differs(family, Name_matched_accepted_family), "differs", "same")]

check(
    !anyNA(tnrs_review$name_check) && !anyNA(tnrs_review$author_check),
    "every name sent has a name and an author comparison"
)

# =============================================================================
# 3. REVIEW THE TNRS SUGGESTIONS (STEP BY STEP)
# =============================================================================
# The Flora de Panamá names have priority, and nothing below changes them
# automatically. Each step prints what TNRS suggests, and a suggestion is
# accepted only by copying its line (printed by print_accept_lines()) into the
# place indicated and running it.
# TNRS is used mainly to correct spelling: the spelling of the name parts (3c)
# and the format of the authority (3d). Differences in content (3b, 3e, 3f)
# are checked by hand. Family and synonymy in the backbone are shown for
# information only (3g, 3h) and are not applied.

# Copy of the names before the review, to record what is changed (3j)
spp_panama <- copy(spp_new)

# print_accept_lines(): print, for each suggestion, the line that accepts it,
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

# --- 3a. Overview ------------------------------------------------------------
tnrs_review[, .N, by = name_check]
tnrs_review[, .N, by = author_check]
tnrs_review[, .N, by = family_check]
tnrs_review[, .N, by = Taxonomic_status]

# --- 3b. Names not matched at their rank -------------------------------------
# TNRS found no name, or only a higher rank (e.g. the genus). Check these names
# in the Flora de Panamá list and the literature; keep them unless you find an
# error, and in that case correct the column by hand, e.g.
# spp_new[mnemonic == "xxxxxx", species := "..."]
tnrs_review[
    name_check == "not matched at its rank",
    .(mnemonic, scientific_name, authority, Name_matched, Name_matched_rank, Overall_score, Unmatched_terms)
]

# Corrections (write the lines here):


# --- 3c. Spelling of genus, species, rank and variety ------------------------
# TNRS matched the name with another spelling. Review each suggestion (score:
# 1 = identical; check the name in Tropicos or IPNI if unsure) and accept only
# true spelling mistakes.
spelling_suggestions <- rbindlist(list(
    tnrs_review[
        name_check == "spelling differs" & differs(genus, Genus_matched),
        .(mnemonic, column = "genus", current = genus, suggested = Genus_matched, score = Genus_score)
    ],
    tnrs_review[
        name_check == "spelling differs" & differs(species, Specific_epithet_matched),
        .(mnemonic, column = "species", current = species, suggested = Specific_epithet_matched, score = Specific_epithet_score)
    ],
    tnrs_review[
        name_check == "spelling differs" & name_type == "infraspecific" & differs(rank, Infraspecific_rank_matched),
        .(mnemonic, column = "rank", current = rank, suggested = Infraspecific_rank_matched, score = NA_real_)
    ],
    tnrs_review[
        name_check == "spelling differs" & name_type == "infraspecific" & differs(variety, Infraspecific_epithet_matched),
        .(mnemonic, column = "variety", current = variety, suggested = Infraspecific_epithet_matched, score = Infraspecific_epithet_score)
    ]
))
spelling_suggestions

# Lines that accept each suggestion: copy the ones you accept below and run them
print_accept_lines(spelling_suggestions)

# Accepted spelling corrections (paste the lines here):


# --- 3d. Authority: spacing only ---------------------------------------------
# Same authority written with other spaces: the backbone writes initials
# without spaces (e.g. "C.DC." for "C. DC."). Accept the ones you want written
# as in the backbone.
author_spacing <- tnrs_review[
    author_check == "spacing only",
    .(mnemonic, scientific_name, column = "authority", current = authority, suggested = Author_matched)
]
author_spacing

# Lines that accept each suggestion: copy the ones you accept below and run them
print_accept_lines(author_spacing)

# Accepted authority formats (paste the lines here):

# --- 3e. Authority: different text -------------------------------------------
# The authority differs in content. The backbone often gives the author of the
# name without the author of the basionym in parentheses, or another variant,
# so keep the Flora de Panamá authority unless it has an error. The comparison
# also helps to spot spelling mistakes in the authority: check these in IPNI or
# Tropicos and correct them by hand.
# `backbone` is the authority as WFO stores it, not a recommendation: WFO
# sometimes keeps the place of publication in it (e.g. "Exell in Pulle", where
# "in Pulle" is the book, not an author), so most of these are not corrections.
author_different <- tnrs_review[
    author_check == "different",
    .(mnemonic, scientific_name, column = "authority", current = authority, backbone = Author_matched, Taxonomic_status)
]
author_different

# Lines that copy the backbone authority as it is (only where the backbone is right)
print_accept_lines(author_different[, .(mnemonic, column, suggested = backbone)])

# Accepted or hand-corrected authorities (paste or write the lines here):


# --- 3f. Authority not compared ----------------------------------------------
# TNRS returned no author for the matched name, usually because it could not
# read the authority (see Unmatched_terms). Check these authorities by eye and
# correct them by hand, e.g.
# spp_new[mnemonic == "xxxxxx", authority := "..."]
tnrs_review[
    author_check == "not compared",
    .(mnemonic, scientific_name, authority, Unmatched_terms)
]

# Corrections (write the lines here):


# --- 3g. Family (information only) -------------------------------------------
# The family of the matched name in the backbone differs from the Flora de
# Panamá family (e.g. Cordiaceae, which WFO includes in Boraginaceae). The
# Flora de Panamá family is kept.
tnrs_review[
    family_check == "differs",
    .(mnemonic, scientific_name, family, Name_matched_accepted_family)
]

# --- 3h. Names that are synonyms in the backbone (information only) ----------
# Flora de Panamá has priority, so these names are kept. The accepted name of
# the backbone is shown for information (e.g. WFO treats the two varieties of
# Swartzia simplex as synonyms of the species).
tnrs_review[
    is.na(Taxonomic_status) | Taxonomic_status != "Accepted",
    .(mnemonic, scientific_name, Taxonomic_status, Accepted_name, Accepted_name_author, Accepted_family)
]

# --- 3i. Morphospecies (not sent to TNRS) ------------------------------------
# Kept as they are; their missing fields are completed in section 4.
spp_new_not_sent[, .(mnemonic, family, genus, species, notes)]

# --- 3j. Changes made to the Flora de Panamá names ---------------------------
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

# Save the review and the changes made
fwrite(tnrs_review, file.path(checks_folder, "species_tnrs_review.csv"))
fwrite(accepted_changes, file.path(checks_folder, "species_tnrs_accepted_changes.csv"))

# =============================================================================
# 4. COMPLETE AND CHECK THE REMAINING FIELDS
# =============================================================================

# --- 4a. Order ---------------------------------------------------------------
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

# --- 4b. Growth form (R. Pérez & S. Aguilar) ---------------------------------
# 2_prepare_viewfulltable.R uses this column (Lifeform_RPerez_SAguilar).
# Rows without one:
spp_new[
    is.na(f_de_vida_r_perez_s_aguilar),
    .(mnemonic, genus, species, f_de_vida_r_foster, pan_forma_de_vida)
]

# apeihy, nects1 and nects3 are trees (as in the previous version of this script)
spp_new[mnemonic == "apeihy", f_de_vida_r_perez_s_aguilar := "árbol"]
spp_new[mnemonic == "nects1", f_de_vida_r_perez_s_aguilar := "árbol"]
spp_new[mnemonic == "nects3", f_de_vida_r_perez_s_aguilar := "árbol"]

# ficuau (Ficus aurea) has no growth form in this column: the list of R. Pérez
# and S. Aguilar gave it to ficutr, now replaced ("árbol grande mayor a 20m"),
# and the Flora de Panamá list says "árbol". Decide and fill it, e.g.
# spp_new[mnemonic == "ficuau", f_de_vida_r_perez_s_aguilar := "..."]

spp_new[is.na(f_de_vida_r_perez_s_aguilar), .(mnemonic, genus, species)]

# --- 4c. Synonyms ------------------------------------------------------------
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
# 5. ASSEMBLE THE SPECIES TABLE AND EXPORT
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
    pan_forma_de_vida = "Lifeform_FloraPanama", # growth form, Flora de Panamá list
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

# count nrows per Lifeform_RPerez_SAguilar
full_out[, .N, by = Lifeform_RPerez_SAguilar]

# --- Export ------------------------------------------------------------------
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
