# =============================================================================
# 0_prepare_species_tables.R
#
# Purpose: Prepare the BCI 50-ha plot species table: the species list (with
#          growth forms) checked against the TNRS name service, completed from
#          the ViewTaxonomy table, with notes on known data issues.
#
# Run from the project root (or from BCI_stem_reconstruction/). The TNRS calls
# need an internet connection.
#
# Inputs (BCI_stem_reconstruction/DATA/RAW/):
#   sp_tables/Lista_bci_mnemonics_formadevida.xlsx   species list ("new" taxonomy;
#                                                    the file must be at this path)
#   ViewFiles_bci_allcensuses/ViewTaxonomy_bci.csv    database taxonomy ("old")
#   ViewFiles_bci_allcensuses/ViewFullTable_bci.csv   census records (mnemonics in use)
# Outputs (BCI_stem_reconstruction/DATA/SPP_TABLE/, which must exist):
#   bci_spptable.txt, bci_spptable.csv (both tab-separated) and
#   bci_spptable.RData (object bci.spptable), one row per mnemonic with
#   Mnemonic, Order, Family, Genus, SpeciesName, InfraspecificRank,
#   InfraspecificEpithet, Authority, Synonyms, Lifeform_RFoster,
#   Lifeform_RPerez_SAguilar, CommonName, Herbarium, Notes.
# Packages: data.table, TNRS, stringr, readxl, inspectdf.
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
# 2_prepare_viewfulltable.R; the ficutr row is removed from this table in Section 4.

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
# Standardise the taxon names and build the strings that will be sent to TNRS.
spp_new <- copy(viewt_tax1_updated_flor_pan[mnemonic != "uniden"])

# Genus must be capitalised for TNRS.
spp_new[, genus := str_to_sentence(genus)]

# Species epithets should be lower-case.
spp_new[, species := str_to_lower(species)]

spp_new[, rank := NA_character_]
spp_new[, variety := NA_character_]


spp_new[codigo %in% inc[1], `:=`(
    rank = subp_in_old[subspmnemonic %in% inc[1]]$rank,
    variety = subp_in_old[subspmnemonic %in% inc[1]]$subspecies
)]

spp_new[codigo %in% inc[2], `:=`(
    rank = subp_in_old[subspmnemonic %in% inc[2]]$rank,
    variety = subp_in_old[subspmnemonic %in% inc[2]]$subspecies
)]

spp_new[codigo %in% inc[1]]
spp_new[codigo %in% inc[2]]

# Strip non-authority annotations (sensu, auct., nom. dub., ined.) from the
# authority field; they confuse the TNRS parser
spp_new[, autoridad := stringr::str_replace_all(
    autoridad,
    regex("sensu.*|auct\\.|nom\\. dub\\.|ined\\.", ignore_case = TRUE),
    ""
)]
spp_new[, autoridad := str_trim(autoridad)]
spp_new[autoridad == "", autoridad := NA_character_]

# TNRS only accepts families ending in -aceae, so other families are omitted
# from the generated name string.
cat("Families not ending in -aceae (will be excluded from name string):\n")
spp_new[
    !is.na(family) & !grepl("aceae$", family, ignore.case = TRUE),
    .(family)
]


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

spp_new[, name_string := build_tnrs_name(
    family = familia,
    genus = genero,
    species = especie,
    infra_rank = rank,
    infra_name = variety,
    authority = autoridad
)]

# Assign a stable row ID used as the join key to merge TNRS results back
spp_new[, ID := .I]

# Spot-check: confirm name strings look correct before sending to TNRS
spp_new[, .(ID, orden, familia, genero, especie, autoridad, name_string)]

# --- 2c. Split by genus availability ----------------------------------------
# Rows without a genus cannot be sent to TNRS; set aside for manual review.

spp_new_no_genero <- spp_new[is.na(genero) | genero == ""]
spp_new_to_check <- spp_new[!is.na(genero) & genero != ""]

cat(nrow(spp_new_no_genero), "rows have no genus and will be skipped\n")
spp_new[is.na(genero)]
# apeihy, nects1 and nects3 have no genus in the new list (Apeiba "hybrida" and two Nectandra morphospecies); their taxonomy is filled from the old table in section 3
cat(nrow(spp_new_to_check), "rows will be sent to TNRS\n")

# Prepare the two-column input expected by TNRS(): ID + name_string
tnrs_input <- spp_new_to_check[, .(ID, name_string)]

# --- 2d. TNRS parse mode: verify parsing before resolving -------------------
# TNRS parse mode checks the name string syntax and component extraction.
# Run this before resolve mode to catch any badly formed names early.

# Boyle, B. L., Matasci, N., Mozzherin, D., Rees, T., Barbosa, G. C., Kumar
# Sajja, R., & Enquist, B. J. (2021). Taxonomic Name Resolution Service, version
# 5.1. In Botanical Information and Ecology Network. https://tnrs.biendata.org/
parsed <- TNRS(
    taxonomic_names = tnrs_input,
    mode = "parse"
)
parsed <- as.data.table(parsed)
parsed

# Display parsed components for review
chk <- parsed[, .(
    ID, Name_submitted, Family, Genus,
    Specific_epithet, Infraspecific_rank,
    Infraspecific_epithet, Author
)]

# Any IDs that were dropped by the parser (should be empty)
chk_missing <- setdiff(min(spp_new_to_check$ID):max(spp_new_to_check$ID), chk$ID)
spp_new_to_check[ID %in% chk_missing] # inspect dropped rows if any

# --- 2e. TNRS resolve mode: retrieve accepted names -------------------------
# Only run after confirming parse output looks correct.
# sources: WFO (World Flora Online) + WCVP (Kew Plants of the World Online)
# matches = "best": return only the single best match per name

results <- TNRS(
    taxonomic_names = tnrs_input,
    # World Flora Online - WFO
    # World Checklist of Vascular Plants - WCVP
    sources         = c("wfo", "wcvp"),
    matches         = "best"
)
results_dt <- as.data.table(results)

# --- 2f. Check for TNRS merged-ID results ------------------------------------
# TNRS can merge two rows that submitted identical name strings into a single
# result with a comma-separated ID (e.g. "1559,1558"). The two Swartzia simplex
# varieties are told apart by the variety added to their name string in
# section 2. The lines below list any merged ID and compare the returned IDs
# with those of the species list; nothing is repaired here (a merged ID would
# become NA in the numeric conversion further down).

# Merged-ID rows, if any
results_dt[grepl(",", ID)]$ID

# IDs returned versus IDs of the species list
setdiff(results_dt$ID, 1:nrow(spp_new)) # IDs in results not in spp_new
setdiff(1:nrow(spp_new), results_dt$ID) # IDs in spp_new not in results

## apeihy, nects1, and nects3 are not in results because they had no genus and were not sent to TNRS
spp_new[ID %in% setdiff(1:nrow(spp_new), results_dt$ID)] # inspect any missing IDs in results

# Convert ID to numeric for joining and sort for readability
results_dt[, ID := as.numeric(as.character(ID))]
setorder(results_dt, ID)

# --- 2g. Diagnose problems ---------------------------------------------------
# Assign a human-readable problem label to each row using TNRS score columns.
# Rules are evaluated in order; fcase() returns the first matching condition.
# Labels are written in Spanish to match the final output format.

results_dt[, problem := fcase(
    # ── MORPHOSPECIES ─────────────────────────────────────────────────────
    grepl("\\bsp\\.?\\s*(\\d+|nov\\.?)?$",
        Name_submitted,
        ignore.case = TRUE
    ),
    "Morfoespecie — solo género, sin resolución a especie posible",
    # ── PERFECT MATCHES ───────────────────────────────────────────────────
    Overall_score == 1 &
        Taxonomic_status == "Accepted" &
        !is.na(Accepted_name),
    "OK",
    Overall_score == 1 & is.na(Accepted_name),
    "Sin resolver en backbone — verificar literatura",
    # ── AUTHORITY FORMAT ONLY ─────────────────────────────────────────────
    Genus_score == 1 &
        Specific_epithet_score == 1 &
        Family_score == 1 &
        Taxonomic_status == "Accepted" &
        !is.na(Accepted_name),
    "OK — solo formato de autoridad",
    # ── FAMILY NOT RECOGNISED ─────────────────────────────────────────────
    Warnings == 4 &
        Genus_score == 1 &
        Specific_epithet_score == 1,
    "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    Warnings == 4 &
        (Genus_score < 1 | Specific_epithet_score < 1),
    "Familia no está en backbone + problema en nombre — verificar manualmente",
    # ── SYNONYMS ──────────────────────────────────────────────────────────
    Taxonomic_status == "Synonym",
    "Sinónimo — reemplazar con nombre_aceptado",
    # ── INFRASPECIFIC ERRORS ──────────────────────────────────────────────
    Genus_score == 1 &
        Specific_epithet_score == 1 &
        Infraspecific_epithet_score < 1,
    "Epíteto infraespecífico mal escrito",
    # ── COMPONENT SPELLING ERRORS ─────────────────────────────────────────
    Genus_score < 1 &
        Specific_epithet_score == 1,
    "Género mal escrito",
    Genus_score == 1 &
        Specific_epithet_score < 1,
    "Epíteto específico mal escrito",
    Genus_score < 1 &
        Specific_epithet_score < 1,
    "Género y epíteto mal escritos",
    Family_score < 1 &
        Genus_score == 1 &
        Specific_epithet_score == 1,
    "Solo discrepancia en familia",
    # ── LOW CONFIDENCE ────────────────────────────────────────────────────
    Overall_score < 0.5,
    "Coincidencia pobre — verificar manualmente",
    is.na(Name_matched),
    "Sin coincidencia encontrada",
    default = "Revisar — caso no clasificado"
)]

# Frequency table of problem categories — useful triage overview
results_dt[, .N, by = problem][order(-N)]

# --- 2h. Select and retain TNRS output columns ------------------------------
# Keep only the columns used in the curation below.

# Inspect any Swartzia varieties that may need special handling.
results_dt[Genus_submitted == "Swartzia"]

# Inspect rows that did not match cleanly.
results_dt[Unmatched_terms != ""]

cols_tnrs_clean <- c(
    # Join key
    "ID",
    "Accepted_name", # full currently-accepted name
    "Accepted_species", # binomial only, no infraspecific rank
    "Accepted_name_author", # standardised authority
    "Accepted_family", # correct family per backbone
    "Infraspecific_rank", # infraspecific rank assigned by TNRS (e.g. subsp., var., f.)
    "Infraspecific_epithet_matched", # infraspecific epithet matched by TNRS (if any)
    # Custom diagnosis
    "problem"
)

results_dt <- results_dt[, ..cols_tnrs_clean]

# --- 2i. Merge TNRS results back to new taxonomy -----------------------------

# Remove the helper column created for the TNRS submission
spp_new[, name_string := NULL]

# Full join on ID (all = TRUE): every row of spp_new is kept, with the TNRS
# columns where a result exists
spp_new <- merge(
    spp_new,
    results_dt,
    by  = "ID",
    all = TRUE
)

# =============================================================================
# 3. APPLY TNRS CORRECTIONS
# =============================================================================
# Use the TNRS diagnosis labels to apply safe fixes.
# - Not sent to TNRS (no genus): accepted here, completed from the old taxonomy below.
# - Reclassified family: update only the family.
# - Morphospecies: accept genus-only names.
# - Authority formatting: standardise authority text.
# After these steps, clean up temporary TNRS helper columns and finalise the table.

# Add a stable row index for reference during curation
spp_new[, I := .I]

# Create a working copy of the diagnosis column; `solution` will be updated
# to "OK" as each issue is resolved, leaving unresolved cases visible.
spp_new[, solution := problem]

# Frequency table before corrections — baseline overview
message("Problem distribution before corrections:")
print(table(spp_new$solution, useNA = "ifany"))

# --- Check 1: Rows that were not sent to TNRS --------------------------------
# The rows without a genus in the new list (apeihy, nects1, nects3) were not
# sent to TNRS and therefore have NA in `problem`. Their family, genus and
# species are filled from the old taxonomy further below; mark them as OK.
spp_new[is.na(solution)]
spp_new[is.na(solution), solution := "OK"]
print(table(spp_new$solution, useNA = "ifany"))

# --- Check 2: Family reclassified by backbone --------------------------------
# Warnings = 4 in TNRS indicates the submitted family is not recognised by
# WFO/WCVP (e.g. Cordiaceae → Boraginaceae, Meliaceae reclassifications).
# The species name is correct; only the family needs updating.
message("Families to update (reclassified in backbone):")
print(
    spp_new[
        solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
        .(familia, Accepted_family)
    ]
)

#       familia familia_aceptada
# 1: Cordiaceae     Boraginaceae
# 2: Cordiaceae     Boraginaceae
# 3: Cordiaceae     Boraginaceae
# 4:  Meliaceae        Meliaceae

# Previous family is Boraginaceae
# We'll keep it like that
sp_bci_raw_input[Family == "Cordiaceae"]
sp_bci_raw_input[Family == "Boraginaceae"]

spp_new[
    solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    familia := Accepted_family # overwrite with backbone-accepted family
]
spp_new[
    solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    solution := "OK"
]
message("After Check 2 (family reclassified):")
print(table(spp_new$solution, useNA = "ifany"))

# --- Check 3: Morphospecies (sp.) -------------------------------------------
# Entries with "sp." cannot be resolved beyond genus level; this is expected
# behaviour, not an error. No name fields need updating.
spp_new[solution == "Morfoespecie — solo género, sin resolución a especie posible"]
spp_new[
    solution == "Morfoespecie — solo género, sin resolución a especie posible",
    solution := "OK"
]
message("After Check 3 (morphospecies):")
print(table(spp_new$solution, useNA = "ifany"))

# --- Check 4: Authority format difference only --------------------------------
# The name and family are correct; only the punctuation/spacing of the
# authority string differs from the backbone standard. Overwrite with the
# standardised backbone authority.
message("Authority differences to standardise:")
print(
    spp_new[
        solution == "OK — solo formato de autoridad",
        .(genero, especie, Accepted_name, autoridad, Accepted_name_author)
    ]
)

spp_new[
    solution == "OK — solo formato de autoridad",
    `:=`(
        autoridad = Accepted_name_author,
        solution  = "OK"
    )
]
message("After Check 4 (authority format):")
print(table(spp_new$solution, useNA = "ifany"))

# --- Final authority sweep ---------------------------------------------------
# Catch any remaining rows where the authority still differs from the backbone
# value (e.g. minor formatting differences not covered by Check 4).
message("Remaining authority mismatches after all checks:")
print(spp_new[autoridad != Accepted_name_author, .(codigo, autoridad, Accepted_name_author)])
spp_new[
    !is.na(Accepted_name_author) & autoridad != Accepted_name_author,
    autoridad := Accepted_name_author
]

spp_new[, texto := NULL] # column of the species list that is not exported
spp_new[, fotos := NULL] # column of the species list that is not exported

spp_new[autoridad != Accepted_name_author]
spp_new[, Accepted_name_author := NULL]
spp_new[, Accepted_name := NULL]
spp_new[, Accepted_species := NULL]

spp_new[familia != Accepted_family]
spp_new[, Accepted_family := NULL]

spp_new[, Infraspecific_rank := ifelse(Infraspecific_rank == "", NA_character_, Infraspecific_rank)]
spp_new[, Infraspecific_epithet_matched := ifelse(Infraspecific_epithet_matched == "", NA_character_, Infraspecific_epithet_matched)]

# Swartzia simplex carries an intraspecific rank (not a true synonym); clear sinonimos for those rows.
spp_new[!is.na(Infraspecific_rank) | !is.na(Infraspecific_epithet_matched)]$codigo
spp_new[codigo %in% spp_new[!is.na(Infraspecific_rank) | !is.na(Infraspecific_epithet_matched)]$codigo, sinonimos := NA_character_]

spp_new[, rank := NULL]
spp_new[, variety := NULL]

table(spp_new$solution)

spp_new[solution == "Sinónimo — reemplazar con nombre_aceptado"]
# Swartzia simplex is flagged as a synonym because of the presence of an intraspecific rank; no fix required.

cols_to_keep <- c(
    "codigo", "orden", "familia", "genero",
    "especie", "Infraspecific_rank", "Infraspecific_epithet_matched",
    "autoridad", "sinonimos", "f_de_vida_r_foster", "f_de_vida_r_perez_s_aguilar",
    "nombre_comum", "herbario", "notes"
)

spp_new <- spp_new[, ..cols_to_keep]

inc <- unique(spp_new[is.na(orden) | is.na(familia) | is.na(genero) | is.na(especie)]$codigo)

to_replace <- spp_old[mnemonic %in% inc]

# Fill missing taxonomy fields from the old ViewTaxonomy (the codes without a genus in the new list)
spp_new[codigo %in% inc[1]]
spp_new[
    codigo %in% inc[1],
    `:=`(
        familia = to_replace[mnemonic %in% inc[1]]$family,
        genero = to_replace[mnemonic %in% inc[1]]$genus,
        especie = to_replace[mnemonic %in% inc[1]]$speciesname
    )
]

spp_new[codigo %in% inc[2]]
spp_new[
    codigo %in% inc[2],
    `:=`(
        familia = to_replace[mnemonic %in% inc[2]]$family,
        genero = to_replace[mnemonic %in% inc[2]]$genus,
        especie = to_replace[mnemonic %in% inc[2]]$speciesname
    )
]

spp_new[codigo %in% inc[3]]
spp_new[
    codigo %in% inc[3],
    `:=`(
        familia = to_replace[mnemonic %in% inc[3]]$family,
        genero = to_replace[mnemonic %in% inc[3]]$genus,
        especie = to_replace[mnemonic %in% inc[3]]$speciesname
    )
]

spp_new[is.na(orden) | is.na(familia) | is.na(genero) | is.na(especie)]

spp_new[is.na(orden) & familia == "Malvaceae", orden := unique(spp_new[!is.na(orden) & familia == "Malvaceae"]$orden)]
spp_new[is.na(orden) & familia == "Lauraceae", orden := unique(spp_new[!is.na(orden) & familia == "Lauraceae"]$orden)]

spp_new[is.na(orden) | is.na(familia) | is.na(genero) | is.na(especie)]

spp_new[codigo == "apeihy"]
spp_new[codigo == "apeihy", f_de_vida_r_perez_s_aguilar := "árbol"]

spp_new[codigo == "nects1"]
spp_new[codigo == "nects1", f_de_vida_r_perez_s_aguilar := "árbol"]

spp_new[codigo == "nects3"]
spp_new[codigo == "nects3", f_de_vida_r_perez_s_aguilar := "árbol"]

spp_new[codigo %in% c(inc, "swars1", "swars2")]

# =============================================================================
# 4. ASSEMBLE FINAL bci.spptable AND EXPORT
# =============================================================================

# Retain only the columns needed for the species table; drop all comparison
# and curation helper columns.
cols_clean <- c(
    "codigo", # BCI mnemonic code (primary key)
    "orden", # accepted order
    "familia", # accepted family
    "genero", # accepted genus
    "especie", # accepted species epithet
    "Infraspecific_rank", # infraspecific rank assigned by TNRS (e.g. subsp., var., f.)
    "Infraspecific_epithet_matched", # infraspecific epithet matched by TNRS (if any)
    "autoridad", # accepted authority (standardised by TNRS)
    "sinonimos", # previous synonymised name where applicable
    "f_de_vida_r_foster", # growth form (from Foster's list)
    "f_de_vida_r_perez_s_aguilar", # growth form (from Pérez & Aguilar's list)
    "nombre_comum", # common name
    "herbario", # herbarium voucher reference
    "notes"
)

full_out <- spp_new[, ..cols_clean]

names(full_out)

# --- Rename columns to English ------------------------------------------

setnames(
    full_out,
    old = c(
        "codigo", "orden", "familia", "genero", "especie", "Infraspecific_rank", "Infraspecific_epithet_matched",
        "autoridad", "sinonimos", "f_de_vida_r_foster", "f_de_vida_r_perez_s_aguilar", "nombre_comum", "herbario", "notes"
    ),
    new = c(
        "Mnemonic", "Order", "Family", "Genus", "SpeciesName", "InfraspecificRank", "InfraspecificEpithet",
        "Authority", "Synonyms", "Lifeform_RFoster", "Lifeform_RPerez_SAguilar", "CommonName", "Herbarium", "Notes"
    )
)

# count nrows per Lifeform_RPerez_SAguilar
full_out[, .N, by = Lifeform_RPerez_SAguilar]

# --- Export ------------------------------------------------------------------
bci.spptable <- full_out

output_path <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "SPP_TABLE")

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
