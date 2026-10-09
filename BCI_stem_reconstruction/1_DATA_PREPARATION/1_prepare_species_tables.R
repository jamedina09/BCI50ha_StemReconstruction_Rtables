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

# Number of threads data.table will use
data.table::getDTthreads()

# =============================================================================
# DEFINE INPUT FILE PATHS
# =============================================================================
workspace_root <- getwd()
if (basename(workspace_root) == "BCI_stem_reconstruction") {
    workspace_root <- dirname(workspace_root)
}

TAXONOMY_NEW <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "sp_tables", "Lista_bci_mnemonics_formadevida.xlsx")
TAXONOMY_OLD <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses", "ViewTaxonomy_bci.csv")
INPUT_FILE <- file.path(workspace_root, "BCI_stem_reconstruction", "DATA", "RAW", "ViewFiles_bci_allcensuses", "ViewFullTable_bci.csv")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

# --- 1a. New taxonomy --------------------------------------------------------
spp_new <- as.data.table(readxl::read_xlsx(TAXONOMY_NEW))

# Remove any leading/trailing whitespace introduced by Excel
colnames(spp_new) <- stringr::str_trim(colnames(spp_new))

# Standardise column names: lower-case, spaces → underscores,
# strip accents (UTF-8 → ASCII transliteration), remove remaining
# non-alphanumeric characters so names are safe for data.table operations.
colnames(spp_new) <- tolower(colnames(spp_new))
colnames(spp_new) <- gsub(" ", "_", colnames(spp_new))
colnames(spp_new) <- iconv(colnames(spp_new), from = "UTF-8", to = "ASCII//TRANSLIT")
colnames(spp_new) <- gsub("[^[:alnum:]_]", "", colnames(spp_new))

# Inspect missing values to understand which columns are incomplete
inspectdf::inspect_na(spp_new)


# --- 1b. Old taxonomy table (ViewTaxonomy): load and clean column names ------
spp_old <- as.data.table(fread(TAXONOMY_OLD))
colnames(spp_old) <- tolower(colnames(spp_old))
colnames(spp_old) <- stringr::str_trim(colnames(spp_old))
# Convert literal "NULL" strings to NA for proper handling
spp_old[spp_old == "NULL"] <- NA

# --- 1c. BCI inventory (ViewFullTable) ---------------------------------------
sp_bci_raw_input <- as.data.table(fread(INPUT_FILE))
# Unique Tag–Mnemonic–CensusID combinations (one row per tree per census)
sp_bci_raw <- unique(sp_bci_raw_input[, .(Tag, Mnemonic, CensusID)])

# =============================================================================
# 1.1. CHECK RAW TAXONOMY DATA
# =============================================================================
# Known data problems in the new taxonomy table: the lines below print the
# rows concerned and record a note (column `notes`) where one applies.
spp_new[, notes := NA_character_]

# Example: Appunia siebertiii is correct; Morinda siebertii is a synonym.
spp_new[especie == "siebertii"]

# Apeiba "hybrida" is no longer present in the plot.
spp_new[codigo %in% "apeihy"]
# present in census 1:3
sp_bci_raw[Mnemonic == "apeihy"]
sp_bci_raw_input[Mnemonic == "apeihy", .(CensusID, ExactDate)]

spp_new[codigo %in% "apeihy", notes := "No longer present; last seen 1990"]
spp_new[codigo %in% "apeihy"]

# Nectandra s1 and Nectandra s3 are Lauraceae morphospecies; one of them died.
# Only one remains alive, and it should be monitored for fertility-based identification.
spp_new[codigo %in% "nects1"]
spp_new[codigo %in% "nects3"]

# present in census 1:9
sp_bci_raw[Mnemonic %in% "nects1"]

# present in census 1:3
sp_bci_raw[Mnemonic %in% "nects3"]
sp_bci_raw_input[Mnemonic == "nects3", .(CensusID, ExactDate)]
spp_new[codigo %in% "nects3", notes := "No longer present; last seen 1990"]

# Pterocarpus officinalis no existe en la parcela, yo personalmente revisé todos los Pterocarpus y
# corresponden a P. rohrii.
spp_new[especie %in% "rohrii"]
spp_new[codigo %in% "pterro"]
spp_new[codigo %in% "pterro", notes := "Previously identified as 'pterof' (Pterocarpus officinalis) are incorrect; those are 'pterro' (Pterocarpus rohrii)"]

# this needs to be replaced with the correct code
unique(sp_bci_raw_input[SpeciesName == "officinalis", .(Mnemonic, Family, Genus, SpeciesName)])

# Correct:
sp_bci_raw_input[Mnemonic == "pterro"]
# Replace with pterro:
sp_bci_raw_input[Mnemonic == "pterof"]

# FIXME: pterof should be replaced by pterro in the BCI inventory (raw data).
# Below, the replacement is made in the copy held in memory, which this script
# only uses to compare mnemonics; 1_prepare_viewfulltable.R.R makes the same
# replacement in the table it prepares.

# Beilschmiedia and Quararibea: compare the codes of the species list with
# those used in the raw inventory.
spp_new[especie %in% "pendula"]
spp_new[especie %in% "asterolepis"]
spp_new[especie %in% "tovarensis"]
spp_new[especie %in% "stenophylla"]

unique(sp_bci_raw_input[Genus == "Beilschmiedia", .(Mnemonic, Family, Genus, SpeciesName)])
unique(sp_bci_raw_input[Genus == "Quararibea", .(Mnemonic, Family, Genus, SpeciesName)])

spp_new[codigo %in% c("beilpe", "quaras")]

# The only inventory fix applied here is the known mnemonic correction
# `pterof -> pterro`.
sp_bci_raw_input[, Mnemonic := ifelse(Mnemonic == "pterof", "pterro", Mnemonic)]
bci_data_mnemonic <- sort(unique(sp_bci_raw_input[, Mnemonic]))

# Find mnemonics present in the BCI inventory that are missing from the new list.
inc <- setdiff(bci_data_mnemonic, spp_new$codigo)
unique(sp_bci_raw_input[Mnemonic %in% inc, .(Mnemonic, Family, Genus, SpeciesName)])

# =============================================================================
# 2. PREPARE NEW TAXONOMY FOR TNRS VALIDATION
# =============================================================================
# Standardise the taxon names and build the strings that will be sent to TNRS.

# Genus must be capitalised for TNRS.
spp_new[, genero := str_to_sentence(genero)]

# Species epithets should be lower-case.
spp_new[, especie := str_to_lower(especie)]

# Preserve variety/subspecies data when it appears in the old taxonomy.
subp_in_old <- spp_old[subspmnemonic %in% spp_new$codigo & !is.na(subspecies)]
subp_in_old <- subp_in_old[rank == "var."]

spp_new[, rank := NA_character_]
spp_new[, variety := NA_character_]
inc <- sort(subp_in_old$subspmnemonic)

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
    !is.na(familia) & !grepl("aceae$", familia, ignore.case = TRUE),
    .(familia)
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
