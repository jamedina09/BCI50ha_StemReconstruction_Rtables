# =============================================================================
# 0_prepare_species_tables.R
#
# Purpose: Validate and unify the new Panamá flora taxonomy (Rolando Pérez,
#          Nov 2025) against the old ForestGEO/CTFS ViewTaxonomy using the TNRS
#          backbone, then reconcile both sources to produce a final BCI 50-ha
#          species table (bci.spptable) with corrected names, families, and
#          authorities.
#
# Workflow overview:
#    1. Load source tables (new taxonomy, old taxonomy, BCI inventory).
#    2. Validate the new taxonomy through TNRS
#       (clean → build name string → parse → resolve → diagnose → export).
#    3. Summarise the BCI inventory by mnemonic (n individuals per species).
#    4. Subset both taxonomy sources to BCI mnemonics only.
#    5. Diagnose coverage discrepancies between sources.
#    6. Prepare old taxonomy columns for merging.
#    7. Merge old and new taxonomy (full outer join).
#    8. Initial export of the merged lookup table.
#    9. Apply TNRS corrections (family reclassifications, synonyms, authority
#       standardisation, morphospecies handling).
#   10. Build a comparison table of old vs. updated name fields.
#   11. Assemble the final bci.spptable and export.
#
# Inputs:
#   DATA/RAW/TAXONOMY_ROLANDO_P/
#       Lista de especies de la Flora de Panamá_nov25 original.xlsx  <- new taxonomy
#   DATA/RAW/ViewFiles_bci_allcensuses/ViewTaxonomy_bci.csv          <- old taxonomy
#   DATA/PROCESSED/1_ViewFullTable_no_missing_tags_census.rds        <- BCI inventory
#
# Outputs:
#   DATA/RAW/TAXONOMY_ROLANDO_P/
#       Lista_de_especies_de_la_Flora_de_Panamá_nov25_original_con_puntajes.xlsx
#       check_bci_mnemonics.xlsx
#       check_bci_mnemonics_to_check.txt
#       check_bci_mnemonics_Rolando_comparison.txt
#   DATA/SPP_TABLE/
#       bci_spptable.txt
#       bci_spptable.RData
# =============================================================================
rm(list = ls())

options(
    datatable.print.class = FALSE,
    datatable.print.keys = TRUE,
    datatable.verbose = FALSE
)

# =============================================================================
# PACKAGES
# =============================================================================
library(data.table) # fast data manipulation
library(here) # project-relative paths
library(TNRS) # Taxonomic Name Resolution Service client
library(stringr) # string helpers (str_to_sentence, str_trim, etc.)

data.table::setDTthreads(1) # use all available CPU threads for data.table
data.table::getDTthreads()

# =============================================================================
# FILE PATHS
# =============================================================================
TAXONOMY_NEW <- here(
    "DATA", "RAW", "TAXONOMY_ROLANDO_P",
    "Lista de especies de la Flora de Panamá_nov25 original.xlsx"
)
TAXONOMY_OLD <- here(
    "DATA", "RAW", "ViewFiles_bci_allcensuses", "ViewTaxonomy_bci.csv"
)
INPUT_FILE <- here(
    "DATA", "PROCESSED", "1_ViewFullTable_no_missing_tags_census.rds"
)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

# --- 1a. New taxonomy (Rolando Pérez, Nov 2025) ------------------------------
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

# --- 1b. Old taxonomy (ForestGEO/CTFS ViewTaxonomy) --------------------------
spp_old <- as.data.table(fread(TAXONOMY_OLD))
colnames(spp_old) <- tolower(colnames(spp_old))
colnames(spp_old) <- stringr::str_trim(colnames(spp_old))
# Replace literal "NULL" strings (CSV artefact from the database export)
# with proper NA values so they behave correctly in logical tests
spp_old[spp_old == "NULL"] <- NA

# --- 1c. BCI inventory (pre-processed; no missing tags/censuses) -------------
sp_bci_raw_input <- as.data.table(readRDS(INPUT_FILE))
# Keep only unique Tag–Mnemonic pairs so each physical tree is counted once;
# census-level rows are not needed at this stage.
sp_bci_raw <- unique(sp_bci_raw_input[, .(Tag, Mnemonic)])

# Quick spot-check for species known to have data-quality issues
sp_bci_raw_input[Mnemonic %in% "apeihy"]
# unique(sp_bci_raw_input[Mnemonic %in% "nects1", .(Tag, CensusID, ExactDate)])
# unique(sp_bci_raw_input[Mnemonic %in% "nects3", .(Tag, CensusID, ExactDate)])

# =============================================================================
# 2. PREPARE NEW TAXONOMY FOR TNRS VALIDATION
# =============================================================================
# TNRS (Taxonomic Name Resolution Service) resolves submitted plant names
# against curated backbones (WFO, WCVP) and returns accepted names, synonyms,
# spelling corrections, and match scores.  The workflow is:
#   Step 2a: Clean and standardise name components in spp_new.
#   Step 2b: Build the TNRS name string from those components.
#   Step 2c: Parse mode  — verify TNRS parses each component correctly.
#   Step 2d: Resolve mode — retrieve accepted names and match scores.
#   Step 2e: Fix any TNRS output artefacts (e.g. merged IDs).
#   Step 2f: Diagnose each row with a human-readable problem label.
#   Step 2g: Subset, rename, and merge results back to spp_new.

# --- 2a. Standardise name components ----------------------------------------

# TNRS requires genus with first letter capitalised
spp_new[, genero := str_to_sentence(genero)]

# Species epithets must be lower-case for correct parsing
spp_new[, especie := str_to_lower(especie)]

# No subspecies column in the source data — create it as NA so the
# build_tnrs_name() helper below has a consistent signature
spp_new[, subespecie := NA_character_]

# Strip non-authority annotations (sensu, auct., nom. dub., ined.) from the
# authority field; they confuse the TNRS parser
spp_new[, autoridad := stringr::str_replace_all(
    autoridad,
    regex("sensu.*|auct\\.|nom\\. dub\\.|ined\\.", ignore_case = TRUE),
    ""
)]
spp_new[, autoridad := str_trim(autoridad)]
spp_new[autoridad == "", autoridad := NA_character_]

# Families that do NOT end in -aceae will be excluded from the name string
# because TNRS only accepts families with the standard -aceae suffix
cat("Families not ending in -aceae (will be excluded from name string):\n")
spp_new[
    !is.na(familia) & !grepl("aceae$", familia, ignore.case = TRUE),
    .(familia)
]

# --- 2b. Build the TNRS name string -----------------------------------------
# TNRS expects a single string per taxon that may include:
#   [Family] Genus [species [subsp. subspecies]] [authority]
# Family is prepended only when it ends in -aceae (TNRS parser requirement).

build_tnrs_name <- function(family, genus, species, subspecies, authority) {
    # Start with genus (minimum required by TNRS)
    name <- genus
    # Append species epithet when available
    name <- ifelse(!is.na(species) & species != "",
        paste(name, species),
        name
    )
    # Append infraspecific rank + epithet when available
    name <- ifelse(!is.na(subspecies) & subspecies != "",
        paste(name, "subsp.", subspecies),
        name
    )
    # Append authority when available
    name <- ifelse(!is.na(authority) & authority != "",
        paste(name, authority),
        name
    )
    # Prepend family only if it ends in -aceae
    name <- ifelse(
        !is.na(family) & family != "" & grepl("aceae$", family, ignore.case = TRUE),
        paste(family, name),
        name
    )
    return(name)
}

spp_new[, name_string := build_tnrs_name(
    family     = familia,
    genus      = genero,
    species    = especie,
    subspecies = subespecie,
    authority  = autoridad
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
cat(nrow(spp_new_to_check), "rows will be sent to TNRS\n")

# Prepare the two-column input expected by TNRS(): ID + name_string
tnrs_input <- spp_new_to_check[, .(ID, name_string)]

# --- 2c. TNRS parse mode: verify parsing before resolving -------------------
# Parse mode does NOT query the backbone; it only shows how TNRS splits the
# name string into components (Family, Genus, epithet, Author).
# Inspect the output to confirm no components are mis-assigned before the
# more expensive resolve call.

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

# --- 2d. TNRS resolve mode: retrieve accepted names -------------------------
# Only run after confirming parse output looks correct.
# sources: WFO (World Flora Online) + WCVP (Kew Plants of the World Online)
# matches = "best": return only the single best match per name

results <- TNRS(
    taxonomic_names = tnrs_input,
    sources         = c("wfo", "wcvp"),
    matches         = "best"
)
results_dt <- as.data.table(results)

# --- 2e. Fix TNRS merged-ID artefact ----------------------------------------
# TNRS occasionally merges two rows that submitted identical name strings into
# a single result with a comma-separated ID (e.g. "1559,1558").
# This happens here for the two Swartzia simplex subspecies.
# Fix: duplicate that result row, assign each original ID separately.

# Identify any merged-ID rows
results_dt[grepl(",", ID)]
# Expected output:
#           ID    Accepted_name     Accepted_species  Accepted_name_author Taxonomic_status Overall_score
# 1: 1559,1558 Swartzia simplex  Swartzia simplex    (Sw.) Spreng.         Accepted         1

double_row <- copy(results_dt[ID == "1559,1558"])

# Create one row per original ID
row_1558 <- copy(double_row)
row_1558[, ID := "1558"]

row_1559 <- copy(double_row)
row_1559[, ID := "1559"]

# Remove the merged row and insert the two corrected rows
results_dt <- results_dt[ID != "1559,1558"]
results_dt <- rbind(results_dt, row_1558, row_1559)

results_dt[ID %in% c("1558", "1559")] # confirm both rows are present with correct IDs

# Verify that the corrected results cover all expected IDs
setdiff(results_dt$ID, 1:nrow(spp_new)) # IDs in results not in spp_new
setdiff(1:nrow(spp_new), results_dt$ID) # IDs in spp_new not in results

# Convert ID to numeric for joining and sort for readability
results_dt[, ID := as.numeric(as.character(ID))]
setorder(results_dt, ID)

# Confirm all consecutive IDs differ by 1 (no gaps or duplicates)
unique(diff(results_dt$ID))

# --- 2f. Diagnose problems ---------------------------------------------------
# Assign a human-readable problem label to each row using TNRS score columns.
# Rules are evaluated in order; fcase() returns the first matching condition.
# Labels are written in Spanish to match the final output format.

results_dt[, problem := fcase(
    # ── PERFECT MATCHES ───────────────────────────────────────────────────────
    # Score = 1, status Accepted, accepted name present → no action needed
    Overall_score == 1 & Taxonomic_status == "Accepted" & !is.na(Accepted_name),
    "OK",
    # Score = 1 but no accepted name → taxon unresolved in the backbone
    Overall_score == 1 & is.na(Accepted_name),
    "Sin resolver en backbone — verificar literatura",
    # ── AUTHORITY FORMAT ONLY ─────────────────────────────────────────────────
    # All name components match (genus, epithet, family) but authority differs
    # in punctuation or spacing → correct name, no taxonomic action needed
    Genus_score == 1 & Specific_epithet_score == 1 &
        Family_score == 1 & Taxonomic_status == "Accepted" & !is.na(Accepted_name),
    "OK — solo formato de autoridad",
    # ── SYNONYMS ──────────────────────────────────────────────────────────────
    # Matched name exists but is outdated → replace with Accepted_name
    Taxonomic_status == "Synonym",
    "Sinónimo — reemplazar con nombre_aceptado",
    # ── FAMILY NOT RECOGNISED BY BACKBONE ────────────────────────────────────
    # Warnings = 4: submitted family was rejected (reclassified / synonymised).
    # The species name itself may still be correct.
    Warnings == 4 & Genus_score == 1 & Specific_epithet_score == 1,
    "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    Warnings == 4 & (Genus_score < 1 | Specific_epithet_score < 1),
    "Familia no está en backbone + problema en nombre — verificar manualmente",
    # ── MORPHOSPECIES (sp.) ───────────────────────────────────────────────────
    # Entry ends in sp. / sp.N / sp.nov. → resolved to genus only; expected
    grepl("\\bsp\\.?\\s*(\\d+|nov\\.?)?$", Name_submitted, ignore.case = TRUE),
    "Morfoespecie — solo género, sin resolución a especie posible",
    # ── COMPONENT SPELLING ERRORS ─────────────────────────────────────────────
    Genus_score < 1 & Specific_epithet_score == 1,
    "Género mal escrito",
    Genus_score == 1 & Specific_epithet_score < 1,
    "Epíteto específico mal escrito",
    Genus_score < 1 & Specific_epithet_score < 1,
    "Género y epíteto mal escritos",
    Infraspecific_epithet_score < 1,
    "Epíteto infraespecífico mal escrito",
    Family_score < 1 & Genus_score == 1 & Specific_epithet_score == 1,
    "Solo discrepancia en familia",
    # ── LOW CONFIDENCE ────────────────────────────────────────────────────────
    Overall_score < 0.5,
    "Coincidencia pobre — verificar manualmente",
    is.na(Name_matched),
    "Sin coincidencia encontrada",
    default = "Revisar — caso no clasificado"
)]

# Frequency table of problem categories — useful triage overview
results_dt[, .N, by = problem][order(-N)]

# --- 2g. Select and retain TNRS output columns ------------------------------
# Keep only the columns that are actionable for downstream curation.
# See the reference block below for a full description of each column.

cols_tnrs_clean <- c(
    # Join key
    "ID",
    # Accepted name components (what to correct TO)
    "Accepted_name", # full currently-accepted name
    "Accepted_species", # binomial only, no infraspecific rank
    "Accepted_name_author", # standardised authority
    # "Accepted_name_rank",   # rank (species / subspecies / variety)
    "Accepted_family", # correct family per backbone
    # Name status and match quality
    "Taxonomic_status", # Accepted / Synonym / No opinion
    "Overall_score", # main quality filter (0–1)
    "Genus_score", # pinpoints genus misspelling
    "Specific_epithet_score", # pinpoints epithet misspelling
    "Family_score", # pinpoints family mismatch (Warnings = 4 cases)
    # What TNRS matched
    "Name_matched", # closest match in backbone (may be a synonym)
    # "Name_matched_rank",    # rank of the matched name
    "Unmatched_terms", # parts TNRS could not parse — non-empty = problem
    # Custom diagnosis
    "problem"
)

results_dt <- results_dt[, ..cols_tnrs_clean]

# ================================================================================
# REFERENCE: TNRS KEPT COLUMNS — DESCRIPTION AND PURPOSE
# ================================================================================
# ── IDENTITY ─────────────────────────────────────────────────────────────────────
# ID
#     Row number linking the TNRS results back to your original data (spp_new).
#     Used as the join key when merging corrections back into your table.
# ── WHAT YOU NEED TO CORRECT TO ──────────────────────────────────────────────────
# Accepted_name
#     The currently accepted full name you should be using in your dataset.
#     This is your primary correction output. If your name was a synonym or
#     misspelled, this column gives you the right name to replace it with.
#     Will be NA if the taxon has unresolved status in the backbone.
# Accepted_species
#     The binomial only (Genus + species epithet) of the accepted name,
#     without infraspecific rank, authority, or any other annotation.
#     Useful when you want a clean species-level name regardless of whether
#     the matched taxon was a subspecies or variety.
# Accepted_name_author
#     The standardized authority of the accepted name as recorded in the
#     backbone (WFO / WCVP). Use this to overwrite your authority column
#     with the correct, standardized format.
# Accepted_name_rank
#     The taxonomic rank of the accepted name: species, subspecies, variety,
#     forma, etc. Important to check when you submitted a species name but
#     the backbone accepted it only at genus level, or vice versa.
# Accepted_family
#     The accepted family for this taxon according to the backbone.
#     Critical for the Warnings = 4 cases, where the family you prepended
#     to the name string is not recognized (e.g. Cordiaceae, Metteniusaceae)
#     because it has been reclassified or synonymized. This column tells you
#     what family the backbone currently places the species in.
# ── WAS YOUR NAME RIGHT OR OUTDATED? ─────────────────────────────────────────────
# Taxonomic_status
#     The status of the name TNRS matched, not of your submitted name.
#     Possible values:
#         Accepted     The matched name is the current accepted name.
#                      If Overall_score = 1, your name is correct.
#         Synonym      The matched name exists in the backbone but is outdated.
#                      You must use Accepted_name instead.
#         No opinion   The backbone has no consensus on this name.
#                      Treat with caution and verify manually.
# ── HOW CONFIDENT IS THE MATCH? ──────────────────────────────────────────────────
# All scores run from 0 to 1. A score of 1 means a perfect match for that
# component. Lower values indicate how far your submitted name deviates
# from the backbone name.
# Overall_score
#     The combined match score across all name components. Your main filter
#     for deciding which names need attention.
#     Interpretation guide:
#         1.00          Perfect match (all components matched exactly)
#         0.90 - 0.99   Minor issue, usually authority formatting only
#         0.80 - 0.89   Moderate issue, inspect the component scores
#         < 0.80        Serious problem, manual check required
# Genus_score
#     Match score for the genus component alone.
#     If < 1 while Specific_epithet_score = 1, the genus is misspelled
#     or does not exist in the backbone.
# Specific_epithet_score
#     Match score for the species epithet component alone.
#     If < 1 while Genus_score = 1, the epithet is misspelled or wrong.
#     If both Genus_score and Specific_epithet_score are < 1, both
#     components have problems.
# Family_score
#     Match score for the family component alone.
#     A value < 1 typically means the family you prepended to the name
#     string is not recognized by the backbone (Warnings = 4 cases).
#     In these cases, check Accepted_family for the correct family.
# ── WHAT DID TNRS ACTUALLY MATCH? ────────────────────────────────────────────────
# Name_matched
#     The closest name found in the backbone for your submitted string.
#     Note: this may still be a synonym. Always use Accepted_name for
#     corrections, not Name_matched.
#     If NA, TNRS found no plausible match at all.
# Name_matched_rank
#     The taxonomic rank at which TNRS was able to match your name:
#     species, subspecies, variety, genus, etc.
#     If your submitted name was a full binomial but this returns "genus",
#     TNRS could only resolve to genus level — the species epithet was
#     not recognized. This is expected for morphospecies (sp.) entries.
# Unmatched_terms
#     Any part of your submitted name string that TNRS could not parse
#     or assign to a known name component.
#     If not empty, something in your name string confused the parser —
#     for example, collection numbers, annotations like "cf." or "aff.",
#     or morphospecies numbers (e.g. "sp. 7").
#     Always inspect rows where this is not empty.
# ── YOUR DIAGNOSTIC COLUMN ───────────────────────────────────────────────────────
# problem
#     Your custom classification of what kind of issue each row has.
#     Use this as your triage column to prioritize corrections.
#     Possible values and recommended actions:
#     OK
#         Name is accepted, score is perfect. No action needed.
#     OK — authority format only
#         Name and family are correct. Score is below 1 only because
#         your authority string differs slightly in punctuation or spacing
#         from the backbone version. No taxonomic action needed. Optionally
#         overwrite your authority with Accepted_name_author.
#     Synonym — replace with Accepted_name
#         Your name exists in the backbone but is outdated. Replace your
#         name with the value in Accepted_name.
#     Unresolved in backbone — check literature
#         TNRS matched the name perfectly (score = 1) but the backbone
#         has no accepted name recorded for it. The taxon may be recently
#         described, have uncertain placement, or not yet processed by
#         WFO/WCVP. Requires manual literature check.
#     Family not in backbone — reclassified (check Accepted_family)
#         The family you prepended is not recognized by WFO/WCVP because
#         it has been synonymized or reclassified (e.g. Cordiaceae is now
#         part of Boraginaceae). The species name itself is likely correct.
#         Check Accepted_family for the current family placement.
#     Family not in backbone + name issue — manual check
#         Same as above but the genus or epithet also has a match problem.
#         Both the family and the name components need to be verified.
#     Morphospecies — genus only, no species resolution possible
#         Your entry contains "sp." indicating an unidentified specimen.
#         TNRS resolved to genus level only. This is expected behaviour,
#         not an error. No correction possible until the specimen is
#         identified to species.
#     Genus misspelled
#         Genus_score < 1 but Specific_epithet_score = 1. The genus name
#         has a spelling error. Check Name_matched and Accepted_name for
#         the correct spelling.
#     Species epithet misspelled
#         Specific_epithet_score < 1 but Genus_score = 1. The species
#         epithet has a spelling error. Check Accepted_name for the
#         correct spelling.
#     Both genus and epithet misspelled
#         Both Genus_score and Specific_epithet_score are < 1. Both
#         components have problems. Use Accepted_name as reference and
#         verify manually.
#     Infraspecific epithet misspelled
#         The subspecies or variety epithet has a spelling error.
#         Check Accepted_name for the correct infraspecific epithet.
#     Family mismatch only
#         Family_score < 1 but genus and epithet are fine. The family
#         you provided does not match what the backbone expects for this
#         taxon. Check Accepted_family.
#     Poor match — check manually
#         Overall_score < 0.5. TNRS could not find a confident match.
#         The name may be too different from any backbone entry, contain
#         unusual annotations, or be genuinely absent from the backbone.
#         Requires manual verification.
#     No match found
#         Name_matched is NA. TNRS found nothing in the backbone
#         resembling your submitted name. Check for major misspellings,
#         invalid names, or names not yet in WFO/WCVP.
#     Review — unclassified
#         Does not fit any of the above rules. Inspect the score columns
#         and Warnings individually to determine the issue.
#     Recommended triage order:
#         1. Synonym — replace with Accepted_name         (easy fix)
#         2. OK — authority format only                   (optional fix)
#         3. Family not in backbone — reclassified        (update family)
#         4. Genus / epithet misspelled                   (fix spelling)
#         5. Unresolved in backbone — check literature    (manual)
#         6. Poor match / No match found                  (manual)
#         7. Morphospecies                                (no fix possible)
# ================================================================================

# --- 2h. Merge TNRS results back to new taxonomy and export -----------------

# Confirm no unexpected values remain in the temporary subspecies column
unique(spp_new$subespecie)

# Remove helper columns created for TNRS submission (no longer needed)
spp_new[, subespecie := NULL]
spp_new[, name_string := NULL]

# Left join: keep all rows of spp_new; attach TNRS scores where available
spp_new <- merge(
    spp_new,
    results_dt,
    by  = "ID",
    all = TRUE
)

# Rename TNRS output columns to Spanish to match the rest of the table
col_translation <- c(
    "Accepted_name"          = "nombre_aceptado",
    "Accepted_species"       = "especie_aceptada",
    "Accepted_name_author"   = "autoridad_aceptada",
    "Accepted_name_rank"     = "rango_aceptado",
    "Accepted_family"        = "familia_aceptada",
    "Taxonomic_status"       = "estado_taxonomico",
    "Overall_score"          = "puntaje_general",
    "Genus_score"            = "puntaje_genero",
    "Specific_epithet_score" = "puntaje_especie",
    "Family_score"           = "puntaje_familia",
    "Name_matched"           = "nombre_coincidente",
    "Name_matched_rank"      = "rango_coincidente",
    "Unmatched_terms"        = "terminos_no_reconocidos",
    "problem"                = "problema"
)

# Only rename TNRS columns that are actually present (safe if some were excluded)
tnrs_cols_present <- intersect(names(col_translation), names(spp_new))
setnames(spp_new, old = tnrs_cols_present, new = col_translation[tnrs_cols_present])

# Export annotated new taxonomy (drop the row-index ID before writing)
writexl::write_xlsx(
    spp_new[, ID := NULL],
    here(
        "DATA", "RAW", "TAXONOMY_ROLANDO_P",
        "Lista_de_especies_de_la_Flora_de_Panamá_nov25_original_con_puntajes.xlsx"
    )
)
message("Exported: TNRS-annotated new taxonomy")

# =============================================================================
# 3. SUMMARISE MNEMONICS IN THE BCI INVENTORY
# =============================================================================
# Count the number of unique tag records per mnemonic (= n individuals observed
# across all censuses for that species).
sp_bci <- sp_bci_raw[, .(ntags_in_bci = .N), by = Mnemonic]
setorder(sp_bci, Mnemonic)

message(sprintf("Distinct mnemonics in BCI inventory: %d", length(unique(sp_bci$Mnemonic))))

# =============================================================================
# 4. SUBSET TAXONOMY FILES TO BCI MNEMONICS
# =============================================================================

# --- 4a. Old taxonomy: subset to BCI mnemonics ------------------------------
# Keep only rows whose mnemonic code appears in the BCI inventory
old_sp <- spp_old[mnemonic %in% sp_bci$Mnemonic]

# Detect mnemonics that map to more than one row in the old taxonomy
# (can happen when the same mnemonic is used for species and subspecies)
old_sp[, duplicated_mnemonic := duplicated(mnemonic) | duplicated(mnemonic, fromLast = TRUE)]
old_sp[duplicated_mnemonic == TRUE] # inspect duplicated rows

# For duplicated mnemonics, replace the mnemonic with the subspecies mnemonic
# so each row gets a unique identifier before re-filtering
old_sp[duplicated_mnemonic == TRUE, mnemonic := subspmnemonic]

# Drop the helper flag column
old_sp[, duplicated_mnemonic := NULL]

# Re-filter after disambiguation to ensure we only retain BCI mnemonics
old_sp <- old_sp[mnemonic %in% sp_bci$Mnemonic]

# --- 4b. New taxonomy: flag BCI species --------------------------------------
# Add a flag column to spp_new indicating whether each species code (codigo)
# is present in the BCI 50-ha plot inventory
spp_new[, BCI50ha := ifelse(codigo %in% sp_bci$Mnemonic, "si", "no")]

# Quick coverage check (commented out; uncomment to re-run)
# setdiff(sp_bci$Mnemonic, old_sp$mnemonic)  # BCI codes absent from old taxonomy
# setdiff(old_sp$mnemonic, sp_bci$Mnemonic)  # old taxonomy codes absent from BCI

# =============================================================================
# 5. DIAGNOSTIC: COVERAGE DISCREPANCIES BETWEEN SOURCES
# =============================================================================
# Compare which mnemonics are present in the BCI inventory, the old taxonomy,
# and the new taxonomy to detect gaps before merging.

# --- 5a. BCI mnemonics missing from the old taxonomy -------------------------
missing_from_old <- setdiff(sp_bci$Mnemonic, old_sp$mnemonic)
message("BCI mnemonics absent from old taxonomy 'mnemonic' column:")
print(missing_from_old)

# Check whether the missing codes appear as subspecies mnemonics instead
# (e.g. swars2 is stored under subspmnemonic, not mnemonic, in the old table)
still_missing <- setdiff(missing_from_old, old_sp$subspmnemonic)
message("Still missing after checking 'subspmnemonic' column:")
print(still_missing)

# Spot-check known subspecies/mnemonic mapping cases
message("Subspecies entry for swars2 in old taxonomy:")
print(old_sp[subspmnemonic == "swars2", .(mnemonic, subspmnemonic)])
# Expected: mnemonic = swars2, subspmnemonic = swars2

message("Subspecies entry for protte in old taxonomy:")
print(old_sp[mnemonic == "protte", .(mnemonic, subspmnemonic)])
# Expected: mnemonic = protte, subspmnemonic = protte

message("Subspecies entry for quaras in old taxonomy:")
print(old_sp[mnemonic == "quaras", .(mnemonic, subspmnemonic)])
# Expected: mnemonic = quaras, subspmnemonic = quaras

# Final symmetric checks after disambiguation
setdiff(sp_bci$Mnemonic, old_sp$mnemonic) # still missing from old taxonomy
setdiff(old_sp$mnemonic, sp_bci$Mnemonic) # old taxonomy codes absent from BCI

# --- 5b. BCI mnemonics missing from the new taxonomy -------------------------
missing_from_new <- setdiff(sp_bci$Mnemonic, spp_new[BCI50ha == "si", codigo])
message("BCI mnemonics absent from new taxonomy 'codigo' column:")
print(missing_from_new)
# Expected: "apeihy" "nects1" "nects3" "uniden"
# (3 species not yet in the new list + 1 unidentified-individual code)

# --- 5c. Old taxonomy codes absent from the new taxonomy --------------------
message("Old taxonomy mnemonics not found in new taxonomy:")
print(setdiff(old_sp$mnemonic, spp_new[BCI50ha == "si", codigo]))
# Expected: "apeihy" "nects1" "nects3" "uniden"

# --- 5d. New taxonomy codes absent from the old taxonomy --------------------
message("New taxonomy codes not found in old taxonomy:")
print(setdiff(spp_new[BCI50ha == "si", codigo], old_sp$mnemonic))

# Summary of what the merged table will look like:
#   apeihy, nects1, nects3, uniden → old-taxonomy columns present, new-taxonomy NAs
#   all other BCI species           → both sources populated

# =============================================================================
# 6. PREPARE OLD TAXONOMY FOR MERGING
# =============================================================================

# Select the relevant taxonomy columns and prefix all with "prev_" to
# distinguish old-taxonomy columns from new-taxonomy columns in the merged table
old_sp_out <- old_sp[, .(
    family, genus, speciesname, subspecies, authority,
    subspauthority, rank, idlevel, mnemonic, subspmnemonic, listofoldnames
)]
colnames(old_sp_out) <- paste("prev", colnames(old_sp_out), sep = "_")

# Check for any remaining duplicated mnemonics in the prepared table
# (should be empty after disambiguation in step 4)
dup_mnemonics <- old_sp_out[
    duplicated(prev_mnemonic) | duplicated(prev_mnemonic, fromLast = TRUE),
    prev_mnemonic
]
message("Duplicated prev_mnemonic in old_sp_out (should be empty):")
print(dup_mnemonics)

# Inspect all mnemonics for reference
print(old_sp_out$prev_mnemonic)

# Attach inventory tag counts so each row carries abundance information
old_sp_out <- merge(
    old_sp_out, sp_bci,
    by.x = "prev_mnemonic", by.y = "Mnemonic",
    all.x = TRUE
)

# Move ntags_in_bci to the first column for easier visual inspection
setcolorder(old_sp_out, c("ntags_in_bci", setdiff(colnames(old_sp_out), "ntags_in_bci")))

# Print the prepared table for a final visual check
print(old_sp_out)

# Sanity check: any old mnemonic that maps to >1 row in the new taxonomy
# would create spurious duplicate rows in the merge — detect and report
dup_check <- old_sp_out[
    prev_mnemonic %in% unique(spp_new[BCI50ha == "si", codigo]),
    .N,
    by = prev_mnemonic
][N > 1]

if (nrow(dup_check) > 0) {
    warning("Old mnemonics with >1 match in new taxonomy — review before merging:")
    print(dup_check)
    # Inspect the problematic rows in both tables for manual resolution
    spp_new[BCI50ha == "si" & (codigo %in% dup_check$prev_mnemonic)]
    old_sp[mnemonic %in% dup_check$prev_mnemonic]
} else {
    message("OK: no old mnemonic matches >1 row in the new taxonomy")
}

# =============================================================================
# 7. MERGE OLD AND NEW TAXONOMY
# =============================================================================
# Full outer join (all = TRUE) ensures species present in only one source are
# retained with NAs for the missing side.
#   Left  (new taxonomy, BCI rows only): join key = codigo
#   Right (old taxonomy):                join key = prev_mnemonic
# Restricting the new taxonomy to BCI50ha == "si" keeps the merged table
# focused on the species actually present in the BCI 50-ha plot.

full_out <- merge(
    spp_new[BCI50ha == "si"], old_sp_out,
    by.x = "codigo", by.y = "prev_mnemonic",
    all = TRUE
)

# Drop the helper BCI flag column (no longer needed after join)
full_out[, BCI50ha := NULL]

# Spot-check a known species to verify merge correctness
full_out[codigo == "quara1"]
old_sp_out[prev_mnemonic == "quara1"]
spp_new[codigo == "quara1"]
sp_bci[Mnemonic == "quara1"]

# =============================================================================
# 8. INITIAL EXPORT OF MERGED LOOKUP TABLE
# =============================================================================
# Exclude the "uniden" (unidentified) code: it carries no taxonomic information
# and would only add noise to the curator review.
# Note: this table was shared with Rolando Pérez and David Roubik for review.

# Remove "uniden" from the in-memory object before any further processing
full_out <- full_out[codigo != "uniden"]

# Excel workbook — primary output for manual curation (shared with curators)
writexl::write_xlsx(
    full_out,
    here("DATA", "RAW", "TAXONOMY_ROLANDO_P", "check_bci_mnemonics.xlsx")
)
message("Exported: check_bci_mnemonics.xlsx")

# Drop annotation columns that are internal to the new taxonomy source file
# and are not relevant for the cross-taxonomy comparison or spptable:
#   mnbc  — internal note field from Rolando's spreadsheet
#   texto — free-text field from Rolando's spreadsheet
#   fotos — photo reference field from Rolando's spreadsheet
full_out[, mnbc := NULL]
full_out[, texto := NULL]
full_out[, fotos := NULL]

# Tab-delimited text file — useful for version control and diff comparisons
fwrite(
    full_out,
    here("DATA", "RAW", "TAXONOMY_ROLANDO_P", "check_bci_mnemonics_to_check.txt"),
    sep = "\t"
)
message("Exported: check_bci_mnemonics_to_check.txt")

# =============================================================================
# 9. APPLY TNRS CORRECTIONS
# =============================================================================
# The `problema` column (populated in section 2f) classifies each species by
# the type of name issue detected. This section resolves each issue category
# in turn, updating the relevant name fields in `full_out` and marking the
# row as "OK" in a working `solution` column.  The checks are applied in
# order of increasing complexity:
#   Check 1 — Species absent from Rolando's list          (mark OK as-is)
#   Check 2 — Family reclassified in backbone             (update familia)
#   Check 3 — Synonym: replace with accepted name         (update all name fields)
#   Check 4 — Morphospecies (sp.)                         (mark OK, no fix needed)
#   Check 5 — Authority format difference only            (update autoridad)
# After all checks, any remaining authority mismatches are standardised.

# Add a stable row index for reference during curation
full_out[, I := .I]

# Create a working copy of the diagnosis column; `solution` will be updated
# to "OK" as each issue is resolved, leaving unresolved cases visible.
full_out[, solution := problema]

# Frequency table before corrections — baseline overview
message("Problem distribution before corrections:")
print(table(full_out$solution, useNA = "ifany"))

# --- Check 1: Species not present in Rolando's list --------------------------
# Some BCI mnemonics (e.g. apeihy, nects1, nects3) were not in the new
# taxonomy and therefore have NA in `problema`. They already carry valid
# name data from the old taxonomy (prev_* columns); mark them as OK.
full_out[is.na(solution), solution := "OK"]
message("After Check 1 (absent from Rolando list):")
print(table(full_out$solution, useNA = "ifany"))

# --- Check 2: Family reclassified by backbone --------------------------------
# Warnings = 4 in TNRS indicates the submitted family is not recognised by
# WFO/WCVP (e.g. Cordiaceae → Boraginaceae, Meliaceae reclassifications).
# The species name is correct; only the family needs updating.
message("Families to update (reclassified in backbone):")
print(
    full_out[
        solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
        .(familia, familia_aceptada)
    ]
)

full_out[
    solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    familia := familia_aceptada # overwrite with backbone-accepted family
]
full_out[
    solution == "Familia no está en backbone — reclasificada (verificar familia_aceptada)",
    solution := "OK"
]
message("After Check 2 (family reclassified):")
print(table(full_out$solution, useNA = "ifany"))

# --- Check 3: Synonym — replace with accepted name ---------------------------
# The submitted name is a synonym; update family, genus, species, and authority
# with the currently accepted values from the TNRS backbone.
# `especie_aceptada` contains the full binomial (Genus + epithet); we split it
# to fill the separate genero and especie fields.
# The original matched name is preserved in `sinonimos` for traceability.
message("Synonyms to update:")
print(full_out[solution == "Sinónimo — reemplazar con nombre_aceptado"])

full_out[
    problema == "Sinónimo — reemplazar con nombre_aceptado",
    `:=`(
        familia   = familia_aceptada,
        genero    = str_split_fixed(especie_aceptada, " ", 2)[, 1], # genus component
        especie   = str_split_fixed(especie_aceptada, " ", 2)[, 2], # epithet component
        autoridad = autoridad_aceptada,
        sinonimos = nombre_coincidente, # record the old (synonymised) name
        solution  = "OK"
    )
]
message("After Check 3 (synonyms replaced):")
print(table(full_out$solution, useNA = "ifany"))

# --- Check 4: Morphospecies (sp.) -------------------------------------------
# Entries with "sp." cannot be resolved beyond genus level; this is expected
# behaviour, not an error. No name fields need updating.
full_out[
    solution == "Morfoespecie — solo género, sin resolución a especie posible",
    solution := "OK"
]
message("After Check 4 (morphospecies):")
print(table(full_out$solution, useNA = "ifany"))

# --- Check 5: Authority format difference only --------------------------------
# The name and family are correct; only the punctuation/spacing of the
# authority string differs from the backbone standard. Overwrite with the
# standardised backbone authority.
message("Authority differences to standardise:")
print(
    full_out[
        solution == "OK — solo formato de autoridad",
        .(autoridad, autoridad_aceptada)
    ]
)

full_out[
    solution == "OK — solo formato de autoridad",
    `:=`(
        autoridad = autoridad_aceptada,
        solution  = "OK"
    )
]
message("After Check 5 (authority format):")
print(table(full_out$solution, useNA = "ifany"))

# --- Final authority sweep ---------------------------------------------------
# Catch any remaining rows where the authority still differs from the backbone
# value (e.g. minor formatting differences not covered by Check 5).
message("Remaining authority mismatches after all checks:")
print(full_out[autoridad != autoridad_aceptada, .(codigo, autoridad, autoridad_aceptada)])
full_out[
    !is.na(autoridad_aceptada) & autoridad != autoridad_aceptada,
    autoridad := autoridad_aceptada
]

# =============================================================================
# 10. BUILD COMPARISON TABLE (OLD vs. UPDATED NAME FIELDS)
# =============================================================================
# Select the columns needed for a side-by-side comparison of the previous
# (old taxonomy) and updated (TNRS-corrected new taxonomy) name fields.
# This table is shared with curators to document what changed for each species.

cols_Rolando_clean <- c(
    "codigo", # BCI mnemonic (join key)
    # "orden",         # taxonomic order — omitted; not used downstream
    "familia", # updated family (corrected by TNRS where needed)
    "genero", # updated genus
    "especie", # updated species epithet
    "autoridad", # updated authority (standardised by TNRS)
    "sinonimos", # old synonymised name (populated only for CHECK 3 rows)
    "forma_de_vida", # growth form from Rolando's list
    "nombre_comum", # common name
    "herbario", # herbarium voucher reference
    # "ntags_in_bci",  # tag count — omitted from curator output
    # Previous (old ForestGEO taxonomy) name fields — all prefixed prev_
    "prev_family", "prev_genus", "prev_speciesname", "prev_subspecies",
    "prev_authority", "prev_subspauthority", "prev_rank", "prev_idlevel",
    "prev_subspmnemonic", "prev_listofoldnames"
)

full_out <- full_out[, ..cols_Rolando_clean]

# Create flag columns to highlight where the updated value differs from the
# previous value — makes it easy to filter changed rows in Excel
full_out[
    ,
    `:=`(
        updated_family    = ifelse(familia != prev_family, "yes", "no"),
        updated_genus     = ifelse(genero != prev_genus, "yes", "no"),
        updated_species   = ifelse(especie != prev_speciesname, "yes", "no"),
        updated_authority = ifelse(autoridad != prev_authority, "yes", "no")
    )
]

# Print changed rows for each field for a quick review
message("Species with updated family:")
print(full_out[updated_family == "yes", .(codigo, prev_family, familia)])
message("Species with updated genus:")
print(full_out[updated_genus == "yes", .(codigo, prev_genus, genero)])
message("Species with updated species epithet:")
print(full_out[updated_species == "yes", .(codigo, prev_speciesname, especie)])
message("Species with updated authority:")
print(full_out[updated_authority == "yes", .(codigo, prev_authority, autoridad)])

# Export comparison table (tab-delimited for easy diff / version control)
fwrite(
    full_out,
    here("DATA", "RAW", "TAXONOMY_ROLANDO_P", "check_bci_mnemonics_Rolando_comparison.txt"),
    sep = "\t"
)
message("Exported: check_bci_mnemonics_Rolando_comparison.txt")

# =============================================================================
# 11. ASSEMBLE FINAL bci.spptable AND EXPORT
# =============================================================================

# --- 11a. Handle subspecies -------------------------------------------------
# Inspect rows that carry subspecies information from the old taxonomy
print(full_out[!is.na(prev_subspecies)])

# Replace empty strings with NA for consistency (empty string = no subspecies)
full_out[prev_subspecies == "", prev_subspecies := NA]

# Promote prev_subspecies to a standalone `subspecies` column in the output
full_out[, subspecies := prev_subspecies]

# For subspecies rows, the synonym field was incorrectly populated with the
# parent species name during Check 3. Clear it to avoid misleading output.
full_out[!is.na(prev_subspecies), sinonimos := NA_character_]

# Inspect remaining NAs across all columns
message("Missing value summary:")
print(data.table(inspectdf::inspect_na(full_out)))

# --- 11b. Fill missing names from old taxonomy -------------------------------
# Species that were absent from Rolando's list (Check 1) have NA in the
# name fields (familia, genero, especie). Populate them from the old
# ForestGEO taxonomy columns so no species row is left without a name.
full_out[
    is.na(especie),
    `:=`(
        familia = prev_family,
        genero  = prev_genus,
        especie = prev_speciesname
    )
]

# --- 11c. Select final columns ----------------------------------------------
# Retain only the columns needed for the species table; drop all comparison
# and curation helper columns.
cols_clean <- c(
    "codigo", # BCI mnemonic code (primary key)
    "familia", # accepted family
    "genero", # accepted genus
    "especie", # accepted species epithet
    "subspecies", # subspecies epithet (NA for full species)
    "autoridad", # accepted authority (standardised by TNRS)
    "sinonimos", # previous synonymised name where applicable
    "forma_de_vida", # growth form (from Rolando's list)
    "nombre_comum", # common name
    "herbario" # herbarium voucher reference
)

full_out <- full_out[, ..cols_clean]

# Final NA check before renaming
print(data.table(inspectdf::inspect_na(full_out)))

# --- 11d. Rename columns to English ------------------------------------------
# ForestGEO/CTFS conventions use English column names; rename to match.
setnames(
    full_out,
    old = c(
        "codigo", "familia", "genero", "especie", "subspecies",
        "autoridad", "sinonimos", "forma_de_vida", "nombre_comum", "herbario"
    ),
    new = c(
        "Mnemonic", "Family", "Genus", "SpeciesName", "Subspecies",
        "Authority", "Synonyms", "Lifeform", "CommonName", "Herbarium"
    )
)

table(full_out$Lifeform, useNA = "ifany") # check lifeform categories and missing values

# translate to English and standardise lifeform categories
full_out[
    Lifeform %in% c("árbol", "arbolito"),
    Lifeform := "tree"
]
full_out[
    Lifeform == "arbusto",
    Lifeform := "shrub"
]
full_out[
    Lifeform == "Helecho arbóreo",
    Lifeform := "tree_fern"
]
full_out[
    Lifeform == "Palma",
    Lifeform := "palm"
]

# if commonname has "abraza palo" and genus is Ficus, 
# change Lifeform to strangler_fig
full_out[
    grepl("abraza palo", CommonName, ignore.case = TRUE) & Genus == "Ficus",
    Lifeform := "strangler_fig"
]

# full_out[is.na(Lifeform), Lifeform := "tree"] # replace any remaining NAs with "unknown"

chk <- sp_bci_raw_input[Mnemonic %in% full_out[is.na(Lifeform)]$Mnemonic, .(Mnemonic, CensusID, DBH)]
chk[, DBH_cm := DBH / 10]

full_out[is.na(Lifeform), Lifeform := "tree"]

table(full_out$Lifeform, useNA = "ifany") # final check of lifeform categories

# --- 11e. Export -------------------------------------------------------------
bci.spptable <- full_out

# Tab-delimited plain text — portable, version-control friendly
fwrite(
    bci.spptable,
    here("DATA", "SPP_TABLE", "bci_spptable.txt"),
    sep = "\t"
)
message("Exported: DATA/SPP_TABLE/bci_spptable.txt")

# R binary format — for direct use in downstream R scripts
save(
    bci.spptable,
    file = here("DATA", "SPP_TABLE", "bci_spptable.RData")
)
message("Exported: DATA/SPP_TABLE/bci_spptable.RData")
