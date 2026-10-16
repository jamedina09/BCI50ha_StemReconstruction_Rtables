################################################################################
# SPECIES TABLE OF THE R TABLES
################################################################################
# Purpose: Build the species table of the BCI R tables (bci.spptable), one
# row per species code of the R tables, with the names of the stage-1 species
# table and notes on how codes and names changed since the last public
# release of the plot data (Condit et al. 2019); and the changes of species
# code since that release, one row per change and one row per tree.
#
# Run from the project root, after 2_create_R_tables_BCI.R (it reads the R
# tables), line by line or with source().
#
# Last public release, the reference of the notes:
#   Condit R, Pérez R, Aguilar S, Lao S, Foster R, Hubbell S (2019) Complete
#   data from the Barro Colorado 50-ha plot: 423617 trees, 35 years. Dryad,
#   doi:10.15146/5xcp-0d46 (published 2019-06-07; version 3, 2020-06-24,
#   changed only the metadata).
#
# Inputs (BCI_stem_reconstruction/DATA/):
#   RTABLES/bci.stem<n>.Rdata             the R tables (2_create_R_tables_BCI.R):
#                                         trees, tags and codes
#   SPP_TABLE/bci_spptable.RData          species table of stage 1
#                                         (1_DATA_PREPARATION/1_prepare_species_tables.R);
#                                         its columns are copied unchanged
#   SPP_TABLE/condit/bci.spptable.rdata   Condit's species table (the file of
#                                         Dryad version 3; SHA-256 checked)
#   RAW/bci_dryad_condit.rds              Condit's plot data: the 8 Dryad stem
#                                         tables (bci.stem1-8) in one table
#
# Outputs:
#   DATA/RTABLES/bci.spptable.rdata (object bci.spptable, a data.frame) and
#   bci.spptable.csv                      the species table: one row per code
#                                         of the R tables, in order of code,
#                                         uniden last
#   DATA/RTABLES/bci.species_code_changes_summary.rdata (object
#   bci.species_code_changes_summary) and .csv
#                                         one row per change of code since
#                                         Condit et al. (2019): spCondit2019,
#                                         LatinCondit2019, sp, Latin, trees
#                                         (number of trees)
#   DATA/RTABLES/bci.species_code_changes.rdata (object
#   bci.species_code_changes) and .csv    one row per tree whose code changed
#                                         since Condit et al. (2019): treeID,
#                                         tag, spCondit2019, LatinCondit2019,
#                                         sp, Latin
#   DATA/CHECKS/check_species_table_bci.txt  failed checks and notes
#
# Columns of bci.spptable (all text, NA when empty):
#   sp ... Herbarium  the columns of the stage-1 table (Mnemonic renamed sp),
#                     unchanged, plus Latin = Genus SpeciesName and, when
#                     there is one, the infraspecific name (InfraspecificRank
#                     InfraspecificEpithet), as in Condit's species table
#   NameCondit2019    name of the code in Condit et al. (2019), when its
#                     genus, epithet or infraspecific name differs from Latin
#   NameChange        kind of change from NameCondit2019 to Latin: genus;
#                     genus and epithet; genus and epithet gender; epithet;
#                     epithet gender; epithet spelling; infraspecific name;
#                     morphospecies label (name_change_kind() in
#                     species_table_functions.R). It describes the change, not
#                     which spelling is correct: the names are those of the
#                     stage-1 table
#   ReidentifiedFrom  codes that the trees now under this code had in Condit
#                     et al. (2019), with the number of trees: "tremsp (84)"
#   ReidentifiedTo    codes that the trees of this code in Condit et al.
#                     (2019) have now, with the number of trees
#   Notes             the note of the stage-1 table, then sentences built from
#                     the columns above: the former name, a code not used in
#                     Condit et al. (2019), the trees of a code no longer used
# Only the codes and the names (genus, epithet, infraspecific name) are
# compared with Condit et al. (2019); authority, family and the other columns
# are not.
#
# uniden (unidentified stems) has no row in the stage-1 table, which lists
# taxa. With INCLUDE_UNIDEN = TRUE (section 1) this script adds a row without
# names as the last row of the table, so that every code of the R tables has
# a row: a merge() of a census table with the species table (an inner join by
# default) keeps the uniden stems, and the 47 trees recorded as ingaun (Inga
# unidentified) in Condit et al. (2019) have a row for their note. Condit's
# species table has the code too ("Unidentified species").
#
# Sections:
#   1. Setup
#   2. Load: R tables, stage-1 species table, Condit's species table and plot data
#   3. Codes of the R tables
#   4. Trees whose code changed since Condit et al. (2019)
#   5. Names compared with Condit et al. (2019)
#   6. Notes
#   7. Assemble and check
#   8. Export
#
# Packages: data.table, digest.
################################################################################

# ========================================================================
# SECTION 1: SETUP
# ========================================================================

# Remove all objects from workspace to start fresh
rm(list = ls())

library(data.table)
# digest (SHA-256 of a file) is called with digest::

# Paths (run from the project root)
main_path <- getwd()
site <- "bci"
DATA_folder <- file.path(main_path, "BCI_stem_reconstruction", "DATA")
RTABLES_folder <- file.path(DATA_folder, "RTABLES")
CHECK_folder <- file.path(DATA_folder, "CHECKS")
SPP_STAGE1_FILE <- file.path(DATA_folder, "SPP_TABLE", "bci_spptable.RData")
CONDIT_SPPTABLE_FILE <- file.path(DATA_folder, "SPP_TABLE", "condit", "bci.spptable.rdata")
CONDIT_PLOT_FILE <- file.path(DATA_folder, "RAW", "bci_dryad_condit.rds")

# Condit et al. (2019): SHA-256 of bci.spptable.rdata in Dryad version 3, as
# listed by https://datadryad.org/api/v2/versions/31849/files (checked on
# 2026-10-15; 50,401 bytes), and the number of trees in the title of the dataset
CONDIT_SPPTABLE_SHA256 <- "254f42986af6e1e28306f25e39910b9920448d28982bbda7124643d998a43955"
CONDIT_N_TREES <- 423617L
CONDIT_REF <- "Condit et al. (2019)"

# Code of the unidentified stems, and whether it gets a row (see the header)
UNIDEN_CODE <- "uniden"
INCLUDE_UNIDEN <- TRUE

# Log of failed checks and notes (an existing one is replaced)
dir.create(CHECK_folder, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(CHECK_folder, paste0("check_species_table_", site, ".txt"))
if (file.exists(log_file)) file.remove(log_file)
invisible(file.create(log_file))

# ── Hard check that stays visible in interactive (line-by-line) runs ────────
# As in 2_create_R_tables_BCI.R: on failure it prints a ❌ line (count +
# examples), raises an immediate warning, writes the same text to the log
# file, and only then stops. On success it prints a ✓ line.
bio_check <- function(ok, msg, examples = NULL, n_bad = NULL) {
  if (isTRUE(all(ok))) {
    cat("✓", msg, "\n")
    return(invisible(TRUE))
  }
  n_txt <- if (!is.null(n_bad)) sprintf(" [%d case(s)]", n_bad) else ""
  ex_txt <- if (length(examples) > 0L) {
    paste0(" | examples: ", paste(head(unique(examples), 10), collapse = ", "))
  } else {
    ""
  }
  full_msg <- paste0("CHECK FAILED: ", msg, n_txt, ex_txt)
  cat("❌", full_msg, "\n")
  warning(full_msg, call. = FALSE, immediate. = TRUE)
  cat(full_msg, file = log_file, append = TRUE, sep = "\n")
  stop(full_msg, call. = FALSE)
}

# Visible note that does not stop (also written to the log file)
note <- function(msg) {
  cat("⚠️ ", msg, "\n")
  cat(paste("NOTE:", msg), file = log_file, append = TRUE, sep = "\n")
}

# Codes as one printable line
show_codes <- function(x) if (length(x)) paste(sort(x), collapse = ", ") else "(none)"

# Name, code and note rules (shared with tests/test_species_table.R)
source(file.path(main_path, "BCI_stem_reconstruction", "3_PREPARE_R_TABLES", "species_table_functions.R"))

inputs <- c(SPP_STAGE1_FILE, CONDIT_SPPTABLE_FILE, CONDIT_PLOT_FILE)
bio_check(all(file.exists(inputs)), "all input files exist", examples = inputs[!file.exists(inputs)])

# ========================================================================
# SECTION 2: LOAD
# ========================================================================

# --- 2a. R tables: trees, tags and codes -----------------------------------------
rtable_files <- list.files(RTABLES_folder, pattern = sprintf("^%s\\.stem[0-9]+\\.Rdata$", site), full.names = TRUE)
rtable_census <- as.integer(sub(sprintf("^%s\\.stem([0-9]+)\\.Rdata$", site), "\\1", basename(rtable_files)))
rtable_files <- rtable_files[order(rtable_census)]
bio_check(length(rtable_files) > 0L, "R tables found in DATA/RTABLES (run 2_create_R_tables_BCI.R first)")
basename(rtable_files) # bci.stem1.Rdata ... bci.stem9.Rdata

# One row per tree and table. treeID and tag are text in the R tables (tags
# keep their leading zeros).
rtable_trees <- rbindlist(lapply(rtable_files, function(f) {
  e <- new.env()
  obj <- load(f, envir = e)
  x <- as.data.table(get(obj[1L], envir = e))
  unique(x[, .(treeID = as.character(treeID), tag = as.character(tag), sp = as.character(sp))])[, table := basename(f)]
}))
bio_check(
  rtable_trees[, !anyNA(treeID) && !anyNA(tag) && !anyNA(sp)],
  "every row of the R tables has a treeID, a tag and a code"
)
bio_check(
  rtable_trees[, .N, by = treeID][, all(N == length(rtable_files))],
  "every R table lists every tree"
)

trees_new <- unique(rtable_trees[, .(treeID, tag, sp)])
rm(rtable_trees)
bio_check(
  !anyDuplicated(trees_new$treeID),
  "each tree of the R tables has one tag and one code in every table",
  examples = trees_new[duplicated(treeID), treeID], n_bad = sum(duplicated(trees_new$treeID))
)
bio_check(!anyDuplicated(trees_new$tag), "one tree per tag in the R tables", examples = trees_new[duplicated(tag), tag])
trees_new[, treeID_int := suppressWarnings(as.integer(treeID))]
bio_check(!anyNA(trees_new$treeID_int), "every treeID of the R tables is a whole number (to match Condit's)")
nrow(trees_new) # 460,211 trees (2026-10-15)

# --- 2b. Species table of stage 1 ----------------------------------------------
# Object bci.spptable, one row per code of the plot except uniden, all text.
e <- new.env()
obj <- load(SPP_STAGE1_FILE, envir = e)
spp_stage1 <- as.data.table(get(obj[1L], envir = e))
rm(e, obj)

stage1_cols <- c(
  "Mnemonic", "Order", "Family", "Genus", "SpeciesName", "InfraspecificRank", "InfraspecificEpithet",
  "Authority", "Synonyms", "Lifeform_RFoster", "Lifeform_RPerez_SAguilar", "CommonName", "Herbarium", "Notes"
)
bio_check(
  identical(names(spp_stage1), stage1_cols),
  "the stage-1 species table has the expected columns, in order",
  examples = c(setdiff(names(spp_stage1), stage1_cols), setdiff(stage1_cols, names(spp_stage1)))
)
bio_check(all(vapply(spp_stage1, is.character, NA)), "every column of the stage-1 species table is text")
bio_check(
  !anyNA(spp_stage1$Mnemonic) && !anyDuplicated(spp_stage1$Mnemonic),
  "one row per code in the stage-1 species table"
)
bio_check(
  spp_stage1[, all(clean_text(Genus) != "" & clean_text(SpeciesName) != "")],
  "every code of the stage-1 species table has a genus and an epithet"
)
nrow(spp_stage1) # 327 codes (2026-10-15)

# --- 2c. Condit's species table -----------------------------------------------
# It lists the species of all of Condit's Panama plots (1,412 codes), so the
# codes used in the 50-ha plot are taken from his plot data (2d).
bio_check(
  digest::digest(file = CONDIT_SPPTABLE_FILE, algo = "sha256") == CONDIT_SPPTABLE_SHA256,
  "Condit's species table is the file of Dryad version 3 (SHA-256)"
)
e <- new.env()
load(CONDIT_SPPTABLE_FILE, envir = e)
condit_spp <- as.data.table(e$bci.spptable)
rm(e)
bio_check(!anyNA(condit_spp$sp) && !anyDuplicated(condit_spp$sp), "one row per code in Condit's species table")
# Latin = Genus + Species + subsp, where subsp holds the infraspecific name
# with its rank ("var. grandiflora"). IDLevel is "species" unless the code is
# identified only to genus or family, or is a group of species.
bio_check(
  condit_spp[, all(full_name(Genus, Species, subsp) == clean_text(Latin))],
  "Condit's Latin = Genus Species subsp for every code"
)

# --- 2d. Condit's plot data ----------------------------------------------------
# The 8 stem tables of the Dryad release (censuses 1 to 8, 1982-2015) in one
# table; each has a row for every stem, so every tree of the release is here.
condit_plot <- as.data.table(readRDS(CONDIT_PLOT_FILE))
bio_check(
  nrow(condit_plot) == 8L * uniqueN(condit_plot$stemID),
  "Condit's plot data: every stem once in each of the 8 censuses"
)
trees_condit <- unique(condit_plot[, .(treeID = as.integer(treeID), tag = as.character(tag), sp = as.character(sp))])
rm(condit_plot)
invisible(gc())
bio_check(
  nrow(trees_condit) == CONDIT_N_TREES && !anyDuplicated(trees_condit$treeID),
  "Condit's plot data: 423,617 trees, each with one tag and one code"
)
bio_check(!anyDuplicated(trees_condit$tag), "Condit's plot data: one tree per tag")
codes_condit <- sort(unique(trees_condit$sp))
bio_check(all(codes_condit %in% condit_spp$sp), "every code of Condit's plot data is in his species table")

# ========================================================================
# SECTION 3: CODES OF THE R TABLES
# ========================================================================
codes_rtables <- sort(unique(trees_new$sp))
length(codes_rtables) # 328 (2026-10-15)

# The stage-1 table holds exactly the codes of the plot, except uniden
codes_without_row <- setdiff(codes_rtables, c(spp_stage1$Mnemonic, UNIDEN_CODE))
bio_check(
  length(codes_without_row) == 0L,
  "every code of the R tables except uniden has a row in the stage-1 species table",
  examples = codes_without_row, n_bad = length(codes_without_row)
)
codes_without_trees <- setdiff(spp_stage1$Mnemonic, codes_rtables)
bio_check(
  length(codes_without_trees) == 0L,
  "every code of the stage-1 species table has trees in the R tables",
  examples = codes_without_trees, n_bad = length(codes_without_trees)
)

# The codes of the species table
codes_table <- if (INCLUDE_UNIDEN) codes_rtables else setdiff(codes_rtables, UNIDEN_CODE)
length(codes_table) # 328 with uniden (2026-10-15)

# Codes against Condit's release: used by his plot data and not by the R
# tables (their trees now have other codes or left the database), and the
# reverse. 2026-10-15: gone ingaun, pterof, tremsp; new cupasc, ingave, pse2sp.
codes_gone <- setdiff(codes_condit, codes_rtables)
codes_new <- setdiff(codes_rtables, codes_condit)
cat("Codes of Condit et al. (2019) no longer used:", show_codes(codes_gone), "\n")
cat("Codes not used in Condit et al. (2019):", show_codes(codes_new), "\n")

# A code missing from Condit's species table has no former name to compare
codes_not_in_condit <- setdiff(codes_table, condit_spp$sp)
if (length(codes_not_in_condit)) {
  note(paste("codes not in Condit's species table (no former name to compare):", show_codes(codes_not_in_condit)))
}

# ========================================================================
# SECTION 4: TREES WHOSE CODE CHANGED SINCE CONDIT ET AL. (2019)
# ========================================================================
# Condit's treeID is the TreeID of the database, so the releases are matched
# on it. Matched on the tag instead they pair the same trees (checked), so a
# tree never changed treeID or tag: the only change of identity of a tree is
# a change of code.
trees <- merge(
  trees_condit[, .(treeID_int = treeID, tag_condit = tag, sp_condit = sp)],
  trees_new[, .(treeID_int, treeID, tag_new = tag, sp_new = sp)],
  by = "treeID_int", all = TRUE
)
trees[, in_condit := !is.na(sp_condit)]
trees[, in_new := !is.na(sp_new)]

bio_check(
  trees[in_condit & in_new, all(tag_condit == tag_new)],
  "same tag in both releases for every tree matched on treeID",
  examples = trees[in_condit & in_new & tag_condit != tag_new, treeID_int]
)
tags_both <- merge(
  trees_condit[, .(tag, treeID_condit = treeID)],
  trees_new[, .(tag, treeID_new = treeID_int)],
  by = "tag"
)
bio_check(
  tags_both[, all(treeID_condit == treeID_new)],
  "every tag in both releases belongs to the same tree (treeID) in both",
  examples = tags_both[treeID_condit != treeID_new, tag], n_bad = tags_both[treeID_condit != treeID_new, .N]
)
rm(tags_both)

# 2026-10-15: 423,610 trees in both releases; 7 of Condit's no longer in the
# database (hassfl, alsebl, rinosy 2, calolo, protco, picrla); 36,601 recruited
# after his release (first measured in 2022-23)
trees[in_condit & in_new, .N]
trees[in_condit & !in_new, .(treeID = treeID_int, tag = tag_condit, sp = sp_condit)]
trees[!in_condit & in_new, .N]

# Trees whose code changed, and the pairs of codes (old -> new)
changed <- trees[in_condit & in_new & sp_condit != sp_new]
pairs <- changed[, .N, by = .(sp_condit, sp_new)]
setorder(pairs, -N, sp_condit, sp_new)
nrow(changed) # 302 trees (2026-10-15)
nrow(pairs) # 111 pairs (2026-10-15)
# 2026-10-15: tremsp -> tremmi 84, ingaun -> uniden 47, swars1 -> swars2 26,
# swars2 -> swars1 10, pterof -> pterro 9 (by 2_prepare_viewfulltable.R), the
# other pairs 1-4 trees re-identified one by one
pairs[N > 4L]

# Trees per code: trees now = trees in Condit's release - moved to another
# code + from another code - left the database + recruited after his release
counts <- data.table(sp = sort(unique(c(codes_condit, codes_rtables))))
count_cols <- c("n_condit", "n_new", "n_out", "n_in", "n_deleted", "n_recruited")
counts[, (count_cols) := 0L]
counts[trees_condit[, .N, by = sp], n_condit := i.N, on = "sp"]
counts[trees_new[, .N, by = sp], n_new := i.N, on = "sp"]
counts[changed[, .N, by = sp_condit], n_out := i.N, on = .(sp = sp_condit)]
counts[changed[, .N, by = sp_new], n_in := i.N, on = .(sp = sp_new)]
counts[trees[in_condit & !in_new, .N, by = sp_condit], n_deleted := i.N, on = .(sp = sp_condit)]
counts[trees[!in_condit & in_new, .N, by = sp_new], n_recruited := i.N, on = .(sp = sp_new)]
bio_check(
  counts[, all(n_new == n_condit - n_out + n_in - n_deleted + n_recruited)],
  "the trees of every code add up between the two releases",
  examples = counts[n_new != n_condit - n_out + n_in - n_deleted + n_recruited, sp]
)

# One row per tree whose code changed (names added in section 7)
species_code_changes <- changed[, .(treeID, tag = tag_new, spCondit2019 = sp_condit, sp = sp_new)]
setorder(species_code_changes, spCondit2019, sp, tag)

# ========================================================================
# SECTION 5: NAMES COMPARED WITH CONDIT ET AL. (2019)
# ========================================================================
spp <- copy(spp_stage1)
setnames(spp, "Mnemonic", "sp")
if (INCLUDE_UNIDEN && !UNIDEN_CODE %in% spp$sp) {
  spp <- rbind(spp, data.table(sp = UNIDEN_CODE), fill = TRUE)
}
if (!INCLUDE_UNIDEN) spp <- spp[sp != UNIDEN_CODE]
# In order of code, with uniden (not a taxon) as the last row
setorder(spp, sp)
spp <- rbind(spp[sp != UNIDEN_CODE], spp[sp == UNIDEN_CODE])
bio_check(
  setequal(spp$sp, codes_table) && !anyDuplicated(spp$sp),
  "one row per code of the R tables in the species table (uniden as set by INCLUDE_UNIDEN)"
)
bio_check(
  spp[sp != UNIDEN_CODE, identical(sp, sort(sp, method = "radix"))] && (!INCLUDE_UNIDEN || spp[.N, sp] == UNIDEN_CODE),
  "the codes are in order, and uniden is the last row"
)

# Latin: genus, epithet and the infraspecific name when there is one
spp[, Latin := full_name(Genus, SpeciesName, InfraspecificRank, InfraspecificEpithet)]
spp[Latin == "", Latin := NA_character_]
bio_check(spp[is.na(Latin), all(sp == UNIDEN_CODE)], "every code except uniden has a Latin name")

# Names matched by code, never by Latin name; uniden has no name to compare
names_cmp <- merge(
  spp[sp != UNIDEN_CODE, .(
    sp,
    genus_new = Genus, epithet_new = SpeciesName,
    infra_new = full_name(InfraspecificRank, InfraspecificEpithet), latin_new = Latin
  )],
  condit_spp[, .(sp, genus_old = Genus, epithet_old = Species, infra_old = subsp, latin_old = clean_text(Latin), IDLevel)],
  by = "sp"
)
names_cmp[, NameChange := name_change_kind(
  genus_old, epithet_old, infra_old, genus_new, epithet_new, infra_new,
  id_level_old = IDLevel
)]
bio_check(
  setequal(names_cmp[!is.na(NameChange), sp], names_cmp[latin_new != latin_old, sp]),
  "the codes whose genus, epithet or infraspecific name changed are those whose full name changed"
)

# 2026-10-15: 49 codes: epithet 29, genus 6, genus and epithet 6, epithet
# spelling 3 (appuse, ardiba, rinosy), infraspecific name 2 (myrcga, protte),
# epithet gender 1 (tremmi), genus and epithet gender 1 (cappfr),
# morphospecies label 1 (sapisp)
names_cmp[!is.na(NameChange), .N]
names_cmp[!is.na(NameChange), .N, by = NameChange][order(-N)]
names_cmp[!is.na(NameChange) & NameChange != "epithet", .(sp, latin_old, latin_new, NameChange)][order(NameChange, sp)]

spp[, c("NameCondit2019", "NameChange") := NA_character_]
spp[names_cmp[!is.na(NameChange)], `:=`(NameCondit2019 = i.latin_old, NameChange = i.NameChange), on = "sp"]

# ========================================================================
# SECTION 6: NOTES
# ========================================================================

# --- 6a. Re-identified trees, both ways, with their counts ----------------------
spp[, c("ReidentifiedFrom", "ReidentifiedTo") := NA_character_]
spp[pairs[, .(x = format_code_counts(sp_condit, N)), by = sp_new], ReidentifiedFrom := i.x, on = .(sp = sp_new)]
spp[pairs[, .(x = format_code_counts(sp_new, N)), by = sp_condit], ReidentifiedTo := i.x, on = .(sp = sp_condit)]

# --- 6b. Trees of the codes no longer used ---------------------------------------
# A code no longer used has no row, so its trees are described in the notes
# of the codes that hold them now. 2026-10-15: tremsp -> tremmi (84),
# ingaun -> uniden (47), pterof -> pterro (9).
absorbed <- pairs[sp_condit %in% codes_gone]
absorbed[condit_spp, name_old := clean_text(i.Latin), on = .(sp_condit = sp)]
absorbed[counts, n_old := i.n_condit, on = .(sp_condit = sp)]
absorbed[, sentence := absorbed_sentence(sp_condit, name_old, N, n_old, reference = CONDIT_REF)]
setorder(absorbed, sp_new, -N, sp_condit)
absorbed[, .(sp_condit, sp_new, N, n_old)]

# Codes no longer used whose trees all left the database, or went to a code
# without a row (uniden when INCLUDE_UNIDEN = FALSE), have no note: they are
# only in the files of the changes of code (none on 2026-10-15)
codes_gone_unplaced <- c(setdiff(codes_gone, absorbed$sp_condit), absorbed[!sp_new %in% spp$sp, sp_condit])
if (length(codes_gone_unplaced)) {
  note(paste("codes no longer used without a note in the species table:", show_codes(unique(codes_gone_unplaced))))
}

# --- 6c. Notes --------------------------------------------------------------------
# The note of the stage-1 table first, then: uniden is not a taxon; the former
# name; a code not used in Condit et al. (2019); the trees of codes no longer
# used.
spp[, note_uniden := fifelse(sp == UNIDEN_CODE, "Stems not identified to species; not a taxon.", "")]
spp[, note_name := fifelse(is.na(NameCondit2019), "", paste0("Named ", NameCondit2019, " in ", CONDIT_REF, "."))]
spp[, note_new := fifelse(sp %in% codes_new, paste0("Code not used in ", CONDIT_REF, "."), "")]
spp[, note_absorbed := ""]
spp[absorbed[, .(x = paste(sentence, collapse = " ")), by = sp_new], note_absorbed := i.x, on = .(sp = sp_new)]
spp[, Notes := join_notes(Notes, note_uniden, note_name, note_new, note_absorbed)]
spp[, c("note_uniden", "note_name", "note_new", "note_absorbed") := NULL]

absorbed_rows <- merge(absorbed[, .(sp = sp_new, sp_condit)], spp[, .(sp, Notes)], by = "sp")
bio_check(
  absorbed_rows[, all(mapply(grepl, paste0("recorded as ", sp_condit, " ("), Notes, fixed = TRUE))],
  "every code no longer used is named in the notes of the codes that hold its trees"
)
rm(absorbed_rows)

# The codes with a note (2026-10-15: 55)
spp[!is.na(Notes), .(sp, Notes)]

# ========================================================================
# SECTION 7: ASSEMBLE AND CHECK
# ========================================================================
out_cols <- c(
  "sp", "Order", "Family", "Genus", "SpeciesName", "Latin", "InfraspecificRank", "InfraspecificEpithet",
  "Authority", "Synonyms", "Lifeform_RFoster", "Lifeform_RPerez_SAguilar", "CommonName", "Herbarium",
  "NameCondit2019", "NameChange", "ReidentifiedFrom", "ReidentifiedTo", "Notes"
)
spptable <- spp[, ..out_cols]

# The columns of the stage-1 table are copied unchanged, and its notes stay at
# the start of Notes
copied_cols <- setdiff(stage1_cols, c("Mnemonic", "Notes"))
copy_cmp <- merge(spp_stage1, spptable, by.x = "Mnemonic", by.y = "sp", suffixes = c("_stage1", ""))
bio_check(
  nrow(copy_cmp) == nrow(spp_stage1) &&
    all(vapply(copied_cols, function(k) identical(copy_cmp[[paste0(k, "_stage1")]], copy_cmp[[k]]), NA)),
  "every code of the stage-1 table is in the species table with its columns unchanged"
)
bio_check(
  copy_cmp[!is.na(Notes_stage1), all(startsWith(Notes, clean_text(Notes_stage1)))],
  "the notes of the stage-1 table open the Notes of their codes"
)
rm(copy_cmp)

# The uniden row has no names: only its re-identified trees and its notes
if (INCLUDE_UNIDEN) {
  bio_check(
    spptable[sp == UNIDEN_CODE, all(is.na(unlist(.SD))), .SDcols = c(copied_cols, "Latin", "NameCondit2019", "NameChange")],
    "the uniden row has no names, only re-identified trees and notes"
  )
}

# The text built here has no empty strings, "NA" or outer spaces
built_cols <- c("Latin", "NameCondit2019", "NameChange", "ReidentifiedFrom", "ReidentifiedTo", "Notes")
built_text <- unlist(spptable[, ..built_cols], use.names = FALSE)
built_text <- built_text[!is.na(built_text)]
bio_check(
  !any(built_text == "" | built_text == "NA" | built_text != trimws(built_text) | grepl("  ", built_text, fixed = TRUE)),
  "no empty text, \"NA\" text, outer or double spaces in the columns built here"
)

# Names of the trees whose code changed: Condit's name of the former code and
# the Latin of the code now
species_code_changes[condit_spp, LatinCondit2019 := clean_text(i.Latin), on = .(spCondit2019 = sp)]
species_code_changes[spptable, Latin := i.Latin, on = "sp"]
setcolorder(species_code_changes, c("treeID", "tag", "spCondit2019", "LatinCondit2019", "sp", "Latin"))
bio_check(
  nrow(species_code_changes) == nrow(changed) && !anyDuplicated(species_code_changes$treeID),
  "one row per tree whose code changed"
)
bio_check(
  !anyNA(species_code_changes$LatinCondit2019),
  "every former code has its name in Condit's species table"
)
bio_check(
  species_code_changes[is.na(Latin), all(sp == UNIDEN_CODE)],
  "every code now has a Latin name, except uniden"
)

# One row per change of code (former code -> code now) with its number of
# trees, the largest changes first. 2026-10-15: 111 changes, 302 trees.
species_code_changes_summary <- species_code_changes[, .(trees = .N), by = .(spCondit2019, LatinCondit2019, sp, Latin)]
setorder(species_code_changes_summary, -trees, spCondit2019, sp)
bio_check(
  nrow(species_code_changes_summary) == nrow(pairs) &&
    !anyDuplicated(species_code_changes_summary[, .(spCondit2019, sp)]) &&
    sum(species_code_changes_summary$trees) == nrow(species_code_changes),
  "one row per change of code, and the trees of the changes add up to the trees whose code changed"
)
species_code_changes_summary[trees > 4L]

# ========================================================================
# SECTION 8: EXPORT
# ========================================================================
# The species table and the changes of code (one row per change, one row per
# tree), as .rdata (a data.frame, the ForestGEO convention) and .csv
spptable_name <- paste0(site, ".spptable")
summary_name <- paste0(site, ".species_code_changes_summary")
changes_name <- paste0(site, ".species_code_changes")
assign(spptable_name, as.data.frame(spptable))
assign(summary_name, as.data.frame(species_code_changes_summary))
assign(changes_name, as.data.frame(species_code_changes))
save(list = spptable_name, file = file.path(RTABLES_folder, paste0(spptable_name, ".rdata")))
save(list = summary_name, file = file.path(RTABLES_folder, paste0(summary_name, ".rdata")))
save(list = changes_name, file = file.path(RTABLES_folder, paste0(changes_name, ".rdata")))
# The .csv files are UTF-8 with a byte order mark (BOM), so that Excel shows
# the Spanish characters (árbol, Pérez) correctly; R reads them with
# fread() or read.csv(fileEncoding = "UTF-8-BOM"). fwrite() never writes a
# header in append mode (its documented behaviour, regardless of col.names),
# so the header is written by a normal fwrite() first, and the BOM is
# prepended to the finished file, header included.
write_csv_utf8 <- function(x, file) {
  x <- copy(x)
  for (col in names(x)[vapply(x, is.character, NA)]) set(x, j = col, value = enc2utf8(x[[col]]))
  tmp <- tempfile()
  fwrite(x, tmp)
  body <- readBin(tmp, "raw", file.size(tmp))
  file.remove(tmp)
  writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), body), file)
}
write_csv_utf8(spptable, file.path(RTABLES_folder, paste0(spptable_name, ".csv")))
write_csv_utf8(species_code_changes_summary, file.path(RTABLES_folder, paste0(summary_name, ".csv")))
write_csv_utf8(species_code_changes, file.path(RTABLES_folder, paste0(changes_name, ".csv")))

# The bytes: BOM first, valid UTF-8, header present, and the Spanish
# characters of the table written as UTF-8 (á = C3 A1), not as Latin-1 (E1)
csv_cols <- list(spptable = out_cols, species_code_changes_summary = names(species_code_changes_summary), species_code_changes = names(species_code_changes))
for (nm in c(spptable_name, summary_name, changes_name)) {
  csv_file <- file.path(RTABLES_folder, paste0(nm, ".csv"))
  bytes <- readBin(csv_file, "raw", file.size(csv_file))
  text <- rawToChar(bytes[-(1:3)])
  Encoding(text) <- "UTF-8"
  header <- strsplit(sub("\n.*", "", text), ",", fixed = TRUE)[[1]]
  bio_check(
    identical(bytes[1:3], as.raw(c(0xEF, 0xBB, 0xBF))) && validUTF8(text),
    paste0(nm, ".csv is UTF-8 with a BOM")
  )
  bio_check(
    identical(header, csv_cols[[sub(paste0("^", site, "\\."), "", nm)]]),
    paste0(nm, ".csv has its column names as the first row")
  )
}
spp_bytes <- readBin(file.path(RTABLES_folder, paste0(spptable_name, ".csv")), "raw", file.size(file.path(RTABLES_folder, paste0(spptable_name, ".csv"))))
spp_text <- rawToChar(spp_bytes[-(1:3)])
Encoding(spp_text) <- "UTF-8"
bio_check(
  !grepl("\u00c3|\ufffd", spp_text) && grepl("\u00e1rbol", spp_text),
  "the Spanish characters of bci.spptable.csv are written correctly (árbol)"
)
rm(bytes, text, header, csv_file, csv_cols, spp_bytes, spp_text)

# Read back: the .rdata files hold the same objects, and the .csv files the
# same text (tags read as text, empty = NA; the tree counts are compared as
# text too)
for (nm in c(spptable_name, summary_name, changes_name)) {
  e <- new.env()
  loaded <- load(file.path(RTABLES_folder, paste0(nm, ".rdata")), envir = e)
  bio_check(
    identical(loaded, nm) && identical(e[[nm]], get(nm)),
    paste0(nm, ".rdata reads back identical")
  )
  csv_back <- fread(file.path(RTABLES_folder, paste0(nm, ".csv")), colClasses = "character", na.strings = "", encoding = "UTF-8")
  as_text <- as.data.frame(lapply(get(nm), as.character), stringsAsFactors = FALSE)
  bio_check(
    isTRUE(all.equal(as.data.frame(csv_back), as_text, check.attributes = FALSE)),
    paste0(nm, ".csv reads back identical (as text)")
  )
}
rm(e, loaded, csv_back, as_text)

cat(
  "\nSummary\n",
  " species table:", nrow(spptable), "codes", if (INCLUDE_UNIDEN) "(uniden included, last row)" else "(uniden left out)", "\n",
  " codes whose name changed since Condit et al. (2019):", spptable[!is.na(NameChange), .N], "\n",
  " codes with re-identified trees:", spptable[!is.na(ReidentifiedFrom) | !is.na(ReidentifiedTo), .N], "\n",
  " codes with a note:", spptable[!is.na(Notes), .N], "\n",
  " changes of code:", nrow(species_code_changes_summary), "| trees whose code changed:", nrow(species_code_changes), "\n",
  " written to", RTABLES_folder, "\n"
)
