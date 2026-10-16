# ========================================================================
# TESTS: species table of the R tables (3_create_species_table_BCI.R)
# ========================================================================
# Run from the project root (or from this folder):
#   testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_species_table.R")
#
#   1. The functions of species_table_functions.R on hand-written cases
#   2. Real data: bci.spptable, bci.species_code_changes_summary and
#      bci.species_code_changes, recomputed from the R tables and Condit's
#      files with separate code (skipped when the outputs are missing)
#
# Option (environment variable):
#   SPECIES_RTABLES_DIR  folder with the R tables and the three outputs
#                        (default DATA/RTABLES)
# ========================================================================
suppressMessages({
  library(testthat)
  library(data.table)
})

find_stage3_file <- function(f) {
  cand <- c(f, file.path("..", f), file.path("BCI_stem_reconstruction", "3_PREPARE_R_TABLES", f))
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("cannot find ", f, ": run from the project root or from the tests folder")
  normalizePath(hit)
}
stage3_dir <- dirname(find_stage3_file("species_table_functions.R"))
project_root <- normalizePath(file.path(stage3_dir, "..", ".."))
source(file.path(stage3_dir, "species_table_functions.R"))

# ------------------------------------------------------------------------
# 1. Functions
# ------------------------------------------------------------------------
test_that("clean_text and full_name ignore NA, outer and repeated spaces", {
  expect_identical(clean_text(c(NA, "  a  b ", "c")), c("", "a b", "c"))
  expect_identical(full_name("Swartzia", "simplex", "var.", "grandiflora"), "Swartzia simplex var. grandiflora")
  expect_identical(full_name("Ficus", " trigonata ", NA, ""), "Ficus trigonata")
  expect_identical(full_name(NA, NA), "")
  expect_identical(full_name(c("A", "B"), c("x", NA)), c("A x", "B"))
})

test_that("name_change_kind: one case per kind, with the real names", {
  kind <- function(old, new, id = "species") {
    name_change_kind(old[1], old[2], old[3], new[1], new[2], new[3], id_level_old = id)
  }
  expect_identical(kind(c("Swartzia", "simplex", "var. grandiflora"), c("Swartzia", "simplex", "var.  grandiflora")), NA_character_)
  expect_identical(kind(c("Sapium", "broadleaf", NA), c("Sapium", "sp.", NA), id = "genus"), "morphospecies label")
  expect_identical(kind(c("Protium", "tenuifolium", "subsp. sessiliflorum"), c("Protium", "tenuifolium", NA)), "infraspecific name")
  expect_identical(kind(c("Clidemia", "dentata", NA), c("Miconia", "dentata", NA)), "genus")
  expect_identical(kind(c("Capparidastrum", "frondosum", NA), c("Morisonia", "frondosa", NA)), "genus and epithet gender")
  expect_identical(kind(c("Laetia", "procera", NA), c("Casearia", "bicolor", NA)), "genus and epithet")
  expect_identical(kind(c("Appunia", "seibertii", NA), c("Appunia", "siebertii", NA)), "epithet spelling")
  expect_identical(kind(c("Ardisia", "bartlettii", NA), c("Ardisia", "bartletii", NA)), "epithet spelling")
  expect_identical(kind(c("Rinorea", "sylvatica", NA), c("Rinorea", "silvatica", NA)), "epithet spelling")
  expect_identical(kind(c("Trema", "micrantha", NA), c("Trema", "micranthum", NA)), "epithet gender")
  expect_identical(kind(c("Mosannona", "garwoodii", NA), c("Mosannona", "garwoodiae", NA)), "epithet gender")
  expect_identical(kind(c("Ficus", "crocata", NA), c("Ficus", "trigonata", NA)), "epithet")
  # another genus and a small change that is not a gender ending
  expect_identical(kind(c("Xus", "seibertii", NA), c("Yus", "siebertii", NA)), "genus and epithet")
})

test_that("latin_stem: epithets that differ only in their gender ending share a stem", {
  same <- function(a, b) latin_stem(a) == latin_stem(b)
  expect_true(same("micrantha", "micranthum"))
  expect_true(same("frondosum", "frondosa"))
  expect_true(same("garwoodii", "garwoodiae"))
  expect_true(same("rubra", "rubrum") && same("ruber", "rubra"))
  expect_true(same("aspera", "asperum") && same("asper", "aspera"))
  expect_true(same("nobilis", "nobile") && same("costaricensis", "costaricense"))
  expect_false(same("seibertii", "siebertii"))
  expect_false(same("bartlettii", "bartletii"))
  expect_false(same("sylvatica", "silvatica"))
  expect_false(same("trinervium", "trinerve"))
})

test_that("name_change_kind is vectorised, recycles IDLevel and takes no rows", {
  k <- name_change_kind(
    c("Trema", "Ficus"), c("micrantha", "crocata"), c(NA, NA),
    c("Trema", "Ficus"), c("micranthum", "crocata"), c(NA, NA)
  )
  expect_identical(k, c("epithet gender", NA))
  expect_identical(name_change_kind(character(0), character(0), character(0), character(0), character(0), character(0)), character(0))
  # A morphospecies whose label did not change has no change
  expect_identical(name_change_kind("Nectandra", "fuzzy", NA, "Nectandra", "fuzzy", NA, id_level_old = "genus"), NA_character_)
})

test_that("format_code_counts: largest group first, ties by code, NA when empty", {
  expect_identical(format_code_counts(c("protpa", "tet2pa", "swars2"), c(2L, 2L, 26L)), "swars2 (26), protpa (2), tet2pa (2)")
  expect_identical(format_code_counts(character(0), integer(0)), NA_character_)
})

test_that("absorbed_sentence: all the trees, one tree, part of the trees", {
  expect_identical(
    absorbed_sentence("tremsp", "Trema unidentified", 84L, 84L),
    "Includes the 84 trees recorded as tremsp (Trema unidentified) in Condit et al. (2019), a code no longer used."
  )
  expect_identical(
    absorbed_sentence("xxxxxx", "Xus yus", 1L, 1L),
    "Includes the tree recorded as xxxxxx (Xus yus) in Condit et al. (2019), a code no longer used."
  )
  expect_identical(
    absorbed_sentence(c("a", "b"), c("A a", "B b"), c(2L, 3L), c(5L, 3L), reference = "R"),
    c(
      "Includes 2 of the 5 trees recorded as a (A a) in R, a code no longer used.",
      "Includes the 3 trees recorded as b (B b) in R, a code no longer used."
    )
  )
})

test_that("join_notes: full stops, empty parts dropped, NA when nothing", {
  out <- join_notes(
    c("No longer present: dead in 1990", NA, "", "Is it?"),
    c("Named X in R.", NA, " ", ""),
    c("", NA, "", "Code new.")
  )
  expect_identical(out, c("No longer present: dead in 1990. Named X in R.", NA, NA, "Is it? Code new."))
})

# ------------------------------------------------------------------------
# 2. Real data
# ------------------------------------------------------------------------
data_dir <- file.path(project_root, "BCI_stem_reconstruction", "DATA")
rtables_dir <- Sys.getenv("SPECIES_RTABLES_DIR", file.path(data_dir, "RTABLES"))
spptable_file <- file.path(rtables_dir, "bci.spptable.rdata")
summary_file <- file.path(rtables_dir, "bci.species_code_changes_summary.rdata")
changes_file <- file.path(rtables_dir, "bci.species_code_changes.rdata")
stem_files <- list.files(rtables_dir, pattern = "^bci\\.stem[0-9]+\\.Rdata$", full.names = TRUE)
real_ok <- file.exists(spptable_file) && file.exists(summary_file) && file.exists(changes_file) && length(stem_files) > 0L &&
  file.exists(file.path(data_dir, "SPP_TABLE", "condit", "bci.spptable.rdata")) &&
  file.exists(file.path(data_dir, "RAW", "bci_dryad_condit.rds")) &&
  file.exists(file.path(data_dir, "SPP_TABLE", "bci_spptable.RData"))

if (real_ok) {
  load_one <- function(f) {
    e <- new.env()
    get(load(f, envir = e)[1L], envir = e)
  }
  spp <- as.data.table(load_one(spptable_file))
  smr <- as.data.table(load_one(summary_file))
  chg <- as.data.table(load_one(changes_file))
  stage1 <- as.data.table(load_one(file.path(data_dir, "SPP_TABLE", "bci_spptable.RData")))
  condit_spp <- as.data.table(load_one(file.path(data_dir, "SPP_TABLE", "condit", "bci.spptable.rdata")))
  # Trees of the R tables (the last table lists every tree) and of Condit's release
  last_table <- as.data.table(load_one(stem_files[which.max(as.integer(gsub("\\D", "", basename(stem_files))))]))
  new_trees <- unique(last_table[, .(treeID = as.integer(treeID), tag = as.character(tag), sp = as.character(sp))])
  rm(last_table)
  condit_trees <- unique(as.data.table(readRDS(file.path(data_dir, "RAW", "bci_dryad_condit.rds")))[
    , .(treeID = as.integer(treeID), tag = as.character(tag), sp = as.character(sp))
  ])
}

test_that("real data: one row per code of the R tables, the stage-1 columns unchanged", {
  skip_if_not(real_ok, "no R tables / species-table outputs found")
  expect_false(anyDuplicated(spp$sp) > 0)
  expect_setequal(spp$sp, unique(new_trees$sp))
  # in order of code, uniden (not a taxon) last
  expect_identical(spp$sp, c(sort(setdiff(spp$sp, "uniden"), method = "radix"), "uniden"))
  m <- spp[match(stage1$Mnemonic, sp)]
  for (k in setdiff(names(stage1), c("Mnemonic", "Notes"))) expect_identical(m[[k]], stage1[[k]], info = k)
  # Latin is genus, epithet and infraspecific name; only uniden has none
  latin <- trimws(gsub("\\s+", " ", paste(
    spp$Genus, spp$SpeciesName,
    fifelse(is.na(spp$InfraspecificRank), "", spp$InfraspecificRank),
    fifelse(is.na(spp$InfraspecificEpithet), "", spp$InfraspecificEpithet)
  )))
  expect_identical(spp[sp != "uniden", Latin], latin[spp$sp != "uniden"])
  expect_identical(spp[is.na(Latin), sp], "uniden")
})

test_that("real data: the trees whose code changed, recomputed tree by tree", {
  skip_if_not(real_ok, "no R tables / species-table outputs found")
  both <- merge(condit_trees, new_trees, by = "treeID", suffixes = c("_old", "_new"))
  expect_true(all(both$tag_old == both$tag_new))
  moved <- both[sp_old != sp_new]
  expect_identical(nrow(chg), nrow(moved))
  m <- moved[match(as.integer(chg$treeID), treeID)]
  expect_false(anyNA(m$treeID))
  expect_identical(chg$tag, m$tag_new)
  expect_identical(chg$spCondit2019, m$sp_old)
  expect_identical(chg$sp, m$sp_new)
  expect_identical(chg$LatinCondit2019, condit_spp[match(chg$spCondit2019, sp), trimws(gsub("\\s+", " ", Latin))])
  expect_identical(chg$Latin, spp[match(chg$sp, sp), Latin])
})

test_that("real data: one row per change of code, with its number of trees", {
  skip_if_not(real_ok, "no R tables / species-table outputs found")
  both <- merge(condit_trees, new_trees, by = "treeID", suffixes = c("_old", "_new"))
  moved <- both[sp_old != sp_new, .(trees = .N), by = .(spCondit2019 = sp_old, sp = sp_new)]
  expect_identical(nrow(smr), nrow(moved))
  expect_false(anyDuplicated(smr[, .(spCondit2019, sp)]) > 0)
  m <- moved[smr, on = .(spCondit2019, sp)]
  expect_identical(m$trees, m$i.trees)
  expect_identical(sum(smr$trees), nrow(chg))
  # the largest changes first, then by code
  expect_identical(smr[order(-trees, spCondit2019, sp, method = "radix"), .(spCondit2019, sp)], smr[, .(spCondit2019, sp)])
  # the names are those of the trees' file
  expect_identical(smr$LatinCondit2019, chg[match(smr$spCondit2019, spCondit2019), LatinCondit2019])
  expect_identical(smr$Latin, spp[match(smr$sp, sp), Latin])
})

test_that("real data: ReidentifiedFrom / ReidentifiedTo equal the counts of the trees", {
  skip_if_not(real_ok, "no R tables / species-table outputs found")
  fmt <- function(codes) {
    if (!length(codes)) return(NA_character_)
    t <- table(codes)
    o <- order(-as.integer(t), names(t))
    paste0(names(t)[o], " (", as.integer(t)[o], ")", collapse = ", ")
  }
  from <- vapply(spp$sp, function(s) fmt(chg[sp == s, spCondit2019]), "")
  to <- vapply(spp$sp, function(s) fmt(chg[spCondit2019 == s, sp]), "")
  expect_identical(spp$ReidentifiedFrom, unname(from))
  expect_identical(spp$ReidentifiedTo, unname(to))
})

test_that("real data: former names and notes", {
  skip_if_not(real_ok, "no R tables / species-table outputs found")
  cl <- function(x) trimws(gsub("\\s+", " ", ifelse(is.na(x), "", x)))
  old <- condit_spp[match(spp$sp, sp), cl(Latin)]
  changed <- spp$sp != "uniden" & !is.na(old) & cl(spp$Latin) != old
  expect_identical(!is.na(spp$NameCondit2019), changed)
  expect_identical(spp$NameCondit2019[changed], old[changed])
  expect_identical(!is.na(spp$NameChange), changed)
  # Every former name, code new to the plot and code no longer used is in the notes
  expect_true(all(mapply(grepl, paste0("Named ", spp$NameCondit2019[changed], " in Condit et al. (2019)."), spp$Notes[changed], fixed = TRUE)))
  codes_new <- setdiff(spp$sp, condit_trees$sp)
  expect_true(all(grepl("Code not used in Condit et al. (2019).", spp[sp %in% codes_new, Notes], fixed = TRUE)))
  gone <- setdiff(unique(condit_trees$sp), unique(new_trees$sp))
  for (g in intersect(gone, chg$spCondit2019)) {
    holders <- unique(chg[spCondit2019 == g, sp])
    expect_true(all(grepl(paste0("recorded as ", g, " ("), spp[sp %in% holders, Notes], fixed = TRUE)), info = g)
  }
  # Stage-1 notes open the Notes
  s1 <- stage1[!is.na(Notes)]
  expect_true(all(startsWith(spp[match(s1$Mnemonic, sp), Notes], s1$Notes)))
})
