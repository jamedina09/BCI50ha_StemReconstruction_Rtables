# Tests for dp_global/R/measurement_rejoin.R
# Run from the project root:
#   testthat::test_file("BCI_stem_reconstruction/2_STEM_IDENTIFICATION/tests/test_measurement_rejoin.R")
suppressMessages({
  library(testthat)
  library(data.table)
})
cand <- c("dp_global/R/measurement_rejoin.R", "../../../dp_global/R/measurement_rejoin.R")
source(cand[file.exists(cand)][1])

# One tree from rows: list(c(census, StemID, recon stem, DBH mm, HOM, status, codes, pin), ...)
# Unless a test says otherwise, the trees carry a sign of a moved point of
# measurement (rule 8): a HOM away from 1.3 m or a buttress code "B".
tree <- function(rows, tag = "T1", growth_form = "tree") {
  x <- rbindlist(lapply(rows, function(r) {
    data.table(
      Tag = tag, CensusID = as.integer(r[[1]]), StemID = as.character(r[[2]]), ReconstructedStemID = as.numeric(r[[3]]),
      DBH = as.numeric(r[[4]]), HOM = as.numeric(r[[5]]), Status = r[[6]], ListOfTSM = r[[7]], TrueStemID = r[[8]]
    )
  }))
  x[, `:=`(
    single_stem_tags = FALSE, dbh_with_best_candidate_taper_corrected = DBH,
    ExactDate = as.Date("1982-01-01") + round((CensusID - 1L) * 365.25 * 3), obs_row_id = seq_len(.N),
    ReconstructionMethod = "dp", growth_form = growth_form
  )]
  x
}
# run(): apply_measurement_rejoin() with the bounds of the BCI stage-2 run
# stems(): number of measured stems in its result
run <- function(x, ...) apply_measurement_rejoin(x, max_shrink = -0.5, max_growth = 5, recruit_max_mm = 260, verbose = FALSE, ...)
stems <- function(r) r$dt[!is.na(DBH), uniqueN(ReconstructedStemID)]

test_that("buttressed trunk restarting above the recruit limit is joined; StemIDs are not used", {
  # different database StemIDs (renumbered): joined
  x <- tree(list(list(1, "s1", 1, 1775, 1.3, "alive", NA, NA), list(2, "s9", 2, 1140, 1.3, "alive", NA, NA), list(3, "s9", 2, 1134, 6.5, "alive", "B", NA)))
  r <- run(x)
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$route, "impossible_recruit")
  expect_equal(stems(r), 1L)
  expect_equal(r$dt[CensusID == 1, ReconstructionMethod], "measurement_rejoin")
  expect_identical(r$dt[, .(DBH, HOM, Status, StemID, TrueStemID)], x[, .(DBH, HOM, Status, StemID, TrueStemID)]) # nothing else changes
  # the same tree with one database StemID gives exactly the same result
  y <- copy(x)[, StemID := "s1"]
  expect_identical(run(y)$dt$ReconstructedStemID, r$dt$ReconstructedStemID)
})

test_that("a new stem below the recruit limit stays split, even with the same database StemID", {
  r <- run(tree(list(list(1, "s1", 1, 300, 1.3, "alive", NA, NA), list(2, "s1", 2, 200, 1.3, "alive", NA, NA))))
  expect_equal(nrow(r$pairs), 0L)
  expect_equal(stems(r), 2L)
})

test_that("an earlier stem below 10 cm stays split", {
  r <- run(tree(list(list(1, "s1", 1, 90, 1.3, "alive", NA, NA), list(2, "s1", 2, 280, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "earlier stem below min_dbh")
})

test_that("broken trunk (size ratio below 0.4) stays split", {
  r <- run(tree(list(list(1, "s1", 1, 1000, 1.3, "alive", NA, NA), list(2, "s1", 2, 300, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "size ratio below ratio_min")
})

test_that("a new stem more than 1.5 times larger stays split (also an upward recording error)", {
  r <- run(tree(list(list(1, "s1", 1, 250, 1.3, "alive", NA, NA), list(2, "s9", 2, 498, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "size ratio above ratio_max")
  r <- run(tree(list(list(1, "s1", 1, 391, 1.3, "alive", NA, NA), list(2, "s1", 2, 642, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "size ratio above ratio_max")
})

test_that("a break or resprout code on the later measurement keeps the stems apart (rule R1)", {
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 2, 600, 1.3, "broken below", NA, NA))))
  expect_equal(r$candidates$why, "break/resprout code on the later measurement")
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 2, 600, 1.3, "alive", "R", NA))))
  expect_equal(r$candidates$why, "break/resprout code on the later measurement")
  # a code on the earlier measurement only (the stem started as a resprout) does not block
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "broken below", "R", NA), list(2, "s1", 2, 600, 1.3, "alive", "B", NA))))
  expect_equal(nrow(r$pairs), 1L)
})

test_that("a split the engine chose inside the growth bounds is left alone", {
  x <- tree(list(list(1, "s1", 1, 300, 1.3, "alive", NA, NA), list(2, "s1", 2, 290, 1.3, "alive", NA, NA)))
  r <- run(x)
  expect_equal(r$candidates$why, "link inside the growth bounds")
  expect_equal(nrow(r$pairs), 0L)
  expect_equal(nrow(r$obs_pairs), 0L) # not joined in the posterior samples either
  expect_equal(stems(r), 2L)
})

test_that("two different pins keep the stems apart; one pin decides the target id", {
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, "p1"), list(2, "s1", 2, 600, 1.3, "alive", "B", "p2"))))
  expect_equal(r$candidates$why, "two different pins")
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 2, 600, 4, "alive", NA, "p2"), list(3, "s1", 2, 610, 4, "alive", NA, "p2"))))
  expect_equal(r$dt[, unique(ReconstructedStemID)], 2) # the pinned stem keeps its id
})

test_that("no join when the ending stem has a later record or the new stem an earlier one", {
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 1, NA, 1.3, "dead", NA, NA), list(2, "s1", 2, 600, 1.3, "alive", NA, NA))))
  expect_equal(nrow(r$pairs), 0L)
})

test_that("two equally good candidates stay split (no clearly best pairing)", {
  r <- run(tree(list(
    list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(1, "s2", 3, 880, 1.3, "alive", NA, NA),
    list(2, "s9", 2, 700, 1.3, "alive", "B", NA)
  )))
  expect_equal(nrow(r$pairs), 0L)
  expect_equal(r$candidates[, sort(unique(why))], "ambiguous (no clearly best pairing)")
})

test_that("a small stem ending in the same census does not block the join", {
  r <- run(tree(list(
    list(1, "s1", 1, 900, 4.2, "alive", NA, NA), list(1, "s2", 2, 120, 1.3, "alive", NA, NA),
    list(2, "s9", 3, 700, 4.2, "alive", NA, NA)
  )))
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs[, c(rs0, rs1)], c("1", "3"))
  expect_equal(r$dt[!is.na(DBH), uniqueN(ReconstructedStemID)], 2L)
})

test_that("two trunks measured at two heights: each pairs with the stem at its own height", {
  r <- run(tree(list(
    list(1, "a", 1, 748, 5.45, "alive", NA, NA), list(1, "b", 2, 830, 1.3, "alive", NA, NA),
    list(2, "c", 3, 697, 5.45, "alive", NA, NA), list(2, "d", 4, 606, 1.3, "alive", NA, NA)
  )))
  expect_equal(nrow(r$pairs), 2L)
  expect_setequal(r$pairs[, paste(rs0, rs1)], c("1 3", "2 4"))
})

test_that("forks after one trunk stay split, and the trunk is not joined across them", {
  # trunk at 1.3 m, two forks, then a trunk at 2.5 m (HOM raised, smaller)
  r <- run(tree(list(
    list(1, "t", 1, 631, 1.3, "alive", NA, NA),
    list(2, "f1", 2, 500, 1.3, "alive", NA, NA), list(2, "f2", 3, 400, 1.3, "alive", NA, NA),
    list(3, "n", 4, 600, 2.5, "alive", NA, NA)
  )))
  expect_equal(r$candidates[rs1 %in% c("2", "3"), unique(why)], "ambiguous (no clearly best pairing)")
  # the trunk was measured in census 2 (as forks): no join across that census
  expect_equal(r$candidates[rs0 == "1" & rs1 == "4", why], "the trunk may be measured in the census between")
  expect_equal(nrow(r$pairs), 0L)
})

test_that("across a census without a measurement: joined only when the HOM was raised and the DBH is smaller", {
  gap_why <- "across a census without a measurement: no recorded HOM rise with a smaller DBH"
  # 778 mm at 1.3 m, then 670 mm at 2.5 m: POM raised, joined
  r <- run(tree(list(list(1, "s1", 1, 778, 1.3, "alive", NA, NA), list(3, "s9", 2, 670, 2.5, "alive", NA, NA))))
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$route, "impossible_recruit_gap")
  expect_equal(r$dt[, uniqueN(ReconstructedStemID)], 1L)
  # 595 mm at 3.2 m, then 271 mm at 1.6 m ten years later: HOM lowered, another stem
  r <- run(tree(list(list(1, "s1", 1, 595, 3.2, "alive", NA, NA), list(3, "s9", 2, 271, 1.6, "alive", NA, NA))))
  expect_equal(r$candidates$why, gap_why)
  # HOM raised but a larger DBH: growth, not a moved POM
  r <- run(tree(list(list(1, "s1", 1, 203, 1.3, "alive", NA, NA), list(3, "s9", 2, 261, 3.9, "alive", NA, NA))))
  expect_equal(r$candidates$why, gap_why)
})

test_that("across a census without a measurement: a trunk only missed in one census is not joined, nor after a death record", {
  # missed in census 2 ("missing", no DBH), ordinary growth, buttressed tree: not a moved POM
  r <- run(tree(list(list(1, "s1", 1, 341, 1.3, "alive", "B", NA), list(2, "s1", 1, NA, 1.3, "missing", NA, NA), list(3, "s9", 2, 357, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "across a census without a measurement: no recorded HOM rise with a smaller DBH")
  expect_equal(nrow(r$pairs), 0L)
  # recorded dead in census 2: the stem that follows is another stem
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 1, NA, NA, "dead", "D", NA), list(3, "s9", 2, 880, 1.3, "alive", NA, NA))))
  expect_equal(r$candidates$why, "death or break record in the census between")
  expect_equal(nrow(r$pairs), 0L)
  # no more than one census in between
  r <- run(tree(list(list(1, "s1", 1, 341, 1.3, "alive", NA, NA), list(4, "s9", 2, 357, 1.3, "alive", NA, NA))))
  expect_equal(nrow(r$candidates), 0L)
})

test_that("a drop without a sign of a moved POM stays split (recording error or shrinking trunk)", {
  # 435 mm then 288 mm, both at 1.3 m, no buttress code, ordinary growth form
  rows <- list(list(1, "s1", 1, 296, 1.3, "alive", NA, NA), list(2, "s1", 1, 435, 1.3, "alive", NA, NA), list(3, "s9", 2, 288, 1.3, "alive", NA, NA))
  r <- run(tree(rows))
  expect_equal(r$candidates$why, "no evidence of a moved point of measurement")
  expect_equal(r$candidates$evidence, "")
  expect_equal(nrow(r$pairs), 0L)
})

test_that("the signs of a moved POM are read on the whole tree, in any census", {
  base <- list(list(1, "s1", 1, 642, 1.3, "alive", NA, NA), list(2, "s9", 2, 420, 1.3, "alive", NA, NA))
  # a buttress code on the later stem, two censuses later
  r <- run(tree(c(base, list(list(4, "s9", 2, 440, 1.3, "alive", "B;cylN", NA)))))
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$evidence, "buttress")
  # another stem of the tree measured at 2 m
  r <- run(tree(c(base, list(list(1, "s5", 3, 150, 2, "alive", NA, NA), list(2, "s5", 3, 155, 2, "alive", NA, NA)))))
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$evidence, "HOM")
  # a code containing the letter B is not the buttress code
  r <- run(tree(c(base, list(list(4, "s9", 2, 440, 1.3, "alive", "BR", NA)))))
  expect_equal(r$candidates$why, "no evidence of a moved point of measurement")
})

test_that("a strangler growth form is a sign when irregular_forms is given", {
  x <- tree(list(list(1, "s1", 1, 750, 1.3, "alive", NA, NA), list(2, "s9", 2, 600, 1.3, "alive", NA, NA)), growth_form = "strangler")
  expect_equal(run(x)$candidates$why, "no evidence of a moved point of measurement")
  r <- run(x, irregular_forms = "strangler|estrangulador")
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$evidence, "growth form")
  # the label in the species list (Spanish) works too; the test is on the whole tree
  x[, growth_form := "Estrangulador"]
  expect_equal(nrow(run(x, irregular_forms = "strangler|estrangulador")$pairs), 1L)
})

test_that("a species that usually has buttresses is a sign only where no buttress codes exist", {
  # 595 mm then 370 mm at 1.3 m in censuses 1-2, no code anywhere; species "lue"
  x <- tree(list(list(1, "s1", 1, 595, 1.3, "alive", NA, NA), list(2, "s9", 2, 370, 1.3, "alive", NA, NA)))
  x[, Mnemonic := "lue"]
  expect_equal(run(x)$candidates$why, "no evidence of a moved point of measurement")
  r <- run(x, buttress_species = "lue", coded_censuses = 3:9)
  expect_equal(nrow(r$pairs), 1L)
  expect_equal(r$pairs$evidence, "species")
  # another species: still split
  expect_equal(run(x, buttress_species = "oth", coded_censuses = 3:9)$candidates$why, "no evidence of a moved point of measurement")
  # the trunk is measured in census 3 (codes exist) without a B code: its own record outranks the species
  y <- rbind(x, tree(list(list(3, "s9", 2, 380, 1.3, "alive", NA, NA)))[, Mnemonic := "lue"])
  expect_equal(run(y, buttress_species = "lue", coded_censuses = 3:9)$candidates$why, "no evidence of a moved point of measurement")
})

test_that("mr_buttress_species lists species by the share of large trees with the buttress code", {
  mk <- function(tag, sp, code, census = 3L) data.table(Tag = tag, CensusID = census, DBH = 400, Mnemonic = sp, ListOfTSM = code)
  d <- rbindlist(c(lapply(1:10, function(i) mk(paste0("a", i), "bt", "B;M")),   # 100%
                   lapply(1:10, function(i) mk(paste0("b", i), "nb", NA_character_)), # 0%
                   lapply(1:9, function(i) mk(paste0("c", i), "few", "B")),      # 9 trees: too few
                   list(mk("d1", "old", "B", census = 1L))))                     # no census with codes
  expect_equal(mr_buttress_species(d, coded_censuses = 3:9), "bt")
})

test_that("a larger DBH beyond the bounds is joined only when a recorded HOM change explains it", {
  rise_why <- "larger DBH not explained by a recorded HOM change"
  # 785 mm at 1.3 m, then 923 mm at 8.8 m: the recorded DBH grew 4.6 cm/yr,
  # the taper-corrected one (1160 mm) 12.5 cm/yr
  x <- tree(list(list(1, "s1", 1, 785, 1.3, "alive", NA, NA), list(2, "s9", 2, 923, 8.8, "alive", NA, NA)))
  x[CensusID == 2, dbh_with_best_candidate_taper_corrected := 1160]
  r <- run(x)
  expect_equal(nrow(r$pairs), 1L)
  # the same rise at an unchanged HOM (buttressed tree): not a moved POM
  y <- tree(list(list(1, "s1", 1, 785, 1.3, "alive", "B", NA), list(2, "s9", 2, 1000, 1.3, "alive", NA, NA)))
  expect_equal(run(y)$candidates$why, rise_why)
  # HOM changed, but the recorded DBH also grew beyond the bound (7.2 cm/yr)
  z <- tree(list(list(1, "s1", 1, 785, 1.3, "alive", NA, NA), list(2, "s9", 2, 1000, 3, "alive", NA, NA)))
  expect_equal(run(z)$candidates$why, rise_why)
})

test_that("idempotent", {
  x <- tree(list(list(1, "s1", 1, 1775, 1.3, "alive", NA, NA), list(2, "s1", 2, 1140, 1.3, "alive", "B", NA)))
  r1 <- run(x)
  expect_equal(nrow(r1$pairs), 1L)
  expect_identical(r1$obs_pairs[, .(obs0, obs1)], data.table(obs0 = 1L, obs1 = 2L)) # the posterior gets the export joins
  r2 <- run(r1$dt)
  expect_equal(nrow(r2$pairs), 0L)
  expect_identical(r2$dt, r1$dt)
})

test_that("posterior samples: the pair is joined only where it is split with a clean end/start and one pin", {
  pairs <- data.table(Tag = "T1", obs0 = 1L, obs1 = 2L, c0 = 1L, c1 = 2L)
  info <- data.table(tag = "T1", obs = 1:3, c = c(1L, 2L, 2L), pin = NA_character_)
  paths <- data.table(
    tag = "T1", treeID = "1", path_sig = c("a", "b", "c"), path_count = c(100L, 60L, 40L), path_prob = c(.5, .3, .2),
    recon = c("1:1;2:2;3:3", "1:1;2:1;3:3", "1:1;2:2;3:1") # split cleanly | already joined | obs 1 continues into obs 3
  )
  np <- apply_measurement_rejoin_to_paths(paths, pairs, info, verbose = FALSE)
  expect_equal(np[, sum(path_count)], 200L)
  expect_equal(np[recon == "1:1;2:1;3:3", path_count], 160L) # the clean split collapses into the joined path
  expect_equal(np[recon == "1:1;2:2;3:1", path_count], 40L) # unchanged
  # a pair across a census without a measurement (c1 = c0 + 2) is joined the same way
  pg <- copy(pairs)[, c1 := 3L]
  ig <- copy(info)[, c := c(1L, 3L, 3L)]
  ng <- apply_measurement_rejoin_to_paths(paths, pg, ig, verbose = FALSE)
  expect_equal(ng[recon == "1:1;2:1;3:3", path_count], 160L)
  info2 <- copy(info)[, pin := c("p1", "p2", NA)]
  np2 <- apply_measurement_rejoin_to_paths(paths, pairs, info2, verbose = FALSE)
  expect_equal(np2[recon == "1:1;2:2;3:3", path_count], 100L) # two pins: not joined
})
