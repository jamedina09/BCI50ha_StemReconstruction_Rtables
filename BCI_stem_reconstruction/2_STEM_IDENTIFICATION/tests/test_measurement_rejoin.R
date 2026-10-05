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
tree <- function(rows, tag = "T1") {
  x <- rbindlist(lapply(rows, function(r) {
    data.table(
      Tag = tag, CensusID = as.integer(r[[1]]), StemID = as.character(r[[2]]), ReconstructedStemID = as.numeric(r[[3]]),
      DBH = as.numeric(r[[4]]), HOM = as.numeric(r[[5]]), Status = r[[6]], ListOfTSM = r[[7]], TrueStemID = r[[8]]
    )
  }))
  x[, `:=`(
    single_stem_tags = FALSE, dbh_with_best_candidate_taper_corrected = DBH,
    ExactDate = as.Date("1982-01-01") + round((CensusID - 1L) * 365.25 * 3), obs_row_id = seq_len(.N),
    ReconstructionMethod = "dp"
  )]
  x
}
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
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "broken below", "R", NA), list(2, "s1", 2, 600, 1.3, "alive", NA, NA))))
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
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, "p1"), list(2, "s1", 2, 600, 1.3, "alive", NA, "p2"))))
  expect_equal(r$candidates$why, "two different pins")
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 2, 600, 4, "alive", NA, "p2"), list(3, "s1", 2, 610, 4, "alive", NA, "p2"))))
  expect_equal(r$dt[, unique(ReconstructedStemID)], 2) # the pinned stem keeps its id
})

test_that("no join when the ending stem has a later record or the new stem an earlier one", {
  r <- run(tree(list(list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(2, "s1", 1, NA, 1.3, "dead", NA, NA), list(2, "s1", 2, 600, 1.3, "alive", NA, NA))))
  expect_equal(nrow(r$pairs), 0L)
})

test_that("a unique candidate is needed: two stems ending at the same census stay split", {
  r <- run(tree(list(
    list(1, "s1", 1, 900, 1.3, "alive", NA, NA), list(1, "s2", 3, 880, 1.3, "alive", NA, NA),
    list(2, "s9", 2, 700, 1.3, "alive", NA, NA)
  )))
  expect_equal(nrow(r$pairs), 0L)
})

test_that("idempotent", {
  x <- tree(list(list(1, "s1", 1, 1775, 1.3, "alive", NA, NA), list(2, "s1", 2, 1140, 1.3, "alive", NA, NA)))
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
  info2 <- copy(info)[, pin := c("p1", "p2", NA)]
  np2 <- apply_measurement_rejoin_to_paths(paths, pairs, info2, verbose = FALSE)
  expect_equal(np2[recon == "1:1;2:2;3:3", path_count], 100L) # two pins: not joined
})
