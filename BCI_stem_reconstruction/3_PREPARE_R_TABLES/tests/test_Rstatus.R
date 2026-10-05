# ========================================================================
# TESTS: Rstatus and dbh rules of 2_create_R_tables_BCI.R
# ========================================================================
# Run from the project root (or from this folder):
#   testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R")
#
#   1. Hand-written fixtures, one per rule and correction case
#   2. Exhaustive raw inputs: rstatus_functions.R (what the script runs) gives
#      exactly the Rstatus and dbh of the reference implementation
#      (reference_rstatus.R), and every invariant holds
#   3. Rstatus outputs: validity of every P/A/G/D sequence, soundness and
#      reachability
#   4. Real data: the exported R tables equal the rules and satisfy every
#      invariant (Rstatus, dbh, DFstatus, ExactDate, date, location, stemID);
#      offending stem / tree IDs are printed
#
# Options (environment variables):
#   RSTATUS_TEST_LEVEL   "full" (default; about 10 min on 16 cores) or "quick" (a few min)
#   RSTATUS_RTABLES_DIR  folder with bci.stem1..9.Rdata (default DATA/RTABLES)
#   RSTATUS_STAGE2_FILE  stage-2 table the R tables were built from (default
#                        DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds)
#   RSTATUS_CORES        cores for the reference implementation (default: all but 2)
# ========================================================================
suppressMessages({
  library(testthat)
  library(data.table)
})

find_test_file <- function(f) {
  cand <- c(f, file.path("tests", f), file.path("BCI_stem_reconstruction", "3_PREPARE_R_TABLES", "tests", f))
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("cannot find ", f, ": run from the project root or from the tests folder")
  normalizePath(hit)
}
tests_dir <- dirname(find_test_file("reference_rstatus.R"))
stage3_dir <- dirname(tests_dir)
project_root <- normalizePath(file.path(stage3_dir, "..", ".."))
source(file.path(stage3_dir, "rstatus_functions.R"))
source(file.path(tests_dir, "reference_rstatus.R"))
source(file.path(tests_dir, "rstatus_validity.R"))

level <- Sys.getenv("RSTATUS_TEST_LEVEL", "full")
cores <- as.integer(Sys.getenv("RSTATUS_CORES", max(1L, parallel::detectCores() - 2L)))
cat(sprintf("\nRstatus tests: level %s, %d cores\n", level, cores))

# Rstatus and dbh of rstatus_functions.R for a table of cells
# (data.table(TreeID, StemID, census, Status, DBH), one row per stem x census).
rules_table <- function(cells) {
  cells <- cells[order(TreeID, StemID, census)]
  n <- max(cells$census)
  stems <- unique(cells[, .(TreeID, StemID)])
  St <- matrix(cells$Status, ncol = n, byrow = TRUE)
  Db <- matrix(cells$DBH, ncol = n, byrow = TRUE)
  o <- compute_rstatus(St, Db, stems$TreeID)
  data.table(
    TreeID = rep(stems$TreeID, n), StemID = rep(stems$StemID, n), census = rep(seq_len(n), each = nrow(stems)),
    Rstatus = as.vector(o$Rstatus), dbh = as.vector(o$dbh)
  )
}
same_num <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
same_chr <- same_num
show_ids <- function(ids, k = 10) paste(head(unique(ids), k), collapse = ", ")

# ========================================================================
# 1. HAND-WRITTEN FIXTURES
# ========================================================================
# One tree from raw records: stems = list of raw Status vectors (NA = no row
# or no status); dbh = list of DBH vectors (default 10 on every "alive").
tree <- function(stems, dbh = NULL) {
  st <- do.call(rbind, stems)
  db <- if (is.null(dbh)) matrix(ifelse(st %in% "alive", 10, NA_real_), nrow(st)) else do.call(rbind, lapply(dbh, as.numeric))
  out <- compute_rstatus(st, db, rep("T", nrow(st)))
  list(R = apply(out$Rstatus, 1, paste, collapse = ""), dbh = unname(out$dbh), ref = ref_tree(st, db))
}
expect_tree <- function(x, R, dbh = NULL) {
  expect_identical(x$R, R)
  expect_identical(apply(x$ref$Rstatus, 1, paste, collapse = ""), R) # the reference agrees
  expect_equal(x$dbh, x$ref$dbh)
  if (!is.null(dbh)) expect_equal(x$dbh, do.call(rbind, lapply(dbh, as.numeric)))
}

test_that("missed census between alive records is A, with no DBH (never imputed)", {
  expect_tree(tree(list(c("alive", NA, "alive"))), "AAA", list(c(10, NA, 10)))
  expect_tree(tree(list(c("alive", "missing", "alive"))), "AAA", list(c(10, NA, 10)))
  expect_tree(tree(list(c("alive", NA, NA, "alive"))), "AAAA", list(c(10, NA, NA, 10)))
})

test_that("false death: a dead record followed by alive is A, and its DBH is kept", {
  expect_tree(tree(list(c("alive", "dead", "alive"))), "AAA")
  expect_tree(tree(list(c("alive", "stem dead", "alive")), list(c(10, 12, 14))), "AAA", list(c(10, 12, 14)))
  expect_tree(tree(list(c("alive", "dead", NA, "dead", "alive"))), "AAAAA")
  expect_tree(tree(list(c("alive", "broken below", "alive"))), "AAA")
})

test_that("dead record with a DBH, never alive again: G/D from that census, DBH kept as recorded", {
  expect_tree(tree(list(c("alive", "dead", NA)), list(c(10, 12, NA))), "ADD", list(c(10, 12, NA)))
  expect_tree(tree(list(c("alive", "stem dead", "dead")), list(c(10, 12, NA))), "ADD", list(c(10, 12, NA)))
  expect_tree(
    tree(list(c("alive", "stem dead", NA), c("alive", "alive", "alive")), list(c(10, 12, NA), c(10, 11, 12))),
    c("AGG", "AAA"), list(c(10, 12, NA), c(10, 11, 12))
  )
})

test_that("request example: G, unrecorded, G beside A, unrecorded, A -> G G G and A A A, never D", {
  # stem 1 dies at census 2; census 3 has no record of either stem
  expect_tree(tree(list(c("alive", "dead", NA, "dead"), c("alive", "alive", NA, "alive"))), c("AGGG", "AAAA"))
  # the same with stem 1 first recorded dead (never alive)
  expect_tree(tree(list(c(NA, "dead", NA, "dead"), c("alive", "alive", NA, "alive"))), c("PGGG", "AAAA"))
})

test_that("unrecorded census between two dead records stays dead", {
  expect_tree(tree(list(c("alive", "dead", NA, "dead"))), "ADDD")
  expect_tree(tree(list(c("alive", "dead", NA, "dead"), c("alive", "alive", "alive", "alive"))), c("AGGG", "AAAA"))
})

test_that("P only before the first record (late recruits)", {
  expect_tree(tree(list(c(NA, NA, "alive"))), "PPA")
  expect_tree(tree(list(c("missing", "alive", "alive"))), "PAA")
  expect_tree(tree(list(c(NA, "alive", "alive"), c("alive", "alive", "alive"))), c("PAA", "AAA"))
})

test_that("a dead record before the stem's first alive record is A (false death)", {
  expect_tree(tree(list(c(NA, "broken below", "alive")), list(c(NA, NA, 11))), "PAA", list(c(NA, NA, 11)))
  expect_tree(tree(list(c("dead", "alive", "alive"))), "AAA")
  expect_tree(tree(list(c("dead", NA, "alive"))), "AAA")
  expect_tree(tree(list(c(NA, "dead", "alive")), list(c(NA, 9, 11))), "PAA", list(c(NA, 9, 11)))
})

test_that("a stem never alive is P until its first record, then G or D", {
  expect_tree(tree(list(c(NA, "broken below", "dead"))), "PDD")
  expect_tree(tree(list(c(NA, "broken below", "dead"), c("alive", "alive", "alive"))), c("PGG", "AAA"))
  expect_tree(tree(list(c("dead", NA, NA))), "DDD")
  expect_tree(tree(list(c(NA, "dead", NA), c("alive", "alive", "dead"))), c("PGD", "AAD"))
})

test_that("a stem never recorded stays P in every census", {
  expect_tree(tree(list(c(NA, "missing", NA))), "PPP")
  expect_tree(tree(list(c(NA, NA, NA), c("alive", "alive", "dead"))), c("PPP", "AAD"))
})

test_that("alive, then no record to the end: dead from the first census without a record", {
  expect_tree(tree(list(c("alive", NA, NA))), "ADD")
  expect_tree(tree(list(c("alive", "alive", "missing"))), "AAD") # missing in the final census
  expect_tree(tree(list(c("alive", NA, NA), c("alive", "alive", "alive"))), c("AGG", "AAA"))
})

test_that("unrecorded between alive and a later dead record: dead from the first unrecorded census", {
  expect_tree(tree(list(c("alive", NA, "dead"))), "ADD")
  expect_tree(tree(list(c("alive", NA, "dead"), c("alive", "alive", "alive"))), c("AGG", "AAA"))
})

test_that("tree status: G while the tree lives, D once every stem is dead, and D stays D", {
  expect_tree(tree(list(c("alive", "dead", "dead", "dead"), c("alive", "alive", "dead", NA))), c("AGDD", "AADD"))
  expect_tree(tree(list(c("alive", "dead", NA), c("alive", "alive", "alive"), c(NA, NA, "alive"))), c("AGG", "AAA", "PPA"))
})

test_that("single-stem trees never get G", {
  expect_tree(tree(list(c("alive", "dead", "dead"))), "ADD")
  expect_tree(tree(list(c("alive", "stem dead", NA))), "ADD")
})

test_that("a tree with no record in a census but alive later is alive there (no D)", {
  expect_tree(tree(list(c("alive", NA, "alive"), c("alive", NA, "dead"))), c("AAA", "AGG"))
  expect_tree(tree(list(c("alive", NA, NA, "alive"))), "AAAA")
})

test_that("broken below with a DBH is a measured resprout (A); without a DBH it is a dead record", {
  expect_tree(tree(list(c("alive", "broken below")), list(c(10, 5))), "AA", list(c(10, 5)))
  expect_tree(tree(list(c("alive", "broken below")), list(c(10, NA))), "AD", list(c(10, NA)))
})

test_that("a missing or status-less record with a DBH is alive, and the DBH is kept", {
  expect_tree(tree(list(c("alive", "missing", "alive")), list(c(10, 11, 12))), "AAA", list(c(10, 11, 12)))
  expect_tree(tree(list(c("alive", NA, "alive")), list(c(10, 11, 12))), "AAA", list(c(10, 11, 12)))
})

test_that("an alive record without a DBH is A with no DBH", {
  expect_tree(tree(list(c("alive", "alive", "alive")), list(c(10, NA, 12))), "AAA", list(c(10, NA, 12)))
})

test_that("an unknown raw status stops", {
  expect_error(rstatus_raw_code("sleeping", NA_real_), "unknown raw status")
  expect_error(ref_evidence("sleeping", NA_real_), "unknown raw status")
})

# ========================================================================
# 2. EXHAUSTIVE RAW INPUTS
# ========================================================================
# Raw states found in the BCI data (Status x DBH present) plus "no row". The
# 12 states fall into 6 behaviour classes (equivalence checked below), so
# larger trees use one representative per class; alive / dead / no record is
# complete for Rstatus.
raw_states <- data.table(
  name = c("alive+dbh", "alive", "bb+dbh", "bb", "dead+dbh", "dead", "stemdead+dbh", "stemdead", "missing+dbh", "missing", "nostatus+dbh", "norow"),
  Status = c("alive", "alive", "broken below", "broken below", "dead", "dead", "stem dead", "stem dead", "missing", "missing", NA, NA),
  dbh = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
  class = c("alive+dbh", "alive", "alive+dbh", "dead", "dead+dbh", "dead", "dead+dbh", "dead", "missing+dbh", "norow", "alive+dbh", "norow")
)
states <- function(nm) raw_states[match(nm, name)]
S6 <- states(c("alive+dbh", "alive", "dead", "dead+dbh", "missing+dbh", "norow"))
S5 <- states(c("alive+dbh", "alive", "dead", "dead+dbh", "norow"))
S4 <- states(c("alive+dbh", "dead", "dead+dbh", "norow"))
S3 <- states(c("alive+dbh", "dead", "norow"))

# Every combination of raw states for trees of S stems x n censuses.
gen_cells <- function(S, n, st) {
  K <- nrow(st)
  k <- 0:(K^(S * n) - 1)
  rbindlist(lapply(seq_len(S), function(s) {
    rbindlist(lapply(seq_len(n), function(c) {
      d <- (k %/% (K^((s - 1L) * n + (c - 1L)))) %% K
      data.table(
        TreeID = as.character(k), StemID = paste0(k, "_", s), census = c, state = st$name[d + 1L],
        Status = st$Status[d + 1L], DBH = ifelse(st$dbh[d + 1L], 100 * s + c, NA_real_)
      )
    }))
  }))
}

check_exhaustive <- function(S, n, st) {
  t0 <- Sys.time()
  cells <- gen_cells(S, n, st)
  x <- rules_table(cells)
  r <- ref_rstatus_table(cells[, .(TreeID, StemID, census, Status, DBH)], cores)
  x <- r[x, on = .(StemID, census)]
  x <- cells[, .(StemID, census, raw_DBH = DBH, code = rstatus_raw_code(Status, DBH))][x, on = .(StemID, census)]
  n_trees <- uniqueN(x$TreeID)
  bad_R <- x[Rstatus != Rstatus_ref, unique(TreeID)]
  bad_D <- x[!same_num(dbh, dbh_ref), unique(TreeID)]
  expect(length(bad_R) == 0L, sprintf("Rstatus differs from the reference in %d trees, e.g. %s", length(bad_R), show_ids(bad_R)))
  expect(length(bad_D) == 0L, sprintf("dbh differs from the reference in %d trees, e.g. %s", length(bad_D), show_ids(bad_D)))
  # DBH invariants: dbh is the raw DBH everywhere; a P cell never has one
  expect_identical(x[Rstatus == "P" & !is.na(dbh), .N], 0L)
  expect_identical(x[!same_num(dbh, raw_DBH), .N], 0L)
  # Rstatus validity: invalid only when a stem has no record at all (all P)
  w <- dcast(x, TreeID + StemID ~ census, value.var = "Rstatus")
  seqs <- do.call(paste0, w[, -(1:2)])
  v <- classify_trees(split(seqs, w$TreeID))
  no_rec_tree <- x[, .(no_rec = all(code == "N")), by = .(TreeID, StemID)][, .(any_no_rec = any(no_rec)), by = TreeID]
  setkey(no_rec_tree, TreeID)
  any_no_rec <- no_rec_tree[names(split(seqs, w$TreeID)), any_no_rec]
  expect_identical(sum(!v$valid & !any_no_rec), 0L)
  expect_true(all(v$reasons[!v$valid] == "all P (stem never recorded)"))
  cat(sprintf(
    "  exhaustive %d stem(s) x %d censuses, %d raw states: %s trees, %d invalid (all of them hold a never-recorded stem), %.0f s\n",
    S, n, nrow(st), format(n_trees, big.mark = ","), sum(!v$valid), as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
  invisible(x)
}

test_that("all 12 raw states, one stem, 3-5 censuses: rules == reference, invariants hold", {
  for (n in 3:5) check_exhaustive(1L, n, raw_states)
})

test_that("the 12 raw states fall into 6 behaviour classes", {
  cells <- gen_cells(1L, 5L, raw_states)
  # every state replaced by its class representative (DBH presence is kept,
  # so the DBH values stay the same: 100 * stem + census)
  cl <- copy(cells)[, state := raw_states$class[match(state, raw_states$name)]]
  cl[, Status := raw_states$Status[match(state, raw_states$name)]]
  expect_identical(!is.na(cl$DBH), raw_states$dbh[match(cl$state, raw_states$name)])
  a <- rules_table(cells)
  b <- rules_table(cl)
  expect_identical(a$Rstatus, b$Rstatus)
  expect_true(all(same_num(a$dbh, b$dbh)))
})

test_that("class representatives, larger trees: rules == reference, invariants hold", {
  check_exhaustive(1L, 6L, S6)
  check_exhaustive(2L, 3L, S6)
  check_exhaustive(3L, 3L, S3)
  if (level == "full") {
    check_exhaustive(2L, 4L, S6)
    check_exhaustive(2L, 5L, S4)
    check_exhaustive(3L, 3L, S5)
    check_exhaustive(3L, 4L, S3)
  }
})

# ========================================================================
# 3. RSTATUS OUTPUTS: VALIDITY, SOUNDNESS, REACHABILITY
# ========================================================================
all_seqs <- function(n) do.call(paste0, rev(expand.grid(rep(list(c("P", "A", "G", "D")), n), stringsAsFactors = FALSE)))
grammar_ok <- function(s) grepl("^P*A*G*D*$", s) & grepl("[AGD]", s)

test_that("every P/A/G/D sequence of 3-6 censuses: valid stems are exactly P*A*G*D* with a non-P", {
  for (n in 3:6) {
    s <- all_seqs(n)
    expect_length(s, 4^n)
    expect_identical(grammar_ok(s), !nzchar(stem_reasons(s)))
    expect_identical(sum(grammar_ok(s)), as.integer(choose(n + 3, 3) - 1))
    cat(sprintf("  n = %d: %d sequences, %d valid stems, %d valid single-stem trees\n", n, length(s), sum(grammar_ok(s)), sum(classify_trees(as.list(s))$valid)))
  }
  # examples from the request
  expect_true(all(grammar_ok(c("AAA", "PAAG", "PAGGD", "ADD"))))
  expect_false(any(grammar_ok(c("AGA", "AP", "DG", "DA", "GA", "PPP"))))
  # PAAG and PAGGD need another stem alive later in the same tree
  expect_false(any(classify_trees(list("PAAG"))$valid, classify_trees(list("PAGGD"))$valid))
  expect_true(classify_trees(list(c("PAAG", "AAAA")))$valid)
})

valid_combos <- function(S, n) {
  ok <- all_seqs(n)
  ok <- ok[grammar_ok(ok)]
  idx <- as.matrix(expand.grid(rep(list(seq_along(ok)), S)))
  idx <- idx[apply(idx, 1, function(z) !is.unsorted(z)), , drop = FALSE] # stems are exchangeable
  trees <- lapply(seq_len(nrow(idx)), function(i) ok[idx[i, ]])
  v <- classify_trees(trees)
  data.table(combo = vapply(trees, function(t) paste(sort(t), collapse = "|"), ""), valid = v$valid, flag = v$flag, reasons = v$reasons, request_invalid = v$request_checks_invalid)
}

test_that("trees of 2 and 3 stems: tree-level validity, and the gap in the two request checks", {
  sizes <- list(c(2, 3), c(2, 4), c(2, 5), c(2, 6), c(3, 3), c(3, 4), c(3, 5))
  if (level != "full") sizes <- sizes[1:5]
  for (z in sizes) {
    V <- valid_combos(z[1], z[2])
    missed <- V[!valid & !request_invalid]
    # the request checks miss only "G with no stem alive then or later"
    expect_true(all(missed$reasons == "G while no stem is A then or later (tree dead)"))
    cat(sprintf(
      "  %d stems x %d censuses: %d combinations of valid stems, %d valid (%d flagged), %d invalid (%d pass the two request checks)\n",
      z[1], z[2], nrow(V), V[valid == TRUE, .N], V[valid & flag, .N], V[valid == FALSE, .N], nrow(missed)
    ))
  }
})

test_that("soundness and reachability: every produced combination is valid, every valid one is produced", {
  sizes <- list(c(1, 3), c(1, 4), c(1, 5), c(1, 6), c(2, 3), c(2, 4), c(2, 5), c(2, 6), c(3, 3), c(3, 4))
  if (level != "full") sizes <- sizes[c(1:6, 9)]
  for (z in sizes) {
    S <- z[1]
    n <- z[2]
    cells <- gen_cells(S, n, S3)
    x <- rules_table(cells)
    w <- dcast(x, TreeID + StemID ~ census, value.var = "Rstatus")
    w[, s := do.call(paste0, .SD), .SDcols = -(1:2)]
    produced <- unique(w[, .(combo = paste(sort(s), collapse = "|")), by = TreeID]$combo)
    pv <- classify_trees(strsplit(produced, "|", fixed = TRUE))
    # soundness: the only invalid output is a stem never recorded (all P)
    expect_true(all(pv$reasons[!pv$valid] == "all P (stem never recorded)"))
    V <- if (S == 1) {
      ok <- all_seqs(n)
      ok <- ok[grammar_ok(ok)]
      v1 <- classify_trees(as.list(ok))
      data.table(combo = ok, valid = v1$valid, flag = v1$flag)
    } else {
      valid_combos(S, n)
    }
    unreached <- V[valid == TRUE & !combo %in% produced, combo]
    expect(length(unreached) == 0L, sprintf("%d valid combinations never produced, e.g. %s", length(unreached), show_ids(unreached, 5)))
    cat(sprintf(
      "  %d stem(s) x %d censuses: %d produced (%d invalid, all with a never-recorded stem); reached %d of %d valid (%d flagged)\n",
      S, n, length(produced), sum(!pv$valid), V[valid == TRUE & combo %in% produced, .N], V[valid == TRUE, .N], V[valid & flag, .N]
    ))
  }
})

# ========================================================================
# 4. REAL DATA: THE EXPORTED R TABLES
# ========================================================================
rtables_dir <- Sys.getenv("RSTATUS_RTABLES_DIR", file.path(project_root, "BCI_stem_reconstruction", "DATA", "RTABLES"))
stage2_file <- Sys.getenv("RSTATUS_STAGE2_FILE", file.path(project_root, "BCI_stem_reconstruction", "DATA", "PROCESSED", "complete_dataset_final_with_reconstructed_stemids.rds"))
real_ok <- file.exists(stage2_file) && length(list.files(rtables_dir, "^bci\\.stem[0-9]+\\.Rdata$")) > 0L

if (real_ok) {
  cat(sprintf("\nReal data: %s\n           %s\n", rtables_dir, stage2_file))
  # Raw input, built the way the script builds it (one row per stem x census)
  raw <- as.data.table(readRDS(stage2_file))
  raw <- raw[, .(
    TreeID = as.character(TreeID),
    StemID = paste(as.character(TreeID), as.character(as.numeric(as.character(ReconstructedStemID))), sep = "_"),
    census = as.integer(as.character(CensusID)), Status = as.character(Status), DBH = as.numeric(DBH),
    ExactDate = as.Date(ExactDate), PX = as.numeric(PX), PY = as.numeric(PY)
  )]
  n_cens <- max(raw$census)
  grid <- unique(raw[, .(TreeID, StemID)])[, .(census = seq_len(n_cens)), by = .(TreeID, StemID)]
  cells <- raw[grid, on = .(TreeID, StemID, census)]
  files <- file.path(rtables_dir, sprintf("bci.stem%d.Rdata", seq_len(n_cens)))
  # The exported stemID is the stem number within its tree; StemID below is
  # the script's internal key "TreeID_stemID" ("TreeID_NA" for stemID NA).
  rt <- rbindlist(lapply(seq_len(n_cens), function(i) {
    e <- new.env()
    load(files[i], envir = e)
    as.data.table(get(ls(e)[1], envir = e))[, .(
      StemID = paste(treeID, stemID, sep = "_"), stemID, TreeID = treeID, census = i, Rstatus, dbh, DFstatus,
      ExactDate, date, gx, gy, quadrat, order = .I
    )]
  }))
}

test_that("real data: one row per stem per census, same stems in the same order in every table", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  expect_identical(rt[, .N, by = .(StemID, census)][N > 1L, .N], 0L)
  expect_identical(nrow(rt), nrow(cells))
  ord <- dcast(rt, order ~ census, value.var = "StemID")
  expect_true(all(vapply(ord[, -1], identical, logical(1), ord[[2]])))
  expect_setequal(unique(rt$StemID), unique(cells$StemID))
})

test_that("real data: exported Rstatus and dbh equal the rules and the reference, cell by cell", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  x <- rules_table(cells[, .(TreeID, StemID, census, Status, DBH)])
  r <- ref_rstatus_table(cells[, .(TreeID, StemID, census, Status, DBH)], cores)
  y <- rt[, .(StemID, census, Rstatus_rt = Rstatus, dbh_rt = dbh)][x, on = .(StemID, census)][r, on = .(StemID, census)]
  b1 <- y[Rstatus_rt != Rstatus, unique(StemID)]
  b2 <- y[Rstatus != Rstatus_ref, unique(StemID)]
  b3 <- y[!same_num(dbh_rt, dbh) | !same_num(dbh, dbh_ref), unique(StemID)]
  expect(length(b1) == 0L, sprintf("exported Rstatus differs from the rules for %d stems, e.g. %s", length(b1), show_ids(b1)))
  expect(length(b2) == 0L, sprintf("rules differ from the reference for %d stems, e.g. %s", length(b2), show_ids(b2)))
  expect(length(b3) == 0L, sprintf("exported dbh differs for %d stems, e.g. %s", length(b3), show_ids(b3)))
})

test_that("real data: Rstatus invariants (stem sequences and tree consistency)", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  w <- dcast(rt, TreeID + StemID ~ census, value.var = "Rstatus")
  seqs <- do.call(paste0, w[, -(1:2)])
  by_tree <- split(seqs, w$TreeID)
  v <- classify_trees(by_tree)
  # the only invalid stems are stems with no record at all (all P)
  codes <- cells[, .(no_rec = all(rstatus_raw_code(Status, DBH) == "N")), by = .(TreeID, StemID)]
  all_P <- w$StemID[grepl("^P+$", seqs)]
  expect_setequal(all_P, codes[no_rec == TRUE, StemID])
  bad_trees <- names(by_tree)[!v$valid]
  other <- bad_trees[v$reasons[!v$valid] != "all P (stem never recorded)"]
  expect(length(other) == 0L, sprintf("%d trees break the Rstatus rules, e.g. %s", length(other), show_ids(other)))
  cat(sprintf(
    "  real trees: %s | invalid only because of a never-recorded stem: %d | flagged: %d (stems dead without ever being A: %d)\n",
    format(length(by_tree), big.mark = ","), length(bad_trees), sum(v$flag), sum(grepl("[GD]", seqs) & !grepl("A", seqs))
  ))
  # a tree with no record in a census but alive later is alive there: none of its stems is D
  rec <- cells[, .(any_rec = any(rstatus_raw_code(Status, DBH) != "N")), by = .(TreeID, census)]
  alive_later <- rt[, .(alive_later = any(Rstatus == "A")), by = .(TreeID, census)][order(TreeID, -census)][, alive_later := cummax(alive_later) > 0, by = TreeID]
  gap <- rec[any_rec == FALSE][alive_later[alive_later == TRUE], on = .(TreeID, census), nomatch = 0L]
  d_in_gap <- rt[gap, on = .(TreeID, census), nomatch = 0L][Rstatus == "D", unique(TreeID)]
  expect(length(d_in_gap) == 0L, sprintf("%d trees have a D stem in an unrecorded census although alive later, e.g. %s", length(d_in_gap), show_ids(d_in_gap)))
})

test_that("real data: dbh is exported exactly as recorded", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  y <- rt[, .(StemID, census, Rstatus, dbh)][cells[, .(StemID, census, Status, raw_DBH = DBH)], on = .(StemID, census)]
  b1 <- y[Rstatus == "P" & !is.na(dbh), unique(StemID)]
  b2 <- y[!same_num(dbh, raw_DBH), unique(StemID)]
  on_dead <- y[Rstatus %in% c("G", "D") & !is.na(dbh)]
  b3 <- on_dead[!Status %in% c("dead", "stem dead"), unique(StemID)]
  expect(length(b1) == 0L, sprintf("dbh on a P cell for %d stems, e.g. %s", length(b1), show_ids(b1)))
  expect(length(b2) == 0L, sprintf("exported dbh differs from the raw DBH for %d stems, e.g. %s", length(b2), show_ids(b2)))
  expect(length(b3) == 0L, sprintf("dbh on a G/D cell that is not a dead / stem dead record for %d stems, e.g. %s", length(b3), show_ids(b3)))
  dead_dbh <- y[!is.na(raw_DBH) & Status %in% c("dead", "stem dead")]
  cat(sprintf(
    "  DBH on dead / stem dead records: %d, all kept | on G/D (never alive later): %d | on A (false death): %d\n",
    nrow(dead_dbh), dead_dbh[Rstatus != "A", .N], dead_dbh[Rstatus == "A", .N]
  ))
  chk <- file.path(dirname(rtables_dir), "CHECKS", "dbh_on_dead_records.csv")
  if (file.exists(chk)) {
    listed <- fread(chk)
    expect_setequal(paste(listed$StemID, listed$census), paste(on_dead$StemID, on_dead$census))
  }
})

test_that("real data: stemID is the reconstructed stem number within its tree", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  expect_type(rt$stemID, "integer")
  expect_true(all(is.na(rt$stemID) | rt$stemID >= 1L))
  # treeID + stemID identify one row per census, and the same stems in every census
  expect_identical(rt[, .N, by = .(TreeID, stemID, census)][N > 1L, .N], 0L)
  # stemID NA only for stems without a DP identity (no status and no DBH ever)
  no_id <- rt[is.na(stemID), unique(StemID)]
  expect_true(all(grepl("_NA$", no_id)))
  expect_identical(rt[is.na(stemID) & (!is.na(dbh) | !is.na(DFstatus)), .N], 0L)
  cat(sprintf(
    "  stemID: %s stems in %s trees | largest stem number %d | stems without a DP identity (stemID NA): %d\n",
    format(uniqueN(rt$StemID), big.mark = ","), format(uniqueN(rt$TreeID), big.mark = ","),
    max(rt$stemID, na.rm = TRUE), length(no_id)
  ))
})

test_that("real data: ExactDate is exactly the recorded date; date fills every row and equals it where recorded", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  y <- rt[, .(StemID, census, ExactDate, date)][cells[, .(StemID, census, raw_date = ExactDate)], on = .(StemID, census)]
  rec_days <- as.numeric(y$raw_date - as.Date("1960-01-01"))
  b1 <- y[!same_num(as.numeric(ExactDate), as.numeric(raw_date)), unique(StemID)]
  b2 <- y[is.na(date), unique(StemID)]
  b3 <- y[!is.na(raw_date) & date != rec_days, unique(StemID)]
  setorder(y, StemID, census)
  b4 <- y[, .(ok = all(diff(date) > 0)), by = StemID][ok == FALSE, StemID]
  expect(length(b1) == 0L, sprintf("ExactDate differs from the recorded date for %d stems, e.g. %s", length(b1), show_ids(b1)))
  expect(length(b2) == 0L, sprintf("a row has no date for %d stems, e.g. %s", length(b2), show_ids(b2)))
  expect(length(b3) == 0L, sprintf("date differs from the recorded date (days since 1960-01-01) for %d stems, e.g. %s", length(b3), show_ids(b3)))
  expect(length(b4) == 0L, sprintf("date does not increase from census to census for %d stems, e.g. %s", length(b4), show_ids(b4)))
  cat(sprintf(
    "  ExactDate: recorded %s, NA %s (exported as recorded) | date: every row (%s filled)\n",
    format(y[!is.na(raw_date), .N], big.mark = ","), format(y[is.na(raw_date), .N], big.mark = ","),
    format(y[is.na(raw_date), .N], big.mark = ",")
  ))
})

test_that("real data: one location per tree, the position most censuses agree on", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  loc <- unique(rt[, .(TreeID, gx, gy, quadrat)])
  multi <- loc[, .N, by = TreeID][N > 1L, TreeID]
  expect(length(multi) == 0L, sprintf("%d trees have more than one (gx, gy, quadrat), e.g. %s", length(multi), show_ids(multi)))
  # Recompute the votes from the raw coordinates: one vote per tree and census
  # (its modal pair; ties -> smallest x, then y). The exported pair must have
  # the most census votes of its tree.
  pool <- cells[!is.na(PX) & !is.na(PY), .N, by = .(TreeID, census, PX, PY)]
  setorder(pool, TreeID, census, -N, PX, PY)
  votes <- pool[pool[, .I[1L], by = .(TreeID, census)]$V1][, .(n_votes = .N), by = .(TreeID, PX, PY)]
  best <- votes[, .(max_votes = max(n_votes)), by = TreeID]
  chosen <- votes[unique(rt[!is.na(gx), .(TreeID, PX = gx, PY = gy)]), on = .(TreeID, PX, PY)]
  chosen <- best[chosen, on = "TreeID"]
  b1 <- chosen[is.na(n_votes) | n_votes < max_votes, TreeID]
  expect(length(b1) == 0L, sprintf("%d trees are not at their most-voted position, e.g. %s", length(b1), show_ids(b1)))
  no_xy <- setdiff(unique(rt$TreeID), chosen$TreeID)
  expect_setequal(no_xy, setdiff(unique(cells$TreeID), unique(pool$TreeID)))
  cat(sprintf(
    "  location: %s trees, one position each | with raw positions that disagree: %s | without coordinates: %d\n",
    format(uniqueN(rt$TreeID), big.mark = ","), format(votes[, .N, by = TreeID][N > 1L, .N], big.mark = ","), length(no_xy)
  ))
})

test_that("real data: DFstatus is exactly the raw Status (never modified)", {
  skip_if_not(real_ok, "no R tables / stage-2 file found")
  y <- rt[, .(StemID, census, DFstatus)][cells[, .(StemID, census, Status)], on = .(StemID, census)]
  bad <- y[!same_chr(as.character(DFstatus), Status), unique(StemID)]
  expect(length(bad) == 0L, sprintf("DFstatus differs from the raw Status for %d stems, e.g. %s", length(bad), show_ids(bad)))
})
