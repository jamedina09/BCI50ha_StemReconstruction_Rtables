# ========================================================================
# REFERENCE IMPLEMENTATION OF THE RSTATUS AND DBH RULES
# ========================================================================
# Written to be read, not to be fast: one tree, one stem and one census at a
# time, straight from the written rules. test_Rstatus.R checks that
# rstatus_functions.R (what 2_create_R_tables_BCI.R runs) and the exported R
# tables give exactly the same Rstatus and dbh.
#
# Raw record -> evidence:
#   alive (with or without DBH)                  -> alive
#   broken below WITH DBH                        -> alive (resprout measured)
#   broken below without DBH, dead, stem dead    -> dead (with or without DBH)
#   missing or no status WITH DBH                -> alive
#   missing or no status without DBH, no row     -> none (not recorded)
# Step 1 (each stem):
#   P before the stem's first record;
#   A from the first record to the last evidence of life (missed censuses and
#     false deaths in between, including a dead first record, are A);
#   dead after the last evidence of life; a stem never alive is dead from its
#     first record; a stem never recorded is P in every census.
# Step 2: the tree is alive at c if some stem is A at c or at a later census.
# Step 3: a dead stem is G when its tree is alive, D when it is dead.
# Step 4: dbh = the raw DBH, as recorded (a P cell never has one).
#
# ref_evidence(status, dbh): one record -> "alive", "dead" or "none"; stops
#   on an unknown raw status.
# ref_tree(status, dbh): one tree -> list(Rstatus, dbh), matrices [stem x census].
# ref_rstatus_table(cells, cores): many trees (see below). Needs data.table
#   (loaded by test_Rstatus.R); cores > 1 uses parallel::mclapply.
# ========================================================================

ref_evidence <- function(status, dbh) {
  has_dbh <- !is.na(dbh)
  if (is.na(status)) {
    return(if (has_dbh) "alive" else "none")
  }
  switch(status,
    "alive" = "alive",
    "broken below" = if (has_dbh) "alive" else "dead",
    "dead" = "dead",
    "stem dead" = "dead",
    "missing" = if (has_dbh) "alive" else "none",
    stop("unknown raw status: ", status)
  )
}

# One tree: status and dbh are matrices [stem x census].
ref_tree <- function(status, dbh) {
  n_stems <- nrow(status)
  n <- ncol(status)
  # Step 1
  life <- matrix("P", n_stems, n)
  for (s in seq_len(n_stems)) {
    evidence <- character(n)
    for (c in seq_len(n)) evidence[c] <- ref_evidence(status[s, c], dbh[s, c])
    records <- which(evidence != "none")
    if (length(records) == 0L) next # never recorded: P in every census
    first_record <- min(records)
    alive_at <- which(evidence == "alive")
    last_alive <- if (length(alive_at) > 0L) max(alive_at) else first_record - 1L
    for (c in seq_len(n)) {
      if (c < first_record) {
        life[s, c] <- "P"
      } else if (c <= last_alive) {
        life[s, c] <- "A"
      } else {
        life[s, c] <- "dead"
      }
    }
  }
  # Step 2
  tree_alive <- logical(n)
  for (c in seq_len(n)) tree_alive[c] <- any(life[, c:n] == "A")
  # Steps 3 and 4
  rstatus <- matrix("P", n_stems, n)
  out_dbh <- matrix(NA_real_, n_stems, n)
  for (s in seq_len(n_stems)) {
    for (c in seq_len(n)) {
      if (life[s, c] == "dead") {
        rstatus[s, c] <- if (tree_alive[c]) "G" else "D"
      } else {
        rstatus[s, c] <- life[s, c]
      }
      out_dbh[s, c] <- dbh[s, c]
    }
  }
  list(Rstatus = rstatus, dbh = out_dbh)
}

# Many trees. cells: data.table(TreeID, StemID, census, Status, DBH) with one
# row per stem x census. Returns data.table(StemID, census, Rstatus_ref, dbh_ref).
ref_rstatus_table <- function(cells, cores = 1L) {
  cells <- cells[order(TreeID, StemID, census)]
  n <- max(cells$census)
  stems <- unique(cells[, .(TreeID, StemID)])
  stopifnot(nrow(cells) == nrow(stems) * n)
  St <- matrix(cells$Status, ncol = n, byrow = TRUE)
  Db <- matrix(cells$DBH, ncol = n, byrow = TRUE)
  rows <- split(seq_len(nrow(stems)), stems$TreeID)
  chunks <- split(seq_along(rows), cut(seq_along(rows), max(1L, min(length(rows), 8L * cores)), labels = FALSE))
  run_chunk <- function(k) {
    lapply(rows[k], function(r) ref_tree(St[r, , drop = FALSE], Db[r, , drop = FALSE]))
  }
  res <- if (cores > 1L) parallel::mclapply(chunks, run_chunk, mc.cores = cores) else lapply(chunks, run_chunk)
  res <- unlist(res, recursive = FALSE)
  R <- matrix("", nrow(stems), n)
  D <- matrix(NA_real_, nrow(stems), n)
  r_all <- unlist(rows[unlist(chunks)], use.names = FALSE)
  R[r_all, ] <- do.call(rbind, lapply(res, `[[`, "Rstatus"))
  D[r_all, ] <- do.call(rbind, lapply(res, `[[`, "dbh"))
  data.table(
    StemID = rep(stems$StemID, n), census = rep(seq_len(n), each = nrow(stems)),
    Rstatus_ref = as.vector(R), dbh_ref = as.vector(D)
  )
}
