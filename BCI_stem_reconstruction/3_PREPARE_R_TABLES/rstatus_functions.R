# ========================================================================
# RSTATUS AND DBH RULES FOR THE BCI R TABLES
# ========================================================================
# Sourced by 2_create_R_tables_BCI.R and by tests/test_Rstatus.R, so the
# script and the tests run exactly the same code. One small function per
# step of the rules (the script calls them in this order, with its logs and
# checks in between; compute_rstatus() chains them for one table):
#   rstatus_raw_code()   raw field record -> A / D / N             (Section 4)
#   rstatus_propagate()  P before the first record; no record after
#                        A or D -> D                                (Section 6)
#   fix_resurrections()  alive later => never dead                  (Section 8)
#   rstatus_tree_dg()    dead stem -> G (tree alive) or D (tree dead) (Section 10)
#
# Rules:
#   Evidence of life : an "alive" record (with or without DBH), or a DBH on a
#                      "broken below", "missing" or status-less record.
#   Dead record      : "dead", "stem dead" (with or without DBH) and
#                      "broken below" without DBH. A DBH on a dead record is
#                      not evidence of life (the stem is A there only if it
#                      is alive later); the DBH itself is kept.
#   No record        : "missing" or no status without DBH, or no row.
#   P                : only before a stem's first record.
#   A                : from the first record to the last evidence of life
#                      (missed censuses and false deaths in between are A).
#   After the last evidence of life, or for a stem never alive: dead from the
#                      first census without evidence of life (G or D).
#   Tree alive at c  : some stem of the tree is A at c or at a later census.
#   G / D            : dead stem with its tree alive / dead at that census.
#   DBH              : exported as recorded on every record (A, G or D); a P
#                      cell never has one; never imputed.
# ========================================================================

# Raw field record -> "A" (evidence of life), "D" (dead record) or "N" (no
# record). status and dbh are vectors of the same length (one per stem x
# census); status NA means no status (or no row).
rstatus_raw_code <- function(status, dbh) {
  known <- c("alive", "dead", "stem dead", "broken below", "missing")
  unknown <- setdiff(unique(status[!is.na(status)]), known)
  if (length(unknown) > 0L) {
    stop("rstatus_raw_code: unknown raw status: ", paste(unknown, collapse = ", "))
  }
  has_dbh <- !is.na(dbh)
  bb <- status %in% "broken below"
  code <- rep("N", length(status))
  code[status %in% "alive"] <- "A"
  code[status %in% c("dead", "stem dead")] <- "D"
  code[bb & has_dbh] <- "A" # resprout measured
  code[bb & !has_dbh] <- "D" # broken stem, nothing measured
  code[(is.na(status) | status %in% "missing") & has_dbh] <- "A"
  code
}

# Section 6 propagation on encounter histories (one string per stem, one
# character per census, codes A / D / N):
#   ^N -> P, PN -> PP  P until the stem's first record
#   AN -> AD           no record after alive: dead from the first census
#                      without a record (undone in fix_resurrections() when
#                      the stem is alive later)
#   DN -> DD           no record after a dead record: stays dead
rstatus_propagate <- function(histories) {
  h <- sub("^N", "P", histories)
  while (any(grepl("PN", h, fixed = TRUE))) h <- gsub("PN", "PP", h, fixed = TRUE)
  while (any(grepl("AN", h, fixed = TRUE))) h <- gsub("AN", "AD", h, fixed = TRUE)
  while (any(grepl("DN", h, fixed = TRUE))) h <- gsub("DN", "DD", h, fixed = TRUE)
  h
}

# ------------------------------------------------------------------------
# HELPER FUNCTION: fix_resurrections
# ------------------------------------------------------------------------
# Section 8: resolve D->A / G->A transitions in encounter histories (see the
# PARAMETERS / RETURNS notes inside). The default log_file reads CHECK_folder,
# a global of 2_create_R_tables_BCI.R; pass log_file = NULL to write no log.
# Stops when the histories differ in length or dbh_matrix has another shape.
fix_resurrections <- function(status_vec,
                              dbh_matrix,
                              stem_ids = NULL,
                              dbh_aware = TRUE,
                              log_file = file.path(CHECK_folder, "log_resurrections.txt"),
                              verbose = TRUE) {
  # PARAMETERS
  #   status_vec : char vector of encounter histories (one per stem)
  #   dbh_matrix : numeric matrix [stem x census], NA = not measured
  #   stem_ids   : optional StemIDs for the log
  #   dbh_aware  : TRUE  -> per-cell decision uses DBH evidence
  #                FALSE -> lifespan rule: every DA/GA -> AA (backfill all);
  #                         the mode used by 2_create_R_tables_BCI.R and by
  #                         compute_rstatus()
  #   log_file   : audit log path (NULL: no log)
  #   verbose    : print a one-line summary
  #
  # RETURNS list(status_vec = corrected histories, n_backfill, n_demote,
  #   flagged = one row per D->A / G->A cell with the action taken).
  #
  # IMPLEMENTATION (single matrix pass + vectorized string ops):
  #   1. Build status matrix once.
  #   2. Find ALL (i, j) cells where smat[i, j-1] in {D,G} and smat[i, j] = A.
  #      Each such cell is an immediate D->A / G->A transition.
  #   3. Classify each cell by DBH at j:
  #        has DBH  -> backfill_to_A   (the A is real; earlier D/G is wrong)
  #        no  DBH  -> demote_to_DG    (the A is a zombie)
  #   4. Apply demotes by direct matrix-index assignment (single vectorized op).
  #   5. Apply backfills via gsub("DA|GA","AA",...) loop on the string vector
  #      (only on rows that actually need backfilling; vectorized C code).
  #   6. Write a single log of all classified cells (no per-pass overhead).
  stopifnot(is.character(status_vec))
  status_lengths <- nchar(status_vec)
  if (length(unique(status_lengths)) != 1L) {
    stop("fix_resurrections: all encounter histories must have the same length.")
  }
  n_censuses <- status_lengths[1]
  n_stems <- length(status_vec)
  if (!is.matrix(dbh_matrix) || nrow(dbh_matrix) != n_stems || ncol(dbh_matrix) != n_censuses) {
    stop(sprintf(
      "fix_resurrections: dbh_matrix must be %d x %d (got %d x %d).",
      n_stems, n_censuses, nrow(dbh_matrix), ncol(dbh_matrix)
    ))
  }
  if (is.null(stem_ids)) stem_ids <- as.character(seq_len(n_stems))

  # ---- single matrix scan -----------------------------------------------
  smat <- do.call(rbind, strsplit(status_vec, "", fixed = TRUE))

  prev <- smat[, 1:(n_censuses - 1), drop = FALSE]
  curr <- smat[, 2:n_censuses, drop = FALSE]
  hits <- which((prev == "D" | prev == "G") & curr == "A", arr.ind = TRUE)

  if (nrow(hits) == 0L) {
    if (verbose) cat("fix_resurrections: no D->A / G->A cells found.\n")
    return(list(
      status_vec = status_vec, n_backfill = 0L, n_demote = 0L,
      flagged = data.table()
    ))
  }

  cell_rows <- hits[, 1]
  cell_cols <- hits[, 2] + 1L # the "A" cell column
  prev_cols <- hits[, 2] # the "D"/"G" cell column
  prev_codes <- smat[cbind(cell_rows, prev_cols)]
  cell_dbh <- dbh_matrix[cbind(cell_rows, cell_cols)]
  has_dbh <- !is.na(cell_dbh)

  if (!dbh_aware) {
    action <- rep("backfill_to_A", length(cell_rows))
  } else {
    action <- ifelse(has_dbh, "backfill_to_A", "demote_to_DG")
  }
  is_back <- action == "backfill_to_A"
  is_dem <- action == "demote_to_DG"
  n_backfill <- sum(is_back)
  n_demote <- sum(is_dem)

  # ---- apply demotes (one vectorized matrix assignment) -----------------
  if (n_demote > 0L) {
    smat[cbind(cell_rows[is_dem], cell_cols[is_dem])] <- prev_codes[is_dem]
  }

  # ---- rebuild string vector ---------------------------------------------
  new_vec <- do.call(paste0, lapply(seq_len(n_censuses), function(j) smat[, j]))

  # ---- apply backfills via gsub on the rows that need it ----------------
  # After demotes, rows with backfill cells still contain "DA" or "GA".
  # gsub propagates DDDA -> DDAA -> DAAA -> AAAA in a small loop on a
  # subset of rows only.
  if (n_backfill > 0L) {
    back_rows <- unique(cell_rows[is_back])
    sub <- new_vec[back_rows]
    while (any(grepl("DA|GA", sub, fixed = FALSE))) {
      sub <- gsub("DA|GA", "AA", sub)
    }
    new_vec[back_rows] <- sub
  }

  # ---- build log --------------------------------------------------------
  flagged <- data.table(
    stem_idx = cell_rows,
    StemID   = stem_ids[cell_rows],
    census   = cell_cols,
    pattern  = paste0(prev_codes, "A"),
    DBH      = cell_dbh,
    action   = action
  )

  if (!is.null(log_file)) {
    con <- file(log_file, open = "w")
    on.exit(close(con), add = TRUE)
    writeLines(c(
      paste0("# log_resurrections.txt - generated ", format(Sys.time())),
      paste0("# mode            : ", if (dbh_aware) "DBH-aware" else "naive (gsub-equivalent)"),
      paste0("# n_stems         : ", n_stems),
      paste0("# n_censuses      : ", n_censuses),
      paste0(
        "# n_backfill_to_A : ", n_backfill,
        if (dbh_aware) {
          "  (D/G->A transitions backfilled; the A cell had a DBH)"
        } else {
          "  (D/G->A transitions backfilled; the whole D/G run before each A becomes A)"
        }
      ),
      paste0("# n_demote_to_DG  : ", n_demote, "  (A rewritten to D/G; cell had no DBH; DBH-aware mode only)"),
      "#"
    ), con)
    if (nrow(flagged) > 0L) {
      writeLines("# per-cell disposition (one row per D->A or G->A cell):", con)
      write.table(flagged, con,
        sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE
      )
    } else {
      writeLines("# no resurrection cells found.", con)
    }
  }

  if (verbose) {
    cat(sprintf(
      "fix_resurrections: mode=%s | backfill_to_A=%d | demote_to_DG=%d\n",
      if (dbh_aware) "DBH-aware" else "naive",
      n_backfill, n_demote
    ))
    if (n_demote > 0L) {
      cat(sprintf(
        "  ⚠ %d zombie A cell(s) demoted (no DBH after D/G). See: %s\n",
        n_demote, log_file
      ))
    }
  }

  list(
    status_vec = new_vec,
    n_backfill = n_backfill,
    n_demote = n_demote,
    flagged = flagged
  )
}

# Section 10 tree-aware D/G. status_matrix: [stem x census] codes A / D / G /
# P; tree_id: one value per stem. A dead cell (D or G) at census c becomes G
# when some stem of its tree is A at c or later (tree alive), D otherwise.
# A and P cells are never changed.
rstatus_tree_dg <- function(status_matrix, tree_id) {
  is_A <- (status_matrix == "A") * 1L
  stem_last_A <- ifelse(rowSums(is_A) > 0L, max.col(is_A, ties.method = "last"), 0L)
  tree_last_A <- ave(stem_last_A, tree_id, FUN = max)
  tree_alive <- col(status_matrix) <= tree_last_A
  dead <- status_matrix == "D" | status_matrix == "G"
  out <- status_matrix
  out[dead & tree_alive] <- "G"
  out[dead & !tree_alive] <- "D"
  out
}

# All steps for one table of stems (used by tests/test_Rstatus.R).
#   status, dbh : matrices [stem x census] of raw Status and DBH
#   tree_id     : one value per stem
# Returns list(Rstatus, dbh): matrices of the same shape; dbh is the raw DBH,
# exported unchanged (no DBH is removed or imputed).
compute_rstatus <- function(status, dbh, tree_id) {
  code <- matrix(rstatus_raw_code(as.vector(status), as.vector(dbh)), nrow = nrow(status))
  h <- do.call(paste0, as.data.frame(code, stringsAsFactors = FALSE))
  h <- rstatus_propagate(h)
  h <- fix_resurrections(h, dbh, dbh_aware = FALSE, log_file = NULL, verbose = FALSE)$status_vec
  R <- rstatus_tree_dg(do.call(rbind, strsplit(h, "", fixed = TRUE)), tree_id)
  list(Rstatus = R, dbh = dbh)
}
