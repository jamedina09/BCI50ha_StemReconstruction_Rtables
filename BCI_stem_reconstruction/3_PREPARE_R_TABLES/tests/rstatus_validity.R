# ========================================================================
# VALIDITY OF RSTATUS SEQUENCES (one stem) AND COMBINATIONS (one tree)
# ========================================================================
# Stem level: valid = P* A* G* D* with at least one non-P census.
# Tree level (Steps 2-3): the tree is alive at census c iff some stem is A at
# c or later; a G needs the tree alive, a D needs it dead. The two checks in
# the request ("any D => all non-P stems D", "any A => no D") are reported
# separately: they miss "G with no stem alive then or later".
# Flags (valid, but reported): a stem dead without ever being A, a G/D in
# census 1, and a P stem while its tree is dead.
# ========================================================================

stem_reasons <- function(s) {
  r <- character(length(s))
  add <- function(r, hit, txt) ifelse(hit, ifelse(nzchar(r), paste(r, txt, sep = "; "), txt), r)
  r <- add(r, !grepl("[AGD]", s), "all P (stem never recorded)")
  r <- add(r, grepl("[AGD]P", s), "P after a record")
  r <- add(r, grepl("[GD].*A", s), "A after G/D (false death not corrected)")
  r <- add(r, grepl("D.*G", s), "G after D (dead tree comes back)")
  r
}

stem_flags <- function(s) {
  r <- character(length(s))
  add <- function(r, hit, txt) ifelse(hit, ifelse(nzchar(r), paste(r, txt, sep = "; "), txt), r)
  r <- add(r, grepl("[GD]", s) & !grepl("A", s), "dead without ever being A")
  r <- add(r, grepl("^[GD]", s), "G/D in census 1")
  r
}

# trees: list of character vectors (one string per stem), all of one length.
classify_trees <- function(trees) {
  s <- unlist(trees, use.names = FALSE)
  tree <- rep(seq_along(trees), lengths(trees))
  n <- nchar(s[1])
  stopifnot(all(nchar(s) == n))
  M <- do.call(rbind, strsplit(s, "", fixed = TRUE))
  sr <- stem_reasons(s)
  sf <- stem_flags(s)
  L <- data.table(tree = rep(tree, n), c = rep(seq_len(n), each = length(s)), R = as.vector(M))
  L[, `:=`(isA = R == "A", isD = R == "D", isG = R == "G", isP = R == "P")]
  last_A <- L[isA == TRUE, .(last_A = max(c)), by = tree]
  L <- last_A[L, on = "tree"]
  L[is.na(last_A), last_A := 0L]
  L[, `:=`(g_dead = isG & c > last_A, d_alive = isD & c <= last_A, p_dead = isP & c > last_A & last_A > 0L)]
  tc <- L[, .(nA = sum(isA), nD = sum(isD), nG = sum(isG), g_dead = sum(g_dead), d_alive = sum(d_alive), p_dead = sum(p_dead)), by = .(tree, c)]
  tc[, `:=`(D_beside_live = nD > 0L & (nA + nG) > 0L, A_and_D = nA > 0L & nD > 0L)]
  tr <- tc[, .(
    G_tree_dead = sum(g_dead) > 0L, D_tree_alive = sum(d_alive) > 0L, P_tree_dead = sum(p_dead) > 0L,
    D_beside_live = sum(D_beside_live) > 0L, A_and_D = sum(A_and_D) > 0L
  ), by = tree]
  setorder(tr, tree)
  st <- data.table(tree = tree, bad = nzchar(sr), txt = sr, flag = nzchar(sf))
  st <- st[, .(bad = sum(bad) > 0L, flag = sum(flag) > 0L, txt = paste(unique(txt[nzchar(txt)]), collapse = "; ")), by = tree]
  setorder(st, tree)
  reasons <- paste0(
    st$txt,
    ifelse(tr$G_tree_dead, "; G while no stem is A then or later (tree dead)", ""),
    ifelse(tr$D_tree_alive, "; D while a stem is A then or later (tree alive)", ""),
    ifelse(tr$D_beside_live, "; [request check] D beside a non-P stem that is not D", ""),
    ifelse(tr$A_and_D, "; [request check] A and D in the same census", "")
  )
  list(
    valid = !st$bad & !tr$G_tree_dead & !tr$D_tree_alive & !tr$D_beside_live & !tr$A_and_D,
    flag = st$flag | tr$P_tree_dead,
    request_checks_invalid = st$bad | tr$D_beside_live | tr$A_and_D,
    reasons = sub("^; ", "", reasons)
  )
}
