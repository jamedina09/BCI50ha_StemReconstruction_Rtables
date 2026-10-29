# =============================================================================
# load-bci-tables.R — loads the files of the BCI 50-ha data package and
# unpacks the posterior samples of the stem reconstruction. Base R only.
#
#   source("code/load-bci-tables.R")
#   bci <- load_bci_tables("data")      # list: bci.stem1 ... bci.stem9, bci.spptable
#   bci.stem9 <- bci$bci.stem9          # one data.frame per census
#   load("data/bci.stem9.Rdata")        # or one table at a time
#   paths <- load_bci_posterior("data") # posterior samples (one row per sample)
#   one <- draw_bci_reconstruction(paths)
#   ind <- bci_individuals(bci)         # state of every individual in every census
#   bci_individual_counts(ind)          # alive, dead, recruits, deaths per census
#
# Text columns stay text, so tags keep their leading zeros and signs
# (e.g. "000006", "-05599"); ExactDate is a Date.
# =============================================================================

.load_object <- function(path) {
    if (!file.exists(path)) stop("File not found: ", path)
    e <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    get(nm[1], envir = e)
}

# The nine census tables (bci.stem1 ... bci.stem9) and the species table
# (bci.spptable), as data.frames named as the ForestGEO R tables.
load_bci_tables <- function(data_dir = "data", censuses = 1:9) {
    out <- lapply(censuses, function(k) .load_object(file.path(data_dir, sprintf("bci.stem%d.Rdata", k))))
    names(out) <- paste0("bci.stem", censuses)
    out$bci.spptable <- .load_object(file.path(data_dir, "bci.spptable.rdata"))
    out
}

# Posterior samples of the stem reconstruction: one row per distinct sampled
# reconstruction of a tree (tag, treeID, path_sig, path_count, path_prob,
# recon). With long = TRUE, one row per record and sample: treeID, path (row
# of the samples), path_count, StemPaths, label.
load_bci_posterior <- function(data_dir = "data", long = FALSE) {
    x <- as.data.frame(readRDS(file.path(data_dir, "posterior_sampled_paths.rds")))
    if (!long) return(x)
    pairs <- strsplit(x$recon, ";", fixed = TRUE)
    n <- lengths(pairs)
    kv <- strsplit(unlist(pairs, use.names = FALSE), ":", fixed = TRUE)
    data.frame(
        treeID = rep(x$treeID, n),
        path = rep(seq_len(nrow(x)), n),
        path_count = rep(x$path_count, n),
        StemPaths = as.integer(vapply(kv, `[`, "", 1L)),
        label = vapply(kv, `[`, "", 2L),
        stringsAsFactors = FALSE
    )
}

# Draws one sampled reconstruction per tree, with probability
# path_count / sum(path_count) within the tree, and returns its records with
# their stem label in that draw (treeID, StemPaths, label). Join to a census
# table on treeID + StemPaths; records with the same treeID and label are one
# stem in this draw. The samples cover each tree's records up to its anchor
# census (2010 for most trees); a later record, fixed by its stem tag, belongs
# to the sampled stem of the last sampled record with the same stemID.
# Repeat the draw to propagate the uncertainty.
draw_bci_reconstruction <- function(posterior) {
    rows <- unlist(lapply(split(seq_len(nrow(posterior)), posterior$treeID), function(i) {
        if (length(i) == 1L) i else sample(i, 1L, prob = posterior$path_count[i])
    }), use.names = FALSE)
    one <- posterior[rows, ]
    pairs <- strsplit(one$recon, ";", fixed = TRUE)
    kv <- strsplit(unlist(pairs, use.names = FALSE), ":", fixed = TRUE)
    data.frame(
        treeID = rep(one$treeID, lengths(pairs)),
        StemPaths = as.integer(vapply(kv, `[`, "", 1L)),
        label = vapply(kv, `[`, "", 2L),
        stringsAsFactors = FALSE
    )
}

# ---- Individuals: alive, broken below or dead in every census ---------------
# An individual is one free-standing woody plant (treeID: a tree, shrub, palm,
# tree fern or strangler fig), whatever its number of stems. bci_individuals()
# puts every individual in one state in every census, from the census tables:
#   "stem at 1.3 m"    alive, with at least one living stem at 1.3 m: Rstatus A
#                      on a record that is not a stump, nor a census just before
#                      a stump (below). These are the individuals counted in a
#                      census.
#   "broken below"     alive, with no living stem at 1.3 m but a stump record:
#                      DFstatus "broken below" without a diameter (the field
#                      crew recorded the plant alive, broken below 1.3 m).
#   "alive, not seen"  alive, but missing or recorded dead in this census: it
#                      has a living stem or a stump record in a later census.
#   "dead"             recorded alive before, with no living stem and no stump
#                      record in this census or any later one (every stem D).
#   "not yet"          before its first living stem or stump record.
#   "never alive"      only dead or missing records: from its first record on.
# Alive = "stem at 1.3 m" + "broken below" + "alive, not seen". A census just
# before a stump of the same stem, with no sign of a stem at 1.3 m there (no
# "alive" record, and no diameter on a record that is not dead), has no living
# stem at 1.3 m: the stem broke between its last measurement and the stump.
# Needs the nine census tables, because a state can depend on later censuses
# (as Rstatus does): for some censuses only, subset the result. Returns a
# data.frame with one row per individual and census: treeID, tag, sp,
# CensusID, state (a factor) and first_census (the first census in which it
# is alive; NA if never). Stops if a state disagrees with Rstatus.
BCI_INDIVIDUAL_STATES <- c("stem at 1.3 m", "broken below", "alive, not seen", "dead", "not yet", "never alive")
bci_individuals <- function(bci) {
    censuses <- 1:9
    tabs <- bci[paste0("bci.stem", censuses)]
    if (any(vapply(tabs, is.null, TRUE))) {
        stop("bci_individuals() needs the nine census tables: load them with load_bci_tables() and its default censuses")
    }
    t1 <- tabs[[1]]
    if (!all(vapply(tabs, function(t) identical(t$treeID, t1$treeID) && identical(t$stemID, t1$stemID), TRUE))) {
        stop("The census tables must hold the same stems in the same order")
    }
    n_c <- length(censuses)
    # [stem x census]
    R <- sapply(tabs, function(t) t$Rstatus)
    DF <- sapply(tabs, function(t) t$DFstatus)
    has_dbh <- sapply(tabs, function(t) !is.na(t$dbh))
    is_in <- function(m, v) matrix(m %in% v, nrow(m))
    stump <- is_in(DF, "broken below") & !has_dbh
    at13 <- is_in(DF, "alive") | (has_dbh & !is_in(DF, c("dead", "stem dead")))
    # the stem's next event after each census: 1 a stump, 2 a sign of a stem at 1.3 m
    nxt <- matrix(0L, nrow(R), n_c)
    if (n_c > 1L) {
        for (j in (n_c - 1L):1L) nxt[, j] <- ifelse(stump[, j + 1L], 1L, ifelse(at13[, j + 1L], 2L, nxt[, j + 1L]))
    }
    stem13 <- R == "A" & !stump & (at13 | nxt != 1L)
    # [individual x census]
    ids <- unique(t1$treeID)
    tree <- match(t1$treeID, ids)
    by_tree <- function(m) rowsum(m * 1L, tree, reorder = FALSE) > 0
    S13 <- by_tree(stem13)
    STU <- by_tree(stump)
    sign <- S13 | STU
    has_sign <- rowSums(sign) > 0
    # alive from its first sign, or its first A stem if that comes earlier
    first <- max.col((sign | by_tree(R == "A")) * 1L, ties.method = "first")
    last <- n_c + 1L - max.col(sign[, n_c:1L, drop = FALSE] * 1L, ties.method = "first")
    first_rec <- max.col(by_tree(!is.na(DF)) * 1L, ties.method = "first")
    cc <- matrix(seq_len(n_c), length(ids), n_c, byrow = TRUE)
    state <- matrix("dead", length(ids), n_c)
    state[cc < last] <- "alive, not seen"
    state[STU] <- "broken below"
    state[S13] <- "stem at 1.3 m"
    state[cc < first] <- "not yet"
    never <- !has_sign
    state[never, ] <- ifelse(cc[never, , drop = FALSE] < first_rec[never], "not yet", "never alive")
    # the states agree with Rstatus: an individual with an A stem is alive; a
    # dead one has every stem D (or P)
    alive <- matrix(state %in% BCI_INDIVIDUAL_STATES[1:3], nrow(state))
    if (any(by_tree(R == "A") & !alive)) stop("An individual with an A stem is not alive")
    if (any(by_tree(R == "A" | R == "G") & state == "dead")) stop("A dead individual has a stem that is not D")
    one <- !duplicated(t1$treeID)
    data.frame(
        treeID = rep(ids, each = n_c),
        tag = rep(t1$tag[one], each = n_c),
        sp = rep(t1$sp[one], each = n_c),
        CensusID = rep(censuses, times = length(ids)),
        state = factor(as.vector(t(state)), levels = BCI_INDIVIDUAL_STATES),
        first_census = rep(ifelse(has_sign, censuses[first], NA_integer_), each = n_c),
        stringsAsFactors = FALSE
    )
}

# Counts per census from bci_individuals(): alive individuals (all, and by
# state), dead individuals (cumulative), recruits (individuals alive for the
# first time; none for the first census), deaths (alive in the previous census,
# dead in this one), individuals not yet alive and never alive, and all
# individuals. Takes bci_individuals() for all censuses or for consecutive
# ones. Stops unless the numbers add up: the states hold every individual in
# every census; alive = stem at 1.3 m + broken below + alive, not seen; dead
# grows by the deaths; alive + dead grows by the recruits.
bci_individual_counts <- function(ind) {
    tab <- unclass(table(ind$CensusID, ind$state))
    k <- as.integer(rownames(tab))
    if (any(diff(k) != 1L)) stop("bci_individual_counts() needs consecutive censuses")
    prev <- ind[order(ind$treeID, ind$CensusID), c("treeID", "CensusID", "state", "first_census")]
    was_alive <- c(FALSE, prev$state[-nrow(prev)] %in% BCI_INDIVIDUAL_STATES[1:3] & prev$treeID[-nrow(prev)] == prev$treeID[-1L])
    died <- tapply(was_alive & prev$state == "dead", prev$CensusID, sum)
    rec <- tapply(!is.na(prev$first_census) & prev$first_census == prev$CensusID, prev$CensusID, sum)
    out <- data.frame(
        CensusID = k,
        alive = tab[, "stem at 1.3 m"] + tab[, "broken below"] + tab[, "alive, not seen"],
        alive_stem_1.3m = tab[, "stem at 1.3 m"],
        alive_broken_below = tab[, "broken below"],
        alive_not_seen = tab[, "alive, not seen"],
        dead = tab[, "dead"],
        recruits = c(NA, rec[-1L]),
        deaths = c(NA, died[-1L]),
        not_yet = tab[, "not yet"],
        never_alive = tab[, "never alive"],
        individuals = rowSums(tab),
        row.names = NULL
    )
    n_ind <- length(unique(ind$treeID))
    stopifnot(
        "Every individual is in one state in every census" = all(out$individuals == n_ind),
        "Dead individuals grow by the deaths" = all(diff(out$dead) == out$deaths[-1L]),
        "Alive + dead individuals grow by the recruits" = all(diff(out$alive + out$dead) == out$recruits[-1L])
    )
    out
}
