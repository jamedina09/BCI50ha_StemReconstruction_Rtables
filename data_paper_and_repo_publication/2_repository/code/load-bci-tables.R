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
