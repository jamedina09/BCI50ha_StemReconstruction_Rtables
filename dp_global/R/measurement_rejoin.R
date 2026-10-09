# ========================================================================
# MEASUREMENT-DISCONTINUITY REJOIN (post-engine repair of stem identity)
# ========================================================================
# Why: when the point of measurement (POM) of a trunk moves, its recorded DBH
# can change by far more than any growth: BCI 1982 diameters of buttressed
# trees were taken at 1.3 m around the buttresses (1982 has no HOM records),
# and from 1985 the POM was raised above them (recorded, or left at 1.3 m by
# mistake). The taper correction cannot remove buttress flare, so the
# taper-corrected DBH still drops 20-50 %. Recording errors do the same (one
# census with a wrong DBH). The engines' pruning removes every link outside
# the hard growth bounds [MAX_SHRINK, MAX_GROWTH] cm/yr before costs are
# compared, while a recruit above the recruit limit is only penalised, so the
# engine books a death plus an impossible recruit for one physical trunk.
#
# What: in trees flagged as multi-stem (single_stem_tags == FALSE), join an
# ended stem s0 (last measurement at census c) with a stem s1 of the same
# tree that starts at c + 1, when
#   1. s1 is an impossible recruit: it starts at or above the recruit limit
#      (no new stem starts that large) in a tree measured before;
#   2. s0 is the only stem of the tree whose measurements end at c (the only
#      candidate);
#   3. the link is outside the hard growth bounds (the engine could not
#      consider it);
#   4. the later measurement has no break / resprout status or code (rule R1:
#      a broken-below record with a DBH starts a new stem, so it is never
#      joined to the stem before it);
#   5. s0 was >= min_dbh_mm and the taper-corrected size ratio d1 / d0 lies in
#      [ratio_min, ratio_max] (BCI: buttress cases reach 0.49, broken trunks
#      and resprouts stay below 0.3; a much larger "continuation" is a
#      different stem);
#   6. clean end/start: s0 has no row after c and s1 none before c + 1 (no
#      census gets two rows), the two stems carry at most one distinct
#      TrueStemID pin, and neither stem takes part in another candidate join.
# The database StemIDs are not used: before the anchor they are the identity
# the reconstruction replaces, and pins are used only to keep two pinned stems
# apart. A split the engine chose inside the growth bounds is left as it is.
# Relabelled rows keep every column except ReconstructedStemID (they take the
# id of the stem carrying the pin, else of the later stem) and
# ReconstructionMethod (= "measurement_rejoin"). Posterior samples: the
# observation pairs joined in the export are joined in each sample where they
# are split with a clean end/start and at most one pin. Idempotent: a second
# call finds nothing to join.
#
# Used by 2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R (export) and
# 3_PREPARE_R_TABLES/1_prepare_posteriors_BCI.R (posterior samples), both in
# BCI_stem_reconstruction/.
#
# Units: DBH and the taper-corrected DBH (dbh_with_best_candidate_taper_corrected)
# are in mm, as are recruit_max_mm and min_dbh_mm; max_shrink and max_growth
# are in cm/yr and are compared with the growth of the taper-corrected DBH.
# ========================================================================

# Visible check (prints ✓ / ❌; on failure also warns, then stops).
#   ok       : TRUE when the check passes
#   msg      : what is checked
#   examples : optional values shown on failure (first 10)
#   n_bad    : optional number of failing cases shown on failure
mr_check <- function(ok, msg, examples = NULL, n_bad = NULL) {
    ok <- isTRUE(ok)
    if (ok) {
        cat("✓", msg, "\n")
        return(invisible(TRUE))
    }
    txt <- paste0(
        "CHECK FAILED: ", msg, if (!is.null(n_bad)) paste0(" [", n_bad, " case(s)]") else "",
        if (length(examples)) paste0(" | examples: ", paste(utils::head(examples, 10), collapse = ", ")) else ""
    )
    cat("❌", txt, "\n")
    warning(txt, call. = FALSE, immediate. = TRUE)
    stop(txt, call. = FALSE)
}

.mr_r_regex <- "\\b(R|RP|RF|RT|QR|OR)\\b"

# Rows of multi-stem trees that carry a ReconstructedStemID (measured or not),
# with what the rule needs. `dt` must hold the row index .mr_row.
.mr_obs <- function(dt) {
    dt[single_stem_tags %in% FALSE & !is.na(ReconstructedStemID), .(
        row = .mr_row, Tag = as.character(Tag), c = as.integer(as.character(CensusID)),
        rs = as.character(ReconstructedStemID), d = DBH,
        tc = dbh_with_best_candidate_taper_corrected, t = as.numeric(ExactDate),
        R = Status %in% "broken below" | grepl(.mr_r_regex, fcoalesce(as.character(ListOfTSM), ""), perl = TRUE),
        pin = as.character(TrueStemID), obs = as.integer(obs_row_id)
    )]
}

# Eligibility of an observation pair (vectorised): returns "join" or the
# reason the pair is kept apart (rules 3-5 of the header, in the order tested).
#   d0 : DBH (mm) of the earlier measurement; ratio : d1 / d0 (taper-corrected)
#   g_tc : annual growth of the taper-corrected DBH (cm/yr)
#   R1 : TRUE when the later measurement has a break / resprout status or code
.mr_eligibility <- function(d0, ratio, g_tc, R1, max_shrink, max_growth, min_dbh_mm, ratio_min, ratio_max) {
    outside <- g_tc < max_shrink | g_tc > max_growth
    fcase(
        R1, "break/resprout code on the later measurement",
        !outside, "link inside the growth bounds",
        d0 < min_dbh_mm, "earlier stem below min_dbh",
        ratio < ratio_min, "size ratio below ratio_min",
        ratio > ratio_max, "size ratio above ratio_max",
        default = "join"
    )
}

# Candidate pairs in the current table, with the reason each is joined or kept apart.
#   dt : stage-2 table with the row index .mr_row added by
#        apply_measurement_rejoin(), and Tag, CensusID, DBH,
#        dbh_with_best_candidate_taper_corrected, ExactDate, Status, ListOfTSM,
#        TrueStemID, obs_row_id, single_stem_tags, ReconstructedStemID
#   max_shrink, max_growth : hard growth bounds of the stage-2 run (cm/yr)
#   recruit_max_mm : recruit limit (mm); min_dbh_mm, ratio_min, ratio_max : rule 5
# Returns one row per candidate pair (rules 1-2) with `why` ("join" or the
# reason) and, for joins, `target` and `source` (the ids kept and replaced).
measurement_rejoin_candidates <- function(dt, max_shrink, max_growth, recruit_max_mm, min_dbh_mm = 100,
                                          ratio_min = 0.4, ratio_max = 1.5) {
    w <- .mr_obs(dt)
    span <- dt[single_stem_tags %in% FALSE & !is.na(ReconstructedStemID),
        .(any_first = min(as.integer(as.character(CensusID))), any_last = max(as.integer(as.character(CensusID)))),
        by = .(Tag = as.character(Tag), rs = as.character(ReconstructedStemID))
    ]
    m <- w[!is.na(d)]
    fm <- m[m[, .I[which.min(c)], by = .(Tag, rs)]$V1]
    lm <- m[m[, .I[which.max(c)], by = .(Tag, rs)]$V1]
    pins <- w[!is.na(pin), .(pins = list(unique(pin))), by = .(Tag, rs)]
    # impossible recruit with exactly one stem of the tree ending just before
    tag_first <- m[, .(tag_first = min(c)), by = Tag]
    s1 <- span[fm[, .(Tag, rs, c1 = c, row1 = row, d1 = d)], on = .(Tag, rs)][tag_first, on = "Tag", nomatch = 0L]
    s1 <- s1[c1 > tag_first & d1 >= recruit_max_mm & any_first == c1]
    s0 <- span[lm[, .(Tag, rs, c0 = c, row0 = row)], on = .(Tag, rs)][any_last == c0]
    B <- s0[, .(Tag, c1 = c0 + 1L, row0)][s1[, .(Tag, c1, row1)], on = .(Tag, c1), nomatch = 0L]
    P <- B[B[, .I[.N == 1L], by = .(Tag, row1)]$V1, .(Tag, row0, row1, route = "impossible_recruit")]
    if (!nrow(P)) return(P)
    key <- w[, .(row, c, rs, d, tc, t, R, obs)]
    k0 <- copy(key); setnames(k0, paste0(names(k0), "0"))
    k1 <- copy(key); setnames(k1, paste0(names(k1), "1"))
    P <- k0[P, on = "row0"]
    P <- k1[P, on = "row1"]
    P <- span[, .(Tag, rs0 = rs, last0 = any_last)][P, on = .(Tag, rs0)]
    P <- span[, .(Tag, rs1 = rs, first1 = any_first)][P, on = .(Tag, rs1)]
    P <- pins[, .(Tag, rs0 = rs, p0 = pins)][P, on = .(Tag, rs0)]
    P <- pins[, .(Tag, rs1 = rs, p1 = pins)][P, on = .(Tag, rs1)]
    P[, `:=`(g_tc = (tc1 - tc0) / 10 / ((t1 - t0) / 365.25), ratio = tc1 / tc0)]
    P[, n_pins := mapply(function(a, b) length(unique(c(unlist(a), unlist(b)))), p0, p1)]
    P[, why := .mr_eligibility(d0, ratio, g_tc, R1, max_shrink, max_growth, min_dbh_mm, ratio_min, ratio_max)]
    P[why == "join" & rs0 == rs1, why := "already joined"]
    P[why == "join" & (last0 != c0 | first1 != c1), why := "not a clean end/start"]
    P[why == "join" & n_pins > 1L, why := "two different pins"]
    P[why == "join", `:=`(n0 = .N), by = .(Tag, rs0)][why == "join", `:=`(n1 = .N), by = .(Tag, rs1)]
    P[why == "join" & (n0 > 1L | n1 > 1L), why := "ambiguous (stem in two joins)"]
    P[why == "join", target := fifelse(lengths(p0) > 0L & lengths(p1) == 0L, rs0, rs1)]
    P[why == "join", source := fifelse(target == rs0, rs1, rs0)]
    P[]
}

# Stage-2 table -> list(dt = repaired table, pairs = joined pairs, candidates = first pass,
# obs_pairs = the joined observation pairs, for the posterior samples).
# The input is not modified. Joins are applied in passes (at most max_iter):
# a join can make a new candidate pair, and the loop ends when a pass joins
# nothing. `verbose` prints the number of joins and the reasons of the first pass.
apply_measurement_rejoin <- function(dt, max_shrink, max_growth, recruit_max_mm, min_dbh_mm = 100,
                                     ratio_min = 0.4, ratio_max = 1.5,
                                     max_iter = 10L, verbose = TRUE) {
    dt <- data.table::copy(data.table::as.data.table(dt))
    dt[, .mr_row := .I]
    rs_class <- class(dt$ReconstructedStemID)
    pairs <- list()
    first <- NULL
    for (it in seq_len(max_iter)) {
        P <- measurement_rejoin_candidates(dt, max_shrink, max_growth, recruit_max_mm, min_dbh_mm, ratio_min, ratio_max)
        if (is.null(first)) first <- P
        J <- if (nrow(P)) P[why == "join"] else P
        if (!nrow(J)) break
        J[, pass := it]
        pairs[[it]] <- J
        for (k in seq_len(nrow(J))) {
            sel <- which(as.character(dt$Tag) == J$Tag[k] & as.character(dt$ReconstructedStemID) == J$source[k])
            data.table::set(dt, i = sel, j = "ReconstructedStemID", value = methods::as(J$target[k], rs_class))
            data.table::set(dt, i = sel, j = "ReconstructionMethod", value = "measurement_rejoin")
        }
    }
    dt[, .mr_row := NULL]
    pairs <- data.table::rbindlist(pairs, fill = TRUE)
    # the posterior samples get the same observation pairs as the export
    obs_pairs <- if (nrow(pairs)) unique(pairs[, .(Tag, obs0, obs1, c0, c1, route)]) else
        data.table::data.table(Tag = character(), obs0 = integer(), obs1 = integer(), c0 = integer(), c1 = integer(), route = character())
    if (isTRUE(verbose)) {
        cat(sprintf(
            "[measurement_rejoin] %d stem pair(s) joined in %d tree(s) (impossible recruit with a single candidate)\n",
            nrow(pairs), data.table::uniqueN(pairs$Tag)
        ))
        if (!is.null(first) && nrow(first)) print(first[, .N, by = .(why)][order(-N)])
    }
    list(dt = dt, pairs = pairs, candidates = first, obs_pairs = obs_pairs)
}

# Posterior samples: join each eligible observation pair in every path where
# it is split with a clean end/start and the two labels carry at most one pin
# (paths hold measured observations only).
#   paths      : data.table(tag, treeID, path_sig, path_count, path_prob, recon)
#   obs_pairs  : data.table(Tag, obs0, obs1, c0, c1), from apply_measurement_rejoin()
#   obs_info   : data.table(tag, obs, c, pin), census and TrueStemID of every measured observation
# Returns a copy of `paths` keyed by tag. In the trees of obs_pairs, paths that
# become identical are collapsed (path_count summed, path_sig of the first
# one) and path_prob is recomputed as path_count / sum(path_count).
apply_measurement_rejoin_to_paths <- function(paths, obs_pairs, obs_info, verbose = TRUE) {
    paths <- data.table::copy(paths)
    if (!nrow(obs_pairs)) return(paths)
    pr <- obs_pairs[, .(tag = as.character(Tag), obs0 = as.integer(obs0), obs1 = as.integer(obs1), c0 = as.integer(c0), c1 = as.integer(c1))]
    setorder(pr, tag, c0)
    hit <- paths[tag %in% pr$tag]
    n_join <- 0L
    out <- lapply(split(hit, by = "tag"), function(tp) {
        prt <- pr[tag == tp$tag[1]]
        oi <- obs_info[tag == tp$tag[1]]
        tp[, recon := vapply(strsplit(recon, ";", fixed = TRUE), function(r) {
            kv <- do.call(rbind, strsplit(r, ":", fixed = TRUE))
            obs <- as.integer(kv[, 1])
            lab <- kv[, 2]
            cc <- oi$c[match(obs, oi$obs)]
            pn <- oi$pin[match(obs, oi$obs)]
            for (q in seq_len(nrow(prt))) {
                i0 <- match(prt$obs0[q], obs)
                i1 <- match(prt$obs1[q], obs)
                if (is.na(i0) || is.na(i1) || lab[i0] == lab[i1]) next
                l0 <- lab == lab[i0]
                l1 <- lab == lab[i1]
                if (any(cc[l0] > prt$c0[q]) || any(cc[l1] < prt$c1[q])) next
                if (length(unique(stats::na.omit(pn[l0 | l1]))) > 1L) next
                lab[l1] <- lab[i0]
                n_join <<- n_join + 1L
            }
            paste(obs, lab, sep = ":", collapse = ";")
        }, character(1))]
        tp
    })
    hit <- data.table::rbindlist(out)
    # identical partitions collapse into one path; sample counts add up
    hit <- hit[, .(path_sig = path_sig[1L], path_count = sum(path_count)), by = .(tag, treeID, recon)]
    hit[, path_prob := path_count / sum(path_count), by = tag]
    paths <- rbind(paths[!tag %in% pr$tag], hit[, names(paths), with = FALSE])
    data.table::setkeyv(paths, "tag")
    if (isTRUE(verbose)) {
        cat(sprintf("[measurement_rejoin] posterior: %d join(s) inside the samples of %d tree(s)\n", n_join, data.table::uniqueN(hit$tag)))
    }
    paths
}
