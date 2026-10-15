# ========================================================================
# MEASUREMENT-DISCONTINUITY REJOIN (post-engine repair of stem identity)
# ========================================================================
# Why: when the point of measurement (POM) of a trunk moves (buttresses
# growing up the trunk, deformities, strangler figs), its recorded DBH can
# change by far more than any growth, whether the new height of measurement
# (HOM) was recorded or not. The taper correction cannot remove buttress
# flare, so the taper-corrected DBH can still drop by half. The engines'
# pruning removes every link outside the hard growth bounds [MAX_SHRINK,
# MAX_GROWTH] cm/yr before costs are compared, while a recruit above the
# recruit limit is only penalised, so the engine books a stem loss plus an
# impossible recruit for one physical trunk. The engines also cannot continue
# a stem across a census in which it has no measurement.
#
# What: in trees flagged as multi-stem (single_stem_tags == FALSE), join a
# stem s1 that starts at census c1 with a stem s0 of the same tree whose
# measurements end at c0 = c1 - 1, or at c0 = c1 - 2 (s0 not measured in the
# census between), only when the split has to be undone and a moved POM
# explains it. All of the following must hold.
#
# The split has to be undone
#   1. s1 is an impossible recruit: its first DBH is at or above the recruit
#      limit of the stage-2 run (BCI: 26 cm) in a tree measured before. Left
#      apart, one trunk counts as a stem loss plus a large new stem.
#   2. Clean handover: s0 has no row from c1 on and s1 none before c1, so no
#      census gets two rows.
#   3. The engine could not link them: between consecutive censuses the
#      taper-corrected growth is outside the hard bounds (a split the engine
#      chose inside the bounds is its own decision and is kept); across a
#      census without a measurement the engine cannot link at all.
#
# It is not a break or a death
#   4. The later measurement has no break / resprout status or code (rule R1:
#      a broken-below record with a DBH starts a new stem), and s0 has no
#      dead, stem-dead or broken-below record in the census between.
#
# The size change is the one a moved POM gives
#   5. A smaller DBH: the taper-corrected ratio d1 / d0 is at least ratio_min
#      (BCI: the joined POM drops reach 0.42; the candidates below 0.4 are
#      broken trunks, 0.13-0.33, and one strangler-fig stem at an unchanged
#      HOM, 0.39), and s0 was >= min_dbh_mm.
#   6. A larger DBH (growth above the bound) only when the recorded HOM
#      changed between the two measurements and the recorded DBH grew inside
#      the bounds: the taper correction at the new height, not the trunk,
#      broke the bound. A move up the trunk that was not recorded gives a
#      smaller DBH, never a larger one. The taper-corrected ratio stays at or
#      below ratio_max (BCI: the two joined rises are 1.47 and 1.48).
#   7. Across a census without a measurement, the recorded HOM rose and the
#      DBH is smaller: sizes alone cannot tell a moved POM from a trunk that
#      was only missed in one census, and a missed census is not a moved POM.
#      No other measurement of the tree in the census between may be the same
#      trunk (a stem not measured at c0 whose taper-corrected size is within
#      [ratio_min, ratio_max] of s0: the trunk was then measured, under
#      another stem or as forks).
#
# The tree shows that its POM moves (at least one)
#   8. a. a measurement of the tree was taken away from the standard height
#         (|HOM - hom_std| > hom_tol): its POM had to be moved;
#      b. a record of the tree carries the buttress code (buttress_code);
#      c. its growth form has no stable POM (irregular_forms; BCI: strangler
#         figs);
#      d. its species usually has buttresses (buttress_species) and neither
#         stem of the pair is measured in a census that recorded buttress
#         codes (coded_censuses; BCI: from 1990, so only 1982 and 1985). Where
#         the codes exist, the trunk's own record outranks its species: a
#         trunk measured then without a B code has no buttress. The species
#         list comes from mr_buttress_species().
#      The signs are read on the whole tree, in any census, because the stem
#      identities are what is being repaired. A census without HOM records is
#      not a sign on its own (BCI 1982: every HOM is 1.3 m): a drop from 1982
#      in a tree with no other sign cannot be told from a recording error or a
#      shrinking trunk (in three of the six BCI cases the 1990 DBH rises back
#      towards the 1982 value).
#
# One identity
#   9. The two stems carry at most one distinct TrueStemID pin.
#  10. The pairing is clearly the best one: among the pairs that pass 1-9, it
#      is the best pair of s1 and the best pair of s0, ranked by (a) same
#      height of measurement (|HOM1 - HOM0| <= hom_tol), (b) size ratio
#      closest to 1 (|log(d1 / d0)|); a pair is clearly best when it wins on
#      (a), or when its |log ratio| is below ratio_margin times that of the
#      next pair. Pairs of consecutive censuses are resolved first; pairs
#      across a census without a measurement are resolved afterwards, among
#      the stems not joined in the first step. Pairs that fail 1-9 (for
#      example a small sprout that ends in the same census) do not compete;
#      when no pair is clearly best (two similar forks after one trunk, two
#      similar trunks before one stem), nothing is joined.
#
# Not joined, by design (BCI: 005827, a species with buttresses but a trunk
# measured from 1990 to 2015 with no B code): recording errors and trunks missed in one census in
# trees without a sign of a moved POM, and any split inside the growth
# bounds. The database StemIDs are not used: before the anchor they are the
# identity the reconstruction replaces, and pins are used only to keep two
# pinned stems apart. Relabelled rows keep every column except
# ReconstructedStemID (they take the id of the stem carrying the pin, else of
# the later stem) and ReconstructionMethod (= "measurement_rejoin").
# Posterior samples: the observation pairs joined in the export are joined in
# each sample where they are split with a clean end/start and at most one
# pin. Idempotent: a second call finds nothing to join.
#
# Used by 2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R (export) and
# 3_PREPARE_R_TABLES/1_prepare_posteriors_BCI.R (posterior samples), both in
# BCI_stem_reconstruction/.
#
# Units: DBH and the taper-corrected DBH (dbh_with_best_candidate_taper_corrected)
# are in mm, as are recruit_max_mm and min_dbh_mm; max_shrink and max_growth
# are in cm/yr and are compared with the growth of the taper-corrected DBH
# (rule 6 also with that of the recorded DBH); HOM, hom_std and hom_tol are
# in m.
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

# Species whose trees usually have buttresses (rule 8 d): the share of trees that
# reached min_dbh_mm and were measured in a census with buttress codes
# (coded_censuses) carrying the code in any census, for species with at least
# min_trees such trees.
#   dt : stage-2 table (all tags) with Tag, CensusID, DBH, ListOfTSM, Mnemonic
mr_buttress_species <- function(dt, coded_censuses, buttress_code = "B", min_dbh_mm = 300,
                                min_trees = 10L, min_share = 0.8) {
    m <- dt[!is.na(DBH), .(Tag = as.character(Tag), CensusID = as.integer(as.character(CensusID)), DBH,
                           sp = as.character(Mnemonic), tsm = fcoalesce(as.character(ListOfTSM), ""))]
    m[, B := vapply(strsplit(tsm, ";", fixed = TRUE), function(v) any(v %in% buttress_code), NA)]
    tr <- m[, .(sp = sp[1L], dmax = max(DBH), coded = any(CensusID %in% coded_censuses), B = any(B)), by = Tag]
    tr <- tr[dmax >= min_dbh_mm & coded]
    s <- tr[, .(trees = .N, share = mean(B)), by = sp]
    s[trees >= min_trees & share >= min_share, sp]
}

# Rows of multi-stem trees that carry a ReconstructedStemID (measured or not),
# with what the rule needs. `dt` must hold the row index .mr_row.
#   buttress_code, irregular_forms : rule 8 b and c (see measurement_rejoin_candidates())
.mr_obs <- function(dt, buttress_code = "B", irregular_forms = NULL) {
    hom <- if ("HOM" %in% names(dt)) "HOM" else NULL
    b_regex <- if (length(buttress_code)) sprintf("\\b(%s)\\b", paste(buttress_code, collapse = "|")) else NULL
    w <- dt[single_stem_tags %in% FALSE & !is.na(ReconstructedStemID), .(
        row = .mr_row, Tag = as.character(Tag), c = as.integer(as.character(CensusID)),
        rs = as.character(ReconstructedStemID), d = DBH,
        tc = dbh_with_best_candidate_taper_corrected, t = as.numeric(ExactDate),
        R = Status %in% "broken below" | grepl(.mr_r_regex, fcoalesce(as.character(ListOfTSM), ""), perl = TRUE),
        ends = Status %in% c("dead", "stem dead", "broken below") | grepl(.mr_r_regex, fcoalesce(as.character(ListOfTSM), ""), perl = TRUE),
        B = if (is.null(b_regex)) FALSE else grepl(b_regex, fcoalesce(as.character(ListOfTSM), ""), perl = TRUE),
        pin = as.character(TrueStemID), obs = as.integer(obs_row_id)
    )]
    w[, h := if (is.null(hom)) NA_real_ else as.numeric(dt[w$row, get(hom)])]
    w[, sp := if ("Mnemonic" %in% names(dt)) as.character(dt$Mnemonic[w$row]) else NA_character_]
    w[, irregular := if (is.null(irregular_forms) || !"growth_form" %in% names(dt)) FALSE else
        grepl(irregular_forms, as.character(dt$growth_form[w$row]), ignore.case = TRUE)]
    w
}

# Eligibility of an observation pair (vectorised): returns "eligible" or the
# reason the pair is kept apart (rules 3-8 of the header, in the order tested).
#   d0, d1 : DBH (mm) of the two measurements; ratio : d1 / d0 (taper-corrected)
#   g_tc, g_raw : annual growth of the taper-corrected and of the recorded DBH (cm/yr)
#   R1 : TRUE when the later measurement has a break / resprout status or code
#   gap : c1 - c0 (1 or 2); ended_between : s0 has a dead / break record in between
#   raised : the recorded HOM of the later measurement is higher; moved : the
#   recorded HOM differs; measured_between : the tree has, in the census
#   between, a measurement that may be the same trunk; evidence : rule 8 holds
.mr_eligibility <- function(d0, d1, ratio, g_tc, g_raw, R1, gap, ended_between, raised, moved, measured_between, evidence,
                            max_shrink, max_growth, min_dbh_mm, ratio_min, ratio_max) {
    outside <- g_tc < max_shrink | g_tc > max_growth
    fcase(
        R1, "break/resprout code on the later measurement",
        ended_between, "death or break record in the census between",
        gap == 1L & !outside, "link inside the growth bounds",
        gap > 1L & !(raised & d1 < d0), "across a census without a measurement: no recorded HOM rise with a smaller DBH",
        gap > 1L & measured_between, "the trunk may be measured in the census between",
        d0 < min_dbh_mm, "earlier stem below min_dbh",
        ratio < ratio_min, "size ratio below ratio_min",
        ratio > ratio_max, "size ratio above ratio_max",
        g_tc > max_growth & !(moved & g_raw <= max_growth), "larger DBH not explained by a recorded HOM change",
        !evidence, "no evidence of a moved point of measurement",
        default = "eligible"
    )
}

# Within one group (the eligible pairs of one stem), TRUE for the pair that is
# clearly best (rule 10), FALSE for the others.
.mr_best <- function(gap, moved, dev, margin) {
    n <- length(gap)
    out <- logical(n)
    o <- order(gap, moved, dev)
    if (n == 1L) {
        out[o[1L]] <- TRUE
        return(out)
    }
    b <- o[1L]
    s <- o[2L]
    out[b] <- gap[b] < gap[s] || moved[b] < moved[s] || dev[b] < margin * dev[s]
    out
}

# Candidate pairs in the current table, with the reason each is joined or kept apart.
#   dt : stage-2 table with the row index .mr_row added by
#        apply_measurement_rejoin(), and Tag, CensusID, DBH, HOM (optional),
#        dbh_with_best_candidate_taper_corrected, ExactDate, Status, ListOfTSM,
#        TrueStemID, obs_row_id, single_stem_tags, ReconstructedStemID,
#        growth_form (optional, rule 8 c)
#   max_shrink, max_growth : hard growth bounds of the stage-2 run (cm/yr)
#   recruit_max_mm : recruit limit (mm); min_dbh_mm, ratio_min, ratio_max : rules 5-7
#   max_gap : largest c1 - c0 (2: one census without a measurement in between)
#   hom_tol (m) : HOM differences up to hom_tol count as the same height (rules 6-8, 10)
#   ratio_margin : rule 10
#   hom_std (m) : standard height of measurement (rule 8 a)
#   buttress_code : code(s) of ListOfTSM marking a buttress (rule 8 b); NULL: not used
#   irregular_forms : regular expression on growth_form for growth forms
#                     without a stable POM (rule 8 c); NULL: not used
#   buttress_species, coded_censuses : rule 8 d (species of mr_buttress_species();
#                     censuses that recorded buttress codes); NULL: not used. dt
#                     needs Mnemonic when buttress_species is given
# Returns one row per candidate pair (rules 1-2) with `why` ("join" or the
# reason), `evidence` (the signs of rule 8 found on the tree: "HOM",
# "buttress", "growth form", "species") and, for joins, `target` and `source` (the ids
# kept and replaced).
measurement_rejoin_candidates <- function(dt, max_shrink, max_growth, recruit_max_mm, min_dbh_mm = 100,
                                          ratio_min = 0.4, ratio_max = 1.5, max_gap = 2L,
                                          hom_tol = 0.05, ratio_margin = 0.5, hom_std = 1.3,
                                          buttress_code = "B", irregular_forms = NULL,
                                          buttress_species = NULL, coded_censuses = NULL) {
    w <- .mr_obs(dt, buttress_code, irregular_forms)
    span <- w[, .(any_first = min(c), any_last = max(c)), by = .(Tag, rs)]
    m <- w[!is.na(d)]
    fm <- m[m[, .I[which.min(c)], by = .(Tag, rs)]$V1]
    lm <- m[m[, .I[which.max(c)], by = .(Tag, rs)]$V1]
    pins <- w[!is.na(pin), .(pins = list(unique(pin))), by = .(Tag, rs)]
    # impossible recruits (rule 1) with a clean start, and every stem of the tree
    # whose measurements end one or max_gap censuses before, with a clean end
    tag_first <- m[, .(tag_first = min(c)), by = Tag]
    s1 <- span[fm[, .(Tag, rs, c1 = c, row1 = row, d1 = d)], on = .(Tag, rs)][tag_first, on = "Tag", nomatch = 0L]
    s1 <- s1[c1 > tag_first & d1 >= recruit_max_mm & any_first == c1]
    s0 <- span[lm[, .(Tag, rs, c0 = c, row0 = row)], on = .(Tag, rs)]
    P <- s0[, .(Tag, rs0 = rs, c0, row0, last0 = any_last)][s1[, .(Tag, rs1 = rs, c1, row1)], on = "Tag", allow.cartesian = TRUE, nomatch = 0L]
    P <- P[rs0 != rs1 & c0 < c1 & c1 - c0 <= max_gap & last0 < c1]
    if (!nrow(P)) {
        return(P[, .(Tag, row0, row1, route = character())])
    }
    P[, route := fifelse(c1 - c0 == 1L, "impossible_recruit", "impossible_recruit_gap")]
    # dead / break records of s0 in the census between (rule 4)
    stops <- w[is.na(d) & ends == TRUE, .(Tag, rs0 = rs, cb = c)]
    P[, ended_between := FALSE]
    if (nrow(stops)) P[stops, on = .(Tag, rs0, c0 < cb, c1 > cb), ended_between := TRUE]
    P <- P[, .(Tag, row0, row1, route, ended_between)]
    key <- w[, .(row, c, rs, d, h, tc, t, R, obs)]
    k0 <- copy(key); setnames(k0, paste0(names(k0), "0"))
    k1 <- copy(key); setnames(k1, paste0(names(k1), "1"))
    P <- k0[P, on = "row0"]
    P <- k1[P, on = "row1"]
    P <- span[, .(Tag, rs0 = rs, last0 = any_last)][P, on = .(Tag, rs0)]
    P <- span[, .(Tag, rs1 = rs, first1 = any_first)][P, on = .(Tag, rs1)]
    P <- pins[, .(Tag, rs0 = rs, p0 = pins)][P, on = .(Tag, rs0)]
    P <- pins[, .(Tag, rs1 = rs, p1 = pins)][P, on = .(Tag, rs1)]
    P[, `:=`(gap = c1 - c0, g_tc = (tc1 - tc0) / 10 / ((t1 - t0) / 365.25), g_raw = (d1 - d0) / 10 / ((t1 - t0) / 365.25), ratio = tc1 / tc0)]
    P[, n_pins := mapply(function(a, b) length(unique(c(unlist(a), unlist(b)))), p0, p1)]
    P[, `:=`(raised = !is.na(h0) & !is.na(h1) & h1 > h0 + hom_tol, moved = !is.na(h0) & !is.na(h1) & abs(h1 - h0) > hom_tol)]
    # rule 8: signs on the tree, in any census, that its POM moves
    ev <- w[, .(ev_hom = any(!is.na(d) & !is.na(h) & abs(h - hom_std) > hom_tol), ev_B = any(B), ev_form = any(irregular)), by = Tag]
    P <- ev[P, on = "Tag"]
    P[, ev_species := FALSE]
    if (length(buttress_species) && length(coded_censuses)) {
        cov <- m[, .(coded = any(c %in% coded_censuses)), by = .(Tag, rs)]
        P <- cov[, .(Tag, rs0 = rs, coded0 = coded)][P, on = .(Tag, rs0)]
        P <- cov[, .(Tag, rs1 = rs, coded1 = coded)][P, on = .(Tag, rs1)]
        spt <- w[, .(sp = sp[1L]), by = Tag]
        P <- spt[P, on = "Tag"]
        P[, ev_species := !coded0 & !coded1 & sp %in% buttress_species]
    }
    P[, evidence := paste0(fifelse(ev_hom, "HOM;", ""), fifelse(ev_B, "buttress;", ""), fifelse(ev_form, "growth form;", ""), fifelse(ev_species, "species;", ""))]
    P[, evidence := sub(";$", "", evidence)]
    # measurements of the tree in the census between, on stems not measured at c0
    P[, measured_between := FALSE]
    if (any(P$gap > 1L)) {
        on_c <- unique(m[, .(Tag, rs, c)])
        mid <- m[d >= min_dbh_mm, .(Tag, rsb = rs, cb = c, tcb = tc)]
        G <- P[gap > 1L, .(Tag, rs0, rs1, c0, c1, tc0)][mid, on = .(Tag), allow.cartesian = TRUE, nomatch = 0L][cb > c0 & cb < c1 & rsb != rs0 & rsb != rs1]
        G <- G[!on_c[, .(Tag, rsb = rs, c0 = c)], on = .(Tag, rsb, c0)]
        G <- G[tcb / tc0 >= ratio_min & tcb / tc0 <= ratio_max, .(Tag, rs0, rs1)]
        if (nrow(G)) P[unique(G), on = .(Tag, rs0, rs1), measured_between := TRUE]
    }
    P[, why := .mr_eligibility(d0, d1, ratio, g_tc, g_raw, R1, gap, ended_between, raised, moved, measured_between, nzchar(evidence),
                               max_shrink, max_growth, min_dbh_mm, ratio_min, ratio_max)]
    P[why == "eligible" & n_pins > 1L, why := "two different pins"]
    # rule 10: the clearly best pairing among the eligible pairs; pairs of
    # consecutive censuses first, then pairs across a census without a
    # measurement, among the stems not joined in the first step
    P[, dev := abs(log(ratio))]
    for (g in sort(unique(P$gap))) {
        if (g > 1L) {
            done <- P[why == "join"]
            P[why == "eligible" & gap == g & (paste(Tag, rs0) %in% paste(done$Tag, done$rs0) | paste(Tag, rs1) %in% paste(done$Tag, done$rs1)),
              why := "another pairing fits better"]
        }
        sel <- P[, why == "eligible" & gap == g]
        if (!any(sel)) next
        P[sel, best1 := .mr_best(gap, moved, dev, ratio_margin), by = .(Tag, rs1)]
        P[sel, best0 := .mr_best(gap, moved, dev, ratio_margin), by = .(Tag, rs0)]
        P[sel, `:=`(clear1 = any(best1)), by = .(Tag, rs1)]
        P[sel, `:=`(clear0 = any(best0)), by = .(Tag, rs0)]
        P[sel, why := fcase(
            best1 & best0, "join",
            (best1 | !clear1) & (best0 | !clear0), "ambiguous (no clearly best pairing)",
            default = "another pairing fits better"
        )]
    }
    P[why == "join", target := fifelse(lengths(p0) > 0L & lengths(p1) == 0L, rs0, rs1)]
    P[why == "join", source := fifelse(target == rs0, rs1, rs0)]
    P[, intersect(c("best1", "best0", "clear1", "clear0"), names(P)) := NULL]
    P[]
}

# Stage-2 table -> list(dt = repaired table, pairs = joined pairs, candidates = first pass,
# obs_pairs = the joined observation pairs, for the posterior samples).
# The input is not modified. Joins are applied in passes (at most max_iter):
# a join can make a new candidate pair, and the loop ends when a pass joins
# nothing. `verbose` prints the number of joins and the reasons of the first pass.
# The arguments after recruit_max_mm are those of measurement_rejoin_candidates().
apply_measurement_rejoin <- function(dt, max_shrink, max_growth, recruit_max_mm, min_dbh_mm = 100,
                                     ratio_min = 0.4, ratio_max = 1.5, max_gap = 2L,
                                     hom_tol = 0.05, ratio_margin = 0.5, hom_std = 1.3,
                                     buttress_code = "B", irregular_forms = NULL,
                                     buttress_species = NULL, coded_censuses = NULL,
                                     max_iter = 10L, verbose = TRUE) {
    dt <- data.table::copy(data.table::as.data.table(dt))
    dt[, .mr_row := .I]
    rs_class <- class(dt$ReconstructedStemID)
    pairs <- list()
    first <- NULL
    for (it in seq_len(max_iter)) {
        P <- measurement_rejoin_candidates(dt, max_shrink, max_growth, recruit_max_mm,
            min_dbh_mm = min_dbh_mm, ratio_min = ratio_min, ratio_max = ratio_max, max_gap = max_gap,
            hom_tol = hom_tol, ratio_margin = ratio_margin, hom_std = hom_std,
            buttress_code = buttress_code, irregular_forms = irregular_forms,
            buttress_species = buttress_species, coded_censuses = coded_censuses
        )
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
            "[measurement_rejoin] %d stem pair(s) joined in %d tree(s) (impossible recruit; %d across a census without a measurement)\n",
            nrow(pairs), data.table::uniqueN(pairs$Tag), if (nrow(pairs)) sum(pairs$route == "impossible_recruit_gap") else 0L
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
