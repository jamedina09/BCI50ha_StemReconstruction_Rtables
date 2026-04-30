# TrueStemID Hard Invariant — Analysis & POST Code Changes

**Status:** Drafted while PRE-change baseline runs are in progress.
**Goal:** Make `TrueStemID` a hard end-to-end invariant: **whenever `TrueStemID`
is non-NA on a row, the final `ReconstructedStemID` for that row MUST equal
`TrueStemID`, with `ReconstructionMethod = "given"`. No exceptions, in either
the full DP or the probabilistic-greedy fallback.**

The DP is only allowed to *infer* `ReconstructedStemID` for rows where
`TrueStemID` is NA. For all other rows it must respect the supplied identity.

---

## 1. What we did so far (Step 1 + Step 2)

In [dp_global/scripts/main_cpp_bci.R](dp_global/scripts/main_cpp_bci.R) (single-tag debug
driver) the input data are pre-processed to populate `TrueStemID` from
biologically certain evidence:

- **Step 1a** — `StemTag` non-NA: physical field tag → `TrueStemID := OriginalStemID`.
- **Step 1b** — `CensusID >= 7` (BCI re-tagging campaign year): the database
  ID is trustworthy → `TrueStemID := OriginalStemID`.
- **Step 2a** — compute `last_dbh_census` per `(Tag, OriginalStemID)`.
- **Step 2b** — direct anchor of terminal-event rows (`Status %in% {dead, stem
  dead, broken below}` or R-family code `R|RP|RF|RT|QR|OR` in `ListOfTSM`) when
  strictly after `last_dbh_census`. No prior anchor required (handles pure
  pre-C7 stems and spans-C7 with gaps before C7).
- **Step 2c** — bidirectional `nafill` LOCF/NOCB inside the post-last-DBH zone
  only, by `(Tag, OriginalStemID)`. Carries 2b/Step-1 anchors forward and a
  later C7+ anchor backward into terminal-phase rows.
- **Step 2d** — drop the temporary `.last_dbh_census` column.

**Pre-last-DBH rows are deliberately left as `NA`** so the DP can resolve them
freely.

### Critical gap (not yet addressed by Step 2)

A stem can carry an `R`/`broken below`/`dead` row that is itself **the
terminal anchor** for a long pre-anchor live trajectory:

```
C1: DBH=23.9  TrueStemID=NA          (live, but not yet anchored)
C2: DBH=NA    TrueStemID=4769  (R)   (death/resprout — anchored by Step 2b)
... [no C3..C6 records for this stem]
C7: stem absent (no re-tag)
```

Step 2 anchors only the **terminal NA-DBH row** (C2). It does **not** push
`TrueStemID` backward onto the live C1 row.

This is the conceptual missing piece: *the terminal state defines the identity
of the entire pre-terminal live trajectory* — but since that propagation
requires a biological-rules judgement (does the live C1 belong to the same
stem as the terminal C2?), the cleanest answer is **delegate it to the DP** by
making sure the DP has the necessary information (the terminal `TrueStemID`)
**in its track set** so its solution can naturally pick it up.

---

## 2. How the DP and probabilistic matcher (mis)handle TrueStemID today

### 2.1 Pin mechanism overview

Both engines build a fixed list of "tracks" (potential stem identities) and
then try to match each observation row to a track. When a row has a
non-`NA` `TrueStemID` and `pin_truestemid=TRUE`, the match index is forced.

### 2.2 BUG #1 — Pre-anchor `TrueStemID` silently dropped from track set

In [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R):

- **L963**: `anchor_obs <- tree_data[CensusID == anchor_start & !is.na(DBH)]`
- **L997**: `anchor_ids <- sort(unique(anchor_obs$TrueStemID))` — *anchor census only*
- **L1258**:
  ```r
  track_ids <- c(anchor_ids,
                 if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra)
                 else integer(0))
  ```

The `track_ids` are built from `anchor_ids` (anchor-census IDs only) plus
synthetic recruit slots — **pre-anchor `TrueStemID` values are never added**.

Then at **L1327** the pin map is built via `match(tsid, track_ids)` and
`NA` results are silently dropped at **L1338+**.

**Consequence**: a pre-anchor row with a perfectly valid `TrueStemID` (e.g.
`4769`) that does *not* exist at the anchor census produces `match(4769,
track_ids) = NA`, so the row is left unpinned and the DP is free to assign it
any track. Final `ReconstructedStemID != TrueStemID`. Invariant violated.

### 2.3 BUG #2 — NA-DBH terminal rows invisible to the DP, never receive `ReconstructedStemID`

The DP filters observations to `!is.na(DBH)` (`obs_row_idx[[p]]` builder):
NA-DBH rows are excluded from the state machine and the back-trace at
[dp_global/R/dp_global_dp.R L2060](dp_global/R/dp_global_dp.R#L2060):

```r
tree_data[obs_idx, ReconstructedStemID := track_ids[sv]]
```

only writes to indexed (non-NA-DBH) rows.

Existing post-anchor handling (L2086–2088) covers post-anchor rows with a
known `TrueStemID`:

```r
.pinned_post_unassigned <- which(!is.na(tree_data$TrueStemID) &
                                 tree_data$CensusID > anchor_start &
                                 is.na(tree_data$ReconstructedStemID))
tree_data[.pinned_post_unassigned, ReconstructedStemID := as.integer(TrueStemID)]
```

…but there is **no symmetric block for pre-anchor NA-DBH rows**. A row with
`TrueStemID=4769`, `CensusID=2`, `DBH=NA`, `Status="dead"` ends up with
`ReconstructedStemID=NA` after the DP. Invariant violated.

### 2.4 BUG #3 — Probabilistic matcher mirrors BUG #1

In [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R):

- **L143–157**: `anchor_ids` built from anchor-census obs only.
- **L184**: `tidx <- match(as.integer(tsid), anchor_ids)` — silent drop for
  pre-anchor TrueStemIDs absent at anchor.

Same failure as BUG #1, on the fallback path.

### 2.5 BUG #4 — `main_cpp.R` and `main_cpp_chunk.R` lack Step 1+2 reconstruction

`grep_search` confirms that **only [main_cpp_bci.R](dp_global/scripts/main_cpp_bci.R) builds `TrueStemID`**.
The simulated-data driver [main_cpp.R](dp_global/scripts/main_cpp.R) and the
chunked production driver [main_cpp_chunk.R](dp_global/scripts/main_cpp_chunk.R)
expect `TrueStemID` to come pre-populated from input.

The simulated CSV has a `TrueStemID` column (column 4 of
`data_simulation/data/simulated_data_1.csv`) so the simulator handles its own
ground truth — main_cpp.R needs no patch for the simulated data.

For BCI runs through the chunked driver (production path), the same Step 1+2
preprocessing is needed; otherwise the DP has nothing to pin against. We can
either (a) port Step 1+2 into main_cpp_chunk.R behind a switch (e.g. `if (any
column == "StemTag")`), or (b) preprocess the BCI RDS once and cache the
result. Recommendation: (a), guarded by the presence of the BCI-specific
columns (`StemTag`, `ListOfTSM`, `Status`).

### 2.6 BUG #5 — Post-anchor block in main_cpp.R/main_cpp_chunk.R requires `!is.na(DBH)`

```r
.post[!is.na(TrueStemID) & !is.na(DBH), `:=`(
    ReconstructedStemID = as.integer(TrueStemID),
    ReconstructionMethod = "given"
)]
```

This is the runner-side counterpart of BUG #2 for post-anchor rows. Drop the
`& !is.na(DBH)` filter to honour the invariant on terminal NA-DBH rows.

---

## 3. Patches to apply (POST phase)

### Patch A — Extend DP `track_ids` with pre-anchor `TrueStemID` (FIXES BUG #1)

**File:** [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R) around L1255–1258

**Replace:**
```r
n_extra <- K - length(anchor_ids)
current_max <- suppressWarnings(max(tree_data$TrueStemID, na.rm = TRUE))
if (!is.finite(current_max)) current_max <- 0
track_ids <- c(anchor_ids, if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0))
```

**With:**
```r
# Hard-invariant extension: every pre-anchor TrueStemID must also appear in
# track_ids so that match(tsid, track_ids) succeeds for non-anchor pinned rows.
# Without this, pre-anchor TrueStemIDs that are absent at the anchor census
# (e.g. a stem that died before the BCI re-tagging campaign at C7) are
# silently dropped from the pin map and the DP is free to violate the
# invariant ReconstructedStemID == TrueStemID.
.pre_anchor_tsids <- sort(unique(
    tree_data$TrueStemID[
        !is.na(tree_data$TrueStemID) &
        tree_data$CensusID < anchor_start &
        !(tree_data$TrueStemID %in% anchor_ids)
    ]
))

# Make sure K covers anchor_ids plus all pre-anchor TrueStemIDs
.required_K <- length(anchor_ids) + length(.pre_anchor_tsids)
if (.required_K > K) {
    K <- .required_K
}

n_extra <- K - length(anchor_ids) - length(.pre_anchor_tsids)
current_max <- suppressWarnings(max(c(tree_data$TrueStemID, anchor_ids, .pre_anchor_tsids), na.rm = TRUE))
if (!is.finite(current_max)) current_max <- 0
track_ids <- c(
    anchor_ids,
    .pre_anchor_tsids,
    if (n_extra > 0L) seq.int(from = current_max + 1L, length.out = n_extra) else integer(0)
)
```

> **Note:** the existing slack/anchor split logic (`.anchor_track_set <-
> which(track_ids %in% anchor_ids)` and `.slack_track_set <- which(!(track_ids
> %in% anchor_ids))`) at L1273–1274 is **preserved** because pre-anchor
> TrueStemIDs are *not* in `anchor_ids` and so they correctly land in the
> slack track set (they are empty at the anchor census, like recruits/deaths).

### Patch B — Assign `ReconstructedStemID` to pre-anchor NA-DBH rows after DP (FIXES BUG #2)

**File:** [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R), after the
post-anchor block at ~L2088, **before** the NA-R barrier post-processing.

**Insert:**
```r
# Hard-invariant pre-anchor sweep: every pre-anchor row with a non-NA
# TrueStemID that the DP did not assign (e.g. NA-DBH terminal rows that were
# filtered out of obs_row_idx) is set to its TrueStemID. This complements
# the .pinned_post_unassigned block above and guarantees that ReconstructedStemID
# == TrueStemID on every row where TrueStemID is non-NA, regardless of DBH.
if (isTRUE(pin_truestemid)) {
    .pinned_pre_unassigned <- which(
        !is.na(tree_data$TrueStemID) &
        tree_data$CensusID < anchor_start &
        is.na(tree_data$ReconstructedStemID)
    )
    if (length(.pinned_pre_unassigned) > 0L) {
        tree_data[.pinned_pre_unassigned, `:=`(
            ReconstructedStemID = as.integer(TrueStemID),
            ReconstructionMethod = "given"
        )]
    }
}
```

> **Ordering note:** must run *before* the NA-R barrier block (~L2099+) so
> that the barrier sever logic sees these rows already labelled with their
> correct identity (the existing `is.na(TrueStemID) | provisional_dp` guards
> already protect pinned rows).

### Patch C — Mirror Patch A in the probabilistic matcher (FIXES BUG #3)

**File:** [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R) after the `anchor_ids` block (~L160), **before** the `K <- max(length(anchor_ids), max_obs)` line.

**Insert:**
```r
# Hard-invariant extension (mirrors dp_global_dp.R Patch A): include
# every pre-anchor TrueStemID in the anchor_ids set used by the pin map,
# so match(as.integer(tsid), anchor_ids) succeeds for non-anchor pinned
# rows whose TrueStemID is absent at the anchor census.
.pre_anchor_tsids_prob <- sort(unique(
    tree_data$TrueStemID[
        !is.na(tree_data$TrueStemID) &
        tree_data$CensusID < anchor_pos_census &  # anchor census variable used in this fn
        !(tree_data$TrueStemID %in% anchor_ids)
    ]
))
if (length(.pre_anchor_tsids_prob) > 0L) {
    anchor_ids <- c(anchor_ids, as.integer(.pre_anchor_tsids_prob))
}
```

> The exact name of the anchor-census variable inside this function is to be
> confirmed at apply time (the snippet I read uses `obs_data[[i]]$census_id`
> for the per-iteration census; the anchor census itself is `census_range[n_census]`
> or similar). Will pin the exact identifier at edit time.

### Patch D — Pre-anchor NA-DBH sweep in the probabilistic matcher (FIXES BUG #3 cont.)

**File:** [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R), at the end, after the existing pin-based "given" labelling (~L397+), before the function returns.

**Insert** the same sweep as Patch B (substituting whatever variable names the
probabilistic matcher uses for the data table — likely `tree_data` here too).

### Patch E — Step 1+2 in the chunked production driver (FIXES BUG #4)

**File:** [dp_global/scripts/main_cpp_chunk.R](dp_global/scripts/main_cpp_chunk.R), inside `run_main_chunked` after the input is loaded but **before** the species/bio columns are attached.

**Logic:** detect BCI inputs by the presence of `StemTag` + `Status` +
`ListOfTSM` columns; if all three are present and `TrueStemID` is missing or
has < 1% non-NA, then run the same Step 1+2 reconstruction as
`main_cpp_bci.R` lines 184–249. Otherwise leave `TrueStemID` untouched.

This change is **deferred** (out of scope for the single-tag BCI debugger
study). It is required for production runs against full BCI but does not
affect the simulated-data comparison.

### Patch F — Drop the `!is.na(DBH)` guard in the post-anchor block (FIXES BUG #5)

**Files:**
- [dp_global/scripts/main_cpp.R](dp_global/scripts/main_cpp.R) L816
- [dp_global/scripts/main_cpp_chunk.R](dp_global/scripts/main_cpp_chunk.R) L829

**Replace:**
```r
.post[!is.na(TrueStemID) & !is.na(DBH), `:=`(
    ReconstructedStemID = as.integer(TrueStemID),
    ReconstructionMethod = "given"
)]
```

**With:**
```r
.post[!is.na(TrueStemID), `:=`(
    ReconstructedStemID = as.integer(TrueStemID),
    ReconstructionMethod = "given"
)]
```

> The DP-engine equivalent in `propagate_post_anchor_given()` (L478–481)
> already does the right thing — only these two outer-driver code paths need
> updating.

---

## 4. Algorithmic invariants the patches must preserve

After applying A–F, the following must hold *by construction*, not by luck:

| Invariant | Where enforced |
|-----------|----------------|
| `TrueStemID` non-NA  =>  `ReconstructedStemID == TrueStemID`, method `given` | Patches A+B for pre-anchor; Patches A+B equivalents already in code for anchor and post-anchor |
| Pre-anchor `TrueStemID` participates in the DP track set | Patch A |
| Pre-anchor NA-DBH `TrueStemID` rows receive their identity in output | Patch B |
| Probabilistic fallback enforces the same invariant | Patches C+D |
| `K` (track count) is large enough to host all pinned identities | Patch A's `.required_K` guard |
| Slack/anchor track partition is preserved (recruits/deaths still get free slots) | Patch A keeps `track_ids %in% anchor_ids` semantics; pre-anchor TrueStemIDs land in slack as desired |
| NA-R barrier sever logic does not destroy a Patch-B pin | Existing guards `is.na(TrueStemID) | provisional_dp` already exclude pinned rows |
| Post-anchor NA-DBH stumps with R/BB code keep `given` | Existing `propagate_post_anchor_given` (L490–497) already does this |

---

## 5. Risks, edge cases, and open questions

1. **State-space inflation.** Patch A grows `K` by the number of distinct
   pre-anchor TrueStemIDs not present at the anchor. For BCI tags this can
   double or triple the state count and may push more tags through the
   `K_too_small` / `state_too_large` fallback to the probabilistic matcher.
   **Mitigation:** that's the correct behaviour — the probabilistic matcher
   is also fixed (Patches C+D) to honour the invariant. The pre/post
   comparison study will quantify how often this occurs.

2. **Anchor-OK gating.** The slack-tracks block (L1208) sets `anchor_ok`
   based on whether at least one anchor-census stem is recruitable. Patch A
   adds tracks that are *not* anchored, but these are pre-anchor TrueStemIDs
   that *predate* the anchor — they should not affect the recruit-DBH check.
   The split `anchor_track_set / slack_track_set` puts them in `slack_track_set`,
   which already has the right semantics (always allowed at any census).

3. **`provisional_dp` interaction.** When the anchor census has missing
   TrueStemIDs and `ALLOW_PROVISIONAL_DP_ANCHOR=TRUE`, the code at L975–989
   fabricates anchor IDs into `tree_data$TrueStemID`. Patch A's
   `.pre_anchor_tsids` query happens *after* this fabrication, so the
   fabricated anchor IDs are not included again (they are in `anchor_ids`).
   Safe.

4. **Duplicate TrueStemID across `(Tag, OriginalStemID)` groups.** Within
   one Tag we should not see the same TrueStemID on two different
   OriginalStemIDs (data invariant from the upstream cleaning). If it
   ever happens, `match()` will pick the first track and the DP will be
   forced into a single-track-multiple-rows configuration that may be
   infeasible. Recommend a one-time `validate_truestemid_uniqueness()` check
   right after Step 2d.

5. **`USE_MEASUREMENT_ERROR` with extended track set.** Posterior bins are
   computed per track; extending the track set increases compute slightly
   (linear in K). No correctness risk.

6. **NA-R barrier sever, "anchor-validated crossings" filter (L2120–2143).**
   This already tests `tree_data$TrueStemID[CensusID > .cc_bar ...]`. With
   Patch B writing `ReconstructedStemID` on pre-anchor NA-DBH rows, this
   filter is unchanged in semantics (it's about post-barrier knowns).
   Safe.

7. **`finalize_out` re-stamps `TrueStemID` for `provisional_dp` rows
   (L519–525).** This restores fabricated anchors. Patch B happens *before*
   `finalize_out` is called, on rows where `TrueStemID` is genuine (Step
   1/2). The re-stamp loop only touches rows with `ReconstructionMethod ==
   "provisional_dp"` and so does not undo Patch B. Safe.

8. **Patch C variable name.** Pending verification of the anchor-census
   variable name inside `match_stems_probabilistic`. Will fix at edit time
   (likely `obs_data[[n_census]]$census_id` or just `anchor_start`).

9. **Patch E (chunked driver) is deferred.** Not required for this study;
   simulated CSV provides `TrueStemID` directly. Production BCI will need
   it before the next chunked production run.

---

## 6. Sequencing for the POST batch

Apply patches in this order (they are independent but B must run after A
in the file):

1. **Patch A** — extend `track_ids` (~L1255–1260 of dp_global_dp.R)
2. **Patch B** — pre-anchor NA-DBH sweep (~L2089 of dp_global_dp.R)
3. **Patch C** — extend `anchor_ids` in probabilistic (~L160 of
   dp_probabilistic_matching.R)
4. **Patch D** — mirror sweep at end of probabilistic
5. **Patch F** — drop `!is.na(DBH)` filter in main_cpp.R + main_cpp_chunk.R

Then re-run the same 30 commands and diff PRE vs POST.

---

## 7. Acceptance criteria for "fully stable code"

After POST runs:

- For every output CSV across all 30 runs (PRE and POST):
  - `nrow(out[!is.na(TrueStemID) & ReconstructedStemID != TrueStemID, ])` must
    be **0** in POST (was likely > 0 in PRE).
  - `nrow(out[!is.na(TrueStemID) & is.na(ReconstructedStemID), ])` must be
    **0** in POST.
  - `ReconstructionMethod` must be `"given"` on all those rows in POST.
- POST may show:
  - More tags falling back from full DP to probabilistic at `DP_MAX_STATES=10000`
    (because Patch A grows K). Check the `DP_FallbackReason` distribution.
  - Identical ground-truth-honouring results between full DP and probabilistic
    at `DP_MAX_STATES=2` (both routed to probabilistic, both honour the pin).

---

*End of analysis.*
