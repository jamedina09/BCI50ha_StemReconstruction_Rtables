## Reverse the direction of `ReconstructedStemID` numbering

### Background and motivation

In the BCI dataset (and in the convention shared by most ForestGEO plots),
`OriginalStemID` was assigned **forward in time as new stems were
encountered**: a stem first measured in census 1 received a smaller integer
than a stem first measured in census 5. The numerical order of
`OriginalStemID` therefore mirrors the chronological order of stem
appearance.

The reconstruction engine works **backward in time** from a late-census
anchor (default `ANCHOR_START_CENSUS = 7`). Stem identities are derived
from the anchor and propagated to earlier censuses.

---

### Where `ReconstructedStemID` values are minted today

There are five independent minting sites, all of which use ascending
`max + 1L` allocation:

**1. `dp_global_dp.R` — provisional anchor IDs**
When the anchor census has DBH rows but some or all `TrueStemID` values
are missing, the engine assigns provisional IDs:

```r
current_max <- max(tree_data$TrueStemID, na.rm = TRUE)
prov_ids <- seq.int(from = current_max + 1L, length.out = n)
```

These IDs go into both `TrueStemID` and `ReconstructedStemID` on the
anchor rows, tagged `ReconstructionMethod = "provisional_dp"`.

**2. `dp_global_dp.R` — extra slack track IDs (`track_ids`)**
The DP needs `K` tracks where `K = n_anchor_stems + n_pre_anchor_pins +
n_slack`. The slack tracks (one per dead/recruit stem that has no anchor
identity) receive fresh IDs:

```r
current_max <- max(c(TrueStemID, anchor_ids, .pre_anchor_tsids), na.rm=TRUE)
track_ids <- c(anchor_ids, .pre_anchor_tsids,
               seq.int(from = current_max + 1L, length.out = n_extra))
```

These track IDs become the `ReconstructedStemID` values for pre-anchor
observations that the DP assigns to those slots.

**3. `dp_global_dp.R` — NA-R (resprout) barrier splitting**
When a track crosses a hard resprout boundary (all live stems gone at
census *c*, NA-DBH R-coded rows present), the pre-barrier segment is
severed and given a new ID:

```r
.cur_max_id <- max(tree_data$ReconstructedStemID, na.rm=TRUE)
for (.old_id in .crossing) {
    .new_id <- .cur_max_id + 1L
    .cur_max_id <- .new_id
    ...
}
```

**4. `dp_global_dp.R` — segment-split offset**
When a tag is too complex to solve as one unit, it is split into
sub-segments. After the pre-segment is solved, its `ReconstructedStemID`
values are shifted upward by the maximum post-segment ID to avoid
collision:

```r
.offset <- .max_post_id
ReconstructedStemID := ReconstructedStemID + .offset
# DP_PosteriorTop{k}ID columns are offset by the same value
```

**5. `dp_global_main.R` — `apply_broken_below_invariants()`**
After the engine returns, broken-below life-cycle invariants are enforced.
Rows that must be split into a new trajectory receive:

```r
.new_id <- .cur_max_id + 1L
.cur_max_id <- .new_id
```

These IDs are tagged `ReconstructionMethod ∈
{bb_split, bb_split_carry, bb_post_terminator_split,
bb_post_terminator_split_carry}`.

**Critical side-effect at minting site 5:** when the invariant enforcement
overrides a pre-stamped `TrueStemID` pin (the BCI driver pre-stamps
`TrueStemID = OriginalStemID` on every BB+DBH row), the function also
updates `TrueStemID` to the newly minted ID so the two columns stay
consistent:

```r
if (!is.na(trueid[j]) && trueid[j] == old_id) trueid[j] <- new_id
...
if (has_true) out[, TrueStemID := trueid]
```

This means that, after `apply_broken_below_invariants()` runs, some rows
have a bb-minted integer in their `TrueStemID` column — **not** a real
BCI database ID. The `known_ids` computation in the renumbering function
must account for this.

**Net effect today:** all five minting sites produce IDs *above* the
maximum known `TrueStemID`. A stem first appearing at census 1 that was
severed by a resprout barrier can carry a `ReconstructedStemID` larger
than a stem first appearing at census 6 — opposite to the chronological
convention in `OriginalStemID`.

#### Additional minting sites detected during 3-pass code review

The five sites above are the ones documented in the user-facing pipeline.
A careful pass over the full engine reveals **four more sites** that
mint engine IDs the same way (`max + 1L` ascending):

1. **`dp_global_dp.R` ~line 706** — NA-R barrier split inside the DP
   fallback path (when the constrained DP cannot solve and routes to a
   simpler enumeration). Same semantics as minting site 3 but in a
   different code branch.
2. **`dp_global_dp.R` ~line 755** — live-R-boundary track sever in the
   same fallback path; severs tracks that continue past a live R-coded
   row.
3. **`dp_probabilistic_matching.R` ~line 124** — anchor IDs for the
   probabilistic matcher (used for `PROB_SPECIES` such as figs, palms,
   `bactma`, `oenoma`, and as DP fallback). Pads anchor rows with
   sequential IDs starting above the max known `TrueStemID`.
8b. **`dp_probabilistic_matching.R` ~line 184** — commits those anchor
   IDs to `ReconstructedStemID`; also `seq_len(.N)` IDs at line 124 for
   the single-census fast-path branch.
4. **`export_probabilistic_posteriors()` in
   `dp_probabilistic_matching.R`** — writes posterior path files for
   the probabilistic matcher in exactly the same format as
   `dp_global_dp.R`, and calls the same `apply_bb_invariants_to_samples()`
   per-sample bb minter.

All of these are caught automatically by the renumbering algorithm
because step 3 (`engine_ids ← setdiff(all_recon, known_ids)`) is a
black-box partition: any ID that is not a real database ID and not a
`given_orphan` StemID is treated as engine-minted regardless of which
branch produced it. The algorithm therefore remains correct in the
presence of the additional sites; no per-site code changes are needed.

---

### Columns that reference `ReconstructedStemID`

Any renumbering plan must update **all** of the following consistently:

| Column / artefact | Where | Notes |
|---|---|---|
| `ReconstructedStemID` | main output table | the primary output |
| `ReconstructedStemID_PreSweep` | main output table | audit snapshot; rename for consistency |
| `DP_PosteriorTop{k}ID` (k = 1…posterior_top_k) | main output table | per-observation top-K candidate IDs |
| `ReconstructedStemID` in posterior path files | `posteriors/tag_*_paths.feather/rds/csv` | per-sample full-path reconstructions; `recon` field encodes `ObsRowID:ReconstructedStemID` pairs; `path_sig` is a `-`-joined string of per-census IDs — see dedicated subsection below |

`DP_PosteriorTop{k}Prob`, `DP_PosteriorEntropy`,
`DP_PosteriorReconstructedProb`, `DP_PosteriorUnlinkedProb`, and
`DP_PosteriorBin` are **probability / bin values**, not IDs — they do not
need renaming.

`TrueStemID` is **never renamed** — it is the external BCI database
reference.

---

### Proposed implementation: post-engine renumber function

Rather than changing the five minting sites (which would require
carefully avoiding negative integers, collisions, and inter-site
dependencies), the cleanest approach is a single **post-engine
renumbering pass** applied once per tag as the last step of
`run_dp_one_group()`, after all post-engine helpers have run.

#### Timing within the per-chunk pipeline

```
run_dp_one_group()
  ├─ match_stems_dp_global_backward_marginals_batch()   ← minting sites 1-4
  │    └─ finalize_out() / inner TrueStemID sweep
  ├─ script-level backstop sweep
  └─ return out

in run_main_chunked(), per chunk:
  ├─ maybe_add_posterior_bins()          ← uses probabilities only; OK before rename
  ├─ apply_carried_terminal_backfill()   ← LOCF, no new IDs minted
  ├─ apply_orphan_stem_backfill()        ← uses StemID (database ID), not engine-minted
  ├─ apply_broken_below_invariants()     ← minting site 5
  └─ renumber_engine_minted_ids()        ← NEW: must run last
```

`renumber_engine_minted_ids()` runs after all helpers because
`apply_broken_below_invariants()` (minting site 5) adds IDs that must
also be renumbered.

#### Algorithm

```
renumber_engine_minted_ids(out, posterior_top_k):

  1. Determine known_ids — real BCI database IDs that must never be renamed.

     Step 1a. Start from TrueStemID values that have NOT been overridden
     by either apply_broken_below_invariants() or DP minting site 1
     (provisional anchor IDs):

       ENGINE_MINTED_INTO_TRUESTEMID <- c(
           "provisional_dp",                # DP minting site 1 — anchor IDs
           "bb_split", "bb_split_carry",    # bb pin-override (minting site 5)
           "bb_post_terminator_split",
           "bb_post_terminator_split_carry"
       )

       db_ids ← unique(out$TrueStemID[
         !is.na(out$TrueStemID) &
         !(out$ReconstructionMethod %in% ENGINE_MINTED_INTO_TRUESTEMID)
       ])

     Why exclude bb_* rows: apply_broken_below_invariants() overwrites
     TrueStemID with the newly minted ID on pin-override rows (to keep
     TrueStemID and ReconstructedStemID consistent). Those overwritten
     TrueStemID values are engine-minted, not real database IDs.

     Why exclude provisional_dp rows: minting site 1 in dp_global_dp.R
     (DP provisional anchor when TrueStemID is missing on anchor rows)
     writes the freshly minted prov_id into BOTH ReconstructedStemID AND
     TrueStemID, with ReconstructionMethod = "provisional_dp". Without
     this exclusion, the prov_id would be classified as a real DB ID via
     TrueStemID, would not be renumbered, and would remain at a large
     positive value while all other engine_ids in the same tag get
     renumbered below min(known_ids) — breaking the chronological
     monotonicity invariant. The probabilistic engine does NOT have this
     problem: when it pads anchor rows with fabricated IDs (line ~184)
     it writes only to ReconstructedStemID (not TrueStemID) and tags the
     row "probabilistic".

     Side-effect of excluding provisional_dp: after renumbering, the
     TrueStemID column on provisional_dp rows still holds the original
     engine-minted prov_id (large value), while ReconstructedStemID
     holds the renumbered value. To keep TrueStemID and
     ReconstructedStemID consistent on these rows (mirroring what
     apply_broken_below_invariants does for bb_* rows), step 8 of the
     algorithm must also apply the mapping to TrueStemID on rows where
     ReconstructionMethod ∈ ENGINE_MINTED_INTO_TRUESTEMID. See step 8
     update below.

     Step 1b. Also protect IDs assigned by apply_orphan_stem_backfill():

       orphan_ids ← unique(out$ReconstructedStemID[
         !is.na(out$ReconstructionMethod) &
         out$ReconstructionMethod == "given_orphan"
       ])

     Why: apply_orphan_stem_backfill() sets ReconstructedStemID = StemID
     for rows where TrueStemID = NA and DBH = NA. StemID is a real BCI
     database ID, not engine-minted. Since TrueStemID = NA for these rows,
     they would not appear in db_ids; without explicit protection they
     would be misclassified as engine-minted and renumbered.

       known_ids ← union(db_ids, orphan_ids)

  2. all_recon ← unique(out$ReconstructedStemID[!is.na(out$ReconstructedStemID)])

  3. engine_ids ← setdiff(all_recon, known_ids)
     # All IDs that are not real database IDs: DP slack tracks, provisional
     # anchor IDs elevated above max TrueStemID, NA-R barrier splits, segment-
     # split offsets, and bb-minted splits (minting sites 1-5).

  4. For each engine_id, compute first_census:
       first_census[id] ← min(out$CensusID[out$ReconstructedStemID == id])

  5. Sort engine_ids by first_census ascending.
     Tie-break by original ReconstructedStemID value ascending for determinism.
     Result: engine_ids_sorted[1] = track with earliest first appearance.

  6. Compute new ID base:
       base ← min(known_ids) - length(engine_ids)
     Assign new IDs:
       new_id[engine_ids_sorted[i]] ← base + (i - 1)
     So engine_ids_sorted[1] (earliest) gets the smallest ID (base),
     engine_ids_sorted[n] (latest) gets base + n - 1, which is still
     strictly below min(known_ids).

  7. Build mapping table:
       mapping ← data.table(Tag = tag, old_id = engine_ids,
                            new_id = new_id[engine_ids])

  8. Apply mapping to all ID columns in out:
       - ReconstructedStemID                       (always, all rows)
       - ReconstructedStemID_PreSweep              (always, all rows;
         stores the pre-sweep engine ID; may differ from
         ReconstructedStemID on rows where the TrueStemID sweep overrode
         the DP assignment; rename for auditing consistency)
       - DP_PosteriorTop{k}ID for k = 1 … posterior_top_k (always, all rows)
       - TrueStemID ONLY on rows where
         ReconstructionMethod ∈ ENGINE_MINTED_INTO_TRUESTEMID
         (keeps TrueStemID == ReconstructedStemID on provisional_dp and
         bb_* rows; other rows' TrueStemID is the real DB ID and must
         not be touched)
     (vectorised lookup: new_val <- mapping$new_id[match(col, mapping$old_id)];
     NA cells are left as NA)

  9. Handle posterior path files — see dedicated subsection below.

  10. Return:
       list(out = renamed_out, mapping = mapping)
```

#### Edge cases

| Case | Handling |
|---|---|
| No engine-minted IDs (all rows are `given` or `given_orphan`) | `engine_ids` is empty; mapping is empty; function is a no-op |
| `known_ids` is empty (tag has no real TrueStemID and no orphan backfill) | use `0L` as the base so new IDs are 1-based sequential negative integers; the `new_id < min(known_ids)` guarantee cannot be made, but the chronological ordering within the tag is preserved; log a warning |
| Two engine-minted tracks with the same `first_census` | tie-break by original `ReconstructedStemID` value ascending; ensures determinism across re-runs |
| `ReconstructedStemID_PreSweep` is NA for a row | leave as NA; do not apply mapping to NA cells |
| `DP_PosteriorTop{k}ID` is NA for a row | leave as NA; only rename non-NA values |
| Posterior path files are absent (`posterior_samples = 0`) | skip step 9 entirely |

---

### Posterior path files — detailed treatment

#### What the path files contain

When `posterior_samples > 0`, the engine writes a
`posteriors/tag_{Tag}_posterior_samples_{ts}_paths.feather` file (or
`.rds` / `.csv`) with one row per unique reconstruction path. Key
columns:

- `path_sig`: a `-`-joined string of `ReconstructedStemID` values in
  `CensusID` order (e.g. `"890123-890124-891000"`). Used to collapse
  identical reconstructions.
- `recon`: a `;`-joined string of `ObsRowID:ReconstructedStemID` pairs
  (e.g. `"42:890123;43:891000"`). Used to re-attach a path to the main
  table via `ObsRowID`.
- `path_prob`, `path_count`: importance weights and counts.

#### Engine routing and posterior filename uniqueness

Both the DP engine
(`match_stems_dp_global_backward_marginals_batch()` in `dp_global_dp.R`)
and the probabilistic engine (`export_probabilistic_posteriors()` in
`dp_probabilistic_matching.R`) write posterior path files with the
**identical filename pattern**
`tag_{Tag}_posterior_samples_{ts}_paths.{feather|rds|csv}`.

Despite the shared pattern, **only one engine writes per tag per run**:

- Tags whose species is in `PROB_SPECIES` (figs, palms, `bactma`,
  `oenoma`, etc.) are routed directly to the probabilistic engine, which
  writes the posterior file.
- Tags that hit `do_fallback()` inside the DP (state-space too large,
  anchor missing, K too small) are also routed to the probabilistic
  engine; the DP returns the prob engine's output without writing its
  own posterior file.
- All other tags are solved by the DP, which writes the posterior file
  directly.
- Segment-split sub-calls inside the DP explicitly disable posteriors
  (`posterior_samples = 0L`) so they never write a second file for the
  same tag.

Consequence for the renumbering plan: **the renumbering function must
be engine-agnostic** — it receives the merged `out` table without
needing to know which engine produced any row, and the test code must
load `tag_{Tag}_posterior_samples_*_paths.*` by glob pattern (not
assume a specific engine).

#### Why path files are harder to rename than the main table

The `ReconstructedStemID` values in path files come from **two
independent sources**:

1. **DP track IDs** (minting sites 1–2): the `track_ids` vector built
   inside `match_stems_dp_global_backward_marginals_batch()`. These are
   the same integers that appear in the `DP_PosteriorTop{k}ID` columns
   of the main table.

2. **Per-sample bb-invariant IDs**: `apply_bb_invariants_to_samples()`
   is called inside the engine, *before* path files are written, to
   enforce the same R1/R2 broken-below contract on each posterior sample.
   It mints fresh IDs locally per sample using `max(sample_rid) + 1L`.
   These per-sample IDs are **independent** of the MAP-level IDs minted
   by the outer `apply_broken_below_invariants()` call (minting site 5)
   that runs after the engine returns.

Consequence: **a simple `old_id → new_id` mapping table is insufficient
for path files**. The mapping covers the DP track IDs (source 1), but
the per-sample bb IDs (source 2) cannot be translated by it because they
were generated independently inside the engine and do not correspond to
any ID in the mapping.

#### Recommended architecture (correct, requires engine change)

Instead of writing path files inside the engine, return `samples_dt`
(the per-sample reconstruction data.table) as part of the engine's
return value, and write path files *after* the renumbering pass:

```
1. match_stems_dp_global_backward_marginals_batch() returns
     list(out = ..., samples_dt = ...)   # samples_dt instead of writing

2. renumber_engine_minted_ids(out, ...) → produces mapping

3. Apply mapping to samples_dt$ReconstructedStemID
     (translates DP track IDs; per-sample bb IDs become stale)

4. Re-run apply_bb_invariants_to_samples(samples_dt, renamed_out)
     (re-derives per-sample bb IDs using renamed track IDs as the base)

5. Recompute path_sig, path_counts, paths_summary

6. Write paths_summary to the posteriors directory
```

This approach keeps all ID columns — MAP table and path files — in the
same renumbered space, and avoids the independent-bb-ID inconsistency.

Changes required in `dp_global_dp.R`:

- Remove the path-file writing block at the end of
  `match_stems_dp_global_backward_marginals_batch()`.
- Instead, attach `samples_dt` (and `sampling_profile`) to the returned
  list as `attr(out, "samples_dt")` and `attr(out, "sampling_profile")`.
- In `1_main_cpp_chunk_bci.R` (and `scripts/main_cpp_chunk.R`), after
  calling `renumber_engine_minted_ids()`, retrieve `samples_dt`, apply
  steps 3–6 above, and write the path files.

#### Fallback (simpler, documented limitation)

If the architectural change is not desired, write a companion mapping
file alongside each path file:

```
post_dir/tag_{Tag}_id_mapping.feather
  columns: Tag, old_id, new_id
```

Downstream users must:

1. Apply `old_id → new_id` to the `recon` field for DP track IDs.
2. Note that per-sample bb IDs that appear in `recon` but are **absent**
   from the mapping are internal per-sample integers and are not directly
   comparable to any ID in the renamed MAP table. Downstream uncertainty
   quantification that joins paths to the MAP table via `ObsRowID` (not
   via `ReconstructedStemID`) is unaffected.

---

### Files to modify

| File | Change |
|---|---|
| `dp_global/R/dp_global_main.R` | Add `renumber_engine_minted_ids()` function in the shared helpers section (section 7 or new section 8) |
| `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R` | In `run_main_chunked()`, call `renumber_engine_minted_ids()` after `apply_broken_below_invariants()` and before writing chunk outputs |
| `dp_global/scripts/main_cpp_chunk.R` | Same call added to the non-BCI chunked runner for parity |
| `dp_global/R/dp_global_dp.R` | If the recommended architecture is adopted: remove path-file writing; return `samples_dt` as `attr(out, "samples_dt")`; otherwise no changes needed for the fallback approach |

---

### Invariants that must hold after renumbering

1. `ReconstructedStemID == TrueStemID` for every row where
   `ReconstructionMethod == "given"` — **unchanged by renaming**
   (known_ids are never touched; the row's `TrueStemID` is a real DB ID
   and stays in `db_ids`).

1b. `ReconstructedStemID == TrueStemID` for every row where
   `ReconstructionMethod ∈ ENGINE_MINTED_INTO_TRUESTEMID`
   (`provisional_dp`, `bb_split`, `bb_split_carry`,
   `bb_post_terminator_split`, `bb_post_terminator_split_carry`) —
   **maintained** because step 8 applies the same mapping to both
   columns on these rows.

2. `new_id < min(known_ids)` for all engine-minted IDs — guaranteed by
   the base computation in step 6. Note: `known_ids` is the union of
   real-DB `db_ids` and `orphan_ids`; engine-minted IDs that were
   written into `TrueStemID` (`provisional_dp`, `bb_*`) are excluded
   from `db_ids` and therefore correctly classified as engine-minted.

3. The chronological ordering property: if track A first appears at an
   earlier census than track B, then `new_id(A) < new_id(B)` —
   guaranteed by the sort in step 5.

4. `DP_PosteriorTop1ID` reflects the DP-assigned MAP track ID for each
   row before NA-R barrier splitting. After renumbering, the mapping is
   applied consistently to both `DP_PosteriorTop1ID` and
   `ReconstructedStemID`. Note: for rows that went through barrier
   splitting (minting site 3), `DP_PosteriorTop{k}ID` was not updated by
   the barrier split — this is a pre-existing condition, not introduced
   by the renumbering.

5. Posterior path files are consistent with the MAP table: either via
   the recommended architecture (path files written after renumbering,
   using re-derived per-sample bb IDs) or via the fallback (companion
   mapping file plus documented limitation on per-sample bb IDs).

---

### What does NOT change

- `TrueStemID` on rows where `ReconstructionMethod == "given"` — these
  hold real BCI database IDs and are never renamed. (TrueStemID *is*
  renamed on `provisional_dp` and `bb_*` rows, but those values were
  already engine-minted, not real DB IDs; see Invariant 1b.)
- `StemID` / `OriginalStemID` columns — not touched.
- Probability and entropy columns (`DP_PosteriorTop{k}Prob`,
  `DP_PosteriorEntropy`, `DP_PosteriorReconstructedProb`,
  `DP_PosteriorUnlinkedProb`, `DP_PosteriorBin`) — numerical values;
  no ID renaming needed.
- `ReconstructionMethod` labels — unchanged.
- Audit flags (`SweepAuditOverride`, `SweepRollbackToPreSweep`) — unchanged.
- The core DP algorithm, cost functions, and minting sites — no changes.

---

### Verification of the plan against the engine code (6-pass review)

The plan was checked against the full engine on six independent passes
(three initial + three additional). Findings:

1. **Pass 1 — minting sites.** The five originally enumerated sites plus
   the four additional ones documented above account for every
   assignment to `ReconstructedStemID` in `dp_global/R/dp_global_dp.R`,
   `dp_global/R/dp_global_main.R`, and
   `dp_global/R/dp_probabilistic_matching.R`. The grep audit found
   49 assignment expressions; all of them either (a) mint a fresh ID
   (covered by the partition rule), (b) copy from `TrueStemID` (rows are
   in `known_ids`), (c) copy from `StemID` via the orphan backfill
   (rows are explicitly added to `known_ids`), or (d) carry an existing
   `ReconstructedStemID` to another row (no new ID introduced).

2. **Pass 2 — column coverage.** The grep audit confirmed the only
   columns holding `ReconstructedStemID`-valued integers are:
   `ReconstructedStemID`, `ReconstructedStemID_PreSweep`,
   `DP_PosteriorTop{k}ID` (k = 1..DP_POSTERIOR_TOP_K), and the
   `ReconstructedStemID` column inside `samples_dt` / posterior path
   files. The plan renames all four.

3. **Pass 3 — ordering of helpers.** The plan places
   `renumber_engine_minted_ids()` strictly after
   `apply_broken_below_invariants()`. Verified against
   `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`
   (`run_main_chunked()` invokes them in the order:
   `maybe_add_posterior_bins` → `apply_carried_terminal_backfill` →
   `apply_orphan_stem_backfill` → `apply_broken_below_invariants`).
   The plan extends this chain by one step.

4. **Pass 3b — single-tag driver.** `dp_global/scripts/main_cpp_bci.R`
   runs the same chain of helpers (its section 9b). The renumbering
   call must be added there too, after `apply_broken_below_invariants()`
   and before the writers in section 10.

5. **Pass 4 — engine routing.** Verified the DP-vs-probabilistic
   routing in `match_stems_dp_global_backward_marginals_batch()`. A tag
   takes exactly one of three paths: (a) DP success, (b) direct route to
   `match_stems_probabilistic()` for `PROB_SPECIES`, or (c) DP failure
   followed by `do_fallback()` → `match_stems_probabilistic()`. In all
   three cases the renumbering function operates on the merged `out`
   table and is engine-agnostic. Segment-split sub-calls inside the DP
   recursively call the DP with `posterior_samples = 0L` (line 1168),
   so no posterior file is written by sub-segments.

6. **Pass 5 — posterior file naming.** Both engines write posteriors
   with the identical filename pattern
   `posteriors/tag_{Tag}_posterior_samples_{ts}_paths.{ext}`.
   `tag_local` / `tag_val` is identical between the two writers and the
   `ts_local` timestamp comes from the shared `BATCH_TS`. Because only
   one engine fires per tag per run, the two never collide. The test
   loader must use a glob (e.g.
   `Sys.glob(file.path(out_dir, "posteriors", paste0("tag_", tag, "_posterior_samples_*_paths.*")))`)
   rather than assuming a specific engine wrote the file.

7. **Pass 6 — segment-split posterior consistency.** When the DP splits
   a tag into sub-segments and the segment-split offset (minting site 4)
   shifts pre-segment IDs upward, the offset is also applied to
   `DP_PosteriorTop{k}ID` columns (`dp_global_dp.R` ~line 1267). The
   renumbering function therefore sees a self-consistent main table at
   the time it runs, regardless of whether the tag was split. The
   partition invariant (Stage 1 test) is the operational guarantee that
   the renumbering preserved this consistency.

No additional changes to the plan were required after the 6-pass
review. The black-box partition strategy was specifically chosen to be
robust against undiscovered minting sites; it has now been validated to
cover the four sites in `dp_global_dp.R` fallback paths, the prob-engine
sites in `dp_probabilistic_matching.R`, and the segment-split offset.

#### Three additional passes — segment splits with mixed engines

A further three passes were performed focused on the user's specific
scenario: tags where the engine splits the process (segment split at an
NA-R / R barrier, or bb-split at the MAP level) and the two halves may
be handled by different engines (DP vs probabilistic).

8. **Pass 7 — mixed-engine segment splits.** The segment-split branch
   at `dp_global_dp.R` ~line 1133 recursively calls
   `match_stems_dp_global_backward_marginals_batch()` on the pre- and
   post-resprout halves independently. Each sub-call is free to either
   solve via DP or hit `do_fallback()` and route to the probabilistic
   engine. Therefore, **for a single tag, the pre-segment and the
   post-segment can be solved by different engines**. The two halves
   are recombined via `rbindlist` (line ~1278) and returned through
   `finalize_out(combined)`. Verification:
   - The renumbering function operates on the merged `out` table per
     tag; the black-box partition `engine_ids = setdiff(all_recon,
     known_ids)` is engine-agnostic, so it catches engine-minted IDs
     from both halves regardless of which engine produced them.
   - The current segment-split offset (minting site 4) pushes
     pre-segment IDs *above* post-segment IDs numerically, but the
     pre-segment has *earlier* censuses by construction. After the
     renumber sort by `first_census` ascending, the pre-segment engine
     IDs receive the smallest renumbered values and the post-segment
     engine IDs receive larger values — chronologically correct, and
     also exactly reversing the current backward convention.
   - Pre-segment provisional anchor IDs (DP minting site 1, inside the
     pre-sub-call) are tagged `provisional_dp` and are correctly
     excluded from `db_ids` by the Pass 7 fix above.
   - Pre-segment prob-fallback anchor IDs are tagged `probabilistic`
     and never written to `TrueStemID`, so they are caught by the
     partition automatically.

9. **Pass 8 — segment-split posterior files.** Confirmed that the
   outer engine call returns at the segment-split path early
   (`return(finalize_out(combined))` at `dp_global_dp.R` line ~1281),
   which is *before* the posterior-writing block at line ~2456.
   Both sub-calls pass `posterior_samples = 0L`, so even if a sub-call
   falls back to the prob engine, the prob engine's writer (gated by
   `if (posterior_samples > 0L)`) is skipped. Net consequence: **no
   posterior path file is written for segment-split tags by any engine
   in any sub-call**. Step 9 of the renumbering algorithm (posterior
   path file handling) is a no-op for segment-split tags. This
   simplifies the recommended architecture: the `samples_dt` attribute
   will simply be `NULL` for segment-split tags, and the post-engine
   path-file writer can skip those.

10. **Pass 9 — bb-split at the MAP level on non-segment-split tags.**
    For a tag that runs as a single engine call (no segment split) and
    then has `apply_broken_below_invariants()` mint `bb_split` /
    `bb_post_terminator_split` IDs post-engine, the situation is:
    - One engine produced the entire `out` table for the tag (either
      DP or prob — never mixed, since there is no segment split).
    - The single posterior path file (if any) was written by that one
      engine and contains track IDs in the engine's original space
      plus per-sample bb IDs minted by `apply_bb_invariants_to_samples`
      *inside* the engine.
    - `apply_broken_below_invariants()` then mints fresh `bb_*` IDs at
      the MAP level, independently of the per-sample IDs.
    - The renumbering function sees one merged table with one set of
      engine_ids (DP/prob tracks + MAP-level bb IDs). The partition
      catches all of them. The recommended-architecture posterior
      handling (re-run `apply_bb_invariants_to_samples` after track-ID
      renaming) re-derives consistent per-sample bb IDs.
    Conclusion: the user's stated concern — that bb-splitting could
    cause one half of a tag to be handled by DP and another by prob —
    can only happen via the segment-split mechanism (covered by Pass 7
    above), not via `apply_broken_below_invariants()`, which runs on
    the already-merged table.

Net result of the three additional passes: **one bug fix needed** —
adding `"provisional_dp"` to the exclusion list in step 1a, with a
corresponding extension to step 8 to apply the mapping to `TrueStemID`
on engine-minted-into-TrueStemID rows. Both have been incorporated
into the Algorithm section above. No other plan changes are required.
The mixed-engine segment-split scenario is fully covered.

---

### Test plan

The renumbering will be tested via a **500-tag before/after regression**
on the hardest cases. The test driver is the single-tag debug script
`dp_global/scripts/main_cpp_bci.R`, which loads
`bci_data/multistem_tags.rds` and runs the full pipeline for one tag.

#### Inputs

- **Data**: `bci_data/multistem_tags.rds` (the canonical multistem
  dataset loaded by `main_cpp_bci.R` via `INPUT_FILE`).
- **Tags**: a stratified sample of **500 tags** drawn from
  `bci_data/check_sweep_audit_override_tags.csv` (~32 349 zero-padded
  tag strings). These are tags that triggered
  `SweepAuditOverride = TRUE` in past runs, meaning the engine produced
  a `ReconstructedStemID` that the script-level sweep had to override
  against `TrueStemID`. They concentrate the hardest cases:
  segment-split tags, NA-R barrier crossings, BB pin overrides, and
  fallback-path tags. If the renumbering survives this set, it will
  survive the easier tags.

#### Stratification of the 500-tag sample

Sample with a fixed seed so the same 500 tags are used before and
after the code change. Suggested stratification (proportional to
frequency in the source list, with a floor per stratum):

| Stratum | Heuristic | Target count |
|---|---|---|
| `PROB_SPECIES` tags (figs / palms / `bactma` / `oenoma`) | species lookup | 50 |
| DP fallback (`PROB_LOOKAHEAD_WEIGHT` indicates `DP_FALLBACK_GROWTH_FORMS`) | growth-form lookup | 50 |
| Segment-split tags (large `n_obs` per tag at C7) | top of `n_obs` sorted desc | 50 |
| BB pin-override tags (any `bb_split*` method in baseline) | from baseline `ReconstructionMethod` | 100 |
| NA-R barrier-crossing tags | tags with R-coded NA-DBH rows at any census | 100 |
| Random remainder | uniform sample | 150 |

Save the resulting list to `bci_data/test_renumber_500_tags.csv` so
both the before and after runs use exactly the same tags.

#### Phase A — Baseline run (BEFORE code changes)

On the current `main` branch (or the parent commit of the renumbering
branch):

```bash
mkdir -p logs out/baseline
while read tag; do
    Rscript dp_global/scripts/main_cpp_bci.R \
        --WHICH_TAG=$tag \
        --POSTERIOR_SAMPLES=50 \
        --POSTERIOR_SAMPLE_SEED=42 \
        --POSTERIOR_SAMPLES_FORMAT=feather \
        >> logs/baseline_run.log 2>&1
done < <(tail -n +2 bci_data/test_renumber_500_tags.csv)
```

For each tag, archive the run's `out_dir` contents to
`out/baseline/tag_<tag>/`:

- `stem_reconstruction_dp_global_rcpp.rds` — the main output table.
- `posteriors/tag_<tag>_posterior_samples_*_paths.feather` — the
  posterior path file (may have been written by either the DP or the
  probabilistic engine; load by glob).
- The `run_started.txt` / `run_finished.txt` markers (timestamps).

#### Phase B — Renumbered run (AFTER code changes)

On the feature branch with `renumber_engine_minted_ids()` added:

```bash
git checkout rename_stemids
mkdir -p out/renumbered
while read tag; do
    Rscript dp_global/scripts/main_cpp_bci.R \
        --WHICH_TAG=$tag \
        --POSTERIOR_SAMPLES=50 \
        --POSTERIOR_SAMPLE_SEED=42 \
        --POSTERIOR_SAMPLES_FORMAT=feather \
        >> logs/renumbered_run.log 2>&1
done < <(tail -n +2 bci_data/test_renumber_500_tags.csv)
```

Archive each tag's `out_dir` to `out/renumbered/tag_<tag>/` with the
same structure as Phase A.

Using the **same** `POSTERIOR_SAMPLE_SEED = 42` ensures the upstream
sampler draws the same paths in both runs; the only differences in the
output should be the renumbered ID values.

#### Phase C — Regression assertions

Write the assertion script as
`dp_global/tests/test_renumber_regression.R`. For each tag in the
500-tag list, load both the baseline and renumbered RDS files and
apply:

1. **Schema preservation**
   - Identical row count, identical `ObsRowID` set.
   - Identical column names (the renumbering must not add or drop
     columns beyond, optionally, the mapping table).

2. **Untouched columns**
   - `TrueStemID` identical between baseline and renumbered **except**
     on rows where `ReconstructionMethod ∈ ENGINE_MINTED_INTO_TRUESTEMID`
     (`provisional_dp`, `bb_split`, `bb_split_carry`,
     `bb_post_terminator_split`, `bb_post_terminator_split_carry`).
     On those rows, verify instead that
     `renumbered$TrueStemID == renumbered$ReconstructedStemID` (Invariant 1b).
   - `ReconstructionMethod` identical.
   - `StemID`, `OriginalStemID`, `DBH`, `Status`, `ListOfTSM`, `Tag`,
     `CensusID`, `Species` identical.

3. **Partition invariance** (the core operational guarantee)
   - For all pairs `(i, j)`:
     `baseline$ReconstructedStemID[i] == baseline$ReconstructedStemID[j]`
     iff `renumbered$ReconstructedStemID[i] == renumbered$ReconstructedStemID[j]`.
   - Implemented efficiently as: build a contingency table
     `table(baseline$ReconstructedStemID, renumbered$ReconstructedStemID)`;
     the table must have exactly one non-zero entry per row and per
     column.

4. **Chronological monotonicity for engine IDs**
   - Build the renumbered table's `(ReconstructedStemID, first_census)`
     pairs; restrict to IDs that are not in `TrueStemID` and not
     `given_orphan`; verify `cor(first_census, ReconstructedStemID,
     method = "spearman") == 1` within each tag.

5. **Engine-minted IDs below all real database IDs**
   - For tags where `known_ids` is non-empty:
     `max(renumbered$ReconstructedStemID[ReconstructionMethod ∈
     engine_methods]) < min(known_ids)`.

6. **Posterior probabilities and entropies unchanged**
   - `all.equal(baseline$DP_PosteriorTop{k}Prob,
     renumbered$DP_PosteriorTop{k}Prob)` for k = 1..DP_POSTERIOR_TOP_K.
   - Same for `DP_PosteriorEntropy`, `DP_PosteriorReconstructedProb`,
     `DP_PosteriorUnlinkedProb`, and `DP_PosteriorBin`.

7. **Posterior top-K IDs partition-equivalent**
   - `DP_PosteriorTopkID` must satisfy the same partition invariance
     property as `ReconstructedStemID` (Assertion 3) within each k.

8. **Sweep audit unchanged**
   - `sum(baseline$SweepAuditOverride, na.rm = TRUE) ==
     sum(renumbered$SweepAuditOverride, na.rm = TRUE)`.
   - `identical(baseline$ReconstructedStemID_PreSweep ≠
     baseline$ReconstructedStemID, renumbered$ReconstructedStemID_PreSweep
     ≠ renumbered$ReconstructedStemID)` — the set of rows where the
     sweep overrode the engine must be identical.

9. **Posterior path file**
   - Locate the file by glob (works for both engines):
     `path <- Sys.glob(file.path(out_dir, "posteriors",
     paste0("tag_", tag, "_posterior_samples_*_paths.feather")))`.
   - Read baseline and renumbered `paths_summary`.
   - `path_prob` and `path_count` columns must match.
   - The `recon` column (a `ObsRowID:ReconstructedStemID;...` string)
     must satisfy partition invariance: parse both into
     `data.table(ObsRowID, RecID)`, then verify that the equivalence
     classes defined by RecID match between baseline and renumbered.

10. **Determinism**
    - Re-run Phase B once more with the same seed; the renumbered
      output for each tag must be byte-identical to the first Phase B
      run.

#### Phase D — Aggregate report

Write results to `dp_global/tests/test_renumber_500_results.csv` with
columns:

```
tag, engine_used, n_obs, n_engine_ids, n_known_ids,
assertion_1_schema, assertion_2_untouched,
assertion_3_partition_main, assertion_4_chronological,
assertion_5_below_known, assertion_6_probs,
assertion_7_partition_topk, assertion_8_audit,
assertion_9_posterior_paths, assertion_10_determinism,
overall_pass
```

`engine_used` is one of `dp`, `prob_species`, `dp_fallback_to_prob`,
inferred from the log files or from `DP_FallbackReason` in the output
table.

#### Exit criteria

- `overall_pass == TRUE` for **all 500 tags**.
- Any failing tag is added to a permanent regression-test fixture under
  `dp_global/tests/fixtures/<tag>/{baseline,renumbered}/` so the
  failure can be reproduced and the unit test added to
  `dp_global/tests/test_renumber_regression.R`.
- Logs are inspected for any unexpected warnings, especially the
  edge-case warning emitted when `known_ids` is empty.
