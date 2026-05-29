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

6. **`dp_global_dp.R` ~line 706** — NA-R barrier split inside the DP
   fallback path (when the constrained DP cannot solve and routes to a
   simpler enumeration). Same semantics as minting site 3 but in a
   different code branch.
7. **`dp_global_dp.R` ~line 755** — live-R-boundary track sever in the
   same fallback path; severs tracks that continue past a live R-coded
   row.
8. **`dp_probabilistic_matching.R` ~line 124** — anchor IDs for the
   probabilistic matcher (used for `PROB_SPECIES` such as figs, palms,
   `bactma`, `oenoma`, and as DP fallback). Pads anchor rows with
   sequential IDs starting above the max known `TrueStemID`.
8b. **`dp_probabilistic_matching.R` ~line 184** — commits those anchor
   IDs to `ReconstructedStemID`; also `seq_len(.N)` IDs at line 124 for
   the single-census fast-path branch.
9. **`export_probabilistic_posteriors()` in
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
     by apply_broken_below_invariants():

       db_ids ← unique(out$TrueStemID[
         !is.na(out$TrueStemID) &
         !(out$ReconstructionMethod %in% c(
             "bb_split", "bb_split_carry",
             "bb_post_terminator_split", "bb_post_terminator_split_carry"
         ))
       ])

     Why exclude bb_* rows: apply_broken_below_invariants() overwrites
     TrueStemID with the newly minted ID on pin-override rows (to keep
     TrueStemID and ReconstructedStemID consistent). Those overwritten
     TrueStemID values are engine-minted, not real database IDs.

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
       - ReconstructedStemID
       - ReconstructedStemID_PreSweep  (stores the pre-sweep engine ID;
         may differ from ReconstructedStemID on rows where the TrueStemID
         sweep overrode the DP assignment; rename for auditing consistency)
       - DP_PosteriorTop{k}ID for k = 1 … posterior_top_k
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
`posteriors/tag_{Tag}_paths.feather` file (or `.rds` / `.csv`) with one
row per unique reconstruction path. Key columns:

- `path_sig`: a `-`-joined string of `ReconstructedStemID` values in
  `CensusID` order (e.g. `"890123-890124-891000"`). Used to collapse
  identical reconstructions.
- `recon`: a `;`-joined string of `ObsRowID:ReconstructedStemID` pairs
  (e.g. `"42:890123;43:891000"`). Used to re-attach a path to the main
  table via `ObsRowID`.
- `path_prob`, `path_count`: importance weights and counts.

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

1. `ReconstructedStemID == TrueStemID` for every row where `TrueStemID`
   is non-NA and `ReconstructionMethod == "given"` — **unchanged by
   renaming** (known_ids are never touched).

2. `new_id < min(TrueStemID)` for all engine-minted IDs — guaranteed by
   the base computation in step 6.

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

- `TrueStemID` column — never renamed.
- `StemID` / `OriginalStemID` columns — not touched.
- Probability and entropy columns (`DP_PosteriorTop{k}Prob`,
  `DP_PosteriorEntropy`, `DP_PosteriorReconstructedProb`,
  `DP_PosteriorUnlinkedProb`, `DP_PosteriorBin`) — numerical values;
  no ID renaming needed.
- `ReconstructionMethod` labels — unchanged.
- Audit flags (`SweepAuditOverride`, `SweepRollbackToPreSweep`) — unchanged.
- The core DP algorithm, cost functions, and minting sites — no changes.

---

### Verification of the plan against the engine code (3-pass review)

The plan was checked against the full engine on three independent passes.
Findings:

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

No additional changes to the plan were required after the 3-pass review.
The black-box partition strategy was specifically chosen to be robust
against undiscovered minting sites; it has now been validated to also
cover the four sites in `dp_global_dp.R` fallback paths and in
`dp_probabilistic_matching.R`.

---

### Test plan

The renumbering will be tested in two stages, using the single-tag
debug driver `dp_global/scripts/main_cpp_bci.R` and the list of
problematic tags in `bci_data/check_sweep_audit_override_tags.csv`.

#### Inputs

- **Data**: `bci_data/multistem_tags.rds` (the canonical multistem
  dataset loaded by `main_cpp_bci.R` via `INPUT_FILE`).
- **Tags**: `bci_data/check_sweep_audit_override_tags.csv` (single column
  `check_tags`, ~32 349 zero-padded tag strings). These are tags that
  triggered `SweepAuditOverride = TRUE` in past runs, meaning the engine
  produced a `ReconstructedStemID` that the script-level sweep had to
  override against `TrueStemID`. They concentrate the hardest cases:
  segment-split tags, NA-R barrier crossings, BB pin overrides, and
  fallback-path tags. If the renumbering survives this set, it will
  survive the easier tags.

#### Stage 1 — Per-tag invariant tests on a representative subset

Goal: confirm the renumbering preserves all behavioural invariants on a
tractable sample.

1. Sample ~100 tags from `check_sweep_audit_override_tags.csv` stratified
   by complexity (a few one-stem tags, several multi-stem, plus 5-10
   tags from each of: BB pin override, segment-split, NA-R barrier,
   probabilistic fallback). Save the sampled tag list to
   `bci_data/test_renumber_tags.csv`.

2. For each tag, run the driver twice with identical settings:

   ```
   # Baseline (current main):
   git checkout main
   Rscript dp_global/scripts/main_cpp_bci.R \
       --WHICH_TAG=<tag> \
       --POSTERIOR_SAMPLES=50 \
       --POSTERIOR_SAMPLE_SEED=42

   # Renumbered (feature branch):
   git checkout rename_stemids
   Rscript dp_global/scripts/main_cpp_bci.R \
       --WHICH_TAG=<tag> \
       --POSTERIOR_SAMPLES=50 \
       --POSTERIOR_SAMPLE_SEED=42
   ```

3. For each tag, load both output RDS files
   (`stem_reconstruction_dp_global_rcpp.rds` in the run's `out_dir`) and
   apply the following assertions in a test script (e.g.
   `dp_global/tests/test_renumber_invariants.R`):

   - **Equal row count and `ObsRowID` set** between baseline and
     renumbered.
   - **Identical `TrueStemID`** column (renumbering must not touch real
     database IDs).
   - **Identical `ReconstructionMethod`** column.
   - **Same partition of rows by reconstructed track**: for every pair
     of rows `(i, j)` with `ObsRowID_i`, `ObsRowID_j`,
     `baseline$ReconstructedStemID[i] == baseline$ReconstructedStemID[j]`
     iff `renumbered$ReconstructedStemID[i] == renumbered$ReconstructedStemID[j]`.
     This is the operationally important invariant: the renumbering must
     be a pure relabelling, not a re-assignment.
   - **Chronological monotonicity for engine IDs**: build the renumbered
     table's `(ReconstructedStemID, first_census)` pairs; restrict to
     IDs that are not in `TrueStemID` and not `given_orphan`; verify
     `cor(first_census, ReconstructedStemID, method = "spearman") == 1`
     within each tag.
   - **`ReconstructedStemID < min(TrueStemID)`** for all engine-minted
     IDs in tags that have at least one real `TrueStemID`.
   - **Posterior consistency**: for each row,
     `renumbered$DP_PosteriorTop1ID == renumbered$ReconstructedStemID`
     iff the same was true in the baseline.
   - **Posterior probabilities unchanged**:
     `all.equal(baseline$DP_PosteriorTop{k}Prob, renumbered$DP_PosteriorTop{k}Prob)`
     for k = 1..DP_POSTERIOR_TOP_K; same for `DP_PosteriorEntropy`,
     `DP_PosteriorReconstructedProb`, `DP_PosteriorUnlinkedProb`.
   - **Posterior path file**: load
     `out_dir/posteriors/tag_<tag>_posterior_samples_<ts>_paths.feather`.
     Apply the renumbering mapping to the `recon` field; verify the
     translated path probabilities and partition match the baseline.
   - **Sweep audit unchanged**:
     `sum(baseline$SweepAuditOverride) == sum(renumbered$SweepAuditOverride)`.

4. Aggregate per-tag pass/fail into a CSV
   (`dp_global/tests/test_renumber_results.csv`) with columns
   `tag, n_engine_ids, n_known_ids, all_invariants_passed,
   failing_invariant`. Investigate any tag where
   `all_invariants_passed == FALSE`.

#### Stage 2 — Full-set regression with posterior sampling on a hard subset

Goal: confirm the renumbering is stable across the full problematic
tag list and across reruns with different random seeds.

1. Run the driver across **all** ~32 k tags in
   `check_sweep_audit_override_tags.csv` with posteriors disabled
   (`POSTERIOR_SAMPLES = 0`) to keep runtime tractable:

   ```
   while read tag; do
       Rscript dp_global/scripts/main_cpp_bci.R \
           --WHICH_TAG=$tag \
           --POSTERIOR_SAMPLES=0 \
           >> logs/renumber_stage2.log 2>&1
   done < <(tail -n +2 bci_data/check_sweep_audit_override_tags.csv)
   ```

   (Or use `RUN_ALL_TAGS=TRUE` once the driver supports an explicit tag
   list — currently it defaults to a single tag.)

2. For each tag, run the same invariant assertions as Stage 1 except
   the posterior-path file check (skipped because
   `POSTERIOR_SAMPLES = 0`).

3. Sample 200 of the most complex tags (largest `n_engine_ids`) and
   re-run with `POSTERIOR_SAMPLES = 200` and two different seeds
   (`POSTERIOR_SAMPLE_SEED = 42` and `= 7`). For each, verify:

   - The MAP `ReconstructedStemID` is identical across seeds (renumbering
     is deterministic given the engine output).
   - The posterior path partition is identical across seeds (different
     samples may draw different paths but the path *labels* should be
     stable under the renumbering applied to the same engine output).

#### Stage 3 — Comparison against the Dryad reference table

Goal: confirm the renumbered output still matches the Dryad
reconstruction at the same rate as the baseline.

1. Run the existing comparison script
   `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/3_initial_comparisson_dryad_and_me.R`
   on the renumbered output. The comparison is partition-based (it
   compares which observations share a reconstructed identity, not the
   integer values), so the agreement rate with Dryad must be **identical**
   between baseline and renumbered.

2. If the agreement rate changes by more than ±0.01 percentage points,
   the renumbering has changed the partition (a bug) — investigate.

#### Exit criteria

All three stages must pass with zero invariant violations before the
renumbering is merged. Any tag in Stage 2 that fails an invariant must
be added to a permanent regression-test fixture under
`dp_global/tests/fixtures/` and a unit test added in
`dp_global/tests/test_renumber_invariants.R`.
