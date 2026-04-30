# Fix: R/OR Code Handling Is Biologically Inverted

**File:** `dp_global/R/dp_global_dp.R`
**Status:** Known bug, not yet fixed
**Priority:** High — Bug 4 (unlabelled rows) is the dominant observable failure and is fixable post-hoc; Bugs 1/2/5 are deeper logic errors requiring DP re-run

---

## Biological Meaning of the R Code

### What the field crew actually records

The R family of codes describes the **current physical state** of a stem at the time of census measurement. They are recorded by field workers when the bole has broken below the standard breast-height measurement point (1.3 m). They are NOT death codes — they drive `DFstatus = "broken below"`, not `"dead"` or `"stem dead"`. The tree individual may still be alive through basal resprouting.

| Code | Full name | Precise meaning | DFstatus |
|------|-----------|-----------------|----------|
| `R`  | Broken | The bole has broken below 1.3 m. Primary trigger code. DBH may be NA (stump too low to measure) or >0 (stump remnant still measurable above break point). | `broken below` |
| `RP` | Broken — resprout | Bole broken below 1.3 m **and** a new vegetative resprout is already visible and recorded at this census. The resprout itself receives a new StemID. | `broken below` |
| `RF` | Broken — fallen | Broken bole has also fallen to the ground. | `broken below` |
| `RT` | Broken — at top | Broken at or above the measurement point, but not high enough for a standard DBH re-measurement. | `broken below` |
| `OR` | Other breakage | Structural damage not captured by R, RP, RF, or RT (e.g. partial splits, hinge breaks). | `broken below` |
| `QR` | *(unverified)* | Present in the code regex throughout `dp_global_dp.R` but **not found in the BCI ForestGEO census data** or in the ForestGEO code reference (`forestgeo_codes_reference.qmd`). May be a legacy or site-specific code. Should be confirmed before the fix is applied. | — |

### What the row represents for identity reconstruction

The R-coded row records the state of the **old, pre-break bole** at the census where the break was first observed. It is the last time that physical bole was measurable. Crucially:

- The R-coded row **IS the old stem** — it must carry the same `ReconstructedStemID` as all prior censuses where that stem was alive.
- Any new bole(s) growing from the base or stump (often coded `RP` or appearing as a new StemID in subsequent censuses) are **genuinely new organisms** with new IDs.
- The break event happened **between** the preceding census and this one. The R row marks the first observation of the post-break state, not the last observation of an intact bole.

DBH on an R-coded row varies by field conditions:
- `DBH = NA` — stump is too low (below 1.3 m) for any measurement. Common for `OR` stump chains.
- `DBH > 0` — a remnant portion of the bole remains above 1.3 m and was still measured, even though it is no longer a growing bole.

| Census | Row    | Code  | Meaning |
|--------|--------|-------|--------|
| C1     | Stem A | —     | Alive, growing normally |
| C2     | Stem A | **R** | Same bole, now broken — **last record of Stem A's physical bole** |
| C2     | Stem B | —     | New basal resprout — **first record of Stem B** |
| C3     | Stem B | —     | Growing |

The R-coded row at C2 **IS Stem A** — it must carry the same `ReconstructedStemID` as C1.
Stem B at C2 is genuinely new and needs a fresh ID.

---

## Observed Manifestations in Real Data

Analysis of ~70 tags with R/OR codes reveals three structural patterns. Bug 4 dominates.

### Pattern A — OR stump-chain (majority, ~45 tags)

Examples: tags 216273, 236424, 252448, 400337, 600060 and many others.

```
Old StemID:  C1–Ck   alive, real DBH  →  RStemID = X           ✓
Stump StemID C(k+1)  R,  DBH=NA       →  RStemID = NA          ✗ Bug 4
Stump StemID C(k+2)  OR, DBH=NA       →  RStemID = NA          ✗ Bug 4 (pre-anchor)
Stump StemID C(a+1)  OR, DBH=NA       →  RStemID = StemID      ✓ post-anchor path
```

The stump StemID represents a dead/inactive physical base with no measurable DBH.
The biological assignment is unambiguous — it should carry its own StemID as RStemID
for all rows, pre- and post-anchor alike. The intra-stem inconsistency (NA pre-anchor,
StemID post-anchor for the same stump StemID) is a direct consequence of Bug 4.

**Pre-break old stem is correctly handled.** The algorithm correctly stops the old
stem's RStemID at its last live census and does not extend it to the stump rows.

**Bugs 1 and 2 do NOT visibly fire** in this pattern because OR-only stumps never
have a live DBH measurement in the post-break period — the DP has nothing to assign
the stump a track to, so the R-recruit constraint and segment-split boundary are
never exercised against a real DBH row.

---

### Pattern B — NA-R gap with recovery (~15 tags)

Examples: tags 017702, 018279, 050727, 116578, 220679, 282902.

```
StemID:  C1–Ck    alive, real DBH  →  RStemID = X   ✓
StemID:  C(k+1)   RP/R, DBH=NA    →  RStemID = NA  ✗ Bug 4 (gap row)
StemID:  C(k+2)   missing/NA      →  RStemID = NA  ✗ Bug 4 (gap row)
StemID:  C(k+3)+  alive, real DBH  →  RStemID = X   ✓ correctly re-linked
```

The same StemID persists across the break; the DP correctly re-links the track once
measurements resume. Only the gap rows (NA DBH, R/RP/missing) lack labels. No
identity mis-assignment occurs — only missing method labels on intermediate rows.

---

### Pattern C — Live R code with real DBH, multi-stem tags (~10 tags)

Examples: tags 005204, 027029, 106471, 412392.

These are the cases where Bug 1 (R-recruit constraint inverted) and Bug 2 (segment
boundary) could fire. However, in the observed data all such cases have `TrueStemID`
set on the recovery censuses (C4+ with real DBH), which allows the DP to pin the
correct track via `pin_truestemid`. The pin overrides the buggy constraint.
The observable failure is again limited to Bug 4: the NA-DBH gap rows between the
break census and the first recovery measurement have `RStemID=NA`.

Bug 1 would cause incorrect identity assignment specifically when:
- A stem has a live R-coded row WITH real DBH (not NA)
- No `TrueStemID` pin is available for that row
- The DP must infer identity across the break from DBH continuity alone

This scenario exists in the data (e.g., tags 027029 row 4 `R, DBH=11`;
116578 row 4 `R, DBH=13`; 282902 row 3 `R, DBH=13`) but is handled correctly
because those rows have `TrueStemID` pinning. The bug would surface in datasets
where `TrueStemID` is unavailable or not yet assigned.

---

## Can These Be Fixed Post-Hoc (After the DP Has Already Run)?

### Short answer

| Bug | Post-hoc fixable? | Confidence |
|-----|-------------------|------------|
| Bug 4 — unlabelled NA rows | **Yes — deterministically** | High |
| OR-chain intra-stem inconsistency | **Yes — deterministically** | High |
| Bug 1 — R-recruit inverted | Partially — only when fallback is detectable | Medium |
| Bug 2 — segment boundary | **No** — requires DP re-run | — |
| Bug 3 — K inflation | N/A — computation only, no output effect | — |
| Bug 5 — pin vs R-recruit conflict | **No** — requires DP re-run | — |

---

### Post-hoc fix for Bug 4 and OR-chain inconsistency

Both problems are fixable with a deterministic post-processing pass on the output
table, using the following rules applied **per StemTag group within each Tag**:

**Rule 1 — Fill RStemID for pre-anchor NA-DBH R/OR rows from the same StemTag chain.**

For every row where `DBH = NA` and `ListOfTSM` contains an R/OR code:
- Find other rows with the same `(Tag, StemTag)` that have a non-NA `RStemID`.
- If a non-NA `RStemID` exists anywhere in the chain (pre- or post-anchor), assign
  that value to all NA-RStemID rows in the chain.
- This resolves both Bug 4 (no label) and the OR-chain inconsistency (NA pre-anchor
  vs StemID post-anchor for the same stump).

**Rule 2 — Stamp `ReconstructionMethod` on newly filled rows.**

Rows updated by Rule 1 should receive `ReconstructionMethod = "dp"` to distinguish
them from rows assigned by the DP's normal track assignment.

**Rule 3 — Handle chains with no recoverable RStemID.**

If the entire StemTag chain has only NA-DBH rows and no post-anchor RStemID is
available (e.g., tag 241413 where the stump StemID has no real DBH anywhere):
- Assign `RStemID = StemID` (the stump's own ForestGEO StemID) as the identity.
- Set `ReconstructionMethod = "broken_stump_no_history"`.
- This is biologically correct: the stump IS the identity, even if it can't be
  measured.

**Sketch of the R post-processing function:**

```r
fix_or_stump_labels <- function(dt, resprout_regex) {
  # dt: output data.table with Tag, StemTag, CensusID, DBH,
  #     ListOfTSM, RStemID (= ReconstructedStemID), ReconstructionMethod

  dt[, RStemID := as.character(RStemID)]

  # Rows that are R/OR-coded with no DBH and no method label
  is_r_na <- is.na(dt$DBH) &
             grepl(resprout_regex, dt$ListOfTSM, perl = TRUE) &
             is.na(dt$ReconstructionMethod)

  # For each (Tag, StemTag) group, propagate the non-NA RStemID
  dt[, .fill_group := {
    known_id <- unique(RStemID[!is.na(RStemID)])
    if (length(known_id) == 1L) {
      # Exactly one identity in this StemTag — fill all NA-R rows
      RStemID[is_r_na & .I_group] <- known_id
      ReconstructionMethod[is_r_na & .I_group & is.na(ReconstructionMethod)] <- "dp"
    } else if (length(known_id) == 0L) {
      # No known identity — assign StemID itself
      RStemID[is_r_na & .I_group] <- StemID[is_r_na & .I_group]
      ReconstructionMethod[is_r_na & .I_group & is.na(ReconstructionMethod)] <-
        "broken_stump_no_history"
    }
    NULL  # by-reference modification
  }, by = .(Tag, StemTag)]

  dt
}
```

> **Note:** This sketch uses `.I_group` as a placeholder for within-group row indexing;
adjust to proper `data.table` by-group idiom (e.g., `.SD` rows).

---

### Post-hoc fix for Bug 1 (partial)

Bug 1 causes the DP to fall back to probabilistic for groups where the R-recruit
constraint makes all transitions infeasible. You can detect affected groups:

```r
# Groups where a live-R row (DBH != NA, R/OR code) exists and the group used
# probabilistic fallback — these may have incorrect identity assignment
affected <- dt[
  !is.na(DBH) &
  grepl(resprout_regex, ListOfTSM, perl = TRUE) &
  grepl("probabilistic", ReconstructionMethod)
, unique(Tag)]
```

For these groups, the only reliable fix is to re-run the DP after correcting Bug 1.
Manual inspection of `RStemID` assignments for the affected stems may reveal
single-stem tags where identity is trivially correct regardless, but for multi-stem
tags no post-hoc rule can substitute for the correct DP transition logic.

---

### Summary: what to do now vs. later

**Do now (post-hoc, no re-run needed):**
1. Apply the `fix_or_stump_labels()` pass above to fill all NA-RStemID R/OR rows
   and resolve the pre/post-anchor inconsistency.
2. Flag tags with live-R rows and probabilistic fallback for manual review.

**Do later (requires DP code fix + re-run):**
1. Fix Bug 1 to eliminate incorrect fallbacks on live-R rows without TrueStemID.
2. Fix Bug 2 to ensure live-R census goes into the correct DP segment.
3. Fix Bug 5 to resolve the pin-conflict fallback.

---

## What the Code Currently Does (Wrong)

### Bug 1 — R-recruit constraint is inverted

**Location:** `dp_global_dp.R` ~line 1807, backward-pass transition loop

The current code treats the R-coded observation as a **new recruit**:

```r
# Current comment (incorrect):
# "An observation with a resprout code ... is a NEW organism going backward
#  in time — the track it occupies at p+1 MUST BE EMPTY at p"
```

This forces R-coded rows to start from empty tracks — the opposite of the biological
truth. The R-coded row should **continue** a pre-existing track (the old stem's identity).
The non-R sibling stems at the same census are the genuine new organisms; standard recruit
cost logic already handles those correctly.

**Consequence:** All transitions where the R-coded stem continues an existing track are
pruned. For tags where this is the only biologically valid assignment, the DP either:
- Finds no feasible states → `no_feasible_edges` fallback to probabilistic, or
- Produces an incorrect pairing (old stem matched to wrong identity).

---

### Bug 2 — Segment split puts the R census in the wrong segment

**Location:** `dp_global_dp.R` ~line 1082, resprout segment split

```r
# Current (wrong):
.post_data <- tree_data[CensusID >= .r_boundary_census]  # R census = start of post-segment
.pre_data  <- tree_data[CensusID <  .r_boundary_census]  # R census absent from pre-segment
```

The R census is the **last census of the old stem** and belongs in the pre-segment.
The post-segment should start at the census **after** R, where only genuinely new boles exist.

**Consequence:** The R-coded row (old stem's final measurement) is fed into
`post_segment_all_recruits = TRUE`, treating the old stem's last observation as a
brand-new recruit. The old stem's census history is reconstructed in two disconnected
sub-problems with no identity continuity across the break event.

---

### Bug 3 — K inflation overcounts needed tracks

**Location:** `dp_global_dp.R` ~line 1172

```r
n_resprout_extra <- sum(...)  # one extra track per R-coded obs
K_base <- K_base + n_resprout_extra
```

Under the correct model, R-coded rows reuse the pre-existing track slot of the old stem —
they do not need fresh slots. Only genuinely new sibling boles need extra tracks, and those
are already counted by the standard `births_needed` computation. The extra K inflation is
unnecessary and inflates the state space, increasing computation time.

---

### Bug 4 — NA-R rows with no prior stem get no method label

**Location:** `dp_global_dp.R` ~line 2107, NA-R barrier loop

When the R census is the earliest census in the DP range there are no pre-barrier stems to
link the NA-R row to. The row gets `ReconstructedStemID = NA` and no `ReconstructionMethod`
— it looks like a missing assignment rather than a deliberate "broken stem with no prior
history in range."

---

### Bug 5 — Pin vs R-recruit constraint conflict (live R rows with TrueStemID)

**Location:** `dp_global_dp.R` ~line 1807 and ~line 2075

When a live R-coded row has a known `TrueStemID` and `pin_truestemid = TRUE`:
- The pin forces the stem onto its known track (correct).
- The R-recruit constraint demands that track be empty at the preceding census (wrong).

These two constraints are mutually exclusive. The R-recruit constraint wins silently by
pruning all transitions involving the pinned track → fallback. The NA-R post-processor
then refuses to re-split because of the pin guard, logging:

```
WARNING: pinned TrueStemID X crosses NA-R barrier — respecting pin (not reassigning)
```

Result: either fallback to probabilistic or identity incorrectly severed/merged.

---

## Examples

### Example 1 — Simple break, one new bole

```
Tag 001 | C1: Stem A (DBH=50)
Tag 001 | C2: Stem A (DBH=51, code=R),  Stem B (DBH=3)
Tag 001 | C3: Stem B (DBH=4)
```

**Expected reconstruction:**
- Stem A: C1 DBH=50 → C2 DBH=51_R  →  StemID=1 (both censuses)
- Stem B: C2 DBH=3  → C3 DBH=4      →  StemID=2

**Current output (wrong):** R-recruit constraint prunes the A→A transition. DP may
pair Stem B (new) with StemID=1 and create a spurious ID for Stem A's C1 record, or
fall back to probabilistic entirely.

---

### Example 2 — R row with known TrueStemID (pin conflict, Bug 5)

```
Tag 002 | C1: Stem A (DBH=30, TrueStemID=5)
Tag 002 | C2: Stem A (DBH=31, code=R, TrueStemID=5),  Stem B (DBH=2)
Tag 002 | C3 (anchor): Stem B (DBH=3, TrueStemID=6)
```

`pin_truestemid=TRUE` forces C2 Stem A onto track 5 (correct).
R-recruit constraint demands track 5 be empty at C1 (wrong).
Constraints are mutually exclusive → all feasible states pruned → fallback to probabilistic.

---

### Example 3 — NA-DBH R row (handled correctly — reference case)

```
Tag 003 | C1: Stem A (DBH=20)
Tag 003 | C2: Stem A (DBH=NA, code=R)   ← break recorded, not measured
Tag 003 | C3: Stem B (DBH=2)
```

The NA-R barrier path severs tracks across the barrier and links the C2 NA-R row to Stem
A's pre-barrier ID. **This case works correctly.** The fix must not break it.

---

### Example 4 — NA-R row with no prior census (Bug 4)

```
Tag 004 | C1: Stem A (DBH=NA, code=R)   ← first census in DP range, no history
Tag 004 | C2: Stem B (DBH=5)
```

`length(.cens_before) == 0L → next`. C1 Stem A gets `ReconstructedStemID = NA`,
no `ReconstructionMethod` label.

---

## Proposed Fix

### Fix 1 — Invert the R-recruit constraint (~line 1817)

Replace "R-coded obs must be on an empty track" with the opposite:
**R-coded obs must continue an occupied track.**

```r
# Identify unpinned R-coded observations (pinned ones are already correct via the pin)
.r_pos_p1_unpinned <- .r_pos_p1
if (isTRUE(pin_truestemid) && !is.null(pin_tidx_at_census[[p + 1L]])) {
    .pins_p1 <- pin_tidx_at_census[[p + 1L]]
    .r_pos_p1_unpinned <- .r_pos_p1[is.na(.pins_p1[.r_pos_p1])]
}

# NEW: R-coded obs at p+1 must occupy a track that WAS occupied at p
.keep_r <- rep(TRUE, n_feasible)
for (.rj in seq_along(.r_pos_p1_unpinned)) {
    .rt <- .r_tracks[, .r_pos_p1_unpinned[.rj]]
    .was_occupied <- matrix(FALSE, nrow = n_feasible, ncol = .n_obs_p)
    for (.pc in seq_len(.n_obs_p)) {
        .was_occupied[, .pc] <- .rt == .p_assign[, .pc]
    }
    # KEEP pairs where the R track was occupied at p (rowSums > 0)
    .keep_r <- .keep_r & (rowSums(.was_occupied) > 0L)
}
```

---

### Fix 2 — Move the R census into the pre-segment (~line 1082)

```r
# Fixed:
.post_data <- tree_data[CensusID >  .r_boundary_census]  # census AFTER R
.pre_data  <- tree_data[CensusID <= .r_boundary_census]  # includes the R census
```

The pre-segment anchor becomes the R census itself (last observed census with DBH for old
stems). `post_segment_all_recruits = TRUE` remains correct for the post-segment because the
first census after R contains only genuinely new boles.

The `do_fallback` live-R boundary severing loop also uses this boundary direction and must
be updated consistently: change `CensusID >= .cc_fb` to `CensusID > .cc_fb` so that the R
census itself stays in the pre-boundary segment.

---

### Fix 3 — Remove K inflation (~line 1172)

Remove the `n_resprout_extra` block entirely:

```r
# Remove this block:
# n_resprout_extra <- sum(...)
# K_base <- K_base + n_resprout_extra
```

R-coded rows reuse existing tracks. The `births_needed` computation already accounts for
any genuine stem-count increases at and after the break census.

---

### Fix 4 — Label unlabelled NA-R rows (~line 2120)

After the NA-R barrier linking loop, stamp a method on barrier rows that remain unlinked:

```r
tree_data[
    CensusID == .cc_bar & is.na(DBH) &
    !is.na(ListOfTSM) & grepl(resprout_regex, ListOfTSM, perl = TRUE) &
    is.na(ReconstructedStemID),
    ReconstructionMethod := "dp"
]
```

Apply the same stamp in the equivalent block inside `do_fallback`.

---

## Affected Code Locations

| Location | File | Line (approx.) | Change |
|----------|------|---------------|--------|
| R-recruit constraint | `dp_global_dp.R` | ~1807–1848 | Invert: must be occupied, not empty; skip pinned obs |
| Segment split boundary | `dp_global_dp.R` | ~1082–1086 | `>=` → `>` for post, `<` → `<=` for pre |
| K inflation (`n_resprout_extra`) | `dp_global_dp.R` | ~1172 | Remove block |
| NA-R barrier, unlabelled rows | `dp_global_dp.R` | ~2120 | Add `ReconstructionMethod := "dp"` |
| `do_fallback` live-R severing loop | `dp_global_dp.R` | ~613 | Boundary `>= .cc_fb` → `> .cc_fb` |
| Error-handler R loops in runners | `main_cpp.R`, `main_cpp_chunk.R` | ~775, ~763 | Same boundary direction fix |

---

## Notes

- The NA-R barrier path (Example 3) is a **separate code path** and appears correct for the
  case where pre-barrier stems exist. Verify it remains correct after the live-R fix.
- `R`, `RP`, `RF`, `RT`, `OR` all mean the same thing for identity reconstruction: **the
  coded row IS the old stem's last record**; it must continue the pre-break track. This is
  confirmed by their shared `DFstatus = "broken below"` and by the BCI census data (all five
  codes appear exclusively in `broken below` rows in `codes_identification.csv`).
- `QR` **does not appear in the BCI census data or the ForestGEO code reference**
  (`forestgeo_codes_reference.qmd`). It is included in all four `resprout_regex` patterns in
  `dp_global_dp.R` but its presence there is unverified. Before applying the fix, confirm
  whether `QR` is a valid site-specific code or a historical artefact, and remove it from the
  regex if it is not present in the target dataset.
- `OR` is "other breakage" (partial splits, hinge breaks), not "other breakage resprout" — it
  does not imply a resprout is present. Its identity semantics are identical to `R`: the row is
  the old stem, not a new organism.
- After fixing the DP, re-run tag 005558 (14 stems at C4, previously memory-crashed) and
  any tag known to have R events before the anchor to validate correctness.
