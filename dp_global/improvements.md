# Improvements

This file collects forward-looking design plans that affect the engine's
output schema or numerical conventions but are **not yet implemented**.
Each plan is self-contained: it states the motivation, the proposed
behaviour, the precise touch-points in the codebase, the expected
backward-compatibility impact, and a validation strategy.

> Implementation status: **plan only — no code changes**.

---

## Plan 1 — Reverse the direction of `ReconstructedStemID` numbering

### 1.1 Motivation and current behaviour

In the BCI dataset (and in the convention shared by most ForestGEO plots),
`OriginalStemID` was assigned **forward in time as new stems were
encountered**: a stem first measured in census 1 received a smaller
integer than a stem first measured in census 5. The numerical order of
`OriginalStemID` therefore mirrors the chronological order of stem
appearance.

The reconstruction engine works **backward in time** from a late-census
anchor (default `ANCHOR_START_CENSUS = 7`). Stem identities are derived
from the anchor and propagated to earlier censuses. Two consequences
follow from how `ReconstructedStemID` values are minted today:

1. **At the anchor census**, `ReconstructedStemID` is initialised from
   `TrueStemID` where available, then padded with sequential integers
   starting one above the maximum known `TrueStemID` — see
   [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R#L150-L162)
   (lines 150–162):

   ```r
   anchor_ids[has_true] <- as.integer(anchor_obs$TrueStemID[has_true])
   next_id <- max(anchor_ids[has_true]) + 1L
   anchor_ids[!has_true] <- seq.int(next_id, length.out = n_missing)
   ```

2. **At NA-R resprout barriers**, when the engine needs to mint a fresh
   identity for a pre-barrier track that the DP severed from its
   post-barrier continuation, it uses ascending `max + 1L` allocation —
   see [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R#L2215-L2280)
   (around lines 2215–2280):

   ```r
   .cur_max_id <- max(tree_data$ReconstructedStemID, na.rm = TRUE)
   for (.old_id in .crossing) {
       .new_id <- as.integer(.cur_max_id) + 1L
       .cur_max_id <- .new_id
       ...
   }
   ```

3. The probabilistic matcher's death-track allocator uses the same
   pattern — see [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R#L1054)
   (line 1054):

   ```r
   death_base <- max(anchor_ids, na.rm = TRUE) + 1L
   ```

The net effect is that **engine-minted IDs are always larger than any
anchor or `OriginalStemID`**. A stem first observed at census 1 and
extended backward through several barriers can carry a `ReconstructedStemID`
that is *larger* than a stem first observed at census 6, which is the
opposite of the chronological convention encoded in `OriginalStemID`.

The user's request is to **make engine-minted IDs follow the same
chronological convention as `OriginalStemID`**: smaller integer ⇒ earlier
appearance (going back in time, **descending** from the anchor's largest
known id).

### 1.2 Target behaviour

After the change, a complete `(Tag, Stem)` history reconstructed by the
engine should satisfy the following invariant for any **single tag**:

> If stem A's earliest non-NA `ReconstructedStemID` row appears at a
> census strictly **earlier** than stem B's earliest non-NA
> `ReconstructedStemID` row (within the same `Tag`), then
> `min(ReconstructedStemID(A)) < min(ReconstructedStemID(B))`.

Equivalently, when the engine needs to mint a brand-new identity for a
pre-anchor track that does not match any anchor stem, it should
**descend** from the smallest currently-used id rather than ascend from
the largest. The descent base is the anchor's smallest non-NA
`TrueStemID`, and each new mint uses `min(...) - 1L`.

This is consistent with the source-data convention because `OriginalStemID`
itself was assigned forward in chronological order: the earliest stems in
the dataset have the smallest ids. When the engine fabricates a stem
identity that the field crew never recorded (because the stem died before
the anchor was assigned), placing it numerically **below** the anchor ids
keeps the chronological-vs-numerical correspondence intact.

### 1.3 Scope: what this changes and what it does not

**Changes:**

- The numerical value (and width / sign convention, see §1.5) of
  `ReconstructedStemID` for rows with `ReconstructionMethod ∈ {"dp",
  "probabilistic", "provisional_dp", "dp_mf_inferred"}` whose chosen
  identity is **not** carried by any input `TrueStemID` /
  `OriginalStemID` / `StemID` value — i.e. only the engine-minted ids
  change.
- The numerical value used for the **NA-R barrier severance** mint inside
  [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R#L2259) and the
  **death-track base** inside
  [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R#L1054).
- The padding for anchor rows that lack `TrueStemID` in
  [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R#L150-L162).

**Does not change:**

- Any row whose `ReconstructionMethod` is `"given"`, `"given_orphan"`, or
  `"carried_terminal"` — those are anchored to source ids
  (`TrueStemID` / `OriginalStemID` / `StemID`) and must retain their
  source numerical value.
- The `OriginalStemID` and `TrueStemID` input columns.
- The structural decisions of the DP / probabilistic engine (which rows
  belong to the same track) — only the **labels** assigned to the tracks
  are renumbered.
- The `ReconstructedStemID_PreSweep` audit column. It mirrors whatever
  the engine assigned, which under this plan will be a descending mint
  rather than an ascending one.

### 1.4 Algorithm: descending allocator

Replace every `max(...) + 1L` mint site with a `min(...) - 1L` allocator
that maintains a per-tag descending counter. Concretely:

1. **Initialisation (per tag, at the start of `run_dp_one_group()`):**
   compute the descent base from the anchor and from `OriginalStemID`:

   ```r
   .ts_at_anchor <- tree_data[CensusID == anchor_start & !is.na(TrueStemID),
                              TrueStemID]
   .src_global    <- tree_data[, OriginalStemID]   # or StemID in production
   .descent_base  <- suppressWarnings(min(c(.ts_at_anchor, .src_global),
                                          na.rm = TRUE))
   if (!is.finite(.descent_base)) .descent_base <- 0L
   ```

   Store `.descent_base` on the tag-local environment so all three mint
   sites share the same counter.

2. **Per-mint allocator** (replaces both `max(...) + 1L` sites):

   ```r
   .cur_min_id <- .descent_base
   .new_id <- as.integer(.cur_min_id) - 1L
   .cur_min_id <- .new_id     # update the counter
   ```

   The **first** mint inside a tag is therefore one less than the smallest
   anchor / `OriginalStemID` in that tag's history; subsequent mints
   continue to descend.

3. **Anchor padding** (the
   [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R#L160)
   `seq.int(next_id, length.out = n_missing)` line):

   ```r
   anchor_ids[!has_true] <- seq.int(.descent_base - 1L,
                                    by = -1L,
                                    length.out = n_missing)
   ```

   The descending sequence starts one below `.descent_base` and steps
   downward (`by = -1L`).

### 1.5 Negative-id and sign considerations

If `OriginalStemID` (or `TrueStemID`) for a tag includes the smallest
possible plot-wide id (often `1L`), a single descending mint will produce
`0L`, and additional mints will become negative. There are three
acceptable resolutions; the implementation should pick exactly one and
document it loudly:

**Option A — allow negative integers.** Cleanest. Engine-minted ids carry
the chronological-order semantics for free, and the sign acts as a
visible audit signal that the row was never anchored to a source id.
Downstream BA accounting is unaffected because BA aggregations group on
`(Tag, ReconstructedStemID)` and never interpret the integer as a count.

**Option B — descend within a per-tag negative offset block.** Allocate
each tag a private id space such as
`-Tag * 1000L - <local descending counter>`. This keeps engine-minted ids
disjoint across tags and easily greppable, but breaks the global
chronological-order promise across tags.

**Option C — descend from a large positive offset reserved for the
engine.** E.g. start at `100000000L - <local descending counter>` per tag.
Avoids negatives but places engine-minted ids *above* `OriginalStemID`,
which contradicts the original motivation.

**Recommendation: Option A.** It is the only option that fulfils the
chronological-order invariant in §1.2 with no caveats. Under Option A,
a downstream sanity assertion — that no row carrying
`ReconstructionMethod ∈ {"given", "given_orphan", "carried_terminal"}`
has a non-positive `ReconstructedStemID` — gives a one-line check that
the source-id rows were never re-minted.

### 1.6 Touch-points (precise)

| File | Lines | Change |
|------|------:|--------|
| [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R) | 2215, 2259 | Replace `max(... ReconstructedStemID, na.rm = TRUE)` → descent base; replace `as.integer(.cur_max_id) + 1L` → `as.integer(.cur_min_id) - 1L`; rename `.cur_max_id` → `.cur_min_id`. |
| [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R) | 150–162 | Replace `next_id <- max(anchor_ids[has_true]) + 1L` and `seq.int(next_id, length.out = n_missing)` with the descending-base allocator from §1.4 step 3. |
| [dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R) | 1054 | Replace `death_base <- max(anchor_ids, na.rm = TRUE) + 1L` with `death_base <- min(c(anchor_ids, OriginalStemID), na.rm = TRUE) - 1L` (or read from the shared per-tag counter if §1.4 step 1 promotes it to a closure-level variable). |
| [dp_global/R/dp_global_main.R](dp_global/R/dp_global_main.R) | new | Add a tiny helper, e.g. `make_descent_allocator(tree_data, anchor_start)`, returning a stateful closure with two methods — `next_id()` and `seq_n(n)` — so all three mint sites use the same counter and the same definition of `.descent_base`. This keeps the rule in one place and avoids drift if a future fourth mint site is added. |
| [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R) | new | At the start of `run_dp_one_group()` (or at the top of the per-tag preprocessing block), construct the allocator and pass it as an extra argument to `match_stems_probabilistic()` and to the NA-R barrier block. |

### 1.7 Backward-compatibility strategy

Since the change alters the **integer label** on engine-minted rows,
every downstream consumer that joins on `ReconstructedStemID` will see
different keys. To avoid silently breaking pinned anchor tests and to let
users run side-by-side comparisons, gate the new behaviour behind a CLI
flag:

| Flag | Default | Effect |
|------|---------|--------|
| `--RECONSTRUCTED_ID_DIRECTION` | `"ascending"` (legacy) | `"ascending"` reproduces today's behaviour. `"descending"` enables Plan 1. |

Plumb the flag into `dp_global/scripts/main_cpp.R`, propagate it through
`run_dp_one_group()` into the allocator builder, and document it in
[dp_global/scripts/README.md](dp_global/scripts/README.md) under the CLI
flag table.

Once `"descending"` has been validated end-to-end on the BCI test set and
the simulated-data regression suite, switch the default and keep
`"ascending"` available as an opt-out for one release cycle, then remove
it.

### 1.8 Validation plan

1. **Unit-level invariant test (per tag):**
   - Run any tag with at least one engine-minted id (e.g. BCI tag
     `258411`, where Fix 2 left `ReconstructedStemID = 995113` on C6
     row 3).
   - Assert that the smallest non-NA `ReconstructedStemID` in the tag's
     output equals `min(c(TrueStemID, OriginalStemID), na.rm = TRUE) - K`
     for some `K >= 0`.
   - Assert that every row with `ReconstructionMethod ∈ {"given",
     "given_orphan", "carried_terminal"}` carries a `ReconstructedStemID`
     equal to its own source id (unchanged from current behaviour).

2. **Chronological-order invariant test (per tag):**
   - For each pair of distinct `ReconstructedStemID` values in the tag,
     compute the earliest census at which each appears.
   - Assert that the order of `min(ReconstructedStemID)` matches the
     order of those earliest censuses (ascending census ⇒ ascending id).

3. **Regression test against the current 26-tag harness:**
   - Re-run `/tmp/run_tag_list.R --TAG_LIST=/tmp/twenty_six_tags.txt`
     under both `--RECONSTRUCTED_ID_DIRECTION=ascending` and
     `--RECONSTRUCTED_ID_DIRECTION=descending`.
   - Build a row-by-row mapping from old `ReconstructedStemID` →
     new `ReconstructedStemID` per `(Tag, CensusID, OriginalStemID,
     StemTag)`. Assert that the mapping is a **bijection within each
     tag** (i.e. no two old ids collapse to one new id, and vice versa).
     This proves the renumbering is purely a relabelling.
   - Re-run the duplicate-Recon-per-(Tag, CensusID) check from this
     session: it must still report **0** duplicates.

4. **Downstream BA-uncertainty test:**
   - Run `dp_global/scripts/basal_area_uncertainty.R` on both outputs.
   - Tag-level `basal_area_tag_census.csv` must be **byte-identical**
     (totals are identity-invariant).
   - `basal_area_tag_change.csv` decompositions (growth / loss / gain)
     must be statistically equivalent — the posterior means and CIs may
     differ only by Monte Carlo noise from re-sampling. Confirm with a
     two-sided Kolmogorov–Smirnov test on the per-interval `DeltaBA`
     samples.

5. **Probabilistic-matcher death-track sanity:**
   - On any tag that exercises the probabilistic fallback (e.g. BCI tags
     with 7+ stems per census), assert that no death-track id collides
     with any anchor id, with any `OriginalStemID`, or with any
     barrier-mint id from §1.4 step 2. This is the descending-allocator
     analogue of today's `death_base = max(...) + 1L` non-collision
     guarantee.

### 1.9 Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Negative `ReconstructedStemID` breaks a downstream consumer that assumes positive ids. | The CLI flag (§1.7) keeps the legacy ascending default until every consumer is audited. The `apply_orphan_stem_backfill()` and `apply_carried_terminal_backfill()` helpers in [dp_global/R/dp_global_main.R](dp_global/R/dp_global_main.R) read `ReconstructedStemID` only via `is.na()` and assignment, so they are unaffected. The basal-area uncertainty script groups by `ReconstructedStemID` without arithmetic, so it is also unaffected. |
| The descent base is computed per tag, but a downstream join keys on `(Plot, ReconstructedStemID)` globally, expecting non-collision across tags. | Today's ascending allocator already collides across tags whenever two tags share an `OriginalStemID` numerical range — joins are always done on `(Tag, ReconstructedStemID)`. Document this explicitly in the new README section so the assumption is not silently inherited. |
| The NA-R barrier severance mint and the probabilistic death-track mint are currently independent counters. Sharing a single descending counter across them changes their numerical sequence in non-trivial ways. | The bijection regression test in §1.8 step 3 is exactly the safety net for this — any non-bijection is a hard failure of the renumbering plan and signals a deeper structural assumption that must be examined before merging. |
| `min(c(TrueStemID, OriginalStemID), na.rm = TRUE)` over a tag may be `Inf` (all NA). | Guard with `if (!is.finite(...)) <- 0L` (matches the current `max` guard in [dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R#L2216)). The first mint then becomes `-1L`. |

### 1.10 Out of scope (explicit non-goals)

- Renumbering pre-existing on-disk run outputs. The plan only changes
  future runs; older `*.rds` / `*.csv` outputs remain valid under the
  legacy convention.
- Changing the column type of `ReconstructedStemID`. It remains
  `integer`; only the value distribution changes (negative values become
  possible under Option A).
- Changing the allocation strategy inside the C++ transition-cost code
  ([dp_global/src/transition_cost_rcpp.cpp](dp_global/src/transition_cost_rcpp.cpp)).
  That code operates on track *indices* (1..K), not on emitted
  `ReconstructedStemID` values, and is unaffected.

---

## Plan 2 — Generalise the BCI Steps 1–3 pre-DP `TrueStemID` propagation

### 2.1 Motivation

The BCI debug driver
[dp_global/scripts/main_cpp_bci.R](dp_global/scripts/main_cpp_bci.R)
inlines ≈250 lines of pre-DP propagation logic (Steps 1, 2, 2b, 2c,
2d, 3, 3a) that build a clean `TrueStemID` column from
`OriginalStemID`, terminal events, and last-DBH census. Step 3a.5
(same-OS continuity, Fix 1) has been **partially ported** to
[dp_global/scripts/main_cpp_chunk.R](dp_global/scripts/main_cpp_chunk.R)
as a guarded no-op block: it runs only when `TrueStemID` already
exists in the input, so it is active for simulated inputs (and
correctly does nothing there because `TrueStemID` is never all-NA per
group), but it is a **no-op for raw BCI input** because Steps 1–3
are still missing — there is no `TrueStemID` column to anchor onto.

Result: a chunked production run on raw BCI data still emits the 5
Class-B born-resprout duplicates that Fix 1 resolves in
`main_cpp_bci.R`. Until Steps 1–3 are also ported, those tags must
be either (a) handled by pre-processing the BCI input through the
same Steps 1–3 logic before invoking `main_cpp_chunk.R`, or (b)
accepted as known residual duplicates and patched in a follow-up.

### 2.2 Proposal

Extract Steps 1, 2, 2b, 2c, 2d, 3, 3a, **and 3a.5** into a single
function `propagate_true_stem_id(dt, status_alive_set,
status_terminal_set, source_id_col)` placed in a new module
`dp_global/R/dp_global_preprocess.R`. Replace both the inline block
in `main_cpp_bci.R` and the partial 3a.5-only block in
`main_cpp_chunk.R` with a single call to the helper, guarded by a
flag `--PROPAGATE_TRUE_STEM_ID` (default `TRUE` for BCI input,
`FALSE` for synthetic data where `TrueStemID` is supplied by the
simulator).

### 2.3 Touch-points

- New file: `dp_global/R/dp_global_preprocess.R`.
- Loader: register in `load_dp_global()` in
  [dp_global/R/dp_global_main.R](dp_global/R/dp_global_main.R).
- Drivers: replace the inline block in `main_cpp_bci.R` and the
  partial Step 3a.5 block in `main_cpp_chunk.R` with the helper
  call; add the same call to `main_cpp.R` for symmetry.

### 2.4 Validation

Run the 26-tag regression harness with both flag values; with the
flag enabled in `main_cpp.R` / `main_cpp_chunk.R`, per-tag output
must be byte-identical to `main_cpp_bci.R` modulo the audit columns
added by the script-level sweep. Then run a full BCI chunked pass
and confirm the 5 Class-B tags (060145, 233660, 606162, 639010,
739002) are duplicate-free.

---

## Plan 3 — Lift duplicate-aware sweep into engine-level `finalize_out`

### 3.1 Motivation

Fix 2 (duplicate-aware pinning) currently lives in the script-level
post-engine sweep in
[dp_global/scripts/main_cpp.R](dp_global/scripts/main_cpp.R) and is
mirrored in
[dp_global/scripts/main_cpp_chunk.R](dp_global/scripts/main_cpp_chunk.R).
The probabilistic matcher
([dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R))
and the inner `finalize_out` sweep in
[dp_global/R/dp_global_dp.R](dp_global/R/dp_global_dp.R) can in
principle still emit duplicates of their own (the script-level sweep
currently catches them only because it runs after both). Defense in
depth would push the duplicate check down to every place a Recon is
written.

### 3.2 Proposal

Factor the rollback-vs-pin decision into a small helper
`commit_pin_or_rollback(out, rows, true_vec)` in
`dp_global/R/dp_global_main.R`, and call it from:

1. `finalize_out()` in `dp_global/R/dp_global_dp.R` (engine sweep).
2. The probabilistic matcher's id-emission step
   ([dp_global/R/dp_probabilistic_matching.R](dp_global/R/dp_probabilistic_matching.R)
   lines 150–162 and 1054).
3. The script-level sweep in both drivers.

### 3.3 Validation

Run the 26-tag regression harness — outcome must be unchanged. Then
remove the script-level sweep block (in a separate follow-up) and
re-run; outcome must still be unchanged. If it changes, the engine
sweep is missing a code path the script-level sweep was silently
patching.

---

## Plan 4 — Post-pass sanity assertion for residual duplicates

### 4.1 Motivation

The current pipeline depends on the layered ordering of Steps 5.5b /
5.5c / posterior_bins / script-level sweep to produce a duplicate-free
output. There is no terminal assertion — a future refactor that
reorders steps could silently re-introduce duplicates and the run
would still complete.

### 4.2 Proposal

Add `assert_no_duplicate_recon(out)` at the very end of `main_cpp.R`,
`main_cpp_chunk.R` (per-chunk), and `main_cpp_bci.R`. The check is

```r
dup <- out[!is.na(ReconstructedStemID),
           .N, by = .(Tag, CensusID, ReconstructedStemID)][N > 1]
if (nrow(dup) > 0L) {
    fwrite(dup, file.path(out_dir, "DUPLICATE_RECON_AUDIT.csv"))
    stop(sprintf("[assert] %d duplicate (Tag,CensusID,ReconstructedStemID) triples; see DUPLICATE_RECON_AUDIT.csv", nrow(dup)))
}
```

Cheap (single `data.table` aggregation), proportional to chunk size.

### 4.3 Validation

Drop into the current pipeline and re-run the 26-tag harness — must
not abort.

---

## Plan 5 — Vectorise the per-row R loop in Fix 2

### 5.1 Motivation

Fix 2 is implemented as a `for (.r in .ts_rows)` loop in
[dp_global/scripts/main_cpp.R](dp_global/scripts/main_cpp.R) (≈lines
907–924) and mirrored in `main_cpp_chunk.R`. The hot path is
`O(length(.ts_rows))` per chunk with R-level `c()` accumulation. For
chunks of ≈10⁵ pinnable rows this becomes the dominant cost.

### 5.2 Proposal

Replace the loop with a `data.table` group-aware operation:

```r
.dt <- data.table(.r = .ts_rows,
                  key = .grp_key[.ts_rows],
                  v   = .true_vec_all[.ts_rows])
# precompute, per (Tag,CensusID), the set of peer Recon values
# (engine-assigned baseline) and check membership in O(N log N)
```

Use `data.table`'s `set()` for the pin/rollback writes and accumulate
indices via `which(...)` rather than `c()`-growth.

### 5.3 Validation

Run the 26-tag harness; output must be byte-identical (including
audit columns). Then benchmark on a single 100k-row chunk and confirm
≥10× speedup.

---

## Plan 6 — Promote the regression harness into the repository

### 6.1 Motivation

The 26-tag regression suite lives in `/tmp/run_tag_list.R`,
`/tmp/twenty_six_tags.txt`, and `/tmp/cmp2.R`. These are essential
to re-validate Fix 1 + Fix 2 after any future change but are outside
the repo and therefore not version-controlled.

### 6.2 Proposal

Move into `tests/regression/`:

- `tests/regression/run_tag_list.R`
- `tests/regression/twenty_six_tags.txt`
- `tests/regression/expected/post_final.csv` (golden output)
- `tests/regression/cmp.R` (comparison + summary)
- `tests/regression/README.md` (one-screen run-and-interpret guide)

Add a `make regression` target in `Makefile` that runs the harness
and exits non-zero on any duplicate or control-row regression.

---

## Plan 7 — Unit test for `apply_orphan_stem_backfill` source-id auto-detect

### 7.1 Motivation

`apply_orphan_stem_backfill()` in
[dp_global/R/dp_global_main.R](dp_global/R/dp_global_main.R) auto-detects
whether the source-id column is `StemID` (production) or
`OriginalStemID` (this test repo). The detection logic is exercised
only implicitly via end-to-end runs, so a column-rename or schema
drift could silently flip the wrong column with no error.

### 7.2 Proposal

Add a focused unit test in `tests/testthat/test-orphan-backfill.R`:

- Synthetic 4-row `data.table` with both `StemID` and `OriginalStemID`
  columns; assert `StemID` wins.
- Synthetic 4-row `data.table` with only `OriginalStemID`; assert
  it is used.
- Synthetic 4-row `data.table` with neither; assert the function
  returns the input unchanged with an informative `message()`.

---

## Plan 8 — Surface empirical justification for Step 3a.5 in the README

### 8.1 Motivation

Step 3a.5 (same-OS continuity anchor — Fix 1) is justified by the
evidence-ratio analysis in
[bci_data/dead_pattern.qmd](bci_data/dead_pattern.qmd) and
[bci_data/broken_below_pattern.qmd](bci_data/broken_below_pattern.qmd)
(≈99.8% of dead transitions and ≈60% of broken-below transitions
follow the same-OS-continuity pattern). These numbers are not
referenced in the engine README, so a future maintainer reading
[dp_global/README.md](dp_global/README.md) sees the rule but not the
evidence.

### 8.2 Proposal

Add a one-paragraph "Empirical justification" subsection to
[dp_global/README.md](dp_global/README.md) immediately under the
Step 3a.5 description, citing the two qmd files and the
99.8% / ~60% evidence ratios. No code changes.

---

*End of Plans 1–8. Future plans should be appended below as new top-level
sections.*
