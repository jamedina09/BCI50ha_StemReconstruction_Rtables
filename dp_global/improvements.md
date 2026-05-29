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
