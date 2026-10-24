# BCI Stem Identification: Chunked Run Instructions

This guide explains how to launch and manage chunked stem-identification runs for BCI data using the provided R scripts. All commands should be run from the **project root** (the directory containing both `dp_global/` and `BCI_stem_reconstruction/`). `caffeinate -i` (macOS) keeps the computer awake during the run; leave it out on other systems.

---

## Launch a New Full Run

To start a new chunked run, use the following command:

```sh
caffeinate -i Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe,ficutr,oreoca" \
  --DP_FALLBACK_GROWTH_FORMS="strangler" \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE \
  --BASE_OUT_DIR=/Users/medinaja/outputs_bci_stem_identification
```

**Flag explanations:**

- `DP_MAX_STATES`: Largest number of assignment states per census that the DP solves; a tag above it goes to the probabilistic matcher (higher = more tags solved exactly, longer runtime).
- `PROB_SPECIES`: Comma-separated species codes (`Mnemonic`) routed to the probabilistic matcher: the clump-forming palms *Oenocarpus mapora* (`oenoma`) and *Bactris major* (`bactma`) and every species with a strangler growth form: eight strangler figs (`ficu*`, including `ficutr`) and *Oreopanax capitatus* (`oreoca`). `DP_FALLBACK_GROWTH_FORMS="strangler"` already sends all strangler trees to the matcher; listing them here as well keeps them there if the growth-form option is changed. Every code must occur in the data; otherwise the run stops at load with `❌ CHECK FAILED`.
- `DP_FALLBACK_GROWTH_FORMS`: `growth_form` labels routed to the probabilistic matcher. The labels are those the script assigns: `tree`, `shrub`, `palm`, `strangler`, `fern` (singular, matched exactly). A value that is not a label (e.g. `palms`) stops the run at load with `❌ CHECK FAILED` instead of silently routing nothing.
- `PROB_BIRTH_DEATH` (default `TRUE`, not needed on the command line): the probabilistic matcher lets every stem continue, die or be recruited in every census pair, scored by the same biological model as the DP. `--PROB_BIRTH_DEATH=FALSE` selects count-based slots, which force every stem to continue when the stem count does not change (before 2010 this gives 0.8–4.4 palm stem replacements, a death and a recruit in one tree, per 100 trees per year, against 4.6–6.9 with the 2010–2023 stem tags). The DP is not affected.
- `RECRUIT_RATE_UNIT` (default `"tree"`, not needed on the command line): the recruitment rate of both engines is the number of new stems per established tree per year (2010–2023 data). `--RECRUIT_RATE_UNIT=slot` selects the rate per empty slot of the estimation grid (~0.08/yr for every species).
- `COVERAGE_RULE` (default `"union"`, not needed on the command line): which species get their own parameter set. A species with enough 2010–2023 data (≥ 20 trees, ≥ 20 growth pairs) needs diameters that fill 4 equal log-DBH bins with ≥ 3 measurements each. `"union"` accepts a species whose bins are filled either from its smallest to its largest diameter or over its middle 95% of diameters, so one out-of-range value (*Beilschmiedia pendula*: a single 2 mm record) does not send a well-sampled species to a pooled set. In BCI it gives *beilpe*, *annoac*, *bactma*, *ast1st*, *hameax* and *pipeco* their own sets (5,887 trees of the engine) and removes none; the run log lists them (`✓ COVERAGE_RULE union …`). In a test against the 2010–2023 stem tags (identities hidden, parameters cross-validated by tree), 65.1% of the links of the affected trees are recovered with `union` against 60.0% with `minmax` (*Bactris major* 43.4% against 20.9%; *Beilschmiedia* 91.1% against 92.8%). `--COVERAGE_RULE=minmax` selects the rule with bins from the smallest to the largest diameter only.
- Settings tested in the sandbox comparison of 2026-10-22 (the defaults are the values kept; not needed on the command line; each can be set back for a comparison run):
  - `DP_ALLOW_SEGMENT_SPLIT=FALSE`: the DP does not end every stem of a tree at its first broken-below record before 2010; only the broken stem is affected.
  - `BB_DBH_RULE="main_ends"`: a broken-below record with a DBH (before 2015, the break of the tree's main stem written on the main-stem record) ends the main stem of the census before, and is matched like any other measurement but cannot continue it. `"r1"` restores rule R1 (the record always starts a new stem).
  - `PIN_BB_REQUIRE_OTHER_RECORD=TRUE`: a measured broken-below record is pinned to its own StemID only when that StemID has another record.
  - `PIN_STOP_AT_GAP=TRUE`: a pin is not carried across a census in which its StemID has no record (e.g. tag 200431).
  - `PIN_STOP_AT_FAST_LINK=TRUE` (Step 3c, added 2026-10-24): a database StemID link made before the field tags can be contradicted by the measurements. Its growth is above the species' 99.9% quantile of tagged growth, while another record of the tree in the earlier census grows into the later record within the 99% quantile. That other record is unpinned, and its StemID has no measured record in the later census. The records of the StemID up to that census then lose their pin, and the engine places them.
    - The database had given the StemID of a small 2005 record to a 2010 stem that continues another 2005 record, which then ended in 2005. Example, tag 170028: the 29 mm stem of 2010 carries the StemID of a 10 mm 2005 record; the main stem measured 29 mm in 2005.
    - The run log reports the count (`✓ Step 3c PIN_STOP_AT_FAST_LINK: …`), and the links are listed in `pin_stop_at_fast_link.csv` in the run folder.
    - Sandbox of 2026-10-24 (1,494 trees; the unchanged code reproduces run 20261022_220659 in every column and posterior sample):
      - 616 links in 616 trees, all 2005 → 2010. In 616 of 616 trees the 2010 stem now continues a 2005 record that fits, and no link above the 99.9% quantile is left.
      - Median 2005 → 2010 growth on these links goes from 6.8 to 0.6 mm/yr.
      - Basal-area loss and gain from 2005 to 2010 in these trees fall from 5.35 / 5.67 m² to 0.30 / 0.63 m² (whole plot).
      - The posterior samples agree in every tree. The 829 other trees, near misses and random controls, are identical.
      - In 7 trees the change moves a break to 2000 → 2005, because a pre-2010 pin of the 2000 stem cannot join the 2010 pin.
      - An alternative tested only against the 99.9% quantile (665 links) also freed 8 ambiguous trees, which then got a new 2010 stem; hence the 99% test.
  - `HARD_PINS=TRUE`: a pin (`TrueStemID`) holds after the engine.
    - The broken-below rules never cut between two rows with the same pin, and never rewrite a pin.
    - Field-tagged stems (`StemTag` from `FIELD_TAG_FIRST_CENSUS = 7`, 2010; pinned before the engine, one pin per tag) therefore stay one stem, their broken-below records without a DBH included.
    - Exception, the stump rule: a stump, death or missing record without a DBH and without a field tag that comes before any life of its own StemID goes to the stem that broke, pinned or not; so do all records of a StemID that never lives. Up to 2010 the database stores a break under the new StemID of the resprout measured later.
  - `PROB_EXPORT_RULE="mea"`: the probabilistic matcher exports its most representative sample. `"fewest_me"` (not adopted) exports, among the samples, one with the fewest links that shrink beyond measurement error. In the blind test against the 2010–2023 stem tags it recovers no more true links than `"mea"` (43.1% vs 43.5%; precision 41.1% in both): the shrinking links it removes are replaced by other wrong links that pass the growth filter. The samples are the same either way.
  - `DP_POSTERIOR_LABEL_FREE=TRUE`: the DP posterior counts each reconstruction once.
    - The problem: a DP track with no stem in 2023 and no pinned measurement is a label without identity.
      - Without this setting, the posterior counts a reconstruction once per labeling of those tracks. The samples then lean toward stems that end before 2023, by an amount that depends on the number of spare tracks.
      - `DP_PosteriorReconstructedProb` of a row on such a track is split across the labels: median 0.5 in run 20261014_220452 where the label-free value is 1.
    - With this setting, those rows get the share of posterior samples in which their stem ends at the same observation as in the export.
    - The exported reconstruction is the same either way.
  - Blind test (diagnostics only): `--BLIND_MIN_CENSUS=7 --ANCHOR_START_CENSUS=9 --BIO_ESTIMATION_START_CENSUS=7` hides the 2010 and 2015 identities (only the 2023 pins are kept) so that the reconstructed links can be scored against the stem tags.
- `POSTERIOR_SAMPLE_SEED`: Integer seed for reproducible sampling.
- `MANUAL_CORES`/`MANUAL_CORES_VALUE`: Enable and set the number of parallel workers.
- `DP_CHUNK_SIZE`: Number of tags processed per parallel chunk (match to core count for efficiency).
- `USE_MEASUREMENT_ERROR`: Whether the DBH measurement-error model is used in the parameter estimation and in the DP growth likelihood. The script default is `TRUE`; the command above turns it off (folder name `NME`).
- `DBH_ROUND_CENSUSES` (default `"1,2"`, not needed on the command line): censuses whose small-stem DBH (< 55 mm) was recorded in 5 mm classes, rounded down (BCI 1982 and 1985). The engines read such a DBH as a size within its class, so a 5 mm class step is not taken as an impossible growth jump that splits one stem into a death and a recruit. The run log shows `✓ DBH rounded down to 5 mm classes …` (the data are checked at run time) and the output folder name carries `NME_R12`. `--DBH_ROUND_CENSUSES=none` turns it off.
- Output folder label: with the defaults above it is `NME_R12_BD_LT_CU_NS_MS_PO_PG_PF_HP_LF` (run 20261022_220659, before `PIN_STOP_AT_FAST_LINK`, has no `PF`):
  - `BD` = birth-death matcher;
  - `LT` = recruitment rate per tree;
  - `CU` = coverage rule `union`;
  - `NS` = no whole-tree resprout split;
  - `MS` = main-stem rule;
  - `PO` and `PG` = the two pin rules of 2026-10-22;
  - `PF` = no pin on the records before a fast database StemID link (`PIN_STOP_AT_FAST_LINK`);
  - `HP` = hard pins;
  - `LF` = label-free DP posterior;
  - `FM` would mark `PROB_EXPORT_RULE="fewest_me"`.
- `BASE_OUT_DIR`: Root directory for all output (a timestamped subdirectory is created for each run).

After all chunks finish, run `2_merge_chunks_to_datatable.R` to merge per-chunk Feather files into the final dataset.

---

## Resume a Partial Run

If a run is interrupted, you can resume it. Chunks with a `_done.txt` marker are skipped automatically.

```sh
caffeinate -i Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --OUT_DIR_OVERRIDE=/Users/medinaja/outputs_bci_stem_identification/<timestamp_run_code> \
  --DP_CHUNK_RESUME=TRUE \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe,ficutr,oreoca" \
  --DP_FALLBACK_GROWTH_FORMS="strangler" \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE \
  [--DP_CHUNK_START=<start>] [--DP_CHUNK_END=<end>]
```

- Replace `<timestamp_run_code>` with the actual directory name from your interrupted run.
- Keep every other setting as in the interrupted run, so that the remaining chunks are processed like the finished ones.
- Optionally specify `DP_CHUNK_START` and/or `DP_CHUNK_END` to process a specific range of chunks.

---

## Utility: Check Maximum DP States for a Tag

To estimate whether a tag's stem count will fit within your `DP_MAX_STATES` budget, run this in R:

```r
n <- 5          # max observed stems in any census for the tag
K <- n + 2      # number of tracks, here one slack track plus one track for a recruit (see dp_global/README.md, "Choosing K")
states <- prod(K:(K - n + 1))   # permutations P(K, n), the worst case for one census
cat("Max DP states for n =", n, ":", states, "\n")
# Example: n=5  →  P(7,5) = 2,520  (well within 10,000)
#          n=6  →  P(8,6) = 20,160 (exceeds 10,000 — the tag goes to the probabilistic matcher)
```

---

## Output Structure

- Each run creates a timestamped output directory inside `BASE_OUT_DIR`.
- Per-chunk Feather files (with `_done.txt` markers), the `posteriors/` folder with one path file per tag, and the run log and parameter files are written to this directory.
- After all chunks finish, merge them using `2_merge_chunks_to_datatable.R`.

---

For more details, see comments in the driver scripts and the main project README files.
