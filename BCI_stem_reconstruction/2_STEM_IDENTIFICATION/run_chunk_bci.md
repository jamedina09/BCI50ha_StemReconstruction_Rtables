# BCI Stem Identification: Chunked Run Instructions

This guide explains how to launch and manage chunked stem-identification runs for BCI data using the provided R scripts. All commands should be run from the **project root** (the directory containing both `dp_global/` and `BCI_stem_reconstruction/`). `caffeinate -i` (macOS) keeps the computer awake during the run; leave it out on other systems.

---

## Launch a New Full Run

To start a new chunked run, use the following command:

```sh
caffeinate -i Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe" \
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
- `PROB_SPECIES`: Comma-separated species codes (`Mnemonic`) routed to the probabilistic matcher: the clump-forming palms *Oenocarpus mapora* (`oenoma`) and *Bactris major* (`bactma`) and seven strangler figs (`ficu*`). Every code must occur in the data; otherwise the run stops at load with `❌ CHECK FAILED`.
- `DP_FALLBACK_GROWTH_FORMS`: `growth_form` labels routed to the probabilistic matcher. The labels are those the script assigns: `tree`, `shrub`, `palm`, `strangler`, `fern` (singular, matched exactly). A value that is not a label (e.g. `palms`) stops the run at load with `❌ CHECK FAILED` instead of silently routing nothing.
- `PROB_BIRTH_DEATH` (default `TRUE`, not needed on the command line): the probabilistic matcher lets every stem continue, die or be recruited in every census pair, scored by the same biological model as the DP. `--PROB_BIRTH_DEATH=FALSE` selects count-based slots, which force every stem to continue when the stem count does not change (before 2010 this gives 0.8–4.4 palm stem replacements, a death and a recruit in one tree, per 100 trees per year, against 4.6–6.9 with the 2010–2023 stem tags). The DP is not affected.
- `RECRUIT_RATE_UNIT` (default `"tree"`, not needed on the command line): the recruitment rate of both engines is the number of new stems per established tree per year (2010–2023 data). `--RECRUIT_RATE_UNIT=slot` selects the rate per empty slot of the estimation grid (~0.08/yr for every species).
- `COVERAGE_RULE` (default `"union"`, not needed on the command line): which species get their own parameter set. A species with enough 2010–2023 data (≥ 20 trees, ≥ 20 growth pairs) needs diameters that fill 4 equal log-DBH bins with ≥ 3 measurements each. `"union"` accepts a species whose bins are filled either from its smallest to its largest diameter or over its middle 95% of diameters, so one out-of-range value (*Beilschmiedia pendula*: a single 2 mm record) does not send a well-sampled species to a pooled set. In BCI it gives *beilpe*, *annoac*, *bactma*, *ast1st*, *hameax* and *pipeco* their own sets (5,887 trees of the engine) and removes none; the run log lists them (`✓ COVERAGE_RULE union …`). In a test against the 2010–2023 stem tags (identities hidden, parameters cross-validated by tree), 65.1% of the links of the affected trees are recovered with `union` against 60.0% with `minmax` (*Bactris major* 43.4% against 20.9%; *Beilschmiedia* 91.1% against 92.8%). `--COVERAGE_RULE=minmax` selects the rule with bins from the smallest to the largest diameter only.
- `POSTERIOR_SAMPLE_SEED`: Integer seed for reproducible sampling.
- `MANUAL_CORES`/`MANUAL_CORES_VALUE`: Enable and set the number of parallel workers.
- `DP_CHUNK_SIZE`: Number of tags processed per parallel chunk (match to core count for efficiency).
- `USE_MEASUREMENT_ERROR`: Whether the DBH measurement-error model is used in the parameter estimation and in the DP growth likelihood. The script default is `TRUE`; the command above turns it off (folder name `NME`).
- `DBH_ROUND_CENSUSES` (default `"1,2"`, not needed on the command line): censuses whose small-stem DBH (< 55 mm) was recorded in 5 mm classes, rounded down (BCI 1982 and 1985). The engines read such a DBH as a size within its class, so a 5 mm class step is not taken as an impossible growth jump that splits one stem into a death and a recruit. The run log shows `✓ DBH rounded down to 5 mm classes …` (the data are checked at run time) and the output folder name carries `NME_R12` (with the defaults above, `NME_R12_BD_LT_CU`: `BD` = birth-death matcher, `LT` = recruitment rate per tree, `CU` = coverage rule `union`). `--DBH_ROUND_CENSUSES=none` turns it off.
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
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe" \
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
