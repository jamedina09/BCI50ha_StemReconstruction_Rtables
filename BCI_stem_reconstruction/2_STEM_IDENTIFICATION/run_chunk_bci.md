# Engine run (BCI)

How to run `1_main_cpp_chunk_bci.R`, the stage-2 engine driver, and what its
settings do. Run every command from the project root. `caffeinate -i` (macOS)
keeps the computer awake during the run; leave it out on other systems.

## Table of contents

- [How to run](#how-to-run)
- [Settings](#settings)
- [Output folder](#output-folder)
- [Resume a run](#resume-a-run)
- [Blind test](#blind-test)
- [Check the DP state count of a tag](#check-the-dp-state-count-of-a-tag)

## How to run

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

The command sets four values that differ from the script defaults:
`DP_FALLBACK_GROWTH_FORMS`, `POSTERIOR_SAMPLE_SEED`, `USE_MEASUREMENT_ERROR`
and `BASE_OUT_DIR`. The other flags restate the defaults. Any setting below
can be given on the command line as `--NAME=value`.

When the run finishes (`run_finished.txt` in the run folder), run
`2_merge_chunks_to_datatable.R` (see `README.md`).

## Settings

Engine routing and posterior:

| Setting | Default | What it does |
| --- | --- | --- |
| `DP_MAX_STATES` | `10000` | Largest number of assignment states of a census that the DP solves; a larger tag goes to the probabilistic matcher |
| `PROB_SPECIES` | the 11 codes of the command | Species (`Mnemonic`) sent to the probabilistic matcher: the clump palms *Oenocarpus mapora* (`oenoma`) and *Bactris major* (`bactma`), and every species with a strangler growth form, eight strangler figs (`ficu*`) and *Oreopanax capitatus* (`oreoca`). `DP_FALLBACK_GROWTH_FORMS = "strangler"` already sends the stranglers there; listing them keeps them there if that option changes. A code not in the data stops the run at load |
| `DP_FALLBACK_GROWTH_FORMS` | none | Growth forms sent to the probabilistic matcher (labels `tree`, `shrub`, `palm`, `strangler`, `fern`, matched exactly); a value that is not a label (e.g. `palms`) stops the run at load instead of routing nothing |
| `PROB_BIRTH_DEATH` | `TRUE` | The matcher lets every stem continue, die or be recruited in every census pair, scored by the DP's biological model; `FALSE` (count-based slots) forces every stem to continue when the stem count does not change, which hides stem replacements in palm clumps. The DP is not affected |
| `PROB_EXPORT_RULE` | `"mea"` | The matcher exports its most representative sample; `"fewest_me"` exports, among the samples, one with the fewest links that shrink beyond measurement error |
| `POSTERIOR_SAMPLES` | `200` | Posterior samples per tag (`0` writes no path files) |
| `POSTERIOR_SAMPLES_FORMAT` | `"feather"` | Format of the posterior path files (`feather`, `rds`, `csv`) |
| `POSTERIOR_SAMPLE_SEED` | `123` | Seed of the posterior samples |
| `DP_POSTERIOR_LABEL_FREE` | `TRUE` | The DP posterior counts each reconstruction once: a track with no stem in the last census and no pin is a label without identity. Otherwise a reconstruction counts once per labeling of those tracks, so the samples lean toward stems that end before the anchor and `DP_PosteriorReconstructedProb` of their rows is split across the labels. With it, those rows get the share of samples in which their stem ends at the same observation as in the export. The exported reconstruction is the same either way |

Biological model:

| Setting | Default | What it does |
| --- | --- | --- |
| `USE_MEASUREMENT_ERROR` | `TRUE` | Measurement-error mixture in the parameter estimation and the DP growth likelihood (`FALSE` gives the folder label `NME`) |
| `DBH_ROUND_CENSUSES` | `"1,2"` | Censuses whose DBH below `DBH_ROUND_MAX_MM` (55 mm) was recorded in `DBH_ROUND_WIDTH_MM` (5 mm) classes, rounded down (1982, 1985); the engines read such a DBH as a size within its class, so a one-class step is not read as an impossible growth jump that splits one stem into a death and a recruit. The driver checks the rounding in the data (`✓ DBH rounded down to 5 mm classes …`); `"none"` turns it off |
| `MAX_GROWTH_FIXED`, `MAX_SHRINK_FIXED` | `5`, `-0.5` | Hard bounds of annual growth (cm/yr) |
| `RECRUIT_MAX_FIXED` | `25.9999` | Largest DBH of a new stem (cm) |
| `PRUNE_BOUND_FACTOR` | `5` | Pruning bounds = factor × the hard bounds and the recruit limit |
| `K_SHRINK_FIXED`, `K_GROWTH_FIXED` | `0`, `0` | Soft shrinkage and growth penalties (off) |
| `RECRUIT_RATE_UNIT` | `"tree"` | Recruitment rate = new stems per established tree per year; `"slot"` = per empty slot of the estimation grid |
| `COVERAGE_RULE` | `"union"` | A species with enough 2010–2023 data (≥ 20 trees, ≥ 20 growth pairs) gets its own parameter set when its diameters fill 4 equal log-DBH bins (≥ 3 measurements each) from its smallest to its largest diameter or over its middle 95%, so one out-of-range diameter (e.g. a single 2 mm record of *Beilschmiedia pendula*) does not send a well-sampled species to a pooled set; the log lists the species added (`✓ COVERAGE_RULE union …`). `"minmax"` uses the full range only |
| `ANCHOR_START_CENSUS` | `7` | First census of known identity (2010); the parameters are estimated from it on |

Pins and breaks:

| Setting | Default | What it does |
| --- | --- | --- |
| `FIELD_TAG_FIRST_CENSUS` | `7` | Field tags are the `StemTag`s from this census on |
| `HARD_PINS` | `TRUE` | Pins hold after the engine: the broken-below rules never cut between two records with the same pin and never rewrite a pin, so each field-tagged stem (pinned before the engine, one pin per tag) stays one stem, its broken-below records without a DBH included. Exception, the stump rule: a stump, death or missing record without a DBH and without a field tag that comes before any life of its own `StemID` goes to the stem that broke, and so do all records of a `StemID` that never lives (up to 2010 the database stores a break under the new `StemID` of the resprout measured later) |
| `PIN_BB_REQUIRE_OTHER_RECORD` | `TRUE` | Step 3a pins a measured broken-below record only when its `StemID` has another record |
| `PIN_STOP_AT_GAP` | `TRUE` | Step 3b does not carry a pin across a census without a record of its `StemID` |
| `PIN_STOP_AT_FAST_LINK` | `TRUE` | Step 3c (below) |
| `BB_DBH_RULE` | `"main_ends"` | A broken-below record with a DBH ends the main stem of the census before and is matched like any other measurement, but cannot continue that main stem; `"r1"` makes it always a new stem |
| `DP_ALLOW_SEGMENT_SPLIT` | `FALSE` | `TRUE` makes the DP end every stem of a tree at its first broken-below record before the anchor |

Step 3c (`PIN_STOP_AT_FAST_LINK`): a pinned link of a `StemID` between two
consecutive censuses before 2010 is checked against the measurements. When its
growth is above the 99.9% quantile of tagged growth of the species, and another
record of the tree in the earlier census grows into the later record within the
99% quantile (that record unpinned, its `StemID` without a measured record in
the later census), the records of the `StemID` up to the earlier census lose
their pin and the engine places them. This catches a database `StemID` carried
from a small 2005 record to the 2010 stem of another 2005 record (tag 170028:
the 29 mm stem of 2010 carries the `StemID` of a 10 mm record, while the main
stem measured 29 mm in 2005). The log reports the count
(`✓ Step 3c PIN_STOP_AT_FAST_LINK: …`) and `pin_stop_at_fast_link.csv` in the
run folder lists the links.

Run control:

| Setting | Default | What it does |
| --- | --- | --- |
| `BASE_OUT_DIR` | `BCI_stem_reconstruction/output` | Folder in which the run folder is created |
| `OUT_DIR_OVERRIDE` | none | Writes into this folder instead (used to resume a run) |
| `DP_CHUNK_SIZE` | `18` | Tags per chunk |
| `MANUAL_CORES`, `MANUAL_CORES_VALUE` | `TRUE`, `18` | Number of parallel workers (`FALSE`: all cores but one) |
| `DP_CHUNK_RESUME` | `TRUE` | Skips chunks that have a `_done.txt` marker |
| `DP_CHUNK_START`, `DP_CHUNK_END` | none | Runs only this range of chunks |
| `TAG_FILTER_FILE` | none | Runs only the tags listed in the first column of this file (the parameters are still estimated from all tags) |
| `BLIND_MIN_CENSUS`, `BIO_ESTIMATION_START_CENSUS` | none | Blind test (below) |

## Output folder

The run folder is `<BASE_OUT_DIR>/<date>_<time>_<label>`. With the defaults and
the command above the label is
`unknown_allT_DP_MB_NME_R12_BD_LT_CU_NS_MS_PO_PG_PF_HP_LF_g5_sm0p5_kg0_ks0_rcpp`:

| Part | Meaning |
| --- | --- |
| `unknown` | no `CONFIG_NAME` |
| `allT` | all tags |
| `DP_MB` | DP with marginals and posterior bins |
| `NME` | no measurement-error model |
| `R12` | DBH classes in censuses 1 and 2 |
| `BD` | birth-death matcher |
| `LT` | recruitment rate per tree |
| `CU` | coverage rule `union` |
| `NS` | no whole-tree split at a break |
| `MS` | main-stem rule |
| `PO`, `PG`, `PF` | pin rules of Steps 3a, 3b and 3c |
| `HP` | hard pins |
| `LF` | label-free DP posterior |
| `g5`, `sm0p5` | growth bounds 5 and −0.5 cm/yr |
| `kg0`, `ks0` | soft penalties 0 |

`FM` would mark `PROB_EXPORT_RULE = "fewest_me"`.

The folder holds `run_started.txt`, `run_finished.txt`, `run_log.txt`,
`run_parameters_full.txt`, `pin_stop_at_fast_link.csv`, one
`stem_reconstruction_dp_global_rcpp_chunk_NNN.feather` with its `_done.txt`
marker per chunk (`_failed.txt` for a chunk that raised an error) and
`posteriors/` with one `tag_<Tag>_posterior_samples__paths.feather` per tag.

## Resume a run

Run the same command with the folder of the interrupted run; chunks with a
`_done.txt` marker are skipped:

```sh
caffeinate -i Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R \
  --OUT_DIR_OVERRIDE=/Users/medinaja/outputs_bci_stem_identification/<run folder> \
  --DP_CHUNK_RESUME=TRUE \
  --DP_MAX_STATES=10000 \
  --PROB_SPECIES="oenoma,bactma,ficuob,ficupo,ficuc2,ficubu,ficuc1,ficuci,ficupe,ficutr,oreoca" \
  --DP_FALLBACK_GROWTH_FORMS="strangler" \
  --POSTERIOR_SAMPLE_SEED=42 \
  --MANUAL_CORES=TRUE \
  --MANUAL_CORES_VALUE=18 \
  --DP_CHUNK_SIZE=18 \
  --USE_MEASUREMENT_ERROR=FALSE
```

Keep every other setting as in the interrupted run. `DP_CHUNK_START` and
`DP_CHUNK_END` limit the run to a range of chunks.

## Blind test

`--BLIND_MIN_CENSUS=7 --ANCHOR_START_CENSUS=9 --BIO_ESTIMATION_START_CENSUS=7`
gives the engine censuses 7–9 only and keeps the pins of census 9: the 2010 and
2015 identities are hidden, so the reconstructed links can be scored against
the stem tags. It is a diagnostic; its output is not used downstream.

## Check the DP state count of a tag

The number of assignment states of a census with `n` measured stems on `K`
tracks is P(K, n) = K! / (K − n)!. A tag whose largest census exceeds
`DP_MAX_STATES` goes to the probabilistic matcher:

```r
n <- 5                         # largest number of measured stems in one census
K <- n + 2                     # tracks: one slack track and one for a new stem
prod(K:(K - n + 1))            # P(7, 5) = 2,520: DP
# n = 6, K = 8: P(8, 6) = 20,160 > 10,000: probabilistic matcher
```
