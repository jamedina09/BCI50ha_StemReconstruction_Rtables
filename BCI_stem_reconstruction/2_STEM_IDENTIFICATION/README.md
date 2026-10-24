# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. This stage runs the
`dp_global` engines on the cleaned ViewFullTable from
`BCI_stem_reconstruction/1_DATA_PREPARATION/` and produces the reconstructed
stem dataset used by `BCI_stem_reconstruction/3_PREPARE_R_TABLES/`. Both
scripts are run from the project root.

## Scripts

### `1_main_cpp_chunk_bci.R`

Chunked DP driver.

- Loads `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds`.
- Checks at load that every growth form and species code named in the routing
  settings (`DP_FALLBACK_GROWTH_FORMS`, `NON_TAPER_CORRECTED_GROWTH_FORMS`,
  `PROB_SPECIES`) occurs in the data; a value that does not (e.g. `palms` for
  the label `palm`) stops the run with `❌ CHECK FAILED`.
- Estimates the biological parameters from the 2010–2023 stems of all tags:
  one set per species with enough data and pooled sets per growth form for
  the others (all trees, all shrubs, palms + tree ferns, strangler figs). The
  growth SD is fitted as a constant when it falls with size, and the
  recruitment rate is new stems per established tree per year
  (`RECRUIT_RATE_UNIT`; see `dp_global/README.md`, *Parameter Estimation*). A
  species gets its own parameter set when its diameters cover the size range
  either from smallest to largest or over their middle 95%
  (`COVERAGE_RULE = "union"`), so a single out-of-range diameter does not
  send a well-sampled species to a pooled set.
- Pins the rows whose identity the database fixes (`TrueStemID`, Steps 1–3 in
  the script: stem tags, 2010+ records, measured break rows and the rows of
  the same StemID); the engines resolve the other rows. A measured break row
  is pinned only when its StemID has another record
  (`PIN_BB_REQUIRE_OTHER_RECORD`), and a pin is not carried across a census in
  which its StemID has no record (`PIN_STOP_AT_GAP`). Some StemID links made
  before 2010 are contradicted by the measurements: the link grows faster than
  99.9% of the species' tagged links, while another record of the tree, free
  and ending there, grows into the later record within 99% of them. The records
  before such a link lose their pin (`PIN_STOP_AT_FAST_LINK`, Step 3c: 616
  links, all 2005 → 2010, where the database gave the StemID of a small 2005
  record to the 2010 stem; see `run_chunk_bci.md`). The database StemID is
  used only to set these pins; the engines never read it. Pins hold after the
  engine (`HARD_PINS`), so every field-tagged stem stays one stem. Field tags
  are the StemTags from 2010 (`FIELD_TAG_FIRST_CENSUS`; a tag the database
  re-keyed gets one pin): a broken-below record without a DBH with a field tag
  is that stem with no measurable DBH at 1.3 m, so it stays with its tag. The
  exception keeps the stump rule for records without a field tag: a stump,
  death or missing record without a DBH that comes before any life of its own
  StemID (up to 2010 the database stores a break under the new StemID of the
  resprout measured later, and shows that StemID's later tag on it) still goes
  to the stem that broke, and so do all records of a StemID that never lives.
  The one other use of the StemID is for unmeasured records left
  without an identity: a death / stump / missing record joins the stem of the
  previous record with its StemID, and a stem never measured keeps its StemID.
- Breaks before 2015 (the BCI convention up to 2010: one broken-below record
  per tree, on the main-stem record, with the DBH of the largest living stem
  or no DBH): the main stem of the census before ends at the break, and the
  broken record is matched like any other measurement but cannot continue
  that main stem (`BB_DBH_RULE = "main_ends"`); the tree's other stems are not
  cut at the break (`DP_ALLOW_SEGMENT_SPLIT = FALSE`). See
  `dp_global/README.md`, *ForestGEO Code Handling*.
- Runs `dp_global` on multi-stem tags in parallel chunks: the exact DP, or the
  probabilistic matcher for palm clumps (*Oenocarpus mapora*, *Bactris major*),
  strangler figs and trees too complex for the DP. The matcher lets every stem
  continue, die or be recruited in every census pair (`PROB_BIRTH_DEATH`), and
  exports its most representative sample (`PROB_EXPORT_RULE = "mea"`; the
  option `"fewest_me"`, the sample with the fewest links that shrink beyond
  measurement error, recovered no more true links in the blind test and is
  not used).
- The DP posterior (marginals, posterior samples, `DP_PosteriorReconstructedProb`)
  counts each reconstruction once. Tracks with no stem in 2023 and no pinned
  measurement are interchangeable labels and are not counted once per labeling
  (`DP_POSTERIOR_LABEL_FREE`; see `dp_global/README.md`, *Choosing K*).
- A blind-test mode checks the engines against the 2010–2023 stem tags
  (`BLIND_MIN_CENSUS`, `BIO_ESTIMATION_START_CENSUS`; diagnostics only).
- Applies the post-engine helper chain to every chunk (pin-track rejoin,
  backfills, terminal-to-host, broken-below invariants, renumbering of each
  tag's stems to 1..N, posterior path files; see `dp_global/README.md`).
- Writes one Feather file per chunk and the posterior path files to
  `BASE_OUT_DIR/<run folder>/` (the folder name starts with the run's time
  stamp).
- Supports resuming interrupted runs.

### `2_merge_chunks_to_datatable.R`

Merges completed chunk outputs into final files.

- Reads Feather chunk outputs from `home_dir/run_code` (both set at the top
  of the script; the run must have finished).
- Converts them to temporary Parquet parts and merges them.
- Writes `merged_output.parquet` and `merged_output.rds` to
  `BCI_stem_reconstruction/DATA/<run_code>/`.
- Compares the merged table with the stage-1 table, restores the recorded DBH
  (mm) and adds the single-stem tags, whose one stem gets
  `ReconstructedStemID = 1` (`ReconstructionMethod =
  "single_stem_tag_no_reconstructed"`).
- Joins stems that the engine split only because the point of measurement
  (POM) of the trunk moved (buttresses, deformities, strangler figs; the new
  height recorded or not), so that a diameter change fell outside its hard
  growth bounds. An ended stem is joined to the stem of the same tree that
  starts in the next census (or two censuses later, when the trunk was not
  measured in between and its recorded height of measurement rose) when the
  new stem would be an impossible recruit (it starts above the recruit
  limit), the size change is the one a moved POM gives, the tree shows a sign
  of a moved POM (a height of measurement away from 1.3 m, a buttress code B,
  a strangler growth form, or, in the censuses without buttress codes, a
  species that usually has buttresses) and the ended stem is clearly its best
  candidate (same height of measurement, then the closest size). Recording
  errors and trunks only missed in one census are not joined; database
  StemIDs are not used (`apply_measurement_rejoin()`; rules and their
  justification in `dp_global/R/measurement_rejoin.R`, summary in
  `dp_global/README.md`; `ReconstructionMethod = "measurement_rejoin"`). The
  joined stem keeps one of its two IDs, so the tree's IDs skip one number per
  join.
- Restores the recorded `ExactDate`. The stage-2 driver fills a missing date
  (the mean date of the tag in its census, else of the census), because the
  engines and the rejoin need one to compute growth; after the rejoin the
  merge sets those dates back to the stage-1 value (NA), as it does for the
  DBH. Checks: stage 2 changed no recorded date, and every stage-1 column of
  the final table equals the stage-1 table (`growth_form` aside, the stage-2
  label on multiple-stem tags). Stage 3 gives every row a date in its own
  column `date`.
- Restores the pins (`TrueStemID`) of the rows relabelled by the broken-below
  rules R1/R2: with hard pins each is the row's database StemID, but in run
  20261022_220659 the stage-2 renumbering turned 127 of them (37 trees) into
  per-tree numbers (fixed since in `renumber_engine_minted_ids()`). Checks that
  every such pin is a database StemID of its tree. The reconstructed stems do
  not change.
- Checks that field-tagged stems stay whole and stops with `❌ CHECK FAILED`
  otherwise: every row with a field tag (2010 on) is pinned, every field tag
  (2010 on, by tag value) is one reconstructed stem, and no reconstructed stem
  holds two field tags. The table of run 20261014_220452 fails it with 703
  field tags on two stems (tagged 2010 stumps on the stem that broke, stems cut
  by rule R2, one re-keyed tag). Rows before 2010 whose copied StemTag lost its
  pin in Step 3c are counted in the log (93 in the 2026-10-24 sandbox).
- Writes the final reconstructed stem table
  `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` (the
  input of stage 3), the joins to `DATA/PROCESSED/measurement_rejoin_audit.csv`
  and the joined observation pairs, for the posterior samples, to
  `DATA/PROCESSED/measurement_rejoin_pairs.csv`.

## Other files

- `run_chunk_bci.md` — run and resume commands for `1_main_cpp_chunk_bci.R`.
- `tests/test_measurement_rejoin.R` — tests of the measurement rejoin
  (`testthat::test_file("BCI_stem_reconstruction/2_STEM_IDENTIFICATION/tests/test_measurement_rejoin.R")`).
- `comparissons/` (not under version control) — comparisons with the
  published BCI stem tables (Condit et al. 2019):
  - `compare_reconstruction_vs_dryad.R`: stem trajectories, link by link, and
    page-by-page PDFs per tag;
  - `compare_records_vs_condit.R`: the final R tables record by record
    (trees, records, dates, locations, measurements, statuses, stems), with
    the report `records_comparison_report.html` and a CSV of every
    difference;
  - `dryad_approach_weaknesses.qmd`.

## Requirements

R packages: `data.table`, `parallel`, `arrow` (Feather / Parquet), and those
loaded by `dp_global/R/dp_global_main.R` (`Rcpp`, `igraph`, `MASS`,
`lpSolve`); `testthat` for the tests.

## Data flow

```text
DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
        │
        ▼
1_main_cpp_chunk_bci.R
        └──▶ BASE_OUT_DIR/<run folder>/  (chunk Feather outputs, posteriors/)

2_merge_chunks_to_datatable.R
        ├──▶ BCI_stem_reconstruction/DATA/<run_code>/merged_output.{parquet,rds}
        └──▶ BCI_stem_reconstruction/DATA/PROCESSED/
               complete_dataset_final_with_reconstructed_stemids.rds
               measurement_rejoin_audit.csv, measurement_rejoin_pairs.csv
```
