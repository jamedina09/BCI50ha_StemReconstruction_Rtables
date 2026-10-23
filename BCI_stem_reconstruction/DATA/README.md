# DATA: status and pending tasks

Last updated 2026-10-21. Run in use: `20261014_220452_unknown_allT_DP_MB_NME_R12_BD_LT_CU_g5_sm0p5_kg0_ks0_rcpp`.

## Folders

| Folder | Content |
|---|---|
| `RAW/` | database exports (`oct_13_2026`, `oct_16_2026`), species sources (`sp_tables/`), Condit et al. (2019) |
| `SPP_TABLE/` | stage-1 species table, Condit's table, comparison helper |
| `PROCESSED/` | stage-1 stem table, engine input |
| `<run folder>/` | stage-2 engine output and merge |
| `POSTERIORS/`, `RTABLES/` | stage-3 posterior preparation and R tables |
| `CHECKS/` | check files; `revision_oct16/` holds the Oct 16 review lists |
| `BACKUPS/` | the run of 2026-10-05 and its side files (zip archives); species-table outputs archived before the rebuild (`species_tables_20261021/`) |

## Known issues (reported, not waited for)

Decision 2026-10-21: the pipeline does not wait for the data manager or the
botanists. Each issue below is handled as stated, and kept here for the
release notes. Details and CSVs: `CHECKS/revision_oct16/`.

Data manager (database export):

- swars1 / swars2: variety labels shifted in the Oct 16 export. The pipeline reads no database names, so the R tables are not affected.
- protte (SubspeciesID 15) and quaras (SubspeciesID 13): duplicate taxonomy rows, ID without a name. Not read.
- 48 codes with old or mixed families (Fabaceae subfamilies, Myrsinaceae, Clusiaceae, turpoc Tapisciaceae). Not read: families come from the Flora de Panamá and the lists of email 1.
- 6 stems with two measurements at the same HOM and date; 3 stems with two StemIDs for one StemTag (trees 29078, 315977, 155532); tree 236797 stem 1a dead beside live stems. Kept as recorded.
- 1,123 stems and 1,267 trees dead and then alive; 34 dead records with DBH; 557 alive without DBH. Handled by the Rstatus rules (alive later means alive); `DFstatus` and `dbh` stay as recorded.
- 174 probable DBH typos (tree 1114: 2375 mm in 1985); 189 stems losing more than half their DBH. Kept as recorded.
- 177 records without date; 98 dates outside their census; tree 73712 stem 5 dated 2013. `ExactDate` stays as recorded (NA without a date); stage 3 fills `date`.
- 86 trees without coordinates; 4 trees with two locations; tree 550132 quadrat vs coordinates. One location per tree (stage 3); trees without coordinates keep NA.
- StemNumber (0 or empty from 2015), PrimaryStem empty from 2015, LargeStem empty, HOM on records without DBH. Not used.
- The partial census 9 (2020) is removed by `2_prepare_viewfulltable.R`.

Botanists (species table):

- Spelling: appuse *siebertii* (Flora de Panamá; IPNI *seibertii*), ardiba *bartletii*, rinosy *silvatica* Kuntze, malmsp *garwoodii*, pipeco *trinervium*: the Flora de Panamá spellings are kept.
- guargr, laetth, swars1, swars2: WFO treats them as synonyms; the Flora de Panamá names are kept (`CHECKS/species_table_vs_wfo.csv`).
- Synonym column: `CHECKS/species_synonyms_vs_wfo.csv` (2026-10-21: 20 synonyms that WFO places in another taxon, 3 not found, swars1 lists its own name).
- Authorities thevah and anacex: the WFO / IPNI authorities are adopted (section 9e of `1_prepare_species_tables.R`).
- Infraspecific rank: protte (subsp. *sessiliflorum*), myrcga ("tip. *gatunensis*"), as in the Flora de Panamá list.
- Morphospecies apeihy, nects1, nects3, sapisp: database labels (hybrid, fuzzy, sp.4_(tiny_leaf); sapisp "sp."), checked against email 1; growth form "árbol".
- Empty fields: common name (19 codes), herbarium (76 codes), left empty.

## To do: pipeline

- [x] Rebuild the species tables (2026-10-21; previous outputs in `BACKUPS/species_tables_20261021/`).
- [x] Export: stay on the Oct 13 export (2026-10-21). With the Oct 16 export `2_prepare_viewfulltable.R` gives the same engine input in every column the pipeline reads; only the database name columns differ (SpeciesName and SpeciesID on 117 rows, Subspecies on 57,687 rows).
- [x] Merge: the raw ExactDate is restored (2026-10-21; merge and stage 3 rerun in production: 16 cells of the final table and 66 of the R tables changed, stage 4 unchanged).
- [x] Stage 4 outputs: current (run 20261014_220452; identical after the date fix).
- [x] Archive the run of 2026-10-05: `BACKUPS/run_20261005_232052.zip` and `BACKUPS/not_run_20261005_20261014.zip` (2026-10-21); its chunk folder is no longer on disk.
- [x] `COMPARE_DATASETS/bci.stem.mine/` is gone.
- [x] Commit and push the species-table changes (2026-10-21: ad6eeca..23bac81) and the date fix (6be0bed).
- [x] Comparison with Condit et al. (2019) on the final R tables (2026-10-21): `COMPARE_DATASETS/` rerun (results identical to 2026-10-15) and `2_STEM_IDENTIFICATION/comparissons/compare_records_vs_condit.R` (records, dates, coordinates, statuses, stems; `records_comparison_report.html`).
- [ ] Publication folder (`data_paper_and_repo_publication/`) from run 20261014_220452.

## Identity fixes for the final engine run (sandbox test 2026-10-22)

**How they were tested:**
- Each fix is a setting of `2_STEM_IDENTIFICATION/1_main_cpp_chunk_bci.R`; the default is the value kept (list in `run_chunk_bci.md`).
- The sandbox ran the production command on 2,414 test tags: the trees each fix affects, 600 DP control trees, 120 stump-only trees and 200 matcher trees without breaks.
- With every setting set back, the sandbox reproduces run 20261014_220452 exactly.
- A fix was kept only when it improved its metric and left the control trees unchanged.
- Code and documents are not committed yet.

- [x] **Whole-tree resprout split off** (`DP_ALLOW_SEGMENT_SPLIT = FALSE`). At the first break before 2010 the DP no longer ends every stem of the tree. Results, with the main-stem rule:
  - the other stems continue 88% of the time (6.8% with the split);
  - +531 links, with the same share shrinking beyond measurement error (1.58%);
  - control trees unchanged.
- [x] **Main-stem rule instead of R1** (`BB_DBH_RULE = "main_ends"`; deterministic, no size threshold).
  - The data: up to 2010 the break is written on the main-stem record, with the DBH of the largest living stem or no DBH.
  - The rule: the main stem of the census before ends at the break. The broken record's DBH is matched like any other measurement of the tree but cannot continue that main stem; with one stem before, it is always a new stem.
  - Results: similar-size stems that end beside a new stem starting at the break record fell from 109 to 23; control trees unchanged.
- [x] **Step 3a pins a broken record only when its StemID has another record** (`PIN_BB_REQUIRE_OTHER_RECORD = TRUE`).
  - Similar-size pairs 2005–2010: 45 → 15.
  - +235 links; beyond measurement error 1.58% → 1.51%.
  - The links it creates are better than those it replaces (1.5% vs 2.1% beyond measurement error).
- [x] **Step 3b does not carry a pin across a census in which its StemID has no record** (`PIN_STOP_AT_GAP = TRUE`; tag 200431).
  - Records joined across such a gap: 57 → 32.
  - Joins with a more than 2-fold size change: 20 → 5.
  - Control trees unchanged.
- [x] **Not kept: matcher export with the fewest links beyond measurement error** (`PROB_EXPORT_RULE = "fewest_me"`).
  - It removes 86% of those links in matcher trees.
  - In the blind test (2010 and 2015 identities hidden, 839 tags) it recovers no more true links: 43.1% vs 43.5%; precision 41.1% in both; per tree, 41 better and 53 worse.
  - The shrinking links it removes are replaced by other wrong links that pass the growth filter.
- [x] **DP posterior counts each reconstruction once** (`DP_POSTERIOR_LABEL_FREE = TRUE`). Found 2026-10-22 while checking the posteriors; the problem is present in production.
  - The problem: a DP track with no stem in 2023 and no pinned measurement is a label without identity.
    - The DP counted each reconstruction once per labeling of those tracks, so its samples leaned toward stems ending before 2023.
    - `DP_PosteriorReconstructedProb` of a row on such a track was split across the labels: median 0.5 in production for unpinned stems ending before 2023, where the label-free value is 1.
  - Results, export identical for every tree and column:
    - those rows now median 1.0, equal to the share of samples that agree;
    - DP samples hold 0.29% fewer links than the export (0.43% without the correction; 0.25% in production); trees whose samples are off by half a link or more are now balanced (4 fewer, 3 more; 23 and 2 before);
    - with 1 vs 3 spare tracks (500 DP trees) the probabilities of stems alive in 2023 change in no row (53 without the correction), and the expected number of posterior links stays the same (2,730 vs 2,728; without the correction 2,725 vs 2,712);
    - the posteriors of 12 of the 600 control trees change: link frequencies by at most 0.095, expected links identical.
- [x] **Pins hold; field tags are kept** (`HARD_PINS = TRUE`, `FIELD_TAG_FIRST_CENSUS = 7`). User rules 2026-10-22: the database StemID only defines the pins; stems with a StemTag stay pinned; "a broken below with no dbh with a stemtag indicates there is a stem, just that there is no measurable dbh at 1.3".
  - The engines never read the StemID; the driver uses it only to set the pins. It is also used for unmeasured records left without an identity: a death / stump / missing record joins the stem of the previous record with its StemID, and a stem never measured keeps its StemID.
  - R1/R2 never cut between rows with the same pin and never rewrite a pin (run 20261014_220452 rewrote 4,724).
  - Field tags are the StemTags from 2010, the first census that tagged every stem, counted by tag value. The one tag the database re-keyed and split (523314 stem 2: StemID 890756 in 2010–2015, 1120187 in 2023) gets one pin (driver Step 1(a')).
  - A broken-below record without a DBH that carries a field tag stays with its tag's stem:
    - before the tag's first measurement: 674 records (666 in 2010), which run 20261014_220452 put on the stem that broke (e.g. tree 630310: the 2010 stump with tag 1 now starts the tag-1 stem measured from 2015, and the 2005 stem ends in 2005);
    - between two measurements: a stem measured, broken below and measured again under its tag is one stem; R2 no longer cuts it (60 stems in the whole data set, DBH after / before median 0.88).
  - The stump rule of 2026-09-30 still holds for records without a field tag: a stump, death or missing record without a DBH that comes before any life of its StemID is the end record of the stem that broke, and a resprout measured later is a new stem.
    - Before 2010 the database stores a break under a new StemID, the one the resprout measured later carries, and shows that StemID's later tag on all its records. Those tags are not field tags, so the 1,896 records that carry one keep the stump rule, as in run 20261014_220452.
    - All records of a StemID that never lives also go to the stem that broke, tagged or not; that stem never carries another field tag (0 of 6,077 records).
    - A first version that let every pin hold put 14,839 stump records on the resprout's stem, and stage 3 read 15,150 stump censuses as alive without a DBH.
  - The merge stops with `❌ CHECK FAILED` if a field tag (2010 on, by tag value) is on two stems or a stem holds two field tags. Run 20261014_220452 fails it with 703 field tags on two stems.
- [x] **Stump records stay with the stem that broke** (`apply_terminal_to_host()`):
  - a "missing" record after a stump on its StemID moves with the stump (run 20261014_220452 left 347 on another stem, e.g. 1985 stumps recorded "missing" in 1990);
  - a whole-tree break of a multi-stem tree up to 2010 (every stem ends together, so the rule "the one stem that ended" finds none): the stump goes to the main stem of the census before, where the field wrote the break. Run 20261014_220452 left 2,639 of these stumps (7,170 records) as stems that never live; in every case the main stem ends there and can take them.
- [x] **Stage 4: a stump is alive but out of the stock** (user decision 2026-10-22; stage 3 and the R statuses are unchanged). An A record without a DBH whose `DFstatus` is "broken below" is a stem with no measurable DBH at 1.3 m: the three stage-4 scripts give it no DBH, do not count it as a stem, and give it no basal area or biomass. Before the stem's first measurement the stem enters the stock at that measurement (a recruit); between two measurements its basal area or biomass leaves the stock and comes back. Without this, stage 4 would fill the tagged stumps above (alive in stage 3, because alive later) with a DBH as missed measurements. To be checked on the final run's R tables.
- **All kept settings together** (sandbox run of the final configuration):
  - links up to 2010 in trees and shrubs: 10,300 → 11,677; beyond measurement error 1.62% → 1.51%; growth above 5 mm/yr 2.63% → 2.47%;
  - stems ending beside a similar-size new stem: 122 → 40;
  - other stems continuing at the first break: 6.8% → 81%;
  - joins across a census without a record: 57 → 32;
  - DP control and stump-only trees: same export;
  - posterior checks pass; export among the sampled paths 98.6% → 99.4%;
  - trees solved by the DP 1,802 → 1,829 (dead ends 46 → 7); run time unchanged.
- [ ] **Residual, not changed: about 14 similar-size pairs at 2005–2010 in the whole data set.** The 46 production cases are all in the test set; 14 remain.
  - Each is a pre-2010 StemID pin against a different 2010 StemID pin (tag 198338: 89 → 90 mm).
  - Fixing them needs pre-2010 pins to be link constraints instead of track labels, in both engines and the post-engine helpers.
- [ ] Probabilistic matcher before 1995: 39 new stems of 26 cm or more on trees already measured (e.g. 004934, 005286). Not caused by the DP resprout split (3 of about 60); review the cases before choosing a fix.
- [ ] Model (optional, needs the blind test): palm growth tail (fast palm links split), link likelihood units (links favoured about 5x, deaths under-detected), mortality of large stems, strangler 3501.
- [ ] Full engine run with these settings, then the merge, stage 3 and stage 4; commit the code and documents.
- [x] Dropped 2026-10-21: a growth check in `apply_pin_track_rejoin()` (the rejoin only restores links the engine itself made).

## Validation against known identities (2026-10-22)

Before 2010 there are no tags, so the fixes were tested where identities are known: the 2010, 2015 and 2023 censuses of 4,115 multiple-stem trees, rewritten so that 2010 and 2015 look like censuses before 2010. The 2,115 test trees were those whose 2010 main stem was broken in 2015 while another stem was measured.

- One broken record per tree, on the main stem, carrying the DBH of the largest surviving stem (whose own record is removed).
- No records of dead or broken secondary stems.
- Database StemIDs renumbered at the rates the database shows before 2010: single-stem trees keep them 89–93%, main stems 8–13%, secondary stems 0–15%, broken records 64%.
- Tags only where the database carried them back.

2023 keeps its tags (the pins). Scripts in the session sandbox: `sim_build.R`, `sim_eval.R`.

Links between consecutive censuses of one stem, against the true stems (recall / precision, %):

| Configuration | All links | Break trees 2010→2015 | Break trees 2015→2023 | Other trees | DP trees |
| --- | --- | --- | --- | --- | --- |
| Database StemIDs only (Condit's identity) | 22.6 / 100 | | | | |
| Production settings (run 20261014_220452) | 59.9 / 74.2 | 20.5 / 38.2 | 44.9 / 75.3 | 83–85 / 79–85 | 65.9 / 86.4 |
| Final settings | **79.2 / 75.9** | **66.2 / 62.5** | **80.3 / 75.2** | same | **89.5 / 86.4** |
| Final, rule R1 instead of the main-stem rule | 69.1 / 74.4 | 20.4 / 37.7 | 80.6 / 75.5 | same | 77.2 / 85.2 |
| Final, whole-tree split on | 70.6 / 76.1 | 66.5 / 62.8 | 45.9 / 75.0 | same | 78.6 / 87.5 |
| Final, Step 3a pins as before | 73.3 / 73.7 | 66.5 / 62.8 | 57.3 / 65.6 | same | 82.0 / 84.0 |

Each fix carries its own part of the gain:
- The main-stem rule (instead of R1) recovers the links at the break.
- No whole-tree split and the Step 3a rule recover the links after the break (the split cuts every stem of the tree; a pin on a broken record's StemID that has no other record keeps the record off its tagged stem).

The table also shows three things:

- **The break.** The broken record's DBH is linked to its true 2010 stem in 84.5% of cases (0% with the production settings, where R1 makes it a new stem) and to its true 2023 stem in 95.4% (58.3%). No broken main stem continues into 2015 (27 with the production settings).
- **Controls.** Trees without a break are unchanged.
- **Posteriors.**
  - Link probabilities in break trees (2010→2015): mean probability of the true links 0.17 → 0.61; Brier score 0.255 → 0.120.
  - Calibration: links given 0.95–1 are true 97.5% of the time, links given 0.4–0.6 47%.
  - `DP_PosteriorReconstructedProb` is calibrated only with the label-free correction: DP rows given 0.2–0.4 are correct 44% of the time, against 87% without it, where the label split understates certainty.

Remaining limits:

- The probabilistic matcher (palm clumps, stranglers) recovers 48% of true links at 45% precision.
- In DP trees the remaining errors are mostly between stems under 2 cm.

**Real data, random sample of 2,000 multiple-stem trees (trees and shrubs, before 2010; no truth, the tagged era 2010–2023 as the yardstick).**

| Before 2010 | Database StemIDs (Condit) | Production engine | Final engine | Tagged era |
| --- | --- | --- | --- | --- |
| Stems that continue to the next census | 0% (6% in 2005–10) | 77–89% | 82–89% | 84% (2010–15) |
| Stem deaths, %/yr | every stem ends | 3.6–5.2 | 3.0–4.0 | 3.6 |
| New stems, %/yr (1990–2010) | every stem is new | 4.5–6.8 | 3.2–5.3 | 4.1–4.2 |
| Similar-size end/new pairs per 1,000 stems | 77–97 | 0–4.7 | 0–2.0 | 0 |

- With the database StemIDs alone, a stem of a multiple-stem tree never continues before 2010: each census every stem ends and a new one starts.
- The final engine's stem turnover matches the tagged era.
- Field tags in this sample (2,018, from 2010, by tag value): the production settings put 9 on two stems (tagged 2010 stumps on the stem that broke, R2 cuts) and rewrote 91 pins. The final rules for field tags (fix list above) were tested on the test set, not on this sample; the turnover figures above depend only on measured records, which those rules do not change.

**Real blind test (2010 and 2015 identities hidden, 839 trees, mostly palm clumps).** The final settings give the same result as the production settings (43.5% of true links at 41.1% precision). The fixes act only on records of the pre-2010 convention, so they leave the post-2010 reconstruction unchanged.

## Fixed rules

- `DFstatus` is never modified; it stays equal to the raw Status.
- The 2005 duplicate records are not a problem; ignore them.
- Flora de Panamá names have priority; WFO is only a spelling check.
