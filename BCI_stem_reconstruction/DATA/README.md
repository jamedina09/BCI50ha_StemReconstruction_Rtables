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

## Proposed identity fixes (from the comparison with Condit, 2026-10-21)

Details: section 7 of `2_STEM_IDENTIFICATION/comparissons/records_comparison_report.html`; candidate lists in `links_failing_growth_filter.csv` and `pins_reused_after_gap.csv` there. Each needs a sandbox engine run on the affected tags, then the full run.

- [ ] Pins across a gap (tag 200431): Step 3b of `1_main_cpp_chunk_bci.R` copies a database StemID pin to every record of the (Tag, StemID) group. Stop it at a census in which the StemID has no record. Rechecked 2026-10-21: of the 422 StemIDs measured again after a gap, 386 have a record in the gap (dead, stem dead, broken below or missing; the database followed the stem, keep them); 36 have none, all from before 2010 to 2010+, all joined in our tables, 17 change more than 2-fold.
- [ ] Clumps in the probabilistic matcher: `select_consensus_trajectory()` exports a noisy sample; export the lowest-cost assignment or repair it with local swaps (385 joins rejected by the growth filter, 82% in palms, mostly *Oenocarpus mapora*).
- [x] ~~`apply_pin_track_rejoin()`: add a growth check~~ Dropped 2026-10-21: the rejoin only restores links the engine itself made (rows of one engine track), so a check would overrule the engine; the small excess (2.3% vs 1.8% rejected without rounding) comes from harder trees.
- [ ] Whole-tree resprout split in the DP (replaces the "2010 re-tagging" item; checked 2026-10-21). Before 2010 the database gives every secondary stem a new StemID at every census and records no stem tags (all stems are tagged from 2010), so 2010 is not special. The 2005-2010 cases come from `allow_segment_split` in `match_stems_dp_global_backward_marginals_batch()`: at the first broken-below / R record with a DBH before 2010, every stem of the tree ends and every stem of the next census starts new. 15,281 multi-stem trees are split; at the split census the other stems continue 2.9% (91.5% in trees without a split); about 5,500 links are cut in trees and shrubs (3-14% of stem deaths and 1-11% of new stems per interval, 1985-2010). Fix: `allow_segment_split = FALSE`, so the DP treats the first broken-below census like the later ones and only the broken stem is split (R1), as in the matcher since 2026-09-30. Also: Step 3a pins a broken-below / R row to its own StemID even when that StemID has no other record; that pin blocks the link to the pinned 2010 stem (61 of the 121 similar-size 2005-2010 cases, e.g. tree 5096: 87 -> 94 (R) -> 98 -> 110 mm). Skip the pin in that case.
- [ ] Probabilistic matcher before 1995: 39 new stems of 26 cm or more on trees already measured (e.g. 004934, 005286). Not caused by the DP resprout split (3 of about 60); review the cases before choosing a fix.
- [ ] Decide: rule R1. Checked 2026-10-21: before 2010 a break is ONE broken-below record per tree, always on the main-stem record. Its DBH fits the tree's next-largest stem of the census before in 91-94% (the main stem alone 0.6-2%, a resprout 5-7%) and continues in the next census, so it is a surviving stem written on the main record. R1 makes it a new stem: about 5,400 false deaths plus new stems in trees and shrubs, 1985-2005. A stump-only record (34%) is a whole-tree break. Option: read the record as "the main stem broke" and let its DBH continue the stem it fits; a new stem only when it fits none.
- [ ] Model (optional, needs the blind test): palm growth tail (fast palm links split), link likelihood units (links favoured about 5x, deaths under-detected), mortality of large stems, strangler 3501.
- [ ] Before the final run: a metrics script for all of the above, a fixed test set (~5,000 tags), and a sandbox A/B (current engine vs candidate); keep only the fixes that improve it.

## Fixed rules

- `DFstatus` is never modified; it stays equal to the raw Status.
- The 2005 duplicate records are not a problem; ignore them.
- Flora de Panamá names have priority; WFO is only a spelling check.
