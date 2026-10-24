# DATA

Inputs and outputs of the pipeline. `RAW/` holds the inputs; every other folder
is written by a stage and is rebuilt by running stages 1–4
(`BCI_stem_reconstruction/README.md`).

## Contents

| Folder | Content | Written by |
| --- | --- | --- |
| `RAW/` | Inputs: database exports (`oct_13_2026/`, `oct_16_2026/`, `ViewFiles_bci_allcensuses/`), species sources (`sp_tables/`), the stem data of Condit et al. (2019) (`bci_dryad_condit.rds`) | — |
| `SPP_TABLE/` | Stage-1 species table (`bci_spptable.*`); Condit's species table (`condit/bci.spptable.rdata`, input); the helper `compare_species_condit_and_new_release.R` and its outputs in `condit/` | stage 1, helper |
| `PROCESSED/` | Stage-1 census table (`ViewFullTable_single_vs_multiple_stem_tags.rds`, the engine input); stage-2 final stem table (`complete_dataset_final_with_reconstructed_stemids.rds`) and the measurement-rejoin files | stages 1 and 2 |
| `<run_code>/` | Merged engine output (`merged_output.parquet`, `.rds`); the engine run folder itself is in `/Users/medinaja/outputs_bci_stem_identification/` | stage 2 |
| `POSTERIORS/` | Posterior samples of every tree (`posterior_sampled_paths.rds`) | stage 3 |
| `RTABLES/` | R tables `bci.stem1`–`bci.stem9`, species table and species code changes | stage 3 |
| `CHECKS/` | Check and review files: `check_bci.txt` (read it first), `status_patterns_bci.csv`, `field_records_to_review_bci.csv`, `reconstruction_to_review_bci.csv`, `check_species_table_bci.txt`, `species_changes_from_flora_panama_bci.csv` (stage 3) and the WFO lists (stage 1) | stages 1 and 3 |
| `BACKUPS/` | Archives of earlier outputs, each with its original paths: `before_final_run_20261024/` (with `README.txt` and an md5 manifest), `species_tables_20261021/`, and zip archives | — |

## Inputs

| File | Content |
| --- | --- |
| `RAW/oct_13_2026/ViewFullTable_bci.csv` | Census records and species codes: the database export read by the pipeline |
| `RAW/oct_16_2026/` | A later export with the same records; not read |
| `RAW/ViewFiles_bci_allcensuses/ViewFullTable_bci.csv` | An earlier export; stage 1 reads its `DBHID` column to check the records |
| `RAW/sp_tables/` | Species lists of the botanists (`email_1/`, `email_2/`), the Flora de Panamá list (`flora_panama_oct_9_2026/`) and the WFO backbone (`wfo_plant_list_2026_06/`) |
| `RAW/bci_dryad_condit.rds` | Stem data of Condit et al. (2019), for the species code changes of stage 3 |

## Status

The output folders are empty; their earlier contents are in
`BACKUPS/before_final_run_20261024/`. Running stages 1–4 rebuilds them.

## Known data issues

How the pipeline treats the issues found in the database export. `DFstatus`,
`dbh` and `ExactDate` are always exported as recorded.

Database export:

| Issue | Treatment |
| --- | --- |
| Species names, families and taxonomy rows of the export (swars1 / swars2 variety labels, duplicate taxonomy rows of protte and quaras, old or mixed families of 48 codes) | Not read: codes come from the export, names from the botanists' lists |
| Stems dead and then alive (1,123 stems, 1,267 trees); dead records with a DBH (34); alive records without a DBH (557) | Rstatus rules: alive later means alive; a dead record with a DBH is listed for review |
| Probable DBH typos (174; tree 1114: 2375 mm in 1985); stems losing more than half their DBH (189) | Kept as recorded |
| Two measurements at the same HOM and date (6 stems); two `StemID`s for one `StemTag` (trees 29078, 315977, 155532); tree 236797 stem 1a dead beside live stems | Kept as recorded |
| Records without a date (177); dates outside their census (98; tree 73712 stem 5 dated 2013) | `ExactDate` kept; stage 3 fills `date` |
| Trees without coordinates (86); trees with two locations (4); quadrat that differs from the coordinates (tree 550132) | One location per tree (stage 3); trees without coordinates keep `NA` |
| `StemNumber` (0 or empty from 2015), `LargeStem` (empty), HOM on records without a DBH | Not used |
| `PrimaryStem` (empty from 2015) | Used up to 2010, where a break is written on the main stem's record |
| Partial census of 2020 (135 of 1,250 quadrats) | Removed in stage 1 |

Species table:

| Issue | Treatment |
| --- | --- |
| Spellings that differ from IPNI: appuse *siebertii* (IPNI *seibertii*), ardiba *bartletii*, rinosy *silvatica* Kuntze, malmsp *garwoodii*, pipeco *trinervium* | Flora de Panamá spellings kept |
| Names that WFO treats as synonyms: guargr, laetth, swars1, swars2 | Flora de Panamá names kept (`CHECKS/species_table_vs_wfo.csv`) |
| Synonym column against WFO | Listed in `CHECKS/species_synonyms_vs_wfo.csv` |
| Authorities of thevah and anacex | WFO / IPNI authorities adopted (section 9e of `1_prepare_species_tables.R`) |
| Infraspecific ranks of protte (subsp. *sessiliflorum*) and myrcga ("tip. *gatunensis*") | As in the Flora de Panamá list |
| Morphospecies apeihy, nects1, nects3, sapisp | Database labels, checked against email 1; growth form "árbol" |
| Empty common names (19 codes) and herbarium fields (76 codes) | Left empty |

## Fixed rules

- `DFstatus` is never modified; it equals the raw `Status`.
- The 2005 records that look like duplicates are kept as recorded.
- Flora de Panamá names have priority; WFO is only a spelling check.
