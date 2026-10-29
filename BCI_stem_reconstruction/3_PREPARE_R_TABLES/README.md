# 3_PREPARE_R_TABLES

Stage 3 of the pipeline. It turns the stem table of stage 2 into the ForestGEO
R tables (one table per census), the species table and the review files, and
collects the posterior samples used by stage 4.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Inputs](#inputs)
- [Outputs](#outputs)
  - [R-table columns](#r-table-columns)
- [What it does](#what-it-does)
  - [`1_prepare_posteriors_BCI.R`](#1_prepare_posteriors_bcir)
  - [`2_create_R_tables_BCI.R`](#2_create_r_tables_bcir)
  - [`3_create_species_table_BCI.R`](#3_create_species_table_bcir)
- [Checks](#checks)
- [Settings](#settings)
- [Tests](#tests)
- [Requirements](#requirements)

## Contents

| File | What it does |
| --- | --- |
| `1_prepare_posteriors_BCI.R` | Collects the posterior path files of the engine run into one table |
| `2_create_R_tables_BCI.R` | Assigns the corrected status (`Rstatus`), one location per tree and a date to every row, and writes the R tables and the check files |
| `3_create_species_table_BCI.R` | Writes the species table of the R tables and the files of the species code changes |
| `rstatus_functions.R` | The `Rstatus` and `dbh` rules, and `tree_state()`, an independent check of the tree status; sourced by script 2 and the tests |
| `species_table_functions.R` | The name, code and note rules of script 3; sourced by script 3 and the tests |
| `tests/` | Tests of the `Rstatus` rules and of the species table |

## How to run

From the project root, after the stage-2 merge. Set `home_dir` and `run_code`
at the top of `1_prepare_posteriors_BCI.R` to the engine run (the same run as
in the merge).

```sh
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/1_prepare_posteriors_BCI.R
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/2_create_R_tables_BCI.R
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/3_create_species_table_BCI.R
```

Scripts 1 and 2 are independent: script 2 does not read the posterior file,
which is used by `4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R`.
Script 3 reads the R tables, so it runs after script 2.

## Inputs

| File | Content |
| --- | --- |
| `<home_dir>/<run_code>/posteriors/*_paths.feather` | Posterior path files of the engine run (script 1) |
| `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds` | Stem table of stage 2 (scripts 1 and 2) |
| `DATA/PROCESSED/measurement_rejoin_pairs.csv` | Observation pairs joined by the stage-2 measurement rejoin (script 1) |
| `DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds` | Census table of stage 1, for the `treeID` of each tag (script 1) |
| `DATA/SPP_TABLE/bci_spptable.RData` | Species table of stage 1 (script 3) |
| `DATA/SPP_TABLE/condit/bci.spptable.rdata`, `DATA/RAW/bci_dryad_condit.rds` | Species table and stem data of Condit et al. (2019) (script 3) |
| `DATA/RAW/sp_tables/flora_panama_oct_9_2026/` | Flora de Panamá list (script 3) |

## Outputs

| File | Content |
| --- | --- |
| `DATA/POSTERIORS/posterior_sampled_paths.rds` | One row per sampled path of a tree: `tag`, `treeID`, `path_sig`, `path_count`, `path_prob`, `recon`; attribute `engine_run` (the run code). Draw a tree's paths with weight `path_count` |
| `DATA/RTABLES/bci.stem1.Rdata` … `bci.stem9.Rdata` (and `.csv`) | One R table per census (columns below) |
| `DATA/RTABLES/bci.spptable.rdata` (object `bci.spptable`, a data.frame) and `.csv` | Species table: one row per species code of the R tables |
| `DATA/RTABLES/bci.species_code_changes_summary.rdata` and `.csv` | One row per change of species code since Condit et al. (2019), with its number of trees |
| `DATA/RTABLES/bci.species_code_changes.rdata` and `.csv` | One row per tree whose species code changed since Condit et al. (2019) |
| `DATA/CHECKS/check_posteriors_bci.txt` | Every check of script 1 |
| `DATA/CHECKS/check_bci.txt` | Every check of script 2 and its summary tables; read it first |
| `DATA/CHECKS/status_patterns_bci.csv` | One row per encounter history (raw codes, `Rstatus` and tree histories, number of stems, one example stem) |
| `DATA/CHECKS/field_records_to_review_bci.csv` | Raw records that need a look, for the data manager (one row per case, column `reason`) |
| `DATA/CHECKS/reconstruction_to_review_bci.csv` | Stems whose reconstruction needs a look (one row per case, column `reason`) |
| `DATA/CHECKS/check_species_table_bci.txt` | Checks and notes of script 3 |
| `DATA/CHECKS/species_changes_from_flora_panama_bci.csv` | Codes whose name or authority differs from the Flora de Panamá list, for the botanists |

The `.csv` files of the species table and of the code changes are UTF-8 with a
byte order mark: `fread()` reads them as they are, `read.csv()` needs
`fileEncoding = "UTF-8-BOM"`.

### R-table columns

| Column | Type | Content |
| --- | --- | --- |
| `CensusID` | number | Census 1–9 |
| `treeID` | text | Tree |
| `stemID` | integer | Reconstructed stem within its tree; `treeID` + `stemID` identify a stem |
| `tag`, `StemTag` | text | Tree tag and stem tag (leading zeros and signs kept, e.g. `000006`, `-05599`) |
| `sp` | text | Species code |
| `quadrat`, `gx`, `gy` | text, number, number | One location per tree |
| `dbh` | number | DBH in mm, as recorded |
| `hom` | number | Height of measurement in m |
| `ExactDate` | date | Field date as recorded (`NA` without a record) |
| `date` | number | Date in days since 1960-01-01, on every row |
| `codes` | text | Field codes |
| `DFstatus` | text | Raw field status |
| `Rstatus` | text | Corrected status (`P`, `A`, `G`, `D`) |
| `StemPaths` | integer | Row key of the posterior paths (`ObsRowID`) |

Every table lists the same stems in the same order, sorted by `treeID` and
then `stemID` as numbers (tree 20 before tree 100). Read the `.csv` with the
text columns as text.

## What it does

### `1_prepare_posteriors_BCI.R`

1. Reads every posterior path file of the run and adds the `treeID` of each
   tag.
2. Applies the measurement-rejoin joins of the merge to every sample in which
   the two stems are split with a clean end and start.
3. Writes `posterior_sampled_paths.rds`.

### `2_create_R_tables_BCI.R`

1. Gives every stem a corrected status in every census:

   | Rstatus | Meaning |
   | --- | --- |
   | `P` | before the stem's first record |
   | `A` | alive: from the first record to the last evidence of life |
   | `G` | dead, while its tree is alive (another stem of the tree is `A` then or later) |
   | `D` | dead, and its whole tree dead |

   - Evidence of life: an `alive` record; a `broken below` record (with a DBH,
     a measured resprout; without one, a stump: the stem broke below 1.3 m and
     the field crew recorded it alive); or a DBH on a `missing` or status-less
     record. Dead record: `dead` or `stem dead` (with or without DBH).
   - Between a stem's first record and its last evidence of life the stem is
     `A`, also in a census without a record or with a dead record (alive later
     means never dead); a dead first record of a stem alive later is `A` too.
     After its last evidence of life it is dead from the first census without
     one.
   - A stem never alive is `P` until its first record, then `G` or `D`. A stem
     with only `missing` records is `P` in every census.
   - Placeholder rows of stage 1 that carry no stem are dropped.
2. Exports `dbh` and `ExactDate` as recorded, never removed or imputed, and
   `DFstatus` equal to the raw status. A `G` or `D` cell carries a DBH only when
   one was recorded on a dead record of a stem never alive later. An `A` cell
   without a DBH is a missed measurement, or a stump: a `broken below` record
   without a DBH (alive, with no stem at 1.3 m; stage 4 keeps it out of the
   stock), or a dead, missing or absent record just before a stump of the
   same stem (alive because the stump comes later; no stem at 1.3 m is known).
   A tree whose only `A` stems are stumps is alive, broken below 1.3 m.
3. Numbers the stems of each tree (`stemID`) as the engine does (1, 2, … in
   order of first appearance), so `treeID` + `stemID` identify a stem. A stem
   joined by the measurement rejoin keeps one of its two numbers, so such a tree
   skips one number per join. Internally stems are keyed on
   `TreeID_ReconstructedStemID`; every exported stem has a `stemID`.
4. Fills `date` on every row, in days since 1960-01-01 (the unit of the
   ForestGEO database's `Date` field; `as.Date(date, origin = "1960-01-01")`
   gives the calendar date):
   1. the recorded date (`ExactDate`) where there is one;
   2. otherwise the most common recorded date of the same tree in that census;
   3. otherwise the most common recorded date of its quadrat in that census;
   4. otherwise the most common recorded date of the census.

   Only recorded dates vote, and ties go to the earliest date. The most common
   date (the mode) is used rather than a mean or median because it is always a
   day on which the field crew was recording; the tree's own date comes first
   because the records of one quadrat can span months (1982). Census intervals
   (growth, recruitment, mortality) are computed from `date`.
5. Gives each tree one location: each census votes for the most common
   (`PX`, `PY`) pair of the tree's stems (ties: smallest x, then y), and the
   tree takes the pair with the most votes (ties: the most recent census).
   `quadrat` (20 m, named `XXYY`) comes from that pair; a tree without
   coordinates keeps `NA` and its most common raw quadrat.
6. Checks the order and the column types of every table, then writes the R
   tables and the check files.

### `3_create_species_table_BCI.R`

1. Writes one row per species code (`sp`) of the R tables, in code order with
   `uniden` last. The names and the other columns come unchanged from the
   stage-1 table (`Mnemonic` renamed `sp`), plus `Latin` = `Genus SpeciesName`,
   followed by the infraspecific name when there is one (`Swartzia simplex
   var. grandiflora`), as in Condit's species table. `uniden` (unidentified
   stems) gets a row without names (`INCLUDE_UNIDEN = TRUE`), so that every
   code of the R tables has a row: `merge()` of a census table with the species
   table is an inner join by default and would drop the unidentified stems
   otherwise.
2. Compares each code with the release of Condit R, Pérez R, Aguilar S, Lao S,
   Foster R, Hubbell S (2019), *Complete data from the Barro Colorado 50-ha
   plot: 423617 trees, 35 years*, Dryad, doi:10.15146/5xcp-0d46 (the last
   public release; its species table is checked by SHA-256). The releases are
   matched tree by tree on `treeID`, which pairs the same trees as matching on
   `tag`. Only codes and names (genus, epithet, infraspecific name) are
   compared, not authority, family or other columns.

   | Column | Content |
   | --- | --- |
   | `NameCondit2019` | Name of the code in Condit et al. (2019), when it differs from `Latin` |
   | `NameChange` | Kind of change: `genus`, `genus and epithet`, `genus and epithet gender`, `epithet`, `epithet gender` (only the gender ending, e.g. *micrantha* → *micranthum*), `epithet spelling` (at most 2 other letters, edit distance), `infraspecific name` or `morphospecies label`. It describes the change, not which spelling is correct |
   | `Notes` | Only what no column says, so each fact is written once: the note of the stage-1 table; that `uniden` is not a taxon; that the code was not used in Condit et al. (2019); and the codes no longer used whose trees are now under this code, with their former name and number of trees, e.g. "Code tremsp, Trema unidentified in Condit et al. (2019), is no longer used; its 84 trees are under this code." |

   A code no longer used has no row: the code that holds its trees names it in
   `Notes`. Trees re-identified from one code to another are not in the species
   table but in the two files of step 3.

3. Writes the changes of species code: `bci.species_code_changes_summary`
   (one row per change: `spCondit2019`, `LatinCondit2019`, `sp`, `Latin`,
   `trees`, largest changes first) and `bci.species_code_changes` (one row per
   tree: `treeID`, `tag`, `spCondit2019`, `LatinCondit2019`, `sp`, `Latin`;
   tags are text with their leading zeros, so read the `.csv` with the tag
   column as text).
4. Writes `species_changes_from_flora_panama_bci.csv`, for the botanists. The
   stage-1 names come from the Flora de Panamá list and change only by the hand
   corrections of `1_DATA_PREPARATION/1_prepare_species_tables.R` (section 9);
   this file lists them: one row per code whose genus, epithet, infraspecific
   name or authority differs from the Flora list, with `sp`, `FloraName` and
   `Name` (full names with authority), `Changed` (the parts that differ) and
   the parts of both names. The authority is compared by its authors only:
   another way of writing the same authors (spaces, periods, capitals, accents:
   `C. DC.` = `C.DC.`) is not a change. Family, synonyms and the other columns
   are not compared. The morphospecies `apeihy`, `nects1` and `nects3` are not
   in the Flora list (no baseline) and are named in the log. The stage-1 file
   `CHECKS/species_wfo_vs_flora_panama.csv` lists, from the WFO side, the WFO
   corrections adopted against the Flora list. It is a review file, not part of
   the R tables.

## Checks

A failed check prints ❌, raises a warning and writes to the script's check
file before it stops.

Each script first checks that it reads the right inputs:

- script 1: the run `home_dir/run_code` finished, `run_code` is the run of the
  merge (`DATA/PROCESSED/engine_run.txt`), every posterior tag is a
  multiple-stem tag of the final table, and every observation of the paths is
  a measured observation of its tag there;
- script 2: the stage-2 table is newer than the stage-1 table (after rerunning
  stage 1, rerun the engine and the merge) and `engine_run.txt` exists;
- script 3: the stage-1 census table, the stage-2 table and the R tables were
  written in this order. The stage-1 species table is not dated: its names can
  be corrected and script 3 rerun without rebuilding the other stages.

Script 2 checks the order and column types of every
table before writing it, and writes every check to `check_bci.txt` (✓ passed,
❌ failed) with the summary tables: raw records alive because of their DBH,
status codes, the biology assessment of every version of the histories, the
evidence behind each stem death, status × DBH by census (stumps apart from
missed measurements), and at the end the number of cases in each review list.
Read it first.

`status_patterns_bci.csv` has one row per encounter history: `raw_history`
(raw codes, one per census: `A` evidence of life, `D` dead record, `N` no
record), `dbh_recorded` (1 where a DBH was recorded), `Rstatus_history`,
`tree_history`, the number of stems and one example stem (`example_treeID`,
`example_stemID`).

The two review lists have one row per case with `reason`, `treeID`, `tag`,
`stemID` (NA for a whole tree), `census`, `dbh`, `DFstatus`, `Rstatus` and
`details`:

| File | Reasons |
| --- | --- |
| `field_records_to_review_bci.csv` (raw records, for the data manager) | `dead_record_of_stem_alive_later` (a dead / stem dead record of a stem recorded alive later, exported `A`), `dbh_on_dead_record`, `stem_dead_but_tree_dead` (raw "stem dead" while no stem of the tree is alive then or later), `duplicate_measurement` (two stems of one tree measured on the same date at heights at least 0.5 m apart, taper-corrected DBHs within 20%: the trunk measured at the old and the new height when its point of measurement was raised), `two_dates_in_census`, `location_spread_over_1m`, `no_coordinates`, `quadrat_differs_from_coordinates` |
| `reconstruction_to_review_bci.csv` (stems) | `identity_break_candidate` (one stem ends and a new stem of similar size starts in the same tree and census), `new_stem_50mm_or_more` (ingrowth or a new tree first measured at 50 mm or more), `stem_never_alive`, `stem_never_recorded`, `stump_of_stem_measured_later` (a stump, exported `A` without a DBH, on a stem recorded at 1.3 m in a later census) |

Script 3 writes its failed checks and notes to `check_species_table_bci.txt`.
The stage-1 species script writes the WFO comparison lists
(`species_*wfo*.csv`, `species_table_vs_wfo.csv`).

## Settings

| Setting | Script | Default | What it does |
| --- | --- | --- | --- |
| `home_dir`, `run_code` | 1 | engine output folder; run folder name | The engine run whose posterior files are read (the same run as in the stage-2 merge) |
| `INCLUDE_UNIDEN` | 3 | `TRUE` | Add the row of `uniden` to the species table |

## Tests

```r
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R")
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_species_table.R")
```

`test_Rstatus.R` checks:

- the rules on hand-written cases;
- every one of the 12 raw states (each `Status` with and without a DBH, and no
  row) against its code (evidence of life, dead record, no record);
- every stem history over the 9 censuses (3^9 code sequences) in every tree
  context against the rules written from scratch (`tests/reference_rstatus.R`).
  `Rstatus` depends only on the stem's codes and the last census in which its
  tree is alive, so this covers every stem;
- the validity of every `Rstatus` sequence (`tests/rstatus_validity.R`);
- `tree_state()` (the script's check I4) against the tree rule in every tree of
  1–4 stems built from the 55 stem histories the rules can produce; three
  stems already give every combination a census can hold, so this covers trees
  of any size;
- the exported R tables cell by cell (`Rstatus`, `dbh`, `DFstatus`,
  `ExactDate`, `date`, location, `stemID`, row order, column types, raw state ×
  `Rstatus`).

Options (environment variables): `RSTATUS_TEST_LEVEL` (`full`, the default,
about 10 minutes on 16 cores; `quick`, the smaller enumerations only, a few
minutes with the real-data tests),
`RSTATUS_RTABLES_DIR` and `RSTATUS_STAGE2_FILE` (other tables to test),
`RSTATUS_CORES` (default: all cores but two). The real-data tests are skipped
when the tables are missing.

`test_species_table.R` checks the functions of `species_table_functions.R` on
hand-written cases, and recomputes the species outputs from the R tables and
Condit's files, and the Flora de Panamá changes from the stage-1 table, with
separate code. `SPECIES_RTABLES_DIR` and `SPECIES_CHECKS_DIR` point it at other
folders; it skips the real-data tests when the outputs are missing.

## Requirements

R packages: `data.table`, `arrow`, `inspectdf`; `digest` and `readxl`
(script 3); `testthat` and `parallel` (tests).
