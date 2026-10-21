# 3_PREPARE_R_TABLES

Stage 3 of the BCI stem-reconstruction pipeline. This stage consumes the
reconstructed stem identities from `BCI_stem_reconstruction/2_STEM_IDENTIFICATION/` and generates the
ForestGEO-format census tables plus QC exports used for downstream analysis.

## Scripts

- `1_prepare_posteriors_BCI.R` — Consolidates `_paths.feather` posterior files
  from a completed stage 2 run (`home_dir` and `run_code` at the top of the
  script; the same run as in the stage-2 merge) into
  `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`,
  applying to every sample the measurement-discontinuity joins of the merge
  step (`DATA/PROCESSED/measurement_rejoin_pairs.csv`). One row per sampled
  path of a tree, with `tag`, `treeID`, `path_sig`, `path_count`, `path_prob`
  and `recon`; sample a path with weight `path_count`.
- `2_create_R_tables_BCI.R` — Builds the final census tables from
  `DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds`.
  It assigns each stem a corrected status (`Rstatus`) in every census,
  gives every tree one location, fills the `date` column, and
  exports ForestGEO-format `.Rdata` (and supporting `.csv`) tables.
- `3_create_species_table_BCI.R` — Builds the species table of the R tables,
  `DATA/RTABLES/bci.spptable.rdata`, from the stage-1 species table, with
  notes on how codes and names changed since the last public release of the
  plot data (Condit et al. 2019), and the changes of species code since then,
  one row per change and one row per tree (see [Species table](#species-table)).
  For the botanists it also lists the codes whose name or authority differs
  from the Flora de Panamá list (`DATA/CHECKS/species_changes_from_flora_panama_bci.csv`).
- `rstatus_functions.R` — The `Rstatus` and `dbh` rules as small functions,
  sourced by `2_create_R_tables_BCI.R` and by the tests.
- `species_table_functions.R` — The name, code and note rules of
  `3_create_species_table_BCI.R`, sourced by it and by the tests.

All scripts are run from the project root. The first two read the output of
the stage-2 merge (`2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R`), so
run it first. They do not depend on each other: `2_create_R_tables_BCI.R` does
not read the posterior file, which is used by
`4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R`.
`3_create_species_table_BCI.R` reads the R tables, so it runs after
`2_create_R_tables_BCI.R`.

## Rstatus and dbh rules

`Rstatus` is the corrected status; `DFstatus` keeps the raw field status.

| Rstatus | Meaning |
| --- | --- |
| `P` | before the stem's first record |
| `A` | alive: from the first record to the last evidence of life |
| `G` | stem dead, its tree alive (some stem of the tree is `A` then or later) |
| `D` | stem dead and its whole tree dead |

- Evidence of life: an `alive` record, or a DBH on a `broken below`,
  `missing` or status-less record.
- Dead record: `dead` / `stem dead` (with or without DBH), and `broken below`
  without DBH.
- Missed censuses and false deaths between a stem's first record and its last
  evidence of life are `A` (alive later ⇒ never dead). A dead first record of
  a stem alive later is `A` too.
- After the last evidence of life a stem is dead from the first census
  without one. A stem never alive is `P` until its first record, then `G`/`D`.
  A stem never recorded (only `missing` records) is `P` in every census.
  Empty placeholder rows without a stem identity (stage 1 adds one for a
  census in which the tree has no record) are dropped, so every exported stem
  has a `stemID`.
- `dbh` is exported exactly as recorded: no DBH is removed or imputed. A `G`
  or `D` cell carries a DBH only when it was recorded on a `dead` / `stem dead`
  record of a stem never alive later (listed in
  `CHECKS/dbh_on_dead_records.csv`); a `P` cell never has one; an `A` without
  DBH is a missed measurement.

## Other exported columns

- `stemID`: the reconstructed stem number within its tree, so `treeID` +
  `stemID` identify a stem. The engine numbers the stems of a tree 1, 2, ...
  in order of first appearance; in a tree changed by the measurement rejoin of
  the stage-2 merge, the joined stem keeps one of its two numbers, so the
  tree skips one number per join. Internally the script keys stems on
  `TreeID_ReconstructedStemID`; a stem without a DP identity (records with no
  status and no DBH) gets `stemID` NA.
- `ExactDate`: the field date exactly as recorded (`NA` where there is no
  record); never imputed, like `dbh` and `DFstatus`.
- `date`: the ForestGEO R-table date, in days since 1960-01-01 (the unit of
  the ForestGEO database's `Date` field and of Condit's tables;
  `as.Date(date, origin = "1960-01-01")` gives the calendar date). Every row
  has one:
  1. the recorded date (`ExactDate`) where there is one;
  2. otherwise the most common recorded date of the same tree in that census;
  3. if the tree has none, the most common recorded date of the same quadrat
     in that census;
  4. if the quadrat has none, the most common recorded date of the census.

  Only recorded dates vote, and ties go to the earliest date. The most common
  date (the mode) is used rather than a mean or median because it is always a
  day on which the field crew was recording. Condit's Dryad tables fill `date`
  with the mean recorded date of the quadrat instead: the two agree within
  about a day in 1985–2015, but in 1982 the records of one quadrat span months,
  so the tree's own date is closer. Census intervals (growth, recruitment,
  mortality) are computed from `date`.
- `gx`, `gy`, `quadrat`: one location per tree, the same for every stem and
  census. Each census casts one vote, the most common (`PX`, `PY`) pair among
  the tree's stems in that census (ties go to the smallest x, then y); the
  tree takes the pair with the most votes (ties go to the most recent census).
  `quadrat` is derived from that pair (20 m quadrats named `XXYY`); a tree
  without coordinates keeps `NA` coordinates and its most common raw quadrat.

## Species table

`3_create_species_table_BCI.R` writes one row per species code (`sp`) of the
R tables, in order of code with `uniden` last. The names and the other columns come unchanged from the stage-1
table (`DATA/SPP_TABLE/bci_spptable.RData`, `Mnemonic` renamed `sp`), plus
`Latin` = `Genus SpeciesName`, followed by the infraspecific name when there is
one (`Swartzia simplex var. grandiflora`), as in Condit's species table.

The new columns compare each code with the last public release of the plot
data: Condit R, Pérez R, Aguilar S, Lao S, Foster R, Hubbell S (2019),
*Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years*,
Dryad, doi:10.15146/5xcp-0d46. Only codes and names (genus, epithet,
infraspecific name) are compared, not authority, family or other columns.

| Column | Content |
| --- | --- |
| `NameCondit2019` | the name of the code in Condit et al. (2019), when it differs from `Latin` |
| `NameChange` | kind of change: `genus`, `genus and epithet`, `genus and epithet gender`, `epithet`, `epithet gender` (only the gender ending, e.g. *micrantha* → *micranthum*), `epithet spelling` (at most 2 other letters), `infraspecific name` or `morphospecies label`. It describes the change; it does not say which spelling is correct |
| `ReidentifiedFrom` | the codes that the trees now under this code had in Condit et al. (2019), with the number of trees, e.g. `tremsp (84)` |
| `ReidentifiedTo` | the codes that the trees of this code in Condit et al. (2019) have now |
| `Notes` | only what no column says, so that each fact is written once: the note of the stage-1 table, then that `uniden` is not a taxon, that the code was not used in Condit et al. (2019), and the codes no longer used whose trees are now under this code, with their former name, e.g. "Code tremsp, Trema unidentified in Condit et al. (2019), is no longer used." Former names are only in `NameCondit2019` and `NameChange`, and numbers of trees only in `ReidentifiedFrom` and `ReidentifiedTo` |

The releases are matched tree by tree on `treeID`, which pairs the same trees
as matching on `tag`. A code no longer used (on 2026-10-15: tremsp, ingaun,
pterof) has no row: the code that holds its trees names it in `Notes`, and
counts the trees in `ReidentifiedFrom`. `uniden` (unidentified stems) gets a row without names, the
last row (`INCLUDE_UNIDEN` at the top of the script), so that every code of the R
tables has a row: `merge()` of a census table with the species table is an
inner join by default and would drop the uniden stems otherwise.

The changes of species code since Condit et al. (2019) are in two files:

- `bci.species_code_changes_summary`: one row per change (former code → code
  now), with `spCondit2019`, `LatinCondit2019`, `sp`, `Latin` and `trees`, the
  number of trees, largest changes first;
- `bci.species_code_changes`: one row per tree, with `treeID`, `tag`,
  `spCondit2019`, `LatinCondit2019`, `sp` and `Latin`. Tags are text with their
  leading zeros; read the `.csv` with the tag column as text.

### Changes from the Flora de Panamá list

The stage-1 names come from the Flora de Panamá list (R. Pérez, 2026-10-09,
`DATA/RAW/sp_tables/flora_panama_oct_9_2026/`) and change only by the hand
corrections of `1_DATA_PREPARATION/1_prepare_species_tables.R` (section 9).
For the botanists to review those corrections, the script writes
`DATA/CHECKS/species_changes_from_flora_panama_bci.csv`: one row per code
whose genus, epithet, infraspecific name or authority differs from the Flora
list, with `sp`, `FloraName` and `Name` (the full names with authority),
`Changed` (the parts that differ: `genus`, `epithet`, `infraspecific name`,
`authority`) and the parts of both names. The authority is compared by its
authors only: another way of writing the same authors (spaces, periods,
capitals, accents: `C. DC.` = `C.DC.`) is not a change. Family, synonyms and
the other columns are not compared. Codes not in the Flora list (the
morphospecies apeihy, nects1, nects3) have no baseline and are named in the
log. On 2026-10-21 the file lists 10 codes, all by the authority (anacex,
ficuc2, heisco, ouralu, pse2sp, tab1ro, thevah, tri2tu, tropra, zantpr): the
WFO authorities adopted in stage 1. They are the same codes as in the
stage-1 file `CHECKS/species_wfo_vs_flora_panama.csv` (WFO against the Flora
de Panamá). It is a review file and is not part of the R tables.

The four `.csv` files are UTF-8 with a byte order mark, so that Excel shows
the Spanish characters (`árbol`, `Pérez`) correctly; `fread()` reads them as
they are, and `read.csv()` needs `fileEncoding = "UTF-8-BOM"`.

## Tests

`tests/test_Rstatus.R` (testthat) checks the rules with hand-written
fixtures, every combination of raw states for small trees against a simple
reference implementation (`tests/reference_rstatus.R`), the validity of every
`Rstatus` sequence (`tests/rstatus_validity.R`), and the exported R tables
cell by cell (`Rstatus`, `dbh`, `DFstatus`, `ExactDate`, `date`, location and
`stemID`). Run from the project root:

```r
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R")
```

`RSTATUS_TEST_LEVEL=quick` runs the smaller enumerations only (a few minutes
with the real-data tests; `full`, the default, takes about 10 minutes on 16
cores). `RSTATUS_RTABLES_DIR`
and `RSTATUS_STAGE2_FILE` point the real-data tests at other tables (default:
`DATA/RTABLES` and `DATA/PROCESSED`); if the tables or the stage-2 file are
missing, the real-data tests are skipped. `RSTATUS_CORES` sets the cores of
the reference implementation (default: all but two).

`tests/test_species_table.R` (testthat) checks the functions of
`species_table_functions.R` on hand-written cases, recomputes the species
outputs from the R tables and Condit's files with separate code, and the
changes from the Flora de Panamá list from the stage-1 table and the Flora
list (skipped when the outputs are missing; `SPECIES_RTABLES_DIR` and
`SPECIES_CHECKS_DIR` point it at other folders, default `DATA/RTABLES` and
`DATA/CHECKS`):

```r
testthat::test_file("BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_species_table.R")
```

## Requirements

R packages: `data.table`, `arrow` and `inspectdf`; `igraph` (transition
diagram in `2_create_R_tables_BCI.R`); `digest` and `readxl` (`3_create_species_table_BCI.R`);
`testthat` and `parallel` for the tests.

## Outputs

- `BCI_stem_reconstruction/DATA/RTABLES/<site>.stemN.Rdata` (and `<site>.stemN.csv`),
  one table per census with the columns `CensusID`, `treeID`, `stemID`, `tag`,
  `StemTag`, `sp`, `quadrat`, `gx`, `gy`, `dbh` (mm), `hom` (m), `ExactDate`,
  `date`, `codes`, `DFstatus`, `Rstatus` and `StemPaths` (the row key of the
  posterior paths); every table lists the same stems in the same order
- `BCI_stem_reconstruction/DATA/RTABLES/bci.spptable.rdata` (object
  `bci.spptable`, a data.frame) and `bci.spptable.csv`: the species table
- `BCI_stem_reconstruction/DATA/RTABLES/bci.species_code_changes_summary.rdata`
  and `.csv`: the changes of species code since Condit et al. (2019), one row
  per change with its number of trees
- `BCI_stem_reconstruction/DATA/RTABLES/bci.species_code_changes.rdata` and
  `.csv`: the trees whose species code changed since Condit et al. (2019)
- `BCI_stem_reconstruction/DATA/POSTERIORS/posterior_sampled_paths.rds`
- QC exports in `BCI_stem_reconstruction/DATA/CHECKS/`, among them
  `dbh_on_dead_records.csv`, `location_conflicts.csv`,
  `subset_stems_never_alive.csv`, `subset_stems_never_recorded.csv` and
  `duplicate_measurements.csv` (likely duplicates: two stems of one tree
  measured on the same date at heights at least 0.5 m apart, taper-corrected
  DBHs within 20 %; the trunk measured at the old and the new height when its
  point of measurement was raised), `species_changes_from_flora_panama_bci.csv`
  (the codes whose name or authority differs from the Flora de Panamá list,
  for the botanists) and `check_species_table_bci.txt` (failed checks and
  notes of `3_create_species_table_BCI.R`)
