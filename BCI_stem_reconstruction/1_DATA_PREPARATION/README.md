# 1_DATA_PREPARATION

Stage 1 of the pipeline. It builds the species table from the botanists' lists
and the census table that the engine reads in stage 2
(`2_STEM_IDENTIFICATION/`). It reconstructs no stems.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [What it does](#what-it-does)
  - [`1_prepare_species_tables.R`](#1_prepare_species_tablesr)
  - [`2_prepare_viewfulltable.R`](#2_prepare_viewfulltabler)
- [Checks](#checks)
- [Requirements](#requirements)

## Contents

| File | What it does |
| --- | --- |
| `1_prepare_species_tables.R` | Builds the species table `bci_spptable` and the WFO check lists |
| `2_prepare_viewfulltable.R` | Builds the census table: one measurement per stem and census, taper-corrected DBH, single- and multiple-stem tags labelled |

## How to run

From the project root, in this order (script 2 reads the species table of
script 1):

```sh
Rscript BCI_stem_reconstruction/1_DATA_PREPARATION/1_prepare_species_tables.R
Rscript BCI_stem_reconstruction/1_DATA_PREPARATION/2_prepare_viewfulltable.R
```

Both scripts are written to be run line by line: besides their outputs they
print inspection tables, and script 2 draws diagnostic plots.

## Inputs

All in `BCI_stem_reconstruction/DATA/RAW/`.

| File | Content |
| --- | --- |
| `oct_13_2026/ViewFullTable_bci.csv` | Census records, database export of 2026-10-13 (tab-separated). Script 1 reads only the species codes and names (`Mnemonic`, `Family`, `Genus`, `SpeciesName`) |
| `ViewFiles_bci_allcensuses/ViewFullTable_bci.csv` | An earlier export. Script 2 reads only its `DBHID` column, to check that both exports hold the same records |
| `sp_tables/email_1/Lista_bci_mnemonics_formadevida.xlsx` | Species list of R. Pérez and S. Aguilar, with growth forms, common names and herbarium collections |
| `sp_tables/email_2/bci_name_changes.xlsx` | Name changes of S. Lao and the reply of R. Pérez |
| `sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx` | Flora de Panamá list (R. Pérez, 2026-10-09) |
| `sp_tables/wfo_plant_list_2026_06/classification.csv` | World Flora Online Plant List, June 2026 (v.2026-06, CC0): file `_DwC_backbone_R.zip` of doi:10.5281/zenodo.20782718, unzipped |

## Outputs

| File | Content |
| --- | --- |
| `DATA/SPP_TABLE/bci_spptable.txt`, `.csv` (tab-separated) and `.RData` (object `bci.spptable`) | One row per species code of the plot, `uniden` (unidentified stems) aside. Columns: `Mnemonic`, `Order`, `Family`, `Genus`, `SpeciesName`, `InfraspecificRank`, `InfraspecificEpithet`, `Authority`, `Synonyms`, `Lifeform_RFoster`, `Lifeform_RPerez_SAguilar`, `CommonName`, `Herbarium`, `Notes` |
| `DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds` | The census table: the `ViewFullTable` columns plus `dbh_with_best_candidate_taper_corrected` (mm), `single_stem_tags`, `Lifeform` and `RowID` |
| `DATA/CHECKS/species_wfo_review.csv` | Each name next to its WFO record |
| `DATA/CHECKS/species_wfo_candidates.csv` | Names with several WFO records |
| `DATA/CHECKS/species_wfo_accepted_changes.csv` | The WFO corrections adopted in section 9 of script 1 |
| `DATA/CHECKS/species_wfo_vs_flora_panama.csv` | Names whose WFO record differs from the Flora de Panamá in family, genus, epithet, infraspecific name or authority, with both versions of each part |
| `DATA/CHECKS/species_table_vs_wfo.csv` | Names of the final species table that still differ from their WFO record (family, genus, epithet, infraspecific name, authority, or only the writing of the authority), or that WFO treats as synonyms, with the accepted name |
| `DATA/CHECKS/species_synonyms_vs_wfo.csv` | The `Synonyms` column checked against WFO |

## What it does

### `1_prepare_species_tables.R`

1. Takes the species codes (`Mnemonic`) from the database export. The
   export's species names are not used, except for the three morphospecies
   that no list names (step 3).
2. Names every code from the botanists' lists: email 1, then the name changes
   of email 2, then the Flora de Panamá list, which has priority.
3. Names the morphospecies Apeiba hybrid (`apeihy`), Nectandra fuzzy
   (`nects1`) and Nectandra sp.4_(tiny_leaf) (`nects3`) from the export, and
   checks them against email 1 (Apeiba "hybrida", Nectandra s1 and s3): the
   genus and the family must be those of the email, and the family that of
   the other species of the genus in the list.
4. Replaces the code `pterof` by `pterro`: every Pterocarpus of the plot is
   *P. rohrii*.
5. Checks every name against the WFO backbone with the `WorldFlora` package.
   WFO only suggests corrections; each correction is adopted by hand in
   section 9. The authority is compared by its authors only (spaces, periods,
   capitals and accents aside).
6. Writes the species table and the check lists.

### `2_prepare_viewfulltable.R`

1. Loads the census records (DBH in mm, HOM in m), removes the partial census
   of 2020 (export census 9: 135 of the 1,250 quadrats), replaces `pterof` by
   `pterro` and numbers the censuses 1 to 9 by date.
2. Adds a placeholder row (no stem, no DBH) for every census between a tag's
   first and last record in which the tag has no row. Stage 3 drops the
   placeholders that stay without a stem.
3. Keeps one measurement per stem and census: the highest HOM, ties going to
   the larger DBH.
4. Flags likely DBH entry errors (a jump in and back out over three
   consecutive measurements) and replaces them, in a separate column, by the
   geometric mean of their neighbours. The recorded `DBH` is kept.
5. Corrects the DBH to 1.3 m with the taper model of Cushman et al. (2014)
   (`dbh_with_best_candidate_taper_corrected`) for species whose growth form
   contains "árbol" or "arbusto" (trees and shrubs, strangler figs included)
   or that have none. Palms and tree ferns keep the recorded value.
6. Labels each tag as single-stem (one `StemID`, no `StemTag`) or
   multiple-stem. Stage 2 reconstructs the multiple-stem tags only.
7. Takes only the code and the growth form (`Lifeform_RPerez_SAguilar`) from
   the species table, so species names cannot change the stem identification.

## Checks

A failed check prints ❌, raises a warning and writes to the script's log in
`DATA/CHECKS/` before the script stops. Script 2 checks, among others, that
both exports hold the same `DBHID`s.

## Requirements

R packages: `data.table`, `WorldFlora`, `stringr`, `janitor`, `readxl`
(script 1); `data.table`, `ggplot2`, `inspectdf` (script 2).
