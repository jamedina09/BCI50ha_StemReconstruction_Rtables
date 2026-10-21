# 1_DATA_PREPARATION

This folder contains the first stage of the BCI stem reconstruction pipeline.
It turns the raw ForestGEO export and the botanists' species lists into the
species table and the cleaned census table read by the stem identification
stage (`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/`).

## Scripts (run in order, from the project root)

Both scripts are written to be run line by line: besides their outputs they
print inspection tables, and script 2 draws diagnostic plots. A failed check
prints ❌, raises a warning and writes to a log in `DATA/CHECKS/` before it
stops.

1. `1_prepare_species_tables.R`
   - Takes the six-letter species codes (`Mnemonic`) from the database
     export. The export's names are out of date and are used only for the
     three codes that no list names (see below).
   - Names every code from the botanists' lists:
     - the list of R. Pérez and S. Aguilar (email 1), which also gives the
       growth forms, common names and herbarium collections;
     - the name changes of S. Lao and the reply of R. Pérez (email 2);
     - the Flora de Panamá list, which has priority.
   - Takes the names of the three morphospecies that no list names from the
     database export. These are Apeiba hybrid (`apeihy`), Nectandra fuzzy
     (`nects1`) and Nectandra sp.4_(tiny_leaf) (`nects3`), the labels of
     Condit et al. (2019). It checks them against email 1, which calls them
     Apeiba "hybrida" and Nectandra s1 and s3:
     - the genus must be the one in the email;
     - the family must be the one in the email (Lauraceae for Nectandra) and
       that of the other species of the genus in the list.
   - Replaces the code `pterof` by `pterro`: all Pterocarpus in the plot are
     P. rohrii.
   - Checks every name against a local copy of the World Flora Online (WFO)
     backbone of June 2026, with the `WorldFlora` package. WFO only suggests
     corrections; each one is adopted by hand, in section 9.
   - Writes the species table `bci.spptable`: one row per code of the plot,
     except `uniden` (unidentified stems).

2. `2_prepare_viewfulltable.R`
   - Loads the raw `ViewFullTable` (DBH in mm, HOM in m). It removes census 9
     (2020, partial), replaces the code `pterof` by `pterro`, and renumbers
     the censuses 1..n by date.
   - Adds a placeholder row for every census between a tag's first and last
     record in which the tag has no row. The stage-3 R tables drop the
     placeholders that stage 2 leaves without a stem.
   - Keeps one measurement per stem and census: the highest HOM, with ties
     going to the larger DBH.
   - Flags likely DBH entry errors (a jump in and back out between three
     consecutive measurements). It replaces them, in a new column, by the
     geometric mean of their neighbours; the recorded `DBH` is kept.
   - Applies the Cushman et al. 2014 taper correction to 1.3 m for species
     with a tree or shrub growth form. Palms and tree ferns keep the
     uncorrected value.
   - Labels each tag as single-stem (one `StemID`, no `StemTag`) or
     multiple-stem. Only multiple-stem tags are reconstructed in stage 2.
   - Reads only the code and the growth form (`Lifeform_RPerez_SAguilar`)
     from the species table, so the species names do not change the stem
     identification.

## Inputs

Under `BCI_stem_reconstruction/DATA/RAW/`:

- `oct_13_2026/ViewFullTable_bci.csv`: census records, the export received
  on 2026-10-13 (tab-separated). Script 1 reads only its species codes and
  names (`Mnemonic`, `Family`, `Genus`, `SpeciesName`), and uses the names
  only for the three morphospecies.
- `ViewFiles_bci_allcensuses/ViewFullTable_bci.csv`: the previous export.
  Script 2 reads only its `DBHID` column, to check that the records are the
  same.
- `sp_tables/email_1/Lista_bci_mnemonics_formadevida.xlsx`: species list of
  R. Pérez and S. Aguilar.
- `sp_tables/email_2/bci_name_changes.xlsx`: name changes of S. Lao.
- `sp_tables/flora_panama_oct_9_2026/Lista de especies de la Flora de Panamá_nov25 original.xlsx`:
  Flora de Panamá list (R. Pérez, 2026-10-09).
- `sp_tables/wfo_plant_list_2026_06/classification.csv`: WFO Plant List,
  June 2026 (v.2026-06, CC0). It is the file `_DwC_backbone_R.zip` of
  doi:10.5281/zenodo.20782718, unzipped.

## Outputs

- `BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.txt` and `.csv` (both
  tab-separated), and `.RData` (object `bci.spptable`). The columns are:
  - `Mnemonic`;
  - `Order`, `Family`, `Genus`, `SpeciesName`, `InfraspecificRank`,
    `InfraspecificEpithet`, `Authority`, `Synonyms`;
  - `Lifeform_RFoster`, `Lifeform_RPerez_SAguilar`;
  - `CommonName`, `Herbarium`, `Notes`.
- `BCI_stem_reconstruction/DATA/CHECKS/`, the review of the WFO check:
  - `species_wfo_review.csv`: each name next to its WFO record;
  - `species_wfo_candidates.csv`: names with several WFO records;
  - `species_wfo_accepted_changes.csv`: the changes adopted in the review;
  - `species_wfo_vs_flora_panama.csv`, for the botanists: the names whose WFO
    record differs from the Flora de Panamá in family, genus, species,
    infraspecific name or authority, with both versions of each part. The
    authority is compared by its authors only (spaces, periods, capitals and
    accents aside);
  - `species_table_vs_wfo.csv`: the names of the final species table that
    still differ from their WFO record (family, genus, species,
    infraspecific name, authority, or only the writing of the authority), or
    that WFO treats as synonyms, with the accepted name.
- `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds`:
  the `ViewFullTable` columns plus:
  - `dbh_with_best_candidate_taper_corrected` (mm);
  - `single_stem_tags`;
  - `Lifeform`;
  - `RowID`.

## Requirements

R packages:

- script 1: `data.table`, `WorldFlora`, `stringr`, `janitor`, `readxl`;
- script 2: `data.table`, `ggplot2`, `inspectdf`.

## Notes

- This stage does not reconstruct stems; it only prepares the inputs for the
  stem identification stage.
