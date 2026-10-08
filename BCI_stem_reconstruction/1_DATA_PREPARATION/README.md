# 1_DATA_PREPARATION

This folder contains the first stage of the BCI stem reconstruction pipeline.
It converts the raw ForestGEO exports into the species table and the cleaned
census table read by the stem identification stage
(`BCI_stem_reconstruction/2_STEM_IDENTIFICATION/`).

## Scripts (run in order, from the project root)

Both scripts are written to be run interactively: besides their outputs they
print inspection tables and draw diagnostic plots.

1. `0_prepare_species_tables.R`
   - Reads the species list with growth forms, the database taxonomy and the
     mnemonics used in the census records.
   - Checks the names against the Taxonomic Name Resolution Service (TNRS;
     needs an internet connection), updates families and authorities, and
     completes the three codes without a genus from the database taxonomy.
   - Writes the species table `bci.spptable` (one row per mnemonic).

2. `1_prepare_viewfulltable.R.R`
   - Loads the raw `ViewFullTable` (DBH in mm, HOM in m), replaces the
     mnemonic `pterof` by `pterro`, and renumbers the censuses 1..n by date.
   - Adds a placeholder row for every census between a tag's first and last
     record in which the tag has no row.
   - Keeps one measurement per stem and census (the highest HOM; ties go to
     the larger DBH).
   - Flags likely DBH entry errors (a jump in and back out between three
     consecutive measurements) and replaces them by the geometric mean of
     their neighbours in a new column; the recorded `DBH` is kept.
   - Applies the Cushman et al. 2014 taper correction to 1.3 m for species
     with a tree or shrub growth form (palms and tree ferns keep the
     uncorrected value).
   - Labels each tag as single-stem (one `StemID`, no `StemTag`) or
     multiple-stem; only multiple-stem tags are reconstructed in stage 2.

## Inputs

Under `BCI_stem_reconstruction/DATA/RAW/`:

- `ViewFiles_bci_allcensuses/ViewFullTable_bci.csv` — census records
  (tab-separated)
- `ViewFiles_bci_allcensuses/ViewTaxonomy_bci.csv` — database taxonomy
- `sp_tables/Lista_bci_mnemonics_formadevida.xlsx` — species list with growth
  forms (the script reads it from this path)

## Outputs

- `BCI_stem_reconstruction/DATA/SPP_TABLE/bci_spptable.txt`, `.csv` (both
  tab-separated) and `.RData` (object `bci.spptable`). The folder must exist.
- `BCI_stem_reconstruction/DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds`
  — the `ViewFullTable` columns plus
  `dbh_with_best_candidate_taper_corrected` (mm), `single_stem_tags`,
  `Lifeform` and `RowID`.

## Requirements

R packages: `data.table`, `TNRS`, `stringr`, `readxl`, `inspectdf`, `ggplot2`.

## Notes

- This stage does not perform stem-level reconstruction; it only prepares the
  inputs for the stem identification stage.
