# BCI stem reconstruction

Four-stage pipeline that reconstructs the stem identities of the multiple-stem
trees of the BCI 50-ha plot and exports the ForestGEO R tables. Stage 1
prepares the census table, stage 2 reconstructs the stems with the `dp_global`
engine, stage 3 writes the R tables and stage 4 gives example analyses. Each
stage folder has its own README.

## Table of contents

- [Contents](#contents)
- [How to run](#how-to-run)
- [Data flow](#data-flow)
- [Conventions](#conventions)
- [Requirements](#requirements)

## Contents

| Folder | What it does |
| --- | --- |
| `1_DATA_PREPARATION/` | Builds the species table and the census table read by the engine |
| `2_STEM_IDENTIFICATION/` | Reconstructs the stems of the multiple-stem trees and merges the result into one table |
| `3_PREPARE_R_TABLES/` | Collects the posterior samples and writes the R tables, the species table and the check files |
| `4_EXAMPLE_STRUCTURE_ASSESSMENT/` | Computes biomass and basal-area stocks and fluxes, their identity uncertainty and plot summaries |
| `DATA/` | Inputs (`RAW/`) and the outputs of every stage; `DATA/README.md` describes each folder |

## How to run

Run every script from the project root (the folder that holds
`BCI_stem_reconstruction/` and `dp_global/`), in this order:

```sh
Rscript BCI_stem_reconstruction/1_DATA_PREPARATION/1_prepare_species_tables.R
Rscript BCI_stem_reconstruction/1_DATA_PREPARATION/2_prepare_viewfulltable.R
# engine run: command in 2_STEM_IDENTIFICATION/run_chunk_bci.md
Rscript BCI_stem_reconstruction/2_STEM_IDENTIFICATION/2_merge_chunks_to_datatable.R
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/1_prepare_posteriors_BCI.R
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/2_create_R_tables_BCI.R
Rscript BCI_stem_reconstruction/3_PREPARE_R_TABLES/3_create_species_table_BCI.R
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/biomass_stocks_fluxes.R
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/basal_area_uncertainty.R
Rscript BCI_stem_reconstruction/4_EXAMPLE_STRUCTURE_ASSESSMENT/general_plot_information.R
```

Before the merge and the posterior preparation, set `home_dir` (the
`BASE_OUT_DIR` of the engine run) and `run_code` (the run folder name) at the
top of `2_merge_chunks_to_datatable.R` and `1_prepare_posteriors_BCI.R`.

The scripts can also be run line by line. A failed check prints ❌, raises a
warning and writes to the log of its script before it stops.

The scripts guard against running out of order or on the wrong run: the engine
stops when stage 1 ran out of order; the merge stops when the run did not
finish, a chunk is missing or failed, the run was a test run (tag filter,
blind test, chunk range), it did not read the current stage-1 table, or its
growth bounds differ from the merge's; stage 3 stops when `run_code` is not the
merged run or a table is older than the one it is built from; stage 4 stops
when the R tables or the posterior samples are older than the stage-2 table or
come from another run.

## Data flow

```text
DATA/RAW/ (database export, species lists)
  → stage 1 → DATA/SPP_TABLE/bci_spptable.*, DATA/PROCESSED/ViewFullTable_single_vs_multiple_stem_tags.rds
  → stage 2 → engine run folder (outside the repository), DATA/<run_code>/merged_output.*,
              DATA/PROCESSED/complete_dataset_final_with_reconstructed_stemids.rds
  → stage 3 → DATA/POSTERIORS/, DATA/RTABLES/, DATA/CHECKS/
  → stage 4 → 4_EXAMPLE_STRUCTURE_ASSESSMENT/outputs/
```

## Conventions

- Censuses are numbered by date: 1 (1982), 2 (1985), 3 (1990), 4 (1995),
  5 (2000), 6 (2005), 7 (2010), 8 (2015), 9 (2022–23).
- DBH is in mm and HOM in m in stage 1 and in the R tables; the engine works in
  cm.
- A tree is a `tag` (`treeID` in the R tables). A stem is identified by
  `treeID` + `stemID`; `stemID` numbers the reconstructed stems of a tree.
- A field tag is a `StemTag` from census 7 (2010) on, the first census that
  tagged every stem. Before 2010 the database shows a stem's later tag on its
  earlier records.
- `DFstatus` is the raw field status and is never modified; `Rstatus` is the
  corrected status of stage 3.

## Requirements

R with the packages listed in each stage README, and a C++ compiler for the
engine (`dp_global/src/`).
