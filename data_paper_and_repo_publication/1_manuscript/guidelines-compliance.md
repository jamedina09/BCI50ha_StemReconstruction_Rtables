# Ecology Data Paper 2.0: what this submission meets and what is open

Checked on 2026-10-25, updated on 2026-10-28, against the files in `data_paper_and_repo_publication/guidelines/`: the Data Paper 2.0 Author Guidelines, the Submission Checklist, and the manuscript, metadata and README templates. Ecology requires this format for every Data Paper received from 1 September 2026.

Status: ✓ met · ✗ not yet met · ◐ met in part, or a choice of the authors.

## The files: one per element of a Data Paper

| Element | File | Goes to | Status |
| --- | --- | --- | --- |
| Manuscript | `1_manuscript/bci-50ha_data-paper_manuscript_v3.docx` (latest version), edited by hand | ScholarOne, as the "Main Document" (the only part typeset in Ecology) | ◐ Draft; what the authors must supply is highlighted |
| Metadata | `2_repository/2_Metadata.md` | Repository | ◐ Template filled; fields marked `[ ]` to complete |
| README | `2_repository/3_README.md` | Repository | ◐ Template filled; ORCIDs, DOIs and license to complete |
| Data | `2_repository/data/`: 11 files, 102.0 MB | Repository | ✓ The stage-3 files, unchanged; objects, columns, types, keys and the numbers of the metadata checked (78 checks, all passed) |
| Code | `2_repository/code/`: `load-bci-tables.R`, `1_tables.R`, `2_data-description.R`, `3_census-9-summary.R` | Repository | ✓ Tested: reproduces every number, table and the figure of the manuscript (outputs in `manuscript_outputs/`) |
| Code of the processing | GitHub repository of the pipeline | Repository or Zenodo | ✗ To archive with a DOI (see Open Research) |

The data files are the stage-3 outputs, copied unchanged (same MD5):

- the census tables `bci.stem1.Rdata` … `bci.stem9.Rdata`;
- the species table `bci.spptable.rdata`;
- the posterior samples `posterior_sampled_paths.rds`.

`bci.species_code_changes` stays internal and is not released.

## Manuscript

| Requirement | Draft | Status |
| --- | --- | --- |
| Title of at most 120 characters, sentence case | 72 characters | ✓ |
| Authors with their affiliation during data collection; present addresses at the bottom of the title page | Six authors, all at STRI. R. Condit, S. P. Hubbell and R. B. Foster are left out until they confirm. ORCIDs of S. J. Davies and J. A. Medina-Vega given; the other four blank | ◐ |
| One corresponding author with one email address | J. A. Medina-Vega, <medinaja@si.edu> | ✓ |
| 6–12 keywords, alphabetical, the same in manuscript, metadata, README and ScholarOne | 7, alphabetical, identical in the three files | ✓ |
| Acknowledgments: every funder | NSF grants to S. P. Hubbell, the Smithsonian Institution (STRI and ForestGEO) since the 2015 census, and the other funders of the project; the 2022–2023 census team named. Grant DEB-00753102 kept as written | ◐ |
| Field station named, with its location, in Methods and in Acknowledgments | Summary of Methods (Site) and Acknowledgments | ✓ |
| Contributor roles (CRediT) | Draft for each author to confirm | ◐ |
| Conflict-of-interest statement | Draft | ◐ |
| Abstract of at most 500 words that synthesizes every section | 350 words; its census-9 summary is also in the Data Description (Census 9) | ✓ |
| Sections in order: Overview, Data Availability, Data Description, Summary of Methods, Limitations and Challenges, Scientific Opportunities, Tables and Figures, References | All present, in order | ✓ |
| Main text of at most 2,000 words, title page, abstract, tables, captions and references excluded | 1,693 words with the section headings (1,680 without) | ✓ |
| At most four tables and figures, strictly descriptive | 3 tables (columns of the census tables, columns of the species table, status codes) and 1 figure (trees through the censuses), from `1_tables.R` | ✓ |
| Data Description strictly descriptive, with no statistical analysis | Counts and extents only | ✓ |
| Data Availability: DOI of the repository and license | Placeholders | ✗ |
| References complete, author–year | Every DOI checked; two book chapters lack page numbers; Muller-Landau and Wright (2024), the book of those chapters, is listed but not cited | ◐ |
| Word file | `.docx`, in the template's format: Arial 11 pt, 1.15 spacing, 0.75 in margins | ✓ |

Note 1: the census-9 basal area (31.35 m² ha⁻¹) uses the diameters as recorded (`use_taper_stocks <- FALSE` in `3_census-9-summary.R`); with the taper correction it would be 34.31 m² ha⁻¹. The script still writes the bootstrap intervals to `census-9_summary.csv`; the manuscript does not give them.

## Metadata (`2_repository/2_Metadata.md`)

✓ Every section and field of the template is present, in order:

- identifiers, authors and status;
- study design and site description;
- sampling design and research methods, including taxonomy and the full data processing;
- data file descriptors and observation information;
- variable information, with a table per file: name, definition, units, type, codes and values;
- archiving and history of usage.

✗ To complete:

- DOIs, license and citation;
- the ORCIDs of D. Mitre, R. Pérez, S. Aguilar and S. Lao;
- instrumentation, permits, data forms, storage of the forms and data-entry verification (field team and data manager);
- the herbaria;
- redundant archival sites;
- the date of the first release.

## README (`2_repository/3_README.md`)

✓ Template sections filled: Identifiers, Authors, Study Design, Data Files (naming scheme, file sizes, file contents). ✗ To complete: DOIs, license, the ORCIDs of four authors.

## Data formatting guidelines

| Guideline | This package | Status |
| --- | --- | --- |
| Tidy data: each variable a column, each observation a row, each kind of observational unit a table (exceptions case by case) | Census tables: one row per stem and census, one column per variable, one table per census (the ForestGEO R-table format). Every table lists every stem, before it appeared (P) and after it died (G, D). The posterior samples pack the labels of one sample into `path_sig` and `recon`; one row per record and sample would make the file many times larger. Both exceptions are explained in the metadata | ◐ |
| Nonproprietary formats such as `.csv` or `.txt` | R data files (`.Rdata`, `.rds`) by the authors' choice: R is free and open source, and the files total 102 MB, against 743 MB as CSV. The metadata gives the format, the loading commands and a reader for other software. The technical review may still ask for plain text | ◐ |
| Software as source code | R scripts | ✓ |
| File names: letters, digits, `-`, `_`, `.` only | All files | ✓ |
| File names of at most 35 characters | Longest: `posterior_sampled_paths.rds` (27) | ✓ |
| Dates in file names as YYYY-MM-DD; words joined by `-`, groupings by `_` | No dates in names. The data files keep the ForestGEO names (`bci.stem1.Rdata`), which use `.`; the code files use `-` | ◐ |
| Unique, descriptive file names | Yes | ✓ |

## Open Research policy

| Requirement | Status |
| --- | --- |
| Data, metadata and README in a suitable repository at submission (published, or private for peer review) | ✗ Choose a repository. Dryad, Zenodo and Figshare are pre-approved; Dryad is free for Ecology Data Papers up to 10 GB but publishes under CC0 only |
| Persistent identifier (DOI), open license, ORCID integration | ✗ With the repository |
| Code that processes or summarizes the data and makes the tables of the paper | ✓ `2_repository/code/` |
| Code of a new method, for peer review | ✗ The stem reconstruction is a new method. The manuscript links the public GitHub repository; archive a release with a DOI (Zenodo) and cite it, or decide otherwise |

## Submission checklist

| Item | Status |
| --- | --- |
| Data formatted according to the data formatting guidelines | ◐ See above |
| Metadata formatted according to the metadata guidelines | ◐ Fields to complete |
| README formatted according to the README guidelines | ◐ Fields to complete |
| Data, metadata and README uploaded to an open-access repository from the list (private for peer review or published) | ✗ |
| Manuscript formatted according to the manuscript guidelines | ◐ Placeholders |
| Link to the data repository in the Data Availability section | ✗ After deposit |
| Manuscript submitted as "Main Document" on ScholarOne | ✗ |

## Open items for the authors

1. **Authors:**
   - R. Condit, S. P. Hubbell and R. B. Foster: left out of the author list, the contributor roles and the citation until they confirm (S. Davies to write to them); if they join, add them back with their affiliations;
   - the ORCIDs of D. Mitre, R. Pérez, S. Aguilar and S. Lao (blank on the title page);
   - the role of D. Mitre.
2. **Acknowledgments:** research permits.
3. **Species:**
   - the note that *Nectandra* sp.4 was merged into `nects3` is not supported by the database (one tree, one code and one label in both exports): leave it out unless the botanists confirm it;
   - check the list against S. Lao's latest list and against H. Muller-Landau's version;
   - choose how to cite the Flora de Panamá list.
4. **License and repository:** both options stay open: Dryad (CC0 only; free for Ecology Data Papers up to 10 GB) or Zenodo (CC BY 4.0). See Open Research.
5. **File names:** the Data Availability section names the metadata `Metadata.md`; the files are now `2_Metadata.md` and `3_README.md`. Rename them for the deposit or change the text.
6. **References:** cite Muller-Landau and Wright (2024) or delete it (the two chapters already name the book).

## How to reproduce the outputs of the manuscript

In R (or with `Rscript`), with the working directory set to `2_repository/`, run:

```r
source("code/1_tables.R")             # the figure and the tables
source("code/2_data-description.R")   # numbers of the Data Description and Limitations
source("code/3_census-9-summary.R")   # census-9 summary of the Abstract and the Data Description
```

The scripts write CSV files to `manuscript_outputs/`:

- `figure-1_tree-dynamics.png`, `.pdf` and `.csv`, and its caption `figure-1_caption.txt`;
- `table_columns-census-tables.csv`, `table_columns-species-table.csv`, `table_status-codes.csv`;
- `data-description_numbers.csv`;
- `census-9_summary.csv`, `census-9_species.csv`, `census-9_families.csv`.

They need the R packages data.table, boot and HDInterval.

## Upload

Upload the contents of `2_repository/`: `3_README.md`, `2_Metadata.md`, `data/` and `code/`. Google Drive puts hidden `Icon` files (`Icon\r`) in every synced folder; leave them, and any `.DS_Store`, out of the upload.
