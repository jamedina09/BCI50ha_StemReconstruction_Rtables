## Identifiers
(data set) title: Tree censuses of the Barro Colorado Island 50-ha plot, Panama, 1982–2023
data identifier DOI: [to be issued by the repository]
data paper DOI: [to be issued by Ecology]
license: [CC0 1.0 or CC BY 4.0, to decide]
keywords: Barro Colorado Island; demography; forest dynamics plot; ForestGEO; long-term monitoring; Panama; tropical moist forest
data paper citation: Mitre, D., R. Pérez, S. Aguilar, S. Lao, S. J. Davies, and J. A. Medina-Vega. [Year]. "Tree censuses of the Barro Colorado Island 50-ha plot, Panama, 1982–2023." Ecology [volume]: [article number]. [DOI]

## Authors (for each author)
author name (req): David Mitre
author ORCID (req): [to add]
author email (req-corresponding; optional-others):
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

author name (req): Rolando Pérez
author ORCID (req): [to add]
author email (req-corresponding; optional-others):
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

author name (req): Salomón Aguilar
author ORCID (req): 0000-0002-8694-6723
author email (req-corresponding; optional-others):
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

author name (req): Suzanne Lao
author ORCID (req): [to add]
author email (req-corresponding; optional-others):
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

author name (req): Stuart J. Davies
author ORCID (req): 0000-0002-8596-7522
author email (req-corresponding; optional-others):
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

author name (req): José A. Medina-Vega
author ORCID (req): 0000-0001-5468-5605
author email (req-corresponding; optional-others): <medinaja@si.edu>
author affiliation (optional): Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama

## Study Design
project title: Barro Colorado Island 50-ha Forest Dynamics Plot, Forest Global Earth Observatory (ForestGEO), Smithsonian Tropical Research Institute
objectives: To record the growth, mortality and recruitment of every free-standing woody stem of at least 1 cm diameter at breast height in 50 ha of old-growth tropical moist forest, about every five years since 1982, as a basis for studies of tropical forest diversity, demography and change. This data set holds the nine censuses (1982–2023) with the identity of every stem reconstructed across all censuses.
funding: U.S. National Science Foundation (grants to S. P. Hubbell: DEB-0640386, DEB-0425651, DEB-0346488, DEB-0129874, DEB-00753102, DEB-9909347, DEB-9615226, DEB-9405933, DEB-9221033, DEB-9100058, DEB-8906869, DEB-8605042, DEB-8206992, DEB-7922197); since the 2015 census, the Smithsonian Institution, through the Smithsonian Tropical Research Institute and the Forest Global Earth Observatory (ForestGEO; formerly the Center for Tropical Forest Science and the Smithsonian Institution Global Earth Observatories), the primary funder; HSBC Climate Partnership; John D. and Catherine T. MacArthur Foundation; Andrew W. Mellon Foundation; Frank H. Levinson Family Foundation; Small World Institute Fund; numerous private individuals; and, in the first two decades of the project, the Smithsonian Scholarly Studies Program, World Wildlife Fund, Earthwatch Center for Field Studies, National Geographic Society, Geraldine R. Dodge Foundation, W. Alton Jones Foundation, Pew Charitable Trusts and Guggenheim Foundation.

## Data Files

naming scheme: data/bci.stem<k>.Rdata is census k (1–9) and holds the ForestGEO R table bci.stem<k>; data/bci.spptable.rdata holds the species table bci.spptable; data/posterior_sampled_paths.rds holds the posterior samples of the stem reconstruction; code/ holds R scripts; Metadata.md describes every file and column.
file size(s): 102.1 MB in total.
  - bci.stem1.Rdata … bci.stem9.Rdata: 6.6–7.7 MB each (646,813 rows × 17 columns each)
  - bci.spptable.rdata: 0.02 MB (328 rows × 17 columns)
  - posterior_sampled_paths.rds: 36.1 MB (1,554,448 rows × 6 columns)
  - code/load-bci-tables.R, code/1_tables.R, code/2_data-description.R, code/3_census-9-summary.R: < 0.01 MB each
file contents:
  - Census tables (bci.stem1 … bci.stem9): one row per stem ever recorded in the plot, the same stems in the same order in all nine tables; columns CensusID, treeID, stemID, tag, StemTag, sp, quadrat, gx, gy, dbh (mm), hom (m), ExactDate, date (days since 1960-01-01), codes, DFstatus (field status), Rstatus (P, A, G, D), StemPaths. A stem is identified by treeID + stemID.
  - Species table (bci.spptable): one row per species code (sp), with taxonomy, authority, synonyms, growth forms, common names, herbarium collections, the 2019 name and the kind of name change.
  - Posterior samples: one row per distinct sampled reconstruction of a tree (tag, treeID, path_sig, path_count, path_prob, recon); join to the census tables on treeID + StemPaths.
  - code/load-bci-tables.R loads the tables and draws reconstructions from the posterior samples. code/1_tables.R (Figure 1, its caption and Tables 1–3 of the Data Paper), code/2_data-description.R (numbers of the Data Description and Limitations) and code/3_census-9-summary.R (census-9 summary of the Abstract and the Data Description) reproduce the tables and numbers of the Data Paper; run them from this folder.
  - Format: R data files. load("data/bci.stem9.Rdata") creates the data.frame bci.stem9; readRDS("data/posterior_sampled_paths.rds") returns the posterior samples. Missing values are NA; text columns are text, so tags keep their leading zeros.
