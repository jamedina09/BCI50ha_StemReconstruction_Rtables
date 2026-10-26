::: {custom-style="Title"}
Tree censuses of the Barro Colorado Island 50-ha plot, Panama, 1982–2023, with reconstructed stem histories
:::

David Mitre^1^, Rolando Pérez^1^, Salomón Aguilar^1^, Suzanne Lao^1^, Richard Condit^2^, Stephen P. Hubbell^3^, Robin B. Foster^4^, Stuart J. Davies^1^, José A. Medina-Vega^1^

^1^ Smithsonian Tropical Research Institute, Apartado 0843-03092, Balboa, Ancón, Republic of Panama\
^2^ [Affiliation to confirm (University of California, Santa Cruz, in Condit et al. 2019a)]{.mark}\
^3^ [Affiliation to confirm (University of California, Los Angeles, in Condit et al. 2019a)]{.mark}\
^4^ [Affiliation to confirm (Field Museum of Natural History, in Condit et al. 2019a)]{.mark}

Corresponding author: [name, email]{.mark}

# Keywords

Barro Colorado Island; demography; forest dynamics plot; ForestGEO; long-term monitoring; mortality; multi-stemmed trees; Panama; recruitment; stem identity reconstruction; tree growth; tropical moist forest

# Acknowledgments

**Funding.** The BCI forest dynamics research project was made possible by National Science Foundation grants to Stephen P. Hubbell: DEB-0640386, DEB-0425651, DEB-0346488, DEB-0129874, DEB-00753102, DEB-9909347, DEB-9615226, DEB-9405933, DEB-9221033, DEB-9100058, DEB-8906869, DEB-8605042, DEB-8206992, DEB-7922197, support from the Forest Global Earth Observatory (ForestGEO), the Smithsonian Tropical Research Institute, the John D. and Catherine T. MacArthur Foundation, the Mellon Foundation, the Small World Institute Fund, and numerous private individuals, and through the hard work of over 100 people from 10 countries over the past four decades. The plot project is part of ForestGEO, a global network of large-scale demographic tree plots. [This is ForestGEO's required wording, with "three decades" changed to "four". ForestGEO lists DEB-9615226 twice and writes DEB-00753102, which has one digit too many; confirm the number. Add the support of the 2010, 2015 and 2022–2023 censuses: the Smithsonian Institution (Office of the Under Secretary for Science and Research), ForestGEO and its donors, including the Frank H. Levinson Family Foundation, with grant numbers.]{.mark}

**Field stations, marine laboratories and related research institutions.** The data were collected in the Barro Colorado Nature Monument, Panama, at the Barro Colorado Island research station of the Smithsonian Tropical Research Institute (9.15° N, 79.85° W). We thank the Smithsonian Tropical Research Institute for logistical support [and the authority that issued the research permits]{.mark}.

**Additional acknowledgments.** We thank the hundreds of field assistants, botanists and data managers who have censused the plot since 1980, and the team of the 2022–2023 census: [names, from D. Mitre]{.mark}. Their names and photographs are on the ForestGEO website ([URL]{.mark}). We thank Frank H. Levinson for his long support of ForestGEO [confirm wording]{.mark}, and [H. C. Muller-Landau and other reviewers of the data package]{.mark}.

# Contributor Roles

[Draft; each author to confirm.]{.mark} David Mitre: investigation, data curation. Rolando Pérez and Salomón Aguilar: investigation, data curation (taxonomy). Suzanne Lao: data curation, investigation. Richard Condit: conceptualization, methodology, investigation, data curation, funding acquisition, project administration. Stephen P. Hubbell: conceptualization, methodology, investigation, funding acquisition, project administration. Robin B. Foster: conceptualization, methodology, investigation. Stuart J. Davies: funding acquisition, project administration, resources, supervision. José A. Medina-Vega: methodology (stem reconstruction), software, formal analysis, data curation (R tables), validation, writing – original draft. All authors: writing – review and editing.

# Conflict of Interest

The authors declare no conflicts of interest. [Confirm.]{.mark}

# Abstract

The 50-ha Forest Dynamics Plot on Barro Colorado Island (BCI), Panama, is the oldest plot of the Forest Global Earth Observatory (ForestGEO) and one of the most studied tropical forests in the world. In each census, every free-standing woody stem at least 1 cm in diameter at breast height (DBH) in the 1000 × 500 m plot of lowland moist forest is tagged, mapped, identified and measured. We present nine censuses (1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015 and 2022–2023) as nine census tables in the ForestGEO R-table format, with a species table. They describe 460,211 trees and 646,813 stems of 327 taxa in 62 families, with 2,325,471 diameter measurements. This release adds the ninth census, updates the species names, and follows stems consistently across all censuses. Before 2010 only trees were tagged, and the database gave the secondary stems of multi-stemmed trees new identifiers in many censuses, so their histories could not be followed; since 2010 every stem has been tagged. For the 132,220 trees with more than one database stem identifier or a stem tag, we reconstructed the history of every stem with a dynamic-programming model of growth, mortality and recruitment estimated from the tagged censuses, keeping every identity fixed by a field tag or by the database, and we provide posterior samples so that the uncertainty of the reconstruction can be propagated. All nine tables list the same stems in the same order. Each row keeps the field record unchanged (DBH, height of measurement, date, codes and status) and adds an analytical status (prior, alive, dead stem of a living tree, dead tree) assigned by explicit rules. In 2022–2023 the plot held 203,489 living trees (4,070 ha^−1^; bootstrap 95% interval 4,023–4,114) with 244,204 living stems, of 300 species, 180 genera and 59 families, and a basal area of 34.31 m^2^ ha^−1^ (95% interval 33.00–35.53). Abundance was uneven: 212 species (70.7%) had fewer than 339 individuals (half the mean abundance), 67 had ten or fewer, and 22 were singletons. *Pombalia prunifolia*, *Faramea occidentalis* and *Trichilia tuberculata* made up 27.7% of the individuals, while *Quararibea stenophylla*, *Alseis blackiana* and *Prioria copaifera* had the largest basal areas. Rubiaceae, Violaceae and Fabaceae had the most individuals, and Malvaceae, Fabaceae and Rubiaceae the largest basal areas. Users should note that stem identities before 2010 are reconstructed, not observed; that diameters below 55 mm were recorded in 5-mm classes in 1982 and 1985; that the height of measurement rises over time on buttressed and damaged trunks; and that census 9 was taken 7.5 years after census 8. The data support studies of tree demography and species coexistence, of multi-stemmed and resprouting trees, of biomass and carbon dynamics, and of forest responses to climate. They supersede the 2019 release (Condit et al. 2019a).

# Overview

Large forest dynamics plots, in which every stem above a small size is mapped, identified and remeasured at regular intervals, underpin much of what is known about tropical forest diversity and dynamics (Condit 1998; Anderson-Teixeira et al. 2015; Davies et al. 2021). The 50-ha plot on Barro Colorado Island (BCI), Panama, was the first. Stephen Hubbell and Robin Foster established it in 1980; its first census, the 1982 census, ran from 1981 to 1983 (Hubbell and Foster 1983; Hubbell et al. 2024). It became the model for the Forest Global Earth Observatory (ForestGEO), now more than 70 plots (Davies et al. 2021). Its censuses have been used to study density dependence and coexistence (Hubbell et al. 1990, 1999; Comita et al. 2010), habitat associations (Harms et al. 2001), immigration (Condit et al. 2012), demographic trade-offs (Rüger et al. 2018, 2020), compositional change and climate (Feeley et al. 2011; Condit et al. 2017), and biomass dynamics (Chave et al. 2003; Rutishauser et al. 2020; Piponiot et al. 2022, 2024). Muller-Landau and Wright (2024) review a century of research on the island.

The last public release (Condit et al. 2019a) held eight censuses, 1982–2015, with stems as the database recorded them. This release adds the ninth census, taken in 2022–2023 after the census begun in 2020 was stopped by the COVID-19 pandemic. It updates the species names and documents the names that changed since 2019. And it makes stems comparable across censuses. Before 2010 the database followed trees, not stems: the tree tag was on the main stem, and the other stems of a tree often received a new database identifier at each census, which cut their histories. Since 2010 every stem has carried its own tag. We reconstructed the earlier stem histories under a model of growth, mortality and recruitment, with their uncertainty, so that stem growth, mortality and recruitment can be computed over four decades.

The data are for ecologists who study tropical forest demography and diversity, for carbon-cycle scientists, and for modellers who calibrate or test vegetation models with plot data.

# Data Availability

The data, metadata (`Metadata.md`), README and code are archived at [repository, https://doi.org/xxxx]{.mark} under [license: CC0 1.0 or CC BY 4.0]{.mark}. The data are R data files (Table 2): the nine census tables and the species table in the ForestGEO R-table format, and the posterior samples of the stem reconstruction. R is free and open source. The package code loads the files and reproduces every table and number of this paper. The code that builds the tables from the export of the ForestGEO database, including the stem reconstruction, is at <https://github.com/jamedina09/BCI50ha_StemReconstruction_Rtables>, archived at [Zenodo DOI]{.mark}. Please cite this paper and the data package, and the publications that ForestGEO asks users of the BCI data to cite (Condit 1998; Hubbell et al. 1999; Condit et al. 2019a, 2019b). Please work from the census tables rather than rebuild them from the database records, so that analyses of the plot remain comparable. Questions go to the corresponding author; copies of publications to ForestGEO (forestgeo@si.edu).

# Data Description

**Extent and resolution.** The plot is a 1000 × 500 m rectangle (50 ha) of 1,250 quadrats of 20 × 20 m. The nine censuses span March 1981 to December 2023 (Table 1); a census took from nine months (1985) to almost three years (1990–1992). Consecutive censuses of a tree were about 5 years apart, except 3.4 years (censuses 1–2) and 7.5 years (8–9). Every record carries its date.

**Taxa.** The censuses include every free-standing woody plant: trees, shrubs, palms, tree ferns and strangler figs; lianas are not included. The 327 taxa of the species table belong to 192 genera and 62 families. Four are morphospecies, and one further code (`uniden`) marks 74 trees never identified.

**Sampling unit.** Each row is a stem in one census. A tree (`treeID`, with its field `tag`) has one or more stems, and `treeID` and `stemID` together identify a stem. Of the 460,211 trees, 104,889 had more than one stem during the study.

**Tables.** The data package holds nine census tables, a species table and the posterior samples of the stem reconstruction (Table 2). The census tables (`bci.stem1` to `bci.stem9`) follow the ForestGEO R-table format, one table per census. Each has 646,813 rows and 17 columns (Table 3): one row for every stem ever recorded, in the same order in all tables (by `treeID`, then `stemID`), including the censuses before the stem appeared and after it died. The species table (328 rows, 17 columns) gives for each species code the taxonomy, authority, synonyms, growth forms, common name, herbarium collections and, for the 49 codes renamed since 2019, the 2019 name and the kind of change; it joins the census tables on `sp`. The posterior file holds one row per distinct reconstruction sampled for each of 131,999 trees (1,554,448 rows) and joins the census tables on `treeID` and `StemPaths`.

**Sample sizes.** The tables hold 2,818,748 field records, 2,325,471 of them with a diameter. Living trees ranged from 244,252 (1990) to 203,489 (2022–2023), and living stems from 287,115 (1990) to 244,204 (2022–2023; Table 1).

**Records and status.** Each row keeps the field record unedited: diameter (`dbh`, mm), height of measurement (`hom`, m), date (`ExactDate`), codes and field status (`DFstatus`); these are empty in a census without a record. `Rstatus` gives every stem a status in every census (Table 4), `date` gives every row a date, and `gx` and `gy` give one location per tree.

# Summary of Methods

**Site.** BCI (9.15° N, 79.85° W) is a 15.6-km^2^ island in Gatun Lake, Panama, isolated when the lake was formed during the construction of the Panama Canal (1910–1914) and protected since 1923. It is part of the Barro Colorado Nature Monument, where the Smithsonian Tropical Research Institute runs a field station. The plot lies on the central plateau of the island, at 120–160 m above sea level, in old-growth lowland tropical moist forest that includes a small area of younger forest (Harms et al. 2001). Rainfall averages about 2,600 mm per year, with a dry season from January to April (Leigh 1999).

**Census protocol.** The censuses follow Condit (1998) and the ForestGEO standards (Condit et al. 2014). Every free-standing woody stem with a DBH of at least 10 mm is mapped within its 20 × 20 m quadrat, identified and measured at 1.3 m or, where the trunk is buttressed or deformed there, above the irregularity, and the height of measurement is recorded. Each tree carries a numbered tag; since 2010 every stem has also carried its own tag, whereas before 2010 the other stems of a tree were measured untagged. Each census revisits every tree recorded before and adds the stems that have reached 10 mm. In 1982 and 1985 diameters below 55 mm were recorded in 5-mm classes. The plot's botanists identify the plants in the field.

**From the database to the tables.** The records were exported from the plot's ForestGEO database on 13 October 2026, without the partial census of 2020 (135 of 1,250 quadrats). A stem measured at two heights in one census keeps the higher measurement. Species codes come from the database and names from the lists of the plot's botanists, with the Flora de Panamá list taking priority, checked against World Flora Online (WFO 2026; Kindt 2020).

**Stem reconstruction.** The 327,991 trees with a single database stem and no stem tag needed no reconstruction. In the other 132,220 trees, the identities fixed in the field were kept: every stem tagged from 2010 on, and the database stem identifiers that link the records of one stem, unless the measurements contradict them. An exact dynamic program then assigned the remaining measurements to stems, finding the most probable stem histories under a model of growth, mortality and recruitment estimated from the stems tagged in 2010–2023. Field codes constrain the solution: a stump (broken, with no measurable diameter) ends the stem that broke, and a broken or resprouting stem measured at 1.3 m starts a new stem. Palm clumps, strangler figs and trees too large for the exact program were solved by a probabilistic matcher under the same model. Stems split only because the point of measurement moved were rejoined. The reconstruction of each tree was sampled 200 times, and the tables carry the most probable one (Part 2 and the code give the full method).

**Status, dates and locations.** `Rstatus` follows explicit rules (Table 4). A row without a recorded date takes the most common date of its tree in that census, else of its quadrat, else of the census, and each tree takes the most common coordinates of its stems. Every step is scripted and stops if a check fails.

# Limitations and Challenges

**Reconstructed stems.** Before 2010 the stems of multi-stemmed trees are inferred, not observed. Tree identity (`treeID`, `tag`) is observed, so analyses that treat a tree as the sum of its stems do not depend on the reconstruction. Stem-level analyses before 2010 should be repeated over draws from the posterior samples, one reconstruction per tree weighted by `path_count`. The samples of the dynamic program come from the exact posterior under the model; those of the matcher (palm clumps, strangler figs, the largest trees) are approximate.

**Diameter classes.** In 1982 and 1985 diameters below 55 mm were recorded in 5-mm classes, so the growth of small stems over the first two intervals is coarse.

**Height of measurement.** On buttressed or damaged trunks the point of measurement was raised over time; the 1.1% of measurements taken away from 1.3 m are not comparable with those at 1.3 m. For basal area and biomass they should be corrected to 1.3 m, for example with the taper model of Cushman et al. (2014), except in palms.

**Missing diameters.** An alive stem without a diameter is either a missed measurement (1,789 cells, which can be interpolated from the stem's own measurements) or a stump (777 cells: alive, broken below, with no stem at 1.3 m, to be left out of stocks).

**Census timing.** Censuses lasted up to 2.8 years, intervals differ among quadrats, census 9 came 7.5 years after census 8, and the interrupted census of 2020 is not included. Intervals should be computed from `date`, not from census years.

**Recruitment and ingrowth.** New stems of established trees (ingrowth) are not new individuals (recruitment) and should be counted separately.

**Taxonomy.** Names may change again; species codes are stable within this release.

**Strangler figs and palms.** The diameter of strangler figs that envelop their hosts, and that of palms, does not relate to biomass as it does in other trees (Rutishauser et al. 2020; Piponiot et al. 2024).

**Locations.** Each tree has one location; 86 trees have no coordinates and carry only their quadrat.

# Scientific Opportunities

Nine censuses over four decades, with stems that can be followed throughout, make stem-level questions open that earlier releases could not address. Multiple stems and resprouting are common in tropical trees but rarely followed over time: the tables allow rates of stem ingrowth, stem death and resprouting after damage to be estimated for some 300 species and compared with tree-level recruitment and mortality. The ninth census, after an interval that spans the 2015–2016 El Niño drought, extends tests of directional change in composition, growth and mortality (Feeley et al. 2011; Condit et al. 2017) and of the stability of biomass (Rutishauser et al. 2020; Piponiot et al. 2024). The species table, with growth forms and documented name changes, supports trait-based and phylogenetic analyses of demographic trade-offs (Rüger et al. 2018, 2020) and coexistence (Hubbell et al. 1999; Comita et al. 2010). The posterior samples let the uncertainty of stem identity be carried into estimates of stocks and fluxes, a source of error rarely quantified in plot data. The tables can also calibrate remote sensing of canopy change and test individual-based and demographic vegetation models. Finally, the reconstruction method applies to other ForestGEO plots whose early censuses did not tag every stem.

# Tables

**Table 1.** Censuses of the BCI 50-ha plot in this release. Records: rows with a field record (`DFstatus`). Measurements: records with a diameter. Living stems: `Rstatus` A; living trees: trees with at least one living stem.

| Census | Year | Field dates | Records | Measurements | Living stems | Living trees |
|:--|:--|:----------|--:|--:|--:|--:|
| 1 | 1982 | 1981-03-11 to 1983-07-30 | 263,934 | 263,896 | 263,896 | 235,336 |
| 2 | 1985 | 1985-01-05 to 1985-10-05 | 299,447 | 273,062 | 273,598 | 242,611 |
| 3 | 1990 | 1990-02-06 to 1992-12-04 | 328,634 | 286,910 | 287,115 | 244,252 |
| 4 | 1995 | 1995-01-02 to 1996-02-15 | 317,354 | 267,448 | 267,625 | 229,213 |
| 5 | 2000 | 2000-01-03 to 2001-03-31 | 298,817 | 249,546 | 249,831 | 214,066 |
| 6 | 2005 | 2005-01-03 to 2006-01-24 | 287,494 | 244,793 | 244,965 | 208,542 |
| 7 | 2010 | 2010-01-18 to 2011-03-28 | 317,226 | 246,898 | 247,731 | 208,079 |
| 8 | 2015 | 2013-09-01 to 2016-02-18 | 335,207 | 248,713 | 249,054 | 207,981 |
| 9 | 2022 | 2022-02-15 to 2023-12-23 | 370,635 | 244,205 | 244,204 | 203,489 |

**Table 2.** Files of the data package (102 MB). The census tables and the species table are R data files (`.Rdata`, one data.frame each, named as in the table); the posterior samples are an R data file written by `saveRDS()` (`.rds`).

| File | Rows × columns | Content |
|:-------------|:-----|:--------------------|
| `data/bci.stem1.Rdata` … `bci.stem9.Rdata` | 646,813 × 17 each | Census tables `bci.stem1` … `bci.stem9`: one row per stem, the same stems in the same order in all nine (Table 3) |
| `data/bci.spptable.rdata` | 328 × 17 | Species table `bci.spptable`: one row per species code, `uniden` included |
| `data/posterior_sampled_paths.rds` | 1,554,448 × 6 | Posterior samples of the stem reconstruction: one row per distinct sampled reconstruction of a tree |
| `code/load-bci-tables.R` | | Loads the tables and draws reconstructions from the posterior samples |
| `code/make-manuscript-tables.R` | | Reproduces the tables and numbers of this paper |
| `Metadata.md`, `README.md` | | Documentation |

**Table 3.** Columns of the census tables `bci.stem1` to `bci.stem9`. Empty cells are `NA`.

| Column | Type | Content |
|:-----|:----|:------------------------|
| `CensusID` | number | Census, 1–9 (Table 1) |
| `treeID` | text | Tree; the same in every census |
| `stemID` | integer | Stem number within its tree (1, 2, …). `treeID` and `stemID` together identify a stem |
| `tag` | text | Field tag of the tree, with leading zeros and sign kept (e.g. `000006`, `-05599`) |
| `StemTag` | text | Field tag of the stem, used since 2010; earlier records may show the tag the stem received later |
| `sp` | text | Species code; joins the species table |
| `quadrat` | text | 20 × 20 m quadrat, `XXYY`: column (`gx`/20) and row (`gy`/20) |
| `gx`, `gy` | number, m | Location of the tree, east of the western edge (0–1000) and north of the southern edge (0–500) |
| `dbh` | number, mm | Diameter as recorded, at height `hom`; empty when not measured |
| `hom` | number, m | Height of measurement as recorded |
| `ExactDate` | date | Date of the field record; empty without a record |
| `date` | number, days | Date of every row, in days since 1 January 1960: the recorded date, else the most common date of the tree, the quadrat or the census |
| `codes` | text | Field codes as recorded, separated by `;` |
| `DFstatus` | text | Field status as recorded: alive, broken below, dead, stem dead, missing; empty without a record |
| `Rstatus` | text | Analytical status: P, A, G or D (Table 4) |
| `StemPaths` | integer | Key of the record in the reconstruction, unique within its tree; the posterior samples list the measured records up to the tree's anchor census (2010 for most trees). Empty for trees not reconstructed |

**Table 4.** Analytical status (`Rstatus`). Evidence of life is an alive record, or a diameter on a broken-below, missing or status-less record. A dead record is dead or stem dead, or broken below without a diameter. Counts over the nine tables: A 2,328,019; P 1,755,546; D 1,377,759; G 359,993.

| Code | Meaning | Rule |
|:--|:----|:--------------------------|
| P | prior | Before the stem's first record. A stem with only missing records is P in every census |
| A | alive | From the first record to the last evidence of life, also in censuses without a record or with a dead record in between. An alive cell without a diameter is a missed measurement or a stump |
| G | dead stem | Dead, while another stem of its tree is A in that census or later |
| D | dead tree | Dead, and every stem of its tree dead |

# References

Anderson-Teixeira, K. J., S. J. Davies, A. C. Bennett, E. B. Gonzalez-Akre, H. C. Muller-Landau, S. J. Wright, K. Abu Salim, et al. 2015. "CTFS-ForestGEO: A Worldwide Network Monitoring Forests in an Era of Global Change." *Global Change Biology* 21: 528–549. <https://doi.org/10.1111/gcb.12712>

Chave, J., R. Condit, S. Lao, J. P. Caspersen, R. B. Foster, and S. P. Hubbell. 2003. "Spatial and Temporal Variation of Biomass in a Tropical Forest: Results from a Large Census Plot in Panama." *Journal of Ecology* 91: 240–252. <https://doi.org/10.1046/j.1365-2745.2003.00757.x>

Comita, L. S., H. C. Muller-Landau, S. Aguilar, and S. P. Hubbell. 2010. "Asymmetric Density Dependence Shapes Species Abundances in a Tropical Tree Community." *Science* 329: 330–332. <https://doi.org/10.1126/science.1190772>

Condit, R. 1998. *Tropical Forest Census Plots: Methods and Results from Barro Colorado Island, Panama and a Comparison with Other Plots*. Berlin: Springer-Verlag, and Georgetown, Texas: R. G. Landes Company. <https://doi.org/10.1007/978-3-662-03664-8>

Condit, R., R. A. Chisholm, and S. P. Hubbell. 2012. "Thirty Years of Forest Census at Barro Colorado and the Importance of Immigration in Maintaining Diversity." *PLoS ONE* 7: e49826. <https://doi.org/10.1371/journal.pone.0049826>

Condit, R., S. Lao, A. Singh, S. Esufali, and S. Dolins. 2014. "Data and Database Standards for Permanent Forest Plots in a Global Network." *Forest Ecology and Management* 316: 21–31. <https://doi.org/10.1016/j.foreco.2013.09.011>

Condit, R., R. Pérez, S. Lao, S. Aguilar, and S. P. Hubbell. 2017. "Demographic Trends and Climate over 35 Years in the Barro Colorado 50 Ha Plot." *Forest Ecosystems* 4: 17. <https://doi.org/10.1186/s40663-017-0103-1>

Condit, R., R. Pérez, S. Aguilar, S. Lao, R. Foster, and S. P. Hubbell. 2019a. "Complete Data from the Barro Colorado 50-Ha Plot: 423617 Trees, 35 Years." Dryad, dataset. <https://doi.org/10.15146/5xcp-0d46>

Condit, R., R. Pérez, S. Aguilar, S. Lao, R. Foster, and S. P. Hubbell. 2019b. "BCI 50-Ha Plot Taxonomy." Dryad, dataset. <https://doi.org/10.15146/R3FH61>

Cushman, K. C., H. C. Muller-Landau, R. S. Condit, and S. P. Hubbell. 2014. "Improving Estimates of Biomass Change in Buttressed Trees Using Tree Taper Models." *Methods in Ecology and Evolution* 5: 573–582. <https://doi.org/10.1111/2041-210X.12187>

Davies, S. J., I. Abiem, K. Abu Salim, S. Aguilar, D. Allen, A. Alonso, K. Anderson-Teixeira, et al. 2021. "ForestGEO: Understanding Forest Diversity and Dynamics through a Global Observatory Network." *Biological Conservation* 253: 108907. <https://doi.org/10.1016/j.biocon.2020.108907>

Feeley, K. J., S. J. Davies, R. Pérez, S. P. Hubbell, and R. B. Foster. 2011. "Directional Changes in the Species Composition of a Tropical Forest." *Ecology* 92: 871–882. <https://doi.org/10.1890/10-0724.1>

Harms, K. E., R. Condit, S. P. Hubbell, and R. B. Foster. 2001. "Habitat Associations of Trees and Shrubs in a 50-Ha Neotropical Forest Plot." *Journal of Ecology* 89: 947–959. <https://doi.org/10.1046/j.0022-0477.2001.00615.x>

Hubbell, S. P., R. Condit, and R. B. Foster. 1990. "Presence and Absence of Density Dependence in a Neotropical Tree Community." *Philosophical Transactions of the Royal Society of London, Series B* 330: 269–281. <https://doi.org/10.1098/rstb.1990.0198>

Hubbell, S. P., and R. B. Foster. 1983. "Diversity of Canopy Trees in a Neotropical Forest and Implications for Conservation." In *Tropical Rain Forest: Ecology and Management*, edited by S. L. Sutton, T. C. Whitmore, and A. C. Chadwick, 25–41. Oxford: Blackwell Scientific Publications.

Hubbell, S. P., R. B. Foster, S. T. O'Brien, K. E. Harms, R. Condit, B. Wechsler, S. J. Wright, and S. Loo de Lao. 1999. "Light-Gap Disturbances, Recruitment Limitation, and Tree Diversity in a Neotropical Forest." *Science* 283: 554–557. <https://doi.org/10.1126/science.283.5401.554>

Hubbell, S. P., R. B. Foster, and I. Rubinoff. 2024. "Early History of the Barro Colorado Island Forest Dynamics Project." In *The First 100 Years of Research on Barro Colorado: Plant and Ecosystem Science*, edited by H. C. Muller-Landau and S. J. Wright. Washington, DC: Smithsonian Institution Scholarly Press. <https://doi.org/10.5479/si.26809612>

Kindt, R. 2020. "WorldFlora: An R Package for Exact and Fuzzy Matching of Plant Names against the World Flora Online Taxonomic Backbone Data." *Applications in Plant Sciences* 8: e11388. <https://doi.org/10.1002/aps3.11388>

Leigh, E. G., Jr. 1999. *Tropical Forest Ecology: A View from Barro Colorado Island*. New York: Oxford University Press.

Muller-Landau, H. C., and S. J. Wright, eds. 2024. *The First 100 Years of Research on Barro Colorado: Plant and Ecosystem Science*. 2 vols. Washington, DC: Smithsonian Institution Scholarly Press. <https://doi.org/10.5479/si.26048527>

Piponiot, C., K. J. Anderson-Teixeira, S. J. Davies, D. Allen, N. A. Bourg, D. F. R. P. Burslem, D. Cárdenas, et al. 2022. "Distribution of Biomass Dynamics in Relation to Tree Size in Forests across the World." *New Phytologist* 234: 1664–1677. <https://doi.org/10.1111/nph.17995>

Piponiot, C., R. Condit, S. P. Hubbell, R. Pérez, S. Lao, S. Aguilar, and H. C. Muller-Landau. 2024. "Woody Biomass Stocks and Fluxes in the Barro Colorado Island 50-ha Plot." In *The First 100 Years of Research on Barro Colorado: Plant and Ecosystem Science*, edited by H. C. Muller-Landau and S. J. Wright. Washington, DC: Smithsonian Institution Scholarly Press. <https://doi.org/10.5479/si.26882404>

Rüger, N., L. S. Comita, R. Condit, D. Purves, B. Rosenbaum, M. D. Visser, S. J. Wright, and C. Wirth. 2018. "Beyond the Fast–Slow Continuum: Demographic Dimensions Structuring a Tropical Tree Community." *Ecology Letters* 21: 1075–1084. <https://doi.org/10.1111/ele.12974>

Rüger, N., R. Condit, D. H. Dent, S. J. DeWalt, S. P. Hubbell, J. W. Lichstein, O. R. Lopez, C. Wirth, and C. E. Farrior. 2020. "Demographic Trade-Offs Predict Tropical Forest Dynamics." *Science* 368: 165–168. <https://doi.org/10.1126/science.aaz4797>

Rutishauser, E., S. J. Wright, R. Condit, S. P. Hubbell, S. J. Davies, and H. C. Muller-Landau. 2020. "Testing for Changes in Biomass Dynamics in Large-Scale Forest Datasets." *Global Change Biology* 26: 1485–1498. <https://doi.org/10.1111/gcb.14833>

WFO (World Flora Online Consortium). 2026. "World Flora Online Plant List, June 2026" (version 2026-06). Zenodo. <https://doi.org/10.5281/zenodo.20782718>
