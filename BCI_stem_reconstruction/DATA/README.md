# DATA: pending tasks

Last updated 2026-10-16. Run in use: `20261014_220452_unknown_allT_DP_MB_NME_R12_BD_LT_CU_g5_sm0p5_kg0_ks0_rcpp`.

## Folders

| Folder | Content |
|---|---|
| `RAW/` | database exports (`oct_13_2026`, `oct_16_2026`), species sources (`sp_tables/`), Condit et al. (2019) |
| `SPP_TABLE/` | stage-1 species table, Condit's table, comparison helper |
| `PROCESSED/` | stage-1 stem table, engine input |
| `<run folder>/` | stage-2 engine output and merge |
| `POSTERIORS/`, `RTABLES/` | stage-3 posterior preparation and R tables |
| `CHECKS/` | check files; `revision_oct16/` holds the Oct 16 review lists |
| `BACKUPS/` | copy of the run of 2026-10-05 |

## To do: data manager (Oct 16 export)

Details and CSVs: `CHECKS/revision_oct16/manejador_de_datos/`.

- [ ] swars1 / swars2: variety labels shifted (swars1 must be *grandiflora*, swars2 *continentalis*).
- [ ] protte (SubspeciesID 15) and quaras (SubspeciesID 13): duplicate taxonomy rows, ID without a name.
- [ ] appuse, pipeco, sapisp: names vs Flora de Panamá (wait for the botanists on appuse and pipeco).
- [ ] 48 codes with old or mixed families (Fabaceae subfamilies, Myrsinaceae, Clusiaceae, turpoc Tapisciaceae).
- [ ] 6 stems with two measurements at the same HOM and date.
- [ ] 3 stems with two StemIDs for one StemTag (trees 29078, 315977, 155532); tree 236797 stem 1a dead beside live stems.
- [ ] 1,123 stems and 1,267 trees dead and then alive; 34 dead records with DBH; 557 alive without DBH.
- [ ] 174 probable DBH typos (tree 1114: 2375 mm in 1985); 189 stems losing more than half their DBH.
- [ ] 177 records without date; 98 dates outside their census; tree 73712 stem 5 dated 2013.
- [ ] 86 trees without coordinates; 4 trees with two locations; tree 550132 quadrat vs coordinates.
- [ ] StemNumber (0 or empty from 2015), PrimaryStem empty from 2015, LargeStem empty, HOM on records without DBH.
- [ ] Question: should the partial census 9 (2020) stay in the export?

## To do: botanists (final species list)

Details: `CHECKS/revision_oct16/botanicos/`.

- [ ] Spelling: appuse *seibertii*, ardiba *bartlettii*, rinosy *sylvatica* (Seem.) Kuntze, malmsp *garwoodiae*; pipeco *trinervium* vs *trinerve*.
- [ ] guargr, laetth, swars1, swars2: WFO treats them as synonyms.
- [ ] Synonym column: swars1 synonym equals its own name; 20 entries WFO assigns to other species (consider a "previous names" column); hassfl *Hasseltia panamensis* not found.
- [ ] Authorities: thevah, anacex.
- [ ] Infraspecific rank: protte (subsp. *sessiliflorum*), myrcga ("tip. *gatunensis*").
- [ ] Morphospecies: apeihy, nects1, nects3, sapisp; uniform labels; lifeform for apeihy, nects1, nects3.
- [ ] Empty fields: common name (19 codes), herbarium (76 codes).

## To do: pipeline

- [ ] Paste the spelling lines into `1_DATA_PREPARATION/1_prepare_species_tables.R` section 9c/9e once the botanists answer (appuse, ardiba, rinosy, malmsp).
- [ ] Decide whether to switch the pipeline to the Oct 16 export (needs path changes in stage-1 scripts 0 to 2 and a relaxed same-codes check in script 1, section 2, because pterof is gone).
- [ ] Merge: restore the raw ExactDate on the 16 records filled by stage 2.
- [ ] Project downstream (stage 4 and the publication folder) from run 20261014_220452.
- [ ] Delete `COMPARE_DATASETS/bci.stem.mine/` (unused).
- [ ] Push the local commits.

## Fixed rules

- `DFstatus` is never modified; it stays equal to the raw Status.
- The 2005 duplicate records are not a problem; ignore them.
- Flora de Panamá names have priority; WFO is only a spelling check.
