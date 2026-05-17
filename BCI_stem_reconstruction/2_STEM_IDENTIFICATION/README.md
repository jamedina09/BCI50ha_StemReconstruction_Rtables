# 2_STEM_IDENTIFICATION

Stage 2 of the BCI stem-reconstruction pipeline. Runs the DP-global stem
matching engine on the cleaned ViewFullTable from `1_DATA_PREPARATION/` and
assembles the per-stem reconstructed identities used to build the final
R tables.

## Scripts (run in order)

- `1_main_cpp_chunk_bci.R` — Driver for the chunked DP run. Splits the cleaned
  ViewFullTable into manageable chunks (by Tag/quadrat), invokes the Rcpp
  transition-cost engine for each chunk, and writes per-chunk outputs.
- `2_merge_chunks_to_datatable.R` — Concatenates the per-chunk outputs into a
  single data.table, deduplicates, and reconciles posterior columns across
  chunks.
- `3_initial_comparisson_dryad_and_me.R` — Diagnostic comparison of the
  reconstructed StemIDs against the published Dryad reference dataset
  (per-method DBH agreement, status agreement, method-type breakdown).

## Other files

- `run_chunk_bci.txt` — Notes / command lines for launching long chunked runs.
- `comparissons/` — Local comparison outputs (gitignored).
