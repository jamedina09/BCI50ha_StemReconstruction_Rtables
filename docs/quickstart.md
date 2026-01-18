Quickstart
==========

1. Install prerequisites:
   - R (with packages required by the project).

2. Run a dry-run test (safe; does not execute R driver when `--DRY_RUN` provided):

   ```bash
   ./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
   ```

3. For a real run, run `./bin/run_dp_future.R` to execute multiple configs concurrently (recommended). Do not use `bin/run_dp_full_cpp.sh` except for legacy serial usage.

   ```bash
   ./bin/run_dp_future.R --workers 4 --cores-per-job 4 --configs "fixed data_hard"
   ```

Notes
-----
- If you run jobs in parallel manually, export a shared `BATCH_TS` so all jobs write into a single timestamped set of directories under `dp_global/output/`.
- Because `dp_global/output/` can be large and OneDrive may add latency, consider excluding that directory from cloud sync where possible.
