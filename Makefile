# Makefile - convenience targets
.PHONY: smoke profile docs

smoke:
	@echo "Running smoke tests (uses run_dp_future.R DRY_RUN)..."
	@./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
	@Rscript dp_global/scripts/run_smoke.R

profile:
	@echo "Profiling (cpp-only)..."
	@Rscript --vanilla dp_global/dev/profiling_code.R