# Makefile - convenience targets
.PHONY: smoke profile docs

smoke:
	@echo "Running smoke tests (uses run_dp_future.R DRY_RUN)..."
	@./bin/run_dp_future.R --workers 1 --cores-per-job 1 --configs "fixed" -- --DRY_RUN
	@Rscript dp_global/scripts/run_smoke.R

profile:
	@echo "Profiling (cpp variant)..."
	@PROFILE_VARIANT=cpp Rscript dp_global/dev/profiling_code.R

docs:
	@echo "Building docs (requires pandoc + MathJax)..."
	@if command -v pandoc >/dev/null; then \
	  pandoc -s --from=markdown+tex_math_dollars+raw_tex --mathjax dp_global/README.md -o docs/dp_global.html; \
	else \
	  echo "pandoc not found; install pandoc to build docs."; \
	fi
