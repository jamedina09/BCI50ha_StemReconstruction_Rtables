# Makefile - convenience targets
.PHONY: smoke test test-bb

smoke:
	@echo "Smoke test: verifying dp_global R modules load without errors..."
	@Rscript -e "source('dp_global/R/dp_global_main.R'); cat('[smoke] All modules loaded OK\n')"
