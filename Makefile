# Makefile - convenience targets
.PHONY: smoke test test-bb

smoke:
	@echo "Smoke test: verifying dp_global R modules load without errors..."
	@Rscript -e "source('dp_global/R/dp_global_main.R'); cat('[smoke] All modules loaded OK\n')"

test-bb:
	@echo "Running broken-below invariant tests..."
	@Rscript dp_global/tests/test_broken_below_invariants.R

test: smoke test-bb