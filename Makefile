# Makefile - convenience targets (run from the project root)
#   make smoke       the dp_global modules load
#   make test-dp     tests of the engine (dp_global/tests/)
#   make test-stage3 tests of the Rstatus rules and the species table; the
#                    Rstatus tests also check the exported R tables of stage 3
#                    (a few minutes)
#   make test        test-dp and test-stage3
# A failing test makes the target fail (run_tests.R).
.PHONY: smoke test test-dp test-stage3

DP_TESTS = dp_global/tests/test_enumerate_states.R \
           dp_global/tests/test_estimate_bio_pars.R \
           dp_global/tests/test_probabilistic_birth_death.R
STAGE3_TESTS = BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_Rstatus.R \
               BCI_stem_reconstruction/3_PREPARE_R_TABLES/tests/test_species_table.R

smoke:
	@echo "Smoke test: verifying dp_global R modules load without errors..."
	@Rscript -e "source('dp_global/R/dp_global_main.R'); cat('[smoke] All modules loaded OK\n')"

test-dp:
	@Rscript run_tests.R $(DP_TESTS)

test-stage3:
	@Rscript run_tests.R $(STAGE3_TESTS)

test: test-dp test-stage3
