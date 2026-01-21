############################################################
# load_dp_functions.R
# Helpers to load dp_global R modules into a single object (environment/list)
# - load_dp_functions(): source project modules into a dedicated environment and
#   return a list containing all functions and metadata.
# - save_dp_functions(): save the returned object to an RDS file for reuse.
# - load_dp_functions_from_rds(): read the RDS back into R.
# - attach_dp_functions(): attach the functions into an environment (e.g., .GlobalEnv)
#
# Usage example:
# loader <- load_dp_functions()
# loader$functions$dp_global_main(...)    # call a function
# save_dp_functions(loader, "dp_global_functions.rds")

load_dp_functions <- function(project_root = here::here(),
                              r_files = c(
                                  "dp_global_utils.R",
                                  "dp_global_bio.R",
                                  "dp_global_states.R",
                                  "dp_global_matchers.R",
                                  "dp_global_dp.R",
                                  "dp_global_diag.R"
                              ),
                              required_symbols = c(
                                  "estimate_bio_pars",
                                  "transition_cost_tracks_bio_components",
                                  "match_stems_dp_global_backward_marginals_batch",
                                  "match_stems_optimal_backward",
                                  "enumerate_states_injective",
                                  "add_dp_posterior_bins"
                              ),
                              parent = baseenv()) {
    # Minimal checks
    if (!requireNamespace("here", quietly = TRUE)) stop("Please install 'here' package")

    files_full <- file.path(project_root, "dp_global", "R", r_files)
    for (fp in files_full) {
        if (!file.exists(fp)) stop(sprintf("Required file not found: %s", fp), call. = FALSE)
    }

    env <- new.env(parent = parent)
    symbol_map <- list() # name -> file

    for (i in seq_along(files_full)) {
        fp <- files_full[[i]]
        existing_before <- ls(envir = env, all.names = TRUE)
        # Use sys.source so top-level definitions go into `env` and the file runs
        # exactly as if executed at top-level of that environment.
        sys.source(fp, envir = env)
        existing_after <- ls(envir = env, all.names = TRUE)
        new_objs <- setdiff(existing_after, existing_before)
        if (length(new_objs) > 0L) {
            # Record which file created each new top-level symbol
            for (nm in new_objs) symbol_map[[nm]] <- basename(fp)
        }
    }

    # Collect functions only
    all_names <- ls(envir = env, all.names = TRUE)
    funcs <- Filter(is.function, mget(all_names, envir = env))

    missing <- setdiff(required_symbols, names(funcs))
    if (length(missing) > 0L) {
        warning(sprintf("Some required symbols were not found after loading: %s", paste(missing, collapse = ", ")))
    }

    list(
        functions = funcs,
        env = env,
        files_loaded = r_files,
        symbol_map = symbol_map,
        required_symbols_missing = missing
    )
}

save_dp_functions <- function(obj, file = "sourced_dp_functions.rds") {
    if (is.list(obj) && !is.null(obj$functions)) {
        saveRDS(obj, file)
        invisible(file)
    } else {
        stop("Object must be the loader list returned by load_dp_functions().")
    }
}

load_dp_functions_from_rds <- function(file) {
    readRDS(file)
}

attach_dp_functions <- function(loader_obj, envir = .GlobalEnv) {
    if (!is.list(loader_obj) || is.null(loader_obj$functions)) stop("loader_obj must be the list returned by load_dp_functions()")
    list2env(loader_obj$functions, envir = envir)
    invisible(envir)
}

# End of file

loader <- load_dp_functions()

save_dp_functions(loader, "./dp_global/R/load_dp_functions.rdata")