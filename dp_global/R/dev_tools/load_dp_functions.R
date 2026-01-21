############################################################
# load_dp_functions.R
# Minimal helper to load all `dp_global` R functions into an object and
# export them as an `.rdata` file that you can load elsewhere.
#
# Behavior:
# - Sources core R modules into a dedicated environment and returns a named
#   list of functions.
# - Attempts to compile and load the C++ acceleration via Rcpp::sourceCpp()
#   if available; leaves a graceful warning otherwise.
# - Writes a file `dp_functions.rdata` (containing `dp_functions`) by default.
#
# Usage (from project root):
# source("dp_global/R/load_dp_functions.R")
# # This script will load modules and write dp_global/R/dp_functions.rdata

load_dp_functions <- function(project_root = here::here(),
                              r_files = c(
                                  "dp_global_utils.R",
                                  "dp_global_bio.R",
                                  "dp_global_states.R",
                                  "dp_global_matchers.R",
                                  "dp_global_dp.R",
                                  "dp_global_diag.R"
                              ),
                              try_compile_cpp = TRUE,
                              out_file = file.path("dp_global", "R", "dp_functions.rdata"),
                              env_parent = baseenv()) {
    if (!requireNamespace("here", quietly = TRUE)) stop("Please install the 'here' package")

    files_full <- file.path(project_root, "dp_global", "R", r_files)
    for (fp in files_full) {
        if (!file.exists(fp)) stop(sprintf("Required file not found: %s", fp), call. = FALSE)
    }

    env <- new.env(parent = env_parent)

    # Source R modules into the dedicated env
    for (fp in files_full) {
        sys.source(fp, envir = env)
    }

    # Rcpp/C++ acceleration: try to source the R wrapper (so it exists in env)
    # and capture the C++ source so it can be saved into the .rdata for later
    # re-compilation when loaded elsewhere.
    cpp_wrapper <- file.path(project_root, "dp_global", "src", "transition_cost_rcpp.R")
    cpp_file <- file.path(project_root, "dp_global", "src", "transition_cost_rcpp.cpp")

    if (file.exists(cpp_wrapper)) {
        sys.source(cpp_wrapper, envir = env)
    }

    cpp_source <- NULL
    cpp_filename <- NULL
    compiled_ok <- FALSE
    if (file.exists(cpp_file)) {
        cpp_source <- readLines(cpp_file, warn = FALSE)
        cpp_filename <- basename(cpp_file)
    }

    if (isTRUE(try_compile_cpp)) {
        if (requireNamespace("Rcpp", quietly = TRUE) && file.exists(cpp_file)) {
            tryCatch({
                Rcpp::sourceCpp(cpp_file)
                compiled_ok <- TRUE
                message("[load_dp_functions] Rcpp compiled and loaded: transition_cost_rcpp.cpp")
            }, error = function(e) {
                warning(sprintf("Rcpp compilation failed: %s", e$message))
            })
        } else if (!requireNamespace("Rcpp", quietly = TRUE)) {
            warning("Rcpp not installed; compiled acceleration will not be available.")
        }
    }

    # Collect functions from the env
    all_names <- ls(envir = env, all.names = TRUE)
    funcs <- Filter(is.function, mget(all_names, envir = env))

    list(
        functions = funcs,
        env = env,
        out_file = out_file,
        cpp_source = cpp_source,
        cpp_filename = cpp_filename,
        compiled_ok = compiled_ok
    )
}

# Save loader functions and C++ source as an .rdata file for easy loading elsewhere
save_dp_functions_rdata <- function(loader_obj, file = NULL) {
    if (!is.list(loader_obj) || is.null(loader_obj$functions)) stop("loader_obj must be the list returned by load_dp_functions()")
    file <- file %||% loader_obj$out_file
    dp_functions <- loader_obj$functions
    cpp_source <- loader_obj$cpp_source
    cpp_filename <- loader_obj$cpp_filename
    compiled_ok <- loader_obj$compiled_ok
    save(dp_functions, cpp_source, cpp_filename, compiled_ok, file = file)
    invisible(file)
}

# Utility: load dp functions .rdata and optionally recompile the included C++ source
load_dp_functions_rdata <- function(file, compile_cpp = TRUE, cleanup = TRUE) {
    if (!file.exists(file)) stop("File not found: ", file)
    e <- new.env()
    load(file = file, envir = e)
    # Expect dp_functions in the file
    if (!exists("dp_functions", envir = e)) stop("dp_functions not found in rdata file")
    dp_functions <- get("dp_functions", envir = e)

    # If C++ source is present and compilation requested, write and compile
    if (isTRUE(compile_cpp) && exists("cpp_source", envir = e) && !is.null(e$cpp_source)) {
        if (!requireNamespace("Rcpp", quietly = TRUE)) {
            warning("Rcpp not installed; cannot compile embedded C++ source.")
        } else {
            tf <- tempfile(pattern = "transition_cost_rcpp_", fileext = ".cpp")
            writeLines(e$cpp_source, con = tf)
            tryCatch({
                Rcpp::sourceCpp(tf)
                message("[load_dp_functions_rdata] Compiled embedded C++ source into session.")
            }, error = function(err) {
                warning(sprintf("Failed to compile embedded C++ source: %s", err$message))
            })
            if (isTRUE(cleanup) && file.exists(tf)) unlink(tf)
        }
    }

    list(functions = dp_functions, cpp_source = if (exists("cpp_source", envir = e)) e$cpp_source else NULL,
         cpp_filename = if (exists("cpp_filename", envir = e)) e$cpp_filename else NULL,
         compiled_ok = if (exists("compiled_ok", envir = e)) e$compiled_ok else FALSE)
}

# Minimal utility: infix default for NULL
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper: load .rdata and attach functions into an environment (optionally compile C++ source)
attach_dp_functions_from_rdata <- function(file = file.path("dp_global", "R", "dp_functions.rdata"),
                                           envir = .GlobalEnv,
                                           compile_cpp = TRUE,
                                           cleanup = TRUE) {
    obj <- load_dp_functions_rdata(file, compile_cpp = compile_cpp, cleanup = cleanup)
    if (!is.null(obj$functions) && length(obj$functions) > 0L) {
        list2env(obj$functions, envir = envir)
        message(sprintf("[attach_dp_functions_from_rdata] Attached %d functions into %s", length(obj$functions), deparse(substitute(envir))))
    } else {
        message("[attach_dp_functions_from_rdata] No functions found in rdata")
    }
    invisible(obj)
}

# Run when sourced as a script: load modules and write the .rdata file into project
loader <- load_dp_functions()
save_dp_functions_rdata(loader)
message(sprintf("[load_dp_functions] Wrote: %s", loader$out_file))


