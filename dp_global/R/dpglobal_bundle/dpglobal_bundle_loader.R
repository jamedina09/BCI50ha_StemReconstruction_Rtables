# dpglobal_bundle_loader.R
# Loads all functions required by main_cpp.R for DP_GLOBAL workflow


# Source all required R scripts into a new environment
library(here)

dpglobal_env <- new.env()
# 1) Try to run dp_global_main.R to collect module manifest (r_files) if present
tryCatch(
    {
        source(here::here("dp_global", "R", "dp_global_main.R"), local = dpglobal_env)
        message("[dpglobal_bundle_loader] Sourced dp_global_main.R into dpglobal_env")
    },
    error = function(e) {
        warning(sprintf("[dpglobal_bundle_loader] Failed to source dp_global_main.R: %s. Will attempt to source modules directly.", e$message))
    }
)

# 2) If dp_global_main.R defined a manifest `r_files`, use it to source modules in order
if (exists("r_files", envir = dpglobal_env)) {
    rfs <- get("r_files", envir = dpglobal_env)
    for (rf in rfs) {
        fp <- here::here("dp_global", "R", rf)
        if (file.exists(fp)) {
            tryCatch(
                {
                    sys.source(fp, envir = dpglobal_env)
                    message(sprintf("[dpglobal_bundle_loader] Sourced: %s", rf))
                },
                error = function(e) {
                    warning(sprintf("[dpglobal_bundle_loader] Error sourcing %s: %s", rf, e$message))
                }
            )
        } else {
            warning(sprintf("[dpglobal_bundle_loader] File listed in r_files not found: %s", fp))
        }
    }
} else {
    # Fallback: source all R scripts in dp_global/R (excluding bundle dir)
    rn <- list.files(here::here("dp_global", "R"), pattern = "\\.[rR]$", full.names = TRUE)
    rn <- rn[!grepl("dpglobal_bundle", rn)]
    for (fp in rn) {
        tryCatch(
            {
                sys.source(fp, envir = dpglobal_env)
                message(sprintf("[dpglobal_bundle_loader] Sourced: %s", basename(fp)))
            },
            error = function(e) {
                warning(sprintf("[dpglobal_bundle_loader] Error sourcing %s: %s", fp, e$message))
            }
        )
    }
}

# 3) Source the R-wrapper for Rcpp (if present) into the bundle environment
rcpp_wrapper <- here::here("dp_global", "src", "transition_cost_rcpp.R")
if (file.exists(rcpp_wrapper)) {
    tryCatch(
        {
            sys.source(rcpp_wrapper, envir = dpglobal_env)
            message("[dpglobal_bundle_loader] Sourced transition_cost_rcpp.R into dpglobal_env")
        },
        error = function(e) {
            warning(sprintf("[dpglobal_bundle_loader] Failed to source transition_cost_rcpp.R: %s", e$message))
        }
    )
} else {
    message("[dpglobal_bundle_loader] No transition_cost_rcpp.R wrapper found; skipping.")
}

# 4) Save all loaded functions and objects to an RData file
save(list = ls(envir = dpglobal_env), file = here::here("dp_global", "R", "dpglobal_bundle", "dpglobal_bundle.RData"), envir = dpglobal_env)
message(sprintf("[dpglobal_bundle_loader] Saved %d objects to dpglobal_bundle.RData", length(ls(envir = dpglobal_env))))

# 5) Record a small manifest with package requirements and C++ status
bundle_manifest <- list(
    r_files = if (exists("r_files", envir = dpglobal_env)) get("r_files", envir = dpglobal_env) else NULL,
    required_pkgs_by_file = if (exists("required_pkgs_by_file", envir = dpglobal_env)) get("required_pkgs_by_file", envir = dpglobal_env) else NULL,
    attach_pkgs = if (exists("attach_pkgs", envir = dpglobal_env)) get("attach_pkgs", envir = dpglobal_env) else NULL,
    compiled_acceleration_available = (
        exists("transition_cost_tracks_bio_batch_rcpp_cpp", envir = dpglobal_env, mode = "function") ||
        exists("transition_cost_tracks_bio_batch_rcpp_cpp", envir = globalenv(), mode = "function")
    )
)
tryCatch(
    {
        saveRDS(bundle_manifest, file = here::here("dp_global", "R", "dpglobal_bundle", "dpglobal_bundle_manifest.rds"))
        message("[dpglobal_bundle_loader] Wrote dpglobal_bundle_manifest.rds")
    },
    error = function(e) {
        warning(sprintf("[dpglobal_bundle_loader] Failed to write manifest: %s", e$message))
    }
)

# Print quick summary to help users of the bundle
pkgs_needed <- NULL
if (!is.null(bundle_manifest$required_pkgs_by_file)) {
    pkgs_needed <- unique(unlist(bundle_manifest$required_pkgs_by_file))
    pkgs_needed <- pkgs_needed[!is.na(pkgs_needed) & pkgs_needed != ""]
}
if (!is.null(pkgs_needed) && length(pkgs_needed) > 0) {
    message(sprintf("[dpglobal_bundle_loader] Required packages (per manifest): %s", paste(pkgs_needed, collapse = ", ")))
} else {
    message("[dpglobal_bundle_loader] No per-file package requirements detected in the manifest.")
}
if (bundle_manifest$compiled_acceleration_available) {
    message("[dpglobal_bundle_loader] NOTE: C++ acceleration was available in this session. On target systems, run 'install_transition_cost_rcpp.R' to compile the C++ sources before using Rcpp-backed functions.")
} else {
    message("[dpglobal_bundle_loader] NOTE: C++ acceleration was NOT available in this session. If you want faster execution on target systems, run 'install_transition_cost_rcpp.R' after installing 'Rcpp' and a build toolchain.")
}

# Usage:
# 1. Run this script to create dpglobal_bundle.RData and dpglobal_bundle_manifest.rds
# 2. On another system, load("dpglobal_bundle.RData") to access functions/objects
# 3. If compiled acceleration is required, run the helper script `install_transition_cost_rcpp.R` in this directory to compile the C++ sources

