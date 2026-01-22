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

# Usage:
# 1. Run this script to create dpglobal_bundle.RData
# 2. On another system, load("dpglobal_bundle.RData") to access all functions/objects
