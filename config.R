# Relative project configuration for the network analyses.
project_root <- normalizePath(
  Sys.getenv("NETWORK_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)
analysis_drug <- Sys.getenv("ANALYSIS_DRUG", unset = "tofacitinib")
input_dir <- file.path(project_root, "data", "input")
output_dir <- file.path(project_root, "results", "network")
n_bootstrap <- 1000L
random_seed <- 42L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

