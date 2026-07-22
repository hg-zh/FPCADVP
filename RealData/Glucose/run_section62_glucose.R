#!/usr/bin/env Rscript

# RStudio-friendly Section 6.2 analysis runner.
# Open this file in RStudio and click Source to rerun the full glucose analysis
# with the fixed illustration seed used in the final workflow.
#
# Outputs are written to:
#   newcode/output/section62_glucose_seed_20/

get_this_script_path <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(path)) return(normalizePath(path, mustWork = TRUE))
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  stop("Cannot determine the script path. Save the script and run it from RStudio or via Rscript.")
}

script_path <- get_this_script_path()
script_dir <- dirname(script_path)
project_dir <- script_dir

source(file.path(script_dir, "transport_fpca_modular.R"))
source(file.path(script_dir, "section62_glucose_helpers.R"))

seed <- 20L
data_path <- file.path(project_dir, "data_1.csv")
output_dir <- file.path(script_dir, "output", sprintf("section62_glucose_seed_%s", seed))

results <- glucose_run_section62(
  data_path = data_path,
  seed = seed,
  output_dir = output_dir,
  latent = TRUE
)

cat(sprintf("Finished Section 6.2 analysis for seed %s\n", results$seed))
cat(sprintf("Figure: %s\n", results$figure_file))
cat(sprintf("Results object: %s\n", results$results_file))
cat(sprintf("Selection report: %s\n", results$report_file))
print(results$selected_summary)
