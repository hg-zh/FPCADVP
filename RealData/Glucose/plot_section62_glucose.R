#!/usr/bin/env Rscript

# RStudio-friendly Section 6.2 plotting script.
# Open this file in RStudio and click Source to recreate the final figure from
# the saved seed-20 results object.
#
# This script reads:
#   newcode/output/section62_glucose_seed_20/section62_results_seed_20.rds
#   newcode/output/section62_glucose_seed_20/baseline_results_seed_20.rds
# and rewrites:
#   newcode/output/section62_glucose_seed_20/ind_gluco_seed_20.pdf

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

source(file.path(script_dir, "transport_fpca_modular.R"))
source(file.path(script_dir, "section62_glucose_helpers.R"))

results_path <- file.path(
  script_dir,
  "output",
  "section62_glucose_seed_20",
  "section62_results_seed_20.rds"
)
baseline_results_path <- file.path(
  script_dir,
  "output",
  "section62_glucose_seed_20",
  "baseline_results_seed_20.rds"
)

if (!file.exists(results_path)) {
  stop(
    paste(
      "Saved results not found at",
      results_path,
      "\nRun newcode/run_section62_glucose_seed20.R first."
    )
  )
}
if (!file.exists(baseline_results_path)) {
  stop(
    paste(
      "Saved baseline results not found at",
      baseline_results_path,
      "\nRun newcode/run_section62_glucose_baselines_seed20.R first."
    )
  )
}

results <- readRDS(results_path)
baseline_results <- readRDS(baseline_results_path)
figure_file <- results$files$figure

# ----------------------------
# Edit these plotting settings
# ----------------------------

PANEL_TITLES <- c(
  "Subject 20 (1636-70-1003), t = 3",
  "Subject 20 (1636-70-1003), t = 7",
  "Subject 21 (1636-70-1005), t = 4",
  "Subject 21 (1636-70-1005), t = 5"
)

METHOD_COLORS <- c(
  Observed = "black",
  Proposed = "#1f77b4",
  LQD = "#d95f02",
  Tangent = "#1b9e77"
)

METHOD_LINETYPES <- c(
  Observed = "solid",
  Proposed = "dashed",
  LQD = "dotdash",
  Tangent = "longdash"
)

LINE_WIDTH <- 1. 
ZERO_LINE_WIDTH <- 0.5

AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12
STRIP_TITLE_SIZE <- 14
LEGEND_TEXT_SIZE <- 14
LEGEND_KEY_WIDTH <- 1.5
LEGEND_KEY_WIDTH_UNITS <- "cm"

SHOW_STRIP_TEXT <- FALSE
SHOW_STRIP_BACKGROUND <- FALSE

STRIP_BACKGROUND_FILL <- "white"
STRIP_BACKGROUND_COLOR <- "black"
STRIP_BACKGROUND_LINEWIDTH <- 0.5

SAVE_WIDTH <- 0.8*900/72
SAVE_HEIGHT <- 0.8*300/72
SAVE_UNITS <- "in"

plot_obj <- glucose_plot_method_comparison(
  selected_summary = results$selected_summary,
  observed_holdouts = results$observed_holdouts,
  method_holdouts = list(
    Proposed = results$predicted_holdouts,
    LQD = baseline_results$lqd_holdouts,
    Tangent = baseline_results$tangent_holdouts
  ),
  support_grid = seq(50, 160, length.out = 201),
  out_file = figure_file,
  panel_titles = PANEL_TITLES,
  method_colors = METHOD_COLORS,
  method_linetypes = METHOD_LINETYPES,
  line_width = LINE_WIDTH,
  zero_line_width = ZERO_LINE_WIDTH,
  axis_title_size = AXIS_TITLE_SIZE,
  axis_text_size = AXIS_TEXT_SIZE,
  strip_title_size = STRIP_TITLE_SIZE,
  legend_text_size = LEGEND_TEXT_SIZE,
  legend_key_width = LEGEND_KEY_WIDTH,
  legend_key_width_units = LEGEND_KEY_WIDTH_UNITS,
  show_strip_text = SHOW_STRIP_TEXT,
  show_strip_background = SHOW_STRIP_BACKGROUND,
  strip_background_fill = STRIP_BACKGROUND_FILL,
  strip_background_color = STRIP_BACKGROUND_COLOR,
  strip_background_linewidth = STRIP_BACKGROUND_LINEWIDTH,
  save_width = SAVE_WIDTH,
  save_height = SAVE_HEIGHT,
  save_units = SAVE_UNITS
)

print(plot_obj)
cat(sprintf("Figure written to: %s\n", figure_file))
