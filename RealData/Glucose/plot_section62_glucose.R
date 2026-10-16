#!/usr/bin/env Rscript

# RStudio-friendly Section 6.2 plotting script.
# Open this file in RStudio and click Source to recreate the final figure from
# the saved seed-20 results object.
#
# This script reads:
#   output/section62_glucose_seed_20/section62_results_seed_20.rds
#   output/section62_glucose_seed_20/baseline_results_seed_20.rds
# and rewrites:
#   output/section62_glucose_seed_20/ind_gluco_seed_20.pdf

get_this_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE), mustWork = TRUE))
  }
  source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_path)) return(normalizePath(source_path, mustWork = TRUE))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(path)) return(normalizePath(path, mustWork = TRUE))
  }
  stop("Run this file with Rscript or Source.")
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
      "\nRun run_section62_glucose.R first."
    )
  )
}
if (!file.exists(baseline_results_path)) {
  stop(
    paste(
      "Saved baseline results not found at",
      baseline_results_path,
      "\nRun run_section62_glucose_baselines.R first."
    )
  )
}

results <- readRDS(results_path)
baseline_results <- readRDS(baseline_results_path)
# Reconstruct the destination after relocation; saved fits may contain old absolute paths.
figure_file <- file.path(dirname(results_path), sprintf("ind_gluco_seed_%s.pdf", results$seed))

# ----------------------------
# Edit these plotting settings
# ----------------------------

# Selection is by the smallest proposed-method mean held-out errors.
# Display the selected pair in subject-index order, as in Figure 6.
selected_for_plot <- results$selected_summary[
  order(results$selected_summary$subject_index), , drop = FALSE
]
PANEL_TITLES <- NULL  # Derive titles from the selected subjects and days.

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
  selected_summary = selected_for_plot,
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

if (interactive()) print(plot_obj)
cat(sprintf("Figure written to: %s\n", figure_file))
