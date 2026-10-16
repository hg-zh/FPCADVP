# Convenience runner for direct-U IN female mortality figures (estimator 13).
#
# This script produces:
# 1. `01_mean_eigenfunctions_female.R` -> `mean_IN_F.pdf`
# 2. `02_country_sign_heatmap_female.R` -> `countries_F.pdf`
# 3. `03_trajectory_html_female.R` -> `esp_F.html`, `blr_F.html` by default
#
# Run this script in RStudio with Source, or run the individual scripts if you
# only want one figure at a time.

# =========================
# Parameters
# =========================

# Resolve adjacent files from the entrypoint, independently of the working directory.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file)) {
  script_file <- sub("^--file=", "", script_file[1])
} else {
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(script_file) && requireNamespace("rstudioapi", quietly = TRUE)) {
    script_file <- rstudioapi::getSourceEditorContext()$path
  }
}
if (is.null(script_file) || !nzchar(script_file)) stop("Run this file with Rscript or Source.")
project_root <- dirname(normalizePath(gsub("~+~", " ", script_file, fixed = TRUE), mustWork = TRUE))
helper_file <- file.path(project_root, "mortality_repro_functions.R")
sex <- "Female"
output_dir <- file.path(dirname(helper_file), "output", "estimator13_IN")

# Which outputs to build
build_mean_eigen_plot <- TRUE
build_country_sign_plot <- TRUE
build_trajectory_html <- TRUE

# Shared workflow toggles
use_cairo_pdf <- TRUE
display_plots <- interactive()
display_widgets <- interactive()
trajectory_file_map <- c(
  ESP = paste0("esp_", if (tolower(sex) == "female") "F" else "M", ".html"),
  BLR = paste0("blr_", if (tolower(sex) == "female") "F" else "M", ".html")
)
legacy_heatmap_order <- c(
  "JPN", "SWE", "ITA", "ESP", "NOR", "ISL", "CHE", "NLD", "GRC",
  "AUS", "FRATNP", "AUT", "LUX", "CAN", "FIN", "SVN", "GBRTENW",
  "PRT", "IRL", "BEL", "NZL_NP", "DNK", "CZE", "TWN", "POL", "SVK",
  "USA", "EST", "BGR", "HUN", "LTU", "BLR", "LVA"
)

# =========================
# Run
# =========================

source(helper_file)

context <- create_mortality_repro_context(project_root = project_root, sex = sex)

if (isTRUE(build_mean_eigen_plot)) {
  run_mean_eigen_figure(
    context = context,
    output_dir = output_dir,
    use_cairo_pdf = use_cairo_pdf,
    display_plot = display_plots
  )
}

if (isTRUE(build_country_sign_plot)) {
  run_country_sign_figure(
    context = context,
    output_dir = output_dir,
    ordered_country_codes = legacy_heatmap_order,
    use_cairo_pdf = use_cairo_pdf,
    display_plot = display_plots
  )
}

if (isTRUE(build_trajectory_html)) {
  run_trajectory_html_figures(
    context = context,
    output_dir = output_dir,
    country_codes = names(trajectory_file_map),
    output_filenames = trajectory_file_map,
    legacy_baseline_postprocess = TRUE,
    display_widgets = display_widgets
  )
}
