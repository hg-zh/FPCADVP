# Main convenience runner for the female mortality figures.
#
# This script produces:
# 1. `01_mean_eigenfunctions_female.R` -> `mean_F.pdf`
# 2. `02_country_sign_heatmap_female.R` -> `countries_F.pdf`
# 3. `03_trajectory_html_female.R` -> `bel_F.html`, `blr_F.html` by default
#
# Run this script in RStudio with Source, or run the individual scripts if you
# only want one figure at a time.

# =========================
# Parameters
# =========================

project_root <- NULL
helper_file <- file.path("repro_rstudio_scripts", "mortality_repro_functions.R")
sex <- "Female"
output_dir <- file.path("repro_rstudio_scripts", "output")

# Which outputs to build
build_mean_eigen_plot <- TRUE
build_country_sign_plot <- TRUE
build_trajectory_html <- TRUE

# Shared workflow toggles
use_cairo_pdf <- TRUE
display_plots <- interactive()
display_widgets <- interactive()
trajectory_file_map <- c(
  ESP = paste0("bel_", if (tolower(sex) == "female") "F" else "M", ".html"),
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
