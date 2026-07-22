# Reproduces the female transport-sign heatmap.
# Legacy reference: old code/Country_plot_new.R -> old figure/countries_F.pdf
# Paper/supplement reference: DP-JRSSB-RR1/Fig/contry_plot_f.pdf

# =========================
# Parameters
# =========================

project_root <- NULL
helper_file <- file.path("mortality_repro_functions.R")
sex <- "Female"
output_dir <- file.path("output")
output_filename <- if (tolower(sex) == "female") "countries_F.pdf" else "countries_M.pdf"

# Use the row order from the saved legacy figure for a closer visual match.
ordered_country_codes <- if (tolower(sex) == "female") {
  c(
    "JPN", "SWE", "ITA", "ESP", "NOR", "ISL", "CHE", "NLD", "GRC",
    "AUS", "FRATNP", "AUT", "LUX", "CAN", "FIN", "SVN", "GBRTENW",
    "PRT", "IRL", "BEL", "NZL_NP", "DNK", "CZE", "TWN", "POL", "SVK",
    "USA", "EST", "BGR", "HUN", "LTU", "BLR", "LVA"
  )
} else {
  NULL
}

# Plot styling
plot_width <- 500 / 72
plot_height <- 600 / 72
use_cairo_pdf <- TRUE
fill_low <- "red"
fill_high <- "blue"
show_legend <- FALSE
x_break_years <- seq(1983, 2018, by = 5)
base_family <- "sans"
axis_text_x_size <- 14
axis_text_y_size <- 14
axis_text_angle <- 90
legend_title_size <- 12
legend_text_size <- 12
legend_position <- "right"

# Workflow toggles
save_plot <- TRUE
display_plot <- interactive()

# =========================
# Run
# =========================

source(helper_file)

context <- create_mortality_repro_context(project_root = project_root, sex = sex)

country_sign_result <- run_country_sign_figure(
  context = context,
  output_dir = output_dir,
  output_filename = output_filename,
  ordered_country_codes = ordered_country_codes,
  fill_low = fill_low,
  fill_high = fill_high,
  show_legend = show_legend,
  x_break_years = x_break_years,
  base_family = base_family,
  axis_text_x_size = axis_text_x_size,
  axis_text_y_size = axis_text_y_size,
  axis_text_angle = axis_text_angle,
  legend_title_size = legend_title_size,
  legend_text_size = legend_text_size,
  legend_position = legend_position,
  plot_width = plot_width,
  plot_height = plot_height,
  save_plot = save_plot,
  use_cairo_pdf = use_cairo_pdf,
  display_plot = display_plot
)

country_sign_plot <- country_sign_result$plot
country_sign_data <- country_sign_result$heatmap_data
country_sign_output_path <- country_sign_result$output_path
