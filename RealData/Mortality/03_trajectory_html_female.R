# Reproduces the female 3D transport trajectory HTML files.
# Legacy reference: old code/Individual_illu.R -> old figure/bel_F.html and old figure/blr_F.html
#
# Note:
# In the legacy script, `bel_F.html` was saved from `Pos_fig[[3]]`, which
# corresponds to Spain (`ESP`) rather than Belgium (`BEL`). The defaults below
# mirror the archived HTML files exactly. Edit `legacy_output_country_map` if
# you want filenames that match the country codes instead.

# =========================
# Parameters
# =========================

project_root <- NULL
helper_file <- file.path("mortality_repro_functions.R")
sex <- "Female"
output_dir <- file.path("output")

# Countries and output files
sex_tag <- if (tolower(sex) == "female") "F" else "M"
legacy_output_country_map <- c(
  ESP = paste0("bel_", sex_tag, ".html"),
  BLR = paste0("blr_", sex_tag, ".html")
)
country_codes <- names(legacy_output_country_map)
output_filenames <- legacy_output_country_map

# Proposed-method settings
kappa <- 500
latent <- FALSE
method_xi <- "CE"
project_baselines <- FALSE
u_cap <- 0.999999
z_cap <- NULL
include_mean <- FALSE
score_method <- "xi_est"
legacy_baseline_postprocess <- TRUE
project_reconstructed_rows <- FALSE

# Plot styling
observed_line_color <- "red"
observed_line_width <- 3
observed_line_opacity <- 0.6
fitted_surface_color <- "blue"
fitted_surface_opacity <- 0.6
reference_surface_color <- "gray"
reference_surface_opacity <- 0.5
camera_center <- c(x = 0, y = 0, z = -0.2)
camera_eye <- c(x = 1.25, y = 1.25, z = 1.25)
x_axis_title <- "Calendar Year"
y_axis_title <- "Age"
z_axis_title <- " "
plot_title <- NULL
selfcontained <- FALSE

# Workflow toggles
save_widgets <- TRUE
display_widgets <- interactive()

# =========================
# Run
# =========================

source(helper_file)

context <- create_mortality_repro_context(project_root = project_root, sex = sex)

trajectory_result <- run_trajectory_html_figures(
  context = context,
  output_dir = output_dir,
  country_codes = country_codes,
  output_filenames = output_filenames,
  kappa = kappa,
  latent = latent,
  method_xi = method_xi,
  project_baselines = project_baselines,
  u_cap = u_cap,
  z_cap = z_cap,
  include_mean = include_mean,
  score_method = score_method,
  legacy_baseline_postprocess = legacy_baseline_postprocess,
  project_reconstructed_rows = project_reconstructed_rows,
  observed_line_color = observed_line_color,
  observed_line_width = observed_line_width,
  observed_line_opacity = observed_line_opacity,
  fitted_surface_color = fitted_surface_color,
  fitted_surface_opacity = fitted_surface_opacity,
  reference_surface_color = reference_surface_color,
  reference_surface_opacity = reference_surface_opacity,
  camera_center = camera_center,
  camera_eye = camera_eye,
  x_axis_title = x_axis_title,
  y_axis_title = y_axis_title,
  z_axis_title = z_axis_title,
  plot_title = plot_title,
  selfcontained = selfcontained,
  save_widgets = save_widgets,
  display_widgets = display_widgets
)

trajectory_fit <- trajectory_result$fit
trajectory_scalar_matrix <- trajectory_result$scalar_matrix
trajectory_widgets <- trajectory_result$widgets
trajectory_output_paths <- trajectory_result$output_paths
