# Direct-U mean/eigenfunctions for estimator (13), with fitted cumulative FVE.

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
output_filename <- if (tolower(sex) == "female") "mean_IN_F.pdf" else "mean_IN_M.pdf"

# Proposed-method settings
kappa <- 500
latent <- FALSE
method_xi <- "IN"
project_baselines <- FALSE
u_cap <- 0.999999
z_cap <- NULL

# Plot styling
plot_width <- 500 / 72
plot_height <- 400 / 72
use_cairo_pdf <- TRUE
use_legacy_fve_labels <- FALSE
legacy_fve_labels <- c("93.76%", "98.21%", "99.53%")
eigen_sign_flip <- c(-1, 1, -1)
theme_name <- "bw"
show_panel_border <- FALSE
show_axis_lines <- FALSE
axis_line_color <- "black"
axis_line_width <- 0.4
y_limits <- c(-0.5, 0.5)
mean_color <- "black"
eigen_colors <- c("green", "red", "blue")
mean_line_width <- 1.2
eigen_line_width <- 1.2
base_family <- "sans"
axis_title_size <- 14
axis_text_size <- 14
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

mean_eigen_result <- run_mean_eigen_figure(
  context = context,
  output_dir = output_dir,
  output_filename = output_filename,
  kappa = kappa,
  latent = latent,
  method_xi = method_xi,
  project_baselines = project_baselines,
  u_cap = u_cap,
  z_cap = z_cap,
  use_legacy_fve_labels = use_legacy_fve_labels,
  legacy_fve_labels = legacy_fve_labels,
  eigen_sign_flip = eigen_sign_flip,
  theme_name = theme_name,
  show_panel_border = show_panel_border,
  show_axis_lines = show_axis_lines,
  axis_line_color = axis_line_color,
  axis_line_width = axis_line_width,
  mean_color = mean_color,
  eigen_colors = eigen_colors,
  mean_line_width = mean_line_width,
  eigen_line_width = eigen_line_width,
  y_limits = y_limits,
  base_family = base_family,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size,
  legend_title_size = legend_title_size,
  legend_text_size = legend_text_size,
  legend_position = legend_position,
  plot_width = plot_width,
  plot_height = plot_height,
  save_plot = save_plot,
  use_cairo_pdf = use_cairo_pdf,
  display_plot = display_plot
)

mean_eigen_plot <- mean_eigen_result$plot
mean_eigen_fit <- mean_eigen_result$fit
mean_eigen_output_path <- mean_eigen_result$output_path
