# Shared helpers for reproducing the mortality figures in RStudio.

`%||%` <- function(x, y) {
  if (!is.null(x)) {
    x
  } else {
    y
  }
}

assert_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Please install the following package(s) before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

find_project_root <- function(project_root = NULL) {
  if (!is.null(project_root)) {
    return(normalizePath(project_root, winslash = "/", mustWork = TRUE))
  }

  required_dirs <- c("Data")
  candidates <- character(0)
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- c(candidates, current)

  for (i in seq_len(5)) {
    current <- dirname(current)
    candidates <- c(candidates, current)
  }

  for (candidate in unique(candidates)) {
    has_required_dirs <- all(file.exists(file.path(candidate, required_dirs)))
    if (isTRUE(has_required_dirs)) {
      return(candidate)
    }
  }

  stop(
    "Could not detect the Mortality project root from the current working directory.\n",
    "Set `project_root` manually near the top of the script.",
    call. = FALSE
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

country_label_lookup <- function() {
  c(
    AUS = "Australia",
    AUT = "Austria",
    BEL = "Belgium",
    BGR = "Bulgaria",
    BLR = "Belarus",
    CAN = "Canada",
    CHE = "Switzerland",
    CZE = "Czechia",
    DNK = "Denmark",
    ESP = "Spain",
    EST = "Estonia",
    FIN = "Finland",
    FRATNP = "France",
    GBRTENW = "England and Wales",
    GRC = "Greece",
    HUN = "Hungary",
    IRL = "Ireland",
    ISL = "Iceland",
    ITA = "Italy",
    JPN = "Japan",
    LTU = "Lithuania",
    LUX = "Luxembourg",
    LVA = "Latvia",
    NLD = "Netherlands",
    NOR = "Norway",
    NZL_NP = "New Zealand",
    POL = "Poland",
    PRT = "Portugal",
    SVK = "Slovakia",
    SVN = "Slovenia",
    SWE = "Sweden",
    TWN = "Taiwan",
    USA = "USA"
  )
}

country_code_to_label <- function(country_codes) {
  lookup <- country_label_lookup()
  labels <- unname(lookup[country_codes])
  labels[is.na(labels)] <- country_codes[is.na(labels)]
  labels
}

normalize_sex <- function(sex) {
  sex_lower <- tolower(sex)

  if (sex_lower %in% c("female", "f")) {
    return("Female")
  }

  if (sex_lower %in% c("male", "m")) {
    return("Male")
  }

  stop("`sex` must be 'Female' or 'Male'.", call. = FALSE)
}

sex_suffix <- function(sex) {
  if (identical(normalize_sex(sex), "Female")) "F" else "M"
}

default_mean_output_filename <- function(sex) {
  paste0("mean_", sex_suffix(sex), ".pdf")
}

default_country_output_filename <- function(sex) {
  paste0("countries_", sex_suffix(sex), ".pdf")
}

default_trajectory_output_filenames <- function(country_codes, sex) {
  stats::setNames(
    object = paste0(tolower(country_codes), "_", sex_suffix(sex), ".html"),
    nm = country_codes
  )
}

default_legacy_trajectory_output_filenames <- function(sex) {
  suffix <- sex_suffix(sex)
  stats::setNames(
    object = c(
      paste0("bel_", suffix, ".html"),
      paste0("blr_", suffix, ".html")
    ),
    nm = c("ESP", "BLR")
  )
}

approx_vec <- function(xout, x, y) {
  stats::approx(x = x, y = y, xout = xout, rule = 2)$y
}

compute_mean_quantiles <- function(data_q) {
  n_subjects <- length(data_q)
  n_times <- length(data_q[[1]])

  vapply(
    X = seq_len(n_times),
    FUN = function(time_index) {
      Reduce(
        f = `+`,
        x = lapply(seq_len(n_subjects), function(subject_index) data_q[[subject_index]][[time_index]]$y)
      ) / n_subjects
    },
    FUN.VALUE = numeric(length(data_q[[1]][[1]]$y))
  )
}

compute_mean_cdf_matrix <- function(mean_quantile_matrix, u_grid) {
  vapply(
    X = seq_len(ncol(mean_quantile_matrix)),
    FUN = function(time_index) {
      time_specific_xgrid <- seq(
        from = min(mean_quantile_matrix[, time_index]),
        to = max(mean_quantile_matrix[, time_index]),
        length.out = nrow(mean_quantile_matrix)
      )
      approx_vec(xout = time_specific_xgrid, x = mean_quantile_matrix[, time_index], y = u_grid)
    },
    FUN.VALUE = numeric(nrow(mean_quantile_matrix))
  )
}

compute_observed_transport_list <- function(data_q, mean_cdf_matrix, xgrid, common_u_grid) {
  n_subjects <- length(data_q)
  n_times <- length(data_q[[1]])

  lapply(
    X = seq_len(n_subjects),
    FUN = function(subject_index) {
      transport_matrix <- matrix(NA_real_, nrow = n_times, ncol = length(xgrid))
      for (time_index in seq_len(n_times)) {
        transport_matrix[time_index, ] <- approx_vec(
          xout = mean_cdf_matrix[, time_index],
          x = common_u_grid,
          y = data_q[[subject_index]][[time_index]]$y
        )
      }
      transport_matrix
    }
  )
}

trapz_vec <- function(x, y) {
  if (length(x) != length(y)) {
    stop("trapz_vec() requires vectors with the same length.", call. = FALSE)
  }
  if (length(x) < 2) {
    return(0)
  }
  sum(diff(x) * (y[-length(y)] + y[-1]) / 2)
}

signed_mass_for_transport <- function(transport_curve, xgrid) {
  integrand <- function(y) {
    approx_vec(xout = y, x = xgrid, y = transport_curve) - y
  }

  stats::integrate(
    f = integrand,
    lower = min(xgrid),
    upper = max(xgrid),
    stop.on.error = FALSE,
    subdivisions = 1000,
    rel.tol = .Machine$double.eps^0.5
  )$value
}

compute_signed_mass_matrix <- function(transport_list, xgrid) {
  do.call(
    what = rbind,
    args = lapply(
      X = transport_list,
      FUN = function(transport_matrix) {
        apply(transport_matrix, 1, signed_mass_for_transport, xgrid = xgrid)
      }
    )
  )
}

legacy_monotone_upper_clamp <- function(transport_curve, upper_bound) {
  adjusted_curve <- transport_curve
  adjusted_curve[adjusted_curve > upper_bound] <- upper_bound
  cummax(adjusted_curve)
}

normalize_signed_mass_matrix <- function(signed_mass_matrix) {
  normalized <- signed_mass_matrix

  for (row_index in seq_len(nrow(normalized))) {
    for (col_index in seq_len(ncol(normalized))) {
      current_range <- range(normalized, na.rm = TRUE)

      if (normalized[row_index, col_index] > 0) {
        normalized[row_index, col_index] <- normalized[row_index, col_index] / abs(current_range[2])
      }

      if (normalized[row_index, col_index] < 0) {
        normalized[row_index, col_index] <- normalized[row_index, col_index] / abs(current_range[1])
      }
    }
  }

  normalized
}

build_sign_heatmap_data <- function(signed_mass_matrix, years, country_codes, ordered_country_codes = NULL) {
  normalized_matrix <- normalize_signed_mass_matrix(signed_mass_matrix)

  if (is.null(ordered_country_codes)) {
    row_totals <- rowSums(normalized_matrix)
    ordered_index <- order(row_totals, decreasing = TRUE)
    ordered_codes <- country_codes[ordered_index]
  } else {
    ordered_codes <- ordered_country_codes
    ordered_index <- match(ordered_codes, country_codes)
    if (any(is.na(ordered_index))) {
      stop("`ordered_country_codes` contains codes that are not in the mortality data.", call. = FALSE)
    }
  }

  ordered_labels <- country_code_to_label(ordered_codes)

  plot_df <- expand.grid(
    Year = as.character(years),
    Country = ordered_labels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  plot_df$value <- as.vector(t(normalized_matrix[ordered_index, , drop = FALSE]))
  plot_df$Country <- factor(plot_df$Country, levels = rev(ordered_labels))

  list(
    plot_df = plot_df,
    ordered_country_codes = ordered_codes,
    ordered_country_labels = ordered_labels,
    normalized_matrix = normalized_matrix[ordered_index, , drop = FALSE]
  )
}

save_ggplot_pdf <- function(plot_object, filename, width, height, use_cairo = TRUE) {
  if (isTRUE(use_cairo)) {
    cairo_file <- tempfile(fileext = ".pdf")
    cairo_saved <- tryCatch(
      {
        suppressWarnings(
          ggplot2::ggsave(
            filename = cairo_file,
            plot = plot_object,
            width = width,
            height = height,
            device = grDevices::cairo_pdf
          )
        )
        file.exists(cairo_file)
      },
      error = function(...) FALSE
    )

    if (isTRUE(cairo_saved)) {
      if (file.exists(filename)) {
        unlink(filename)
      }
      file.rename(cairo_file, filename)
      return(invisible(filename))
    }
  }

  ggplot2::ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    device = grDevices::pdf
  )

  invisible(filename)
}

source_transport_module <- function(module_path) {
  module_env <- new.env(parent = globalenv())
  sys.source(module_path, envir = module_env)

  required_functions <- c(
    "fit_transport_fpca",
    "project_transport"
  )

  missing_functions <- required_functions[!vapply(required_functions, exists, logical(1), envir = module_env, inherits = FALSE)]
  if (length(missing_functions) > 0) {
    stop(
      "The proposed-method file is missing required function(s): ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }

  module_env
}

create_mortality_repro_context <- function(project_root = NULL,
                                           sex = "Female",
                                           year_start = 1983,
                                           transport_module_path = NULL) {
  sex <- normalize_sex(sex)
  root <- find_project_root(project_root)
  transport_module_path <- transport_module_path %||%
    file.path(root,  "transport_fpca_modular.R")

  data_env <- new.env(parent = emptyenv())
  raw_env <- new.env(parent = emptyenv())

  load(file.path(root, "Data", paste0("Data_", sex, ".RData")), envir = data_env)
  load(file.path(root, "Data", paste0("Data_raw_", tolower(sex), ".Rdata")), envir = raw_env)

  data_q <- data_env$Data_q
  raw_object_name <- paste0("Data_", tolower(sex))
  raw_data <- raw_env[[raw_object_name]]

  if (is.null(data_q) || is.null(raw_data)) {
    stop("Could not load the mortality data objects for sex = ", sex, call. = FALSE)
  }

  country_codes <- unique(raw_data$country)
  n_times <- length(data_q[[1]])
  years <- seq.int(from = year_start, length.out = n_times)
  u_grid <- data_q[[1]][[1]]$x
  xgrid <- seq(min(data_q[[1]][[1]]$y), max(data_q[[1]][[1]]$y), length.out = length(data_q[[1]][[1]]$y))

  mean_quantile_matrix <- compute_mean_quantiles(data_q)
  mean_cdf_matrix <- compute_mean_cdf_matrix(mean_quantile_matrix = mean_quantile_matrix, u_grid = u_grid)
  transport_list <- compute_observed_transport_list(
    data_q = data_q,
    mean_cdf_matrix = mean_cdf_matrix,
    xgrid = xgrid,
    common_u_grid = u_grid
  )

  list(
    project_root = root,
    sex = sex,
    years = years,
    data_q = data_q,
    raw_data = raw_data,
    country_codes = country_codes,
    country_labels = country_code_to_label(country_codes),
    u_grid = u_grid,
    xgrid = xgrid,
    mean_quantile_matrix = mean_quantile_matrix,
    mean_cdf_matrix = mean_cdf_matrix,
    transport_list = transport_list,
    signed_mass_matrix = compute_signed_mass_matrix(transport_list, xgrid = xgrid),
    transport_module_path = transport_module_path,
    transport_module = source_transport_module(transport_module_path)
  )
}

fit_mortality_transport_fpca <- function(context,
                                         kappa,
                                         latent,
                                         method_xi,
                                         project_baselines = FALSE,
                                         u_cap = 0.999999,
                                         z_cap = NULL,
                                         fpca_options = list(methodMuCovEst = "smooth", usergrid = TRUE)) {
  time_list <- replicate(length(context$transport_list), context$years, simplify = FALSE)

  context$transport_module$fit_transport_fpca(
    T_list = context$transport_list,
    t_list = time_list,
    xgrid = context$xgrid,
    kappa = kappa,
    latent = latent,
    methodXi = method_xi,
    fpca_optns = fpca_options,
    project = project_baselines,
    u_cap = u_cap,
    Z_cap = z_cap
  )
}

compute_cumulative_fve_labels <- function(lambda, digits = 2) {
  lambda <- lambda[is.finite(lambda)]
  cumulative_fve <- cumsum(lambda) / sum(lambda)
  sprintf("%0.*f%%", digits, 100 * cumulative_fve)
}

build_mean_eigen_plot <- function(fit,
                                  use_legacy_fve_labels = TRUE,
                                  legacy_fve_labels = c("93.76%", "98.21%", "99.53%"),
                                  eigen_sign_flip = c(-1, 1, -1),
                                  theme_name = "bw",
                                  show_panel_border = FALSE,
                                  show_axis_lines = TRUE,
                                  axis_line_color = "black",
                                  axis_line_width = 0.4,
                                  mean_color = "black",
                                  eigen_colors = c("green", "red", "blue"),
                                  mean_line_width = 0.7,
                                  eigen_line_width = 0.7,
                                  y_limits = c(-0.5, 0.5),
                                  x_label = "Calendar Year",
                                  y_label = " ",
                                  base_family = "sans",
                                  axis_title_size = 14,
                                  axis_text_size = 14,
                                  legend_title_size = 12,
                                  legend_text_size = 12,
                                  legend_position = "right") {
  assert_packages("ggplot2")

  n_components <- 3
  if (length(eigen_sign_flip) < n_components) {
    stop("Please provide three sign-flip values in `eigen_sign_flip`.", call. = FALSE)
  }

  fve_labels <- compute_cumulative_fve_labels(fit$fpca_fit$lambda, digits = 2)[seq_len(n_components)]
  if (isTRUE(use_legacy_fve_labels)) {
    fve_labels <- legacy_fve_labels
  }

  base_theme <- switch(
    EXPR = theme_name,
    bw = ggplot2::theme_bw(base_family = base_family),
    gray = ggplot2::theme_gray(base_family = base_family),
    grey = ggplot2::theme_gray(base_family = base_family),
    classic = ggplot2::theme_classic(base_family = base_family),
    minimal = ggplot2::theme_minimal(base_family = base_family),
    stop("Unsupported `theme_name`: ", theme_name, call. = FALSE)
  )

  years <- fit$workGrid
  mean_df <- data.frame(Year = years, Value = as.numeric(fit$fpca_fit$mu))
  eigen_df <- do.call(
    what = rbind,
    args = lapply(seq_len(n_components), function(component_index) {
      data.frame(
        Year = years,
        Value = eigen_sign_flip[component_index] * fit$fpca_fit$phi[, component_index],
        Label = fve_labels[component_index],
        stringsAsFactors = FALSE
      )
    })
  )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = mean_df,
      mapping = ggplot2::aes(x = Year, y = Value),
      color = mean_color,
      linewidth = mean_line_width
    ) +
    ggplot2::geom_line(
      data = eigen_df,
      mapping = ggplot2::aes(x = Year, y = Value, color = Label),
      linewidth = eigen_line_width
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(x = x_label, y = y_label, color = "FVE") +
    ggplot2::scale_color_manual(values = stats::setNames(eigen_colors, fve_labels)) +
    base_theme +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      axis.text.x = ggplot2::element_text(size = axis_text_size),
      axis.text.y = ggplot2::element_text(size = axis_text_size),
      panel.border = if (isTRUE(show_panel_border)) {
        ggplot2::element_rect(color = "black", fill = NA)
      } else {
        ggplot2::element_blank()
      },
      axis.line = if (isTRUE(show_axis_lines)) {
        ggplot2::element_line(color = axis_line_color, linewidth = axis_line_width)
      } else {
        ggplot2::element_blank()
      },
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key = ggplot2::element_rect(fill = "white")
    )
}

build_country_sign_heatmap <- function(context,
                                       ordered_country_codes = NULL,
                                       fill_low = "red",
                                       fill_high = "blue",
                                       show_legend = FALSE,
                                       x_break_years = seq(1983, 2018, by = 5),
                                       x_label = "Calendar Year",
                                       y_label = "",
                                       base_family = "sans",
                                       axis_text_x_size = 14,
                                       axis_text_y_size = 14,
                                       axis_text_angle = 90,
                                       legend_title_size = 12,
                                       legend_text_size = 12,
                                       legend_position = "right") {
  assert_packages("ggplot2")

  heatmap_data <- build_sign_heatmap_data(
    signed_mass_matrix = context$signed_mass_matrix,
    years = context$years,
    country_codes = context$country_codes,
    ordered_country_codes = ordered_country_codes
  )

  plot_object <- ggplot2::ggplot(
    data = heatmap_data$plot_df,
    mapping = ggplot2::aes(x = Year, y = Country, fill = value)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = fill_low, high = fill_high, limits = c(-1, 1)) +
    ggplot2::scale_x_discrete(breaks = as.character(x_break_years)) +
    ggplot2::labs(x = x_label, y = y_label, fill = NULL) +
    ggplot2::theme_gray(base_family = base_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = axis_text_angle,
        vjust = 0.5,
        hjust = 1,
        size = axis_text_x_size
      ),
      axis.text.y = ggplot2::element_text(size = axis_text_y_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.position = if (isTRUE(show_legend)) legend_position else "none"
    )

  list(
    plot = plot_object,
    heatmap_data = heatmap_data
  )
}

interpolate_eigenfunctions_to_observed_times <- function(fpca_fit, observed_times, component_index) {
  approx_vec(
    xout = observed_times,
    x = fpca_fit$workGrid,
    y = fpca_fit$phi[, component_index]
  )
}

reconstruct_scalar_trajectories <- function(fit,
                                            include_mean = TRUE,
                                            score_method = c("xi_est", "mc_projection"),
                                            observed_scalar_list = NULL) {
  score_method <- match.arg(score_method)

  if (identical(score_method, "xi_est")) {
    score_matrix <- fit$fpca_fit$xiEst
  } else {
    if (is.null(observed_scalar_list)) {
      stop("`observed_scalar_list` is required when `score_method = 'mc_projection'`.", call. = FALSE)
    }

    n_subjects <- length(observed_scalar_list)
    n_components <- fit$fpca_fit$selectK
    score_matrix <- matrix(0, nrow = n_subjects, ncol = n_components)

    for (subject_index in seq_len(n_subjects)) {
      for (component_index in seq_len(n_components)) {
        phi_at_observed_times <- interpolate_eigenfunctions_to_observed_times(
          fpca_fit = fit$fpca_fit,
          observed_times = fit$t_list[[subject_index]],
          component_index = component_index
        )
        score_matrix[subject_index, component_index] <- sum(
          observed_scalar_list[[subject_index]] * phi_at_observed_times
        )
      }
    }
  }

  reconstructed <- score_matrix %*% t(fit$fpca_fit$phi)
  if (isTRUE(include_mean)) {
    reconstructed <- sweep(reconstructed, 2, fit$fpca_fit$mu, `+`)
  }
  reconstructed
}

reconstruct_transports_from_scalar <- function(fit,
                                               scalar_matrix,
                                               project_rows = FALSE,
                                               project_function = NULL) {
  n_subjects <- nrow(scalar_matrix)
  n_times <- ncol(scalar_matrix)
  xgrid <- fit$xgrid
  lower_bound <- min(xgrid)
  upper_bound <- max(xgrid)

  lapply(seq_len(n_subjects), function(subject_index) {
    subject_matrix <- matrix(NA_real_, nrow = n_times, ncol = length(xgrid))
    for (time_index in seq_len(n_times)) {
      scalar_value <- scalar_matrix[subject_index, time_index]
      if (scalar_value >= 0) {
        subject_matrix[time_index, ] <- xgrid + scalar_value * (fit$T0_plus[subject_index, ] - xgrid)
      } else {
        subject_matrix[time_index, ] <- xgrid + abs(scalar_value) * (fit$T0_minus[subject_index, ] - xgrid)
      }

      if (isTRUE(project_rows)) {
        if (is.null(project_function)) {
          stop("`project_function` is required when `project_rows = TRUE`.", call. = FALSE)
        }

        subject_matrix[time_index, ] <- project_function(
          subject_matrix[time_index, ],
          lower = lower_bound,
          upper = upper_bound
        )
      }
    }
    subject_matrix
  })
}

default_trajectory_filename <- function(country_code, sex = "Female") {
  paste0(tolower(country_code), "_", sex_suffix(sex), ".html")
}

build_trajectory_widget <- function(context,
                                    reconstructed_transport_list,
                                    country_code,
                                    observed_line_color = "red",
                                    observed_line_width = 3,
                                    observed_line_opacity = 0.6,
                                    fitted_surface_color = "blue",
                                    fitted_surface_opacity = 0.6,
                                    reference_surface_color = "gray",
                                    reference_surface_opacity = 0.5,
                                    camera_center = c(x = 0, y = 0, z = -0.2),
                                    camera_eye = c(x = 1.25, y = 1.25, z = 1.25),
                                    x_axis_title = "Calendar Year",
                                    y_axis_title = "Age",
                                    z_axis_title = " ",
                                    plot_title = NULL) {
  assert_packages(c("plotly", "htmlwidgets"))

  country_index <- match(country_code, context$country_codes)
  if (is.na(country_index)) {
    stop("Unknown country code: ", country_code, call. = FALSE)
  }

  years <- context$years
  xgrid <- context$xgrid
  observed_matrix <- context$transport_list[[country_index]]
  fitted_matrix <- reconstructed_transport_list[[country_index]]

  observed_z <- sweep(observed_matrix, 2, xgrid, `-`)
  fitted_z <- sweep(t(fitted_matrix), 1, xgrid, `-`)
  zero_surface <- matrix(0, nrow = length(xgrid), ncol = length(years))

  figure <- plotly::plot_ly()
  for (time_index in seq_along(years)) {
    figure <- plotly::add_trace(
      p = figure,
      x = years[time_index],
      y = xgrid,
      z = observed_z[time_index, ],
      type = "scatter3d",
      mode = "lines",
      opacity = observed_line_opacity,
      line = list(width = observed_line_width, color = observed_line_color)
    )
  }

  scene_settings <- list(
    camera = list(
      center = as.list(camera_center),
      eye = as.list(camera_eye)
    ),
    xaxis = list(title = x_axis_title),
    yaxis = list(title = y_axis_title),
    zaxis = list(title = z_axis_title)
  )

  figure <- plotly::add_surface(
    p = figure,
    x = years,
    y = xgrid,
    z = fitted_z,
    opacity = fitted_surface_opacity,
    colorscale = list(c(0, 1), c(fitted_surface_color, fitted_surface_color))
  )
  figure <- plotly::add_surface(
    p = figure,
    x = years,
    y = xgrid,
    z = zero_surface,
    opacity = reference_surface_opacity,
    colorscale = list(c(0, 1), c(reference_surface_color, reference_surface_color))
  )
  figure <- plotly::layout(
    p = figure,
    title = plot_title,
    scene = scene_settings,
    showlegend = FALSE
  )
  figure <- plotly::hide_colorbar(figure)

  figure
}

save_widget_file <- function(widget,
                             filename,
                             selfcontained = FALSE,
                             title = "plotly") {
  assert_packages("htmlwidgets")
  htmlwidgets::saveWidget(
    widget = widget,
    file = filename,
    selfcontained = selfcontained,
    title = title
  )
}

run_mean_eigen_figure <- function(context = NULL,
                                  project_root = NULL,
                                  output_dir = NULL,
                                  output_filename = NULL,
                                  kappa = 2000,
                                  latent = TRUE,
                                  method_xi = "IN",
                                  project_baselines = FALSE,
                                  u_cap = 0.999999,
                                  z_cap = NULL,
                                  use_legacy_fve_labels = TRUE,
                                  legacy_fve_labels = c("93.76%", "98.21%", "99.53%"),
                                  eigen_sign_flip = c(-1, 1, -1),
                                  theme_name = "bw",
                                  show_panel_border = FALSE,
                                  show_axis_lines = TRUE,
                                  axis_line_color = "black",
                                  axis_line_width = 0.4,
                                  mean_color = "black",
                                  eigen_colors = c("green", "red", "blue"),
                                  mean_line_width = 0.7,
                                  eigen_line_width = 0.7,
                                  y_limits = c(-0.5, 0.5),
                                  x_label = "Calendar Year",
                                  y_label = " ",
                                  base_family = "sans",
                                  axis_title_size = 14,
                                  axis_text_size = 14,
                                  legend_title_size = 12,
                                  legend_text_size = 12,
                                  legend_position = "right",
                                  plot_width = 500 / 72,
                                  plot_height = 400 / 72,
                                  save_plot = TRUE,
                                  use_cairo_pdf = TRUE,
                                  display_plot = interactive()) {
  if (is.null(context)) {
    context <- create_mortality_repro_context(project_root = project_root, sex = "Female")
  }

  fit <- fit_mortality_transport_fpca(
    context = context,
    kappa = kappa,
    latent = latent,
    method_xi = method_xi,
    project_baselines = project_baselines,
    u_cap = u_cap,
    z_cap = z_cap
  )

  plot_object <- build_mean_eigen_plot(
    fit = fit,
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
    x_label = x_label,
    y_label = y_label,
    base_family = base_family,
    axis_title_size = axis_title_size,
    axis_text_size = axis_text_size,
    legend_title_size = legend_title_size,
    legend_text_size = legend_text_size,
    legend_position = legend_position
  )

  output_dir <- output_dir %||% file.path(context$project_root,  "output")
  output_dir <- ensure_dir(output_dir)
  output_filename <- output_filename %||% default_mean_output_filename(context$sex)
  output_path <- file.path(output_dir, output_filename)

  if (isTRUE(save_plot)) {
    save_ggplot_pdf(
      plot_object = plot_object,
      filename = output_path,
      width = plot_width,
      height = plot_height,
      use_cairo = use_cairo_pdf
    )
  }

  if (isTRUE(display_plot)) {
    print(plot_object)
  }

  invisible(list(
    plot = plot_object,
    fit = fit,
    output_path = output_path
  ))
}

run_country_sign_figure <- function(context = NULL,
                                    project_root = NULL,
                                    output_dir = NULL,
                                    output_filename = NULL,
                                    ordered_country_codes = NULL,
                                    fill_low = "red",
                                    fill_high = "blue",
                                    show_legend = FALSE,
                                    x_break_years = seq(1983, 2018, by = 5),
                                    x_label = "Calendar Year",
                                    y_label = "",
                                    base_family = "sans",
                                    axis_text_x_size = 14,
                                    axis_text_y_size = 14,
                                    axis_text_angle = 90,
                                    legend_title_size = 12,
                                    legend_text_size = 12,
                                    legend_position = "right",
                                    plot_width = 500 / 72,
                                    plot_height = 600 / 72,
                                    save_plot = TRUE,
                                    use_cairo_pdf = TRUE,
                                    display_plot = interactive()) {
  if (is.null(context)) {
    context <- create_mortality_repro_context(project_root = project_root, sex = "Female")
  }

  heatmap_result <- build_country_sign_heatmap(
    context = context,
    ordered_country_codes = ordered_country_codes,
    fill_low = fill_low,
    fill_high = fill_high,
    show_legend = show_legend,
    x_break_years = x_break_years,
    x_label = x_label,
    y_label = y_label,
    base_family = base_family,
    axis_text_x_size = axis_text_x_size,
    axis_text_y_size = axis_text_y_size,
    axis_text_angle = axis_text_angle,
    legend_title_size = legend_title_size,
    legend_text_size = legend_text_size,
    legend_position = legend_position
  )

  output_dir <- output_dir %||% file.path(context$project_root,  "output")
  output_dir <- ensure_dir(output_dir)
  output_filename <- output_filename %||% default_country_output_filename(context$sex)
  output_path <- file.path(output_dir, output_filename)

  if (isTRUE(save_plot)) {
    save_ggplot_pdf(
      plot_object = heatmap_result$plot,
      filename = output_path,
      width = plot_width,
      height = plot_height,
      use_cairo = use_cairo_pdf
    )
  }

  if (isTRUE(display_plot)) {
    print(heatmap_result$plot)
  }

  invisible(list(
    plot = heatmap_result$plot,
    heatmap_data = heatmap_result$heatmap_data,
    output_path = output_path
  ))
}

run_trajectory_html_figures <- function(context = NULL,
                                        project_root = NULL,
                                        output_dir = NULL,
                                        country_codes = NULL,
                                        output_filenames = NULL,
                                        kappa = 500,
                                        latent = FALSE,
                                        method_xi = "CE",
                                        project_baselines = FALSE,
                                        u_cap = 0.999999,
                                        z_cap = NULL,
                                        include_mean = FALSE,
                                        score_method = c("xi_est", "mc_projection"),
                                        legacy_baseline_postprocess = TRUE,
                                        project_reconstructed_rows = FALSE,
                                        observed_line_color = "red",
                                        observed_line_width = 3,
                                        observed_line_opacity = 0.6,
                                        fitted_surface_color = "blue",
                                        fitted_surface_opacity = 0.6,
                                        reference_surface_color = "gray",
                                        reference_surface_opacity = 0.5,
                                        camera_center = c(x = 0, y = 0, z = -0.2),
                                        camera_eye = c(x = 1.25, y = 1.25, z = 1.25),
                                        x_axis_title = "Calendar Year",
                                        y_axis_title = "Age",
                                        z_axis_title = " ",
                                        plot_title = NULL,
                                        selfcontained = FALSE,
                                        save_widgets = TRUE,
                                        display_widgets = interactive()) {
  score_method <- match.arg(score_method)

  if (is.null(context)) {
    context <- create_mortality_repro_context(project_root = project_root, sex = "Female")
  }

  fit <- fit_mortality_transport_fpca(
    context = context,
    kappa = kappa,
    latent = latent,
    method_xi = method_xi,
    project_baselines = project_baselines,
    u_cap = u_cap,
    z_cap = z_cap
  )

  if (isTRUE(legacy_baseline_postprocess)) {
    upper_bound <- max(context$xgrid)
    fit$T0_plus <- t(apply(fit$T0_plus, 1, legacy_monotone_upper_clamp, upper_bound = upper_bound))
    fit$T0_minus <- t(apply(fit$T0_minus, 1, legacy_monotone_upper_clamp, upper_bound = upper_bound))
  }

  observed_scalar_list <- NULL
  if (identical(score_method, "mc_projection")) {
    if (isTRUE(latent)) {
      observed_scalar_list <- lapply(fit$U_hat_obs, fit$g_inv)
    } else {
      observed_scalar_list <- fit$U_hat_obs
    }
  }

  scalar_matrix <- reconstruct_scalar_trajectories(
    fit = fit,
    include_mean = include_mean,
    score_method = score_method,
    observed_scalar_list = observed_scalar_list
  )

  reconstructed_transport_list <- reconstruct_transports_from_scalar(
    fit = fit,
    scalar_matrix = scalar_matrix,
    project_rows = project_reconstructed_rows,
    project_function = context$transport_module$project_transport
  )

  output_dir <- output_dir %||% file.path(context$project_root,  "output")
  output_dir <- ensure_dir(output_dir)

  if (is.null(output_filenames) && is.null(country_codes)) {
    output_filenames <- default_legacy_trajectory_output_filenames(context$sex)
    country_codes <- names(output_filenames)
  } else {
    country_codes <- country_codes %||% names(output_filenames)
    output_filenames <- output_filenames %||% default_trajectory_output_filenames(country_codes, context$sex)
  }

  widgets <- list()
  output_paths <- character(0)

  for (country_code in country_codes) {
    widget <- build_trajectory_widget(
      context = context,
      reconstructed_transport_list = reconstructed_transport_list,
      country_code = country_code,
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
      plot_title = plot_title
    )

    widgets[[country_code]] <- widget

    output_filename <- output_filenames[[country_code]] %||% default_trajectory_filename(country_code, sex = context$sex)
    output_path <- file.path(output_dir, output_filename)
    output_paths[country_code] <- output_path

    if (isTRUE(save_widgets)) {
      save_widget_file(
        widget = widget,
        filename = output_path,
        selfcontained = selfcontained,
        title = "plotly"
      )
    }

    if (isTRUE(display_widgets)) {
      print(widget)
    }
  }

  invisible(list(
    fit = fit,
    scalar_matrix = scalar_matrix,
    reconstructed_transport_list = reconstructed_transport_list,
    widgets = widgets,
    output_paths = output_paths
  ))
}
